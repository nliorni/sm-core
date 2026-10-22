#!/usr/bin/env bash
## ============================================================================
## Downloads the HUMAN tutorial dataset for the rnaseq pipeline: SEQC/MAQC-III.
##
## Reads: the published raw FASTQ of the SEQC consortium's Illumina HiSeq 2000
## runs at BGI (SRA study SRP025982, GEO GSE47792), 2x100 bp paired-end:
##
##   A = Universal Human Reference RNA (UHRR, 10 cancer cell lines) + ERCC Mix 1
##   B = Human Brain Reference RNA     (HBRR)                       + ERCC Mix 2
##
## Three independent libraries of each (A_1..A_3, B_1..B_3), each sequenced on
## two flowcells -- so every sample arrives as two FASTQ pairs, the way real
## multi-lane projects do, and the pipeline has to merge them.
##
## Two published truth sets come along for the validation step in README.md:
##   - the ERCC concentration table: 92 synthetic spike-ins added at KNOWN
##     Mix 1 / Mix 2 ratios (4, 1, 2/3, 1/2), and
##   - MAQC-I TaqMan qPCR (GEO GSE5350, platform GPL4097): ~1,000 genes measured
##     in A and B, 4 replicates each -- the classic reference for A-vs-B fold
##     changes.
##
## And the reference: GENCODE's GRCh38 primary assembly with the matching
## GENCODE primary-assembly annotation, plus the ERCC92 sequences and
## annotation, which the tutorial config adds to both as spike-ins
## (reference.spikein_fasta / spikein_gtf).
##
## Reads come from one of two mirrors of the same runs:
##   aws  SRA's open-data bucket on AWS, as .sra files converted with
##        fasterq-dump -- needs sra-tools (fasterq-dump, vdb-validate) on PATH;
##        typically ~100x faster than ENA from Europe and the US.
##   ena  ENA's FASTQ, checked against ENA's MD5s -- needs only curl.
## SOURCE=auto (default) uses aws when sra-tools is installed, ena otherwise.
## pigz is used for compression when present, gzip otherwise.
##
## ~12 GB of FASTQ + ~1 GB of reference.
##
##   ./download_tutorial_data.sh                    # everything
##   SKIP_REFERENCE=1 ./download_tutorial_data.sh   # reads + truth only
##   N_PAIRS=2000000 ./download_tutorial_data.sh    # first 2M pairs per run
##   SOURCE=ena ./download_tutorial_data.sh         # force the ENA mirror
## ============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_DIR="${SCRIPT_DIR}/fastq"
TRUTH_DIR="${SCRIPT_DIR}/truth"
REF_DIR="${SCRIPT_DIR}/reference"

ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq"
AWS="https://sra-pub-run-odp.s3.amazonaws.com/sra"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
GENCODE_RELEASE="${GENCODE_RELEASE:-50}"
# Plain FTP: EBI serves the same files ~10x faster over ftp:// than https://.
GENCODE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
ERCC_ZIP="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"
ERCC_CONC="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt"
TAQMAN="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4097/soft/GPL4097_family.soft.gz"

# 0 = whole runs (default; ~12M pairs per sample once both flowcells are
# merged). A positive value keeps only the first N pairs of each run: reads sit
# in sequencer order, so a leading slice is an unbiased sample of the library,
# and R1/R2 stay paired because both are cut at the same count.
N_PAIRS="${N_PAIRS:-0}"
SKIP_REFERENCE="${SKIP_REFERENCE:-0}"
SOURCE="${SOURCE:-auto}"
THREADS="${THREADS:-8}"

if [[ "${SOURCE}" == auto ]]; then
  if command -v fasterq-dump > /dev/null && command -v vdb-validate > /dev/null; then
    SOURCE=aws
  else
    SOURCE=ena
  fi
fi
GZIP_CMD="gzip -c"
command -v pigz > /dev/null && GZIP_CMD="pigz -c -p ${THREADS}"

# sample  run (flowcell AC0AYTACXX, lane 1)  run (flowcell BC0AYPACXX, lane 1)
RUNS=(
  "SEQC_A_1 SRR896663 SRR896664"
  "SEQC_A_2 SRR896679 SRR896680"
  "SEQC_A_3 SRR896695 SRR896696"
  "SEQC_B_1 SRR896743 SRR896744"
  "SEQC_B_2 SRR896759 SRR896760"
  "SEQC_B_3 SRR896775 SRR896776"
)

mkdir -p "${FASTQ_DIR}" "${TRUTH_DIR}" "${REF_DIR}"

# Anything that fails is recorded and the script exits non-zero at the end. A
# download that half-works must not look like success: the pipeline run that
# follows would fail later and further from the cause.
FAILURES=()

ena_field() {
  # One field of a run's ENA record, e.g. fastq_md5 ("md5_1;md5_2") or
  # read_count (the number of read pairs).
  curl -s "${ENA_API}?accession=$1&result=read_run&fields=$2&format=tsv" \
    | awk 'NR == 2 { print $2 }'
}

fetch_aws() {
  # Download one run's .sra from AWS, check it, convert to a gzipped FASTQ
  # pair. Complete means vdb-validate passes and both mates hold exactly
  # ENA's read count (or N_PAIRS).
  local run="$1" tmp
  tmp="$(mktemp -d "${FASTQ_DIR}/.${run}.XXXX")"
  echo "  fetching ${run}.sra from AWS..."
  if ! curl -sfL --retry 3 --connect-timeout 30 -o "${tmp}/${run}.sra" "${AWS}/${run}/${run}" \
     || ! vdb-validate "${tmp}/${run}.sra" > /dev/null 2>&1; then
    rm -rf "${tmp}"
    return 1
  fi
  fasterq-dump --split-files --threads "${THREADS}" --temp "${tmp}" -O "${tmp}" "${tmp}/${run}.sra" \
    > /dev/null 2>&1 || { rm -rf "${tmp}"; return 1; }
  local expected n1 n2 mate
  expected="$(ena_field "${run}" read_count)"
  (( N_PAIRS > 0 )) && expected="${N_PAIRS}"
  for mate in 1 2; do
    (( N_PAIRS > 0 )) && head -n $(( N_PAIRS * 4 )) "${tmp}/${run}_${mate}.fastq" > "${tmp}/cut" \
      && mv "${tmp}/cut" "${tmp}/${run}_${mate}.fastq"
  done
  n1=$(( $(wc -l < "${tmp}/${run}_1.fastq") / 4 ))
  n2=$(( $(wc -l < "${tmp}/${run}_2.fastq") / 4 ))
  if [[ "${n1}" != "${expected}" || "${n2}" != "${expected}" ]]; then
    echo "  ${run}: ${n1}/${n2} pairs, expected ${expected}" >&2
    rm -rf "${tmp}"
    return 1
  fi
  for mate in 1 2; do
    ${GZIP_CMD} "${tmp}/${run}_${mate}.fastq" > "${tmp}/${run}_${mate}.fastq.gz"
    mv "${tmp}/${run}_${mate}.fastq.gz" "${FASTQ_DIR}/${run}_${mate}.fastq.gz"
  done
  rm -rf "${tmp}"
}

fastq_path() {
  # SRR896663 -> SRR896/SRR896663 (ENA's layout for 9-character run IDs).
  local run="$1"
  echo "${run:0:6}/${run}"
}

fetch_whole() {
  local url="$1" out="$2" md5="$3"
  local attempt
  for attempt in 1 2 3; do
    echo "  fetching $(basename "${out}") (attempt ${attempt})..."
    if curl -sfL --retry 3 --connect-timeout 30 -C - -o "${out}.part" "${url}" \
       && [[ "$(md5sum < "${out}.part" | cut -d' ' -f1)" == "${md5}" ]]; then
      mv "${out}.part" "${out}"
      return 0
    fi
    # A resumed download with a wrong checksum is corrupt from the start.
    rm -f "${out}.part"
  done
  return 1
}

fetch_slice() {
  local url="$1" out="$2"
  # head closes the pipe once it has enough, so curl takes SIGPIPE: expected,
  # hence pipefail off for this one pipeline. Complete means exactly N records.
  set +o pipefail
  curl -sL --connect-timeout 30 "${url}" | zcat 2>/dev/null \
    | head -n $(( N_PAIRS * 4 )) | gzip -c > "${out}.part"
  set -o pipefail
  local n=$(( $(zcat "${out}.part" | wc -l) / 4 ))
  if (( n == N_PAIRS )); then
    mv "${out}.part" "${out}"
    return 0
  fi
  echo "  ${out}: got ${n} of ${N_PAIRS} pairs" >&2
  rm -f "${out}.part"
  return 1
}

echo "== Reads (SEQC, BGI, HiSeq 2000, 2x100 bp) from ${SOURCE} =="
for row in "${RUNS[@]}"; do
  read -r sample run1 run2 <<< "${row}"
  for run in "${run1}" "${run2}"; do
    if [[ -s "${FASTQ_DIR}/${run}_1.fastq.gz" && -s "${FASTQ_DIR}/${run}_2.fastq.gz" ]]; then
      echo "  already present: ${run}"
      continue
    fi
    if [[ "${SOURCE}" == aws ]]; then
      fetch_aws "${run}" || FAILURES+=("${run}")
      continue
    fi
    md5s=""
    (( N_PAIRS == 0 )) && md5s="$(ena_field "${run}" fastq_md5)"
    for mate in 1 2; do
      out="${FASTQ_DIR}/${run}_${mate}.fastq.gz"
      [[ -s "${out}" ]] && continue
      url="${ENA}/$(fastq_path "${run}")/${run}_${mate}.fastq.gz"
      if (( N_PAIRS == 0 )); then
        md5="$(cut -d';' -f"${mate}" <<< "${md5s}")"
        if [[ -z "${md5}" ]]; then
          FAILURES+=("${run}_${mate}: no MD5 from ENA")
          continue
        fi
        fetch_whole "${url}" "${out}" "${md5}" || FAILURES+=("${run}_${mate}")
      else
        fetch_slice "${url}" "${out}" || FAILURES+=("${run}_${mate}")
      fi
    done
  done
done

echo "== Truth sets =="
if [[ ! -s "${TRUTH_DIR}/ERCC_Controls_Analysis.txt" ]]; then
  # ThermoFisher ships this with Windows line endings.
  curl -sfL "${ERCC_CONC}" | tr -d '\r' > "${TRUTH_DIR}/ERCC_Controls_Analysis.txt" \
    || FAILURES+=("ERCC concentration table")
fi
if [[ ! -s "${TRUTH_DIR}/GPL4097_family.soft.gz" ]]; then
  curl -sfL -o "${TRUTH_DIR}/GPL4097_family.soft.gz" "${TAQMAN}" \
    || FAILURES+=("MAQC TaqMan (GPL4097)")
fi

echo "== Reference: GENCODE ${GENCODE_RELEASE} GRCh38 primary assembly + ERCC92 =="
if [[ ! -s "${REF_DIR}/ERCC92.fa" ]]; then
  tmp="$(mktemp -d)"
  if curl -sfL -o "${tmp}/ERCC92.zip" "${ERCC_ZIP}" && unzip -q -o "${tmp}/ERCC92.zip" -d "${tmp}"; then
    mv "${tmp}/ERCC92.fa" "${tmp}/ERCC92.gtf" "${REF_DIR}/"
  else
    FAILURES+=("ERCC92.zip")
  fi
  rm -rf "${tmp}"
fi
if [[ "${SKIP_REFERENCE}" != 1 ]]; then
  for f in GRCh38.primary_assembly.genome.fa.gz "gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"; do
    if [[ -s "${REF_DIR}/${f}" ]]; then
      echo "  already present: ${f}"
      continue
    fi
    echo "  fetching ${f}..."
    if curl -sfL --retry 3 -C - -o "${REF_DIR}/${f}.part" "${GENCODE}/${f}" && gzip -t "${REF_DIR}/${f}.part"; then
      mv "${REF_DIR}/${f}.part" "${REF_DIR}/${f}"
    else
      rm -f "${REF_DIR}/${f}.part"
      FAILURES+=("${f}")
    fi
  done
fi

if (( ${#FAILURES[@]} )); then
  echo >&2
  echo "FAILED: ${FAILURES[*]}" >&2
  echo "Re-run this script: files already downloaded are kept." >&2
  exit 1
fi
echo
echo "Done. Next, from pipelines/rnaseq/:"
echo "  ./run.py -w all -c tests/tutorial/seqc/config_tutorial.yaml -q 32"
