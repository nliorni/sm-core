#!/usr/bin/env bash
## ============================================================================
## Download a matched genome + annotation for sm-core/rnaseq, checksum-verified,
## and print the `reference:` block for the project config.
##
##   get_reference.sh --dest /refs/GRCh38_gencode50                      # human, GENCODE 50
##   get_reference.sh --species mouse --dest /refs/GRCm39_gencodeM39     # mouse, GENCODE M39
##   get_reference.sh --source ensembl --species danio_rerio --dest /refs/GRCz11_ens116
##   get_reference.sh --spikein ercc --dest ...                          # + ERCC92 spike-ins
##   get_reference.sh --check ...                                        # only check the URLs
##
## Options
##   --source   gencode (default; human and mouse) | ensembl (any Ensembl species)
##   --species  human | mouse          for gencode (default human)
##              <ensembl name>         for ensembl, e.g. homo_sapiens, danio_rerio
##   --release  GENCODE release (default 50 human, M39 mouse) or Ensembl release
##              (default 116). Check https://www.gencodegenes.org/ for newer ones.
##   --spikein  ercc                   also fetch the ERCC92 sequences + concentrations
##   --dest     directory to download into (created; required)
##   --check    resolve every URL without downloading
##   --https    use https:// for EBI and Ensembl (automatic when a proxy is set)
##
## What it fetches, and why those files (docs/resources.md has the full story):
##   - the PRIMARY ASSEMBLY genome: chromosomes + unplaced scaffolds, no alt
##     haplotypes, so no reads are lost to multimapping between a chromosome and
##     its alternates;
##   - the COMPREHENSIVE annotation on that same assembly, from the same
##     release, so every sequence name and transcript model matches the genome.
##
## Every file is checked against the publisher's own checksums (GENCODE MD5SUMS,
## Ensembl CHECKSUMS) and gzip-tested; a file that fails is deleted, not kept.
## Re-running resumes: files already present and verified are skipped.
## A MANIFEST.tsv (file, URL, checksum, date) records exactly what was fetched.
##
## Needs curl, gzip, md5sum, sum, unzip. EBI and Ensembl are used over ftp://,
## which EBI serves about ten times faster than https://. Behind an HTTP proxy
## (http_proxy / https_proxy set, as on many clusters) ftp:// directory listings
## do not work, so https:// is used instead.
## ============================================================================
set -euo pipefail

SOURCE=gencode
SPECIES=human
RELEASE=""
SPIKEIN=""
DEST=""
CHECK=0
PROTO=ftp
[[ -n "${https_proxy:-}${HTTPS_PROXY:-}${http_proxy:-}${HTTP_PROXY:-}${ftp_proxy:-}${FTP_PROXY:-}" ]] && PROTO=https

die() { echo "ERROR: $*" >&2; exit 1; }

while (( $# )); do
  case "$1" in
    --source)  SOURCE="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --release) RELEASE="$2"; shift 2 ;;
    --spikein) SPIKEIN="$2"; shift 2 ;;
    --dest)    DEST="$2"; shift 2 ;;
    --check)   CHECK=1; shift ;;
    --https)   PROTO=https; shift ;;
    -h|--help) sed -n '3,/^## ====/{/^## ====/!p}' "$0" | sed 's/^## \{0,1\}//'; exit 0 ;;
    *) die "unknown option $1 (see --help)" ;;
  esac
done
[[ -n "${DEST}" || "${CHECK}" == 1 ]] || die "--dest is required (see --help)"
[[ -z "${SPIKEIN}" || "${SPIKEIN}" == ercc ]] || die "--spikein: only 'ercc' is known"

ERCC_ZIP="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"
ERCC_CONC="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt"

list_dir() {
  # File names in a remote directory: an FTP name listing, or the links of an
  # https:// index page.
  if [[ "${PROTO}" == ftp ]]; then
    curl -sfl --max-time 60 "$1/"
  else
    curl -sf --max-time 60 "$1/" | grep -oE 'href="[^"?/][^"]*"' | sed 's/^href="//; s/"$//'
  fi
}

## ---------------------------------------------------------------------------
## Resolve the two URLs and how to verify them
## ---------------------------------------------------------------------------
case "${SOURCE}" in
  gencode)
    case "${SPECIES}" in
      human) RELEASE="${RELEASE:-50}"; ASM=GRCh38; DIR=Gencode_human; TAG="v${RELEASE}" ;;
      mouse) RELEASE="${RELEASE:-M39}"; ASM=GRCm39; DIR=Gencode_mouse; TAG="v${RELEASE}" ;;
      *) die "--source gencode covers human and mouse; use --source ensembl for ${SPECIES}" ;;
    esac
    BASE="${PROTO}://ftp.ebi.ac.uk/pub/databases/gencode/${DIR}/release_${RELEASE}"
    GENOME="${ASM}.primary_assembly.genome.fa.gz"
    GTF="gencode.${TAG}.primary_assembly.annotation.gtf.gz"
    GENOME_URL="${BASE}/${GENOME}"
    GTF_URL="${BASE}/${GTF}"
    SUMS_URL="${BASE}/MD5SUMS"
    SUM_KIND=md5
    ;;
  ensembl)
    RELEASE="${RELEASE:-116}"
    BASE="${PROTO}://ftp.ensembl.org/pub/release-${RELEASE}"
    FASTA_DIR="${BASE}/fasta/${SPECIES}/dna"
    GTF_DIR="${BASE}/gtf/${SPECIES}"
    listing="$(list_dir "${FASTA_DIR}")" || die "no Ensembl ${RELEASE} FASTA for '${SPECIES}' (${FASTA_DIR})"
    # Unmasked primary assembly if the species has one (it has alt haplotypes
    # to leave out), toplevel otherwise. Never the soft/hard-masked dna_sm/dna_rm.
    GENOME="$(grep -E '\.dna\.primary_assembly\.fa\.gz$' <<< "${listing}" || true)"
    [[ -n "${GENOME}" ]] || GENOME="$(grep -E '\.dna\.toplevel\.fa\.gz$' <<< "${listing}" || true)"
    [[ -n "${GENOME}" ]] || die "no unmasked genome FASTA in ${FASTA_DIR}"
    GTF="$(list_dir "${GTF_DIR}" | grep -E "\.${RELEASE}\.gtf\.gz$" | grep -vE '\.(chr|abinitio|chr_patch_hapl_scaff)\.' || true)"
    [[ -n "${GTF}" ]] || die "no Ensembl ${RELEASE} GTF for '${SPECIES}' (${GTF_DIR})"
    GENOME_URL="${FASTA_DIR}/${GENOME}"
    GTF_URL="${GTF_DIR}/${GTF}"
    SUMS_URL="${FASTA_DIR}/CHECKSUMS ${GTF_DIR}/CHECKSUMS"
    SUM_KIND=sum
    ;;
  *) die "--source must be gencode or ensembl" ;;
esac

if (( CHECK )); then
  for u in "${GENOME_URL}" "${GTF_URL}" ${SPIKEIN:+"${ERCC_ZIP}"}; do
    size="$(curl -sIL --max-time 60 "${u}" | tr -d '\r' | awk 'tolower($1) == "content-length:" { s = $2 } END { print s }')"
    [[ -n "${size}" ]] || die "not found: ${u}"
    if (( size >= 1000000 )); then printf '%6d MB  %s\n' $(( size / 1000000 )) "${u}"
    else printf '%6d KB  %s\n' $(( size / 1000 )) "${u}"; fi
  done
  exit 0
fi

mkdir -p "${DEST}"
DEST="$(cd "${DEST}" && pwd)"
cd "${DEST}"

## ---------------------------------------------------------------------------
## Download + verify
## ---------------------------------------------------------------------------
expected_sum() {
  # The publisher's checksum line for one file.
  local file="$1" url
  for url in ${SUMS_URL}; do
    curl -sfL --max-time 60 "${url}" 2>/dev/null | awk -v f="${file}" '$NF == f'
  done | head -1
}

verified() {
  local file="$1" line
  [[ -s "${file}" ]] || return 1
  line="$(expected_sum "${file}")"
  [[ -n "${line}" ]] || { echo "  WARNING: no published checksum for ${file}; gzip test only" >&2; gzip -t "${file}"; return; }
  if [[ "${SUM_KIND}" == md5 ]]; then
    [[ "$(md5sum < "${file}" | cut -d' ' -f1)" == "$(awk '{ print $1 }' <<< "${line}")" ]]
  else
    # Ensembl CHECKSUMS hold BSD `sum` output: checksum, blocks, name.
    [[ "$(sum "${file}" | awk '{ print $1, $2 }')" == "$(awk '{ print $1, $2 }' <<< "${line}")" ]]
  fi
}

fetch() {
  local url="$1" file
  file="$(basename "${url}")"
  if verified "${file}" 2>/dev/null; then
    echo "  already present and verified: ${file}"
    return
  fi
  echo "  fetching ${file}"
  curl -sfL --retry 3 --connect-timeout 30 -C - -o "${file}.part" "${url}" || { rm -f "${file}.part"; die "download failed: ${url}"; }
  mv "${file}.part" "${file}"
  if ! verified "${file}" || ! gzip -t "${file}"; then
    rm -f "${file}"
    die "${file} failed its checksum; deleted. Re-run to fetch it again."
  fi
  printf '%s\t%s\t%s\t%s\n' "${file}" "${url}" "$(md5sum < "${file}" | cut -d' ' -f1)" "$(date +%F)" >> MANIFEST.tsv
}

[[ -s MANIFEST.tsv ]] || printf 'file\turl\tmd5\tdownloaded\n' > MANIFEST.tsv
echo "== ${SOURCE} ${SPECIES} release ${RELEASE} -> ${DEST}"
fetch "${GENOME_URL}"
fetch "${GTF_URL}"

if [[ "${SPIKEIN}" == ercc ]]; then
  if [[ ! -s ERCC92.fa ]]; then
    echo "  fetching ERCC92 spike-ins"
    tmp="$(mktemp -d)"
    curl -sfL --retry 3 -o "${tmp}/ERCC92.zip" "${ERCC_ZIP}" || die "download failed: ${ERCC_ZIP}"
    unzip -q -o "${tmp}/ERCC92.zip" ERCC92.fa -d "${tmp}"
    [[ "$(grep -c '^>' "${tmp}/ERCC92.fa")" == 92 ]] || die "ERCC92.fa does not hold 92 sequences"
    mv "${tmp}/ERCC92.fa" .
    rm -rf "${tmp}"
    printf 'ERCC92.fa\t%s\t%s\t%s\n' "${ERCC_ZIP}" "$(md5sum < ERCC92.fa | cut -d' ' -f1)" "$(date +%F)" >> MANIFEST.tsv
  fi
  # Thermo Fisher ships the concentration table with Windows line endings.
  [[ -s ERCC_Controls_Analysis.txt ]] || curl -sfL "${ERCC_CONC}" | tr -d '\r' > ERCC_Controls_Analysis.txt
fi

## ---------------------------------------------------------------------------
## Sanity check: the GTF's sequence names must be the genome's
## ---------------------------------------------------------------------------
genome_seqs="$(zcat "$(basename "${GENOME_URL}")" | grep '^>' | sed 's/^>//; s/[[:space:]].*//' | sort -u)"
gtf_seqs="$(zcat "$(basename "${GTF_URL}")" | grep -v '^#' | cut -f1 | sort -u)"
missing="$(comm -13 <(echo "${genome_seqs}") <(echo "${gtf_seqs}") | head -5)"
[[ -z "${missing}" ]] || die "GTF sequences absent from the genome (first few): ${missing//$'\n'/ }"
echo "  checked: all $(wc -l <<< "${gtf_seqs}") GTF sequence names are in the genome ($(wc -l <<< "${genome_seqs}") sequences)"

cat <<EOF

Done. Put this in the project config (SECTION 2):

reference:
  genome: "${DEST}/$(basename "${GENOME_URL}")"
  gtf: "${DEST}/$(basename "${GTF_URL}")"
EOF
[[ "${SPIKEIN}" == ercc ]] && echo "  spikein_fasta: \"${DEST}/ERCC92.fa\""
cat <<EOF
  # Build the indices once, here, and share them between projects:
  dir: "${DEST}/sm-core-rnaseq"

Then build the indices once:  ./run.py -w reference -c <config> -q 16
EOF
