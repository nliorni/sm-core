#!/usr/bin/env bash
## ============================================================================
## Downloads the YEAST smoke-test dataset for the rnaseq pipeline: nf-core's
## own rnaseq test data (S. cerevisiae, GEO GSE110004), so a run here can be
## set beside an nf-core/rnaseq test run of the same reads.
##
## Three paired-end samples, 50,000 read pairs per run, against a 230 kb slice
## of the genome (chromosome I) and its annotation -- which also carries a
## record on a contig that is not in the FASTA (`not_in_genome`), on purpose:
## the pipeline must drop it rather than fail. WT_REP1 was sequenced in two
## runs, so lane merging is exercised too.
##
## The ERCC92 spike-in sequences are added to the reference too (there are no
## spike-ins in these reads, so they count zero): it puts the spike-in path --
## genome + GTF + both indices -- under test on every CI run.
##
## This is a smoke test, not a validation: it proves every rule runs and
## connects, in minutes. The human tutorial (../seqc) is the validation.
##
## Needs only curl and unzip.
## ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_DIR="${SCRIPT_DIR}/fastq"
REF_DIR="${SCRIPT_DIR}/reference"

DATA="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004"
# The reference nf-core/rnaseq's own `test` profile uses, pinned to a commit.
REF="https://raw.githubusercontent.com/nf-core/test-datasets/7f1614baeb0ddf66e60be78c3d9fa55440465ac8/reference"

mkdir -p "${FASTQ_DIR}" "${REF_DIR}"

fetch() {
  local url="$1" out="$2"
  if [[ -s "${out}" ]]; then
    echo "  already present: $(basename "${out}")"
    return
  fi
  echo "  fetching $(basename "${out}")"
  curl -sfL --retry 3 -o "${out}.part" "${url}"
  mv "${out}.part" "${out}"
}

for run in SRR6357070 SRR6357071 SRR6357072 SRR6357076; do
  for mate in 1 2; do
    fetch "${DATA}/${run}_${mate}.fastq.gz" "${FASTQ_DIR}/${run}_${mate}.fastq.gz"
  done
done
fetch "${REF}/genome.fasta" "${REF_DIR}/genome.fasta"
fetch "${REF}/genes.gtf" "${REF_DIR}/genes.gtf"
if [[ ! -s "${REF_DIR}/ERCC92.fa" ]]; then
  fetch "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip" "${REF_DIR}/ERCC92.zip"
  unzip -q -o "${REF_DIR}/ERCC92.zip" ERCC92.fa -d "${REF_DIR}"
  rm -f "${REF_DIR}/ERCC92.zip"
fi

echo
echo "Done. Next, from pipelines/rnaseq/:"
echo "  ./run.py -w all -c tests/tutorial/yeast/config_tutorial.yaml -q 8"
