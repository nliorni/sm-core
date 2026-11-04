#!/usr/bin/env bash
## ============================================================================
## Creates the empty placeholder files that tests/ci/config_ci*.yaml point at.
##
## A Snakemake dry-run checks that every input it cannot build itself *exists*;
## it never opens them. Empty files are therefore enough to resolve the whole
## DAG of every workflow target, which is what CI checks on each push -- no
## reference, no reads, no downloads.
## ============================================================================
set -euo pipefail

CI_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${CI_DIR}/refs" "${CI_DIR}/fastq" "${CI_DIR}/generated"

touch "${CI_DIR}/refs/genome.fa.gz" "${CI_DIR}/refs/genes.gtf.gz" "${CI_DIR}/refs/blacklist.bed" \
      "${CI_DIR}/refs/sites.vcf.gz" "${CI_DIR}/refs/rna.bam"
for f in C1_L1 C1_L2 C2 T1 T2; do
  touch "${CI_DIR}/fastq/${f}_R1.fastq.gz" "${CI_DIR}/fastq/${f}_R2.fastq.gz"
done

# Variants derived from the main CI config, so they cannot drift from it.
# Identity check, with a BAM of another assay.
cp "${CI_DIR}/config_ci.yaml" "${CI_DIR}/generated/config_ci_identity.yaml"
cat >> "${CI_DIR}/generated/config_ci_identity.yaml" <<'YAML'

qc:
  somalier_sites: "tests/ci/refs/sites.vcf.gz"
  identity_bams: ["tests/ci/refs/rna.bam"]
  identity_prefix: "rna_"
YAML

# No blacklist, no groups, no trimming, no differential accessibility.
sed -e '/^da:/,$d' -e '/blacklist:/d' -e 's/\[condition\]/[]/' "${CI_DIR}/config_ci.yaml" \
    > "${CI_DIR}/generated/config_ci_minimal.yaml"
printf '\ntrimming:\n  enabled: false\n' >> "${CI_DIR}/generated/config_ci_minimal.yaml"

# A group column the samplesheet lacks: must be refused at start-up.
sed -e 's/consensus_group_by: \[condition\]/consensus_group_by: [cell_line]/' \
    "${CI_DIR}/config_ci.yaml" > "${CI_DIR}/generated/config_ci_bad.yaml"

# Every consensus group a single sample, with the default p-value cut-off:
# must warn.
sed -e 's/consensus_group_by: \[condition\]/consensus_group_by: [condition, replicate]/' \
    "${CI_DIR}/config_ci.yaml" > "${CI_DIR}/generated/config_ci_unreplicated.yaml"

# A contrast level no sample has: must be refused at start-up.
sed -e 's/\[condition, treat, ctrl\]/[condition, treated, ctrl]/' \
    "${CI_DIR}/config_ci.yaml" > "${CI_DIR}/generated/config_ci_bad_da.yaml"

echo "CI placeholder inputs created under ${CI_DIR}"
