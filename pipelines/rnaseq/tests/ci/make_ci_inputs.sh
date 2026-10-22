#!/usr/bin/env bash
## ============================================================================
## Creates the empty placeholder files that tests/ci/config_ci*.yaml point at.
##
## A Snakemake dry-run checks that every input it cannot build itself *exists*;
## it never opens them. Empty files are therefore enough to resolve the whole
## DAG of every workflow target, which is what CI checks on each push -- no
## reference, no reads, no downloads. The real end-to-end check is the yeast
## tutorial (tests/tutorial/yeast), which CI also runs.
## ============================================================================
set -euo pipefail

CI_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${CI_DIR}/refs" "${CI_DIR}/fastq" "${CI_DIR}/generated"

touch "${CI_DIR}/refs/genome.fa.gz" "${CI_DIR}/refs/genes.gtf.gz" "${CI_DIR}/refs/spikein.fa"
for f in S1_L1 S1_L2 S2 S3; do
  touch "${CI_DIR}/fastq/${f}_R1.fastq.gz" "${CI_DIR}/fastq/${f}_R2.fastq.gz"
done

# Variants derived from the main CI config, so they cannot drift from it.
sed -e 's|^output_dir: .*|output_dir: "tests/ci/output_notrim/"|' \
    "${CI_DIR}/config_ci.yaml" > "${CI_DIR}/generated/config_ci_notrim.yaml"
printf '\ntrimming:\n  enabled: false\n' >> "${CI_DIR}/generated/config_ci_notrim.yaml"

sed -e 's|^output_dir: .*|output_dir: "tests/ci/output_prefixes/"|' \
    -e '/^samplesheet:/d' \
    "${CI_DIR}/config_ci.yaml" > "${CI_DIR}/generated/config_ci_prefixes.yaml"
cat >> "${CI_DIR}/generated/config_ci_prefixes.yaml" <<'YAML'

unit: ["R1", "R2"]
extension: "fastq.gz"
samples:
  S1: ["tests/ci/fastq/S1_L1", "tests/ci/fastq/S1_L2"]
  S2: "tests/ci/fastq/S2"
YAML

echo "CI placeholder inputs created under ${CI_DIR}"
