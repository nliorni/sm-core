#!/usr/bin/env bash
## ============================================================================
## Checks the yeast tutorial's deliverables after
##   ./run.py -w all -c tests/tutorial/yeast/config_tutorial.yaml
## Run from pipelines/rnaseq/. Used by CI; exits non-zero on the first failure.
##
## These are structural checks plus the facts this dataset is known for:
## GSE110004 is a reverse-stranded (dUTP) library, so both branches must infer
## `reverse`, agree with each other, and agree with the samples declared
## `reverse` in the samplesheet.
## ============================================================================
set -euo pipefail

OUT="tests/tutorial/yeast/output/results"
SAMPLES=(WT_REP1 WT_REP2 RAP1_IAA_30M_REP1)
fail() { echo "FAIL: $*" >&2; exit 1; }

for f in counts/star/gene_counts.tsv counts/star/gene_counts.rds counts/star/gene_counts_summary.tsv \
         counts/salmon/gene_counts.tsv counts/salmon/gene_counts_length_scaled.tsv \
         counts/salmon/gene_tpm.tsv counts/salmon/gene_lengths.tsv \
         counts/salmon/transcript_counts.tsv counts/salmon/gene_counts.rds \
         qc/star/multiqc_report.html qc/salmon/multiqc_report.html qc/strandedness.tsv; do
  [[ -s "${OUT}/${f}" ]] || fail "missing or empty: ${f}"
done
for s in "${SAMPLES[@]}"; do
  [[ -s "${OUT}/star/${s}/${s}.bam.bai" ]] || fail "no BAM index for ${s}"
  [[ -s "${OUT}/salmon/${s}/quant.sf" ]] || fail "no quant.sf for ${s}"
done

for branch in star salmon; do
  header="$(head -1 "${OUT}/counts/${branch}/gene_counts.tsv")"
  [[ "${header}" == $'gene_id\tgene_name\tWT_REP1\tWT_REP2\tRAP1_IAA_30M_REP1' ]] \
    || fail "${branch} matrix header: ${header}"
  # Every sample must have counts: a wrong strandedness column gives ~0.
  awk -F'\t' -v b="${branch}" 'NR > 1 { for (i = 3; i <= NF; i++) s[i] += $i }
    END { for (i = 3; i <= 5; i++) if (s[i] < 1000) { print b ": sample column " i " sums to " s[i]; exit 1 } }' \
    "${OUT}/counts/${branch}/gene_counts.tsv" || fail "${branch} counts too low"
done

awk -F'\t' 'NR > 1 && ($3 != "reverse" || $4 != "reverse" || $7 != "yes") { print; bad = 1 }
  END { exit bad }' "${OUT}/qc/strandedness.tsv" || fail "strandedness should be reverse in both branches"

grep -q $'\tMISMATCH\t' "${OUT}"/qc/*/strandedness.tsv && fail "unexpected strandedness MISMATCH"

# Every section of each report made it in: a malformed custom-content header
# makes MultiQC drop that section silently, without failing.
for branch in star salmon; do
  data="${OUT}/qc/${branch}/multiqc_report_data"
  for module in fastp general_stats; do
    [[ -s "${data}/multiqc_${module}.txt" ]] || fail "${branch} report lacks ${module}"
  done
  for section in "strandedness_${branch}" "biotypes_${branch}"; do
    grep -q "${section}" "${data}/multiqc_data.json" || fail "${branch} report lacks the ${section} section"
  done
done
for f in star rna_seqc; do
  [[ -s "${OUT}/qc/star/multiqc_report_data/multiqc_${f}.txt" ]] || fail "star report lacks ${f}"
done
[[ -s "${OUT}/qc/salmon/multiqc_report_data/multiqc_salmon.txt" ]] || fail "salmon report lacks salmon"

# The 92 ERCC spike-ins are genes of both matrices (with zero counts here).
for branch in star salmon; do
  n=$(grep -c '^ERCC-' "${OUT}/counts/${branch}/gene_counts.tsv" || true)
  [[ "${n}" == 92 ]] || fail "${branch} matrix has ${n} ERCC spike-ins, expected 92"
done

# The record on `not_in_genome` must have been dropped, not failed on.
grep -q "not_in_genome" tests/tutorial/yeast/output/logs/reference/prepare_annotation.log \
  || fail "GTF filtering did not report the off-genome record"

echo "yeast tutorial: all checks passed"
