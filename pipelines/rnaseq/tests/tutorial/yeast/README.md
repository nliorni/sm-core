# Yeast smoke test

The fastest way to see every rule of both branches run: a few minutes once the conda
environments exist. CI runs this tutorial on every push.

It is a **smoke test, not a validation**. It proves that the pipeline runs end to end
and that its outputs connect. The [human SEQC tutorial](../seqc/README.md) is the
validation.

## Dataset

nf-core/rnaseq's own test data, so a run here can be set beside an nf-core/rnaseq
`-profile test` run of the same reads: *S. cerevisiae*, GEO
[GSE110004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110004), 50,000 read
pairs per run, against a 230 kb slice of the genome (chromosome I) and its annotation.

| Sample | Runs | Declared strandedness | Exercises |
|---|---|---|---|
| WT_REP1 | SRR6357070 + SRR6357071 | `auto` | lane merging; inference |
| WT_REP2 | SRR6357072 | `reverse` | declared value checked against the data |
| RAP1_IAA_30M_REP1 | SRR6357076 | `reverse` | |

The 92 ERCC spike-in sequences are added to the reference (`reference.spikein_fasta`).
There are no spike-ins in these reads, so they count zero, but the spike-in path through
both indices is tested on every CI run.

The annotation contains, on purpose, a record on a sequence that is not in the genome
(`not_in_genome`). The pipeline must drop it and say so. nf-core's test data carries the
same trap.

nf-core subsampled the reads to chromosome I, so mapping rates stay high (about 80%)
even with a one-chromosome reference.

## Run

```bash
cd pipelines/rnaseq
tests/tutorial/yeast/download_tutorial_data.sh      # ~20 MB, curl only
./run.py -w all -c tests/tutorial/yeast/config_tutorial.yaml -q 4
tests/tutorial/yeast/check_outputs.sh
```

`check_outputs.sh` checks that every deliverable of both branches exists and that both
matrices have counts for every sample. It also checks what is known about this library:
it is reverse-stranded, so both branches must infer `reverse`, agree with each other, and
raise no `MISMATCH` for the samples declared `reverse`.

## Expected results

Measured on this pipeline's own run (2026-10-21), for comparison with yours:

| Sample | Read pairs | STAR uniquely mapped | Salmon mapping rate | STAR counted | Salmon counts | STAR vs Salmon, per gene (Pearson r, log2) |
|---|---|---|---|---|---|---|
| WT_REP1 (2 runs) | 95,473 | 88.9% | 80.0% | 79,311 | 76,361 | 0.94 |
| WT_REP2 | 47,961 | 89.4% | 79.6% | 40,079 | 38,163 | 0.99 |
| RAP1_IAA_30M_REP1 | 47,839 | 90.3% | 80.7% | 40,682 | 38,613 | 0.99 |

- **Strandedness:** both branches infer `reverse` for every sample, independently. The
  forward share of fragments is 0.10–0.15 for STAR and 0.03–0.10 for Salmon, and the
  samples declared `reverse` pass. STAR's share is higher because STAR counts per gene,
  and a read antisense to one gene can still fall inside a gene on the other strand.
- **Lane merging:** WT_REP1 has twice the reads of the others, because both its runs went
  in.
- **RNA-SeQC:** 84–86% exonic, 8–9% intergenic. That is plausible with only chromosome I
  annotated.
- **Annotation:** 124 transcripts kept, 1 record dropped (`not_in_genome`), listed in
  `logs/reference/prepare_annotation.log`.

The counts are small, so the numbers move with tool versions in the last digits. The
strandedness calls and the agreement between branches should not.
