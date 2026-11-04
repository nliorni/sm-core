# Usage

The launcher's flags (local vs cluster runs, retries, locking, resuming, DAG export) are
shared by every sm-core pipeline and documented once, in [the suite's docs/usage.md](../../../docs/usage.md).
This page covers what is specific to atacseq.

## Targets

```bash
./run.py --list-workflows
```

| Target | Runs | Delivers |
|---|---|---|
| `all` | everything | BAMs, peaks, consensus peaks and counts, bigWigs, QC, MultiQC (and the identity check, with `qc.somalier_sites`; and `da`, when `da.analyses` is set) |
| `da` | PyDESeq2 per contrast, plus whatever it needs upstream | `results/da/` ([outputs.md](outputs.md#differential-accessibility--w-da)) |
| `reference` | reference preparation and the bowtie2 index | `reference.dir` |
| `trim` | fastp | trimmed reads (normally temporary) |

## A typical project

```bash
cp config/config.yaml /projects/P123/config_P123_atac.yaml
#   edit: output_dir, samplesheet, reference.* (and reference.dir, to share the index),
#         peaks.consensus_group_by and bigwig.group_by

./run.py -w all -c /projects/P123/config_P123_atac.yaml -q 32 -n      # check
./run.py -w all -c /projects/P123/config_P123_atac.yaml -q 32
```

Then open `results/qc/qc_metrics.png` and `results/qc/multiqc_report.html`
([outputs.md](outputs.md#qc-what-to-look-at-first)). Once the libraries pass, add the
contrasts to the config's `da` section and run the differential accessibility; only the
`da` jobs run, the rest is reused:

```bash
./run.py -w da -c /projects/P123/config_P123_atac.yaml -q 16
```

## Samplesheets and groups

```csv
sample,fastq_1,fastq_2,cell_line,condition,replicate
LINE1_CTRL_1,/data/L1C1_L001_R1.fastq.gz,/data/L1C1_L001_R2.fastq.gz,LINE1,CTRL,1
LINE1_CTRL_1,/data/L1C1_L002_R1.fastq.gz,/data/L1C1_L002_R2.fastq.gz,LINE1,CTRL,1
LINE1_TREAT_1,/data/L1T1_R1.fastq.gz,/data/L1T1_R2.fastq.gz,LINE1,TREAT,1
```

- Rows with the same `sample` are lanes or runs of **one library** and are merged.
- Name the samples as in the project's RNA-seq samplesheet, if there is one. The identity
  check can then pair each ATAC library with its RNA library by name.
- `peaks.consensus_group_by: [cell_line, condition]` makes one replicate group per line
  and condition. With 3 replicates per group, `min_replicates: 2` keeps a region called in
  2 of the 3 replicates of any group.

## Checking sample identity across assays

When the same cells were also profiled by RNA-seq on the same genome assembly, add the
RNA BAMs:

```yaml
qc:
  somalier_sites: /refs/GRCh38_gencode50/sites.hg38.vcf.gz
  identity_bams:
    - /projects/P123/rnaseq/results/star/LINE1_CTRL_1/LINE1_CTRL_1.bam
    # ...
  identity_prefix: "rna_"
  somalier_extra: "--min-ab 0.2"
```

In `results/qc/identity/somalier.html`, `LINE1_CTRL_1` and `rna_LINE1_CTRL_1` should have
relatedness near 1, and every pair from different cell lines near 0.

`results/qc/identity/identity_matches.tsv` does not rely on names: for every library, its
most related ATAC library (`best`, `relatedness`, `sites`) and its most related RNA library
(`best_other`, ...). `match` / `match_other` name the best one when the relatedness is at
least 0.8 (same cells), and are empty otherwise. An ATAC library whose RNA match is a
library of another cell line is a swap on one side. The other assay's samples tell which:
check how they match among themselves and against the ATAC lines. Matches on fewer than
100 sites are noisy; RNA BAMs give fewer informative sites than ATAC ones, as only
expressed SNPs count.

## Resources

Defaults suit human and mouse on a typical HPC node. `benchmarks/<tool>/<sample>.tsv` has
the wall time and peak memory of each job.

| Rule | Threads | Memory (default) |
|---|---|---|
| `bowtie2_index` | `reference.threads` (16) | `reference.bowtie2_index_mem_mb` (32 GB) |
| `fastp` | `trimming.threads` (8) | 4 GB |
| `bowtie2_align` | `align.threads` (16) | `align.mem_mb` (16 GB) |
| `markdup` | 1 | `align.markdup_mem_mb` (16 GB) |
| `atac_qc` | 1 | `qc.mem_mb` (12 GB) |
| `macs3` | 1 | `peaks.mem_mb` (16 GB) |
| `bigwig` | `bigwig.threads` (8) | `bigwig.mem_mb` (8 GB) |

## Troubleshooting

**"reference.blacklist shares no sequence name with the genome"**: a UCSC-named blacklist
(`chr1`) with an Ensembl genome (`1`). Use `resources/get_atac_resources.sh --ensembl`.

**"'<column>' is not a samplesheet column"**: a `consensus_group_by` or `group_by` column
is missing from the samplesheet (names are case-sensitive).

**"no consensus peak"**: every per-sample peak file is empty or tiny, or
`min_replicates` exceeds what the groups can meet. Check `results/qc/qc_metrics.png`
first: a library with almost no signal has almost no peaks.

**Picard killed / "OutOfMemoryError"**: raise `align.markdup_mem_mb`. The JVM gets 80% of
it.

**`atac_qc` killed**: it keeps every fragment position of one library in memory, about
24 bytes per fragment. 100 M fragments need about 3 GB, plus Python's overhead. Raise
`qc.mem_mb`.
