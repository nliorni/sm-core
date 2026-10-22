# Usage

The launcher's flags (local vs PBS runs, retries, locking, resuming, DAG export) are shared
by every sm-core pipeline and documented once, in [the suite's docs/usage.md](../../../docs/usage.md).
This page covers what is specific to rnaseq.

## Targets

```bash
./run.py --list-workflows
```

| Target | Runs | Delivers |
|---|---|---|
| `all` | both branches | everything below, plus `results/qc/strandedness.tsv` |
| `star` | fastp → STAR → GeneCounts, RNA-SeQC | `results/counts/star/`, BAMs, `results/qc/star/multiqc_report.html` |
| `salmon` | fastp → Salmon → tximport | `results/counts/salmon/`, `results/qc/salmon/multiqc_report.html` |
| `reference` | reference preparation and both indices | `reference.dir` |
| `trim` | fastp | trimmed reads (normally temporary) |

Targets combine: `-w star salmon` is `all` without the comparison table. A target builds
whatever it needs that does not exist yet. For example, `-w salmon` on a fresh project
builds the Salmon index but not STAR's.

## A typical project

```bash
cp config/config.yaml /projects/P123/config_P123.yaml
#   edit: output_dir, samplesheet, reference.*  (and reference.dir, to share indices)

./run.py -w all -c /projects/P123/config_P123.yaml -q 32 -n      # check
./run.py -w all -c /projects/P123/config_P123.yaml -cl -qu workq -j 30
```

Then open `results/qc/star/multiqc_report.html` and `results/qc/salmon/multiqc_report.html`.
Check the strandedness tables first ([outputs.md](outputs.md#qc-what-to-look-at-first)).
After that, load a matrix:

```r
library(DESeq2)
se  <- readRDS("results/counts/salmon/gene_counts.rds")
dds <- DESeqDataSet(se, design = ~ condition)     # condition: a samplesheet column
dds <- DESeq(dds)
```

## Samplesheets

```csv
sample,fastq_1,fastq_2,strandedness,condition,batch
CTRL_1,/data/run1/CTRL_1_R1.fastq.gz,/data/run1/CTRL_1_R2.fastq.gz,auto,control,1
CTRL_1,/data/run2/CTRL_1_R1.fastq.gz,/data/run2/CTRL_1_R2.fastq.gz,auto,control,1
TREAT_1,/data/run1/TREAT_1_R1.fastq.gz,/data/run1/TREAT_1_R2.fastq.gz,auto,treated,1
```

- Rows with the same `sample` are merged, in the order they appear: use this for lanes,
  flowcells and top-up runs of the **same library**. Different libraries of the same
  biological sample are technical replicates. Give them separate IDs and decide later
  whether to collapse them (`DESeq2::collapseReplicates`).
- Metadata columns travel into every matrix's `colData`. Name them as you will want them
  in the design formula.
- nf-core/rnaseq samplesheets work unchanged. Rows without `fastq_2` (single-end) are
  refused: this version is paired-end only.

## Adding samples to a finished project

Add the rows and rerun the same command. Per-sample jobs run only for the new samples.
Every matrix and report is rebuilt, because each one spans all samples. If the pipeline
code has changed since the first run, Snakemake would also rerun the old samples;
`--rerun-triggers mtime` prevents that ([the suite's docs/usage.md](../../../docs/usage.md#resuming-and-re-running)).

## Resources

Defaults suit human and mouse on a typical HPC node. Use `benchmarks/<tool>/<sample>.tsv`
(wall time and peak memory, `max_rss`, per job) to tune `threads` and `mem_mb` per tool in
the config. On PBS, a job killed for memory is retried automatically with 2× and then 3×
the request.

| Rule | Threads | Memory (default) |
|---|---|---|
| `star_index` | `reference.threads` (16) | `reference.star_index_mem_mb` (48 GB) |
| `salmon_index` | `reference.threads` (16) | `reference.salmon_index_mem_mb` (32 GB) |
| `fastp` | `trimming.threads` (8) | 4 GB |
| `star_align` | `star.threads` (16) | `star.mem_mb` (40 GB) |
| `salmon_quant` | `salmon.threads` (8) | `salmon.mem_mb` (32 GB) |
| `rnaseqc` | 1 | 8 GB |

## Troubleshooting

**"declare samples with exactly one of `samplesheet:` or `samples:`"**: the config sets
both, or neither.

**"FASTQ not found"**: relative paths in a samplesheet are relative to the samplesheet
file. Relative paths in the config are relative to `pipelines/rnaseq/`.

**"No transcripts left: do the GTF's sequence names match the genome's?"**: a GENCODE
GTF (`chr1`) with an Ensembl genome (`1`), or the reverse. See
[resources.md](resources.md#the-one-rule-genome-and-gtf-from-the-same-release).

**A strandedness `MISMATCH`**: the declared value is wrong, almost always. Set the sample to
`auto`, or to the inferred value. For the STAR branch, the counts were read from the
declared column, so rerun after fixing it.

**STAR killed / "std::bad_alloc"**: not enough memory for the index. Raise `star.mem_mb`, or
on PBS let the retry do it (`-rt`).

**Very low Salmon mapping rate with good STAR mapping**: usually adapter read-through in
untrimmed reads (`trimming.enabled: false`), or a declared library type that is wrong. Check
`results/salmon/<sample>/lib_format_counts.json`.
