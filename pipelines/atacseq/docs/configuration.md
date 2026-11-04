# Configuration

One YAML file per project, validated against `workflow/schemas/config.schema.yaml` at
start-up: a wrong type, an unknown value or a missing required key stops the run before
any job starts, naming the key. Every key not set takes its schema default.
`config/config.yaml` is the annotated template; this page explains each key.

Relative paths are relative to the pipeline directory (`pipelines/atacseq/`); use absolute
paths for real projects.

## Samples

Exactly one of `samplesheet:` or `samples:`.

### `samplesheet`

An nf-core-style CSV (or TSV) with the columns `sample,fastq_1,fastq_2`, plus any number of
metadata columns:

```
sample,fastq_1,fastq_2,cell_line,condition,replicate
LINE1_CTRL_1,/data/L1C1_L001_R1.fastq.gz,/data/L1C1_L001_R2.fastq.gz,LINE1,CTRL,1
LINE1_CTRL_1,/data/L1C1_L002_R1.fastq.gz,/data/L1C1_L002_R2.fastq.gz,LINE1,CTRL,1
```

Rows that share a `sample` are lanes or runs of one library, and are merged before
trimming. The metadata columns define the groups of `peaks.consensus_group_by` and
`bigwig.group_by`, and are written next to the count matrix
(`results/counts/samples.tsv`). Relative FASTQ paths are read relative to the samplesheet.
A `strandedness` column, if present (rnaseq sheets have one), is ignored.

### `samples`, `unit`, `extension`

The sm-core FASTQ-prefix convention, as in rnaseq: `samples: {ID: prefix or [prefixes]}`,
with reads at `{prefix}_{unit[0]}.{extension}` and `{prefix}_{unit[1]}.{extension}`. There
are no metadata columns this way, so groups cannot be used.

### `output_dir` (required)

`results/`, `logs/` and `benchmarks/` are created inside it, and `reference/` too unless
`reference.dir` is set.

### `notify_email` (optional)

An e-mail is sent through the system `mail` command on success or failure.

## `reference`

| Key | Default | Meaning |
|---|---|---|
| `genome` | required | Genome FASTA, plain or gzipped. GENCODE's primary assembly is recommended. |
| `gtf` | required | The matching GTF. It is used for TSSs (TSS enrichment) and for peak annotation. |
| `blacklist` | `null` | ENCODE blacklist BED, with the genome's sequence names. Peaks overlapping it are dropped, and it is left out of the bigWigs. Strongly recommended. A blacklist that shares no sequence name with the genome (e.g. `1` vs `chr1`) stops the run. |
| `dir` | `<output_dir>/reference` | Where the derived reference is built. Point several projects at one shared directory to build the bowtie2 index once, but never let two runs build it at the same time. |
| `bowtie2_index` | `null` | A prebuilt bowtie2 index **prefix**, from exactly this genome. |
| `mito_name` | `auto` | The mitochondrial sequence. `auto` takes `chrM`, `MT`, `chrMT` or `M`, whichever the genome has. |
| `effective_genome_size` | `null` | MACS3 `-g`. By default, the genome's non-N length. |
| `threads`, `bowtie2_index_mem_mb` | 16, 32000 | bowtie2-build resources. |

## `conda_envs` (optional)

`{name: prefix or YAML}` to use pre-built environments instead of the bundled ones. The
names are fastp, bowtie2, samtools, picard, macs3, subread, deeptools, somalier, atacqc,
multiqc and python.

## `trimming`

Shared with rnaseq (`smcore/rules/reads.smk`). fastp always runs for its read QC.

| Key | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Write trimmed reads and align those. With `false`, bowtie2 reads the raw FASTQ. |
| `min_length` | 20 | Drop a pair when either mate is shorter after trimming. |
| `extra` | `""` | Extra fastp options. |
| `threads`, `mem_mb` | 8, 4000 | |

The adapter is found from the overlap of the two mates (`--detect_adapter_for_pe`), which
works for the Nextera adapter of ATAC libraries without configuration.

## `align`

| Key | Default | Meaning |
|---|---|---|
| `max_fragment` | 2000 | bowtie2 `-X`: the longest fragment that still counts as a proper pair. |
| `min_mapq` | 30 | Pairs below it are dropped (ENCODE). |
| `remove_duplicates` | `true` | Remove Picard-marked duplicates from the final BAM. With `false` they are kept, and also counted in peaks and signal. |
| `extra` | `""` | Extra bowtie2 options. |
| `threads`, `mem_mb` | 16, 16000 | bowtie2 (+ samtools sort). |
| `markdup_mem_mb` | 16000 | Picard MarkDuplicates. 80% of it becomes the JVM heap. |

## `peaks`

| Key | Default | Meaning |
|---|---|---|
| `cutoff` | `pvalue` | MACS3 significance cut-off per sample: `pvalue` (`-p`) or `qvalue` (`-q`). |
| `threshold` | 0.01 | Its value. |
| `extra` | `""` | Extra MACS3 callpeak options, e.g. `--call-summits`. |
| `mem_mb` | 16000 | |
| `consensus_group_by` | `[]` | Samplesheet columns whose combined values make a replicate group, e.g. `[cell_line, condition]`. `[]` makes all samples one group. |
| `min_replicates` | 2 | A merged region is kept when at least this many samples of one group called a peak in it. Capped at each group's size, so a group of one still counts. |

The default cut-off, p < 0.01 per sample, is the ENCODE ATAC-seq pipeline's. Each
sample's peaks are permissive, and the replicate consensus is the reproducibility
filter: a region must be called in `min_replicates` samples of one group. On the
[tutorial](../tests/tutorial/encode/README.md) this finds 93% and 91% of ENCODE's
reproducible (IDR) peaks; q < 0.01 found about 81% and 73%.

When a group has a single sample, nothing filters its peaks, and the run warns. For such
designs set `cutoff: qvalue`, `threshold: 0.05` (MACS3's default), or group samples so
that every group has replicates.

## `bigwig`

| Key | Default | Meaning |
|---|---|---|
| `bin_size` | 10 | Bin size in bp. |
| `group_by` | `[]` | Samplesheet columns. One bigWig, averaged over its samples, per group. `[]` gives per-sample tracks only. |
| `threads`, `mem_mb` | 8, 8000 | |

## `annotation`

| Key | Default | Meaning |
|---|---|---|
| `gene_types` | `[protein_coding, lncRNA]` | Gene types (GENCODE `gene_type` / Ensembl `gene_biotype`) that peaks are annotated to. |
| `promoter_upstream`, `promoter_downstream` | 1000, 1000 | A peak is a "promoter" peak of every gene whose window around the TSS (in the gene's direction) it overlaps. |

The order of `gene_types` also breaks ties: when a peak is equally close to the TSSs of several genes, the earlier type wins, then a named gene over one known only by its ID.

## `qc`

| Key | Default | Meaning |
|---|---|---|
| `threads`, `mem_mb` | 2, 12000 | `atac_qc`, which holds every fragment position of one library in memory: about 8 bytes × 3 per fragment. |
| `somalier_sites` | `null` | The somalier sites VCF for the genome build. Setting it turns on the identity check. |
| `identity_bams` | `[]` | BAMs of other assays of the same cells, on the same assembly, added to the identity check (e.g. rnaseq's `results/star/<sample>/<sample>.bam`). |
| `identity_prefix` | `"ext_"` | Prefix of those BAMs' sample names in the report. Their name is the BAM's read-group `SM` (for rnaseq, the sample ID), or the file name without `.bam` when the BAM has no read group. A pair whose names are equal up to the prefix is expected to be the same cells, and somalier flags it if not. |
| `somalier_extra` | `""` | Options for `somalier relate`. `--min-ab 0.2` is somalier's advice when RNA-seq BAMs are included. |
| `multiqc_extra` | `""` | Extra MultiQC options. |

## `da` — differential accessibility (`-w da`)

PyDESeq2 on the consensus peak counts, for every contrast listed. Nothing runs until
`analyses` has an entry; once it has, `-w all` includes it. The `analyses` format is the
same as rnaseq's `de`, so the ATAC and RNA analyses of a project can share their designs.
Methods are in [methods.md](methods.md#differential-accessibility), files in
[outputs.md](outputs.md#differential-accessibility--w-da).

| Key | Default | Meaning |
|---|---|---|
| `normalization` | `median_of_ratios` | Size factors: `median_of_ratios` (DESeq2's, over the peaks: assumes most peaks do not change) or `library_size` (fragments in each final BAM: robust to a global opening or closing of chromatin, but shifted by differences in FRiP). Both are computed and compared in every run; see [methods.md](methods.md#normalisation). |
| `alpha` | `0.05` | Adjusted p-value cut-off: independent filtering and the significant lists. |
| `min_abs_log2fc` | `1.0` | Minimum \|log2 fold change\| for a peak to enter the significant lists. The Wald test itself is against 0. |
| `min_count` | `10` | Pre-filter: keep a peak with at least this many fragments in at least as many samples as the smallest group of the contrast's factor. |
| `shrink` | `true` | Add shrunken log2 fold changes (PyDESeq2 `lfc_shrink`), used for the MA plot. |
| `label_top` | `15` | Peaks labelled with their nearest gene on the volcano plot. |
| `heatmap_top` | `50` | Peaks in the heatmap of each contrast. |
| `plot_formats` | `[png, pdf]` | Every figure in each format. |
| `threads`, `mem_mb` | `4`, `16000` | Per job. |
| `analyses` | `[]` | The models, below. |

Each analysis:

| Key | Meaning |
|---|---|
| `name` | Unique; becomes a folder name. |
| `design` | A formula over samplesheet columns: `"~ condition"`, or `"~ replicate + condition"` for a paired design. |
| `contrasts` | `[[factor, test, reference], ...]`: log2 fold change = test / reference. |
| `split_by` | Optional samplesheet column: one model per value (e.g. per cell line), so each value's dispersions are estimated on its own samples. Without it, all samples share one model and the design can adjust for the column instead (`"~ cell_line + condition"`). |
| `exclude` | Optional list of samples left out. |

```yaml
da:
  analyses:
    - name: per_line
      split_by: cell_line
      design: "~ condition"
      contrasts:
        - [condition, treated, control]
```

Mistakes stop the run before any job starts, naming the analysis: a column the samplesheet
lacks, a level no sample has (in a split), a contrast factor outside the design, or a model
with no replicates left to estimate dispersions.
