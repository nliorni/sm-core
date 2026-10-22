# Configuration

One YAML file per project; `config/config.yaml` is the annotated template. Every key is
validated at start-up against `workflow/schemas/config.schema.yaml`, which also holds each
key's default — a config only has to state what differs. A wrong type, an unknown
strandedness or a missing FASTQ stops the run before anything executes, with a message
naming the key.

Relative paths are relative to the pipeline directory (`pipelines/rnaseq/`), which
`run.py` uses as the working directory wherever it is launched from. Absolute paths are
recommended for real projects.

## Samples

Exactly one of `samplesheet` or `samples` is required.

### `samplesheet`

Path to a CSV (or `.tsv`) in the nf-core/rnaseq format, so an existing nf-core samplesheet
works unchanged:

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample ID; becomes part of every output path, so no spaces or `/`. |
| `fastq_1`, `fastq_2` | yes | R1 and R2. Relative paths are relative to the samplesheet file. Both are required: the pipeline is paired-end only. |
| `strandedness` | no | `auto`, `forward`, `reverse` or `unstranded`; default: the top-level `strandedness`. |
| anything else | no | Sample metadata (condition, batch, donor…). Carried into `samples.tsv` and the `colData` of every `.rds`, so a DESeq2 design formula can use it directly. |

Several rows with the same `sample` are lanes or runs of one library: they are
concatenated in row order before trimming. They must agree on `strandedness` and on every
metadata column.

### `samples`, `unit`, `extension`

The sm-core family convention, for FASTQ named `{prefix}_{unit[0]}.{extension}` and
`{prefix}_{unit[1]}.{extension}`:

```yaml
unit: ["R1_001", "R2_001"]
extension: "fastq.gz"
samples:
  CTRL_1: "/data/CTRL_1_S1_L001"                              # one pair
  CTRL_2: ["/data/CTRL_2_S2_L001", "/data/CTRL_2_S2_L002"]    # lanes merged
```

No metadata columns this way; use a samplesheet when the matrices should carry them.

### `strandedness` — default `auto`

Library strandedness for samples without a samplesheet value. `auto` lets each branch infer
it from its own data and is recommended; declaring it is useful as a check — a declaration
the data contradicts is reported as `MISMATCH` — and for the Salmon branch, which then
quantifies with that library type instead of detecting one. See
[methods.md](methods.md#strandedness).

| Value | Typical kits |
|---|---|
| `reverse` | dUTP/second-strand-marking: Illumina TruSeq Stranded mRNA/Total RNA, Illumina Stranded mRNA, NEBNext Ultra II Directional, KAPA Stranded |
| `forward` | Takara SMARTer Stranded Pico v1, some ligation-based kits |
| `unstranded` | Illumina TruSeq (non-stranded), most SMART-seq / full-length single-cell-style protocols |

### `output_dir` — required

All outputs, logs and benchmarks; also the default reference location.

### `notify_email` — optional

An e-mail at the end of the run, success or failure, via the system `mail`/`mailx`.

## `reference`

| Key | Default | Meaning |
|---|---|---|
| `genome` | — (required) | Genome FASTA, plain or gzipped. See [resources.md](resources.md). |
| `gtf` | — (required) | Annotation GTF from the same source and release as the genome. |
| `spikein_fasta` | `null` | Spike-in sequences (e.g. `ERCC92.fa`), appended to the genome with one gene each. |
| `dir` | `<output_dir>/reference/` | Where the derived reference and indices are built. Share it between projects that use the same genome + GTF. |
| `star_index` | `null` | Use a STAR index built elsewhere — from exactly this genome + GTF (+ spike-ins) and STAR 2.7.11b. |
| `salmon_index` | `null` | Same, for Salmon 1.10.3. |
| `star_sjdb_overhang` | `100` | STAR's `--sjdbOverhang`: ideally read length − 1. |
| `salmon_kmer` | `31` | Salmon index k-mer size. Lower it (e.g. 23) only for reads shorter than ~75 bp. |
| `threads` | `16` | Threads for index building. |
| `star_index_mem_mb` | `48000` | Memory requested for the STAR index build (human: peaks at 38 GB). |
| `salmon_index_mem_mb` | `32000` | Memory requested for the Salmon index build. |

## `conda_envs` — optional

Name → environment prefix or YAML, to use site-maintained environments instead of the
bundled ones. Names: `fastp`, `star`, `salmon`, `samtools`, `gffread`, `rnaseqc`,
`multiqc`, `python`, `bioconductor`. Unlisted names fall back to
`workflow/envs/<name>.yaml`.

## `trimming`

| Key | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Trim with fastp and give both branches the trimmed reads. With `false`, fastp still runs for read QC but writes no reads. |
| `min_length` | `20` | Drop pairs with either mate shorter than this after trimming. |
| `extra` | `""` | Extra fastp arguments, e.g. `--trim_front1 3 --trim_front2 3`. |
| `threads` | `8` | |
| `mem_mb` | `4000` | |

## `star`

| Key | Default | Meaning |
|---|---|---|
| `threads` | `16` | |
| `mem_mb` | `40000` | Human: peaks at 30 GB, the index loaded in memory. |
| `two_pass` | `true` | `--twopassMode Basic`: re-align with the junctions found in a first pass. |
| `sort_threads` | `4` | samtools sort threads for the BAM. |
| `sort_mem_per_thread` | `"1G"` | samtools sort memory per thread. |
| `extra` | `""` | Appended last to the STAR command, so it overrides the defaults in `workflow/rules/star.smk`. |

## `salmon`

| Key | Default | Meaning |
|---|---|---|
| `threads` | `8` | |
| `mem_mb` | `32000` | Human: peaks at 27 GB, since the decoy-aware index holds the genome. |
| `extra` | `"--gcBias --seqBias"` | Extra `salmon quant` arguments. Drop the bias flags for 3′-tag libraries (QuantSeq), where they do not apply. |

## `qc`

| Key | Default | Meaning |
|---|---|---|
| `stranded_threshold` | `0.8` | Minimum share of fragments on one strand to call a library stranded. |
| `unstranded_threshold` | `0.1` | Maximum forward–reverse difference to call it unstranded. |
| `fail_on_strandedness_mismatch` | `false` | Stop the run when a declared strandedness is contradicted, instead of flagging it. |
| `rnaseqc` | `true` | Run RNA-SeQC 2 in the STAR branch. |
| `rnaseqc_extra` | `""` | Extra RNA-SeQC arguments. |
| `multiqc_extra` | `""` | Extra MultiQC arguments. |
