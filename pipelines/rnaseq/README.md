<h1 align="center">sm-core/rnaseq</h1>

<p align="center">
  <em>Paired-end RNA-seq, from raw reads to count matrices — twice, independently.</em>
</p>

<p align="center">
  <img alt="Snakemake" src="https://img.shields.io/badge/snakemake-9.16.0-3C4B9A">
  <img alt="STAR" src="https://img.shields.io/badge/STAR-2.7.11b-E4572E">
  <img alt="Salmon" src="https://img.shields.io/badge/Salmon-1.10.3-E4572E">
  <img alt="Reference" src="https://img.shields.io/badge/reference-GRCh38%20%7C%20GENCODE%2050-3C4B9A">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-green">
</p>

---

**sm-core/rnaseq** is the RNA-seq pipeline of the [sm-core](../../README.md) suite. It takes
paired-end short-read FASTQ files to gene- and transcript-level count matrices, ready for
DESeq2, edgeR or limma. Two independent quantification branches run from the same trimmed
reads: **STAR** aligns to the genome and counts reads per gene, and **Salmon** quantifies
transcripts with bias correction, which tximport then summarises to genes. Each branch
delivers its own count matrix and its own QC report. Neither reads the other's output, so
where the two disagree, you can see it.

## Features

| Target | What it delivers | Main tools |
|---|---|---|
| `star` | Sorted, indexed genome BAM; **integer gene count matrix**; RNA-SeQC metrics; MultiQC report | fastp, STAR, samtools, RNA-SeQC 2, MultiQC |
| `salmon` | Transcript quantification; **gene and transcript count matrices** (counts, length-scaled counts, TPM, length offsets); MultiQC report | fastp, Salmon, tximport, MultiQC |
| `all` | Both, plus a STAR-vs-Salmon comparison of the strandedness each inferred | |
| `reference` | Filtered annotation, transcriptome, STAR and decoy-aware Salmon indices, built once | gffread, STAR, Salmon |
| `trim` | Adapter- and quality-trimmed reads only | fastp |

- **Two matrices, for two purposes.** STAR's counts are integers, traceable fragment by
  fragment to the BAM. Salmon's estimates also use multimapping reads, correct for GC and
  priming biases, and come with the transcript-length offsets that benchmarks favour for
  differential expression. [methods.md](docs/methods.md#why-two-branches) explains when each
  one wins.
- **DESeq2-ready.** Each branch writes a `SummarizedExperiment` (`.rds`) carrying the
  samplesheet's metadata as `colData` and gene names, biotypes and chromosomes as
  `rowData`. `DESeqDataSet(se, ~ condition)` is the next line of your analysis.
- **Strandedness is checked, not assumed.** Each branch infers it from its own evidence.
  A declared value that the data contradicts is flagged, or stops the run if you prefer.
- **nf-core-compatible samplesheets**, with lanes and runs merged and any extra column
  carried through as metadata. The sm-core FASTQ-prefix convention works too.
- **The reference is built from what you give it.** The annotation is filtered to the
  genome, spike-ins are added as genes, and Salmon's transcriptome is cut from the same
  genome by the same GTF. Every index is built by the pipeline, and can be shared across
  projects.
- **QC that answers the usual questions**: read quality, mapping, rRNA and mitochondrial
  share, genomic-DNA contamination (intergenic rate), RNA degradation (3′ bias),
  strandedness. One MultiQC report per branch.
- **Laptop to HPC.** Runs locally, or submits one cluster job per rule (`-cl`), with automatic
  retries that request more memory and walltime on each attempt.
- **A human tutorial with a published truth set**: SEQC/MAQC-III reference RNAs, with ERCC
  spike-ins at known ratios and TaqMan qPCR, validating both matrices. A yeast smoke test
  runs in minutes.

## Quick start

```bash
# 1. Launcher environment (Snakemake; analysis software is built per rule)
conda env create -f ../../environment.yml
conda activate sm-core

# 2. A genome + matching annotation, checksum-verified (human GENCODE 50 by default;
#    --species mouse, --source ensembl --species <name>, --spikein ercc).
#    It prints the `reference:` block for the config.
resources/get_reference.sh --dest /path/to/references/GRCh38_gencode50

# 3. Copy the template; set output_dir, samplesheet, and paste the reference block
cp config/config.yaml config_my_project.yaml

# 4. Dry-run, then run both branches (the first run also builds the indices, ~2-3 h human)
./run.py -w all -c config_my_project.yaml -q 32 -n
./run.py -w all -c config_my_project.yaml -q 32

# ...or only one branch
./run.py -w salmon -c config_my_project.yaml -q 32

# ...or one cluster job per rule
./run.py -w all -c config_my_project.yaml -cl -qu workq -j 20
```

From the repository root, `./run.py rnaseq -w all -c …` is the same command.

## How it works

```
paired-end FASTQ (lanes/runs merged)
  │
  └─ fastp ......... adapter + quality trimming, read QC
        │
        ├─ star ...... STAR two-pass ──→ <sample>.bam ──→ RNA-SeQC 2
        │              └─ GeneCounts (strand from STAR's own counts)
        │                    └─→ counts/star/gene_counts.{tsv,rds}
        │
        └─ salmon .... Salmon, decoy-aware, --gcBias --seqBias
                       └─ tximport ──→ counts/salmon/gene_counts.{tsv,rds}
                                       gene_counts_length_scaled.tsv, gene_tpm.tsv,
                                       transcript_counts.{tsv,rds}, ...

reference (once): genome (+ spike-ins) · GTF filtered to the genome · gffread transcriptome
                  · STAR index · Salmon decoy-aware index · collapsed GTF for RNA-SeQC
```

## Documentation

| Document | Contents |
|---|---|
| [docs/installation.md](docs/installation.md) | Requirements, the launcher environment, per-rule environments, first run |
| [docs/resources.md](docs/resources.md) | Getting the genome, annotation and spike-ins (`resources/get_reference.sh`); what gets built; disk and memory; checklist |
| [docs/usage.md](docs/usage.md) | Targets, samplesheets, local vs cluster runs, adding samples, troubleshooting |
| [docs/configuration.md](docs/configuration.md) | Key-by-key config reference |
| [docs/outputs.md](docs/outputs.md) | Every output file, which matrix to use for what, what to check in QC |
| [docs/methods.md](docs/methods.md) | Tools, versions, parameters, rationale and citations |
| [../../docs/usage.md](../../docs/usage.md) | The launcher flags shared by every sm-core pipeline |

### Tutorials

| Tutorial | Organism | Data | Exercises |
|---|---|---|---|
| [seqc](tests/tutorial/seqc/README.md) | human | SEQC/MAQC-III: 3× UHRR + 3× HBRR, 2×100 bp, two flowcells each, ERCC spike-ins | both branches on a full human reference; validation against ERCC ratios and TaqMan qPCR |
| [yeast](tests/tutorial/yeast/README.md) | *S. cerevisiae* | nf-core/rnaseq's test data | every rule in minutes; lane merging; GTF records off the genome; CI |

[examples/](examples/README.md) holds excerpts of the human tutorial's real outputs:
count matrices, QC tables and the MultiQC reports.

## Requirements

- Linux, conda (or mamba), and Snakemake via the suite's `environment.yml`.
- Human: about 40 GB of RAM for STAR (index build and alignment) and 30 GB for Salmon, and
  about 56 GB of disk for the reference and both indices. See
  [docs/resources.md](docs/resources.md).
- PBS or SLURM for `-cl` (see [../../docs/usage.md](../../docs/usage.md#cluster-runs-pbs-or-slurm)).

## Citation

Cite the pipeline via the repository's [CITATION.cff](../../CITATION.cff), together with the
tools listed in [docs/methods.md](docs/methods.md#citations).

## License

MIT. See [LICENSE](../../LICENSE).
