# Changelog

All notable changes to sm-core/rnaseq are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the pipeline uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-10-20

First release: paired-end short-read RNA-seq, FASTQ to count matrices.

### Added
- Two independent quantification branches, each its own workflow target with its own count
  matrix and MultiQC report:
  - `star`: STAR 2.7.11b two-pass alignment (ENCODE parameters) to a sorted, indexed BAM;
    per-gene fragment counts from `--quantMode GeneCounts`, read from the column that
    matches the sample's strandedness; RNA-SeQC 2.4.2 on a GTEx-style collapsed annotation.
  - `salmon`: Salmon 1.10.3 selective alignment against a decoy-aware index, with GC and
    sequence-bias correction; tximport 1.38.2 gene and transcript matrices (counts,
    `lengthScaledTPM` counts, TPM, average transcript lengths).
  - Both write a `SummarizedExperiment` (`gene_counts.rds`) with samplesheet metadata as
    `colData` and gene name, biotype and chromosome as `rowData`.
- fastp 1.3.7 trimming (overlap-based adapter detection) and read QC.
- Strandedness inference per branch from its own evidence, checked against the declared
  value (nf-core thresholds), and compared between branches with `-w all`.
- Count-composition QC (protein-coding, lncRNA, rRNA, mitochondrial, spike-in shares).
- Reference preparation from a genome + GTF: the GTF is filtered to the genome's
  sequences, spike-ins are added as one gene each, the transcriptome is cut with
  gffread 0.12.7, and both indices are built, all shareable across projects
  (`reference.dir`).
- nf-core-compatible samplesheets (lanes merged, metadata columns carried through), or
  the sm-core FASTQ-prefix convention.
- Config validation and defaults from a JSON schema.
- Tutorials: human SEQC/MAQC-III (ERCC spike-ins and TaqMan qPCR as truth sets) and the
  nf-core yeast test data (with ERCC92 added, so CI covers the spike-in path); CI dry-runs
  every target and runs the yeast tutorial.
- Validated on SEQC: ERCC ratio groups recovered within 0.1 log2 (r = 0.992, both
  branches); TaqMan fold changes r = 0.964 (STAR) / 0.968 (Salmon), 665 genes; STAR vs
  Salmon r = 0.995. `examples/` holds excerpts and both MultiQC reports.

### Fixed during validation
- STAR pinned to bioconda build `2.7.11b=h5ca1c30_4`: later builds are AVX2-only and die
  with SIGILL on pre-AVX2 CPUs. Index generation calls `STAR-plain`, the only build linked
  with OpenMP (single-threaded otherwise: ~6 h instead of under 2 h for human).
- Salmon gentrome: spike-ins were placed after the decoys, which Salmon rejects.
- Default memory raised from what the human run measured: `salmon.mem_mb` 16 → 32 GB
  (peak 27 GB), `reference.star_index_mem_mb` 40 → 48 GB (peak 38 GB).
- MultiQC custom-content headers quoted for YAML; tximport skips inferential replicates.
