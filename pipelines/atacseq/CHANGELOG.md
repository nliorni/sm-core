# Changelog

All notable changes to sm-core/atacseq. Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning: [SemVer](https://semver.org/).

## [0.1.0] - Unreleased

First release.

### Added
- Paired-end ATAC-seq from FASTQ: fastp, bowtie2 (`--very-sensitive -X 2000`), ENCODE
  filtering (proper pairs, primary, MAPQ ≥ 30, orphan removal), Picard duplicate marking,
  duplicate and mitochondrial removal.
- MACS3 peaks per sample on fragments (BAMPE), with blacklisted peaks removed; FRiP.
  ENCODE's cut-off, p < 0.01, by default (`peaks.cutoff`, `peaks.threshold`), with a warning
  when a consensus group has no replicate.
- Consensus peak set, keeping regions called in enough replicates of one group
  (`peaks.consensus_group_by`, `peaks.min_replicates`); annotation with every gene whose
  promoter a peak covers, one nearest gene (canonical transcript TSSs, ties to named
  protein-coding genes), TSS distance and promoter / genic / distal class; fragment count matrix and sample table.
- CPM bigWigs per sample and averaged per group (`bigwig.group_by`).
- ATAC QC: library complexity (NRF, PBC1, PBC2), mitochondrial fraction,
  fragment-size distribution, TSS enrichment, chrX / chrY fragments. Every library is
  scored against the ENCODE standards in a table, three figures and MultiQC.
- Optional sample identity check with somalier, which can include BAMs of other assays
  (`qc.identity_bams`); `identity_matches.tsv` names each library's best genetic match in
  each assay, without relying on sample names.
- `resources/get_atac_resources.sh`: ENCODE blacklist v2 and somalier sites.
- **`-w da`**: differential accessibility of the consensus peaks with PyDESeq2, for every
  contrast in `da.analyses` (the analyses format of rnaseq's `de`: designs, `split_by`,
  `exclude`; checked at start-up). Median-of-ratios or library-size normalisation
  (`da.normalization`), with both always computed and compared. Per contrast: results
  joined with the peak annotation, significant peaks as TSV and BED, a per-gene summary,
  volcano, MA (with the other normalisation's zero), heatmap and peak-class figures. Per
  model: size factors, VST, normalised counts, PCA, sample distances. A summary table and
  figure. Shared with rnaseq: `smcore/de.py`, `smcore/rules/de_design.smk`.
- Tutorial `encode`: GM12878 and K562 at full depth on GRCh38, validated against ENCODE's
  own processing of the same reads.
