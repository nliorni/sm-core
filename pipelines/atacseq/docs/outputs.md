# Outputs

Everything is under `<output_dir>`:

```
results/
├── fastp/                      <sample>.fastp.{json,html}   read QC (trimmed reads are temporary)
├── bam/                        <sample>.bam(.bai)           final BAM: filtered, deduplicated, no chrM
├── peaks/
│   ├── <sample>/               <sample>_peaks.narrowPeak, _summits.bed, _peaks.xls   MACS3
│   └── consensus/
│       ├── consensus_peaks.bed             the consensus peak set (BED6, peak_id in genome order)
│       ├── consensus_peaks.saf             the same, for featureCounts
│       ├── consensus_peaks.support.tsv     how many samples of each group called each peak
│       └── consensus_peaks.annotated.tsv   nearest gene, distance to TSS, class
├── counts/
│   ├── peak_counts.tsv                     fragments per consensus peak per sample
│   ├── samples.tsv                         sample metadata (the samplesheet's extra columns)
│   └── featurecounts.txt.summary           assigned / unassigned fragments per sample
├── bigwig/                     <sample>.bw, groups/<group>.bw   CPM-normalised coverage
├── da/                         differential accessibility (-w da, when da.analyses is set)
└── qc/
    ├── atac_qc_summary.tsv                 every metric of every sample, with ENCODE tiers
    ├── qc_metrics.png/.pdf                 the same, as bars against the ENCODE thresholds
    ├── fragment_sizes.png/.pdf             fragment-length distributions (nucleosome ladder)
    ├── tss_enrichment.png/.pdf             TSS profiles
    ├── multiqc_report.html
    ├── samples/<sample>/                   per-sample QC files (Picard, samtools, FRiP, atac_qc)
    └── identity/                           somalier.{html,pairs.tsv,samples.tsv}, identity_matches.tsv   (qc.somalier_sites)
logs/        one log per job, by tool
benchmarks/  time and memory per job
```

## The count matrix

`counts/peak_counts.tsv` has one row per consensus peak (`peak_id, chrom, start, end`,
0-based half-open as in BED) and one integer column per sample. Each value is the number of
fragments (read pairs) in the final BAM that overlap the peak. A fragment overlapping two
neighbouring peaks counts in both (featureCounts `-O`). `counts/samples.tsv` holds the
sample metadata in the same order. Together they are the input of the differential
accessibility below (`-w da`), or of any other count-based analysis. Join
`consensus_peaks.annotated.tsv` on `peak_id` for the gene names.

### `consensus_peaks.annotated.tsv`

| Column | Meaning |
|---|---|
| `peak_id`, `chrom`, `start`, `end`, `width` | the peak |
| `n_samples` | samples with a peak in this region |
| `nearest_gene_id`, `nearest_gene_name`, `nearest_gene_type`, `nearest_gene_strand` | one gene per peak (of `annotation.gene_types`): its best promoter gene, or else the gene whose TSS is nearest |
| `distance_to_tss` | 0 if the peak covers the TSS; otherwise the distance from the TSS to the nearest edge of the peak, negative upstream and positive downstream (in the gene's direction) |
| `annotation` | `promoter` (overlaps the promoter window of at least one gene), `genic` (inside a gene body), `distal` |
| `promoter_gene_ids`, `promoter_gene_names` | every gene whose promoter window the peak overlaps, comma-separated, best first; empty for non-promoter peaks |

## Differential accessibility (`-w da`)

```
results/da/
├── summary.tsv, summary.png/.pdf         peaks up / down (and at promoters) per contrast
└── <analysis>/<split>/                   <split> is "all" without split_by
    ├── size_factors.tsv, size_factors.png/.pdf   both normalisations, library sizes, fraction in peaks
    ├── vst.tsv, normalized_counts.tsv            variance-stabilised and normalised counts
    ├── pca.png/.pdf, sample_distances.png/.pdf
    └── <test>_vs_<reference>/
        ├── da_results.tsv                       every tested peak, most significant first
        ├── significant_up.tsv / _down.tsv       the significant peaks
        ├── significant_up.bed / _down.bed       the same as BED6, for a browser or motif tools
        ├── genes.tsv                            per gene: its peaks, best peak, best promoter peak
        └── volcano, ma, heatmap, peak_classes   .png/.pdf
```

What to look at, in order:

1. **`pca` and `sample_distances`**: replicates together, groups apart. An outlier or a swap
   shows here first.
2. **`size_factors`**: the two normalisations side by side, with each library's fraction of
   fragments in consensus peaks. Close together: the choice does not matter. Apart, and
   tracking the fraction in peaks: the libraries differ in signal to noise, which
   median-of-ratios corrects and library size does not. Apart between groups with similar
   fractions in peaks: a global change of accessibility, which only library size keeps.
3. **`ma`**: the gold dashed line is where the other normalisation puts zero. Peaks between
   the two lines are called differently by the two methods.
4. **`volcano`, `heatmap`, `peak_classes`**: the size and direction of the response, and
   whether it lies at promoters or at distal (enhancer) sites.

### `da_results.tsv`

| Column | Meaning |
|---|---|
| `peak_id`, `chrom`, `start`, `end`, `width` | the consensus peak (0-based, half-open) |
| `nearest_gene_id`, `nearest_gene_name`, `nearest_gene_type`, `distance_to_tss`, `annotation`, `promoter_gene_names`, `promoter_gene_ids` | from `consensus_peaks.annotated.tsv` |
| `baseMean` | mean normalised count over the model's samples |
| `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` | PyDESeq2 Wald test, test / reference |
| `log2FoldChange_shrunken` | shrunken estimate (`da.shrink`) |
| `mean_normalized_<level>` | mean normalised count in each of the two levels |
| `regulation` | `up` / `down` (padj < `alpha` and \|log2FC\| ≥ `min_abs_log2fc`) or `ns` |

`padj` is empty for peaks removed by independent filtering or flagged as Cook's outliers.

### `genes.tsv`

One row per gene with a tested peak, its peaks being those it is the nearest gene of plus
those covering its promoter: `gene_name`, `gene_type`, `peaks_tested`,
`peaks_up`, `peaks_down`; its most significant peak (`best_peak`, with class, log2FC, padj,
regulation); and its most significant promoter peak (`promoter_peak`, ...), empty when it
has none. This is the table to set beside the RNA-seq results of the same contrast:
`promoter_log2FoldChange` against the gene's expression change.

## QC: what to look at first

1. **`qc/qc_metrics.png`**: every library against ENCODE's thresholds. Green is ideal,
   gold is acceptable, red is below acceptable.
2. **`qc/tss_enrichment.png`**: a sharp peak at the TSS, well above 1. A flat or low
   profile means little signal: dead cells, over- or under-digestion, or a failed library.
3. **`qc/fragment_sizes.png`**: a nucleosome-free peak below 147 bp, then mono- (~200 bp)
   and di-nucleosome (~400 bp) peaks. No ladder means degraded chromatin or a non-Tn5
   library.
4. **The mitochondrial fraction**: 20–60% is common in standard ATAC. Omni-ATAC
   typically brings it under 10%. A high fraction only costs depth.
5. **Identity** (`qc/identity/somalier.html`), when enabled: replicates of one line should
   have relatedness near 1, different lines near 0 (or 0.5 for related individuals). Check
   that the inferred sex matches the expected one. `qc/identity/identity_matches.tsv` names
   each library's best genetic match among the ATAC libraries and among the
   `qc.identity_bams` (usage.md, *Checking sample identity across assays*).

### `qc/atac_qc_summary.tsv`

| Column | Meaning | ENCODE ideal / acceptable |
|---|---|---|
| `aligned_rate` | bowtie2 overall alignment rate | 95% / 80% |
| `fragments_total` | properly paired MAPQ-filtered fragments, before deduplication | |
| `mito_fraction` | share of those on the mitochondrion | |
| `duplicate_fraction` | Picard's duplicate share of read pairs | |
| `fragments_final` | fragments in the final BAM | ≥ 25 M |
| `NRF` | distinct / total fragment positions (nuclear) | 0.9 / 0.7 |
| `PBC1` | positions seen once / distinct positions | 0.9 / 0.7 |
| `PBC2` | positions seen once / positions seen twice | 3 / 1 |
| `FRiP` | fraction of the sample's fragments in its own peaks | 0.3 / 0.2 |
| `peaks` | MACS3 peaks (after the blacklist) | |
| `tss_enrichment` | TSS enrichment score ([methods](methods.md#tss-enrichment)) | 7 / 5 |
| `nucleosome_free_fraction`, `mononucleosome_fraction` | fragments < 147 bp; 147–294 bp | |
| `median_fragment_length` | | |
| `chrX_per_million`, `chrY_per_million` | fragments on X / Y per million final fragments: a quick sex check | |
| `<metric>_status` | `ideal`, `acceptable` or `low` | |
| `below_acceptable` | the metrics below the acceptable tier | |

ENCODE's TSS thresholds were set with its own TSS list and scaling. sm-core's score uses
the GTF's canonical protein-coding TSSs and comes out on a somewhat different scale; the
ENCODE tutorial compares the two directly.

## Logs and benchmarks

Every job logs to `logs/<tool>/`, and the heavy ones record their time and peak memory in
`benchmarks/<tool>/`. The bowtie2 log (`logs/bowtie2/<sample>.log`) holds bowtie2's
alignment summary.
