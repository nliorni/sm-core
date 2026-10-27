# Outputs

Everything is written under `output_dir`. The two branches write side by side and never
into each other's directories.

```
<output_dir>/
├── results/
│   ├── counts/
│   │   ├── star/                       STAR branch count matrix  (-w star)
│   │   │   ├── gene_counts.tsv
│   │   │   ├── gene_counts_summary.tsv
│   │   │   ├── gene_counts.rds
│   │   │   └── samples.tsv
│   │   └── salmon/                     Salmon branch count matrices  (-w salmon)
│   │       ├── gene_counts.tsv
│   │       ├── gene_counts_length_scaled.tsv
│   │       ├── gene_tpm.tsv
│   │       ├── gene_lengths.tsv
│   │       ├── transcript_counts.tsv
│   │       ├── transcript_tpm.tsv
│   │       ├── gene_counts.rds
│   │       ├── transcript_counts.rds
│   │       └── samples.tsv
│   ├── de/<matrix>/                    differential expression  (-w de; see below)
│   ├── qc/
│   │   ├── star/multiqc_report.html    one report per branch
│   │   ├── salmon/multiqc_report.html
│   │   ├── <branch>/strandedness.tsv   declared vs. inferred, per sample
│   │   ├── <branch>/biotypes_mqc.tsv   count composition by gene class
│   │   └── strandedness.tsv            STAR vs. Salmon calls (-w all)
│   ├── fastp/<sample>.fastp.{html,json}
│   ├── star/<sample>/                  <sample>.bam(.bai), Log.final.out, SJ.out.tab, ReadsPerGene.out.tab
│   ├── salmon/<sample>/                quant.sf, lib_format_counts.json, aux_info/, libParams/
│   └── rnaseqc/<sample>/               RNA-SeQC metrics and per-gene tables
├── reference/                          built once (unless reference.dir points elsewhere)
├── logs/<tool>/                        one log per rule and sample
└── benchmarks/<tool>/                  runtime and peak memory per job
```

Trimmed reads and merged lanes are `temp()`: deleted once both branches have read them.

## Which matrix do I use?

| Downstream use | Matrix |
|---|---|
| Differential expression, done by the pipeline | `-w de` with `de.analyses` set: PyDESeq2, ORA and GSEA for each contrast (`results/de/`) |
| Differential expression with DESeq2 | `counts/salmon/gene_counts.rds` → `DESeqDataSet(se, ~ condition)` (uses the length offsets) |
| Differential expression with edgeR / limma-voom from a plain table | `counts/salmon/gene_counts_length_scaled.tsv` |
| Integer read counts, traceable to the BAM; tools that require them | `counts/star/gene_counts.tsv` (or `.rds`) |
| Expression levels for plots, clustering, between-gene comparison | `counts/salmon/gene_tpm.tsv` (TPM: within-sample normalised; do not use for DE) |
| Isoform-level analysis (DTU/DTE) | `counts/salmon/transcript_counts.rds` |

## STAR branch — `-w star`

| File | Contents |
|---|---|
| `counts/star/gene_counts.tsv` | Integer fragment counts: `gene_id`, `gene_name`, one column per sample. Uniquely mapped fragments overlapping one gene only (htseq-count union model), from the strand-appropriate STAR column. |
| `counts/star/gene_counts_summary.tsv` | Per sample: the strandedness used, fragments assigned to genes, and those not counted — no feature, ambiguous (several genes), multimapping, unmapped. A large `N_noFeature` points at intronic/intergenic reads (gDNA, pre-mRNA) or a wrong strandedness. |
| `counts/star/gene_counts.rds` | The same as a `SummarizedExperiment`: assay `counts` (integer); `rowData` gene_name, gene_type, chrom; `colData` `samples.tsv`. |
| `counts/star/samples.tsv` | Sample metadata: every samplesheet column beyond the FASTQ paths, the number of FASTQ pairs merged, and the strandedness declared, inferred by STAR, and used. |
| `star/<sample>/<sample>.bam` (`.bai`) | Coordinate-sorted alignment of every read, unmapped included; read group = sample. |
| `star/<sample>/<sample>.Log.final.out` | STAR's mapping statistics (unique, multi, too short, chimeric…). |
| `star/<sample>/<sample>.SJ.out.tab` | Splice junctions: coordinates, motif, annotated or novel, unique/multi read support. |
| `star/<sample>/<sample>.ReadsPerGene.out.tab` | STAR's raw per-gene counts, all three strandedness columns. |
| `rnaseqc/<sample>/<sample>.metrics.tsv` | RNA-SeQC 2 metrics; see below. |
| `qc/star/multiqc_report.html` | fastp, STAR, RNA-SeQC, strandedness, count composition. |

## Salmon branch — `-w salmon`

| File | Contents |
|---|---|
| `counts/salmon/gene_counts.tsv` | Estimated counts per gene (sum of its transcripts' estimates); not integers. Pair with `gene_lengths.tsv` as offsets for DE. |
| `counts/salmon/gene_lengths.tsv` | Per sample, each gene's expression-weighted average transcript length: the offset that corrects `gene_counts.tsv` for isoform switching. |
| `counts/salmon/gene_counts_length_scaled.tsv` | Counts already scaled for transcript length and library size (`lengthScaledTPM`) — usable directly where offsets are not supported. |
| `counts/salmon/gene_tpm.tsv` | Transcripts per million, per gene. |
| `counts/salmon/transcript_{counts,tpm}.tsv` | Salmon's estimates per transcript. |
| `counts/salmon/gene_counts.rds` | `SummarizedExperiment` with assays `counts`, `abundance` (TPM), `avgTxLength`, `counts_length_scaled`; `rowData` and `colData` as in the STAR branch. `DESeq2::DESeqDataSet(se, design)` picks up `avgTxLength` as offsets automatically. |
| `counts/salmon/transcript_counts.rds` | Transcript-level `SummarizedExperiment` (`counts`, `abundance`, `length`; `rowData` gene_id, gene_name). |
| `counts/salmon/samples.tsv` | As for STAR, with Salmon's strandedness call. |
| `salmon/<sample>/quant.sf` | Salmon's per-transcript output. |
| `salmon/<sample>/lib_format_counts.json` | Fragments per observed library orientation: the strandedness evidence. |
| `salmon/<sample>/aux_info/meta_info.json` | Mapping rate, fragments processed and mapped, bias-model settings. |
| `qc/salmon/multiqc_report.html` | fastp, Salmon, strandedness, count composition. |

## Differential expression — `-w de`

One folder per count matrix, analysis, split (`all` without `split_by`) and contrast:

```
results/de/<matrix>/
├── summary.tsv, summary.{png,pdf}         DEGs up/down and GSEA sets per contrast, all analyses
└── <analysis>/<split>/
    ├── pca.{png,pdf}                      PCA of the model's samples (VST, 500 most variable genes)
    ├── sample_distances.{png,pdf}         Euclidean distances between samples, clustered
    ├── vst.tsv                            variance-stabilised expression (blind), genes x samples
    ├── normalized_counts.tsv              counts / size factors
    ├── size_factors.tsv
    └── <test>_vs_<reference>/
        ├── deseq2_results.tsv             every tested gene (see columns below)
        ├── significant_up.tsv             padj < alpha and log2FC >= min_abs_log2fc
        ├── significant_down.tsv           padj < alpha and log2FC <= -min_abs_log2fc
        ├── volcano.{png,pdf}              log2FC vs -log10 padj; strongest genes named
        ├── ma.{png,pdf}                   mean expression vs shrunken log2FC
        ├── heatmap.{png,pdf}              top genes, z-scored VST, all samples of the model
        └── enrichment/
            ├── ora_up.tsv, ora_down.tsv   over-representation of each list against the tested genes
            ├── ora.{png,pdf}
            ├── gsea.tsv                   preranked GSEA on all tested genes
            ├── gsea.{png,pdf}             strongest sets each way, NES with FDR
            └── gsea_curves.{png,pdf}      running enrichment score, top 3 sets each way
```

`deseq2_results.tsv` columns: `gene_id`, `gene_name`, `gene_type`, `chrom`, `baseMean`,
`log2FoldChange` (test / reference), `lfcSE`, `stat` (Wald), `pvalue`, `padj`
(Benjamini–Hochberg, NA where independent filtering removed the gene),
`log2FoldChange_shrunken`, `mean_normalized_<test>`, `mean_normalized_<reference>`, and
`regulation` (`up` / `down` / `ns`). Genes that failed the `min_count` pre-filter are not
listed. ERCC spike-ins are never tested.

`ora_*.tsv`: `collection`, `term`, `overlap` (significant genes in the set), `set_size`
(tested genes in the set), `pvalue`, `padj` (within the collection), `odds_ratio`,
`genes`. `gsea.tsv`: `collection`, `term`, `size`, `es`, `nes`, `pvalue`, `fdr`, `fwer`,
`leading_edge`. A positive NES means the set sits towards the genes higher in the test
level.

**Reading them.** Look at `pca` first: replicates should group, and a sample in the
wrong place points to a swap or an outlier (`exclude` it, and say so). Then the volcano
and the summary. In ORA and GSEA, pale bars did not pass the cut-off. A GSEA q-value
of `q < 0.001` (with 1,000 permutations) means no permutation reached the observed
score: it is the resolution of the test, not an exact value.

## QC: what to look at first

1. **Strandedness** (`qc/<branch>/strandedness.tsv`, top of each report). Any `MISMATCH` means
   a declared strandedness is wrong, and the STAR counts of that sample were taken from the
   wrong column. `undetermined` means degraded RNA, genomic DNA or a mixed library.
2. **Mapping.** STAR "Uniquely mapped reads %" below ~70% (human), or Salmon mapping rate
   below ~60%, deserves an explanation: rRNA, contamination, another species, adapter
   dimers.
3. **Count composition** (`biotypes_mqc.tsv`). rRNA above a few percent means failed
   rRNA depletion / poly-A selection; mitochondrial above ~20–30% (tissue-dependent) means
   stressed or dying cells.
4. **RNA-SeQC.** *Intergenic rate* above ~5–10% suggests genomic-DNA contamination (check
   for DNase treatment); a high *intronic rate* is expected for total RNA, not for poly-A;
   *3′ bias* well above 0.5–0.6 indicates degraded RNA (low RIN).
5. **Duplication** (fastp) is expected to be high in RNA-seq — highly expressed genes are
   sampled many times — and is not by itself a problem without UMIs.

## Logs and benchmarks

`logs/<tool>/<sample>.log` holds each job's tool output; on failure, start there.
`benchmarks/<tool>/<sample>.tsv` records wall time, CPU time and peak memory (`max_rss`) of
each heavy job — use them to set the `*_mem_mb` and `threads` keys for your data.
