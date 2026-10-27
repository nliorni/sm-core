# Methods

What each step does, with which tool and version, and why that tool. Versions are pinned
exactly in `workflow/envs/*.yaml`; Snakemake builds each environment on first use.

## Tools

| Step | Tool | Version | Environment |
|---|---|---|---|
| Read trimming and read QC | fastp | 1.3.7 | `fastp` |
| Genome alignment, gene counts | STAR | 2.7.11b (build h5ca1c30_4, see `workflow/envs/star.yaml`) | `star` |
| BAM sorting and indexing | samtools | 1.23.1 | `star`, `samtools` |
| Transcript sequences from GTF | gffread | 0.12.7 | `gffread` |
| Transcript quantification | Salmon | 1.10.3 | `salmon` |
| Gene-level aggregation | tximport | 1.38.2 | `bioconductor` |
| Count containers (`.rds`) | SummarizedExperiment | 1.40.0 | `bioconductor` |
| BAM QC | RNA-SeQC | 2.4.2 | `rnaseqc` |
| QC report | MultiQC | 1.35 | `multiqc` |
| Differential expression | PyDESeq2 | 0.5.4 | `de` |
| ORA and GSEA | GSEApy | 1.3.1 | `de` |
| Gene sets | MSigDB (Hallmark, Reactome, GO:BP) | 2026.1 | downloaded once |
| Plots | matplotlib, adjustText | 3.11.1, 1.4.0 | `de` |
| Workflow | Snakemake | 9.16.0 | launcher (`environment.yml`) |

## Pipeline

**Samples.** Lanes or runs of the same library are concatenated (gzip members, no
recompression) before anything else, so every downstream step sees one FASTQ pair per
sample.

**Trimming — fastp.** Adapters are detected from the overlap of the two mates
(`--detect_adapter_for_pe`), which needs no adapter sequence and works for any kit; poly-G
tails from two-colour instruments are trimmed automatically; pairs with a mate shorter than
`trimming.min_length` (20) are dropped. fastp's before/after statistics are the read QC in
both MultiQC reports. With `trimming.enabled: false`, fastp still runs for QC but writes no
reads, and both branches read the raw FASTQ.

**Reference.** The annotation is restricted to sequences present in the genome and
transcripts without exons are dropped; spike-in sequences, if any, are appended to the
genome with one single-exon gene each. Transcript sequences are cut from that genome by
that annotation with gffread, so Salmon's transcriptome and STAR's gene models come from
one annotation with identical IDs. A transcript → gene table (gene name, biotype,
chromosome) built from the same GTF annotates every count matrix.

### STAR branch

**Alignment — STAR** two-pass mode (`--twopassMode Basic`), with the ENCODE long-RNA
parameters also used by nf-core/rnaseq: up to 20 loci per read
(`--outFilterMultimapNmax 20`), mismatches at most 4% of the mapped length, introns of
20 bp – 1 Mb, 8 bp overhang for novel junctions and 1 bp for annotated ones,
`--outFilterType BySJout`. `--runRNGseed 0` makes the placement of multimappers
reproducible. Unmapped reads are kept in the BAM (`--outSAMunmapped Within`). The BAM is
coordinate-sorted by samtools.

**Counts — STAR `--quantMode GeneCounts`.** Fragments are counted per gene during
alignment, with the same rules as `htseq-count` in its default union mode: a fragment is
counted for a gene when it overlaps that gene's exons and no other gene's; fragments mapping
to several loci, to no gene, or to several genes are not counted (they are reported in
`gene_counts_summary.tsv`). STAR writes the counts for unstranded, forward and reverse
libraries at once; the pipeline keeps the column matching the sample's strandedness.

**QC — RNA-SeQC 2** on the BAM, against a collapsed gene model (each gene's isoforms
merged into one, regions shared by two genes on the same strand removed — what GTEx does
for RNA-SeQC): exonic, intronic, intergenic and rRNA rates; 3′/5′ coverage bias; genes
detected; insert size.

### Salmon branch

**Quantification — Salmon** selective alignment against a decoy-aware index: the whole
genome is indexed as decoy sequence alongside the transcripts (k = 31), so a fragment that
matches an unannotated or intronic locus better than any transcript is discarded instead
of being forced onto the most similar transcript. `--gcBias` and `--seqBias` model
fragment-GC and random-hexamer priming biases, which otherwise distort the allocation of
reads among isoforms. The library type is `A` (detected by Salmon) unless a strandedness
is declared.

**Gene level — tximport.** Transcript estimates are summed per gene. Two count matrices
are written, for two kinds of downstream use: `gene_counts.tsv` (estimated counts,
`countsFromAbundance = "no"`) with `gene_lengths.tsv`, the per-sample average transcript
length of each gene, which DESeq2 and edgeR use as an offset to correct for changes in
isoform usage between samples; and `gene_counts_length_scaled.tsv` (`lengthScaledTPM`),
already corrected, for tools that take a plain matrix (limma-voom). `gene_counts.rds`
carries both, plus the offsets, as a SummarizedExperiment.

### Why two branches

They answer the same question with different evidence, and each has a use the other lacks.

- STAR's counts are integers, fragment by fragment traceable to the BAM, and do not depend
  on a statistical model; but they discard every multimapping fragment and every fragment
  overlapping two genes, which biases genes with paralogs or overlapping neighbours
  downwards.
- Salmon's estimates use those fragments, apportioned by likelihood, and correct for
  sequence biases and isoform length; benchmarks favour them for differential expression
  (Soneson et al. 2015; Love et al. 2018). But they are estimates, not observations.

Running both on the same reads makes disagreements visible: a gene that is
differentially expressed in one matrix and not the other is usually one of those
paralogous or overlapping cases, and worth a look before it is reported.

## Strandedness

Each branch infers strandedness from its own evidence: Salmon from the orientation of every
mapped fragment (`lib_format_counts.json`), STAR from which of its three count columns
collects the counts. ≥ 80% on one strand is stranded (forward or reverse); forward and
reverse within 10 points of each other is unstranded; anything else is undetermined — the
thresholds nf-core/rnaseq uses. A declared strandedness the data contradicts is flagged in
`results/qc/<branch>/strandedness.tsv` and the MultiQC report, and stops the run with
`qc.fail_on_strandedness_mismatch: true`. With `-w all`, `results/qc/strandedness.tsv`
compares the two branches' calls.

## Differential expression and enrichment

Run with `-w de` for the analyses in `de.analyses`
([configuration.md](configuration.md#de--differential-expression--w-de)), separately for
each count matrix in `de.matrices`.

**Input.** The STAR matrix is used as is. For Salmon, tximport's length-scaled counts
(`countsFromAbundance = "lengthScaledTPM"`), rounded to integers: the input tximport's
authors recommend for count-based tools that cannot take the per-sample gene-length
offsets, which PyDESeq2 cannot (Soneson et al. 2015). ERCC spike-ins are removed. Each
model uses only its own samples: one `split_by` value, minus `exclude`.

**Pre-filter.** A gene is kept with at least `min_count` (10) reads in at least as many
samples as the smallest group of the contrast's factor, as in the DESeq2 vignette.

**Model.** PyDESeq2, a Python reimplementation of DESeq2 (Muzellec et al. 2023; Love et
al. 2014): median-of-ratios size factors, gene-wise dispersions shrunk towards a fitted
trend, a negative binomial GLM with the configured design, Cook's distance outliers
refitted, and a Wald test of the contrast. Every design column is categorical. The
reference level is set to the contrast's reference, so the coefficient tested is test vs
reference. P-values are adjusted with Benjamini–Hochberg after independent filtering at
`alpha`. Shrunken fold changes come from PyDESeq2's `lfc_shrink` (an apeglm-style prior,
Zhu et al. 2019) and are used for the MA plot and as an optional GSEA ranking; the
significance calls use the unshrunken estimate. A gene is *up* or *down* when padj <
`alpha` and |log2 fold change| ≥ `min_abs_log2fc`.

PyDESeq2 and DESeq2 agree closely but not to the last digit. Report which one you used.

**Overview.** For each model: variance-stabilising transformation (blind to the design),
PCA of the 500 most variable genes (as DESeq2's `plotPCA`), and Euclidean distances
between samples with average-linkage clustering.

**ORA.** For the up and down lists separately, a hypergeometric test of each gene set
(GSEApy `enrich`). The background is the genes that were tested (padj not NA), not the
genome: genes that could never be called would otherwise make every set look enriched.
Sets are kept when `min_size`–`max_size` of their genes were tested. Adjusted p-values
(Benjamini–Hochberg) are within each collection.

**GSEA.** Preranked GSEA (Subramanian et al. 2005; GSEApy `prerank`) of all tested genes,
ranked by the Wald statistic (or the shrunken log2 fold change), with
`gsea.permutations` gene-set permutations and a fixed seed; FDR within each collection.
Where several gene IDs share a symbol, the one with the largest |statistic| is kept.

**Gene sets.** MSigDB (Liberzon et al. 2015) release `msigdb_version`, human (Hs) or
mouse (Mm) gene symbols: Hallmark, Reactome canonical pathways (Milacic et al. 2024)
and GO Biological Process (Gene Ontology Consortium 2023). Downloaded once into
`<reference.dir>/gene_sets/` and reused.

**Validation.** On the SEQC tutorial (brain vs UHRR, 3 vs 3), the log2 fold changes agree
with MAQC TaqMan qPCR at Pearson r = 0.96 in both matrices. Of the genes with a clear
qPCR change (≥ 2-fold), 89–91% are called significant, and 99% of those calls go the
qPCR's way. Genes up in brain are enriched for synaptic signalling and the neuronal
system; genes up in UHRR, a pool of cancer cell lines, for E2F, MYC and G2/M targets.
See [the tutorial](../tests/tutorial/seqc/README.md).

**Figures.** Every plot uses the shared sm-core style (`smcore/plotting`): one colour
meaning across the suite (red up, blue down, grey not significant), Arial-like type, PNG
for viewing and PDF with editable text. The style is adapted from Chen Liu's
[figures4papers](https://github.com/ChenLiu-1996/figures4papers).

## Citations

If you use this pipeline, cite it through the repository's `CITATION.cff`, together with:

- **fastp** — Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
  preprocessor. *Bioinformatics* 34:i884–i890 (2018). doi:10.1093/bioinformatics/bty560
- **STAR** — Dobin A et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*
  29:15–21 (2013). doi:10.1093/bioinformatics/bts635
- **Salmon** — Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast
  and bias-aware quantification of transcript expression. *Nature Methods* 14:417–419
  (2017). doi:10.1038/nmeth.4197
- **Salmon decoys / selective alignment** — Srivastava A et al. Alignment and mapping
  methodology influence transcript abundance estimation. *Genome Biology* 21:239 (2020).
  doi:10.1186/s13059-020-02151-8
- **tximport** — Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq:
  transcript-level estimates improve gene-level inferences. *F1000Research* 4:1521 (2015).
  doi:10.12688/f1000research.7563.2
- **Swimming downstream** — Love MI et al. Swimming downstream: statistical analysis of
  differential transcript usage following Salmon quantification. *F1000Research* 7:952
  (2018). doi:10.12688/f1000research.15398.3
- **RNA-SeQC 2** — Graubert A, Aguet F, Ravi A, Ardlie KG, Getz G. RNA-SeQC 2: efficient
  RNA-seq quality control and quantification for large cohorts. *Bioinformatics*
  37:3048–3050 (2021). doi:10.1093/bioinformatics/btab135
- **gffread** — Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare.
  *F1000Research* 9:304 (2020). doi:10.12688/f1000research.23297.2
- **SAMtools** — Danecek P et al. Twelve years of SAMtools and BCFtools. *GigaScience*
  10:giab008 (2021). doi:10.1093/gigascience/giab008
- **SummarizedExperiment** — Morgan M et al. SummarizedExperiment (Bioconductor).
  doi:10.18129/B9.bioc.SummarizedExperiment
- **MultiQC** — Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis
  results for multiple tools and samples in a single report. *Bioinformatics*
  32:3047–3048 (2016). doi:10.1093/bioinformatics/btw354
- **PyDESeq2** — Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python
  package for bulk RNA-seq differential expression analysis. *Bioinformatics*
  39:btad547 (2023). doi:10.1093/bioinformatics/btad547
- **DESeq2** — Love MI, Huber W, Anders S. Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550 (2014).
  doi:10.1186/s13059-014-0550-8
- **apeglm** — Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence
  count data: removing the noise and preserving large differences. *Bioinformatics*
  35:2084–2092 (2019). doi:10.1093/bioinformatics/bty895
- **GSEApy** — Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing
  gene set enrichment analysis in Python. *Bioinformatics* 39:btac757 (2023).
  doi:10.1093/bioinformatics/btac757
- **GSEA** — Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach
  for interpreting genome-wide expression profiles. *PNAS* 102:15545–15550 (2005).
  doi:10.1073/pnas.0506580102
- **MSigDB Hallmark** — Liberzon A et al. The Molecular Signatures Database hallmark gene
  set collection. *Cell Systems* 1:417–425 (2015). doi:10.1016/j.cels.2015.12.004
- **Reactome** — Milacic M et al. The Reactome Pathway Knowledgebase 2024. *Nucleic Acids
  Research* 52:D672–D678 (2024). doi:10.1093/nar/gkad1025
- **Gene Ontology** — Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023.
  *Genetics* 224:iyad031 (2023). doi:10.1093/genetics/iyad031
- **Benjamini–Hochberg** — Benjamini Y, Hochberg Y. Controlling the false discovery rate.
  *JRSS B* 57:289–300 (1995). doi:10.1111/j.2517-6161.1995.tb02031.x
- **Snakemake** — Mölder F et al. Sustainable data analysis with Snakemake.
  *F1000Research* 10:33 (2021). doi:10.12688/f1000research.29032.2

Tutorial data: SEQC/MAQC-III Consortium. A comprehensive assessment of RNA-seq accuracy,
reproducibility and information content by the Sequencing Quality Control Consortium.
*Nature Biotechnology* 32:903–914 (2014). doi:10.1038/nbt.2957 — and, for the TaqMan
truth set, MAQC Consortium. *Nature Biotechnology* 24:1151–1161 (2006).
doi:10.1038/nbt1239
