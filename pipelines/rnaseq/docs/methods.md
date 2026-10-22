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
- **Snakemake** — Mölder F et al. Sustainable data analysis with Snakemake.
  *F1000Research* 10:33 (2021). doi:10.12688/f1000research.29032.2

Tutorial data: SEQC/MAQC-III Consortium. A comprehensive assessment of RNA-seq accuracy,
reproducibility and information content by the Sequencing Quality Control Consortium.
*Nature Biotechnology* 32:903–914 (2014). doi:10.1038/nbt.2957 — and, for the TaqMan
truth set, MAQC Consortium. *Nature Biotechnology* 24:1151–1161 (2006).
doi:10.1038/nbt1239
