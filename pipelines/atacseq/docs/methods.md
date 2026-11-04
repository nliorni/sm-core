# Methods

What sm-core/atacseq does, with which tool and why. The design follows the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) for alignment,
filtering and QC, and nf-core/atacseq for peak calling on fragments and the consensus peak
set.

## Tools

| Tool | Version | Used for |
|---|---|---|
| fastp | 1.3.7 | adapter and quality trimming, read QC |
| bowtie2 | 2.5.5 | alignment |
| samtools | 1.23.1 | filtering, sorting, statistics |
| Picard | 3.5.0 | duplicate marking |
| MACS3 | 3.1.0 | peak calling |
| subread (featureCounts) | 2.1.1 | peak counts, FRiP |
| deepTools | 4.0.0 | bigWigs |
| somalier | 0.3.6 | sample identity |
| pysam / NumPy | 0.24.1 / 2.3 | ATAC QC metrics |
| MultiQC | 1.35 | report |
| PyDESeq2 | 0.5.4 | differential accessibility |

Every tool runs in its own conda environment, pinned in `workflow/envs/`.

## Reference

The genome is decompressed and indexed with `samtools faidx`, and bowtie2-build builds the
index. From the GTF, `gtf_annotation.py` writes every gene, with its strand and TSS, and the
TSS set for the enrichment score. A gene's TSS is that of its `Ensembl_canonical`
transcript when the GTF tags one, otherwise the gene's 5′ end: the 5′ end is often a minor
transcript far upstream of the main promoter (for BLK, 7.5 kb). These are the TSSs of GENCODE's
`Ensembl_canonical` protein-coding transcripts, one per gene, when the GTF tags them;
otherwise those of all protein-coding transcripts. The blacklist is restricted to the
genome's sequences. The mitochondrial sequence is found by name, and the effective genome
size for MACS3 is the genome's non-N length unless configured.

## Reads and alignment

Lanes are concatenated. fastp trims adapters (detected from mate overlap, which needs no
adapter sequence) and low-quality bases, and drops pairs with a mate shorter than 20 bp.
bowtie2 aligns with `--very-sensitive -X 2000`, as ENCODE and the Omni-ATAC protocol do.
Fragments up to 2 kb are proper pairs, which keeps di- and tri-nucleosome fragments.

## Filtering

As in ENCODE:

1. Reads that are unmapped, have an unmapped mate, are secondary, supplementary or
   QC-failed, are not properly paired, or have MAPQ < 30 are removed.
2. `samtools fixmate -r` and a second filter remove the mates orphaned by step 1.
3. Picard MarkDuplicates marks duplicates on the coordinate-sorted pairs.
4. The final BAM drops the duplicates (`align.remove_duplicates`) and every read on the
   mitochondrion. A proper pair has both mates on one sequence, so dropping a whole
   sequence never splits a pair.

The blacklist is applied to peaks and to bigWigs rather than to reads, as ENCODE does: a
blacklisted region's reads then never reach a peak, a count or a track.

There is no separate Tn5-shifted BAM. Peaks and counts use whole fragments, which a shift
of 4–5 bp at each end does not change. The QC applies the shift where it matters: the cut
sites of the TSS profile. Footprinting tools apply their own.

## Peaks

MACS3 `callpeak -f BAMPE --keep-dup all -p 0.01`, per sample. With BAMPE, MACS3 piles up
each pair's actual fragment, so no shift or extension model is needed (nf-core/atacseq).
Peaks overlapping the blacklist are removed. The permissive per-sample cut-off is the
ENCODE ATAC-seq pipeline's, which then keeps peaks reproduced across replicates (IDR).
Here the replicate consensus below plays that part (`peaks.cutoff`, `peaks.threshold`).

### Consensus peak set

The peaks of all samples are merged where they overlap or touch. A merged region is kept
when at least `peaks.min_replicates` samples of one group called a peak in it, where the
groups come from `peaks.consensus_group_by` and the requirement is capped at each group's
size. nf-core/atacseq applies the replicate requirement over all samples instead. Applied
per group, a region open only in one condition (or one cell line) survives, which is the
point of a differential analysis, while a peak in a single replicate does not.

### Annotation

For each consensus peak, every gene of `annotation.gene_types` whose promoter window
[TSS − 1 kb, TSS + 1 kb] (in the gene's direction) the peak overlaps is a promoter gene of
that peak. A merged peak often covers several: bidirectional promoters, antisense and
readthrough genes. The peak's single nearest gene is its best promoter gene, or else the
gene with the nearest TSS. Ties at the same distance go to the earlier type in
`gene_types`, then to a named gene over one known only by its ID. Without this rule,
GENCODE's new unnamed lncRNA and readthrough genes take the promoters of TAL1 and CD79B.
The distance is 0 when the peak covers the TSS; otherwise it is the distance to the peak's
nearest edge, signed by the gene's direction. A peak is `promoter` when it has a promoter
gene, `genic` when it lies inside any gene body, and `distal` otherwise.

### Counts

featureCounts counts fragments (`-p --countReadPairs -B -C`: both mates mapped, on one
chromosome) overlapping each consensus peak, with `-O` so a fragment spanning two nearby
peaks counts in both.

## Signal tracks

deepTools bamCoverage draws each fragment whole (`--extendReads`) in 10 bp bins,
normalised to counts per million fragments, with the blacklist left out. Group tracks are
the mean of the group's samples (`bigwigAverage`).

## QC

### Library complexity

From the duplicate-marked BAM, over nuclear fragments, duplicates included. Each
fragment's position is its start and length: NRF = distinct positions / fragments, PBC1 =
positions seen once / distinct positions, PBC2 = positions seen once / positions seen
twice. These are ENCODE's definitions, computed on fragments as ENCODE does for
paired-end data.

### TSS enrichment

Each kept fragment contributes two Tn5 cut sites, its start + 4 and its end − 5 (the 9 bp
duplication of the Tn5 insertion). Cut sites are summed per base over ±2 kb around every
TSS, in the gene's direction. The profile is divided by the mean of its outer 100 bp on
each side (the background) and smoothed with a 51 bp moving average. The score is the
maximum within 250 bp of the TSS. This is the ENCODE definition in outline. ENCODE's
implementation differs in its TSS list and smoothing, so the scores are on a related but
not identical scale (see the ENCODE tutorial).

### FRiP

Fragments of the final BAM overlapping the sample's own peaks, over all its fragments
(featureCounts).

### Identity

With `qc.somalier_sites`, somalier extracts allele counts at about 17,000 common,
well-behaved SNPs from every final BAM (and from `qc.identity_bams`), then computes the
relatedness of every pair and each sample's sex from X heterozygosity and Y depth
(Pedersen et al. 2020). ATAC reads cover the sites unevenly, but enough of them fall in
open chromatin to separate individuals.

### ENCODE standards

`qc_summary.py` scores each library against the
[ENCODE ATAC-seq data standards](https://www.encodeproject.org/atac-seq/): the ideal and
acceptable tiers of [outputs.md](outputs.md#qcatac_qc_summarytsv).

## Differential accessibility

The consensus peak counts are tested with PyDESeq2, the Python implementation of DESeq2
(Muzellec et al. 2023), as in rnaseq's differential expression: the code is shared
(`smcore/de.py`). For each analysis and split (one model per value of `split_by`, or all
samples):

1. Peaks are kept when at least as many samples as the smallest group of the contrast's
   factor have ≥ `min_count` fragments (the DESeq2 vignette's pre-filter).
2. Size factors are computed (below), dispersions are fitted and shrunk to their trend, and
   the Wald test compares the two levels of each contrast, with Cook's-distance outlier
   refitting and independent filtering at `alpha`. Benjamini–Hochberg adjusts the p-values
   over the peaks of the model.
3. Shrunken log2 fold changes (PyDESeq2's apeGLM-style prior) are added for the MA plot.
   The significance calls use the unshrunken estimate and the Wald p-value.

A peak is `up` or `down` when padj < `alpha` and |log2FC| ≥ `min_abs_log2fc`. The PCA (500
most variable peaks), sample distances and heatmaps use the blind variance-stabilising
transformation.

### Normalisation

How to normalise ATAC-seq counts is not settled (Reske et al. 2020), because two
assumptions compete:

- **median of ratios** (DESeq2) assumes most peaks do not change. It corrects for
  differences in signal to noise between libraries (FRiP), which are common and technical,
  but it also removes a genuine global opening or closing of chromatin.
- **library size** (DiffBind's default since version 3) scales by the fragments in each
  final BAM. A global change survives, but so does a FRiP difference, as if it were
  biology.

`da.normalization` chooses one. Both are always computed: `size_factors.tsv` and its
figure set them side by side with each library's fraction of fragments in peaks, and the
MA plot draws where the other method puts zero, which is a constant shift: with size
factors s (used) and t (other), a peak unchanged under t has

  log2FC_used = mean_test(log2 t/s) − mean_reference(log2 t/s).

PyDESeq2 has no option for fixed size factors, so for `library_size` sm-core sets them on
the dataset and keeps PyDESeq2 from recomputing them. Given PyDESeq2's own factors, this
reproduces its results exactly (tested).

### Per-gene summary

A gene's peaks are those it is the nearest gene of, plus every peak covering its promoter.
`genes.tsv` gives, per gene, the counts of tested, opened and closed peaks, its most
significant peak and its most significant promoter peak: the gene-level view to compare
with RNA-seq.

## Citations

- **ENCODE ATAC-seq pipeline and standards**: ENCODE Project Consortium (2020), *Nature*
  583:699–710; https://www.encodeproject.org/atac-seq/.
- **ATAC-seq**: Buenrostro JD et al. (2013) *Nat Methods* 10:1213–1218. **Omni-ATAC**:
  Corces MR et al. (2017) *Nat Methods* 14:959–962.
- **nf-core/atacseq**: Patel H et al., doi:10.5281/zenodo.2634132; Ewels PA et al. (2020)
  *Nat Biotechnol* 38:276–278.
- **fastp**: Chen S et al. (2018) *Bioinformatics* 34:i884–i890.
- **bowtie2**: Langmead B & Salzberg SL (2012) *Nat Methods* 9:357–359.
- **samtools**: Danecek P et al. (2021) *GigaScience* 10:giab008.
- **Picard**: Broad Institute, https://broadinstitute.github.io/picard/.
- **MACS**: Zhang Y et al. (2008) *Genome Biol* 9:R137.
- **featureCounts**: Liao Y et al. (2014) *Bioinformatics* 30:923–930.
- **deepTools**: Ramírez F et al. (2016) *Nucleic Acids Res* 44:W160–W165.
- **ENCODE blacklist**: Amemiya HM, Kundaje A & Boyle AP (2019) *Sci Rep* 9:9354.
- **somalier**: Pedersen BS et al. (2020) *Genome Med* 12:62.
- **MultiQC**: Ewels P et al. (2016) *Bioinformatics* 32:3047–3048.
- **DESeq2**: Love MI, Huber W & Anders S (2014) *Genome Biol* 15:550. **PyDESeq2**:
  Muzellec B et al. (2023) *Bioinformatics* 39:btad547.
- **ATAC-seq normalisation**: Reske JJ, Wilson MR & Chandler RL (2020) *Epigenetics &
  Chromatin* 13:22. **DiffBind**: Ross-Innes CS et al. (2012) *Nature* 481:389–393.
- **Snakemake**: Mölder F et al. (2021) *F1000Research* 10:33.
