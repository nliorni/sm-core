# Human tutorial: SEQC/MAQC-III

A real human dataset with two published truth sets. Both count matrices are validated
against them, not only run.

## Dataset

The Sequencing Quality Control (SEQC/MAQC-III) consortium sequenced two commercial
reference RNAs at many sites, both spiked with the ERCC synthetic RNAs:

| Sample | RNA | Spike-in |
|---|---|---|
| **A** | Universal Human Reference RNA (UHRR, Agilent): pooled RNA of 10 cancer cell lines | ERCC Mix 1 |
| **B** | Human Brain Reference RNA (HBRR, Ambion): pooled adult brain | ERCC Mix 2 |

This tutorial takes the BGI site's Illumina HiSeq 2000 runs (SRA SRP025982, GEO
GSE47792): three independent libraries of A and three of B, 2×100 bp. Each library was
sequenced on two flowcells, so each sample is **two FASTQ pairs**, which the pipeline
merges, as for any multi-lane project. That gives about 11–15 million read pairs per
sample.

| Sample | Runs (flowcell AC0AYTACXX, BC0AYPACXX; lane 1) |
|---|---|
| SEQC_A_1 | SRR896663, SRR896664 |
| SEQC_A_2 | SRR896679, SRR896680 |
| SEQC_A_3 | SRR896695, SRR896696 |
| SEQC_B_1 | SRR896743, SRR896744 |
| SEQC_B_2 | SRR896759, SRR896760 |
| SEQC_B_3 | SRR896775, SRR896776 |

Reference: GENCODE release 50, GRCh38 primary assembly, plus the 92 ERCC sequences as
spike-ins (`reference.spikein_fasta`). Strandedness is left to `auto` on purpose, so each
branch has to work it out.

### Truth sets

1. **ERCC spike-ins.** 92 synthetic transcripts in four subgroups of 23, at Mix 1 / Mix 2
   ratios of **4, 1, 2/3 and 1/2**. A carries Mix 1 and B carries Mix 2, so each
   spike-in's A/B ratio is known exactly. UHRR and HBRR hold different amounts of mRNA per
   microgram of total RNA, which shifts every spike-in's ratio by the same constant. The
   validation therefore centres on the ratio-1 subgroup, as the ERCC dashboard does.
2. **MAQC TaqMan qPCR** (GEO GPL4097, MAQC Consortium 2006): about 1,000 genes measured in
   A and B by TaqMan, four replicates each. The 830 genes detected in all eight assays
   give a qPCR log2(B/A) to compare each matrix's log2(B/A) against. This is the reference
   that RNA-seq quantification benchmarks have used for a decade.

## Run

About 12 GB of reads and 1 GB of reference. The download is about 100× faster with
sra-tools installed (`conda install -c bioconda sra-tools`), which lets the script fetch
from SRA's AWS mirror; otherwise it falls back to ENA.

```bash
cd pipelines/rnaseq
tests/tutorial/seqc/download_tutorial_data.sh

./run.py -w all -c tests/tutorial/seqc/config_tutorial.yaml -q 32 -n     # check
./run.py -w all -c tests/tutorial/seqc/config_tutorial.yaml -q 32
# or one cluster job per rule:
./run.py -w all -c tests/tutorial/seqc/config_tutorial.yaml -cl --scheduler slurm -qu <partition> -j 20

python3 tests/tutorial/seqc/validate.py
```

The first run builds the human STAR and Salmon indices (see "Cost" below). After that,
they sit in `output/reference/` and are reused.

## Validation

`validate.py` normalises each branch's `gene_counts.tsv` (DESeq2 median-of-ratios), averages
the three A and the three B libraries, and computes log2(B/A) per gene:

Measured on this pipeline's own run (2026-10-21; GENCODE 50, `strandedness: auto`):

```
                                   STAR branch   Salmon branch   expected
ERCC subgroup A  (Mix1/Mix2 = 4)      +1.91          +1.90         +2.00
ERCC subgroup B  (1)                   0.00           0.00          0.00   (centring group)
ERCC subgroup C  (2/3)                -0.59          -0.59         -0.58
ERCC subgroup D  (1/2)                -0.99          -1.00         -1.00
ERCC, observed vs expected, r          0.992          0.992        (43 / 44 spike-ins detected)
TaqMan log2(B/A), Pearson r            0.964          0.968        (665 genes)
TaqMan log2(B/A), Spearman rho         0.959          0.967
STAR vs Salmon log2(B/A), r                 0.995               (14,216 genes; 72 differ > 2x)
```

How to read it:

- **Both matrices recover the spike-in ratios.** The medians of the four ERCC subgroups fall
  within 0.1 log2 units of their designed values in both branches. The mRNA-content shift
  between UHRR and HBRR (−0.6 log2) is the known effect of the two RNAs carrying different
  mRNA fractions. It's the reason for centring on subgroup B.
- **Both agree closely with qPCR** (r = 0.96–0.97 over 665 genes). Correlations around 0.9
  are commonly reported for RNA-seq against MAQC TaqMan. Salmon is ahead by a small,
  consistent margin, in line with Soneson et al. 2015: it uses the multimapping and
  gene-overlapping fragments that STAR's counts discard.
- **The two branches agree with each other** (r = 0.995). Look first at the 72 genes whose
  fold changes differ by more than 2× if a result depends on the choice of matrix. In 60
  of them, Salmon has more than twice STAR's counts (median 7×): fragments STAR discarded
  as multimapping or ambiguous. They're genes with processed pseudogenes (RPS23, NME2,
  UQCRHL, PFDN5) and genes overlapping an antisense or divergent lncRNA (NDUFA6-DT,
  ENTPD1-AS1). In an unstranded library like this one, STAR cannot tell such genes apart.
  This is the trade-off described in [methods.md](../../../docs/methods.md#why-two-branches).

Other figures from the same run:

| Sample | Read pairs | STAR uniquely mapped | Salmon mapping rate | Exonic / intronic / intergenic (RNA-SeQC) | Genes detected | Mitochondrial share |
|---|---|---|---|---|---|---|
| SEQC_A_1 | 11.2 M | 88.3% | 86.2% | 82.1 / 12.4 / 1.8% | 23,860 | 6.3% |
| SEQC_A_2 | 9.3 M | 90.1% | 87.5% | 83.6 / 11.3 / 1.5% | 23,109 | 6.7% |
| SEQC_A_3 | 11.2 M | 92.2% | 89.1% | 85.6 / 9.6 / 1.3% | 24,019 | 6.6% |
| SEQC_B_1 | 13.5 M | 92.5% | 87.6% | 85.1 / 10.5 / 1.5% | 24,540 | 12.2% |
| SEQC_B_2 | 13.2 M | 92.1% | 87.1% | 84.5 / 10.9 / 1.5% | 24,469 | 12.3% |
| SEQC_B_3 | 14.0 M | 93.4% | 88.2% | 85.9 / 9.8 / 1.4% | 24,735 | 12.2% |

- **Strandedness:** both branches call every library **unstranded**. The forward share is
  0.500 for STAR and 0.499 for Salmon. That is right: SEQC's Illumina sites used the
  non-stranded TruSeq RNA kit. Declaring these libraries `reverse`, a common assumption
  for modern data, would have been flagged as a `MISMATCH`.
- **Brain has twice UHRR's mitochondrial share** (12% vs 6.5%). That is biology, not a
  QC failure: neurons are mitochondria-rich.
- rRNA is below 0.05% everywhere (poly-A selected), and spike-ins are 1–1.6% of counts.

### Differential expression (`-w de`)

The tutorial config also defines one analysis, brain (B) against UHRR (A), run in both
matrices with `-w de` (about 7 minutes). `validate.py` then checks PyDESeq2's results
against TaqMan:

```
                                            STAR branch   Salmon branch
genes called up / down (padj<0.05, 2x)    5,506 / 4,729   5,692 / 4,886
PyDESeq2 log2FC vs TaqMan, Pearson r            0.962           0.960
  shrunken log2FC vs TaqMan                     0.961           0.956
TaqMan genes changed >= 2x: called significant   91%             89%
  ... of those, in the qPCR direction        364 of 366      360 of 363
```

The enrichment results are what the two RNAs predict. Genes up in brain over-represent
synaptic signalling, the neuronal system and axon development. Genes up in UHRR, a pool
of cancer cell lines, over-represent E2F targets, MYC targets, the G2/M checkpoint and
ribosome biogenesis: proliferation. GSEA agrees (Hallmark MYC targets V1 NES -3.4, E2F
targets -3.4). The volcano plot names textbook markers on each side: GFAP, SNAP25 and
CAMK2A for brain; ALB, APOB and AFP (UHRR includes a liver cancer line), immunoglobulin
light chains (a B-cell line) and PMEL (a melanoma line) for UHRR.

`validate.py` uses no random sampling, and GSEA uses a fixed seed, so a rerun on the
same reads with the same tool versions reproduces these figures exactly.

## Cost

On a 2011-era AMD Opteron 6272 node with `-q 32 --resources mem_mb=80000` (a current CPU
should take about half the time):

| Step | Wall time | Peak memory |
|---|---|---|
| STAR index, GRCh38 + GENCODE 50 + ERCC (16 threads), once | 1 h 54 min | 38 GB |
| Salmon decoy-aware index (16 threads), once | 1 h 11 min | 20 GB |
| fastp, per sample | ~4 min | 1.3 GB |
| STAR alignment, per sample (16 threads) | 15–19 min | 30 GB |
| Salmon quantification, per sample (16 threads) | 17–20 min | 27 GB |
| RNA-SeQC, per sample | 6–9 min | 0.4 GB |
| **Whole tutorial, indices included** | **3 h 30 min** | |

Disk: 13 GB of FASTQ, 56 GB of reference and indices, 11 GB of BAM. Everything else is small.
