# Tutorial: ENCODE GM12878 and K562

A complete sm-core/atacseq project on real data, done the way you would do your own: two
human cell lines, two biological replicates each, at full depth, against GRCh38, from the
raw reads to differential accessibility. ENCODE processed the same reads with its own
ATAC-seq pipeline and published the QC metrics and peak sets, so at the end every result
is checked against an independent processing of identical input.

| Library | ENCODE experiment | Replicate | Read pairs |
|---|---|---|---|
| `GM12878_1` | [ENCSR095QNB](https://www.encodeproject.org/experiments/ENCSR095QNB/) | 1 | 38.2 M |
| `GM12878_2` | ENCSR095QNB | 2 | 34.7 M |
| `K562_1` | [ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/) | 1 | 42.0 M |
| `K562_2` | ENCSR483RKN | 2 | 39.4 M |

GM12878 is a B-lymphoblastoid line, K562 a chronic myeloid leukaemia line with erythroid
features. Both experiments come from the same lab (Snyder, Stanford), protocol and
sequencer (2×100 bp), and pass ENCODE's audits.

**What it takes**: 18 GB of reads, 7 GB for the reference (the bowtie2 index included) and
11 GB of results; 16 GB of RAM (Picard MarkDuplicates peaks at 12.6 GB, bowtie2 at
8.4 GB); 11 hours on 48 cores of an older cluster node, of which aligning took 4–8 hours
per library (four in parallel, 16 threads each) and the index 1 hour, built once. On a
current node, expect about half.

## 1. Data

```bash
cd pipelines/atacseq
tests/tutorial/encode/download_tutorial_data.sh
# with a GENCODE GRCh38 genome + GTF already on disk (e.g. the rnaseq reference):
tests/tutorial/encode/download_tutorial_data.sh \
    --genome /refs/GRCh38.primary_assembly.genome.fa.gz \
    --gtf /refs/gencode.v50.primary_assembly.annotation.gtf.gz
```

This writes:

- `fastq/`: the eight FASTQ files, checked against ENCODE's MD5 sums (a rerun resumes and
  skips what is complete);
- `reference/`: links to the genome and GTF, the ENCODE blacklist v2 and the somalier
  sites (`resources/get_atac_resources.sh`, as for a real project);
- `truth/`: what `validate.py` compares with: ENCODE's conservative IDR and
  pseudoreplicated peaks for each line, and `encode_qc.tsv`, the per-replicate QC from
  ENCODE's pipeline.

## 2. Samplesheet and config

[`samplesheet.csv`](samplesheet.csv) is an nf-core-style sheet with two metadata columns:

```csv
sample,fastq_1,fastq_2,cell_line,replicate
GM12878_1,fastq/GM12878_1_R1.fastq.gz,fastq/GM12878_1_R2.fastq.gz,GM12878,1
...
```

Any column after `fastq_2` is metadata: it can define groups and designs. A real project
would add `condition` (and keep `cell_line`, `replicate`).

[`config_tutorial.yaml`](config_tutorial.yaml) sets only what differs from the defaults of
[`config/config.yaml`](../../../config/config.yaml), the template a real project copies:

| Setting | Here | Why |
|---|---|---|
| `reference.genome`, `.gtf`, `.blacklist` | GENCODE 50, GRCh38 primary assembly; blacklist v2 | the same genome and GTF as rnaseq, so the assays of a project line up |
| `peaks.consensus_group_by`, `min_replicates` | `[cell_line]`, 2 | a region is kept when both replicates of one line called it, so a region open in only one line survives. With the default per-sample cut-off (MACS3 p < 0.01, as ENCODE), this replicate requirement is the reproducibility filter |
| `bigwig.group_by` | `[cell_line]` | one averaged track per line, besides the per-sample ones |
| `qc.somalier_sites` | the hg38 sites | turns on the identity check |
| `da.analyses` | `~ cell_line`, K562 vs GM12878 | the contrast; log2 fold change > 0 = more open in K562 |

For a project with treatments in several cell lines, the `da` section would rather be one
model per line, as in rnaseq:

```yaml
da:
  analyses:
    - name: per_line
      split_by: cell_line
      design: "~ condition"
      contrasts:
        - [condition, treated, control]
```

## 3. Run

```bash
./run.py -w all -c tests/tutorial/encode/config_tutorial.yaml -q 32 -n   # dry-run: checks the config, lists the jobs
./run.py -w all -c tests/tutorial/encode/config_tutorial.yaml -q 32
```

The dry-run catches configuration mistakes (a column the samplesheet lacks, a contrast level
no sample has) before anything runs. The first real run also builds the reference: the
bowtie2 index takes 1–2 hours and is built once per `reference.dir`.

On a SLURM cluster, either run the launcher inside one allocation:

```bash
#!/bin/bash
#SBATCH --job-name=atac_encode
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
cd /path/to/sm-core/pipelines/atacseq
conda activate sm-core
./run.py -w all -c tests/tutorial/encode/config_tutorial.yaml -q 32
```

or let it submit one job per rule (`-cl --scheduler slurm -qu <partition> -j 20`, see the
suite's [usage guide](../../../../../docs/usage.md#cluster-runs-pbs-or-slurm)).

## 4. Check the QC

In the order of [outputs.md](../../../docs/outputs.md#qc-what-to-look-at-first):

1. `results/qc/qc_metrics.png`: every library against ENCODE's ideal and acceptable tiers;
2. `results/qc/tss_enrichment.png` and `fragment_sizes.png`: a sharp TSS peak and a
   nucleosome ladder;
3. `results/qc/identity/somalier.html`: replicates related, the two lines not;
4. `results/qc/multiqc_report.html`: everything per tool.

## 5. Differential accessibility

`-w all` includes it because `da.analyses` is set; after changing only the `da` section,
`./run.py -w da ...` reruns just that part. In `results/da/cell_lines/all/`:

1. `pca`, `sample_distances`: the replicates together, the lines apart;
2. `size_factors`: the two normalisations, and each library's fraction of fragments in
   peaks ([outputs.md](../../../docs/outputs.md#differential-accessibility--w-da));
3. `K562_vs_GM12878/`: `volcano`, `ma`, `heatmap`, `peak_classes`; `da_results.tsv`, the
   significant peaks as BED for a genome browser, and `genes.tsv`, the per-gene table to set
   beside RNA-seq.

## 6. Validate against ENCODE

```bash
python3 tests/tutorial/encode/validate.py
```

1. **QC, the same reads processed twice.** The mitochondrial fraction, duplicate fraction,
   NRF and PBC1 must agree with ENCODE's within a small tolerance, since the filters are
   the same. FRiP and TSS enrichment are reported side by side: ENCODE computes them on
   its own peaks and with its own TSS list, so they differ in level.
2. **Peaks.** For each line, the share of ENCODE's conservative IDR regions (the
   reproducible core) found by sm-core's consensus peaks for that line (recall ≥ 0.90),
   and the share of sm-core's consensus peaks inside ENCODE's pseudoreplicated set
   (precision ≥ 0.80). ENCODE's file lists a region once per summit, so distinct regions
   are counted, and the ~600 per line inside the ENCODE blacklist, which sm-core removes
   by design, are left out.
3. **Biology.** Promoters of B-cell genes (CD19, CD79A/B, PAX5, MS4A1, BLK) more open in
   GM12878, and of erythroid genes (GATA1, TAL1, KLF1, HBG1/2, GYPA) more open in K562.
4. **Identity.** Replicates of one line are the same individual, and the two lines are not.
5. **Differential accessibility.** The lineage genes of check 3 as significant promoter
   peaks of the K562 vs GM12878 contrast: erythroid opened in K562, B-cell closed.

It prints a report and exits non-zero if a check fails.

## Results

Run of 4 November 2026: GRCh38 primary assembly, GENCODE 50, sm-core/atacseq 0.1.0.
`validate.py`: all checks passed.

**QC** (sm-core / ENCODE, the same reads):

| Library | Mito fraction | Duplicates | NRF | PBC1 | FRiP | TSS enrichment | Fragments after filtering |
|---|---|---|---|---|---|---|---|
| GM12878_1 | 0.125 / 0.139 | 0.277 / 0.295 | 0.798 / 0.797 | 0.803 / 0.805 | 0.66 / 0.62 | 33.9 / 29.6 | 22.6 M |
| GM12878_2 | 0.089 / 0.103 | 0.169 / 0.184 | 0.888 / 0.886 | 0.892 / 0.892 | 0.59 / 0.53 | 25.7 / 23.1 | 23.9 M |
| K562_1 | 0.095 / 0.121 | 0.230 / 0.255 | 0.824 / 0.817 | 0.830 / 0.826 | 0.46 / 0.39 | 25.7 / 21.9 | 24.6 M |
| K562_2 | 0.105 / 0.131 | 0.230 / 0.259 | 0.831 / 0.824 | 0.837 / 0.833 | 0.43 / 0.35 | 20.9 / 18.5 | 22.7 M |

The complexity metrics agree to the third decimal, and the mitochondrial and duplicate
fractions within 0.03. FRiP and TSS enrichment are higher than ENCODE's (different peaks
and TSS set) but rank the libraries identically. Every library is in ENCODE's ideal tier
for TSS enrichment and FRiP; all four are just under the 25 M fragments ENCODE asks for,
which `qc_metrics.png` shows in red.

**Peaks.** 91–119 k peaks per library; 118,444 consensus peaks, 76,434 of them reproduced
in GM12878 and 72,343 in K562.

| Line | Recall of ENCODE IDR regions | Precision vs ENCODE pseudoreplicated |
|---|---|---|
| GM12878 | 0.926 | 0.976 |
| K562 | 0.909 | 0.885 |

A stricter cut-off (MACS3 q < 0.01) gave a recall of about 0.81 and 0.73: the replicate requirement,
not the per-sample threshold, is what should do the filtering.

**Identity.** Replicates 0.99 related, the two lines −0.16 to −0.30. With ENCODE's
RNA-seq BAMs of the two lines added (`qc.identity_bams`: ENCFF028IUE, K562;
ENCFF450ZKQ, GM12878), each ATAC library matched its own line's RNA at 0.98–1.00, on
~220–310 shared sites, and the other line's at ≤ 0 (`identity_matches.tsv`).

**Differential accessibility**, K562 vs GM12878: 118,256 peaks tested, 40,474 more open in
K562 and 39,633 in GM12878 (padj < 0.05, |log2FC| ≥ 1); 6,707 genes with a promoter peak
opened in K562 and 4,794 closed. Two cell lines of different lineages differ at most of
their regulatory elements, so half the peaks changing is the expected scale here, not a
treatment's. All lineage promoters with a peak go the expected way, e.g. GATA1 +5.1, TAL1
+4.6, HBG1 +6.0 (erythroid, up); MS4A1 −6.8, PAX5 −4.3, BLK −5.7 (B cell, down). Among
the top-ranked peaks: BCR (amplified as BCR-ABL1 in K562), PRAME, and CDKN2B-AS1 (the
CDKN2A/B locus, deleted in K562). Library-size normalisation would shift every log2 fold
change by −0.42 (the gold line on the MA plot): GM12878 has more of its fragments in peaks
(0.59–0.68 vs 0.43–0.46), which median of ratios treats as depth. One model took 10 minutes
and 1.5 GB.
