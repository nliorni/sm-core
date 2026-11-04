<h1 align="center">sm-core/atacseq</h1>

<p align="center">
  <em>Paired-end bulk ATAC-seq, from raw reads to consensus peaks and differential accessibility, checked against the ENCODE standards.</em>
</p>

<p align="center">
  <a href="https://github.com/nliorni/sm-core/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/nliorni/sm-core/actions/workflows/ci.yml/badge.svg?branch=main"></a>
  <img alt="Snakemake" src="https://img.shields.io/badge/snakemake-9.16.0-3C4B9A">
  <img alt="bowtie2" src="https://img.shields.io/badge/bowtie2-2.5.5-E4572E">
  <img alt="MACS3" src="https://img.shields.io/badge/MACS3-3.1.0-E4572E">
  <img alt="Reference" src="https://img.shields.io/badge/reference-GRCh38%20%7C%20GENCODE%2050-3C4B9A">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-green">
</p>

---

**sm-core/atacseq** is the ATAC-seq pipeline of the [sm-core](../../README.md) suite. It
takes paired-end FASTQ files to filtered BAMs, per-sample peaks, a **consensus peak set
with a fragment count matrix**, **differential accessibility** per contrast, signal tracks
for a genome browser, and a QC report that scores every library against the
[ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/). It follows ENCODE's
filtering and the nf-core/atacseq design, and shares its reference, samplesheet and
conventions with [sm-core/rnaseq](../rnaseq/README.md), so the two assays of one project
line up.

## Features

| Target | What it delivers | Main tools |
|---|---|---|
| `all` | Filtered BAMs; MACS3 peaks per sample; **consensus peaks, annotated to genes**; **peak count matrix**; bigWigs per sample and per group; ENCODE QC table and figures; sample identity check; MultiQC | fastp, bowtie2, Picard, samtools, MACS3, featureCounts, deepTools, somalier, MultiQC |
| `da` | **Differential accessibility** of the consensus peaks per contrast: results with nearest genes, significant peaks as BED, a per-gene summary, volcano / MA / heatmap / PCA; two normalisations compared | PyDESeq2 |
| `reference` | Genome, bowtie2 index, gene and TSS annotation, blacklist, built once | bowtie2, samtools |
| `trim` | Adapter- and quality-trimmed reads only | fastp |

- **ENCODE filtering.** Pairs that are properly paired, primary and MAPQ ≥ 30; Picard
  duplicates removed; mitochondrial reads removed after being counted; blacklisted peaks
  dropped.
- **QC scored, not just plotted.** TSS enrichment, FRiP, library complexity (NRF, PBC1,
  PBC2), alignment rate and depth are compared with ENCODE's ideal and acceptable tiers,
  per library, in `results/qc/atac_qc_summary.tsv` and in MultiQC. The fragment-size
  distribution (the nucleosome ladder) and the TSS profile are plotted for every sample.
- **A consensus peak set that respects the design.** A region is kept when enough
  replicates of *one* group called it (`peaks.consensus_group_by`, e.g. cell line ×
  condition), so a region open only under one treatment survives and a single-replicate
  peak does not. Peaks are called as ENCODE does (permissive per sample, reproducibility
  from replicates). Each peak is annotated with the genes whose promoter it covers, its
  nearest gene and the distance to its TSS, and promoter / genic / distal.
- **Differential accessibility with the design you give it.** The same analyses format as
  rnaseq's differential expression (one model per cell line or a shared one, paired
  designs, exclusions), on the same PyDESeq2 code. ATAC-seq has no settled normalisation,
  so both standard choices are computed every time and their disagreement is drawn.
- **Sample identity.** With `qc.somalier_sites`, somalier genotypes every library at
  ~17,000 common SNPs and reports the relatedness of every pair and the inferred sex.
  BAMs from other assays of the same cells (`qc.identity_bams`, e.g. rnaseq's STAR BAMs)
  can join the check, so a swapped or mislabelled sample shows up across assays.
- **nf-core-compatible samplesheets**, with lanes and runs merged and any extra column
  carried through as metadata, as in rnaseq.
- **The same reference as rnaseq**: GENCODE's primary-assembly genome and GTF, plus the
  ENCODE blacklist (`resources/get_atac_resources.sh`).
- **Laptop to HPC**: local runs or one cluster job per rule (`-cl`), with retries that ask
  for more memory.

## Quick start

```bash
# 1. Launcher environment (shared by every sm-core pipeline)
conda env create -f ../../environment.yml && conda activate sm-core

# 2. Genome + GTF (the same files as rnaseq), then the blacklist and somalier sites
../rnaseq/resources/get_reference.sh --dest /refs/GRCh38_gencode50
resources/get_atac_resources.sh --build hg38 --dest /refs/GRCh38_gencode50

# 3. Copy the template; set output_dir, samplesheet, the reference block and the groups
cp config/config.yaml config_my_project.yaml

# 4. Dry-run, then run (the first run also builds the bowtie2 index, ~1-2 h for human)
./run.py -w all -c config_my_project.yaml -q 32 -n
./run.py -w all -c config_my_project.yaml -q 32
```

## How it works

```
paired-end FASTQ (lanes/runs merged)
  │
  └─ fastp ............ adapter + quality trimming (Nextera found by mate overlap)
       └─ bowtie2 ...... --very-sensitive -X 2000
            └─ filter .. proper pairs, primary, MAPQ >= 30 (ENCODE)
                 └─ Picard MarkDuplicates ──→ complexity, mito fraction, fragment sizes,
                      │                       TSS enrichment, chrX/chrY (atac_qc)
                      └─ remove duplicates + chrM ──→ bam/<sample>.bam
                           ├─ MACS3 -f BAMPE ──→ peaks/<sample>/  (blacklist removed) ──→ FRiP
                           │     └─ consensus (>= N replicates of one group)
                           │           ├─ annotation: nearest gene, TSS distance, class
                           │           └─ featureCounts ──→ counts/peak_counts.tsv
                           │                 └─ PyDESeq2 (-w da) ──→ da/<analysis>/<split>/<contrast>/
                           ├─ deepTools ──→ bigwig/<sample>.bw, bigwig/groups/<group>.bw
                           └─ somalier ───→ qc/identity/  (optional)

reference (once): genome · bowtie2 index · genes + TSSs from the GTF · blacklist
```

## Documentation

| Document | Contents |
|---|---|
| [docs/installation.md](docs/installation.md) | Requirements and first run |
| [docs/resources.md](docs/resources.md) | Genome, annotation, blacklist, somalier sites; what gets built; disk and memory |
| [docs/usage.md](docs/usage.md) | Targets, samplesheets, groups, troubleshooting |
| [docs/configuration.md](docs/configuration.md) | Key-by-key config reference |
| [docs/outputs.md](docs/outputs.md) | Every output file and what to check in QC |
| [docs/methods.md](docs/methods.md) | Tools, versions, parameters, rationale and citations |
| [../../docs/usage.md](../../docs/usage.md) | The launcher flags shared by every sm-core pipeline |

### Tutorials

| Tutorial | Data | Exercises |
|---|---|---|
| [encode](tests/tutorial/encode/README.md) | ENCODE GM12878 and K562, two replicates each, full depth, GRCh38 | validation against ENCODE's own processing of the same reads: QC metrics, reproducible peaks, lineage promoters, identity |

## Requirements

- Linux, conda, and Snakemake via the suite's `environment.yml`.
- Human: about 16 GB of RAM for bowtie2 and Picard, about 4 GB of disk for the bowtie2
  index, and roughly 1 GB of final BAM per 10 million fragments. See
  [docs/resources.md](docs/resources.md).

## Citation

Cite the pipeline via the repository's [CITATION.cff](../../CITATION.cff), together with the
tools listed in [docs/methods.md](docs/methods.md#citations).

## License

MIT. See [LICENSE](../../LICENSE).
