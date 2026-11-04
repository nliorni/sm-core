# Installation

## Requirements

- Linux (x86-64).
- [conda](https://conda-forge.org/download/), Miniforge recommended. Snakemake uses it to
  build one environment per tool.
- Memory: about 16 GB for bowtie2, Picard and MACS3 on human, and 32 GB to build the
  bowtie2 index ([resources.md](resources.md)).
- Disk: about 4 GB for the human bowtie2 index, and roughly 1 GB of final BAM per
  10 million fragments.

## 1. The launcher environment

From the repository root, once for every sm-core pipeline:

```bash
conda env create -f environment.yml
conda activate sm-core
```

## 2. Conda for the per-rule environments

As for rnaseq ([rnaseq's installation.md](../../rnaseq/docs/installation.md#2-configure-conda-for-the-per-rule-environments)):
set strict channel priority and keep `~/.condarc` to conda-forge and bioconda. The
environments take 10–30 minutes to build the first time.

## 3. Check the installation

The [ENCODE tutorial](../tests/tutorial/encode/README.md) runs every rule on four real
libraries (two cell lines, two replicates each) at full depth against GRCh38, and checks the results against ENCODE's own
processing of the same reads:

```bash
cd pipelines/atacseq
tests/tutorial/encode/download_tutorial_data.sh          # ~19 GB of reads, plus the genome
./run.py -w all -c tests/tutorial/encode/config_tutorial.yaml -q 32
python3 tests/tutorial/encode/validate.py
```

It takes several hours, most of it building the bowtie2 index and aligning. A dry-run
(`-n`) after the download checks the configuration in seconds.

## 4. A reference

```bash
../rnaseq/resources/get_reference.sh --dest /refs/GRCh38_gencode50        # genome + GTF
resources/get_atac_resources.sh --build hg38 --dest /refs/GRCh38_gencode50 # blacklist + somalier sites
```

Then build the index once with `-w reference`. [resources.md](resources.md) has the
details.
