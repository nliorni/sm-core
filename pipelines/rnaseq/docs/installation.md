# Installation

## Requirements

- Linux (x86-64). macOS works for the Salmon branch and the yeast tutorial; STAR's human
  index needs more memory than most laptops have.
- [conda](https://conda-forge.org/download/) — Miniforge recommended. Snakemake uses it to
  build one environment per tool.
- Memory: about 40 GB for STAR on human or mouse (index build and alignment), and about
  30 GB for Salmon. [resources.md](resources.md#memory-and-disk-at-a-glance-human) has the
  measured figures.
- Disk: about 56 GB for a human reference with both indices, and about 1.5 GB of BAM per
  10 million read pairs.
- A PBS or SLURM scheduler, only for `-cl` (`--scheduler slurm` for SLURM).

No container runtime is needed: every tool is installed from bioconda/conda-forge at a
pinned version.

## 1. The launcher environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate sm-core
```

This installs Snakemake 9.16.0 and the launcher's dependencies only, and is shared by every
pipeline of the suite.

## 2. Configure conda for the per-rule environments

Snakemake builds each tool's environment (`workflow/envs/*.yaml`) the first time a rule
needs it. They only list conda-forge and bioconda, but conda also adds every channel in
your `~/.condarc`. Extra channels make the solves slow, and can make them wrong. Before the
first run:

```bash
conda config --set channel_priority strict
conda config --show channels      # ideally: conda-forge, bioconda
```

The environments take 10–30 minutes to build the first time, and are reused after that.
To build them once for several checkouts or users, keep them in one place:

```bash
./run.py -w all -c my.yaml -q 16 --conda-prefix /shared/conda/sm-core-rnaseq
```

Sites with their own pre-built environments can point the pipeline at them per tool
(`conda_envs:` in the config; see [configuration.md](configuration.md#conda_envs--optional)).

## 3. Check the installation

The yeast tutorial runs every rule of both branches in a few minutes, once the
environments exist:

```bash
cd pipelines/rnaseq
tests/tutorial/yeast/download_tutorial_data.sh
./run.py -w all -c tests/tutorial/yeast/config_tutorial.yaml -q 4
tests/tutorial/yeast/check_outputs.sh
```

## 4. A reference

```bash
resources/get_reference.sh --dest /path/to/references/GRCh38_gencode50     # human, GENCODE 50
```

This downloads a matched genome + annotation, verifies the checksums and prints the
`reference:` block for the config. Then build the indices once with `-w reference`.
[resources.md](resources.md) covers mouse, other species, spike-ins and the reasons for
each choice.
