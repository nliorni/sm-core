# sm-core

A suite of Snakemake pipelines for computational multi-omics. It aims to do for Snakemake
what nf-core does for Nextflow: each pipeline is state-of-the-art, documented key by key,
validated on public data against a published truth set, and runs the same way as every
other pipeline in the suite.

## Pipelines

| Pipeline | What it does | Status |
|---|---|---|
| [rnaseq](pipelines/rnaseq/README.md) | Paired-end RNA-seq: FASTQ → STAR and Salmon count matrices | 0.1.0 |

## Quick start

```bash
git clone https://github.com/nliorni/sm-core.git && cd sm-core
conda env create -f environment.yml          # Snakemake 9.16.0 + launcher; once
conda activate sm-core

./run.py --list                              # the pipelines
./run.py rnaseq --list-workflows             # one pipeline's targets
./run.py rnaseq -w all -c my_config.yaml -q 32 -n     # dry-run
./run.py rnaseq -w all -c my_config.yaml -cl -qu workq -j 20   # one cluster job per rule (PBS; --scheduler slurm for SLURM)
```

`./run.py rnaseq …` and `pipelines/rnaseq/run.py …` are the same command.

## What every pipeline shares

- **One launcher** (`smcore/launcher.py`): local runs (`-q N`) or one cluster job per rule
  (`-cl`), dry-runs, DAG export, automatic retries with more memory (`-rt`), and per-project
  locking, so several projects can run from one checkout. Documented once, in
  [docs/usage.md](docs/usage.md).
- **One set of Snakemake helpers** (`smcore/common.smk`): paired-end inputs from an
  nf-core-style samplesheet or from FASTQ prefixes, with lanes merged; per-tool conda
  environments that a site can substitute; retry-scaled resources; completion e-mails.
- **One shape**: `run.py`, `Snakefile`, `config/config.yaml` with a validating schema,
  `workflow/{rules,envs,scripts,schemas}`, `tests/{ci,tutorial}`, `docs/`, `examples/`.
  See [docs/conventions.md](docs/conventions.md).

## Repository layout

```
run.py              suite dispatcher
bootstrap.py        start a new pipeline from templates/pipeline/
environment.yml     the launcher environment, shared by all pipelines
smcore/             shared launcher and Snakemake helpers
pipelines/<name>/   one self-contained pipeline per directory
templates/pipeline/ the skeleton bootstrap.py copies; runs as a toy demo
docs/               suite-level docs: conventions, launcher usage
```

## Adding a pipeline

```bash
./bootstrap.py chipseq "ChIP-seq: FASTQ to peaks"
```

This creates `pipelines/chipseq/`, which already runs through the shared launcher. The
checklist for turning it into a real pipeline is in [docs/conventions.md](docs/conventions.md#adding-a-pipeline).

## Requirements

Linux, conda, and PBS or SLURM for `-cl` (see
[docs/usage.md](docs/usage.md#cluster-runs-pbs-or-slurm)). Each pipeline's own `docs/installation.md`
lists its memory and disk needs.

## Citation

Cite the suite through [CITATION.cff](CITATION.cff), together with the tools listed in the
pipeline's `docs/methods.md`.

## License

MIT. See [LICENSE](LICENSE).
