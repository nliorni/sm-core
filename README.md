# sm-core

> ### 📣 News!
>
> - **27 October 2026**: **rnaseq 0.2.0** adds differential expression in Python: PyDESeq2
>   for the contrasts in your config, with ORA and GSEA (MSigDB), validated against TaqMan
>   qPCR. All sm-core figures now share one publication style.
> - **22 October 2026**: the first sm-core pipeline is available:
>   **[rnaseq](pipelines/rnaseq/README.md)**, paired-end RNA-seq from FASTQ to STAR and
>   Salmon count matrices, validated on human SEQC/MAQC data. Please have a look!

A suite of Snakemake pipelines for computational multi-omics. It aims to do for Snakemake
what nf-core does for Nextflow: each pipeline is state-of-the-art, documented key by key,
validated on public data against a published truth set, and runs the same way as every
other pipeline in the suite.

## Available pipelines

sm-core is a collection: each pipeline covers one kind of data and lives in its own
directory under `pipelines/`, with its own documentation, tutorial and version. One
checkout and one environment run all of them.

| Pipeline | What it does | Version |
|---|---|---|
| [**rnaseq**](pipelines/rnaseq/README.md) | Paired-end RNA-seq: FASTQ → STAR and Salmon count matrices → differential expression, ORA and GSEA | 0.2.0 |

More pipelines will be added to this table as they are released; [CHANGELOG.md](CHANGELOG.md)
records each one.

## Using a pipeline

Every pipeline is used the same way. Once, set up the suite:

```bash
git clone https://github.com/nliorni/sm-core.git && cd sm-core
conda env create -f environment.yml          # Snakemake 9.16.0 + launcher, shared by all pipelines
conda activate sm-core
./run.py --list                              # the available pipelines
```

Then, for the pipeline you want (here `rnaseq`):

1. **Read its README**, [`pipelines/<name>/README.md`](pipelines/rnaseq/README.md). It has
   that pipeline's quick start, and links to its own `docs/`: installation, resources,
   configuration, outputs and methods.
2. **Get its resources** (a reference genome and annotation, for example). Each
   pipeline's `docs/resources.md` explains what it needs and where to get it.
3. **Copy its config template** and fill it in:
   `cp pipelines/rnaseq/config/config.yaml config_my_project.yaml`.
   Use absolute paths inside it: relative ones are read from the pipeline's directory,
   wherever you launch from.
4. **Pick a target and run it**: dry-run first (`-n`), then for real, locally or on a
   cluster:

```bash
./run.py rnaseq --list-workflows                                # its targets
./run.py rnaseq -w all -c config_my_project.yaml -q 32 -n       # dry-run
./run.py rnaseq -w all -c config_my_project.yaml -q 32          # run locally, 32 cores
./run.py rnaseq -w all -c config_my_project.yaml -cl -qu workq -j 20   # one cluster job per rule (PBS; --scheduler slurm for SLURM)
```

`./run.py rnaseq …` and `pipelines/rnaseq/run.py …` are the same command. The flags are
the same for every pipeline, documented once in [docs/usage.md](docs/usage.md). To check
your setup first, run the pipeline's small tutorial in `pipelines/<name>/tests/tutorial/`.

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
