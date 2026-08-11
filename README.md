# sm-core

A clone-and-customize scaffold for Snakemake pipelines: a PBS launcher CLI
(`run.py`) plus directory conventions for `workflow/{rules,scripts,envs}`,
a `config_template.yaml`, and a `tests/tutorial/` demo.

## Quickstart

```bash
git clone https://github.com/nliorni/sm-core.git my-new-pipeline
cd my-new-pipeline

conda env create -f environment.yml
conda activate sm-core

# Prove it works out of the box:
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 -n   # dry-run
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2      # for real
cat tests/tutorial/output/summary.tsv

# Then make it yours:
./bootstrap.py MyNewPipeline
rm -rf .git && git init
```

See [PIPELINE.md](PIPELINE.md) for what `run.py` actually does (PBS
submission, per-config locking, retry-scaled resources, dry-run/DAG export)
and the conventions the scaffold assumes. See
[tests/tutorial/README.md](tests/tutorial/README.md) for the demo itself.

## What's in here

- **`run.py`** -- launcher CLI. `-q N` runs locally with N cores; `-cl -qu
  <queue> -j <jobs>` submits each rule as its own PBS job instead. Dry-run,
  DAG export, per-config file-locking, conda/Singularity toggles, retry with
  scaled resources. `./run.py --help` / `--list-workflows` for the rest.
- **`Snakefile`** + **`workflow/`** -- a minimal skeleton (`scaled()` retry
  helper, paired-input resolution convention, one demo rule chain) to delete
  and replace with real rules.
- **`config_template.yaml`** -- annotated config to copy per project.
- **`tests/tutorial/`** -- the runnable demo mentioned above.
- **`bootstrap.py`** -- fills in the pipeline name after cloning.
- **`environment.yml`** -- conda env with the exact Snakemake version `run.py` needs.

## Requirements

`conda env create -f environment.yml` installs everything: Python 3.11,
`colorama`, `pyyaml`, and Snakemake pinned to **7.32.4**.

That pin is not optional -- `run.py` drives the workflow through the
`snakemake.snakemake(...)` Python API, which Snakemake 8 removed in favor of
`snakemake.api`/`snakemake.cli`. A newer Snakemake will `import snakemake`
fine and only fail at run time, with `AttributeError: module 'snakemake' has
no attribute 'snakemake'`. If you're installing manually instead of via
`environment.yml` (e.g. an HPC cluster with several Snakemake environments
already around), make sure the one `run.py` runs from is 7.x.

A PBS/OpenPBS cluster if you use `-cl` (for other schedulers, see
`build_cluster_cmd` in `run.py`).

## License

MIT -- see [LICENSE](LICENSE).
