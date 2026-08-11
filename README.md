# sm-core

A clone-and-customize scaffold for Snakemake pipelines: a battle-tested PBS
launcher CLI (`run.py`) plus the directory conventions (`workflow/{rules,
scripts,envs}`, `config_template.yaml`, `tests/tutorial/`, `PIPELINE.md`)
that grew out of copy-pasting the same launcher across several independent
genomics pipelines. This repo pulls that shared 85-90% back into one place so
the next pipeline starts from it instead of another copy-paste.

It ships with a tiny synthetic demo (`tests/tutorial/`) so the template runs
end-to-end right after cloning, with zero bioinformatics dependencies.

## Quickstart

```bash
git clone https://github.com/nliorni/sm-core.git my-new-pipeline
cd my-new-pipeline
pip install "snakemake<8" colorama pyyaml   # run.py needs Snakemake 7.x -- see Requirements below

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

## Requirements

Python 3.9+, `colorama`, `pyyaml`, and **Snakemake 7.x specifically**
(`pip install "snakemake<8"`). `run.py` drives the workflow through the
`snakemake.snakemake(...)` Python API, which Snakemake 8 removed in favor of
`snakemake.api`/`snakemake.cli` -- installing a newer Snakemake will import
fine but fail at run time with `AttributeError: module 'snakemake' has no
attribute 'snakemake'`. If you're on an HPC cluster with several Snakemake
conda envs around, make sure the one `run.py` runs from is a 7.x one.

A PBS/OpenPBS cluster if you use `-cl` (for other schedulers, see
`build_cluster_cmd` in `run.py`).

## License

MIT -- see [LICENSE](LICENSE).
