# Changelog

Changes to the sm-core framework (`smcore/`, the launcher, templates). Each pipeline keeps
its own `pipelines/<name>/CHANGELOG.md`.

## [0.2.0] - 2026-10-20

### Changed
- **From template to suite.** The clone-and-customise scaffold is now a monorepo: shared
  code in `smcore/`, one directory per pipeline in `pipelines/`, the scaffold in
  `templates/pipeline/`. `bootstrap.py <name> "<description>"` creates
  `pipelines/<name>/` instead of renaming a clone.
- The launcher (`run.py`) moved to `smcore/launcher.py` and takes a `Pipeline`
  description, so every pipeline's `run.py` is a ~30-line stub and launcher fixes reach
  all pipelines at once. Behaviour is unchanged: local or PBS execution, per-config
  locking, retry-scaled resources, dry-run, DAG export.
- The Snakemake working directory defaults to the pipeline's own directory, so relative
  config paths resolve the same way from anywhere.
- `smcore/common.smk` holds the shared Snakemake helpers: `scaled()`, `env()` with
  per-site `conda_envs:` overrides, `load_samples()` (FASTQ prefixes or an
  nf-core-style samplesheet, lanes merged), completion e-mails.

### Added
- `./run.py <pipeline> …` dispatcher at the repository root; `--version`.
- `pipelines/rnaseq` 0.1.0.
- Suite docs: `docs/conventions.md`, `docs/usage.md`.

## [0.1.0]

Initial scaffold: PBS launcher (`run.py`) on the Snakemake 9 API, directory conventions,
toy tutorial.
