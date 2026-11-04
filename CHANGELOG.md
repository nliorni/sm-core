# Changelog

Changes to the sm-core framework (`smcore/`, the launcher, templates). Each pipeline keeps
its own `pipelines/<name>/CHANGELOG.md`.

## [Unreleased]

### Added
- `smcore/rules/reads.smk`: lane merging and fastp, moved out of rnaseq so that atacseq
  uses the same rules. rnaseq's outputs and its rerun behaviour are unchanged.
- `smcore/plotting`: one figure style for every pipeline of the suite (semantic colours:
  red up, blue down, grey not significant; Arial-like type found on the system; PNG and
  PDF with editable text), adapted from Chen Liu's figures4papers. Also a matplotlib
  style file, `smcore/plotting/smcore.mplstyle`, for notebooks.
- `SMCORE_ROOT` in `smcore/common.smk`: lets scripts running in a rule's conda
  environment import shared sm-core code.
- `smcore/de.py` and `smcore/rules/de_design.smk`: the differential-analysis code of
  rnaseq's `-w de` (PyDESeq2 models, design checks, figures), moved out so that atacseq's
  `-w da` uses the same. `fit()` and `vst()` also take fixed size factors, which PyDESeq2
  0.5 has no option for; given PyDESeq2's own factors, the results are identical. rnaseq's
  DE tables are byte-identical before and after the move (yeast and seqc tutorials).

### Fixed
- The launcher makes sure Snakemake can run `conda`. `conda activate` provides a shell
  function, not a program, and some shell setups leave the base installation off `PATH`:
  Snakemake then failed with "conda: command not found" / "Error running conda info". The
  launcher now finds the base installation (from `CONDA_EXE`, or from the environment it
  runs in) and adds it to `PATH`, or stops at once with instructions.

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
