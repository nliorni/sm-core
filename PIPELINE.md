# PIPELINE.md

## What this is

`sm-core` is a clone-and-customize scaffold for Snakemake pipelines, factored
out of a run.py launcher + directory layout that was copy-pasted (and
independently drifted) across eight separate pipelines. This file documents
the conventions the scaffold assumes, so a new pipeline built from it starts
consistent with the family instead of reinventing each piece.

Replace this whole file with your real pipeline's own PIPELINE.md once you've
built something -- keep whichever conventions below still apply, drop the
rest.

## Layout

```
run.py                    Launcher CLI (see below) -- rarely needs edits
                           beyond WORKFLOWS and the {{PIPELINE_NAME}} strings
bootstrap.py               One-time script that fills those placeholders
Snakefile                  Top-level: output_dir, SAMPLES, scaled(), includes
config_template.yaml       Annotated config to copy per project
workflow/
  rules/*.smk               Rule definitions, one file per logical stage
  scripts/*.py               Rule `script:` targets (get `snakemake.*` injected)
  envs/*.yaml                 Conda env definitions referenced by `conda:`
tests/tutorial/             A runnable, self-contained demo (see its README)
logs/pbs/                   qsub stdout/stderr when run.py -cl is used (gitignored)
.smk_locks/                 Per-config lockfiles (gitignored, see below)
```

## The `run.py` launcher

A single argparse CLI wrapping `snakemake.snakemake(...)` directly (not the
`snakemake` shell command), covering:

> **Requires Snakemake 7.x.** The `snakemake.snakemake(...)` Python API
> `run.py` calls (`configfiles=`, `forceall=`, `dryrun=`, etc.) was removed
> in Snakemake 8 in favor of `snakemake.api`/`snakemake.cli`. A newer
> Snakemake will `import snakemake` fine and fail at the call site with
> `AttributeError: module 'snakemake' has no attribute 'snakemake'`. Install
> with `pip install "snakemake<8"`, or if you're on a cluster with several
> Snakemake conda envs already around, make sure `run.py` is launched from a
> 7.x one specifically.

- **Local or PBS execution** from the same invocation: `-q N` runs locally
  with N cores; adding `-cl -qu <queue> -j <max-concurrent>` instead submits
  each rule instance as its own `qsub` job. `-cl` needs its own PBS
  submit-command template (`build_cluster_cmd`) -- written for OpenPBS/PBS
  Pro syntax; if your cluster runs SLURM/LSF/etc., that's the one function to
  rewrite.
- **`WORKFLOWS`**: a name -> description dict at the top of the file, both
  documenting and validating `-w`/`--workflow` targets (`--allow-custom-target`
  bypasses validation for one-off rule names). Keep every real workflow you
  add wired to an actual rule/file target in the Snakefile, the way the
  `example` demo target is.
- **Per-config locking** (`acquire_project_lock`): Snakemake's own lock is
  directory-wide, so two runs against *different* configs in the same
  checkout would serialize through it. `run.py` instead takes an flock keyed
  on the config file's absolute path, under `.smk_locks/`, and passes
  `lock=False` to Snakemake. Held for the process's lifetime, released
  automatically on any exit (clean, crashed, or killed) -- no `--unlock`
  cleanup step needed for this lock specifically (Snakemake's own
  `--unlock` flag still exists for its other state).
- **Retry-scaled resources**: `-rt/--restart-times N` (default 2) lets
  Snakemake retry a job killed for exceeding its resources, with the
  Snakefile's `scaled()` helper multiplying `mem_mb`/`runtime` by the attempt
  number so each retry asks for more.
- Dry-run (`-n`), DAG export (`--dag`), Singularity bind paths (`--bind`,
  repeatable), custom resource limits (`--resources key=value ...`),
  colorized run-summary banner before every real invocation.

## Paired-input resolution convention

`samples:` in the config maps a sample name to a file *prefix*; a rule needs
`{prefix}_{unit[N]}{file_ext}`, where `unit` is a 2-element list of tokens
(e.g. `["R1_001", "R2_001"]`) and `file_ext` is the real extension in use.
Centralize this as one pair of helper functions in the Snakefile (`input_1`/
`input_2` in the shipped skeleton) rather than re-deriving the path inline in
each rule -- keeps the naming convention in exactly one place if it ever
needs to change.

## Bootstrapping a new pipeline from this template

```bash
git clone https://github.com/nliorni/sm-core.git my-new-pipeline
cd my-new-pipeline
./bootstrap.py MyNewPipeline
rm -rf .git && git init   # start this pipeline's own history
```

Then, in order:

1. Edit `WORKFLOWS` in `run.py` with your real workflow targets.
2. Replace `workflow/rules/example.smk` (and its `scripts/`/`envs/`) with
   real rules, reusing the `scaled()` and `input_1`/`input_2` conventions
   above.
3. Edit `config_template.yaml` for your project's actual config shape.
4. Update `tests/tutorial/` -- either keep it a fast synthetic smoke test
   (recommended: no tool/conda dependency, safe to run in CI) or swap in a
   `download_tutorial_data.sh` that fetches a small real dataset on demand
   (the pattern used by this family's DNA-seq pipelines, which don't vendor
   real data in git either).
5. Delete `bootstrap.py` once you no longer need it (optional).
