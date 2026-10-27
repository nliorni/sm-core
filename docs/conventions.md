# Conventions

How the suite is put together, and what a pipeline must do to belong to it. The aim is
the one that made nf-core work: a user who has run one pipeline of the suite can run any
other, and a developer who has read one can read any other.

## Layout

```
sm-core/
├── run.py                   suite dispatcher: ./run.py <pipeline> <args>
├── bootstrap.py             new pipeline from templates/pipeline/
├── environment.yml          the one launcher environment, for every pipeline
├── smcore/                  shared code, versioned with the suite
│   ├── launcher.py          the CLI every pipeline's run.py hands over to
│   └── common.smk           Snakemake helpers every Snakefile includes first
├── templates/pipeline/      the skeleton bootstrap.py copies (runnable toy demo)
├── pipelines/
│   └── <name>/              one directory per pipeline -- self-contained:
│       ├── run.py           ~30 lines: name, version, targets -> smcore.launcher
│       ├── Snakefile        includes, then the workflow targets
│       ├── config/          config.yaml: the annotated per-project template
│       ├── workflow/
│       │   ├── rules/       one file per stage; common.smk = pipeline helpers
│       │   ├── envs/        one conda env per tool, exact versions pinned
│       │   ├── scripts/     script: targets and helpers
│       │   └── schemas/     config.schema.yaml: validation + defaults
│       ├── tests/
│       │   ├── ci/          placeholder inputs for DAG dry-runs in CI
│       │   └── tutorial/    real datasets + truth sets, downloaded on demand
│       ├── docs/            installation, resources, usage, configuration, outputs, methods
│       ├── examples/        excerpts of real outputs from the tutorials
│       ├── README.md, CHANGELOG.md
└── docs/                    suite-level: these conventions, launcher usage
```

## Why a monorepo with a shared core

sm-core began as a clone-and-customise template. Its launcher had been copied into eight
pipelines and had drifted in each: a fix made in one never reached the others. Copying a
better template does not fix that; sharing code does. So:

- **What every pipeline does identically lives once**, in `smcore/`: launching, PBS
  submission, locking, retries, environment overrides, sample-input resolution,
  notifications. A pipeline's `run.py` only *describes* the pipeline.
- **What is specific to a pipeline lives only in its directory**, which stays readable on
  its own, as gVAP is.
- **Pipelines are versioned independently** (their `run.py` and `CHANGELOG.md`) but tested
  together: a change to `smcore/` runs every pipeline's CI.

nf-core splits the same three concerns across repositories (`nf-core/tools`,
`nf-core/modules`, one repo per pipeline). With a handful of pipelines maintained by one
group, one repository gives the same separation without cross-repository releases. If a
pipeline later needs its own release cycle, its directory moves out unchanged except for
the two `include`/`sys.path` lines that point at `smcore/`.

### Shared rules: not yet

nf-core also shares *rules* (`nf-core/modules`). sm-core deliberately does not, yet: a rule
is shared when a second pipeline needs it, not before, because a shared rule must be
parameterised for every caller and is harder to read than a local one. When that happens,
it goes in `smcore/modules/<tool>.smk` and is used with Snakemake's `module`/`use rule …
with:` mechanism.

## Rules of the house

These are the conventions gVAP established, applied suite-wide.

1. **Every path, threshold and resource comes from the config**, never hard-coded. The
   schema (`workflow/schemas/config.schema.yaml`) declares each key's type and default;
   `config/config.yaml` documents each one in place, in numbered `SECTION`s, with
   `>>> EDIT PER PROJECT <<<` on what a user must change.
2. **Software by name**: `conda: env("<tool>")`, never a literal path, so sites can
   substitute environments (`conda_envs:`). One env per tool, exact versions.
3. **Scripts by `{SCRIPTS}/<name>.py`** (or `script:` relative to the rule file), so rules
   work from any working directory.
4. **Resources with `scaled()`**: `mem_mb=scaled(N)`, `runtime=scaled(M)`, so retries ask
   for more.
5. **Every rule logs** to `<output_dir>/logs/<tool>/`; heavy rules also write a
   `benchmark:` to `<output_dir>/benchmarks/`.
6. **Inputs from one place**: sample FASTQs come from `load_samples()`
   (`smcore/common.smk`), never re-derived inside a rule.
7. **Optional stages** are chosen in small input functions in `common.smk`, not by making
   rules conditional, so the DAG stays readable and `-n` always tells the truth. Checks that
   depend on the requested target (e.g. "this target needs X") live in the target's input
   function, so they fire only when that target is asked for.
8. **Fail early, name the key**: config problems raise at parse time with a message that
   says which key to change.
9. **Comment why, not what.**
10. **Prefer `script:` for a rule's own Python/R code.** Snakemake re-runs a job when the
    file named in `script:` changes, but not when a script called from a `shell:` command
    changes. After editing one of those, delete the rule's outputs, *and* those of the
    target that collects them: Snakemake does not re-run a job just because an
    intermediate input was deleted.
11. **Figures through `smcore.plotting`.** Every plot a pipeline draws calls
    `apply_style()`, takes colours from it (`UP`, `DOWN`, `NS`, `CATEGORICAL`,
    `DIVERGING`) and is written with `save()`, as PNG and PDF. The colour roles are the
    same everywhere: red up, blue down, grey not significant, and blue then red for the
    first two groups. Scripts import it after `sys.path.insert(0, snakemake.params.smcore)`,
    with the rule passing `smcore=SMCORE_ROOT`; their environment needs only matplotlib.
    QC that MultiQC already plots stays in MultiQC.

## What a pipeline ships

- `run.py`, `Snakefile`, `config/config.yaml`, `workflow/schemas/config.schema.yaml`.
- `tests/ci/`: a `make_ci_inputs.sh` that creates empty placeholders and a `config_ci.yaml`
  — CI dry-runs every target on them, and asserts the structural promises the pipeline
  makes (e.g. "`-w salmon` never runs STAR").
- `tests/tutorial/<name>/`: at least one **real** dataset, from published raw reads, with a
  `download_tutorial_data.sh`, a ready config, and a README that validates the output
  against a **published truth set** — not only "it ran".
- `docs/` in the gVAP shape: `installation`, `resources` (where every reference comes
  from, sizes, memory), `usage`, `configuration` (key by key), `outputs` (every file, and
  which to use), `methods` (tools, versions, citations).
- `examples/`: excerpts of real tutorial outputs, so the repository shows what the
  pipeline produces without running it.
- `CHANGELOG.md` (Keep a Changelog), a version in `run.py`.

## Adding a pipeline

```bash
./bootstrap.py chipseq "ChIP-seq: FASTQ to peaks"
cd pipelines/chipseq
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 --no-conda   # runs already
```

Then replace the toy rules, write the schema and config, add the CI job to
`.github/workflows/ci.yml` (copy `rnaseq-checks`), and add the pipeline to the
"Available pipelines" table in the top-level `README.md`. When it is released, add a
dated line at the top of that README's News! box, newest first.

## Snakemake version

`environment.yml` pins Snakemake **9.16.0**. The launcher drives Snakemake through its
Python API (`snakemake.api`), which only exists from 8.0 and changed shape again in 9; pin
exactly. One known 9.16.0 issue is worked around: `DAGApi.printdag()` returns `None`
instead of the DOT graph, so `--dag` shells out to the `snakemake` CLI instead
(`run_dag_export()` in `smcore/launcher.py`).
