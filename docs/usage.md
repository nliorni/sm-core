# Running a pipeline

Every pipeline of the suite is launched the same way, by the same code
(`smcore/launcher.py`), so everything here applies to all of them. Two equivalent
spellings:

```bash
pipelines/<name>/run.py -w <target> -c <config> -q <cores>      # the pipeline's own launcher
./run.py <name> -w <target> -c <config> -q <cores>              # the suite dispatcher, from the repo root
```

`./run.py --list` lists the pipelines; `<launcher> --list-workflows` a pipeline's targets;
`<launcher> --help` every flag.

## The launcher environment

```bash
conda env create -f environment.yml     # once
conda activate sm-core
```

It holds only Snakemake 9.16.0 and the launcher's own dependencies. Analysis software is
not in it: each rule declares its own conda environment (`workflow/envs/*.yaml`, exact
versions pinned), which Snakemake builds on first use.

## Local runs

```bash
./run.py -w all -c config_my_project.yaml -q 32 -n     # dry-run: what would run, and why
./run.py -w all -c config_my_project.yaml -q 32        # run on 32 cores of this machine
```

Always dry-run first: it resolves the whole job graph, checks that every input the pipeline
cannot build exists, and validates the config — in seconds, without running anything.
`-p` prints each shell command; `--dag | dot -Tsvg > dag.svg` draws the graph.

## Cluster runs (PBS or SLURM)

```bash
./run.py -w all -c config_my_project.yaml -cl -qu workq -j 20                       # PBS: qsub
./run.py -w all -c config_my_project.yaml -cl --scheduler slurm -qu normal -j 20    # SLURM: sbatch
```

`-cl` submits every job (one rule on one sample) as its own cluster job instead of
running it locally: `-qu` is the PBS queue or SLURM partition, `-j` the maximum number of
jobs queued or running at once. The launcher itself stays running on the login node and
must survive until the end — start it inside `tmux`/`screen`, or submit the launcher
itself as a long, 1-core job.

Each job requests the threads, `mem_mb` and `runtime` its rule declares (pipelines set
these from their config); rules that declare nothing get 8 GB and 2 h. Scheduler stdout and
stderr go to `pipelines/<name>/logs/cluster/<rule>.<jobid>.{out,err}`; the tool's own log
is under `<output_dir>/logs/`.

| `--scheduler` | Submits with | Resources |
|---|---|---|
| `pbs` (default) | `qsub -q <queue>` | `-l select=1:ncpus=N:mem=Mmb -l walltime=MIN:00` |
| `slurm` | `sbatch --partition=<queue>` | `--cpus-per-task=N --mem=M --time=MIN` |

**On a SLURM cluster, use `--scheduler slurm`** even if a `qsub` command exists: that is
usually SLURM's PBS-compatibility wrapper (a Perl script — check with `file $(which qsub)`),
which translates `ncpus=N` into `--mincpus=N`, a node-selection constraint, not
`--cpus-per-task`. Under cgroup enforcement a 16-thread job can then be confined to fewer
CPUs than it runs threads. Another scheduler is one more template in
`build_cluster_cmd()` (`smcore/launcher.py`), for the whole suite.

**Under the hood.** Snakemake 8 moved cluster submission out of its core into executor
plugins; `-cl` uses `snakemake-executor-plugin-cluster-generic` (in `environment.yml`)
with the same submit template Snakemake 7's `--cluster` took. Job completion is detected
from each job's output files appearing; `--latency-wait` (60 s) allows for slow shared
filesystems.

## Retries with more memory

`-rt N` (default 2) retries a failed job up to N times. Rules that declare their resources
with `scaled()` ask for 2×, then 3× their memory and walltime on each retry, so a job
killed by the scheduler for exceeding its request comes back with more instead of failing
the run.

## Several projects at once

Run as many projects as you like from one checkout, each with its own config. Snakemake's
own lock would serialise them through one directory-wide lock; the launcher replaces it
with a lock per config file (`.smk_locks/`), so only two runs of the *same* config are
refused. That lock is released automatically however the run ends — no `--unlock` needed
after a crash (`-u` still exists, for Snakemake's own bookkeeping).

## Resuming and re-running

Re-running the same command resumes: finished outputs are kept. After a crash or a killed
job, add `-ri` (`--rerun-incomplete`) to redo jobs whose outputs were left half-written.

By default Snakemake also re-runs a job when its code, parameters, inputs or software
environment changed since it ran. To add samples to a finished project without re-running
everything after an unrelated change, `--rerun-triggers mtime` restricts this to
modification times — but then a genuine fix to a rule is not re-applied to old outputs
either.

## Software environments

Environments are built in `pipelines/<name>/.snakemake/conda/` by default. To build them
once for several checkouts, or on a faster disk, use `--conda-prefix /path/to/envs`.
Sites with their own pre-built environments can substitute them per tool in the config
(`conda_envs:`), without editing any rule. `--no-conda` uses whatever is on `PATH`.

For fast, reproducible solves, configure conda with strict channel priority and only
conda-forge and bioconda:

```bash
conda config --set channel_priority strict
```

## Working directory

The Snakemake working directory is the pipeline's own directory, not the shell's: relative
paths in a config (output, references, samplesheet) mean the same thing wherever the
launcher is started from. `-di` overrides this.

## Notifications

`notify_email:` in a config sends one e-mail at the end of the run, on success or failure,
through the system `mail` command.
