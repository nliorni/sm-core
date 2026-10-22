"""
The sm-core launcher: one CLI shared by every pipeline in the suite.

Each pipeline's own run.py is a few lines that describe the pipeline (name,
version, workflow targets) and hand over to main() here, so local runs, PBS
submission, locking, retries and dry-runs behave identically across the suite
and are fixed in one place.

    from smcore.launcher import Pipeline, main
    PIPELINE = Pipeline(name="rnaseq", version="0.1.0", ...)
    raise SystemExit(main(PIPELINE))
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from colorama import Fore, Style, init
from snakemake.api import (
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    OutputSettings,
    ResourceSettings,
    SnakemakeApi,
    WorkflowSettings,
)
from snakemake.resources import DefaultResources
from snakemake.settings.enums import RerunTrigger

from smcore import __version__ as SMCORE_VERSION

try:
    import yaml
except ImportError:
    yaml = None


init(autoreset=True)

# PBS resource fallback for rules that don't declare their own `resources:
# mem_mb`/`runtime` -- only used when --cluster is passed; local (non-cluster)
# runs are unaffected. Tune for your own cluster's typical job size.
CLUSTER_DEFAULT_MEM_MB = 8000
CLUSTER_DEFAULT_RUNTIME_MIN = 120


@dataclass(frozen=True)
class Pipeline:
    """What the launcher needs to know about one pipeline of the suite."""

    name: str
    version: str
    description: str
    # High-level workflow target -> one-line description, shown by
    # --list-workflows and validated against by default (pass
    # --allow-custom-target to bypass). Every key must be a real rule in the
    # pipeline's Snakefile.
    workflows: Dict[str, str]
    # The pipeline's own directory: holds its Snakefile, and is the default
    # Snakemake working directory (see run_snakemake()).
    directory: Path


def color(text: str, colour: str) -> str:
    return f"{colour}{text}{Style.RESET_ALL}"


def die(message: str, exit_code: int = 1) -> None:
    sys.stderr.write(color(f"ERROR: {message}\n", Fore.RED))
    raise SystemExit(exit_code)


def warn(message: str) -> None:
    sys.stderr.write(color(f"WARNING: {message}\n", Fore.YELLOW))


def info(message: str) -> None:
    sys.stderr.write(color(f"{message}\n", Fore.CYAN))


def success(message: str) -> None:
    sys.stderr.write(color(f"{message}\n", Fore.GREEN))


def logo(pipeline: Pipeline) -> str:
    title = f"sm-core/{pipeline.name} v{pipeline.version}"
    width = max(len(title), len(pipeline.description)) + 4
    return (
        f"\n{Fore.CYAN}    ┌{'─' * width}┐\n"
        f"    │  {Fore.GREEN}{title.ljust(width - 2)}{Fore.CYAN}│\n"
        f"    │  {Fore.YELLOW}{pipeline.description.ljust(width - 2)}{Fore.CYAN}│\n"
        f"    └{'─' * width}┘{Style.RESET_ALL}\n"
    )


def acquire_project_lock(pipeline: Pipeline, configfile: Path):
    # Snakemake's own default locking locks the whole *working directory*
    # (.snakemake/locks/), not just the project a given invocation actually
    # touches -- every project config of a pipeline shares that one directory,
    # so two concurrent runs against DIFFERENT configs would serialize every
    # scheduling step through that one shared lock.
    #
    # Fix: take our own lock keyed on the *config file's* absolute path
    # (i.e. per-project, not per-directory) via flock, then pass lock=False to
    # bypass Snakemake's own directory-wide one entirely. flock is held for
    # this process's entire lifetime and is released automatically by the
    # kernel on any exit (clean, crashed, or killed -- no stale-lock cleanup
    # step needed, unlike Snakemake's own file-based lock which requires an
    # explicit --unlock after a killed process).
    #
    # Still refuses to run the SAME config twice concurrently -- that
    # protection is preserved, just scoped to the one project actually at
    # risk instead of the whole pipeline.
    lock_dir = pipeline.directory / ".smk_locks"
    lock_dir.mkdir(exist_ok=True)
    key = hashlib.sha1(str(configfile.resolve()).encode()).hexdigest()[:16]
    lock_path = lock_dir / f"{key}.lock"
    fh = open(lock_path, "w")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        die(
            f"another run.py invocation is already active for config "
            f"'{configfile}' (lock: {lock_path}). Refusing to start a second "
            "one against the same project -- two concurrent invocations of "
            "the same config can race to declare/clean up the same output "
            "files."
        )
    fh.write(str(os.getpid()))
    fh.flush()
    return fh  # caller must keep this referenced for the process's lifetime


def load_config(configfile: Path) -> Dict[str, Any]:
    if yaml is None:
        warn("PyYAML is not installed; config summary will be skipped.")
        return {}

    try:
        with configfile.open("r") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        warn(f"Could not parse config file for summary: {exc}")
        return {}

    return data or {}


def count_samples(config: Dict[str, Any]) -> str:
    samples = config.get("samples")
    if isinstance(samples, (dict, list)):
        return str(len(samples))
    if config.get("samplesheet"):
        return f"from {config['samplesheet']}"
    return "0"


def print_workflows(pipeline: Pipeline) -> None:
    print(logo(pipeline))
    print(color("Available workflow targets:", Fore.GREEN))
    print()
    for name, description in pipeline.workflows.items():
        print(f"  {color(name, Fore.YELLOW):<20} {description}")
    print()


def parse_resources(resource_args: Optional[List[str]]) -> Optional[Dict[str, int]]:
    if not resource_args:
        return None

    resources: Dict[str, int] = {}

    for item in resource_args:
        if "=" not in item:
            die(f"Invalid resource format: {item}. Use key=value, e.g. igv=1 mem_mb=64000")

        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()

        if not key:
            die(f"Invalid resource key in: {item}")

        try:
            resources[key] = int(value)
        except ValueError:
            die(f"Resource value must be an integer in: {item}")

    return resources


def build_cluster_cmd(scheduler: str, queue: str, logdir: Path) -> str:
    # Braces doubled ({{...}}) are Snakemake's own per-job placeholders and
    # must survive this .format() call untouched; only {queue}/{logdir} are
    # substituted here. The template is handed to
    # snakemake-executor-plugin-cluster-generic (build_cluster_executor_settings):
    # Snakemake 8+ moved generic shell-command cluster submission out of core
    # into that plugin, with the same placeholders as Snakemake 7's --cluster.
    # Another scheduler (LSF, SGE) is one more template here.
    if scheduler == "slurm":
        # Native sbatch. --cpus-per-task is what SLURM binds the job's threads
        # to; SLURM's PBS-compatibility `qsub` wrapper translates ncpus= into
        # --mincpus instead, a node-selection constraint, which under cgroup
        # enforcement can leave a multi-threaded job fewer CPUs than it runs
        # threads. --mem is in MB, --time in minutes.
        template = (
            "sbatch --parsable "
            "--partition={queue} "
            "--nodes=1 --ntasks=1 --cpus-per-task={{threads}} "
            "--mem={{resources.mem_mb}} "
            "--time={{resources.runtime}} "
            "--job-name=smk.{{rule}}.{{jobid}} "
            "--output={logdir}/{{rule}}.{{jobid}}.out "
            "--error={logdir}/{{rule}}.{{jobid}}.err"
        )
    else:
        # OpenPBS / PBS Pro: `-l select=1:ncpus=X:mem=Ymb`, walltime as
        # MM:SS -- `{runtime}:00` is `runtime` minutes.
        template = (
            "qsub -q {queue} "
            "-l select=1:ncpus={{threads}}:mem={{resources.mem_mb}}mb "
            "-l walltime={{resources.runtime}}:00 "
            "-N smk.{{rule}}.{{jobid}} "
            "-o {logdir}/{{rule}}.{{jobid}}.out "
            "-e {logdir}/{{rule}}.{{jobid}}.err"
        )
    return template.format(queue=queue, logdir=logdir)


def build_singularity_args(args: argparse.Namespace) -> str:
    pieces: List[str] = []

    if args.singularity_args:
        pieces.append(args.singularity_args.strip())

    for bind_path in args.bind:
        pieces.append(f"--bind {bind_path}")

    return " ".join(pieces).strip()


def print_run_summary(
    pipeline: Pipeline,
    args: argparse.Namespace,
    snakefile: Path,
    configfile: Path,
    workdir: Path,
    config: Dict[str, Any],
    singularity_args: str,
    resources: Optional[Dict[str, int]],
) -> None:
    output_dir = config.get("output_dir", "not found in config")

    print(logo(pipeline))
    print(color("Run summary", Fore.GREEN))
    print(color("-----------", Fore.GREEN))
    print(f"{color('Snakefile:', Fore.YELLOW)}       {snakefile}")
    print(f"{color('Config:', Fore.YELLOW)}          {configfile}")
    print(f"{color('Targets:', Fore.YELLOW)}         {', '.join(args.workflow or [])}")
    print(f"{color('Cores:', Fore.YELLOW)}           {args.cores}")
    print(f"{color('Workdir:', Fore.YELLOW)}         {workdir}")
    print(f"{color('Output dir:', Fore.YELLOW)}      {output_dir}")
    print(f"{color('Samples:', Fore.YELLOW)}         {count_samples(config)}")

    print(f"{color('Dry run:', Fore.YELLOW)}         {args.dry_run}")
    print(f"{color('Use conda:', Fore.YELLOW)}       {not args.no_conda}")
    print(f"{color('Use Singularity:', Fore.YELLOW)} {not args.no_singularity}")

    if singularity_args:
        print(f"{color('Singularity args:', Fore.YELLOW)} {singularity_args}")

    if resources:
        print(f"{color('Resources:', Fore.YELLOW)}       {resources}")

    if args.forceall:
        print(color("Force all jobs:    yes", Fore.YELLOW))

    if args.rerun_incomplete:
        print(color("Rerun incomplete:  yes", Fore.YELLOW))

    if args.keep_going:
        print(color("Keep going:        yes", Fore.YELLOW))

    if args.cluster:
        print(f"{color('Cluster:', Fore.YELLOW)}         {args.scheduler} (queue={args.queue}, max concurrent jobs={args.jobs})")
        print(f"{color('Restart times:', Fore.YELLOW)}   {args.restart_times}")

    if args.rerun_triggers:
        print(f"{color('Rerun triggers:', Fore.YELLOW)}  {', '.join(args.rerun_triggers)} (default: mtime, params, input, software-env, code)")

    print()


def validate_args(pipeline: Pipeline, args: argparse.Namespace, snakefile: Path, configfile: Path) -> None:
    if not snakefile.exists():
        die(f"Cannot find Snakefile: {snakefile}")

    if not configfile.exists():
        die(f"Cannot find config file: {configfile}")

    if not args.unlock and not args.workflow:
        die(f"No workflow target provided. Use -w {', '.join(pipeline.workflows.keys())}.")

    if args.workflow and not args.allow_custom_target:
        invalid = [target for target in args.workflow if target not in pipeline.workflows]
        if invalid:
            valid = ", ".join(pipeline.workflows.keys())
            die(
                f"Unknown workflow target(s): {', '.join(invalid)}. "
                f"Valid targets are: {valid}. "
                "Use --allow-custom-target to run any rule name manually."
            )


RERUN_TRIGGER_BY_NAME = {
    "mtime": RerunTrigger.MTIME,
    "params": RerunTrigger.PARAMS,
    "input": RerunTrigger.INPUT,
    "software-env": RerunTrigger.SOFTWARE_ENV,
    "code": RerunTrigger.CODE,
}


def rerun_triggers_from_names(names: Optional[List[str]]) -> Optional[Set[RerunTrigger]]:
    if not names:
        return None
    return {RERUN_TRIGGER_BY_NAME[name] for name in names}


def build_deployment_method(args: argparse.Namespace) -> Set[DeploymentMethod]:
    methods: Set[DeploymentMethod] = set()
    if not args.no_conda:
        methods.add(DeploymentMethod.CONDA)
    if not args.no_singularity:
        methods.add(DeploymentMethod.APPTAINER)
    return methods


def build_cluster_executor_settings(pipeline: Pipeline, args: argparse.Namespace):
    # Snakemake 8+ pulled generic shell-command cluster submission out of
    # core into its own plugin -- this is the direct successor to the old
    # `cluster=` kwarg, not a new/different mechanism: same submit-command
    # template (build_cluster_cmd), same placeholder substitution, same
    # "watch for output files" completion detection (no status_cmd is set
    # below, matching Snakemake 7's behavior here).
    try:
        from snakemake_executor_plugin_cluster_generic import (
            ExecutorSettings as ClusterGenericExecutorSettings,
        )
    except ImportError:
        die(
            "The 'snakemake-executor-plugin-cluster-generic' package is required "
            "for -cl/--cluster (see environment.yml). Install it, or drop -cl to "
            "run locally instead."
        )

    cluster_logdir = pipeline.directory / "logs" / "cluster"
    cluster_logdir.mkdir(parents=True, exist_ok=True)
    return ClusterGenericExecutorSettings(
        submit_cmd=build_cluster_cmd(args.scheduler, args.queue, cluster_logdir)
    )


def run_dag_export(snakefile: Path, configfile: Path, workdir: Path, args: argparse.Namespace) -> int:
    # dag_api.printdag() -- the direct Python-API equivalent of --dag --
    # hits a reproducible upstream bug in Snakemake 9.16.0: dag.py's
    # Dag.__str__ compares the PrintDag enum against str(PrintDag.DOT),
    # which is never equal, so it silently returns None instead of the DOT
    # graph and print() blows up. Shelling out to the snakemake CLI itself
    # sidesteps it (the CLI's own --dag hits the same underlying code but
    # through a working path).
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(configfile),
        "--directory",
        str(workdir),
        "--cores",
        str(args.cores),
        "--dag",
        "dot",
        *(args.workflow or []),
    ]
    return subprocess.run(cmd).returncode


def run_snakemake(pipeline: Pipeline, args: argparse.Namespace) -> int:
    snakefile = Path(args.snakefile).resolve() if args.snakefile else pipeline.directory / "Snakefile"
    configfile = Path(args.configfile).resolve()
    # The pipeline's own directory, not the caller's: relative paths in a
    # project config (output_dir, references, samplesheet) then mean the same
    # thing wherever ./run.py is launched from, and each pipeline keeps its
    # Snakemake state (.snakemake/, conda envs) to itself. Absolute paths in
    # the config, as the templates recommend, are unaffected.
    workdir = Path(args.directory).resolve() if args.directory else pipeline.directory

    validate_args(pipeline, args, snakefile, configfile)

    config = load_config(configfile)
    resources = parse_resources(args.resources)
    singularity_args = build_singularity_args(args)

    # Keep a reference for the life of this process -- see
    # acquire_project_lock()'s comment for why. Bypasses Snakemake's own
    # directory-wide lock (lock=False below) in favor of this per-config one.
    _project_lock = acquire_project_lock(pipeline, configfile)

    print_run_summary(
        pipeline=pipeline,
        args=args,
        snakefile=snakefile,
        configfile=configfile,
        workdir=workdir,
        config=config,
        singularity_args=singularity_args,
        resources=resources,
    )

    if args.dag:
        return run_dag_export(snakefile, configfile, workdir, args)

    resource_settings = ResourceSettings(
        # Snakemake caps every rule's `threads:` to min(rule.threads,
        # cores) before rendering `{threads}` into the --cluster submit
        # template -- regardless of execution mode. A too-low cores= here
        # silently caps every job's requested ncpus even though -q's value
        # is irrelevant once --cluster is set, so use a sentinel well above
        # any rule's declared threads instead.
        cores=999 if args.cluster else args.cores,
        nodes=args.jobs if args.cluster else None,
        resources=resources or {},
        default_resources=(
            DefaultResources(
                [
                    f"mem_mb={CLUSTER_DEFAULT_MEM_MB}",
                    f"runtime={CLUSTER_DEFAULT_RUNTIME_MIN}",
                ]
            )
            if args.cluster
            else None
        ),
    )
    config_settings = ConfigSettings(configfiles=[configfile])
    deployment_settings = DeploymentSettings(
        deployment_method=build_deployment_method(args),
        apptainer_args=singularity_args,
        conda_prefix=Path(args.conda_prefix).resolve() if args.conda_prefix else None,
    )
    workflow_settings = WorkflowSettings(wrapper_prefix=args.wrapper_prefix)
    output_settings = OutputSettings(printshellcmds=args.printshellcmds)
    execution_settings = ExecutionSettings(
        latency_wait=args.latency_wait,
        keep_going=args.keep_going,
        lock=False,  # see acquire_project_lock()
        retries=args.restart_times,
    )
    rerun_triggers = rerun_triggers_from_names(args.rerun_triggers)
    dag_settings = DAGSettings(
        targets=set(args.workflow or []),
        forceall=args.forceall,
        force_incomplete=args.rerun_incomplete,
        **({"rerun_triggers": rerun_triggers} if rerun_triggers else {}),
    )

    try:
        with SnakemakeApi(output_settings) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                resource_settings=resource_settings,
                config_settings=config_settings,
                deployment_settings=deployment_settings,
                workflow_settings=workflow_settings,
                snakefile=snakefile,
                workdir=workdir,
            )

            if args.lint is not None:
                workflow_api.lint(json=args.lint == "json")
                return 0

            dag_api = workflow_api.dag(dag_settings=dag_settings)

            if args.unlock:
                dag_api.unlock()
                success("Working directory unlocked.")
                return 0

            executor_settings = (
                build_cluster_executor_settings(pipeline, args) if args.cluster else None
            )
            dag_api.execute_workflow(
                executor="dryrun" if args.dry_run else ("cluster-generic" if args.cluster else "local"),
                execution_settings=execution_settings,
                executor_settings=executor_settings,
            )
    except SystemExit:
        raise
    except Exception as exc:
        die(f"Workflow failed: {exc}")

    success("Workflow finished successfully.")
    return 0


def build_parser(pipeline: Pipeline) -> argparse.ArgumentParser:
    first = next(iter(pipeline.workflows), "all")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=color(f"Launcher for sm-core/{pipeline.name}: {pipeline.description}", Fore.GREEN),
        epilog=f"""
Examples:

  {Fore.CYAN}./run.py -w {first} -c config_MyProject.yaml -q 16 -n{Style.RESET_ALL}
      Dry-run: show what would run, without running it.

  {Fore.CYAN}./run.py -w {first} -c config_MyProject.yaml -q 16{Style.RESET_ALL}
      Run it locally on 16 cores.

  {Fore.CYAN}./run.py -w {first} -c config_MyProject.yaml -cl -qu workq -j 20{Style.RESET_ALL}
      Submit each job as its own PBS job via qsub instead of running
      everything locally in one allocation (--scheduler slurm: sbatch).

  {Fore.CYAN}./run.py --list-workflows{Style.RESET_ALL}
      Show available workflow targets.

  {Fore.CYAN}./run.py -w {first} -c config_MyProject.yaml -q 4 --dag | dot -Tsvg > dag.svg{Style.RESET_ALL}
      Export the DAG as SVG.
""",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"sm-core/{pipeline.name} {pipeline.version} (sm-core {SMCORE_VERSION})",
    )

    required = parser.add_argument_group("Required for normal execution")
    required.add_argument(
        "-w",
        "--workflow",
        nargs="+",
        help=f"Workflow target(s) to run. Main targets: {', '.join(pipeline.workflows.keys())}.",
    )
    required.add_argument(
        "-c",
        "--configfile",
        help="Path to the YAML/JSON Snakemake config file.",
    )
    required.add_argument(
        "-q",
        "--cores",
        type=int,
        default=1,
        help="Number of CPU cores available to Snakemake.",
    )

    execution = parser.add_argument_group("Execution mode")
    execution.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Perform a dry-run without executing jobs.",
    )
    execution.add_argument(
        "-f",
        "--forceall",
        action="store_true",
        help="Force execution of all jobs even if outputs already exist.",
    )
    execution.add_argument(
        "-ri",
        "--rerun-incomplete",
        action="store_true",
        help="Rerun jobs with incomplete outputs.",
    )
    execution.add_argument(
        "-u",
        "--unlock",
        action="store_true",
        help="Unlock the working directory after an interrupted run.",
    )
    execution.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue independent jobs after one job fails.",
    )
    execution.add_argument(
        "--latency-wait",
        type=int,
        default=60,
        help="Seconds to wait for output files to appear on slow filesystems.",
    )

    cluster = parser.add_argument_group("Cluster (PBS or SLURM)")
    cluster.add_argument(
        "-cl",
        "--cluster",
        action="store_true",
        help="Submit each rule instance as its own cluster job, instead of "
        "running everything locally in one allocation (default off; plain -q N "
        "still runs everything locally with N cores, unchanged).",
    )
    cluster.add_argument(
        "--scheduler",
        choices=["pbs", "slurm"],
        default="pbs",
        help="How --cluster submits: qsub (OpenPBS/PBS Pro) or sbatch (SLURM). "
        "On a SLURM cluster whose `qsub` is SLURM's PBS wrapper, use slurm.",
    )
    cluster.add_argument(
        "-qu",
        "--queue",
        type=str,
        default="workq",
        help="PBS queue, or SLURM partition, to submit to when --cluster is set.",
    )
    cluster.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=20,
        help="Max number of concurrent cluster jobs when --cluster is set.",
    )
    cluster.add_argument(
        "-rt",
        "--restart-times",
        type=int,
        default=2,
        help="Number of times to automatically retry a job that fails (e.g. "
        "killed for exceeding its resources), each retry with progressively "
        "more mem_mb/runtime for rules whose resources: use the scaled() "
        "helper (i.e. up to N+1 total attempts).",
    )

    reporting = parser.add_argument_group("Reporting and debugging")
    reporting.add_argument(
        "--list-workflows",
        "--show-workflows",
        dest="show_workflows",
        action="store_true",
        help="Print available high-level workflow targets and exit.",
    )
    reporting.add_argument(
        "-d",
        "--dag",
        action="store_true",
        help=f"Print DAG in DOT format. Example: ./run.py -w {first} ... --dag | dot -Tsvg > dag.svg",
    )
    reporting.add_argument(
        "--printshellcmds",
        "-p",
        action="store_true",
        help="Print shell commands before executing them.",
    )

    reporting.add_argument(
        "-l",
        "--lint",
        default=None,
        help='Run Snakemake lint. Common values: "text" or "json".',
    )

    environment = parser.add_argument_group("Environment and paths")
    environment.add_argument(
        "-di",
        "--directory",
        type=str,
        default=None,
        help="Working directory for Snakemake execution. Default: the pipeline's "
        "own directory, so relative paths in a config resolve the same way "
        "wherever run.py is launched from.",
    )
    environment.add_argument(
        "--snakefile",
        type=str,
        default=None,
        help="Custom Snakefile path. Default: the pipeline's own Snakefile.",
    )
    environment.add_argument(
        "--no-conda",
        action="store_true",
        help="Disable Snakemake conda integration.",
    )
    environment.add_argument(
        "--conda-prefix",
        type=str,
        default=None,
        help="Directory to build the per-rule conda environments in (default: "
        ".snakemake/conda inside the working directory). Point several "
        "checkouts at one shared prefix to build each environment once.",
    )
    environment.add_argument(
        "--no-singularity",
        action="store_true",
        help="Disable Snakemake Singularity/Apptainer integration.",
    )
    environment.add_argument(
        "--bind",
        action="append",
        default=[],
        help="Singularity bind path, e.g. --bind /data. Can be used multiple times.",
    )
    environment.add_argument(
        "--singularity-args",
        type=str,
        default="",
        help="Additional raw Singularity arguments.",
    )
    environment.add_argument(
        "-wr",
        "--wrapper-prefix",
        type=str,
        default=None,
        help="Prefix for a local Snakemake wrapper repository.",
    )

    advanced = parser.add_argument_group("Advanced")
    advanced.add_argument(
        "--resources",
        nargs="+",
        default=None,
        help="Custom Snakemake resources as key=value pairs, e.g. --resources mem_mb=64000.",
    )
    advanced.add_argument(
        "--allow-custom-target",
        action="store_true",
        help="Allow running targets not listed among the main workflow names.",
    )
    advanced.add_argument(
        "--rerun-triggers",
        nargs="+",
        default=None,
        choices=["mtime", "params", "input", "software-env", "code"],
        help="Override Snakemake's rerun-trigger set (default when omitted: mtime "
        "params input software-env code -- Snakemake's own default, reruns a job "
        "if ANY of those changed since last time, even if its output already "
        "exists). Pass --rerun-triggers mtime to fall back to classic make-style "
        "mtime-only comparison -- useful when adding new samples to an existing "
        "project config and an unrelated historical code/params edit elsewhere in "
        "the pipeline would otherwise force recomputing already-finished samples. "
        "Only use mtime-only when you're sure existing outputs are still valid "
        "under the current rule code: it will NOT detect a genuine bug fix to a "
        "rule's shell command/script as a reason to rerun, so a real fix could "
        "silently leave old, wrong outputs in place.",
    )

    return parser


def main(pipeline: Pipeline, argv: Optional[List[str]] = None) -> int:
    parser = build_parser(pipeline)
    args = parser.parse_args(argv)

    if args.show_workflows:
        print_workflows(pipeline)
        return 0

    if not args.configfile:
        die("Missing config file. Use -c config_MyProject.yaml")

    return run_snakemake(pipeline, args)
