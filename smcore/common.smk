###############################################################################
# sm-core shared Snakemake helpers
#
# Included first by every pipeline's Snakefile:
#
#     include: "../../smcore/common.smk"
#
# What lives here is what every pipeline of the suite does the same way:
# the output directory, retry-scaled resources, per-rule software
# environments with site overrides, the paired-end FASTQ input convention
# (config `samples:` or an nf-core-style `samplesheet:`), and completion
# e-mails. Anything specific to one pipeline belongs in that pipeline's own
# workflow/rules/common.smk instead.
###############################################################################

import csv
from pathlib import Path

from snakemake.exceptions import WorkflowError


PIPELINE_DIR = Path(workflow.basedir)
# The sm-core checkout itself: scripts running in a rule's conda environment
# put it on sys.path to import shared code (smcore.plotting), which needs
# nothing beyond matplotlib.
SMCORE_ROOT = str(PIPELINE_DIR.parents[1])
SCRIPTS = str(PIPELINE_DIR / "workflow" / "scripts")
ENVS = PIPELINE_DIR / "workflow" / "envs"


def config_error(message):
    """Stop at parse time with a message naming the config key at fault,
    rather than deep inside a rule with a tool's stack trace."""
    raise WorkflowError(f"config: {message}")


if "output_dir" not in config:
    config_error("'output_dir' is required.")
outputDir = str(config["output_dir"])
if not outputDir.endswith("/"):
    outputDir += "/"

# A sample ID is one path component; without this, a sample whose ID is a
# prefix of another's could match the other's files.
wildcard_constraints:
    sample="[^/]+",


## ============================================================================
## Retry-scaled resources.
## ============================================================================

def scaled(base):
    """Make a resource a function of Snakemake's retry counter.

    A rule declaring `mem_mb=scaled(8000)` asks for 8000 MB on its first
    attempt, 16000 on the first retry, 24000 on the next. Combined with
    run.py's -rt/--restart-times (default 2), a job killed by the scheduler
    for exceeding its declared memory or walltime comes back asking for more
    instead of failing the run.
    """
    def _scaled(wildcards, attempt):
        return int(base) * attempt
    return _scaled


## ============================================================================
## Software environments.
## ============================================================================

def env(name):
    """Resolve a rule's software environment.

    By default this is the conda definition bundled in the pipeline's
    workflow/envs/<name>.yaml, which Snakemake builds on first use. Sites that
    already maintain shared, pre-built environments skip that build by naming
    them in the project config:

        conda_envs:
          star: /shared/condaEnvs/star

    A value can be either an environment prefix (a directory) or a path to
    another environment YAML. Names not listed fall back to the bundled file.
    """
    override = (config.get("conda_envs") or {}).get(name)
    return override if override else str(ENVS / f"{name}.yaml")


## ============================================================================
## Paired-end FASTQ inputs.
##
## Two ways to declare samples, resolved here into one table so no rule ever
## re-derives a path:
##
## 1. `samples:` + `unit:` + `extension:` -- the sm-core family convention.
##    Each value is a FASTQ *prefix*; the reads are
##    {prefix}_{unit[0]}.{extension} and {prefix}_{unit[1]}.{extension}. A
##    value may also be a list of prefixes: the sample was sequenced over
##    several lanes or runs, and they are merged.
##
## 2. `samplesheet:` -- an nf-core-compatible CSV (or TSV) with columns
##    sample,fastq_1,fastq_2 and optionally strandedness. Rows sharing a
##    `sample` are lanes/runs of one library and are merged. Any further
##    column is sample metadata (condition, batch, ...) carried through to
##    the pipeline's outputs. Relative FASTQ paths are read relative to the
##    samplesheet itself, so a sheet stays valid when moved with its data.
## ============================================================================

STRANDEDNESS = ("auto", "unstranded", "forward", "reverse")
SAMPLESHEET_COLUMNS = ("sample", "fastq_1", "fastq_2", "strandedness")


def _check_sample_id(sample):
    if not sample or any(c in sample for c in "/ \t") or sample.startswith("."):
        config_error(
            f"sample ID '{sample}' must be a non-empty name without spaces or '/', "
            "since it becomes part of every output path."
        )


def _samples_from_prefixes():
    unit = config.get("unit")
    extension = config.get("extension")
    if not (isinstance(unit, list) and len(unit) == 2) or not extension:
        config_error(
            "with `samples:`, set `unit:` to the two read tokens (e.g. [\"R1\", \"R2\"]) "
            "and `extension:` (e.g. \"fastq.gz\")."
        )
    table = {}
    for sample, prefixes in config["samples"].items():
        sample = str(sample)
        _check_sample_id(sample)
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        if not prefixes:
            config_error(f"sample '{sample}' has no FASTQ prefix.")
        table[sample] = {
            "units": [
                (f"{p}_{unit[0]}.{extension}", f"{p}_{unit[1]}.{extension}")
                for p in prefixes
            ],
            "strandedness": config.get("strandedness", "auto"),
            "meta": {},
        }
    return table


def _samples_from_samplesheet(path):
    sheet = Path(path)
    if not sheet.exists():
        config_error(f"samplesheet '{path}' does not exist.")
    delimiter = "\t" if sheet.suffix in (".tsv", ".txt") else ","
    with sheet.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    if not rows or "sample" not in rows[0] or "fastq_1" not in rows[0]:
        config_error(f"samplesheet '{path}' needs at least the columns sample,fastq_1,fastq_2.")

    def resolve(p):
        p = p.strip()
        return p if Path(p).is_absolute() else str((sheet.parent / p).resolve())

    table = {}
    for n, row in enumerate(rows, start=2):
        sample = (row.get("sample") or "").strip()
        _check_sample_id(sample)
        r1, r2 = (row.get("fastq_1") or "").strip(), (row.get("fastq_2") or "").strip()
        if not r1 or not r2:
            config_error(
                f"samplesheet line {n} ({sample}): both fastq_1 and fastq_2 are required -- "
                "this pipeline is paired-end only."
            )
        strand = (row.get("strandedness") or config.get("strandedness", "auto")).strip()
        meta = {
            k: (v or "").strip() for k, v in row.items()
            if k not in SAMPLESHEET_COLUMNS and k is not None
        }
        entry = table.setdefault(sample, {"units": [], "strandedness": strand, "meta": meta})
        if entry["strandedness"] != strand:
            config_error(f"samplesheet: sample '{sample}' has rows with different strandedness.")
        if entry["meta"] != meta:
            config_error(f"samplesheet: sample '{sample}' has rows with different metadata.")
        entry["units"].append((resolve(r1), resolve(r2)))
    return table


def load_samples():
    """{sample: {"units": [(r1, r2), ...], "strandedness": str, "meta": {...}}}"""
    has_sheet, has_samples = bool(config.get("samplesheet")), bool(config.get("samples"))
    if has_sheet == has_samples:
        config_error("declare samples with exactly one of `samplesheet:` or `samples:`.")
    table = _samples_from_samplesheet(config["samplesheet"]) if has_sheet else _samples_from_prefixes()
    for sample, entry in table.items():
        if entry["strandedness"] not in STRANDEDNESS:
            config_error(
                f"sample '{sample}': strandedness '{entry['strandedness']}' must be one of "
                f"{', '.join(STRANDEDNESS)}."
            )
        missing = [f for pair in entry["units"] for f in pair if not Path(f).exists()]
        if missing:
            config_error(f"sample '{sample}': FASTQ not found: {', '.join(missing)}")
    return table


## ============================================================================
## Completion notification -- opt-in via `notify_email:` in the project config.
## Omit it and these silently no-op. Uses the system mail/mailx command with a
## 15s timeout so a broken mail relay can't hang the pipeline's completion.
## ============================================================================

def _notify_email(subject, body):
    email = config.get("notify_email")
    if not email:
        return
    import shutil
    import subprocess
    mail_bin = shutil.which("mail") or shutil.which("mailx")
    if not mail_bin:
        return
    try:
        subprocess.run([mail_bin, "-s", subject, email], input=body.encode(), check=False, timeout=15)
    except Exception:
        pass


onsuccess:
    _notify_email(
        f"[sm-core/{PIPELINE_DIR.name}] SUCCESS: {outputDir}",
        f"Pipeline run completed successfully.\nOutput directory: {outputDir}\n",
    )

onerror:
    _notify_email(
        f"[sm-core/{PIPELINE_DIR.name}] FAILED: {outputDir}",
        f"Pipeline run failed. Check the Snakemake log for details:\n{log}\n"
        f"Output directory: {outputDir}\n",
    )
