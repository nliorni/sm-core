outputDir = config["output_dir"]
if not outputDir.endswith("/"):
    outputDir += "/"

# Defaults to {} rather than a bare config["samples"] subscript so the
# Snakefile still parses for a config that legitimately has no `samples:`
# (e.g. while scaffolding a new project before wiring up real inputs).
SAMPLES = config.get("samples", {})


## ============================================================================
## Retry-scaled resources.
##
## Wrap a rule's mem_mb/runtime in scaled() so a job killed for exceeding its
## request (OOM or walltime) is automatically retried by Snakemake with
## progressively more -- 1x on the first attempt, 2x on the first retry, 3x
## on the second, etc. Retries are enabled via run.py's -rt/--restart-times
## (default 2, i.e. up to 3 total attempts).
## ============================================================================

def scaled(base):
    def _scaled(wildcards, attempt):
        return base * attempt
    return _scaled


## ============================================================================
## Shared paired-input resolution.
##
## `samples` in the config maps a sample name to a file *prefix*; the actual
## input is "{prefix}_{unit}{file_ext}", where `unit` is a project-specific
## pair of tokens (e.g. ["R1_001", "R2_001"] for FASTQ) and `file_ext` is
## whatever extension your real inputs use. This is the pattern every rule
## should use to resolve a sample's inputs, rather than re-deriving the path
## inline -- keeps the naming convention in one place. The demo below applies
## it to plain text files; swap in your real file type when you replace
## workflow/rules/example.smk.
## ============================================================================

def input_1(wildcards):
    ext = config.get("file_ext", ".txt")
    return f"{config['samples'][wildcards.sample]}_{config['unit'][0]}{ext}"

def input_2(wildcards):
    ext = config.get("file_ext", ".txt")
    return f"{config['samples'][wildcards.sample]}_{config['unit'][1]}{ext}"


include: "workflow/rules/example.smk"


rule all:
    input:
        outputDir + "summary.tsv"

# Alias so `./run.py -w example` (see run.py's WORKFLOWS dict) resolves to
# the same target as `all` -- keep every entry in WORKFLOWS pointing at a
# real rule/file target like this one.
rule example:
    input:
        outputDir + "summary.tsv"
