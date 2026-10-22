###############################################################################
# rnaseq: pipeline-specific helpers
#
# Included right after smcore/common.smk (outputDir, scaled(), env(),
# load_samples()), so everything defined here is available to every rule
# module: the validated config, the sample table, output locations, and the
# small input-picking functions that let one rule read from different
# upstream files depending on the config.
###############################################################################

import sys

from snakemake.utils import validate

# Fills every `default:` from the schema into `config`, so rules can read
# config["star"]["threads"] directly instead of nesting .get() calls.
validate(config, "../schemas/config.schema.yaml")

SAMPLE_TABLE = load_samples()
SAMPLES = list(SAMPLE_TABLE)

RESULTS = outputDir + "results/"
COUNTS_DIR = RESULTS + "counts/"
QC_DIR = RESULTS + "qc/"
LOGS = outputDir + "logs/"
BENCH = outputDir + "benchmarks/"


## ============================================================================
## Branches
##
## Two independent quantification branches start from the same (trimmed)
## reads and never read each other's output. Each is its own workflow target:
##
##   -w star     STAR genome alignment -> BAM + STAR's per-gene read counts
##   -w salmon   Salmon selective alignment (decoy-aware) -> tximport matrices
##   -w all      both, plus a comparison of their strandedness calls
##
## Each infers strandedness from its own data, and builds only its own index.
## ============================================================================

BRANCHES = ("star", "salmon")


def trimming_enabled():
    return config["trimming"]["enabled"]


## ============================================================================
## Reference locations
##
## Everything derived from the configured genome + GTF is built once under
## reference.dir (default <output_dir>/reference/) by workflow/rules/reference.smk.
## Point reference.dir at a shared location to reuse one build across
## projects that use the same genome and annotation.
## ============================================================================

REF_DIR = config["reference"]["dir"] or (outputDir + "reference")
if not REF_DIR.endswith("/"):
    REF_DIR += "/"

GENOME_FA = REF_DIR + "genome.fa"
GENOME_FAI = REF_DIR + "genome.fa.fai"
ANNOTATION_GTF = REF_DIR + "genes.gtf"
TRANSCRIPTS_FA = REF_DIR + "transcripts.fa"
TX2GENE = REF_DIR + "tx2gene.tsv"


COLLAPSED_GTF = REF_DIR + "genes.collapsed.gtf"


def star_index():
    """A prebuilt index from reference.star_index, or the one built here.

    A prebuilt index must come from exactly the genome and GTF configured
    above: STAR's gene counts are defined by the annotation the index was
    built with, and tx2gene.tsv (gene names, biotypes) by the configured GTF."""
    return config["reference"]["star_index"] or (REF_DIR + "star")


def salmon_index():
    """Same caveat as star_index(): Salmon's transcript IDs must be those of
    the configured GTF, or tximport cannot map them to genes."""
    return config["reference"]["salmon_index"] or (REF_DIR + "salmon")


## ============================================================================
## Reads
## ============================================================================

def raw_reads(wildcards):
    """[R1, R2] for a sample: its single FASTQ pair as given, or the lanes /
    runs concatenated into one pair by merge_lanes when there are several."""
    units = SAMPLE_TABLE[wildcards.sample]["units"]
    if len(units) == 1:
        return list(units[0])
    return [
        RESULTS + f"fastq/{wildcards.sample}_R1.merged.fastq.gz",
        RESULTS + f"fastq/{wildcards.sample}_R2.merged.fastq.gz",
    ]


def analysis_reads(wildcards):
    """What STAR and Salmon read: fastp's trimmed reads, or the raw reads with
    trimming.enabled: false."""
    if trimming_enabled():
        return [
            RESULTS + f"fastp/{wildcards.sample}_R1.trimmed.fastq.gz",
            RESULTS + f"fastp/{wildcards.sample}_R2.trimmed.fastq.gz",
        ]
    return raw_reads(wildcards)


def decompress_cmd(path):
    """STAR's --readFilesCommand for a read file."""
    return "zcat" if str(path).endswith(".gz") else "-"


## ============================================================================
## Strandedness
## ============================================================================

# Salmon library types for paired-end reads: I = inward-facing mates,
# S = stranded, F/R = read 1 comes from the forward/reverse strand of the
# transcript, U = unstranded. A = let Salmon detect it.
SALMON_LIBTYPE = {"auto": "A", "unstranded": "IU", "forward": "ISF", "reverse": "ISR"}


def salmon_libtype(wildcards):
    return SALMON_LIBTYPE[SAMPLE_TABLE[wildcards.sample]["strandedness"]]
