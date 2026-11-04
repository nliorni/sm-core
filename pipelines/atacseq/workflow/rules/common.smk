###############################################################################
# atacseq: pipeline-specific helpers
#
# Included right after smcore/common.smk (outputDir, scaled(), env(),
# load_samples()): the validated config, the sample table, output and
# reference locations, and the sample groups used for consensus peaks and
# averaged signal tracks.
###############################################################################

import re
import sys

from snakemake.utils import validate

# Fills every `default:` from the schema into `config`.
validate(config, "../schemas/config.schema.yaml")

SAMPLE_TABLE = load_samples()
SAMPLES = list(SAMPLE_TABLE)

RESULTS = outputDir + "results/"
BAM_DIR = RESULTS + "bam/"
PEAKS_DIR = RESULTS + "peaks/"
COUNTS_DIR = RESULTS + "counts/"
BIGWIG_DIR = RESULTS + "bigwig/"
QC_DIR = RESULTS + "qc/"
LOGS = outputDir + "logs/"
BENCH = outputDir + "benchmarks/"


## ============================================================================
## Reference locations
##
## Built once under reference.dir (default <output_dir>/reference/) by
## workflow/rules/reference.smk; point several projects on the same genome +
## GTF at one shared directory to build the bowtie2 index once.
## ============================================================================

REF_DIR = config["reference"]["dir"] or (outputDir + "reference")
if not REF_DIR.endswith("/"):
    REF_DIR += "/"

GENOME_FA = REF_DIR + "genome.fa"
GENOME_FAI = REF_DIR + "genome.fa.fai"
CHROM_SIZES = REF_DIR + "chrom.sizes"
# The genome's sequences minus the mitochondrion, as BED: what the final BAMs
# keep (samtools view -L).
NUCLEAR_BED = REF_DIR + "nuclear.bed"
MITO_NAME = REF_DIR + "mito_name.txt"
GENOME_SIZE = REF_DIR + "effective_genome_size.txt"
GENES_TSV = REF_DIR + "genes.tsv"
TSS_BED = REF_DIR + "tss.bed"
BLACKLIST_BED = REF_DIR + "blacklist.bed"
HAS_BLACKLIST = bool(config["reference"]["blacklist"])

BOWTIE2_SUFFIXES = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")


def bowtie2_prefix():
    """A prebuilt index from reference.bowtie2_index, or the one built here."""
    return config["reference"]["bowtie2_index"] or (REF_DIR + "bowtie2/genome")


def bowtie2_index_files():
    prefix = bowtie2_prefix()
    if config["reference"]["bowtie2_index"]:
        # Large genomes (> 4 Gbp) get .bt2l files; accept either.
        large = Path(prefix + ".1.bt2l").exists()
        return [prefix + s + ("l" if large else "") for s in BOWTIE2_SUFFIXES]
    return [prefix + s for s in BOWTIE2_SUFFIXES]


## ============================================================================
## Sample groups
##
## peaks.consensus_group_by and bigwig.group_by name samplesheet columns; a
## group is one combination of their values (e.g. cell_line + condition), and
## its name joins them with "_".
## ============================================================================

def _meta_columns():
    return sorted({k for s in SAMPLES for k in SAMPLE_TABLE[s]["meta"]})


def sample_groups(columns, key):
    """{group name: [samples]} for a list of samplesheet columns."""
    if not columns:
        return {"all": list(SAMPLES)}
    known = _meta_columns()
    for col in columns:
        if col not in known:
            config_error(f"{key}: '{col}' is not a samplesheet column (have: {', '.join(known) or 'none'}).")
    groups = {}
    for s in SAMPLES:
        name = "_".join(SAMPLE_TABLE[s]["meta"][c] for c in columns)
        if not re.fullmatch(r"[A-Za-z0-9_.+-]+", name):
            config_error(f"{key}: group '{name}' (sample {s}) may only contain letters, digits and _ . + -.")
        groups.setdefault(name, []).append(s)
    return groups


CONSENSUS_GROUPS = sample_groups(config["peaks"]["consensus_group_by"], "peaks.consensus_group_by")
# With a p-value cut-off the replicate requirement is what filters noise; a
# group with fewer samples than min_replicates keeps peaks seen once.
_unreplicated = [g for g, members in CONSENSUS_GROUPS.items()
                 if min(config["peaks"]["min_replicates"], len(members)) < 2]
if config["peaks"]["cutoff"] == "pvalue" and _unreplicated:
    print(f"WARNING: consensus group(s) {', '.join(_unreplicated)} keep peaks called in a single sample, "
          f"at a permissive p < {config['peaks']['threshold']}; consider peaks.cutoff: qvalue, "
          "threshold: 0.05 (docs/configuration.md#peaks).", file=sys.stderr)
BIGWIG_GROUPS = (
    sample_groups(config["bigwig"]["group_by"], "bigwig.group_by") if config["bigwig"]["group_by"] else {}
)


## ============================================================================
## Identity check (somalier) -- optional
## ============================================================================

IDENTITY = bool(config["qc"]["somalier_sites"])
IDENTITY_BAMS = list(config["qc"]["identity_bams"])
if IDENTITY_BAMS and not IDENTITY:
    config_error("qc.identity_bams needs qc.somalier_sites (the sites the identity check genotypes).")
for _bam in IDENTITY_BAMS:
    if not Path(_bam).exists():
        config_error(f"qc.identity_bams: {_bam} does not exist.")
