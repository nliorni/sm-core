###############################################################################
# Differential accessibility -- PyDESeq2 on the consensus peak counts
#
# For every analysis in config `da.analyses`, every value of its split_by
# column (or all samples together), and every contrast:
#
#   da_overview   the model's samples: size factors of both normalisations,
#                 VST, PCA, sample distances
#   da_contrast   Wald test, shrunken log2 fold changes, the significant
#                 peaks (TSV and BED), a per-gene summary, volcano / MA /
#                 heatmap, and where the significant peaks lie (promoter /
#                 genic / distal)
#   da_summary    up / down peak counts of every analysis, split and contrast
#
# The models, statistics and figures are shared with rnaseq's differential
# expression (smcore/de.py, smcore/rules/de_design.smk); this file adds what
# is specific to peaks: the normalisation choice, the annotation, BED output.
#
# Results: results/da/<analysis>/<split>/<contrast>/
###############################################################################

DA = config["da"]
DA_DIR = RESULTS + "da/"
DA_FORMATS = DA["plot_formats"]
DA_ENV = env("de")

DA_UNITS = de_units(DA["analyses"], "da")


def da_figures(prefix):
    return [f"{prefix}.{fmt}" for fmt in DA_FORMATS]


def da_outputs(wildcards=None):
    if not DA_UNITS:
        raise WorkflowError("The `da` target needs at least one entry in da.analyses (see docs/configuration.md).")
    out = da_figures(DA_DIR + "summary") + [DA_DIR + "summary.tsv"]
    for analysis, split in sorted({(a, s) for a, s, _ in DA_UNITS}):
        d = f"{DA_DIR}{analysis}/{split}/"
        out += da_figures(d + "pca") + [d + "vst.tsv", d + "size_factors.tsv"]
    for analysis, split, (_, test, reference) in DA_UNITS:
        out.append(f"{DA_DIR}{analysis}/{split}/{test}_vs_{reference}/da_results.tsv")
    return out


def da_params(wildcards):
    return de_model_params(DA["analyses"], DA_UNITS, wildcards)


wildcard_constraints:
    analysis="[A-Za-z0-9_.-]+",
    split="[A-Za-z0-9_.+-]+",
    contrast="[A-Za-z0-9_.+-]+_vs_[A-Za-z0-9_.+-]+",


DA_INPUTS = dict(
    counts=COUNTS_DIR + "peak_counts.tsv",
    samples=COUNTS_DIR + "samples.tsv",
    # every fragment of each final BAM, counted or not: the library sizes
    fc_summary=COUNTS_DIR + "featurecounts.txt.summary",
    annotation=CONSENSUS_DIR + "consensus_peaks.annotated.tsv",
    genes=GENES_TSV,
)


rule da_overview:
    input:
        **DA_INPUTS,
    output:
        size_factors=DA_DIR + "{analysis}/{split}/size_factors.tsv",
        vst=DA_DIR + "{analysis}/{split}/vst.tsv",
        normalized=DA_DIR + "{analysis}/{split}/normalized_counts.tsv",
        plots=da_figures(DA_DIR + "{analysis}/{split}/pca")
        + da_figures(DA_DIR + "{analysis}/{split}/sample_distances")
        + da_figures(DA_DIR + "{analysis}/{split}/size_factors"),
    log:
        LOGS + "da/{analysis}/{split}/overview.log",
    conda:
        DA_ENV
    threads: DA["threads"]
    params:
        model=da_params,
        da=DA,
        bams={s: BAM_DIR + f"{s}.bam" for s in SAMPLES},
        smcore=SMCORE_ROOT,
    resources:
        mem_mb=scaled(DA["mem_mb"]),
        runtime=scaled(60),
    script:
        "../scripts/da_overview.py"


rule da_contrast:
    input:
        **DA_INPUTS,
    output:
        results=DA_DIR + "{analysis}/{split}/{contrast}/da_results.tsv",
        up=DA_DIR + "{analysis}/{split}/{contrast}/significant_up.tsv",
        down=DA_DIR + "{analysis}/{split}/{contrast}/significant_down.tsv",
        up_bed=DA_DIR + "{analysis}/{split}/{contrast}/significant_up.bed",
        down_bed=DA_DIR + "{analysis}/{split}/{contrast}/significant_down.bed",
        genes=DA_DIR + "{analysis}/{split}/{contrast}/genes.tsv",
        plots=da_figures(DA_DIR + "{analysis}/{split}/{contrast}/volcano")
        + da_figures(DA_DIR + "{analysis}/{split}/{contrast}/ma")
        + da_figures(DA_DIR + "{analysis}/{split}/{contrast}/heatmap")
        + da_figures(DA_DIR + "{analysis}/{split}/{contrast}/peak_classes"),
    log:
        LOGS + "da/{analysis}/{split}/{contrast}.log",
    benchmark:
        BENCH + "da/{analysis}/{split}/{contrast}.tsv"
    conda:
        DA_ENV
    threads: DA["threads"]
    params:
        model=da_params,
        da=DA,
        bams={s: BAM_DIR + f"{s}.bam" for s in SAMPLES},
        smcore=SMCORE_ROOT,
    resources:
        mem_mb=scaled(DA["mem_mb"]),
        runtime=scaled(120),
    script:
        "../scripts/da_contrast.py"


rule da_summary:
    input:
        [f"{DA_DIR}{a}/{s}/{c[1]}_vs_{c[2]}/da_results.tsv" for a, s, c in DA_UNITS],
    output:
        tsv=DA_DIR + "summary.tsv",
        plots=da_figures(DA_DIR + "summary"),
    conda:
        DA_ENV
    params:
        units=[(a, s, f"{c[1]}_vs_{c[2]}") for a, s, c in DA_UNITS],
        da=DA,
        da_dir=DA_DIR,
        smcore=SMCORE_ROOT,
    script:
        "../scripts/da_summary.py"
