###############################################################################
# Differential expression -- PyDESeq2, then ORA and GSEA
#
# For every analysis in config `de.analyses`, every value of its split_by
# column (or all samples together), and every contrast:
#
#   de_contrast     PyDESeq2 Wald test, shrunken log2 fold changes, the
#                   significant up/down lists, volcano / MA / heatmap
#   de_overview     the model's samples: VST, PCA, sample distances
#   de_enrichment   ORA of the up and down lists against the genes tested,
#                   and preranked GSEA on all of them
#   de_summary      DEG counts for every analysis, split and contrast
#
# A split is a separate model per value of a samplesheet column (one per
# cell line, say): dispersions are then estimated within each value, which is
# what to do when the groups differ a lot in variability. Without split_by,
# all samples share one model and the design can adjust for that column
# instead (e.g. "~ cell_line + condition").
#
# Results: results/de/<matrix>/<analysis>/<split>/<contrast>/
###############################################################################

DE = config["de"]
DE_DIR = RESULTS + "de/"
DE_FORMATS = DE["plot_formats"]
DE_ENV = env("de")

MSIGDB_FILES = {
    # collection: (human file prefix, mouse file prefix)
    "hallmark": ("h.all", "mh.all"),
    "reactome": ("c2.cp.reactome", "m2.cp.reactome"),
    "go_bp": ("c5.go.bp", "m5.go.bp"),
}
_GS = DE["gene_sets"]
GENE_SET_DIR = REF_DIR + f"gene_sets/msigdb_{_GS['msigdb_version']}_{_GS['species']}/"


def gene_set_files():
    files = {c: GENE_SET_DIR + f"{c}.gmt" for c in _GS["collections"]}
    for name, path in _GS["custom"].items():
        if name in files:
            config_error(f"de.gene_sets.custom: '{name}' is also an MSigDB collection name.")
        path = Path(path) if Path(path).is_absolute() else PIPELINE_DIR / path
        if not path.exists():
            config_error(f"de.gene_sets.custom.{name}: {path} does not exist.")
        files[name] = str(path)
    return files


GENE_SET_FILES = gene_set_files()


# The analyses, checked against the samplesheet (smcore/rules/de_design.smk).
DE_UNITS = de_units(DE["analyses"], "de")


def de_split_dirs():
    return sorted({(a, s) for a, s, _ in DE_UNITS})


def figures(prefix):
    return [f"{prefix}.{fmt}" for fmt in DE_FORMATS]


def de_counts(wildcards):
    if wildcards.branch == "salmon":
        return COUNTS_DIR + "salmon/gene_counts_length_scaled.tsv"
    return COUNTS_DIR + "star/gene_counts.tsv"


def de_outputs(wildcards=None):
    if not DE_UNITS:
        raise WorkflowError("The `de` target needs at least one entry in de.analyses (see docs/configuration.md).")
    out = []
    for branch in DE["matrices"]:
        base = DE_DIR + branch + "/"
        out += figures(base + "summary") + [base + "summary.tsv"]
        for analysis, split in de_split_dirs():
            out += figures(f"{base}{analysis}/{split}/pca") + [f"{base}{analysis}/{split}/vst.tsv"]
        for analysis, split, (_, test, reference) in DE_UNITS:
            d = f"{base}{analysis}/{split}/{test}_vs_{reference}/"
            out += [d + "deseq2_results.tsv", d + "enrichment/gsea.tsv"]
    return out


wildcard_constraints:
    analysis="[A-Za-z0-9_.-]+",
    split="[A-Za-z0-9_.+-]+",
    contrast="[A-Za-z0-9_.+-]+_vs_[A-Za-z0-9_.+-]+",


def de_params(wildcards):
    return de_model_params(DE["analyses"], DE_UNITS, wildcards)


localrules: de_gene_sets


rule de_gene_sets:
    # MSigDB gene sets (gene symbols), downloaded once next to the reference
    # so every project on that reference shares the same pinned release.
    output:
        GENE_SET_DIR + "{collection}.gmt",
    wildcard_constraints:
        collection="|".join(MSIGDB_FILES),
    params:
        version=_GS["msigdb_version"],
        species=_GS["species"],
    run:
        import urllib.request

        human, mouse = MSIGDB_FILES[wildcards.collection]
        tag = "Hs" if params.species == "human" else "Mm"
        prefix = human if params.species == "human" else mouse
        url = (
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
            f"{params.version}.{tag}/{prefix}.v{params.version}.{tag}.symbols.gmt"
        )
        tmp = Path(output[0] + ".part")
        with urllib.request.urlopen(url, timeout=300) as response:
            tmp.write_bytes(response.read())
        lines = tmp.read_text().splitlines()
        if not lines or any(len(line.split("\t")) < 3 for line in lines):
            tmp.unlink()
            raise WorkflowError(f"{url} did not return a GMT file.")
        tmp.rename(output[0])


rule de_contrast:
    input:
        counts=de_counts,
        samples=COUNTS_DIR + "{branch}/samples.tsv",
        tx2gene=TX2GENE,
    output:
        results=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/deseq2_results.tsv",
        up=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/significant_up.tsv",
        down=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/significant_down.tsv",
        plots=figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/volcano")
        + figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/ma")
        + figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/heatmap"),
    log:
        LOGS + "de/{branch}/{analysis}/{split}/{contrast}.log",
    benchmark:
        BENCH + "de/{branch}/{analysis}/{split}/{contrast}.tsv"
    conda:
        DE_ENV
    threads: DE["threads"]
    params:
        model=de_params,
        de=DE,
        smcore=SMCORE_ROOT,
    resources:
        mem_mb=scaled(DE["mem_mb"]),
        runtime=scaled(60),
    script:
        "../scripts/de_contrast.py"


rule de_overview:
    input:
        counts=de_counts,
        samples=COUNTS_DIR + "{branch}/samples.tsv",
        tx2gene=TX2GENE,
    output:
        vst=DE_DIR + "{branch}/{analysis}/{split}/vst.tsv",
        normalized=DE_DIR + "{branch}/{analysis}/{split}/normalized_counts.tsv",
        size_factors=DE_DIR + "{branch}/{analysis}/{split}/size_factors.tsv",
        plots=figures(DE_DIR + "{branch}/{analysis}/{split}/pca")
        + figures(DE_DIR + "{branch}/{analysis}/{split}/sample_distances"),
    log:
        LOGS + "de/{branch}/{analysis}/{split}/overview.log",
    conda:
        DE_ENV
    threads: DE["threads"]
    params:
        model=de_params,
        de=DE,
        smcore=SMCORE_ROOT,
    resources:
        mem_mb=scaled(DE["mem_mb"]),
        runtime=scaled(60),
    script:
        "../scripts/de_overview.py"


rule de_enrichment:
    input:
        results=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/deseq2_results.tsv",
        gene_sets=list(GENE_SET_FILES.values()),
    output:
        ora_up=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/ora_up.tsv",
        ora_down=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/ora_down.tsv",
        gsea=DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/gsea.tsv",
        plots=figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/ora")
        + figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/gsea")
        + figures(DE_DIR + "{branch}/{analysis}/{split}/{contrast}/enrichment/gsea_curves"),
    log:
        LOGS + "de/{branch}/{analysis}/{split}/{contrast}.enrichment.log",
    benchmark:
        BENCH + "de/{branch}/{analysis}/{split}/{contrast}.enrichment.tsv"
    conda:
        DE_ENV
    threads: DE["threads"]
    params:
        collections=list(GENE_SET_FILES),
        de=DE,
        smcore=SMCORE_ROOT,
    resources:
        mem_mb=scaled(DE["mem_mb"]),
        runtime=scaled(120),
    script:
        "../scripts/de_enrichment.py"


def de_summary_inputs(wildcards):
    out = []
    for analysis, split, (_, test, reference) in DE_UNITS:
        d = f"{DE_DIR}{wildcards.branch}/{analysis}/{split}/{test}_vs_{reference}/"
        out += [d + "deseq2_results.tsv", d + "enrichment/gsea.tsv"]
    return out


rule de_summary:
    input:
        de_summary_inputs,
    output:
        tsv=DE_DIR + "{branch}/summary.tsv",
        plots=figures(DE_DIR + "{branch}/summary"),
    conda:
        DE_ENV
    params:
        units=[(a, s, f"{c[1]}_vs_{c[2]}") for a, s, c in DE_UNITS],
        de=DE,
        de_dir=DE_DIR,
        smcore=SMCORE_ROOT,
    script:
        "../scripts/de_summary.py"
