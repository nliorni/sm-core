"""Shared by the de_* scripts: reading a branch's counts and sample table,
selecting one model's samples, and fitting it with PyDESeq2."""

from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd


def log_to(path):
    """Send this script's output (and PyDESeq2's) to the rule's log."""
    handle = open(path, "w")
    sys.stdout = sys.stderr = handle
    return handle


def design_factors(design: str) -> list[str]:
    names = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", design.split("~", 1)[1])
    return [n for n in dict.fromkeys(names) if n not in ("C", "Q")]


def read_counts(counts_tsv: str, tx2gene_tsv: str | None = None):
    """Genes x samples integer counts (spike-ins removed) and a gene table
    (gene_name, gene_type, chrom)."""
    counts = pd.read_csv(counts_tsv, sep="\t", index_col="gene_id")
    names = counts.pop("gene_name")
    genes = pd.DataFrame({"gene_name": names})
    if tx2gene_tsv:
        tx = pd.read_csv(tx2gene_tsv, sep="\t").drop_duplicates("gene_id").set_index("gene_id")
        genes = genes.join(tx[["gene_type", "chrom"]])
        # ERCC spike-ins are controls, not genes: they would distort size
        # factors and dispersion trends.
        keep = genes["gene_type"].ne("spikein")
        counts, genes = counts[keep], genes[keep]
    # Salmon's length-scaled counts are estimates; DESeq2's model needs integers.
    return counts.round().astype(int), genes


def model_samples(samples_tsv: str, model: dict, split: str) -> pd.DataFrame:
    """The sample table of one model: its split, minus excluded samples, design
    factors as categorical (a numeric-looking column such as `replicate` must
    not become a continuous covariate)."""
    meta = pd.read_csv(samples_tsv, sep="\t", dtype=str, index_col="sample").fillna("")
    if model["split_by"]:
        meta = meta[meta[model["split_by"]] == split]
    meta = meta.drop(index=[s for s in model["exclude"] if s in meta.index])
    factors = design_factors(model["design"])
    return meta[factors].astype("category")


def with_reference(meta: pd.DataFrame, factor: str, reference: str) -> pd.DataFrame:
    """Put `reference` first among the factor's levels, so the model's
    coefficient for the test level is test vs reference (needed for shrinkage)."""
    meta = meta.copy()
    levels = [reference] + [l for l in sorted(meta[factor].unique()) if l != reference]
    meta[factor] = pd.Categorical(meta[factor], categories=levels)
    for col in meta.columns:
        if col != factor:
            meta[col] = meta[col].cat.remove_unused_categories()
    return meta


def prefilter(counts: pd.DataFrame, meta: pd.DataFrame, factor: str, min_count: int) -> pd.DataFrame:
    """Keep genes with >= min_count reads in at least as many samples as the
    smallest group of `factor` (DESeq2 vignette's pre-filter)."""
    smallest = int(meta[factor].value_counts().min())
    keep = (counts[meta.index] >= min_count).sum(axis=1) >= smallest
    return counts.loc[keep, meta.index]


def fit(counts: pd.DataFrame, meta: pd.DataFrame, design: str, threads: int):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    inference = DefaultInference(n_cpus=threads)
    dds = DeseqDataSet(
        counts=counts.T,  # PyDESeq2 wants samples x genes
        metadata=meta,
        design=design,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    return dds, inference


def vst(counts: pd.DataFrame, meta: pd.DataFrame, design: str, threads: int) -> pd.DataFrame:
    """Blind variance-stabilised expression (genes x samples) for QC plots."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    dds = DeseqDataSet(
        counts=counts.T, metadata=meta, design=design,
        inference=DefaultInference(n_cpus=threads), quiet=True,
    )
    dds.vst(use_design=False)
    return pd.DataFrame(dds.layers["vst_counts"].T, index=counts.index, columns=counts.columns), dds


def neg_log10(p: pd.Series) -> pd.Series:
    """-log10(p) with p = 0 (below double precision) set just past the smallest
    non-zero value, so those genes stay on the plot."""
    p = p.astype(float)
    floor = p[p > 0].min() if (p > 0).any() else 1e-300
    return -np.log10(p.clip(lower=floor / 10))


def short_term(term: str, width: int = 48) -> str:
    """'HALLMARK_TNFA_SIGNALING_VIA_NFKB' -> 'TNFA signaling via NFKB'."""
    t = re.sub(r"^(HALLMARK|REACTOME|GOBP|GOCC|GOMF|KEGG|WP)_", "", term)
    words = t.split("_")
    # Gene and complex names stay upper case; ordinary words do not.
    t = " ".join(CASED.get(w, w if (w in KEEP_UPPER or any(c.isdigit() for c in w)) else w.lower()) for w in words)
    t = t.replace("TGF \u03b2", "TGF-\u03b2")
    if words and words[0] not in CASED:
        t = t[:1].upper() + t[1:]
    return t if len(t) <= width else t[: width - 1] + "…"


KEEP_UPPER = {"DNA", "RNA", "MRNA", "TNFA", "NFKB", "TGF", "IL2", "IL6", "STAT3", "STAT5", "JAK", "MTORC1",
              "KRAS", "P53", "MYC", "E2F", "G2M", "UV", "PI3K", "AKT", "MTOR", "WNT", "ATP", "GTP", "NAD",
              "TCA", "ER", "UPR", "BCL2", "MCL1", "BH3", "TP53", "HIF", "ROS", "NK", "MHC", "II", "I", "III",
              "ECM", "EMT", "IFN", "GPCR", "MAPK", "ERK", "NOTCH", "HEDGEHOG", "KEGG", "SUMO", "RHO", "RAS"}

# Words that are neither upper nor lower case.
CASED = {"RRNA": "rRNA", "MRNA": "mRNA", "TRNA": "tRNA", "NCRNA": "ncRNA", "SNRNA": "snRNA", "MIRNA": "miRNA",
         "LNCRNA": "lncRNA", "SNORNA": "snoRNA", "DSDNA": "dsDNA", "SSDNA": "ssDNA", "NMD": "NMD", "SRP": "SRP",
         "GTPASE": "GTPase", "ATPASE": "ATPase", "KINASE": "kinase", "NFKB": "NF-κB", "TNFA": "TNFα",
         "IFNG": "IFNγ", "IFNA": "IFNα", "TGF": "TGF", "BETA": "β", "ALPHA": "α", "GAMMA": "γ"}
