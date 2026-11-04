"""sm-core differential analysis: PyDESeq2 models and their figures, shared by
rnaseq's differential expression (`-w de`) and atacseq's differential
accessibility (`-w da`).

A "feature" is a gene in rnaseq and a consensus peak in atacseq. Everything
here works on a features x samples count table and a sample table; what a
feature is, how it is labelled and what is joined to the results stays in the
pipeline's own scripts.

Rule scripts run in the pipeline's DE conda environment (PyDESeq2, pandas,
matplotlib, adjustText) and import this module by putting the sm-core root on
sys.path, as for smcore.plotting.
"""

from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

# ============================================================================
# Models
# ============================================================================


def log_to(path):
    """Send this script's output (and PyDESeq2's) to the rule's log."""
    handle = open(path, "w")
    sys.stdout = sys.stderr = handle
    return handle


def design_factors(design: str) -> list[str]:
    """Column names in a formula: identifiers, minus formula keywords."""
    names = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", design.split("~", 1)[1])
    return [n for n in dict.fromkeys(names) if n not in ("C", "Q")]


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
    """Keep features with >= min_count reads in at least as many samples as the
    smallest group of `factor` (DESeq2 vignette's pre-filter)."""
    smallest = int(meta[factor].value_counts().min())
    keep = (counts[meta.index] >= min_count).sum(axis=1) >= smallest
    return counts.loc[keep, meta.index]


def library_size_factors(library_sizes: pd.Series) -> pd.Series:
    """Size factors proportional to library size, scaled to a geometric mean of
    1 like DESeq2's (the scale cancels in every fold change)."""
    sizes = library_sizes.astype(float)
    return sizes / np.exp(np.log(sizes).mean())


def _set_size_factors(dds, size_factors: pd.Series) -> None:
    """Give a DeseqDataSet fixed size factors.

    PyDESeq2 0.5 has no option for this: deseq2(), vst_fit() and the
    dispersion initialisation all (re)compute median-of-ratios factors through
    fit_size_factors(). Setting the factors and the normalised counts, then
    turning that method into a no-op on this object, makes every later step
    use the given factors. Checked in tests: given PyDESeq2's own
    median-of-ratios factors, the results equal a plain deseq2() run."""
    sf = size_factors.reindex(dds.obs_names).to_numpy(dtype=float)
    if np.isnan(sf).any() or (sf <= 0).any():
        raise ValueError("size factors must be positive and given for every sample")
    # The geometric means of the features, which vst_fit wants to exist, as
    # fit_size_factors would leave them.
    counts = np.asarray(dds.X, dtype=float)
    with np.errstate(divide="ignore"):
        logmeans = np.log(counts).mean(axis=0)
    dds.logmeans = logmeans
    dds.filtered_genes = ~np.isinf(logmeans)
    dds.obs["size_factors"] = sf
    dds.layers["normed_counts"] = counts / sf[:, None]
    dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
    dds.fit_size_factors = lambda *args, **kwargs: None


def fit(counts: pd.DataFrame, meta: pd.DataFrame, design: str, threads: int,
        size_factors: pd.Series | None = None):
    """PyDESeq2 fit of one model. `size_factors` (one per sample) replaces
    DESeq2's median-of-ratios normalisation when given."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    inference = DefaultInference(n_cpus=threads)
    dds = DeseqDataSet(
        counts=counts.T,  # PyDESeq2 wants samples x features
        metadata=meta,
        design=design,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    if size_factors is not None:
        _set_size_factors(dds, size_factors)
    dds.deseq2()
    return dds, inference


def vst(counts: pd.DataFrame, meta: pd.DataFrame, design: str, threads: int,
        size_factors: pd.Series | None = None):
    """Blind variance-stabilised values (features x samples) for QC plots."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    dds = DeseqDataSet(
        counts=counts.T, metadata=meta, design=design,
        inference=DefaultInference(n_cpus=threads), quiet=True,
    )
    if size_factors is not None:
        _set_size_factors(dds, size_factors)
    dds.vst(use_design=False)
    return pd.DataFrame(dds.layers["vst_counts"].T, index=counts.index, columns=counts.columns), dds


def regulation(res: pd.DataFrame, alpha: float, min_abs_log2fc: float) -> np.ndarray:
    """'up' / 'down' / 'ns' from the adjusted p-value and the (unshrunken)
    log2 fold change."""
    significant = (res["padj"] < alpha) & (res["log2FoldChange"].abs() >= min_abs_log2fc)
    return np.where(significant & (res["log2FoldChange"] > 0), "up", np.where(significant, "down", "ns"))


def neg_log10(p: pd.Series) -> pd.Series:
    """-log10(p) with p = 0 (below double precision) set just past the smallest
    non-zero value, so those features stay on the plot."""
    p = p.astype(float)
    floor = p[p > 0].min() if (p > 0).any() else 1e-300
    return -np.log10(p.clip(lower=floor / 10))


def bh(p, n_tests: int) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values for `p`, counting `n_tests` tests in all.

    `n_tests` can exceed len(p): tests left out of `p` (e.g. gene sets with no
    overlap, which gseapy enrich does not return) are counted with p = 1."""
    p = np.asarray(p, dtype=float)
    if p.size == 0:
        return p
    n_tests = max(n_tests, p.size)
    order = np.argsort(p)
    ranked = p[order] * n_tests / np.arange(1, p.size + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def short_term(term: str, width: int = 48) -> str:
    """'HALLMARK_TNFA_SIGNALING_VIA_NFKB' -> 'TNFA signaling via NFKB'."""
    t = re.sub(r"^(HALLMARK|REACTOME|GOBP|GOCC|GOMF|KEGG|WP)_", "", term)
    words = t.split("_")
    # Gene and complex names stay upper case; ordinary words do not.
    t = " ".join(CASED.get(w, w if (w in KEEP_UPPER or any(c.isdigit() for c in w)) else w.lower()) for w in words)
    t = t.replace("TGF β", "TGF-β")
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


# ============================================================================
# Figures (sm-core style)
# ============================================================================


def _counts_note(ax, n_up, n_down):
    # Above the axes, so no data point or label can collide with it.
    from smcore import plotting as sp

    kw = dict(transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")
    ax.text(1.0, 1.02, f"↑ {n_up:,}", color=sp.UP, ha="right", **kw)
    ax.text(0.80, 1.02, f"↓ {n_down:,}", color=sp.DOWN, ha="right", **kw)


def volcano(res: pd.DataFrame, path: str, formats, *, label_col: str, alpha: float, min_abs_log2fc: float,
            label_top: int, title: str, xlabel: str):
    """Volcano of one contrast. The x axis is the log2 fold change the
    significance calls use. Adjusted p-values below double precision (reported
    as 0) are drawn as triangles at the top instead of piling onto one line.
    The labelled features are those both strongly changed and strongly
    significant, half on each side."""
    import matplotlib.pyplot as plt
    from smcore import plotting as sp

    sp.apply_style(font_size=14)
    colors = res["regulation"].map({"up": sp.UP, "down": sp.DOWN, "ns": sp.NS})
    n_up, n_down = int((res["regulation"] == "up").sum()), int((res["regulation"] == "down").sum())
    tested = res.dropna(subset=["padj"])
    capped = tested["padj"] <= 0
    y = neg_log10(tested["padj"])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for is_capped, marker, size in ((False, "o", 10), (True, "^", 18)):
        part = tested[capped == is_capped]
        order = part["regulation"].map({"ns": 0, "down": 1, "up": 1}).sort_values().index
        ax.scatter(part.loc[order, "log2FoldChange"], y[order], c=colors[order], s=size, marker=marker, alpha=0.8,
                   linewidths=0, rasterized=True)
    for x0 in (-min_abs_log2fc, min_abs_log2fc):
        if x0:
            ax.axvline(x0, color=sp.PALETTE["grey_dark"], lw=1, ls="--", zorder=0)
    ax.axhline(-np.log10(alpha), color=sp.PALETTE["grey_dark"], lw=1, ls="--", zorder=0)
    sig = tested[tested["regulation"] != "ns"].copy()
    # Significance capped at its 90th percentile, so among the very significant
    # features the fold change decides and labels spread across the plot.
    y_eff = y[sig.index].clip(upper=y[sig.index].quantile(0.9)) if len(sig) else y[sig.index]
    sig["score"] = sig["log2FoldChange"].abs() * y_eff
    half = max(1, label_top // 2)
    # Only features with a label compete for the labels.
    sig = sig[sig[label_col].fillna("").astype(str) != ""]
    top = (pd.concat([sig[sig["regulation"] == d].nlargest(half, "score") for d in ("up", "down")])
           if label_top else sig.iloc[:0])
    if len(top):
        from adjustText import adjust_text

        texts = [ax.text(r["log2FoldChange"], y[g], r[label_col], fontsize=10, color=sp.PALETTE["ink"])
                 for g, r in top.iterrows()]
        adjust_text(texts, ax=ax, expand=(1.3, 1.6),
                    arrowprops=dict(arrowstyle="-", color=sp.PALETTE["grey_dark"], lw=0.6))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("$-$log$_{10}$ adjusted p")
    ax.set_title(title)
    _counts_note(ax, n_up, n_down)
    if capped.any():
        ax.text(0, -0.16, "▲ adjusted p below double precision (< 1e-300), drawn at the top",
                transform=ax.transAxes, fontsize=9, va="top", color=sp.PALETTE["grey_dark"])
    sp.save(fig, path, formats)


def ma(res: pd.DataFrame, path: str, formats, *, lfc_col: str, ylabel: str, xlabel: str, title: str,
       reference_line: tuple[float, str] | None = None):
    """MA plot: mean normalised count against the (shrunken) log2 fold change.
    `reference_line` = (y, label) draws one more horizontal line, e.g. where
    another normalisation would put zero."""
    import matplotlib.pyplot as plt
    from smcore import plotting as sp

    sp.apply_style(font_size=14)
    colors = res["regulation"].map({"up": sp.UP, "down": sp.DOWN, "ns": sp.NS})
    n_up, n_down = int((res["regulation"] == "up").sum()), int((res["regulation"] == "down").sum())
    fig, ax = plt.subplots(figsize=(7, 5))
    data = res.dropna(subset=[lfc_col])
    order = data["regulation"].map({"ns": 0, "down": 1, "up": 1}).sort_values().index
    ax.scatter(np.log10(data.loc[order, "baseMean"] + 1), data.loc[order, lfc_col], c=colors[order], s=10,
               alpha=0.8, linewidths=0, rasterized=True)
    ax.axhline(0, color=sp.PALETTE["ink"], lw=1.2)
    if reference_line is not None:
        y0, label = reference_line
        ax.axhline(y0, color=sp.PALETTE["gold"], lw=1.4, ls="--")
        ax.text(0.99, y0, label, transform=ax.get_yaxis_transform(), fontsize=9, ha="right",
                va="bottom" if y0 >= 0 else "top", color=sp.PALETTE["grey_dark"],
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _counts_note(ax, n_up, n_down)
    sp.save(fig, path, formats)


def heatmap(expr: pd.DataFrame, labels: list[str], meta: pd.DataFrame, factor: str, levels: list[str],
            path: str, formats, *, title: str):
    """Z-scored variance-stabilised values of the given features (rows of
    `expr`, already selected and ordered) across the model's samples, grouped
    by `factor` in the order of `levels`."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from smcore import plotting as sp

    sp.apply_style(font_size=14)
    cols = [s for level in levels for s in meta.index[meta[factor] == level]]
    z = expr[cols]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    height = max(4, 0.2 * len(z) + 2.2)
    width = max(6, 0.35 * len(cols) + 3.5)
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(z.values, aspect="auto", cmap=sp.DIVERGING, vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(labels, fontsize=8 if len(z) > 30 else 10)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=9)
    ax.tick_params(length=0)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    # Group strip above the columns, colour bar to the right: both in their own
    # inset axes, so neither can overlap the heatmap.
    level_color = sp.colors_for(levels)
    strip = ax.inset_axes([0, 1.01, 1, 0.025 * 50 / max(len(z), 10)])
    strip.imshow([[to_rgb(level_color[meta.loc[s, factor]]) for s in cols]], aspect="auto")
    strip.set_axis_off()
    cax = ax.inset_axes([1.03, 0.62, 0.035 * 8 / max(len(cols), 8) + 0.02, 0.36])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("z-score (VST)", fontsize=11)
    cbar.outline.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=level_color[l]) for l in levels]
    ax.legend(handles, levels, loc="upper left", bbox_to_anchor=(1.03, 0.56), fontsize=10, title=factor,
              title_fontsize=10, handlelength=1, borderaxespad=0, alignment="left")
    ax.set_title(title, pad=14 + 600 / max(len(z), 10))
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def overview_levels(model: dict, meta: pd.DataFrame):
    """The factor to colour samples by (the one the contrasts test), a second
    design factor for the marker shape (e.g. replicate in a paired design),
    and the colour factor's levels, references first."""
    factors = design_factors(model["design"])
    colour_by = model["contrasts"][0][0]
    shape_by = next((f for f in factors if f != colour_by), None)
    levels = list(dict.fromkeys(l for c in model["contrasts"] for l in (c[2], c[1]) if c[0] == colour_by))
    levels += [l for l in sorted(meta[colour_by].unique()) if l not in levels]
    return colour_by, shape_by, levels


def pca(expr: pd.DataFrame, meta: pd.DataFrame, model: dict, path: str, formats, *, title: str,
        n_top: int = 500, what: str = "features"):
    """PCA on the `n_top` most variable rows of `expr`, as DESeq2's plotPCA."""
    import matplotlib.pyplot as plt
    from smcore import plotting as sp

    sp.apply_style(font_size=14)
    colour_by, shape_by, levels = overview_levels(model, meta)
    level_color = sp.colors_for(levels)
    top = expr.var(axis=1).sort_values(ascending=False).index[:n_top]
    x = expr.loc[top].T
    x = x - x.mean()
    u, s, _ = np.linalg.svd(x.values, full_matrices=False)
    pcs = u * s
    explained = s**2 / (s**2).sum() * 100
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]
    shapes = ({v: markers[i % len(markers)] for i, v in enumerate(sorted(meta[shape_by].unique()))}
              if shape_by else {})
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for i, sample in enumerate(meta.index):
        ax.scatter(pcs[i, 0], pcs[i, 1], s=110, color=level_color[meta.loc[sample, colour_by]],
                   marker=shapes.get(meta.loc[sample, shape_by], "o") if shape_by else "o",
                   edgecolors="white", linewidths=1.2, zorder=3)
    ax.margins(x=0.15, y=0.12)
    if len(meta) <= 40:
        from adjustText import adjust_text

        texts = [ax.text(pcs[i, 0], pcs[i, 1], sample, fontsize=8, color=sp.PALETTE["grey_dark"])
                 for i, sample in enumerate(meta.index)]
        adjust_text(texts, x=pcs[:, 0], y=pcs[:, 1], ax=ax, expand=(1.3, 1.8),
                    arrowprops=dict(arrowstyle="-", color=sp.PALETTE["grey"], lw=0.5))
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)" if len(explained) > 1 else "PC2")
    ax.set_title(title)
    handles = [plt.Line2D([], [], marker="o", ls="", color=level_color[l], markersize=10, label=l) for l in levels]
    if shape_by:
        handles += [plt.Line2D([], [], ls="", label="")]
        handles += [plt.Line2D([], [], marker=m, ls="", color=sp.PALETTE["grey_dark"], markersize=9,
                               label=f"{shape_by} {v}") for v, m in shapes.items()]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.04, 0.5), title=colour_by, alignment="left")
    sp.save(fig, path, formats)


def sample_distances(expr: pd.DataFrame, meta: pd.DataFrame, model: dict, path: str, formats, *, title: str):
    """Euclidean distances between samples on all variance-stabilised rows,
    clustered (average linkage)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist, squareform
    from smcore import plotting as sp

    sp.apply_style(font_size=14)
    colour_by, _, levels = overview_levels(model, meta)
    level_color = sp.colors_for(levels)
    dist = squareform(pdist(expr.T.values))
    order = leaves_list(linkage(pdist(expr.T.values), method="average"))
    names = meta.index[order]
    n = len(names)
    size = max(5, 0.3 * n + 3)
    fig, ax = plt.subplots(figsize=(size + 1.5, size))
    im = ax.imshow(dist[np.ix_(order, order)], cmap=sp.SEQUENTIAL_BLUE.reversed(), interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(names, fontsize=9)
    ax.tick_params(length=0, pad=14)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    # Group strips along both edges, in inset axes so they line up with the cells.
    colours = [to_rgb(level_color[meta.loc[s, colour_by]]) for s in names]
    left = ax.inset_axes([-0.035, 0, 0.025, 1])
    left.imshow([[c] for c in colours], aspect="auto")
    left.set_axis_off()
    bottom = ax.inset_axes([0, -0.035, 1, 0.025])
    bottom.imshow([colours], aspect="auto")
    bottom.set_axis_off()
    cax = ax.inset_axes([1.03, 0.6, 0.03, 0.4])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("distance (VST)", fontsize=11)
    cbar.outline.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=level_color[l]) for l in levels]
    ax.legend(handles, levels, ncol=min(len(levels), 6), loc="lower left", bbox_to_anchor=(0, 1.02), fontsize=10,
              handlelength=1, columnspacing=1, borderaxespad=0)
    ax.set_title(title, pad=30)
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def summary_bars(summary: pd.DataFrame, path: str, formats, *, xlabel: str, title: str):
    """Up (right) and down (left) counts for every analysis, split and contrast."""
    import matplotlib.pyplot as plt
    from smcore import plotting as sp

    sp.apply_style(font_size=13)
    labels = [
        c.replace("_vs_", " vs ") + ("" if s == "all" else f"  · {s}")
        + ("" if len(set(summary["analysis"])) == 1 else f"  [{a}]")
        for a, s, c in zip(summary["analysis"], summary["split"], summary["contrast"])
    ]
    n = len(summary)
    fig, ax = plt.subplots(figsize=(10, 0.45 * n + 2))
    y = range(n)[::-1]
    up = ax.barh(list(y), summary["up"], color=sp.UP, height=0.7, label="up")
    down = ax.barh(list(y), -summary["down"], color=sp.DOWN, height=0.7, label="down")
    limit = max(summary["up"].max(), summary["down"].max(), 1) * 1.3
    ax.set_xlim(-limit, limit)
    sp.label_bars(ax, up, [f"{v:,}" for v in summary["up"]])
    sp.label_bars(ax, down, [f"{v:,}" for v in summary["down"]])
    ax.axvline(0, color=sp.PALETTE["ink"], lw=1.2)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=11)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(int(v)):,}"))
    ax.set_xlabel(xlabel)
    ax.legend(loc="lower right", bbox_to_anchor=(1, 1.01), ncol=2, borderaxespad=0)
    ax.set_title(title)
    sp.save(fig, path, formats)
