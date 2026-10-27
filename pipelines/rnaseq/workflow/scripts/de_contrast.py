"""One contrast of one model: PyDESeq2 Wald test, shrunken log2 fold changes,
significant up/down lists, and the volcano, MA and top-gene heatmap."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from de_common import (  # noqa: E402
    fit, log_to, model_samples, neg_log10, prefilter, read_counts, vst, with_reference,
)
from smcore import plotting as sp  # noqa: E402

log_to(sm.log[0])
de = sm.params.de
model = sm.params.model
factor, test, reference = model["contrast"]
split = sm.wildcards.split
formats = de["plot_formats"]

counts, genes = read_counts(sm.input.counts, sm.input.tx2gene)
meta = with_reference(model_samples(sm.input.samples, model, split), factor, reference)
counts = prefilter(counts, meta, factor, de["min_count"])
print(f"{split}: {len(meta)} samples, {counts.shape[0]} genes after the pre-filter, design {model['design']}")

dds, inference = fit(counts, meta, model["design"], sm.threads)

from pydeseq2.ds import DeseqStats  # noqa: E402

stats = DeseqStats(dds, contrast=[factor, test, reference], alpha=de["alpha"], inference=inference, quiet=True)
stats.summary()
res = stats.results_df.copy()
if de["shrink"]:
    coeff = f"{factor}[T.{test}]"
    stats.lfc_shrink(coeff=coeff)
    res["log2FoldChange_shrunken"] = stats.results_df["log2FoldChange"]

normed = pd.DataFrame(dds.layers["normed_counts"].T, index=counts.index, columns=meta.index)
for level in (test, reference):
    res[f"mean_normalized_{level}"] = normed[meta.index[meta[factor] == level]].mean(axis=1)

significant = (res["padj"] < de["alpha"]) & (res["log2FoldChange"].abs() >= de["min_abs_log2fc"])
res["regulation"] = np.where(significant & (res["log2FoldChange"] > 0), "up",
                             np.where(significant, "down", "ns"))
res = genes.join(res, how="inner").sort_values(["padj", "pvalue"], na_position="last")
res.index.name = "gene_id"
res.to_csv(sm.output.results, sep="\t", float_format="%.6g")
res[res["regulation"] == "up"].to_csv(sm.output.up, sep="\t", float_format="%.6g")
res[res["regulation"] == "down"].to_csv(sm.output.down, sep="\t", float_format="%.6g")
n_up, n_down = int((res["regulation"] == "up").sum()), int((res["regulation"] == "down").sum())
print(f"{test} vs {reference}: {n_up} up, {n_down} down (padj < {de['alpha']}, |log2FC| >= {de['min_abs_log2fc']})")

# ---------------------------------------------------------------- plots
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

sp.apply_style(font_size=14)
title = f"{test} vs {reference}" + ("" if split == "all" else f" · {split}")
colors = res["regulation"].map({"up": sp.UP, "down": sp.DOWN, "ns": sp.NS})
base = sm.output.results.rsplit("/", 1)[0] + "/"


def counts_note(ax):
    # Above the axes, so no data point or label can collide with it.
    kw = dict(transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")
    ax.text(1.0, 1.02, f"\u2191 {n_up:,}", color=sp.UP, ha="right", **kw)
    ax.text(0.80, 1.02, f"\u2193 {n_down:,}", color=sp.DOWN, ha="right", **kw)


# Volcano: the x axis is the log2 fold change the significance calls use.
# Adjusted p-values below double precision (reported as 0) are drawn as
# triangles at the top instead of all piling onto one line.
tested = res.dropna(subset=["padj"])
capped = tested["padj"] <= 0
y = neg_log10(tested["padj"])
fig, ax = plt.subplots(figsize=(7.5, 6))
for is_capped, marker, size in ((False, "o", 10), (True, "^", 18)):
    part = tested[capped == is_capped]
    order = part["regulation"].map({"ns": 0, "down": 1, "up": 1}).sort_values().index
    ax.scatter(part.loc[order, "log2FoldChange"], y[order], c=colors[order], s=size, marker=marker, alpha=0.8,
               linewidths=0, rasterized=True)
for x0 in (-de["min_abs_log2fc"], de["min_abs_log2fc"]):
    if x0:
        ax.axvline(x0, color=sp.PALETTE["grey_dark"], lw=1, ls="--", zorder=0)
ax.axhline(-np.log10(de["alpha"]), color=sp.PALETTE["grey_dark"], lw=1, ls="--", zorder=0)
# Label the genes that are both strongly changed and strongly significant,
# half of them on each side, rather than only the top of the p-value list
# (which on a large effect is a single crowded line at the cap).
sig = tested[tested["regulation"] != "ns"].copy()
# Significance capped at its 90th percentile, so among the very significant
# genes the fold change decides and labels spread across the plot.
y_eff = y[sig.index].clip(upper=y[sig.index].quantile(0.9)) if len(sig) else y[sig.index]
sig["score"] = sig["log2FoldChange"].abs() * y_eff
half = max(1, de["label_top"] // 2)
top = pd.concat([sig[sig["regulation"] == d].nlargest(half, "score") for d in ("up", "down")]) if de["label_top"] else sig.iloc[:0]
if len(top):
    from adjustText import adjust_text

    texts = [ax.text(r.log2FoldChange, y[g], r.gene_name, fontsize=10, color=sp.PALETTE["ink"])
             for g, r in top.iterrows()]
    adjust_text(texts, ax=ax, expand=(1.3, 1.6),
                arrowprops=dict(arrowstyle="-", color=sp.PALETTE["grey_dark"], lw=0.6))
ax.set_xlabel(f"log$_2$ fold change ({test} / {reference})")
ax.set_ylabel("$-$log$_{10}$ adjusted p")
ax.set_title(title)
counts_note(ax)
if capped.any():
    ax.text(0, -0.16, "\u25b2 adjusted p below double precision (< 1e-300), drawn at the top",
            transform=ax.transAxes, fontsize=9, va="top", color=sp.PALETTE["grey_dark"])
sp.save(fig, base + "volcano", formats)

# MA: mean expression against the shrunken fold change, the estimate that is
# not inflated for weakly expressed genes.
lfc = "log2FoldChange_shrunken" if de["shrink"] else "log2FoldChange"
fig, ax = plt.subplots(figsize=(7, 5))
ma = res.dropna(subset=[lfc])
order = ma["regulation"].map({"ns": 0, "down": 1, "up": 1}).sort_values().index
ax.scatter(np.log10(ma.loc[order, "baseMean"] + 1), ma.loc[order, lfc], c=colors[order], s=10, alpha=0.8,
           linewidths=0, rasterized=True)
ax.axhline(0, color=sp.PALETTE["ink"], lw=1.2)
ax.set_xlabel("log$_{10}$ mean normalised count")
ax.set_ylabel("log$_2$ fold change" + (" (shrunken)" if de["shrink"] else ""))
ax.set_title(title)
counts_note(ax)
sp.save(fig, base + "ma", formats)

# Heatmap: the most significant genes across every sample of the model,
# z-scored variance-stabilised expression, samples grouped by the factor.
expr, _ = vst(counts, meta, sm.params.model["design"], sm.threads)
picked = res[res["regulation"] != "ns"].head(de["heatmap_top"])
if len(picked) < 2:
    picked = res.dropna(subset=["padj"]).head(de["heatmap_top"])
levels = [reference, test] + [l for l in meta[factor].cat.categories if l not in (reference, test)]
cols = [s for level in levels for s in meta.index[meta[factor] == level]]
z = expr.loc[picked.index, cols]
z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
height = max(4, 0.2 * len(z) + 2.2)
width = max(6, 0.35 * len(cols) + 3.5)
fig, ax = plt.subplots(figsize=(width, height))
im = ax.imshow(z.values, aspect="auto", cmap=sp.DIVERGING, vmin=-2.5, vmax=2.5, interpolation="nearest")
ax.set_yticks(range(len(z)))
ax.set_yticklabels(picked["gene_name"].fillna("").tolist(), fontsize=8 if len(z) > 30 else 10)
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
ax.set_title(f"Top {len(z)} genes \u00b7 {title}", pad=14 + 600 / max(len(z), 10))
for fmt in formats:
    fig.savefig(f"{base}heatmap.{fmt}", dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)
