"""One model's samples at a glance: blind VST, normalised counts, size factors,
PCA and sample-to-sample distances. The first thing to look at before any
contrast: replicates should sit together, and an outlier or a swapped sample
shows here."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from de_common import design_factors, log_to, model_samples, read_counts, vst  # noqa: E402
from smcore import plotting as sp  # noqa: E402

log_to(sm.log[0])
de = sm.params.de
model = sm.params.model
split = sm.wildcards.split
formats = de["plot_formats"]

counts, genes = read_counts(sm.input.counts, sm.input.tx2gene)
meta = model_samples(sm.input.samples, model, split)
counts = counts.loc[counts[meta.index].sum(axis=1) > 0, meta.index]

expr, dds = vst(counts, meta, model["design"], sm.threads)
expr.insert(0, "gene_name", genes.loc[expr.index, "gene_name"])
expr.to_csv(sm.output.vst, sep="\t", float_format="%.5g")
expr = expr.drop(columns="gene_name")
size_factors = dds.obs["size_factors"].rename("size_factor")
size_factors.to_csv(sm.output.size_factors, sep="\t", float_format="%.5g")
normed = counts.div(size_factors, axis=1)
normed.insert(0, "gene_name", genes.loc[normed.index, "gene_name"])
normed.to_csv(sm.output.normalized, sep="\t", float_format="%.5g")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402

sp.apply_style(font_size=14)
base = sm.output.vst.rsplit("/", 1)[0] + "/"
factors = design_factors(model["design"])
# Colour by the factor the contrasts test; a second design factor, if any,
# becomes the marker shape (e.g. replicate / flask in a paired design).
colour_by = model["contrasts"][0][0]
shape_by = next((f for f in factors if f != colour_by), None)
levels = list(dict.fromkeys(l for c in model["contrasts"] for l in (c[2], c[1]) if c[0] == colour_by))
levels += [l for l in sorted(meta[colour_by].unique()) if l not in levels]
level_color = sp.colors_for(levels)
title = "all samples" if split == "all" else split

# PCA on the 500 most variable genes, as DESeq2's plotPCA.
top = expr.var(axis=1).sort_values(ascending=False).index[:500]
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
ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")
ax.set_title(f"PCA · {title}")
handles = [plt.Line2D([], [], marker="o", ls="", color=level_color[l], markersize=10, label=l) for l in levels]
if shape_by:
    handles += [plt.Line2D([], [], ls="", label="")]
    handles += [plt.Line2D([], [], marker=m, ls="", color=sp.PALETTE["grey_dark"], markersize=9, label=f"{shape_by} {v}")
                for v, m in shapes.items()]
ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.04, 0.5), title=colour_by, alignment="left")
sp.save(fig, base + "pca", formats)

# Sample distances: Euclidean on all variance-stabilised genes, clustered.
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
ax.set_title(f"Sample distances \u00b7 {title}", pad=30)
for fmt in formats:
    fig.savefig(f"{base}sample_distances.{fmt}", dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)
