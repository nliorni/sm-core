"""One model's samples at a glance: the size factors of both normalisations
with each library's size and fraction of fragments in consensus peaks; blind
VST; PCA and sample-to-sample distances. Replicates should sit together; an
outlier or a swapped sample shows here."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from da_common import LABEL, METHODS, library_sizes, read_peak_counts, size_factor_table  # noqa: E402
from smcore import de as sde  # noqa: E402

sde.log_to(sm.log[0])
da = sm.params.da
model = sm.params.model
split = sm.wildcards.split
formats = da["plot_formats"]
used = da["normalization"]

counts, coords = read_peak_counts(sm.input.counts)
meta = sde.model_samples(sm.input.samples, model, split)
counts = counts.loc[counts[meta.index].sum(axis=1) > 0, meta.index]
lib = library_sizes(sm.input.fc_summary, sm.params.bams)

sf = size_factor_table(counts, lib, used)
sf.to_csv(sm.output.size_factors, sep="\t", float_format="%.6g")
print(f"{split}: {len(meta)} samples, {len(counts)} peaks; normalisation {used}")
print(sf.to_string())

expr, _ = sde.vst(counts, meta, model["design"], sm.threads, size_factors=sf["size_factor"])
expr.to_csv(sm.output.vst, sep="\t", float_format="%.5g")
counts.div(sf["size_factor"], axis=1).to_csv(sm.output.normalized, sep="\t", float_format="%.5g")

base = sm.output.vst.rsplit("/", 1)[0] + "/"
title = "all samples" if split == "all" else split
sde.pca(expr, meta, model, base + "pca", formats, title=f"PCA · {title}")
sde.sample_distances(expr, meta, model, base + "sample_distances", formats,
                     title=f"Sample distances · {title}")

# Size factors of both methods side by side, with the fraction of fragments
# in peaks: where the two disagree, it is usually because that fraction
# (signal to noise) differs between libraries.
import matplotlib.pyplot as plt  # noqa: E402
from smcore import plotting as sp  # noqa: E402

sp.apply_style(font_size=13)
n = len(sf)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(max(9, 0.5 * n + 7), max(4, 0.32 * n + 1.8)), sharey=True,
                              gridspec_kw={"width_ratios": [3, 2]})
y = np.arange(n)[::-1]
colours = {"median_of_ratios": sp.PALETTE["blue"], "library_size": sp.PALETTE["red"]}
for m in METHODS:
    ax.scatter(np.log2(sf[f"sf_{m}"]), y, s=70 if m == used else 45, color=colours[m],
               marker="o" if m == used else "D", zorder=3,
               label=LABEL[m] + (" (used)" if m == used else ""))
ax.axvline(0, color=sp.PALETTE["grey_dark"], lw=1)
ax.set_yticks(y)
ax.set_yticklabels(sf.index, fontsize=10)
ax.set_xlabel("log$_2$ size factor")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, borderaxespad=0, fontsize=10)
bars = ax2.barh(y, sf["fraction_in_peaks"], color=sp.PALETTE["grey"], height=0.7)
sp.label_bars(ax2, bars, [f"{v:.2f}" for v in sf["fraction_in_peaks"]])
ax2.set_xlim(0, min(1, sf["fraction_in_peaks"].max() * 1.35))
ax2.set_xlabel("fraction of fragments in consensus peaks")
ax2.tick_params(axis="y", length=0)
ax.set_title(f"Normalisation · {title}", loc="left", pad=36)
sp.save(fig, base + "size_factors", formats)
