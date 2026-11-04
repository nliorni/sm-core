"""One contrast of one model: PyDESeq2 Wald test on the consensus peaks,
shrunken log2 fold changes, the significant peaks (TSV and BED), a per-gene
summary, and the volcano, MA, top-peak heatmap and peak-class figures."""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from da_common import LABEL, METHODS, library_sizes, read_peak_counts, size_factor_table, zero_shift  # noqa: E402
from smcore import de as sde  # noqa: E402

sde.log_to(sm.log[0])
da = sm.params.da
model = sm.params.model
factor, test, reference = model["contrast"]
split = sm.wildcards.split
formats = da["plot_formats"]
used = da["normalization"]
other = next(m for m in METHODS if m != used)

counts, coords = read_peak_counts(sm.input.counts)
ann = pd.read_csv(sm.input.annotation, sep="\t", index_col="peak_id",
                  dtype={"promoter_gene_ids": str, "promoter_gene_names": str})
meta = sde.with_reference(sde.model_samples(sm.input.samples, model, split), factor, reference)
counts = sde.prefilter(counts, meta, factor, da["min_count"])
lib = library_sizes(sm.input.fc_summary, sm.params.bams)
sf = size_factor_table(counts, lib, used)
print(f"{split}: {len(meta)} samples, {counts.shape[0]} peaks after the pre-filter, design {model['design']}, "
      f"normalisation {used}")

# median_of_ratios: PyDESeq2's own (identical to sf_median_of_ratios);
# library_size: fixed factors.
dds, inference = sde.fit(counts, meta, model["design"], sm.threads,
                         size_factors=sf["size_factor"] if used == "library_size" else None)

from pydeseq2.ds import DeseqStats  # noqa: E402

stats = DeseqStats(dds, contrast=[factor, test, reference], alpha=da["alpha"], inference=inference, quiet=True)
stats.summary()
res = stats.results_df.copy()
if da["shrink"]:
    stats.lfc_shrink(coeff=f"{factor}[T.{test}]")
    res["log2FoldChange_shrunken"] = stats.results_df["log2FoldChange"]

normed = pd.DataFrame(dds.layers["normed_counts"].T, index=counts.index, columns=meta.index)
test_samples = meta.index[meta[factor] == test]
reference_samples = meta.index[meta[factor] == reference]
res[f"mean_normalized_{test}"] = normed[test_samples].mean(axis=1)
res[f"mean_normalized_{reference}"] = normed[reference_samples].mean(axis=1)
res["regulation"] = sde.regulation(res, da["alpha"], da["min_abs_log2fc"])

ann_cols = ["width", "nearest_gene_id", "nearest_gene_name", "nearest_gene_type", "distance_to_tss", "annotation",
            "promoter_gene_names", "promoter_gene_ids"]
res = coords.join(ann[ann_cols]).join(res, how="inner").sort_values(["padj", "pvalue"], na_position="last")
res.index.name = "peak_id"
res.to_csv(sm.output.results, sep="\t", float_format="%.6g")
for direction, path, bed in (("up", sm.output.up, sm.output.up_bed), ("down", sm.output.down, sm.output.down_bed)):
    hits = res[res["regulation"] == direction]
    hits.to_csv(path, sep="\t", float_format="%.6g")
    # BED6 in genome order, for a browser or a motif tool; the score is
    # -10 log10 padj, capped at 1000.
    score = (-10 * np.log10(hits["padj"].clip(lower=1e-100))).clip(upper=1000).round().astype(int)
    bed_df = pd.DataFrame({"chrom": hits["chrom"], "start": hits["start"], "end": hits["end"],
                           "name": hits.index, "score": score, "strand": "."})
    bed_df.sort_values(["chrom", "start"]).to_csv(bed, sep="\t", header=False, index=False)

n_up, n_down = int((res["regulation"] == "up").sum()), int((res["regulation"] == "down").sum())
fitted = pd.Series(np.asarray(dds.obs["size_factors"], dtype=float), index=dds.obs_names)
shift = zero_shift(fitted, sf[f"sf_{other}"], test_samples, reference_samples)
print(f"{test} vs {reference}: {n_up} peaks up, {n_down} down "
      f"(padj < {da['alpha']}, |log2FC| >= {da['min_abs_log2fc']})")
print(f"A peak unchanged under {LABEL[other]} normalisation has log2FC {shift:+.3f} under {LABEL[used]}.")

# ---------------------------------------------------------------- genes
# Per gene, its peaks: those it is the nearest gene of, and those covering its
# promoter (a peak can cover the promoters of several genes). How many were
# tested and changed, its most significant promoter peak, and its most
# significant peak of any class. The table to set next to the RNA-seq
# results of the same contrast.
tested = res.dropna(subset=["padj"])
pairs = pd.concat([
    tested["nearest_gene_id"].rename("gene_id").reset_index(),
    tested["promoter_gene_ids"].dropna().str.split(",").explode().rename("gene_id").reset_index(),
])
pairs = pairs[pairs["gene_id"].notna() & (pairs["gene_id"] != "")].drop_duplicates()
long = pairs.join(tested[["annotation", "promoter_gene_ids", "log2FoldChange", "padj", "regulation"]], on="peak_id")
long["promoter"] = [g in str(ids).split(",") for g, ids in zip(long["gene_id"], long["promoter_gene_ids"])]
groups = long.groupby("gene_id", sort=False)
info = pd.read_csv(sm.input.genes, sep="\t", index_col="gene_id", usecols=["gene_id", "gene_name", "gene_type"])
genes = pd.DataFrame({
    "peaks_tested": groups.size(),
    "peaks_up": groups["regulation"].agg(lambda r: int((r == "up").sum())),
    "peaks_down": groups["regulation"].agg(lambda r: int((r == "down").sum())),
})
genes = info.reindex(genes.index).join(genes)
for prefix, subset in (("best_peak", long), ("promoter", long[long["promoter"]])):
    # gene -> its peak with the smallest padj (within the subset)
    rows = subset.sort_values(["padj", "peak_id"]).drop_duplicates("gene_id").set_index("gene_id")
    genes["best_peak" if prefix == "best_peak" else "promoter_peak"] = rows["peak_id"]
    if prefix == "best_peak":
        genes["best_peak_class"] = rows["annotation"]
    for col in ("log2FoldChange", "padj", "regulation"):
        genes[f"{prefix}_{col}"] = rows[col]
genes.index.name = "gene_id"
genes = genes.sort_values(["promoter_padj", "best_peak_padj"], na_position="last")
genes.to_csv(sm.output.genes, sep="\t", float_format="%.6g")
print(f"{len(genes)} genes with a tested peak; "
      f"{int(genes['promoter_regulation'].isin(['up', 'down']).sum())} with a significant promoter peak")

# ---------------------------------------------------------------- plots
title = f"{test} vs {reference}" + ("" if split == "all" else f" · {split}")
base = sm.output.results.rsplit("/", 1)[0] + "/"
# Label: the peak's first named promoter gene, else its nearest gene if
# named; genes known only by their ID (ENSG...) tell a reader nothing.
named = lambda n: n if n and not re.match(r"ENS[A-Z]*G\d", n) else ""  # noqa: E731
res["label"] = [next((n for n in str(p).split(",") if named(n)), "") or named(str(g))
                for p, g in zip(res["promoter_gene_names"].fillna(""), res["nearest_gene_name"].fillna(""))]

sde.volcano(res, base + "volcano", formats, label_col="label", alpha=da["alpha"],
            min_abs_log2fc=da["min_abs_log2fc"], label_top=da["label_top"], title=title,
            xlabel=f"log$_2$ fold change ({test} / {reference})")

lfc = "log2FoldChange_shrunken" if da["shrink"] else "log2FoldChange"
sde.ma(res, base + "ma", formats, lfc_col=lfc, title=title,
       xlabel="log$_{10}$ mean normalised count",
       ylabel="log$_2$ fold change" + (" (shrunken)" if da["shrink"] else ""),
       reference_line=(shift, f"0 under {LABEL[other]} normalisation"))

expr, _ = sde.vst(counts, meta, model["design"], sm.threads,
                  size_factors=sf["size_factor"] if used == "library_size" else None)
picked = res[res["regulation"] != "ns"].head(da["heatmap_top"])
if len(picked) < 2:
    picked = res.dropna(subset=["padj"]).head(da["heatmap_top"])
levels = [reference, test] + [l for l in meta[factor].cat.categories if l not in (reference, test)]
labels = [f"{r.label or r.Index} · {r.annotation}" for r in picked.itertuples()]
sde.heatmap(expr.loc[picked.index], labels, meta, factor, levels, base + "heatmap", formats,
            title=f"Top {len(picked)} peaks · {title}")

# Where the peaks lie: promoter / genic / distal among the tested, the
# opened and the closed peaks. A shift towards distal peaks among the changed
# ones is typical of enhancer-driven responses.
import matplotlib.pyplot as plt  # noqa: E402
from smcore import plotting as sp  # noqa: E402

sp.apply_style(font_size=13)
classes = ["promoter", "genic", "distal"]
rows = {"tested": tested, f"up ({n_up:,})": res[res["regulation"] == "up"],
        f"down ({n_down:,})": res[res["regulation"] == "down"]}
share = pd.DataFrame({k: v["annotation"].value_counts(normalize=True).reindex(classes).fillna(0)
                      for k, v in rows.items()}).T
colours = dict(zip(classes, [sp.PALETTE["gold"], sp.PALETTE["teal"], sp.PALETTE["grey_dark"]]))
fig, ax = plt.subplots(figsize=(8, 3.2))
left = np.zeros(len(share))
y = np.arange(len(share))[::-1]
for c in classes:
    bars = ax.barh(y, share[c], left=left, color=colours[c], height=0.65, label=c)
    for yi, l, w in zip(y, left, share[c]):
        if w >= 0.06:
            ax.text(l + w / 2, yi, f"{w:.0%}", ha="center", va="center", fontsize=10, color="white")
    left += share[c].to_numpy()
ax.set_yticks(y)
ax.set_yticklabels(share.index)
ax.tick_params(axis="y", length=0)
ax.set_xlim(0, 1)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=3, borderaxespad=0, fontsize=10)
ax.set_title(title, loc="right", fontsize=12)
sp.save(fig, base + "peak_classes", formats)
