"""One contrast of one model: PyDESeq2 Wald test, shrunken log2 fold changes,
significant up/down lists, and the volcano, MA and top-gene heatmap."""

import sys
from pathlib import Path

import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from de_common import (  # noqa: E402
    fit, log_to, model_samples, prefilter, read_counts, regulation, vst, with_reference,
)
from smcore import de as sde  # noqa: E402

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

res["regulation"] = regulation(res, de["alpha"], de["min_abs_log2fc"])
res = genes.join(res, how="inner").sort_values(["padj", "pvalue"], na_position="last")
res.index.name = "gene_id"
res.to_csv(sm.output.results, sep="\t", float_format="%.6g")
res[res["regulation"] == "up"].to_csv(sm.output.up, sep="\t", float_format="%.6g")
res[res["regulation"] == "down"].to_csv(sm.output.down, sep="\t", float_format="%.6g")
n_up, n_down = int((res["regulation"] == "up").sum()), int((res["regulation"] == "down").sum())
print(f"{test} vs {reference}: {n_up} up, {n_down} down (padj < {de['alpha']}, |log2FC| >= {de['min_abs_log2fc']})")

# ---------------------------------------------------------------- plots
title = f"{test} vs {reference}" + ("" if split == "all" else f" · {split}")
base = sm.output.results.rsplit("/", 1)[0] + "/"

sde.volcano(res, base + "volcano", formats, label_col="gene_name", alpha=de["alpha"],
            min_abs_log2fc=de["min_abs_log2fc"], label_top=de["label_top"], title=title,
            xlabel=f"log$_2$ fold change ({test} / {reference})")

# MA: mean expression against the shrunken fold change, the estimate that is
# not inflated for weakly expressed genes.
lfc = "log2FoldChange_shrunken" if de["shrink"] else "log2FoldChange"
sde.ma(res, base + "ma", formats, lfc_col=lfc, title=title,
       xlabel="log$_{10}$ mean normalised count",
       ylabel="log$_2$ fold change" + (" (shrunken)" if de["shrink"] else ""))

# Heatmap: the most significant genes across every sample of the model,
# z-scored variance-stabilised expression, samples grouped by the factor.
expr, _ = vst(counts, meta, sm.params.model["design"], sm.threads)
picked = res[res["regulation"] != "ns"].head(de["heatmap_top"])
if len(picked) < 2:
    picked = res.dropna(subset=["padj"]).head(de["heatmap_top"])
levels = [reference, test] + [l for l in meta[factor].cat.categories if l not in (reference, test)]
sde.heatmap(expr.loc[picked.index], picked["gene_name"].fillna("").tolist(), meta, factor, levels,
            base + "heatmap", formats, title=f"Top {len(picked)} genes \u00b7 {title}")
