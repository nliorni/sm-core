"""One model's samples at a glance: blind VST, normalised counts, size factors,
PCA and sample-to-sample distances. The first thing to look at before any
contrast: replicates should sit together, and an outlier or a swapped sample
shows here."""

import sys
from pathlib import Path

import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from de_common import log_to, model_samples, read_counts, vst  # noqa: E402
from smcore import de as sde  # noqa: E402

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

base = sm.output.vst.rsplit("/", 1)[0] + "/"
title = "all samples" if split == "all" else split
sde.pca(expr, meta, model, base + "pca", formats, title=f"PCA \u00b7 {title}")
sde.sample_distances(expr, meta, model, base + "sample_distances", formats,
                     title=f"Sample distances \u00b7 {title}")
