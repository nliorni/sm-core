"""rnaseq's side of the differential expression scripts: reading a branch's
gene counts. The models, statistics and figures are shared with atacseq and
live in smcore/de.py (re-exported here, so the de_* scripts import from one
place). The scripts put the sm-core root on sys.path before importing this."""

from __future__ import annotations

import pandas as pd

from smcore.de import (  # noqa: F401  (re-exported)
    bh, design_factors, fit, log_to, model_samples, neg_log10, prefilter, regulation, short_term, vst,
    with_reference,
)


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
