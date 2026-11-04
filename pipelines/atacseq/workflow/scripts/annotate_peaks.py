"""Annotate consensus peaks with genes.

For each peak:

  promoter_gene_ids /     every gene (of the configured types) whose promoter
  promoter_gene_names     window [TSS - promoter_upstream, TSS + promoter_downstream],
                          in the gene's direction, the peak overlaps; comma-separated,
                          best first. A merged peak can cover the promoters of
                          several genes (bidirectional promoters, readthrough
                          and antisense genes), and each of them gets it.
  nearest_gene_*          one gene per peak: the best promoter gene if there
                          is one, else the gene whose TSS is nearest
  distance_to_tss         signed distance from that gene's TSS to the peak: 0
                          when the peak covers it, negative upstream, positive
                          downstream (gene's direction)
  annotation              promoter (any promoter gene), genic (inside the body
                          of any gene of those types), distal (neither)

"Best" among genes at the same distance: the order of annotation.gene_types
(protein_coding first by default), then a named gene over one known only by
its ID (GENCODE's new lncRNA and readthrough genes), then the gene ID. TSS:
the gene's canonical transcript (gtf_annotation.py).
"""

import sys

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")
cfg = sm.params.annotation
up, down = cfg["promoter_upstream"], cfg["promoter_downstream"]
reach = max(up, down)

peaks = pd.read_csv(sm.input.bed, sep="\t", header=None,
                    names=["chrom", "start", "end", "peak_id", "n_samples", "strand"])
genes = pd.read_csv(sm.input.genes, sep="\t", dtype={"chrom": str, "gene_name": str})
genes = genes[genes["gene_type"].isin(cfg["gene_types"])].reset_index(drop=True)
if genes.empty:
    sys.exit(f"no gene of types {cfg['gene_types']} in the annotation (annotation.gene_types).")
print(f"{len(peaks)} peaks; {len(genes)} genes of types {', '.join(cfg['gene_types'])}")

type_rank = {t: i for i, t in enumerate(cfg["gene_types"])}
gene_name = genes["gene_name"].fillna("").astype(str)
unversioned = genes["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
unnamed = (gene_name == "") | (gene_name == genes["gene_id"]) | (gene_name == unversioned)
# Sort key for ties, smaller = better.
genes["_rank"] = genes["gene_type"].map(type_rank) * 2 + unnamed.astype(int)


def gene_bodies(g):
    """Merged gene bodies of one chromosome, as sorted start / end arrays."""
    merged = []
    for a, b in g.sort_values("start")[["start", "end"]].to_numpy():
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return np.array(merged).T


rows = {}
for chrom, p in peaks.groupby("chrom", sort=False):
    g = genes[genes["chrom"] == chrom].sort_values("tss")
    if g.empty:
        for i in p.index:
            rows[i] = (-1, np.nan, False, [])
        continue
    tss, strand = g["tss"].to_numpy(), g["strand"].to_numpy()
    rank, gid, idx = g["_rank"].to_numpy(), g["gene_id"].to_numpy(), g.index.to_numpy()
    ms, me = gene_bodies(g)
    for i, s, e in zip(p.index, p["start"].to_numpy(), p["end"].to_numpy()):
        # Candidates: every TSS within `reach` of the peak, or else the two
        # flanking TSSs (the nearest is one of them).
        lo, hi = np.searchsorted(tss, s - reach), np.searchsorted(tss, e + reach)
        if lo == hi:
            lo, hi = max(lo - 1, 0), min(hi + 1, len(tss))
        k = np.arange(lo, hi)
        t = tss[k]
        dist = np.where(t < s, s - t, np.where(t >= e, t - e + 1, 0))
        # Signed: where the peak lies relative to the TSS, in the gene's direction.
        signed = dist * np.where(t < s, 1, np.where(t >= e, -1, 0)) * np.where(strand[k] == "+", 1, -1)
        in_promoter = (signed >= -up) & (signed <= down)
        pool = k[in_promoter] if in_promoter.any() else k
        pool_d = signed[in_promoter] if in_promoter.any() else signed
        order = np.lexsort((gid[pool], rank[pool], np.abs(pool_d)))
        best = pool[order[0]]
        promoter_genes = [idx[j] for j in pool[order]] if in_promoter.any() else []
        j = np.searchsorted(ms, e) - 1  # last merged body starting before the peak ends
        genic = j >= 0 and me[j] > s
        rows[i] = (idx[best], int(pool_d[order[0]]), genic, promoter_genes)

ann = pd.DataFrame.from_dict(rows, orient="index", columns=["gene_index", "distance", "genic", "promoter"])
ann = ann.loc[peaks.index]
has_gene = ann["gene_index"] >= 0
hit = genes.reindex(ann.loc[has_gene, "gene_index"])

res = peaks[["peak_id", "chrom", "start", "end"]].copy()
res["width"] = res["end"] - res["start"]
res["n_samples"] = peaks["n_samples"]
for col, name in (("gene_id", "nearest_gene_id"), ("gene_name", "nearest_gene_name"),
                  ("gene_type", "nearest_gene_type"), ("strand", "nearest_gene_strand")):
    res[name] = ""
    res.loc[has_gene, name] = hit[col].to_numpy()
res["distance_to_tss"] = ann["distance"].astype("Int64")
is_promoter = ann["promoter"].map(len) > 0
res["annotation"] = np.where(is_promoter, "promoter", np.where(ann["genic"], "genic", "distal"))
res["promoter_gene_ids"] = ann["promoter"].map(lambda js: ",".join(genes.loc[js, "gene_id"]))
res["promoter_gene_names"] = ann["promoter"].map(lambda js: ",".join(gene_name[js]))
res.to_csv(sm.output[0], sep="\t", index=False)

print(res["annotation"].value_counts().to_string())
multi = int((ann["promoter"].map(len) > 1).sum())
print(f"{multi} promoter peaks cover the promoters of more than one gene; "
      f"{len({j for js in ann['promoter'] for j in js})} genes have a promoter peak")
