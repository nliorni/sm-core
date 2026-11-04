"""Best genetic match of every sample, from somalier's pairwise relatedness.

One row per sample (ATAC, and other assays' BAMs from qc.identity_bams): its
most related ATAC sample and its most related other-assay sample, with the
relatedness and the number of sites compared, and `match` / `match_other`:
the best sample when the relatedness reaches MATCH_MIN (the same cells or the
same individual give ~1, unrelated ones ~0 or below), else empty.

The question it answers without relying on sample names: which RNA library
is this ATAC library's cell line? A sample whose best match is not the one
its label says is a swap.
"""

import sys

import pandas as pd

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")

MATCH_MIN = 0.8   # relatedness of the same cells or individual
LOW_SITES = 100   # fewer shared informative sites: a noisy estimate

prefix = sm.params.prefix
has_other = sm.params.has_other
pairs = pd.read_csv(sm.input.pairs, sep="\t").rename(columns={"#sample_a": "sample_a"})
both = pd.concat([
    pairs[["sample_a", "sample_b", "relatedness", "n"]],
    pairs.rename(columns={"sample_a": "sample_b", "sample_b": "sample_a"})[["sample_a", "sample_b", "relatedness", "n"]],
])


def assay(name):
    return "other" if has_other and name.startswith(prefix) else "atac"


both["assay_a"] = both["sample_a"].map(assay)
both["assay_b"] = both["sample_b"].map(assay)

rows = []
for sample, sub in both.groupby("sample_a", sort=False):
    row = {"sample": sample, "assay": assay(sample)}
    for kind, suffix in (("atac", ""), ("other", "_other")):
        if kind == "other" and not has_other:
            continue
        cand = sub[sub["assay_b"] == kind]
        if cand.empty:
            row.update({f"best{suffix}": "", f"relatedness{suffix}": float("nan"),
                        f"sites{suffix}": 0, f"match{suffix}": ""})
            continue
        best = cand.loc[cand["relatedness"].idxmax()]
        row.update({
            f"best{suffix}": best["sample_b"],
            f"relatedness{suffix}": round(float(best["relatedness"]), 3),
            f"sites{suffix}": int(best["n"]),
            f"match{suffix}": best["sample_b"] if best["relatedness"] >= MATCH_MIN else "",
        })
    rows.append(row)

out = pd.DataFrame(rows)
if out.empty:  # a single sample: nothing to compare
    out = pd.DataFrame(columns=["sample", "assay", "best", "relatedness", "sites", "match"])
    out.to_csv(sm.output[0], sep="\t", index=False)
    sys.exit(0)
out["assay"] = pd.Categorical(out["assay"], ["atac", "other"])
out = out.sort_values(["assay", "sample"])
out.to_csv(sm.output[0], sep="\t", index=False)
print(out.to_string(index=False))
low = [c for c in ("sites", "sites_other") if c in out]
if any((out[c] < LOW_SITES).any() for c in low):
    print(f"\nNote: some best matches rest on fewer than {LOW_SITES} sites; treat them with care.")
