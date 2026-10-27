"""Gene-set enrichment for one contrast, two complementary ways:

ORA   the significant genes, up and down separately, against the genes that
      were tested (not the whole genome: a whole-genome background inflates
      every result, since unexpressed genes can never be called).
GSEA  every tested gene, ranked by the Wald statistic (or the shrunken log2
      fold change), so coordinated shifts too small to pass the cut-offs
      still count.

Each collection (Hallmark, Reactome, GO:BP, custom) is tested on its own, and
adjusted p-values / FDR are within a collection."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
sys.path.insert(0, str(Path(sm.scriptdir)))

from de_common import log_to, neg_log10, short_term  # noqa: E402
from smcore import plotting as sp  # noqa: E402

log_to(sm.log[0])
de = sm.params.de
gs_cfg, gsea_cfg = de["gene_sets"], de["gsea"]
test, reference = sm.wildcards.contrast.split("_vs_", 1)
split = sm.wildcards.split
formats = de["plot_formats"]
out_dir = Path(sm.output.gsea).parent

import gseapy as gp  # noqa: E402

res = pd.read_csv(sm.input.results, sep="\t", index_col="gene_id")
res = res.dropna(subset=["gene_name", "stat"])
rank_col = "stat" if gsea_cfg["rank_by"] == "stat" else "log2FoldChange_shrunken"
# Gene sets use symbols. A symbol on two gene IDs (PAR copies, readthroughs)
# keeps the copy with the strongest evidence.
res = res.reindex(res[rank_col].abs().sort_values(ascending=False).index).drop_duplicates("gene_name")
tested = res.dropna(subset=["padj"])
universe = set(tested["gene_name"])
up = tested.loc[tested["regulation"] == "up", "gene_name"].tolist()
down = tested.loc[tested["regulation"] == "down", "gene_name"].tolist()
ranking = res.set_index("gene_name")[rank_col].sort_values(ascending=False)
print(f"{len(universe)} genes tested, {len(up)} up, {len(down)} down; GSEA ranks {len(ranking)} genes by {rank_col}")


def read_gmt(path, genes):
    sets = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        members = sorted(set(fields[2:]) & genes)
        if gs_cfg["min_size"] <= len(members) <= gs_cfg["max_size"]:
            sets[fields[0]] = members
    return sets


ORA_COLUMNS = ["collection", "term", "overlap", "set_size", "pvalue", "padj", "odds_ratio", "genes"]
GSEA_COLUMNS = ["collection", "term", "size", "es", "nes", "pvalue", "fdr", "fwer", "leading_edge"]
ora = {"up": [], "down": []}
gsea, curves = [], {}

for collection, path in zip(sm.params.collections, sm.input.gene_sets):
    sets = read_gmt(path, universe)
    print(f"{collection}: {len(sets)} sets of {gs_cfg['min_size']}-{gs_cfg['max_size']} tested genes")
    if not sets:
        continue
    for direction, genes in (("up", up), ("down", down)):
        if len(genes) < 3:
            continue
        r = gp.enrich(gene_list=genes, gene_sets=sets, background=sorted(universe), outdir=None, cutoff=1.0,
                      no_plot=True, verbose=False).results
        if r is None or r.empty:
            continue
        overlap = r["Overlap"].str.split("/", expand=True).astype(int)
        ora[direction].append(pd.DataFrame({
            "collection": collection,
            "term": r["Term"],
            "overlap": overlap[0],
            "set_size": overlap[1],
            "pvalue": r["P-value"],
            "padj": r["Adjusted P-value"],
            "odds_ratio": r["Odds Ratio"],
            "genes": r["Genes"].str.replace(";", ","),
        }))
    pre = gp.prerank(rnk=ranking, gene_sets=sets, min_size=gs_cfg["min_size"], max_size=gs_cfg["max_size"],
                     permutation_num=gsea_cfg["permutations"], seed=gsea_cfg["seed"], threads=sm.threads,
                     outdir=None, no_plot=True, verbose=False)
    r = pre.res2d
    gsea.append(pd.DataFrame({
        "collection": collection,
        "term": r["Term"],
        "size": [len(pre.results[t]["matched_genes"].split(";")) for t in r["Term"]],
        "es": r["ES"].astype(float),
        "nes": r["NES"].astype(float),
        "pvalue": r["NOM p-val"].astype(float),
        "fdr": r["FDR q-val"].astype(float),
        "fwer": r["FWER p-val"].astype(float),
        "leading_edge": r["Lead_genes"].str.replace(";", ","),
    }))
    for t in r["Term"]:
        curves[(collection, t)] = (np.asarray(pre.results[t]["RES"]), np.asarray(pre.results[t]["hits"]))


def table(parts, columns, sort):
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    return df[columns].sort_values(sort) if len(df) else df


ora_tables = {d: table(ora[d], ORA_COLUMNS, ["padj", "pvalue"]) for d in ora}
ora_tables["up"].to_csv(sm.output.ora_up, sep="\t", index=False, float_format="%.6g")
ora_tables["down"].to_csv(sm.output.ora_down, sep="\t", index=False, float_format="%.6g")
gsea_table = table(gsea, GSEA_COLUMNS, ["fdr", "pvalue"])
gsea_table.to_csv(sm.output.gsea, sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- plots
import matplotlib.pyplot as plt  # noqa: E402

sp.apply_style(font_size=13)
where = "" if split == "all" else f" · {split}"
TOP = 12


def label(row):
    return f"{short_term(row.term)}" + ("" if len(sm.params.collections) == 1 else f"  [{row.collection}]")


def empty(ax, text):
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=sp.PALETTE["grey_dark"])
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)


# ORA: one panel per direction, bar length = -log10 adjusted p, overlap
# written at the bar's end; pale bars did not pass alpha.
fig, axes = plt.subplots(1, 2, figsize=(15, 0.42 * TOP + 2))
for ax, direction, colour, who in ((axes[0], "up", sp.UP, test), (axes[1], "down", sp.DOWN, reference)):
    top = ora_tables[direction].head(TOP).iloc[::-1]
    ax.set_title(f"Higher in {who} ({len(up if direction == 'up' else down):,} genes)", color=colour)
    if top.empty:
        empty(ax, "too few significant genes" if len(up if direction == "up" else down) < 3 else "no gene set")
        continue
    values = neg_log10(top["padj"])
    shades = [colour if p < de["alpha"] else sp.lighten(colour, 0.65) for p in top["padj"]]
    bars = ax.barh(range(len(top)), values, color=shades, height=0.72)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([label(r) for r in top.itertuples()], fontsize=10)
    ax.set_xlabel("$-$log$_{10}$ adjusted p")
    ax.set_xlim(0, values.max() * 1.25)
    ax.axvline(-np.log10(de["alpha"]), color=sp.PALETTE["grey_dark"], lw=1, ls="--", zorder=0)
    sp.label_bars(ax, bars, [f"{r.overlap}/{r.set_size}" for r in top.itertuples()])
    ax.tick_params(axis="y", length=0)
fig.suptitle(f"Over-representation · {test} vs {reference}{where}", x=0.01, ha="left", fontweight="bold")
sp.save(fig, out_dir / "ora", formats)

# GSEA: normalised enrichment scores of the strongest sets each way, FDR on
# the bar; pale bars are above the FDR cut-off.
fig, ax = plt.subplots(figsize=(9, 0.38 * 2 * TOP + 2))
if gsea_table.empty:
    empty(ax, "no gene set")
else:
    pos = gsea_table[gsea_table["nes"] > 0].sort_values(["fdr", "nes"], ascending=[True, False]).head(TOP)
    neg = gsea_table[gsea_table["nes"] < 0].sort_values(["fdr", "nes"], ascending=[True, True]).head(TOP)
    top = pd.concat([neg.sort_values("nes", ascending=False), pos.sort_values("nes")])
    shades = [(sp.UP if n > 0 else sp.DOWN) if q < gsea_cfg["fdr"] else sp.lighten(sp.UP if n > 0 else sp.DOWN, 0.65)
              for n, q in zip(top["nes"], top["fdr"])]
    bars = ax.barh(range(len(top)), top["nes"], color=shades, height=0.72)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([label(r) for r in top.itertuples()], fontsize=10)
    ax.axvline(0, color=sp.PALETTE["ink"], lw=1.2)
    limit = top["nes"].abs().max() * 1.35
    ax.set_xlim(-limit, limit)
    floor = 1 / gsea_cfg["permutations"]  # q = 0 means below what the permutations can resolve
    sp.label_bars(ax, bars, [f"q<{floor:g}" if q < floor else f"q={q:.2g}" for q in top["fdr"]])
    ax.set_xlabel(f"NES  (← {reference}  |  {test} →)")
    ax.tick_params(axis="y", length=0)
ax.set_title(f"GSEA · {test} vs {reference}{where}")
sp.save(fig, out_dir / "gsea", formats)

# Running enrichment score for the top three sets each way.
picked = []
if not gsea_table.empty:
    for sign in (1, -1):
        part = gsea_table[np.sign(gsea_table["nes"]) == sign].sort_values(["fdr", "pvalue"]).head(3)
        picked += [(r.collection, r.term, r.nes, r.fdr) for r in part.itertuples()]
fig, axes = plt.subplots(2, 3, figsize=(15, 8), squeeze=False)
for ax, item in zip(axes.ravel(), picked + [None] * (6 - len(picked))):
    if item is None:
        ax.set_visible(False)
        continue
    collection, term, nes, fdr = item
    res_curve, hits = curves[(collection, term)]
    colour = sp.UP if nes > 0 else sp.DOWN
    x = np.arange(len(res_curve))
    ax.plot(x, res_curve, color=colour, lw=2.5)
    ax.axhline(0, color=sp.PALETTE["ink"], lw=1)
    low = min(res_curve.min(), 0)
    span = res_curve.max() - low
    ax.vlines(hits, low - 0.16 * span, low - 0.04 * span, color=sp.PALETTE["ink"], lw=0.6)
    ax.set_ylim(low - 0.18 * span, res_curve.max() + 0.08 * span)
    ax.set_xlim(0, len(res_curve))
    floor = 1 / gsea_cfg["permutations"]
    q = f"q < {floor:g}" if fdr < floor else f"q = {fdr:.2g}"
    ax.set_title(f"{short_term(term, 40)}\nNES {nes:.2f} \u00b7 {q}", fontsize=12)
    ax.set_xlabel(f"rank ({test} → {reference})", fontsize=11)
    ax.set_ylabel("running ES", fontsize=11)
if not picked:
    axes[0, 0].set_visible(True)
    empty(axes[0, 0], "no gene set")
fig.suptitle(f"GSEA running scores · {test} vs {reference}{where}", x=0.01, ha="left", fontweight="bold")
sp.save(fig, out_dir / "gsea_curves", formats)
