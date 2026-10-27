"""Every analysis, split and contrast of one count matrix in a table and one
figure: how many genes go up and down, and how many gene sets GSEA calls."""

import sys
from pathlib import Path

import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
from smcore import plotting as sp  # noqa: E402

de = sm.params.de
base = Path(sm.params.de_dir) / sm.wildcards.branch
rows = []
for analysis, split, contrast in sm.params.units:
    d = base / analysis / split / contrast
    res = pd.read_csv(d / "deseq2_results.tsv", sep="\t", usecols=["padj", "regulation"])
    gsea = pd.read_csv(d / "enrichment" / "gsea.tsv", sep="\t", usecols=["nes", "fdr"])
    hit = gsea["fdr"] < de["gsea"]["fdr"]
    rows.append({
        "analysis": analysis, "split": split, "contrast": contrast,
        "genes_tested": int(res["padj"].notna().sum()),
        "up": int((res["regulation"] == "up").sum()),
        "down": int((res["regulation"] == "down").sum()),
        "gsea_sets_up": int((hit & (gsea["nes"] > 0)).sum()),
        "gsea_sets_down": int((hit & (gsea["nes"] < 0)).sum()),
    })
summary = pd.DataFrame(rows)
summary.to_csv(sm.output.tsv, sep="\t", index=False)

import matplotlib.pyplot as plt  # noqa: E402

sp.apply_style(font_size=13)
labels = [
    c.replace("_vs_", " vs ") + ("" if s == "all" else f"  · {s}") + ("" if len(set(summary["analysis"])) == 1 else f"  [{a}]")
    for a, s, c in zip(summary["analysis"], summary["split"], summary["contrast"])
]
n = len(summary)
fig, ax = plt.subplots(figsize=(10, 0.45 * n + 2))
y = range(n)[::-1]
up = ax.barh(list(y), summary["up"], color=sp.UP, height=0.7, label="up")
down = ax.barh(list(y), -summary["down"], color=sp.DOWN, height=0.7, label="down")
limit = max(summary["up"].max(), summary["down"].max(), 1) * 1.3
ax.set_xlim(-limit, limit)
sp.label_bars(ax, up, [f"{v:,}" for v in summary["up"]])
sp.label_bars(ax, down, [f"{v:,}" for v in summary["down"]])
ax.axvline(0, color=sp.PALETTE["ink"], lw=1.2)
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=11)
ax.tick_params(axis="y", length=0)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(int(v)):,}"))
ax.set_xlabel(f"genes  (padj < {de['alpha']}, |log$_2$FC| ≥ {de['min_abs_log2fc']})")
ax.legend(loc="lower right", bbox_to_anchor=(1, 1.01), ncol=2, borderaxespad=0)
ax.set_title(f"Differential expression · {sm.wildcards.branch}")
sp.save(fig, Path(sm.output.tsv).with_suffix(""), de["plot_formats"])
