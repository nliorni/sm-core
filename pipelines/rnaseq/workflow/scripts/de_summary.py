"""Every analysis, split and contrast of one count matrix in a table and one
figure: how many genes go up and down, and how many gene sets GSEA calls."""

import sys
from pathlib import Path

import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
from smcore import de as sde  # noqa: E402

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

sde.summary_bars(summary, str(Path(sm.output.tsv).with_suffix("")), de["plot_formats"],
                 xlabel=f"genes  (padj < {de['alpha']}, |log$_2$FC| ≥ {de['min_abs_log2fc']})",
                 title=f"Differential expression · {sm.wildcards.branch}")
