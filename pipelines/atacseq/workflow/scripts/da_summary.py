"""Every analysis, split and contrast in one table and one figure: how many
peaks open and close, and how many of those are at promoters."""

import sys
from pathlib import Path

import pandas as pd

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.smcore)
from smcore import de as sde  # noqa: E402

da = sm.params.da
base = Path(sm.params.da_dir)
rows = []
for analysis, split, contrast in sm.params.units:
    res = pd.read_csv(base / analysis / split / contrast / "da_results.tsv", sep="\t",
                      usecols=["padj", "regulation", "annotation"])
    up, down = res["regulation"] == "up", res["regulation"] == "down"
    promoter = res["annotation"] == "promoter"
    rows.append({
        "analysis": analysis, "split": split, "contrast": contrast,
        "peaks_tested": int(res["padj"].notna().sum()),
        "up": int(up.sum()), "down": int(down.sum()),
        "up_promoter": int((up & promoter).sum()), "down_promoter": int((down & promoter).sum()),
    })
summary = pd.DataFrame(rows)
summary.to_csv(sm.output.tsv, sep="\t", index=False)

sde.summary_bars(summary, str(Path(sm.output.tsv).with_suffix("")), da["plot_formats"],
                 xlabel=f"peaks  (padj < {da['alpha']}, |log$_2$FC| ≥ {da['min_abs_log2fc']})",
                 title=f"Differential accessibility · {da['normalization'].replace('_', ' ')}")
