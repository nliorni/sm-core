"""Every sample's QC against the ENCODE ATAC-seq data standards: a table
(results/qc/atac_qc_summary.tsv), three figures, and the same as MultiQC
custom content.

ENCODE's tiers (https://www.encodeproject.org/atac-seq/), ideal / acceptable
(a value on a boundary counts as the higher tier):
  alignment rate      95% / 80%
  final fragments     25 M non-duplicate, non-mitochondrial pairs
  NRF                 0.9 / 0.7
  PBC1                0.9 / 0.7
  PBC2                3   / 1
  FRiP                0.3 / 0.2
  TSS enrichment      7   / 5   (GRCh38, ENCODE's TSS set; ours is the GTF's
                                 canonical protein-coding TSSs)
Anything below the acceptable tier is "low". The mitochondrial and
duplicate fractions have no ENCODE threshold and are reported as they are.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")
sys.path.insert(0, sm.params.smcore)
from smcore import plotting as sp  # noqa: E402

samples = sm.params.samples
group_of = {s: g for g, members in sm.params.groups.items() for s in members}

# metric: (ideal, acceptable or None, label, format)
STANDARDS = {
    "aligned_rate": (0.95, 0.80, "alignment rate", "{:.1%}"),
    "fragments_final": (25e6, None, "final fragments", "{:,.0f}"),
    "NRF": (0.9, 0.7, "NRF", "{:.3f}"),
    "PBC1": (0.9, 0.7, "PBC1", "{:.3f}"),
    "PBC2": (3.0, 1.0, "PBC2", "{:.2f}"),
    "FRiP": (0.3, 0.2, "FRiP", "{:.3f}"),
    "tss_enrichment": (7.0, 5.0, "TSS enrichment", "{:.2f}"),
}


def status(value, ideal, acceptable):
    if value >= ideal:
        return "ideal"
    if acceptable is not None and value >= acceptable:
        return "acceptable"
    return "low"


def bowtie2_rate(path):
    text = Path(path).read_text()
    m = re.search(r"([\d.]+)% overall alignment rate", text)
    return float(m.group(1)) / 100 if m else np.nan


def picard_duplication(path):
    lines = Path(path).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("## METRICS CLASS")) + 1
    header = lines[start].split("\t")
    rows = []
    for line in lines[start + 1:]:
        if not line.strip():
            break
        rows.append(dict(zip(header, line.split("\t"))))
    dup = sum(float(r["READ_PAIR_DUPLICATES"]) for r in rows)
    pairs = sum(float(r["READ_PAIRS_EXAMINED"]) for r in rows)
    return dup / pairs if pairs else np.nan


def frip(path):
    table = pd.read_csv(path, sep="\t", index_col=0)
    col = table.columns[0]
    total = table[col].sum()
    return table.loc["Assigned", col] / total if total else np.nan


rows = []
for i, s in enumerate(samples):
    m = json.loads(Path(sm.input.metrics[i]).read_text())
    with open(sm.input.peaks[i]) as handle:
        n_peaks = sum(1 for _ in handle)
    rows.append({
        "sample": s,
        "group": group_of.get(s, ""),
        "aligned_rate": bowtie2_rate(sm.input.bowtie2[i]),
        "fragments_total": m["fragments_total"],
        "mito_fraction": m["mito_fraction"],
        "duplicate_fraction": picard_duplication(sm.input.markdup[i]),
        "fragments_final": m["fragments_final"],
        "NRF": m["NRF"],
        "PBC1": m["PBC1"],
        "PBC2": m["PBC2"],
        "FRiP": frip(sm.input.frip[i]),
        "peaks": n_peaks,
        "tss_enrichment": m["tss_enrichment"],
        "nucleosome_free_fraction": m["nucleosome_free_fraction"],
        "mononucleosome_fraction": m["mononucleosome_fraction"],
        "median_fragment_length": m["median_fragment_length"],
        "chrX_per_million": m.get("chrX_per_million", np.nan),
        "chrY_per_million": m.get("chrY_per_million", np.nan),
    })
qc = pd.DataFrame(rows).set_index("sample")
for metric, (ideal, acceptable, _, _) in STANDARDS.items():
    qc[f"{metric}_status"] = [status(v, ideal, acceptable) for v in qc[metric]]
status_cols = [f"{m}_status" for m in STANDARDS]
qc["below_acceptable"] = [
    ",".join(m for m in STANDARDS if qc.loc[s, f"{m}_status"] == "low") for s in qc.index
]
qc.to_csv(sm.output.table, sep="\t", float_format="%.6g")
print(qc[[c for c in qc.columns if not c.endswith("_status")]].to_string())

# ---------------------------------------------------------------- MultiQC
mqc_cols = {
    "aligned_rate": ("Aligned", "{:.1%}"),
    "mito_fraction": ("Mito", "{:.1%}"),
    "fragments_final": ("Final frags (M)", None),
    "NRF": ("NRF", None), "PBC1": ("PBC1", None), "PBC2": ("PBC2", None),
    "FRiP": ("FRiP", None), "peaks": ("Peaks", None),
    "tss_enrichment": ("TSS enr.", None),
    "nucleosome_free_fraction": ("NFR frac.", None),
    "chrY_per_million": ("chrY / M", None),
}
with open(sm.output.mqc_table, "w") as out:
    out.write("# id: 'atac_qc'\n# section_name: 'ATAC-seq QC'\n"
              "# description: 'Per-sample metrics against the ENCODE ATAC-seq standards "
              "(ideal / acceptable): alignment > 95% / 80%, final fragments >= 25 M, NRF > 0.9 / 0.7, "
              "PBC1 > 0.9 / 0.7, PBC2 > 3 / 1, FRiP > 0.3 / 0.2, TSS enrichment > 7 / 5. "
              "The column below_acceptable lists the metrics under the acceptable tier.'\n"
              "# plot_type: 'table'\n")
    table = qc[list(mqc_cols)].copy()
    table["fragments_final"] = table["fragments_final"] / 1e6
    table.columns = [mqc_cols[c][0] for c in mqc_cols]
    table["below_acceptable"] = qc["below_acceptable"].replace("", "-")
    table.index.name = "Sample"
    table.to_csv(out, sep="\t", float_format="%.4g")

fragments = {s: pd.read_csv(p, sep="\t") for s, p in zip(samples, sm.input.fragments)}
tss = {s: pd.read_csv(p, sep="\t") for s, p in zip(samples, sm.input.tss)}


def linegraph(id_, name, description, xlab, ylab, data):
    return {
        "id": id_, "section_name": name, "description": description, "plot_type": "linegraph",
        "pconfig": {"id": id_ + "_plot", "title": f"sm-core/atacseq: {name}", "xlab": xlab, "ylab": ylab},
        "data": data,
    }


frag_data = {}
for s, f in fragments.items():
    f = f[(f["length"] > 0) & (f["length"] <= 1000)]
    total = max(f["count"].sum(), 1)
    frag_data[s] = {int(l): round(c / total * 100, 5) for l, c in zip(f["length"], f["count"])}
Path(sm.output.mqc_fragments).write_text(json.dumps(linegraph(
    "atac_fragment_sizes", "Fragment sizes",
    "Fragment length distribution of the final BAM. A good library shows a nucleosome ladder: a "
    "nucleosome-free peak below 147 bp, then mono- (~200 bp) and di-nucleosome (~400 bp) peaks.",
    "fragment length (bp)", "% of fragments", frag_data)))
tss_data = {s: {int(p): round(v, 4) for p, v in zip(t["position"][::10], t["smoothed"][::10])}
            for s, t in tss.items()}
Path(sm.output.mqc_tss).write_text(json.dumps(linegraph(
    "atac_tss_enrichment", "TSS enrichment",
    "Tn5 cut sites around transcription start sites, relative to the background 1.9-2 kb away. "
    "The peak height is the TSS enrichment score.",
    "distance to TSS (bp)", "enrichment over background", tss_data)))

# ---------------------------------------------------------------- figures
import matplotlib.pyplot as plt  # noqa: E402

sp.apply_style(font_size=12)
groups = list(dict.fromkeys(group_of.get(s, "") for s in samples))
group_color = sp.colors_for(groups)
STATUS_COLOR = {"ideal": sp.PALETTE["green"], "acceptable": sp.PALETTE["gold"], "low": sp.PALETTE["red"]}
panels = list(STANDARDS) + ["mito_fraction"]
n = len(samples)
fig, axes = plt.subplots(2, 4, figsize=(18, max(6, 0.5 * n + 3) * 2 / 1.6), sharey=True)
y = np.arange(n)[::-1]
for ax, metric in zip(axes.ravel(), panels):
    values = qc[metric].replace([np.inf], np.nan)
    if metric == "fragments_final":
        values = values / 1e6
    if metric in STANDARDS:
        ideal, acceptable, label, _ = STANDARDS[metric]
        colors = [STATUS_COLOR[v] for v in qc[f"{metric}_status"]]
        scale = 1e6 if metric == "fragments_final" else 1
        for level, ls in ((ideal, "--"), (acceptable, ":")):
            if level is not None:
                ax.axvline(level / scale, color=sp.PALETTE["grey_dark"], lw=1.2, ls=ls, zorder=0)
    else:
        label, colors = "mitochondrial fraction", [sp.PALETTE["grey_dark"]] * n
    ax.barh(y, values.fillna(0), color=colors, height=0.7)
    ax.set_title(label + (" (M)" if metric == "fragments_final" else ""), fontsize=12)
    ax.tick_params(axis="y", length=0)
    if metric == "PBC2" and values.max() > 20:
        ax.set_xscale("log")
axes[0, 0].set_yticks(y)
axes[0, 0].set_yticklabels(samples, fontsize=9)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in STATUS_COLOR.values()]
handles += [plt.Line2D([], [], color=sp.PALETTE["grey_dark"], ls=ls) for ls in ("--", ":")]
fig.legend(handles, list(STATUS_COLOR) + ["ENCODE ideal", "ENCODE acceptable"], loc="upper right",
           ncol=5, frameon=False, fontsize=10)
fig.suptitle("ATAC-seq QC vs ENCODE standards", x=0.01, ha="left", fontweight="bold")
sp.save(fig, Path(sm.output.fig_metrics).with_suffix(""))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for s, f in fragments.items():
    f = f[(f["length"] > 0) & (f["length"] <= 1000)]
    frac = f["count"] / max(f["count"].sum(), 1)
    for ax in axes:
        ax.plot(f["length"], frac, color=group_color[group_of.get(s, "")], lw=1.2, alpha=0.8)
axes[1].set_yscale("log")
for ax, title in zip(axes, ("linear", "log scale")):
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fraction of fragments")
    ax.set_title(f"Fragment sizes ({title})")
    for x0 in (147, 294):
        ax.axvline(x0, color=sp.PALETTE["grey"], lw=1, ls="--", zorder=0)
if len(groups) > 1:
    axes[1].legend([plt.Line2D([], [], color=group_color[g], lw=2) for g in groups], groups,
                   loc="upper right", fontsize=9)
sp.save(fig, Path(sm.output.fig_fragments).with_suffix(""))

fig, ax = plt.subplots(figsize=(7.5, 5.5))
for s, t in tss.items():
    ax.plot(t["position"], t["smoothed"], color=group_color[group_of.get(s, "")], lw=1.2, alpha=0.8)
ax.axhline(1, color=sp.PALETTE["grey_dark"], lw=1, ls=":")
ax.set_xlabel("distance to TSS (bp)")
ax.set_ylabel("Tn5 cuts / background")
ax.set_title("TSS enrichment")
if len(groups) > 1:
    ax.legend([plt.Line2D([], [], color=group_color[g], lw=2) for g in groups], groups,
              loc="upper right", fontsize=9)
sp.save(fig, Path(sm.output.fig_tss).with_suffix(""))
