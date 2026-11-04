"""featureCounts' table as the consensus-peak count matrix: one row per peak
(0-based BED coordinates, as consensus_peaks.bed), one column per sample; and
the sample metadata table next to it."""

import sys

import pandas as pd

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")

raw = pd.read_csv(sm.input.raw, sep="\t", comment="#")
raw = raw.rename(columns=dict(zip(sm.params.bams, sm.params.samples)))
counts = pd.DataFrame({
    "peak_id": raw["Geneid"],
    "chrom": raw["Chr"],
    "start": raw["Start"] - 1,
    "end": raw["End"],
})
counts = pd.concat([counts, raw[sm.params.samples]], axis=1)
counts.to_csv(sm.output.counts, sep="\t", index=False)

meta = pd.DataFrame.from_dict(sm.params.meta, orient="index")
meta.index.name = "sample"
meta.reindex(sm.params.samples).to_csv(sm.output.samples, sep="\t")

total = counts[sm.params.samples].sum()
print(f"{len(counts)} peaks x {len(sm.params.samples)} samples; fragments in peaks per sample: "
      + ", ".join(f"{s} {n:,}" for s, n in total.items()))
