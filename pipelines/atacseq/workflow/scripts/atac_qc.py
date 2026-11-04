"""Per-sample ATAC-seq QC from the duplicate-marked BAM, in one pass.

The marked BAM still holds duplicates and mitochondrial pairs, which two of
these metrics need. Each properly paired fragment is read once, from its
first mate:

library complexity (ENCODE)
    over nuclear fragments, duplicates included, keyed by position (start,
    length): NRF = distinct / total, PBC1 = positions seen once / distinct,
    PBC2 = positions seen once / positions seen twice
mitochondrial fraction
    mitochondrial fragments / all fragments
fragment sizes
    the length distribution of the fragments kept in the final BAM (nuclear,
    and without duplicates when align.remove_duplicates): the nucleosome
    ladder, with fractions below 147 bp (nucleosome-free) and 147-294 bp
    (mono-nucleosome)
TSS enrichment
    Tn5 cut sites (fragment ends shifted +4 / -5 bp) aggregated per base over
    +-2 kb around every TSS, in the gene's direction; divided by the mean of
    the outer 100 bp at each end, smoothed over 51 bp; the score is the
    maximum within 250 bp of the TSS
sex chromosomes
    fragments on X and Y per million nuclear fragments
"""

import json
import sys
from array import array

import numpy as np
import pysam

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")

WINDOW, FLANK, SMOOTH, CENTRE = 2000, 100, 51, 250
MAX_LEN = 2000  # bowtie2 -X default; longer ones go in the last bin
remove_dups = sm.params.remove_duplicates
mito = open(sm.input.mito).read().strip()

bam = pysam.AlignmentFile(sm.input.bam)
names = list(bam.references)
mito_id = names.index(mito) if mito in names else -1
# X and Y under either naming, reported as chrX / chrY.
sex_ids = {i: "chr" + n.removeprefix("chr") for i, n in enumerate(names) if n in ("chrX", "X", "chrY", "Y")}

keys = [array("q") for _ in names]   # complexity: start * 4096 + length
cuts = [array("q") for _ in names]   # Tn5 cut sites of kept fragments
lengths = np.zeros(MAX_LEN + 1, dtype=np.int64)
n_all = n_mito = n_kept = 0
sex = {"chrX": 0, "chrY": 0}

for r in bam.fetch(until_eof=True):
    flag = r.flag
    # First mate of a proper pair, primary: one record per fragment.
    if flag & 0x942 != 0x42:
        continue
    n_all += 1
    rid = r.reference_id
    if rid == mito_id:
        n_mito += 1
        continue
    length = abs(r.template_length)
    start = min(r.reference_start, r.next_reference_start)
    keys[rid].append(start * 4096 + min(length, 4095))
    if remove_dups and flag & 0x400:
        continue
    n_kept += 1
    lengths[min(length, MAX_LEN)] += 1
    cuts[rid].append(start + 4)
    cuts[rid].append(start + length - 5)
    if rid in sex_ids:
        sex[sex_ids[rid]] += 1

# ---------------------------------------------------------------- complexity
distinct = once = twice = total = 0
for k in keys:
    if not len(k):
        continue
    _, counts = np.unique(np.frombuffer(k, dtype=np.int64), return_counts=True)
    total += int(counts.sum())
    distinct += len(counts)
    once += int((counts == 1).sum())
    twice += int((counts == 2).sum())
del keys

# ---------------------------------------------------------------- TSS enrichment
tss = {}
with open(sm.input.tss) as handle:
    for line in handle:
        f = line.split("\t")
        tss.setdefault(f[0], ([], []))
        tss[f[0]][0].append(int(f[1]))
        tss[f[0]][1].append(f[5].strip() == "-")
profile = np.zeros(2 * WINDOW + 1, dtype=np.int64)
n_tss = 0
for rid, name in enumerate(names):
    if name not in tss or not len(cuts[rid]):
        continue
    c = np.sort(np.frombuffer(cuts[rid], dtype=np.int64))
    pos = np.asarray(tss[name][0])
    minus = np.asarray(tss[name][1])
    lo = np.searchsorted(c, pos - WINDOW)
    hi = np.searchsorted(c, pos + WINDOW, side="right")
    n_tss += len(pos)
    rel = [np.where(m, p - c[a:b], c[a:b] - p) for p, m, a, b in zip(pos, minus, lo, hi) if b > a]
    if rel:
        profile += np.bincount(np.concatenate(rel) + WINDOW, minlength=2 * WINDOW + 1)
background = np.concatenate([profile[:FLANK], profile[-FLANK:]]).mean()
normalised = profile / background if background > 0 else np.zeros(len(profile))
smoothed = np.convolve(normalised, np.ones(SMOOTH) / SMOOTH, mode="same")
tss_score = float(smoothed[WINDOW - CENTRE:WINDOW + CENTRE + 1].max())

# ---------------------------------------------------------------- outputs
kept = max(n_kept, 1)
nuclear = n_all - n_mito
sizes = np.arange(MAX_LEN + 1)
metrics = {
    "sample": sm.wildcards.sample,
    "fragments_total": n_all,
    "fragments_mito": n_mito,
    "mito_fraction": n_mito / n_all if n_all else 0.0,
    "fragments_nuclear": nuclear,
    "fragments_final": n_kept,
    "NRF": distinct / total if total else 0.0,
    "PBC1": once / distinct if distinct else 0.0,
    "PBC2": once / twice if twice else float("inf"),
    "nucleosome_free_fraction": float(lengths[sizes < 147].sum() / kept),
    "mononucleosome_fraction": float(lengths[(sizes >= 147) & (sizes <= 294)].sum() / kept),
    "median_fragment_length": int(np.searchsorted(np.cumsum(lengths), kept / 2)),
    "tss_enrichment": tss_score,
    "tss_count": n_tss,
}
for name, n in sex.items():
    metrics[f"{name}_per_million"] = n / kept * 1e6
with open(sm.output.metrics, "w") as out:
    json.dump(metrics, out, indent=1)
with open(sm.output.fragments, "w") as out:
    out.write("length\tcount\n")
    for length, n in enumerate(lengths):
        out.write(f"{length}\t{n}\n")
with open(sm.output.tss, "w") as out:
    out.write("position\tcuts\tnormalised\tsmoothed\n")
    for i in range(len(profile)):
        out.write(f"{i - WINDOW}\t{profile[i]}\t{normalised[i]:.4f}\t{smoothed[i]:.4f}\n")
print(json.dumps(metrics, indent=1))
