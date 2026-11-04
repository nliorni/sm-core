#!/usr/bin/env python3
"""Validate the ENCODE tutorial against ENCODE's own processing of the same reads.

    python3 tests/tutorial/encode/validate.py [--dir tests/tutorial/encode]

1. QC. Per library, sm-core's numbers next to the ENCODE ATAC-seq pipeline's
   (truth/encode_qc.tsv). The same reads and the same filters (MAPQ 30,
   proper pairs, Picard duplicates, chrM removed) should give the same
   mitochondrial fraction, duplicate fraction and library complexity (NRF,
   PBC1, PBC2): checked within a tolerance. FRiP and TSS enrichment are
   reported only: ENCODE computes FRiP on its own peak set and TSS enrichment
   with its own TSS list and scaling, so the values differ in level, not in
   ranking.

2. Peaks. For each cell line, the consensus peaks sm-core called in that
   line (consensus_peaks.support.tsv, groups_passing) against ENCODE's
   peaks for it:
     recall     the fraction of ENCODE's conservative IDR peaks (the
                reproducible core) overlapped by an sm-core peak; counted
                over distinct regions, as ENCODE's file repeats a region
                once per summit, and leaving out the regions in the ENCODE
                blacklist v2, which sm-core removes by design
     precision  the fraction of sm-core's peaks overlapped by one of ENCODE's
                pseudoreplicated peaks (ENCODE's full default set)

3. Biology. Normalised accessibility (CPM) of the promoter peaks of
   lineage genes: B-cell genes open in GM12878 (a B-lymphoblastoid line),
   erythroid genes open in K562 (an erythroleukaemia line).

4. Identity (somalier): libraries of one line are the same individual
   (relatedness near 1), the two lines are not (near 0).

5. Differential accessibility (-w da, K562 vs GM12878): the lineage genes of
   check 3 again, now as significant promoter peaks (genes.tsv): erythroid
   genes opened in K562, B-cell genes closed.

The tutorial runs replicate 1 of each line; with more libraries in the
samplesheet, every check covers them too.

Standard library only. Prints a report and exits non-zero if a check fails.
"""

import argparse
import csv
import gzip
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

TOLERANCE = {"mito_fraction": 0.03, "duplicate_fraction": 0.05, "NRF": 0.05, "PBC1": 0.05}
RECALL_MIN, PRECISION_MIN = 0.90, 0.80
MARKERS = {
    "GM12878": ["CD19", "CD79A", "CD79B", "PAX5", "MS4A1", "BLK"],
    "K562": ["GATA1", "TAL1", "KLF1", "HBG1", "HBG2", "GYPA"],
}
failures = []


def check(ok, message):
    print(("  PASS  " if ok else "  FAIL  ") + message)
    if not ok:
        failures.append(message)


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def intervals(rows):
    """{chrom: (sorted starts, matching ends)} with overlapping intervals merged."""
    by_chrom = defaultdict(list)
    for chrom, start, end in rows:
        by_chrom[chrom].append((start, end))
    merged = {}
    for chrom, ivs in by_chrom.items():
        out = []
        for s, e in sorted(ivs):
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = ([s for s, _ in out], [e for _, e in out])
    return merged


def overlaps(chrom, start, end, index):
    if chrom not in index:
        return False
    starts, ends = index[chrom]
    i = bisect_right(starts, end - 1) - 1
    return i >= 0 and ends[i] > start


def read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in handle if l.strip())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()
    base = Path(args.dir)
    results = base / "output" / "results"

    # ------------------------------------------------------------ 1. QC
    print("1. QC vs the ENCODE ATAC-seq pipeline (same reads)")
    ours = {r["sample"]: r for r in read_tsv(results / "qc" / "atac_qc_summary.tsv")}
    theirs = {r["sample"]: r for r in read_tsv(base / "truth" / "encode_qc.tsv")}
    cols = ["aligned_rate", "mito_fraction", "duplicate_fraction", "NRF", "PBC1", "PBC2", "FRiP", "tss_enrichment"]
    print(f"  {'sample':<11}" + "".join(f"{c[:14]:>16}" for c in cols))
    for s in ours:
        if s not in theirs:
            check(False, f"{s}: not in ENCODE's QC (truth/encode_qc.tsv)")
    run = [s for s in ours if s in theirs]
    for s in run:
        line = f"  {s:<11}"
        for c in cols:
            a, b = ours[s].get(c, ""), theirs[s].get(c, "")
            line += f"{(f'{float(a):.3f}' if a else '-'):>8}/{(f'{float(b):.3f}' if b else '-'):<7}"
        print(line)
    print("  (sm-core / ENCODE)")
    for s in run:
        for metric, tol in TOLERANCE.items():
            if theirs[s].get(metric):
                a, b = float(ours[s][metric]), float(theirs[s][metric])
                check(abs(a - b) <= tol, f"{s} {metric}: {a:.3f} vs ENCODE {b:.3f} (tolerance {tol})")
    # TSS enrichment and FRiP: same ranking of libraries, not the same scale.
    for metric in ("tss_enrichment", "FRiP"):
        both = [s for s in run if theirs[s].get(metric)]
        if len(both) > 2:
            rank = lambda d: {s: i for i, s in enumerate(sorted(both, key=lambda s: float(d[s][metric])))}
            ro, rt = rank(ours), rank(theirs)
            agree = sum(ro[s] == rt[s] for s in both)
            print(f"  info  {metric}: libraries ranked identically in {agree} of {len(both)} positions")

    # ------------------------------------------------------------ 2. peaks
    print("\n2. Consensus peaks vs ENCODE peaks")
    support = read_tsv(results / "peaks" / "consensus" / "consensus_peaks.support.tsv")
    blacklist = intervals(read_bed(base / "reference" / "hg38-blacklist.v2.bed.gz"))
    for cell in ("GM12878", "K562"):
        mine = [(r["chrom"], int(r["start"]), int(r["end"])) for r in support
                if cell in r["groups_passing"].split(",")]
        idr_all = sorted(set(read_bed(base / "truth" / f"{cell}_conservative_idr.bed.gz")))
        idr = [r for r in idr_all if not overlaps(*r, blacklist)]
        pseudo = intervals(read_bed(base / "truth" / f"{cell}_pseudoreplicated.bed.gz"))
        mine_index = intervals(mine)
        recall = sum(overlaps(c, s, e, mine_index) for c, s, e in idr) / len(idr)
        precision = sum(overlaps(c, s, e, pseudo) for c, s, e in mine) / max(len(mine), 1)
        print(f"  {cell}: {len(mine):,} sm-core peaks; ENCODE {len(idr):,} conservative IDR regions "
              f"(+ {len(idr_all) - len(idr):,} in the blacklist, left out)")
        check(recall >= RECALL_MIN, f"{cell} recall of ENCODE IDR peaks {recall:.3f} (>= {RECALL_MIN})")
        check(precision >= PRECISION_MIN, f"{cell} precision vs ENCODE pseudoreplicated peaks {precision:.3f} (>= {PRECISION_MIN})")

    # ------------------------------------------------------------ 3. biology
    print("\n3. Promoter accessibility of lineage genes (CPM, mean of the line's libraries)")
    counts = read_tsv(results / "counts" / "peak_counts.tsv")
    samples = [c for c in counts[0] if c not in ("peak_id", "chrom", "start", "end")]
    totals = {s: sum(int(r[s]) for r in counts) for s in samples}
    cpm = {r["peak_id"]: {s: int(r[s]) / totals[s] * 1e6 for s in samples} for r in counts}
    ann = read_tsv(results / "peaks" / "consensus" / "consensus_peaks.annotated.tsv")
    promoter = defaultdict(list)
    for r in ann:
        for name in (r.get("promoter_gene_names") or "").split(","):
            if name:
                promoter[name].append(r["peak_id"])
    of_line = {c: [s for s in samples if s.split("_")[0] == c] for c in MARKERS}
    for cell, genes in MARKERS.items():
        other = "K562" if cell == "GM12878" else "GM12878"
        higher = 0
        for g in genes:
            peaks = promoter.get(g, [])
            mean = {c: sum(cpm[p][s] for p in peaks for s in of_line[c]) / len(of_line[c]) for c in MARKERS}
            higher += mean[cell] > mean[other]
            print(f"  {g:<6} {cell:>8} {mean[cell]:8.1f}   {other:>8} {mean[other]:8.1f}   ({len(peaks)} peaks)")
        check(higher >= len(genes) - 1, f"{cell} lineage promoters more open in {cell}: {higher} of {len(genes)}")

    # ------------------------------------------------------------ 4. identity
    pairs_path = results / "qc" / "identity" / "somalier.pairs.tsv"
    if pairs_path.exists():
        print("\n4. Identity (somalier relatedness)")
        for r in read_tsv(pairs_path):
            a, b = r["#sample_a"], r["sample_b"]
            rel = float(r["relatedness"])
            same = a.split("_")[0] == b.split("_")[0]
            check(rel > 0.8 if same else rel < 0.2,
                  f"{a} vs {b}: relatedness {rel:.2f} ({'same line' if same else 'different lines'})")

    # ------------------------------------------------------------ 5. DA
    da_dir = results / "da" / "cell_lines" / "all" / "K562_vs_GM12878"
    if (da_dir / "genes.tsv").exists():
        print("\n5. Differential accessibility, K562 vs GM12878 (promoter peaks, genes.tsv)")
        by_name = {r["gene_name"]: r for r in read_tsv(da_dir / "genes.tsv")}
        summary = read_tsv(results / "da" / "summary.tsv")[0]
        print(f"  {int(summary['peaks_tested']):,} peaks tested: {int(summary['up']):,} more open in K562, "
              f"{int(summary['down']):,} in GM12878")
        for cell, expected in (("K562", "up"), ("GM12878", "down")):
            hits = 0
            for g in MARKERS[cell]:
                r = by_name.get(g, {})
                reg = r.get("promoter_regulation") or "no promoter peak"
                lfc = r.get("promoter_log2FoldChange") or ""
                hits += reg == expected
                print(f"  {g:<6} {reg:>16} {('log2FC ' + format(float(lfc), '+.2f')) if lfc else ''}")
            check(hits >= len(MARKERS[cell]) - 1,
                  f"{cell} lineage promoters significantly {expected} in K562 vs GM12878: {hits} of {len(MARKERS[cell])}")

    print("\n" + ("All checks passed." if not failures else f"{len(failures)} check(s) failed."))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
