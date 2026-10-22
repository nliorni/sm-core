#!/usr/bin/env python3
"""Validate the SEQC tutorial's count matrices against two published truth sets.

    python3 tests/tutorial/seqc/validate.py [--output-dir tests/tutorial/seqc/output]

Both branches (counts/star/gene_counts.tsv, counts/salmon/gene_counts.tsv) are
scored the same way: library sizes normalised by DESeq2's median-of-ratios
method, then the mean of the three A and the three B libraries per gene, and
log2(B/A).

1. ERCC spike-ins (92 synthetic transcripts). Sample A carries Mix 1, sample B
   Mix 2, and the two mixes differ by known ratios in four subgroups of 23:
   Mix1/Mix2 = 4 (group A), 1 (B), 2/3 (C), 1/2 (D). UHRR and HBRR do not hold
   the same mRNA mass per microgram of total RNA, which shifts every spike-in
   ratio by one constant, so ratios are centred on subgroup B (expected 1)
   before comparison -- as in the ERCC dashboard (Munro et al. 2014).
   Reported: the median observed log2(A/B) per subgroup vs. its expected value,
   and the correlation of observed and expected over all detected spike-ins.

2. MAQC TaqMan qPCR (GEO GPL4097): ~1,000 genes measured in A and B, four
   replicates each, normalised to POLR2A. For genes detected in all eight
   (Ct <= 35) whose symbol names exactly one gene in the annotation, the qPCR
   log2(B/A) is compared with each matrix's. POLR2A normalisation shifts every
   qPCR ratio by one constant, so the correlation is the measure, not the
   absolute difference.

Finally the two branches are compared with each other over every gene with a
mean of at least 10 normalised counts in both.

Standard library only. Prints a report; writes nothing.
"""

import argparse
import csv
import gzip
import math
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_matrix(path):
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[2:]
        genes, names, counts = [], {}, {s: [] for s in samples}
        for row in reader:
            genes.append(row[0])
            names[row[0]] = row[1]
            for s, v in zip(samples, row[2:]):
                counts[s].append(float(v))
    return genes, names, samples, counts


def size_factors(samples, counts):
    """DESeq2's median of ratios to the per-gene geometric mean."""
    n = len(counts[samples[0]])
    logmeans = []
    for i in range(n):
        values = [counts[s][i] for s in samples]
        logmeans.append(sum(math.log(v) for v in values) / len(values) if min(values) > 0 else None)
    factors = {}
    for s in samples:
        ratios = [math.log(counts[s][i]) - lm for i, lm in enumerate(logmeans) if lm is not None]
        factors[s] = math.exp(statistics.median(ratios))
    return factors


def condition_means(path):
    genes, names, samples, counts = read_matrix(path)
    sf = size_factors(samples, counts)
    a = [s for s in samples if "_A_" in s]
    b = [s for s in samples if "_B_" in s]
    mean_a, mean_b = {}, {}
    for i, g in enumerate(genes):
        mean_a[g] = statistics.mean(counts[s][i] / sf[s] for s in a)
        mean_b[g] = statistics.mean(counts[s][i] / sf[s] for s in b)
    return names, mean_a, mean_b, sf


def pearson(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def read_ercc(path):
    expected, group = {}, {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected[row["ERCC ID"]] = float(row["log2(Mix 1/Mix 2)"])
            group[row["ERCC ID"]] = row["subgroup"]
    return expected, group


def read_taqman(path):
    """{gene symbol: qPCR log2(B/A)} for assays detected in all A and B replicates."""
    symbol, values, current, in_table = {}, defaultdict(dict), None, None
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table, header = "platform", None
                continue
            if line.startswith("!platform_table_end") or line.startswith("!sample_table_end"):
                in_table = None
                continue
            if line.startswith("!Sample_title"):
                current = line.split("=", 1)[1].strip()          # MAQC_TAQ_1_A1
                continue
            if line.startswith("!sample_table_begin"):
                in_table, header = "sample", None
                continue
            if in_table is None:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields))
            if in_table == "platform":
                symbol[row["ID"]] = row["ORF"]
            elif current and current.split("_")[-1][0] in "AB":
                values[row["ID_REF"]][current.split("_")[-1]] = (float(row["VALUE"]), row["Flag_Detection"])
    out = {}
    for assay, reps in values.items():
        a = [v for k, (v, f) in reps.items() if k.startswith("A") and f == "P"]
        b = [v for k, (v, f) in reps.items() if k.startswith("B") and f == "P"]
        if len(a) == 4 and len(b) == 4 and min(a + b) > 0:
            gene = symbol.get(assay)
            if gene:
                out[gene] = statistics.mean(math.log2(v) for v in b) - statistics.mean(math.log2(v) for v in a)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", default=str(HERE / "output"))
    parser.add_argument("--truth-dir", default=str(HERE / "truth"))
    args = parser.parse_args()

    counts_dir = Path(args.output_dir) / "results" / "counts"
    ercc_expected, ercc_group = read_ercc(Path(args.truth_dir) / "ERCC_Controls_Analysis.txt")
    taqman = read_taqman(Path(args.truth_dir) / "GPL4097_family.soft.gz")

    fold = {}
    for branch in ("star", "salmon"):
        path = counts_dir / branch / "gene_counts.tsv"
        if not path.exists():
            print(f"[{branch}] {path} not found -- skipped")
            continue
        names, mean_a, mean_b, sf = condition_means(path)
        print(f"\n=== {branch} ===")
        print("size factors: " + ", ".join(f"{s} {f:.2f}" for s, f in sf.items()))

        # 1. ERCC
        obs = {e: math.log2((mean_a[e] + 0.5) / (mean_b[e] + 0.5))
               for e in ercc_expected if e in mean_a and min(mean_a[e], mean_b[e]) >= 10}
        if obs:
            shift = statistics.median(v for e, v in obs.items() if ercc_group[e] == "B")
            print(f"\nERCC spike-ins detected (>= 10 normalised counts in A and B): {len(obs)} of 92")
            print(f"  mRNA-content shift (median log2 A/B of ratio-1 group): {shift:+.2f}")
            print("  subgroup  expected log2(A/B)  observed (median, centred)  n")
            for g, exp in (("A", 2.0), ("B", 0.0), ("C", math.log2(2 / 3)), ("D", -1.0)):
                vals = [v - shift for e, v in obs.items() if ercc_group[e] == g]
                if vals:
                    print(f"  {g:8}  {exp:+18.2f}  {statistics.median(vals):+26.2f}  {len(vals):2d}")
            e_list = sorted(obs)
            print(f"  Pearson r (observed vs expected, n={len(e_list)}): "
                  f"{pearson([obs[e] for e in e_list], [ercc_expected[e] for e in e_list]):.3f}")

        # 2. TaqMan
        by_name = defaultdict(list)
        for g, n in names.items():
            by_name[n].append(g)
        pairs = [(taqman[sym], math.log2((mean_b[gs[0]] + 0.5) / (mean_a[gs[0]] + 0.5)))
                 for sym, gs in by_name.items() if sym in taqman and len(gs) == 1]
        if pairs:
            q, r = zip(*pairs)
            print(f"\nTaqMan qPCR genes matched (detected in all 8 assays, unique symbol): {len(pairs)}")
            print(f"  log2(B/A) RNA-seq vs qPCR: Pearson r = {pearson(q, r):.3f}, Spearman rho = {spearman(q, r):.3f}")

        fold[branch] = {g: math.log2((mean_b[g] + 0.5) / (mean_a[g] + 0.5))
                        for g in mean_a if min(mean_a[g], mean_b[g]) >= 10}

    if len(fold) == 2:
        common = sorted(set(fold["star"]) & set(fold["salmon"]))
        x = [fold["star"][g] for g in common]
        y = [fold["salmon"][g] for g in common]
        print(f"\n=== STAR vs Salmon ===\nlog2(B/A) over {len(common)} genes expressed in both: "
              f"Pearson r = {pearson(x, y):.3f}")
        big = sum(1 for a, b in zip(x, y) if abs(a - b) > 1)
        print(f"genes whose fold changes differ by more than 2x between branches: {big}")


if __name__ == "__main__":
    main()
