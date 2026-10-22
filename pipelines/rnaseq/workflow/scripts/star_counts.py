#!/usr/bin/env python3
"""Merge STAR's per-sample ReadsPerGene.out.tab into one gene x sample matrix.

STAR writes three count columns per gene, one per possible library type:

    column 2  unstranded
    column 3  read 1 on the gene's strand          (forward-stranded library)
    column 4  read 1 on the opposite strand        (reverse-stranded, dUTP)

Picking the wrong one silently loses most of the counts, so the column is
chosen per sample from its effective strandedness (strandedness.tsv: the
declared value, or Salmon's inference with strandedness: auto).

The first four rows of each file are STAR's own bookkeeping (N_unmapped,
N_multimapping, N_noFeature, N_ambiguous); they go to the summary table, not
the matrix. Multimapping reads are never counted here -- that is what the
Salmon matrices are for.
"""

import argparse
import csv

COLUMN = {"unstranded": 1, "forward": 2, "reverse": 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--strandedness", required=True)
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--out-counts", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    with open(args.strandedness) as handle:
        strand = {row["sample"]: row["effective"] for row in csv.DictReader(handle, delimiter="\t")}
    names = {}
    with open(args.tx2gene) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            names.setdefault(row["gene_id"], row["gene_name"])

    genes, matrix, summary = [], {}, {}
    for sample, path in zip(args.samples, args.counts):
        col = COLUMN[strand[sample]]
        counts, stats = {}, {}
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if fields[0].startswith("N_"):
                    stats[fields[0]] = int(fields[col])
                else:
                    counts[fields[0]] = int(fields[col])
        if not genes:
            genes = list(counts)
        elif list(counts) != genes:
            raise SystemExit(f"{path}: gene list differs from the other samples (different STAR index?)")
        stats["N_assigned"] = sum(counts.values())
        stats["strandedness"] = strand[sample]
        matrix[sample], summary[sample] = counts, stats

    with open(args.out_counts, "w") as out:
        out.write("\t".join(["gene_id", "gene_name", *args.samples]) + "\n")
        for gene in genes:
            out.write("\t".join([gene, names.get(gene, gene), *(str(matrix[s][gene]) for s in args.samples)]) + "\n")

    keys = ["strandedness", "N_assigned", "N_noFeature", "N_ambiguous", "N_multimapping", "N_unmapped"]
    with open(args.out_summary, "w") as out:
        out.write("\t".join(["sample", *keys]) + "\n")
        for s in args.samples:
            out.write("\t".join([s, *(str(summary[s][k]) for k in keys)]) + "\n")


if __name__ == "__main__":
    main()
