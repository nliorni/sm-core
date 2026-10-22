#!/usr/bin/env python3
"""Share of each sample's counts by gene class, as MultiQC custom content.

Classes, first match wins:
    mitochondrial  any gene on chrM / MT / M
    rRNA           gene_type rRNA, rRNA_pseudogene, Mt_rRNA (nuclear rRNA)
    spikein        gene_type spikein (reference.spikein_fasta)
    protein_coding
    lncRNA         lncRNA, or the older lincRNA / antisense / processed_transcript
    other          everything else (pseudogenes, small RNAs, ...)

A high mitochondrial share points at stressed or dying cells, or at a
cytoplasm-enriched prep; a high rRNA share at failed rRNA depletion or poly-A
selection. Both eat sequencing depth that the rest of the matrix then lacks.
"""

import argparse
import csv
import sys

MITO = {"chrM", "MT", "M", "chrMT"}
LNC = {"lncRNA", "lincRNA", "antisense", "processed_transcript", "sense_intronic", "sense_overlapping"}
CLASSES = ["protein_coding", "lncRNA", "rRNA", "mitochondrial", "spikein", "other"]


def classify(gene_type, chrom):
    if chrom in MITO:
        return "mitochondrial"
    if "rRNA" in gene_type:
        return "rRNA"
    if gene_type == "spikein":
        return "spikein"
    if gene_type == "protein_coding":
        return "protein_coding"
    if gene_type in LNC:
        return "lncRNA"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--branch", required=True)
    args = parser.parse_args()

    gene_class = {}
    with open(args.tx2gene) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_class.setdefault(row["gene_id"], classify(row["gene_type"], row["chrom"]))

    with open(args.counts) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[2:]
        totals = {s: dict.fromkeys(CLASSES, 0.0) for s in samples}
        for row in reader:
            cls = gene_class.get(row[0], "other")
            for s, value in zip(samples, row[2:]):
                totals[s][cls] += float(value)

    out = sys.stdout
    out.write(f"# id: 'biotypes_{args.branch}'\n")
    out.write(f"# section_name: 'Count composition ({args.branch})'\n")
    # Double quotes: MultiQC parses these header lines as YAML, and an
    # apostrophe inside a single-quoted value ends it early.
    out.write('# description: "Percentage of each sample\'s counts by gene class. '
              'High rRNA or mitochondrial shares mean less depth for everything else."\n')
    out.write("# plot_type: 'bargraph'\n")
    out.write("# pconfig:\n#    id: 'biotypes_%s_plot'\n#    cpswitch_c_active: False\n" % args.branch)
    out.write("\t".join(["Sample", *CLASSES]) + "\n")
    for s in samples:
        total = sum(totals[s].values()) or 1.0
        out.write("\t".join([s, *(f"{100 * totals[s][c] / total:.3f}" for c in CLASSES)]) + "\n")


if __name__ == "__main__":
    main()
