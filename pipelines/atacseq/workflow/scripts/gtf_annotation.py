#!/usr/bin/env python3
"""Genes and TSSs from a GTF, for ATAC-seq. Standard library only.

    gtf_annotation.py --gtf genes.gtf[.gz] --fai genome.fa.fai \
        [--out-genes genes.tsv] [--out-tss tss.bed]

genes.tsv   gene_id, gene_name, gene_type, chrom, start, end, strand, tss
            (0-based start, half-open; tss = 0-based TSS of the gene's
            Ensembl_canonical transcript when the GTF tags one, else of the
            gene's 5' end) -- what peaks are annotated against. The gene's 5'
            end is often a minor transcript far upstream of the main
            promoter (BLK: 7.5 kb), hence the canonical one.
tss.bed     one line per distinct TSS of a protein-coding transcript, with
            strand -- what the TSS enrichment score aggregates over. GENCODE's
            Ensembl_canonical transcripts are used when the GTF tags them
            (one well-supported TSS per gene, as ENCODE's TSS set), all
            protein-coding transcripts otherwise.

Only records on sequences of the genome are kept. Gene names and types are
read from GENCODE (gene_name, gene_type) or Ensembl (gene_name,
gene_biotype) attributes.
"""

import argparse
import gzip
import re
import sys

ATTR = re.compile(r'(\S+) "([^"]*)"')


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-genes")
    ap.add_argument("--out-tss")
    args = ap.parse_args()
    if not (args.out_genes or args.out_tss):
        ap.error("give --out-genes, --out-tss or both")

    with open(args.fai) as handle:
        seqs = {line.split("\t")[0] for line in handle}

    genes = []
    canonical, coding = set(), set()
    gene_tss = {}  # gene_id -> TSS of its Ensembl_canonical transcript
    for line in open_text(args.gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in seqs or f[2] not in ("gene", "transcript"):
            continue
        attrs = dict(ATTR.findall(f[8]))
        start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
        tss = start if strand == "+" else end - 1
        if f[2] == "gene":
            gene_id = attrs.get("gene_id", "")
            genes.append((
                gene_id,
                attrs.get("gene_name", gene_id),
                attrs.get("gene_type", attrs.get("gene_biotype", "")),
                f[0], start, end, strand, tss,
            ))
            continue
        is_canonical = "Ensembl_canonical" in re.findall(r'tag "([^"]*)"', f[8])
        if is_canonical:
            gene_tss[attrs.get("gene_id", "")] = tss
        ttype = attrs.get("transcript_type", attrs.get("transcript_biotype", ""))
        if ttype != "protein_coding":
            continue
        site = (f[0], tss, strand)
        coding.add(site)
        if is_canonical:
            canonical.add(site)

    if args.out_genes:
        with open(args.out_genes, "w") as out:
            out.write("gene_id\tgene_name\tgene_type\tchrom\tstart\tend\tstrand\ttss\n")
            for g in sorted(genes, key=lambda g: (g[3], g[4])):
                g = g[:7] + (gene_tss.get(g[0], g[7]),)
                out.write("\t".join(map(str, g)) + "\n")
        print(f"{len(genes)} genes ({sum(g[0] in gene_tss for g in genes)} with a canonical transcript's TSS)",
              file=sys.stderr)

    if args.out_tss:
        sites = canonical or coding
        if not sites:
            sys.exit("no protein-coding transcripts found: the TSS enrichment score needs them.")
        with open(args.out_tss, "w") as out:
            for chrom, pos, strand in sorted(sites):
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t.\t0\t{strand}\n")
        print(f"{len(sites)} TSSs from {'Ensembl_canonical' if canonical else 'all protein-coding'} transcripts",
              file=sys.stderr)


if __name__ == "__main__":
    main()
