#!/usr/bin/env python3
"""Prepare the pipeline's annotation from a user GTF. Standard library only.

    gtf_tools.py prepare --gtf genes.gtf[.gz] --fai genome.fa.fai \
        [--spikein-fasta ercc.fa] --out-gtf genes.gtf --out-tx2gene tx2gene.tsv

1. Keeps only records on sequences present in the genome. A GTF carrying
   patches, alt contigs or another assembly's names would otherwise make
   gffread and STAR fail or warn about every one of them.
2. Drops transcripts left without any exon, which have no sequence.
3. Appends one gene/transcript/exon per spike-in sequence (gene_type
   "spikein"), spanning the whole sequence, so spike-ins are counted like
   genes. ThermoFisher's own ERCC92.gtf is not needed for this.
4. Writes tx2gene.tsv: transcript_id, gene_id, gene_name, gene_type, chrom --
   what tximport aggregates by and what annotates every count matrix.

Gene names and biotypes are read from GENCODE (gene_name, gene_type) or
Ensembl (gene_name, gene_biotype) attributes; a gene without a name is
named by its ID.
"""

import argparse
import gzip
import re
import sys
from collections import Counter

ATTR = re.compile(r'(\S+) "([^"]*)"')


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fai(path):
    with open(path) as handle:
        return {line.split("\t")[0]: int(line.split("\t")[1]) for line in handle}


def read_fasta_names(path):
    lengths, name = {}, None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def prepare(args):
    sequences = read_fai(args.fai)
    dropped = Counter()
    kept_lines = []
    tx_has_exon = set()
    transcripts = {}

    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                sys.exit(f"Malformed GTF line (expected 9 tab-separated fields): {line[:200]}")
            if fields[0] not in sequences:
                dropped[fields[0]] += 1
                continue
            attrs = dict(ATTR.findall(fields[8]))
            tx = attrs.get("transcript_id")
            if fields[2] == "exon":
                if not tx or "gene_id" not in attrs:
                    sys.exit(f"Exon without transcript_id/gene_id: {line[:200]}")
                tx_has_exon.add(tx)
            if tx and tx not in transcripts:
                gene_id = attrs.get("gene_id", "")
                transcripts[tx] = (
                    gene_id,
                    attrs.get("gene_name") or gene_id,
                    attrs.get("gene_type") or attrs.get("gene_biotype") or "unknown",
                    fields[0],
                )
            kept_lines.append((tx, line if line.endswith("\n") else line + "\n"))

    no_exon = set(transcripts) - tx_has_exon
    spikeins = read_fasta_names(args.spikein_fasta) if args.spikein_fasta else {}
    clash = [name for name in spikeins if name in transcripts]
    if clash:
        sys.exit(f"Spike-in names already used by annotated transcripts: {', '.join(clash[:5])}")

    with open(args.out_gtf, "w") as out:
        for tx, line in kept_lines:
            if tx in no_exon:
                continue
            out.write(line)
        for name, length in spikeins.items():
            if name not in sequences:
                sys.exit(f"Spike-in {name} is not in the genome FASTA; it must be appended to it.")
            attrs = (f'gene_id "{name}"; transcript_id "{name}"; gene_name "{name}"; '
                     f'gene_type "spikein"; transcript_type "spikein";')
            for feature in ("gene", "transcript", "exon"):
                a = attrs if feature != "gene" else f'gene_id "{name}"; gene_name "{name}"; gene_type "spikein";'
                out.write(f"{name}\tspikein\t{feature}\t1\t{length}\t.\t+\t.\t{a}\n")
            transcripts[name] = (name, name, "spikein", name)

    with open(args.out_tx2gene, "w") as out:
        out.write("transcript_id\tgene_id\tgene_name\tgene_type\tchrom\n")
        for tx, (gene_id, gene_name, gene_type, chrom) in transcripts.items():
            if tx in no_exon:
                continue
            out.write(f"{tx}\t{gene_id}\t{gene_name}\t{gene_type}\t{chrom}\n")

    kept = len(transcripts) - len(no_exon)
    print(f"transcripts kept: {kept} (incl. {len(spikeins)} spike-ins)", file=sys.stderr)
    print(f"transcripts dropped without exons: {len(no_exon)}", file=sys.stderr)
    for seq, n in dropped.most_common():
        print(f"records dropped, sequence not in genome: {seq}\t{n}", file=sys.stderr)
    if kept == 0:
        sys.exit("No transcripts left: do the GTF's sequence names match the genome's (chr1 vs 1)?")


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def _subtract(intervals, holes):
    """intervals minus holes; both sorted, merged, 1-based closed."""
    out, j = [], 0
    for start, end in intervals:
        cur = start
        while j < len(holes) and holes[j][1] < cur:
            j += 1
        k = j
        while k < len(holes) and holes[k][0] <= end:
            if holes[k][0] > cur:
                out.append([cur, holes[k][0] - 1])
            cur = max(cur, holes[k][1] + 1)
            k += 1
        if cur <= end:
            out.append([cur, end])
    return out


def collapse(args):
    """One transcript per gene, for RNA-SeQC -- the GTEx/ENCODE gene model.

    Each gene's exons from all its isoforms are merged into one exon set, so a
    read in any isoform's exon counts as exonic for the gene. Regions shared by
    exons of two or more genes on the same strand are then removed from all of
    them, so RNA-SeQC's gene-level metrics are not inflated by reads it cannot
    attribute. This is what GTEx's collapse_annotation.py does for its RNA-SeQC
    runs; the count matrices themselves never use this file.
    """
    genes = {}
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            g = genes.setdefault(a["gene_id"], {
                "chrom": f[0], "strand": f[6], "exons": [],
                "name": a.get("gene_name") or a["gene_id"],
                "type": a.get("gene_type") or a.get("gene_biotype") or "unknown",
            })
            g["exons"].append((int(f[3]), int(f[4])))

    for g in genes.values():
        g["exons"] = _merge(g["exons"])

    # Positions covered by exons of >= 2 genes, per chromosome and strand.
    by_strand = {}
    for gene_id, g in genes.items():
        by_strand.setdefault((g["chrom"], g["strand"]), []).extend(
            (s, e, gene_id) for s, e in g["exons"])
    shared = {}
    for key, ivs in by_strand.items():
        events = sorted([(s, 1) for s, e, _ in ivs] + [(e + 1, -1) for s, e, _ in ivs])
        depth, start, holes = 0, None, []
        for pos, delta in events:
            if depth >= 2 and start is not None and pos > start:
                holes.append([start, pos - 1])
            depth += delta
            start = pos
        shared[key] = _merge(holes)

    dropped = 0
    with open(args.out_gtf, "w") as out:
        for gene_id, g in genes.items():
            exons = _subtract(g["exons"], shared[(g["chrom"], g["strand"])])
            if not exons:
                dropped += 1
                continue
            base = (f'gene_id "{gene_id}"; transcript_id "{gene_id}"; '
                    f'gene_name "{g["name"]}"; gene_type "{g["type"]}";')
            lo, hi = exons[0][0], exons[-1][1]
            out.write(f'{g["chrom"]}\tcollapsed\tgene\t{lo}\t{hi}\t.\t{g["strand"]}\t.\t'
                      f'gene_id "{gene_id}"; gene_name "{g["name"]}"; gene_type "{g["type"]}";\n')
            out.write(f'{g["chrom"]}\tcollapsed\ttranscript\t{lo}\t{hi}\t.\t{g["strand"]}\t.\t{base}\n')
            for n, (s, e) in enumerate(exons, start=1):
                out.write(f'{g["chrom"]}\tcollapsed\texon\t{s}\t{e}\t.\t{g["strand"]}\t.\t'
                          f'{base} exon_id "{gene_id}_{n}"; exon_number "{n}";\n')
    print(f"genes: {len(genes)}, dropped (entirely shared with other genes): {dropped}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--gtf", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--spikein-fasta")
    p.add_argument("--out-gtf", required=True)
    p.add_argument("--out-tx2gene", required=True)
    c = sub.add_parser("collapse", help=collapse.__doc__.splitlines()[0])
    c.add_argument("--gtf", required=True)
    c.add_argument("--out-gtf", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    else:
        collapse(args)


if __name__ == "__main__":
    main()
