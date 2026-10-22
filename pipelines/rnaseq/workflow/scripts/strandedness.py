#!/usr/bin/env python3
"""Infer a sample's library strandedness from one branch's own evidence.

    strandedness.py salmon --lib-format-counts lib_format_counts.json ...
    strandedness.py star   --reads-per-gene    ReadsPerGene.out.tab   ...

salmon: Salmon counts, for every mapped fragment, the orientation it was
  observed in -- ISF (read 1 on the transcript's strand) or ISR (read 1 on the
  opposite strand), SF/SR for orphaned mates -- whatever --libType it was
  given.
star: STAR's ReadsPerGene.out.tab counts every gene three ways: column 3
  assuming read 1 is on the gene's strand, column 4 assuming the opposite.
  In a stranded library one of the two collects almost everything; in an
  unstranded one they split evenly.

Either way the forward share says how the library was made:

    forward share >= stranded threshold   -> forward  (e.g. Ligation, SMARTer Stranded)
    reverse share >= stranded threshold   -> reverse  (dUTP: TruSeq Stranded, NEB Directional)
    |forward - reverse| <= unstranded thr -> unstranded
    anything in between                   -> undetermined (degraded, gDNA-rich, mixed)

The thresholds are nf-core/rnaseq's (0.8 / 0.1). Writes a header line and one
row to stdout.
"""

import argparse
import json


def salmon_evidence(path):
    with open(path) as handle:
        lib = json.load(handle)
    forward = sum(lib.get(k, 0) for k in ("ISF", "SF"))
    reverse = sum(lib.get(k, 0) for k in ("ISR", "SR"))
    return forward, reverse


def star_evidence(path):
    forward = reverse = 0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0].startswith("N_"):
                continue
            forward += int(fields[2])
            reverse += int(fields[3])
    return forward, reverse


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("branch", choices=["salmon", "star"])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--lib-format-counts")
    parser.add_argument("--reads-per-gene")
    parser.add_argument("--declared", required=True, choices=["auto", "unstranded", "forward", "reverse"])
    parser.add_argument("--stranded-threshold", type=float, default=0.8)
    parser.add_argument("--unstranded-threshold", type=float, default=0.1)
    args = parser.parse_args()

    if args.branch == "salmon":
        forward, reverse = salmon_evidence(args.lib_format_counts)
    else:
        forward, reverse = star_evidence(args.reads_per_gene)

    total = forward + reverse
    fwd = forward / total if total else 0.0
    rev = reverse / total if total else 0.0

    if total == 0:
        inferred = "undetermined"
    elif fwd >= args.stranded_threshold:
        inferred = "forward"
    elif rev >= args.stranded_threshold:
        inferred = "reverse"
    elif abs(fwd - rev) <= args.unstranded_threshold:
        inferred = "unstranded"
    else:
        inferred = "undetermined"

    if args.declared == "auto":
        status = "inferred" if inferred != "undetermined" else "undetermined"
    elif inferred == "undetermined":
        status = "undetermined"
    else:
        status = "ok" if inferred == args.declared else "MISMATCH"

    # What this branch's downstream steps assume: the declared value when
    # there is one -- for Salmon, that is also what it quantified with -- and
    # otherwise what the data says, falling back to unstranded, which never
    # discards reads, when the data is ambiguous.
    if args.declared != "auto":
        effective = args.declared
    elif inferred != "undetermined":
        effective = inferred
    else:
        effective = "unstranded"

    header = ["sample", "branch", "declared", "inferred", "effective", "status",
              "forward_fraction", "reverse_fraction", "evidence"]
    row = [args.sample, args.branch, args.declared, inferred, effective, status,
           f"{fwd:.4f}", f"{rev:.4f}", str(total)]
    print("\t".join(header))
    print("\t".join(row))


if __name__ == "__main__":
    main()
