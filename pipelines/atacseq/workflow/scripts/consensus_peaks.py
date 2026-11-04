"""Consensus peak set: every sample's peaks merged, kept where enough
replicates of one group agree.

Overlapping (or touching) peaks of all samples are merged into regions. A
region is kept when at least `min_replicates` samples of one group called a
peak in it (capped at the group's size, so a group of one still counts), so a
region open in only one condition survives, while a peak seen in a single
replicate does not. This is nf-core/atacseq's consensus, with the replicate
requirement per group instead of over all samples.

Outputs: BED (peak_id in genome order), SAF for featureCounts, and a support
table (how many samples of each group called each region)."""

import sys
from collections import defaultdict

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")

samples = sm.params.samples
groups = sm.params.groups
min_rep = sm.params.min_replicates
group_of = {s: g for g, members in groups.items() for s in members}
need = {g: min(min_rep, len(members)) for g, members in groups.items()}

with open(sm.input.fai) as handle:
    order = [line.split("\t")[0] for line in handle]
rank = {c: i for i, c in enumerate(order)}

peaks = defaultdict(list)
for sample, path in zip(samples, sm.input.peaks):
    n = 0
    with open(path) as handle:
        for line in handle:
            f = line.split("\t")
            peaks[f[0]].append((int(f[1]), int(f[2]), sample))
            n += 1
    print(f"{sample}: {n} peaks")

regions = []
for chrom in sorted(peaks, key=lambda c: rank.get(c, len(rank))):
    current = None
    for start, end, sample in sorted(peaks[chrom]):
        if current and start <= current[2]:
            current[2] = max(current[2], end)
            current[3].add(sample)
        else:
            if current:
                regions.append(current)
            current = [chrom, start, end, {sample}]
    if current:
        regions.append(current)

group_names = list(groups)
kept = []
for chrom, start, end, members in regions:
    per_group = {g: 0 for g in group_names}
    for s in members:
        per_group[group_of[s]] += 1
    passing = [g for g in group_names if per_group[g] >= need[g]]
    if passing:
        kept.append((chrom, start, end, len(members), per_group, passing))

width = max(6, len(str(len(kept))))
with open(sm.output.bed, "w") as bed, open(sm.output.saf, "w") as saf, open(sm.output.support, "w") as sup:
    saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    sup.write("\t".join(["peak_id", "chrom", "start", "end", "n_samples", "groups_passing"] + group_names) + "\n")
    for i, (chrom, start, end, n, per_group, passing) in enumerate(kept, start=1):
        pid = f"peak_{i:0{width}d}"
        bed.write(f"{chrom}\t{start}\t{end}\t{pid}\t{n}\t.\n")
        saf.write(f"{pid}\t{chrom}\t{start + 1}\t{end}\t+\n")
        sup.write("\t".join([pid, chrom, str(start), str(end), str(n), ",".join(passing)]
                            + [str(per_group[g]) for g in group_names]) + "\n")

print(f"{len(regions)} merged regions; {len(kept)} kept "
      f"(>= {min_rep} replicates of one of {len(groups)} group(s): "
      + ", ".join(f"{g} n={len(m)}" for g, m in groups.items()) + ")")
if not kept:
    sys.exit("no consensus peak: check the per-sample peaks and peaks.min_replicates.")
