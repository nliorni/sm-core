"""Toy demo script: counts lines in the two 'read' files for one sample.

Replace with a real per-sample processing step (typically a wrapped tool
call) when adapting this template to an actual workflow.
"""

from pathlib import Path

r1 = Path(snakemake.input.r1)
r2 = Path(snakemake.input.r2)
out = Path(snakemake.output[0])

n1 = sum(1 for _ in r1.open())
n2 = sum(1 for _ in r2.open())

out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(f"{snakemake.wildcards.sample}\t{n1}\t{n2}\n")
