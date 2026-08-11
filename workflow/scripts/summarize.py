"""Toy demo script: aggregates all samples' per-sample counts into one table."""

from pathlib import Path

out = Path(snakemake.output[0])
rows = [Path(f).read_text().strip() for f in snakemake.input]

out.write_text("sample\tr1_lines\tr2_lines\n" + "\n".join(rows) + "\n")
