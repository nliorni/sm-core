#!/usr/bin/env python3

"""
Launcher for sm-core/rnaseq: paired-end short-read RNA-seq, FASTQ to count
matrices.

Examples
--------
Build the reference indices once:
    ./run.py -w reference -c config/config.yaml -q 16

Dry-run, then run, both branches:
    ./run.py -w all -c config_MyProject.yaml -q 32 -n
    ./run.py -w all -c config_MyProject.yaml -q 32

Submit one PBS job per rule instead of running locally:
    ./run.py -w all -c config_MyProject.yaml -cl -qu workq -j 20

All flags: ./run.py --help (shared by every sm-core pipeline; see
smcore/launcher.py).
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # the sm-core repository root

from smcore.launcher import Pipeline, main  # noqa: E402

PIPELINE = Pipeline(
    name="rnaseq",
    version="0.1.0",
    description="paired-end RNA-seq: FASTQ to count matrices",
    directory=HERE,
    workflows={
        "all": "Both branches below, plus a comparison of their strandedness calls.",
        "star": "STAR branch: genome BAM, STAR gene count matrix, RNA-SeQC, MultiQC.",
        "salmon": "Salmon branch: Salmon quant, tximport gene/transcript matrices, MultiQC.",
        "reference": "Build the reference once: GTF filtering, transcriptome, STAR/Salmon indices.",
        "trim": "fastp adapter/quality trimming only.",
    },
)

if __name__ == "__main__":
    raise SystemExit(main(PIPELINE))
