#!/usr/bin/env python3

"""
Launcher for sm-core/atacseq: paired-end bulk ATAC-seq, FASTQ to filtered
BAMs, peaks, a consensus peak count matrix, signal tracks and ENCODE QC.

Examples
--------
Build the reference once (bowtie2 index, TSS annotation, blacklist):
    ./run.py -w reference -c config/config.yaml -q 16

Dry-run, then run:
    ./run.py -w all -c config_MyProject.yaml -q 32 -n
    ./run.py -w all -c config_MyProject.yaml -q 32

Differential accessibility only (contrasts in the config's da section):
    ./run.py -w da -c config_MyProject.yaml -q 16

All flags: ./run.py --help (shared by every sm-core pipeline; see
smcore/launcher.py).
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # the sm-core repository root

from smcore.launcher import Pipeline, main  # noqa: E402

PIPELINE = Pipeline(
    name="atacseq",
    version="0.1.0",
    description="paired-end ATAC-seq: FASTQ to peaks, a consensus count matrix, signal tracks and QC",
    directory=HERE,
    workflows={
        "all": "Everything: filtered BAMs, per-sample and consensus peaks, peak counts, bigWigs, QC, MultiQC "
        "(and differential accessibility, when da.analyses is set).",
        "da": "Differential accessibility (PyDESeq2) of the consensus peaks per configured contrast.",
        "reference": "Build the reference once: genome, bowtie2 index, TSS and gene annotation, blacklist.",
        "trim": "fastp adapter/quality trimming only.",
    },
)

if __name__ == "__main__":
    raise SystemExit(main(PIPELINE))
