"""atacseq's side of the differential accessibility scripts: the peak count
matrix, library sizes and the two normalisations. The models, statistics and
figures are shared with rnaseq (smcore/de.py). The scripts put the sm-core
root on sys.path before importing this."""

from __future__ import annotations

import numpy as np
import pandas as pd

from smcore.de import library_size_factors

METHODS = ("median_of_ratios", "library_size")
LABEL = {"median_of_ratios": "median of ratios", "library_size": "library size"}


def read_peak_counts(counts_tsv: str):
    """Peaks x samples integer counts, and the peaks' coordinates."""
    table = pd.read_csv(counts_tsv, sep="\t", index_col="peak_id")
    coords = table[["chrom", "start", "end"]]
    return table.drop(columns=["chrom", "start", "end"]).astype(int), coords


def library_sizes(fc_summary: str, bams: dict) -> pd.Series:
    """Fragments in each final BAM: every featureCounts status row (assigned or
    not), i.e. all the read pairs it read. `bams` maps sample -> BAM path, the
    summary's column names."""
    summary = pd.read_csv(fc_summary, sep="\t", index_col=0)
    by_path = summary.sum(axis=0)
    return pd.Series({s: int(by_path[path]) for s, path in bams.items() if path in by_path.index})


def median_of_ratios(counts: pd.DataFrame) -> pd.Series:
    """DESeq2's size factors over the peaks with a fragment in every sample,
    exactly as DESeq2 and PyDESeq2 compute them (no rescaling)."""
    log = np.log(counts[(counts > 0).all(axis=1)].astype(float))
    if log.empty:
        raise ValueError("no peak has a fragment in every sample: median-of-ratios is undefined")
    return np.exp((log.sub(log.mean(axis=1), axis=0)).median(axis=0))


def size_factor_table(counts: pd.DataFrame, lib: pd.Series, used: str) -> pd.DataFrame:
    """Per sample: library size, fragments in consensus peaks, their ratio, and
    the size factors of both methods (the used one also as `size_factor`)."""
    lib = lib[counts.columns]
    table = pd.DataFrame({
        "library_size": lib,
        "fragments_in_peaks": counts.sum(axis=0),
    })
    table["fraction_in_peaks"] = table["fragments_in_peaks"] / table["library_size"]
    table["sf_median_of_ratios"] = median_of_ratios(counts)
    table["sf_library_size"] = library_size_factors(lib)
    table["size_factor"] = table[f"sf_{used}"]
    table.index.name = "sample"
    return table


def zero_shift(sf_used: pd.Series, sf_other: pd.Series, test_samples, reference_samples) -> float:
    """The log2 fold change, under the normalisation used, of a peak that the
    other normalisation calls unchanged: log2FC_used = log2FC_other + shift."""
    r = np.log2(sf_other / sf_used)
    return float(r[list(test_samples)].mean() - r[list(reference_samples)].mean())
