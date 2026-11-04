###############################################################################
# sm-core shared: differential analysis designs (rnaseq `de`, atacseq `da`)
#
#     include: "../../smcore/rules/de_design.smk"
#
# Checks a config section's `analyses` against the samplesheet at parse time,
# so that a wrong column, level or design stops the run before any job starts,
# and lists the models and contrasts to run. The statistics and figures are in
# smcore/de.py.
#
# Needs from the pipeline: SAMPLE_TABLE, SAMPLES (smcore/common.smk's
# load_samples) and config_error.
#
# An analysis:
#   name       unique, becomes a folder name
#   design     a formula over samplesheet columns, e.g. "~ condition"
#   contrasts  [[factor, test level, reference level], ...]
#   split_by   optional column: one model per value (e.g. per cell line)
#   exclude    optional list of samples left out
###############################################################################

import re


def de_design_factors(design):
    # Column names in a formula: identifiers, minus formula keywords.
    names = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", design.split("~", 1)[1])
    return [n for n in dict.fromkeys(names) if n not in ("C", "Q")]


def _de_meta(sample, column):
    return SAMPLE_TABLE[sample]["meta"].get(column, "")


def _de_check_name(value, what):
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", str(value)):
        config_error(f"{what} '{value}' may only contain letters, digits and _ . + - (it becomes a folder name).")


def de_units(analyses, section):
    """[(analysis, split, (factor, test, reference)), ...] for the analyses
    of config section `section` (for messages), checked against the
    samplesheet."""
    columns = sorted({k for s in SAMPLES for k in SAMPLE_TABLE[s]["meta"]})
    names = [a["name"] for a in analyses]
    if len(set(names)) != len(names):
        config_error(f"{section}.analyses: every analysis needs a unique name.")
    units = []
    for a in analyses:
        name = a["name"]
        _de_check_name(name, f"{section}.analyses: name")
        where = f"{section}.analyses '{name}'"
        excluded = set(a.get("exclude") or [])
        unknown = excluded - set(SAMPLES)
        if unknown:
            config_error(f"{where}: exclude lists unknown samples {sorted(unknown)}.")
        samples = [s for s in SAMPLES if s not in excluded]
        factors = de_design_factors(a["design"])
        split_by = a.get("split_by")
        for col in factors + ([split_by] if split_by else []):
            if col not in columns:
                config_error(f"{where}: '{col}' is not a samplesheet column (have: {', '.join(columns) or 'none'}).")
        if split_by in factors:
            config_error(f"{where}: split_by '{split_by}' cannot also be in the design.")
        groups = {}
        for s in samples:
            groups.setdefault(_de_meta(s, split_by) if split_by else "all", []).append(s)
        seen = set()
        for factor, test, reference in a["contrasts"]:
            if factor not in factors:
                config_error(f"{where}: contrast factor '{factor}' is not in the design {a['design']}.")
            _de_check_name(test, f"{where}: level")
            _de_check_name(reference, f"{where}: level")
            if f"{test}_vs_{reference}" in seen:
                config_error(f"{where}: contrast {test}_vs_{reference} appears twice.")
            seen.add(f"{test}_vs_{reference}")
        for split, members in groups.items():
            _de_check_name(split, f"{where}: split value")
            for factor, test, reference in a["contrasts"]:
                levels = {_de_meta(s, factor) for s in members}
                for level in (test, reference):
                    if level not in levels:
                        config_error(
                            f"{where}: no sample with {factor} = '{level}'"
                            + (f" where {split_by} = '{split}'" if split_by else "")
                            + f" (levels: {', '.join(sorted(levels))})."
                        )
                units.append((name, split, (factor, test, reference)))
            # Residual degrees of freedom of an additive design: samples minus
            # coefficients. With none left, dispersions cannot be estimated.
            # (Designs with interactions are left to PyDESeq2 to judge.)
            if not re.search(r"[:*]", a["design"]):
                coefficients = 1 + sum(len({_de_meta(s, f) for s in members}) - 1 for f in factors)
                if len(members) - coefficients < 1:
                    config_error(
                        f"{where}: {len(members)} samples"
                        + (f" in '{split}'" if split_by else "")
                        + f" cannot fit {a['design']} ({coefficients} coefficients): no replicates left "
                        "to estimate dispersions."
                    )
    return units


def de_model_params(analyses, units, wildcards):
    """The `model` param of a rule: one analysis x split, and the contrast of
    the job if it has one."""
    a = next(x for x in analyses if x["name"] == wildcards.analysis)
    mine = [c for (n, s, c) in units if n == wildcards.analysis and s == wildcards.split]
    contrast = next(c for c in mine
                    if f"{c[1]}_vs_{c[2]}" == getattr(wildcards, "contrast", f"{c[1]}_vs_{c[2]}"))
    return {
        "design": a["design"],
        "split_by": a.get("split_by"),
        "exclude": a.get("exclude") or [],
        "contrast": list(contrast),
        "contrasts": [list(c) for c in mine],
    }
