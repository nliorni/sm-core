## ============================================================================
## Toy demo rule chain.
##
## Counts lines in each sample's two "read" files and writes a per-sample
## count, then aggregates all samples into one summary table. Deliberately
## dependency-free (stdlib only) so tests/tutorial/ runs anywhere Snakemake
## itself runs.
##
## Delete and replace once you write real rules -- it exists to show the
## conventions: inputs resolved from SAMPLE_TABLE (smcore/common.smk's
## load_samples(), never re-derived inline), retry-scaled resources via
## scaled(), software via env("<name>").
## ============================================================================


def reads(wildcards):
    # The first FASTQ pair of the sample; a real pipeline merges all units
    # (see pipelines/rnaseq's merge_lanes).
    return SAMPLE_TABLE[wildcards.sample]["units"][0]


rule count_lines:
    input:
        reads,
    output:
        outputDir + "counts/{sample}.linecount.tsv",
    threads: 1
    resources:
        mem_mb=scaled(500),
        runtime=scaled(5),
    conda:
        env("example")
    script:
        "../scripts/count_lines.py"


rule summarize:
    input:
        expand(outputDir + "counts/{sample}.linecount.tsv", sample=SAMPLES),
    output:
        outputDir + "summary.tsv",
    threads: 1
    resources:
        mem_mb=scaled(200),
        runtime=scaled(2),
    script:
        "../scripts/summarize.py"
