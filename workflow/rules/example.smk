## ============================================================================
## Toy demo rule chain.
##
## Counts lines in each sample's two "read" files and writes a per-sample
## count, then aggregates all samples into one summary table. Deliberately
## dependency-free (stdlib only) so `tests/tutorial/` runs anywhere Snakemake
## itself runs, with no conda/Singularity/cluster tooling required.
##
## This whole file is meant to be deleted and replaced once you adapt this
## template to a real workflow -- it exists to demonstrate the conventions
## (paired-input resolution via input_1/input_2, retry-scaled resources via
## scaled(), a conda: env reference) that the rest of this template's docs
## and run.py flags assume.
## ============================================================================

EXAMPLE_ENV = "../envs/example.yaml"


rule count_lines:
    input:
        r1=input_1,
        r2=input_2,
    output:
        outputDir + "counts/{sample}.linecount.tsv",
    threads: 1
    resources:
        mem_mb=scaled(500),
        runtime=scaled(5),
    conda:
        EXAMPLE_ENV
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
