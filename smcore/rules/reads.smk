###############################################################################
# sm-core shared rules: reads -- lane merging and fastp
#
# Included by every short-read pipeline after its own workflow/rules/common.smk:
#
#     include: "../../smcore/rules/reads.smk"
#
# Needs from the pipeline: RESULTS, LOGS, BENCH, SAMPLE_TABLE (load_samples())
# and a `trimming:` config section (enabled, min_length, extra, threads,
# mem_mb) in its schema. Provides raw_reads() and analysis_reads(), the reads
# a pipeline's first real step takes.
###############################################################################


def trimming_enabled():
    return config["trimming"]["enabled"]


## ============================================================================
## Reads
## ============================================================================

def raw_reads(wildcards):
    """[R1, R2] for a sample: its single FASTQ pair as given, or the lanes /
    runs concatenated into one pair by merge_lanes when there are several."""
    units = SAMPLE_TABLE[wildcards.sample]["units"]
    if len(units) == 1:
        return list(units[0])
    return [
        RESULTS + f"fastq/{wildcards.sample}_R1.merged.fastq.gz",
        RESULTS + f"fastq/{wildcards.sample}_R2.merged.fastq.gz",
    ]


def analysis_reads(wildcards):
    """What the pipeline's aligner reads: fastp's trimmed reads, or the raw
    reads with trimming.enabled: false."""
    if trimming_enabled():
        return [
            RESULTS + f"fastp/{wildcards.sample}_R1.trimmed.fastq.gz",
            RESULTS + f"fastp/{wildcards.sample}_R2.trimmed.fastq.gz",
        ]
    return raw_reads(wildcards)


## ============================================================================
## Lane merging and fastp
## ============================================================================

def _units(sample, mate):
    return [pair[mate] for pair in SAMPLE_TABLE[sample]["units"]]


rule merge_lanes:
    # Only reached for samples declared with several FASTQ pairs (lanes,
    # flowcells, top-up runs). Concatenated gzip members are themselves a
    # valid gzip stream, so this is a plain `cat`: no recompression. R1 and R2
    # are concatenated in the same unit order, which keeps the mates paired.
    input:
        r1=lambda wildcards: _units(wildcards.sample, 0),
        r2=lambda wildcards: _units(wildcards.sample, 1),
    output:
        r1=temp(RESULTS + "fastq/{sample}_R1.merged.fastq.gz"),
        r2=temp(RESULTS + "fastq/{sample}_R2.merged.fastq.gz"),
    resources:
        mem_mb=scaled(500),
        runtime=scaled(60),
    run:
        for mate in ("r1", "r2"):
            sources = input[mate]
            if not all(f.endswith(".gz") for f in sources):
                raise WorkflowError(
                    f"{wildcards.sample}: merging lanes needs gzipped FASTQ; got {', '.join(sources)}"
                )
            shell("cat {sources} > {out}", sources=" ".join(sources), out=output[mate])


def _fastp_outputs():
    reports = {
        "json": RESULTS + "fastp/{sample}.fastp.json",
        "html": RESULTS + "fastp/{sample}.fastp.html",
    }
    if trimming_enabled():
        # temp(): the trimmed reads only exist to be aligned/quantified.
        reports["r1"] = temp(RESULTS + "fastp/{sample}_R1.trimmed.fastq.gz")
        reports["r2"] = temp(RESULTS + "fastp/{sample}_R2.trimmed.fastq.gz")
    return reports


rule fastp:
    # Always runs, for its before/after read QC (MultiQC reads the JSON). With
    # trimming.enabled it also writes the trimmed reads the aligner uses;
    # without it, no reads are written and the aligner reads the raw FASTQ.
    #
    # Adapters are found by mate overlap (--detect_adapter_for_pe), which
    # needs no adapter sequence and works for TruSeq, Nextera and MGI alike.
    # Poly-G tails from two-colour instruments (NovaSeq/NextSeq) are trimmed
    # automatically when fastp recognises those read headers.
    input:
        reads=raw_reads,
    output:
        **_fastp_outputs(),
    log:
        LOGS + "fastp/{sample}.log",
    benchmark:
        BENCH + "fastp/{sample}.tsv"
    conda:
        env("fastp")
    threads: config["trimming"]["threads"]
    resources:
        mem_mb=scaled(config["trimming"]["mem_mb"]),
        runtime=scaled(120),
    params:
        out=lambda wildcards, output: (
            f"-o {output.r1} -O {output.r2} --detect_adapter_for_pe "
            f"--length_required {config['trimming']['min_length']} {config['trimming']['extra']}"
            if trimming_enabled()
            else "--disable_adapter_trimming --disable_quality_filtering --disable_length_filtering"
        ),
    shell:
        """
        fastp \
            -i {input.reads[0]} -I {input.reads[1]} \
            {params.out} \
            --thread {threads} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample} \
            2> {log}
        """
