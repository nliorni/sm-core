###############################################################################
# Salmon branch: transcript quantification straight from the reads
#
# Selective alignment against the decoy-aware index (reference.smk), with
# --gcBias and --seqBias (salmon.extra) correcting fragment-GC and
# random-hexamer priming biases, which otherwise distort isoform-level
# estimates and, through them, gene-level counts (Patro et al. 2017; Love et
# al. 2016). The library type comes from each sample's strandedness: A
# (auto-detect) unless the samplesheet or config states it.
#
# Independent of the STAR branch: nothing here reads its output.
###############################################################################

SALMON_DIR = RESULTS + "salmon/{sample}"


rule salmon_quant:
    input:
        reads=analysis_reads,
        index=salmon_index(),
    output:
        quant=SALMON_DIR + "/quant.sf",
        lib=SALMON_DIR + "/lib_format_counts.json",
        meta=SALMON_DIR + "/aux_info/meta_info.json",
    log:
        LOGS + "salmon/{sample}.log",
    benchmark:
        BENCH + "salmon/{sample}.tsv"
    conda:
        env("salmon")
    threads: config["salmon"]["threads"]
    resources:
        mem_mb=scaled(config["salmon"]["mem_mb"]),
        runtime=scaled(120),
    params:
        outdir=SALMON_DIR,
        libtype=salmon_libtype,
        extra=config["salmon"]["extra"],
    shell:
        """
        salmon quant \
            --index {input.index} \
            --libType {params.libtype} \
            --mates1 {input.reads[0]} \
            --mates2 {input.reads[1]} \
            --threads {threads} \
            {params.extra} \
            --output {params.outdir} \
            > {log} 2>&1
        """


rule strandedness_salmon:
    # Salmon tallies every fragment's orientation whatever --libType it was
    # given, so this also catches a wrongly *declared* strandedness -- which
    # otherwise silently costs most of the reads, as Salmon discards fragments
    # incompatible with the declared type.
    input:
        SALMON_DIR + "/lib_format_counts.json",
    output:
        SALMON_DIR + "/strandedness.tsv",
    conda:
        PYTHON_ENV
    params:
        declared=lambda wildcards: SAMPLE_TABLE[wildcards.sample]["strandedness"],
        stranded=config["qc"]["stranded_threshold"],
        unstranded=config["qc"]["unstranded_threshold"],
    shell:
        """
        python {SCRIPTS}/strandedness.py salmon \
            --sample {wildcards.sample} \
            --lib-format-counts {input} \
            --declared {params.declared} \
            --stranded-threshold {params.stranded} \
            --unstranded-threshold {params.unstranded} \
            > {output}
        """
