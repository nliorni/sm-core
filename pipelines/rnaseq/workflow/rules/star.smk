###############################################################################
# STAR branch: genome alignment and per-gene read counts
#
#   <sample>.bam             coordinate-sorted genome alignment of every read
#                            (unmapped included): IGV, RNA-SeQC, splicing,
#                            later variant work
#   ReadsPerGene.out.tab     per-gene fragment counts, htseq-count union model,
#                            for all three strandedness options at once
#   SJ.out.tab               splice junctions found
#   Log.final.out            mapping statistics (MultiQC)
#
# Independent of the Salmon branch: nothing here reads its output, and the
# strandedness used to pick the count column is inferred from STAR's own
# counts (strandedness_star).
###############################################################################

STAR_DIR = RESULTS + "star/{sample}/"


rule star_align:
    input:
        reads=analysis_reads,
        index=star_index(),
    output:
        bam=STAR_DIR + "{sample}.bam",
        gene_counts=STAR_DIR + "{sample}.ReadsPerGene.out.tab",
        junctions=STAR_DIR + "{sample}.SJ.out.tab",
        log_final=STAR_DIR + "{sample}.Log.final.out",
    log:
        LOGS + "star/{sample}.log",
    benchmark:
        BENCH + "star/{sample}.tsv"
    conda:
        env("star")
    threads: config["star"]["threads"]
    resources:
        mem_mb=scaled(config["star"]["mem_mb"]),
        runtime=scaled(240),
    params:
        prefix=STAR_DIR + "{sample}.",
        tmp=STAR_DIR + "_STARtmp",
        read_cmd=lambda wildcards, input: decompress_cmd(input.reads[0]),
        two_pass="--twopassMode Basic" if config["star"]["two_pass"] else "",
        sort_threads=config["star"]["sort_threads"],
        sort_mem=config["star"]["sort_mem_per_thread"],
        extra=config["star"]["extra"],
    shell:
        # ENCODE long-RNA settings (also nf-core/rnaseq's): at most 20 loci per
        # read, mismatches capped at 4% of the read pair, introns 20 bp-1 Mb,
        # junctions not in the annotation need 8 bp of overhang and annotated
        # ones 1. Two-pass mode re-maps every read with the junctions found in
        # pass one, which recovers reads spanning junctions the GTF lacks.
        #
        # --runRNGseed 0 makes multimapper placement reproducible. The genome
        # BAM is streamed unsorted to samtools for sorting rather than sorted
        # by STAR, whose in-memory sort is the usual cause of its
        # out-of-memory failures.
        r"""
        rm -rf {params.tmp}
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.reads[0]} {input.reads[1]} \
            --readFilesCommand {params.read_cmd} \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp} \
            --outSAMtype BAM Unsorted \
            --outStd BAM_Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS NM MD \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            --outSAMstrandField intronMotif \
            --quantMode GeneCounts \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --runRNGseed 0 \
            {params.two_pass} \
            {params.extra} \
            2> {log} \
          | samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_mem} \
            -T {params.prefix}sort \
            -o {output.bam} \
            - 2>> {log}
        # STAR already wrote ReadsPerGene.out.tab, SJ.out.tab and Log.final.out
        # under the prefix, i.e. at the declared output paths; only its
        # scratch files are left to clean up.
        cat {params.prefix}Log.out >> {log}
        rm -rf {params.tmp} {params.prefix}Log.out {params.prefix}Log.progress.out \
               {params.prefix}_STARgenome {params.prefix}_STARpass1
        """


rule samtools_index:
    input:
        STAR_DIR + "{sample}.bam",
    output:
        STAR_DIR + "{sample}.bam.bai",
    log:
        LOGS + "samtools_index/{sample}.log",
    conda:
        SAMTOOLS_ENV
    threads: 4
    resources:
        mem_mb=scaled(1000),
        runtime=scaled(30),
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule strandedness_star:
    # STAR writes counts for all three library types at once; which column
    # holds most of the counts says how the library was made.
    input:
        STAR_DIR + "{sample}.ReadsPerGene.out.tab",
    output:
        STAR_DIR + "{sample}.strandedness.tsv",
    conda:
        PYTHON_ENV
    params:
        declared=lambda wildcards: SAMPLE_TABLE[wildcards.sample]["strandedness"],
        stranded=config["qc"]["stranded_threshold"],
        unstranded=config["qc"]["unstranded_threshold"],
    shell:
        """
        python {SCRIPTS}/strandedness.py star \
            --sample {wildcards.sample} \
            --reads-per-gene {input} \
            --declared {params.declared} \
            --stranded-threshold {params.stranded} \
            --unstranded-threshold {params.unstranded} \
            > {output}
        """
