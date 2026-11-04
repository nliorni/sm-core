###############################################################################
# Alignment and filtering
#
#   bowtie2_align   bowtie2 -> pairs that are properly paired, primary and
#                   above the MAPQ cut-off (ENCODE's filter), coordinate-sorted
#   markdup         Picard MarkDuplicates (flags; nothing removed yet)
#   final_bam       duplicates removed (align.remove_duplicates) and the
#                   mitochondrion dropped: the BAM every later step reads
#   bam_stats       samtools flagstat / idxstats / stats, for MultiQC
#
# The marked BAM, with its duplicates and mitochondrial reads still in, is
# what the library-complexity and mitochondrial-fraction QC needs; it is
# deleted once that QC has read it.
###############################################################################

# SAM flags: 4 unmapped, 8 mate unmapped, 256 secondary, 512 QC fail,
# 1024 duplicate, 2048 supplementary.
FLAGS_BAD = 4 + 8 + 256 + 512 + 2048          # 2828
FLAGS_BAD_DUP = FLAGS_BAD + 1024              # 3852


rule bowtie2_align:
    input:
        reads=analysis_reads,
        index=bowtie2_index_files(),
    output:
        bam=temp(BAM_DIR + "work/{sample}.filtered.bam"),
    log:
        LOGS + "bowtie2/{sample}.log",
    benchmark:
        BENCH + "bowtie2/{sample}.tsv"
    conda:
        env("bowtie2")
    threads: config["align"]["threads"]
    resources:
        mem_mb=scaled(config["align"]["mem_mb"]),
        runtime=scaled(480),
    params:
        index=bowtie2_prefix(),
        max_fragment=config["align"]["max_fragment"],
        mapq=config["align"]["min_mapq"],
        bad=FLAGS_BAD,
        extra=config["align"]["extra"],
        tmp=BAM_DIR + "work/{sample}.sort",
    shell:
        # --very-sensitive -X 2000: ENCODE's and the Omni-ATAC protocol's
        # settings; fragments up to 2 kb count as proper pairs, which keeps
        # the di- and tri-nucleosome fragments. bowtie2 writes the two mates
        # of a pair next to each other, which is what fixmate needs. The
        # first filter can drop one mate of a pair (say, a low-MAPQ one);
        # fixmate -r then clears its partner's proper-pair flag and the
        # second filter drops the orphan.
        r"""
        bowtie2 --very-sensitive -X {params.max_fragment} --threads {threads} \
            -x {params.index} -1 {input.reads[0]} -2 {input.reads[1]} \
            --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
            {params.extra} 2> {log} \
          | samtools view -u -f 2 -F {params.bad} -q {params.mapq} - \
          | samtools fixmate -r - - \
          | samtools view -u -f 2 -F {params.bad} - \
          | samtools sort -@ 4 -m 1G -T {params.tmp} -o {output.bam} - 2>> {log}
        """


rule markdup:
    input:
        BAM_DIR + "work/{sample}.filtered.bam",
    output:
        bam=temp(BAM_DIR + "work/{sample}.marked.bam"),
        metrics=QC_DIR + "samples/{sample}/{sample}.markdup.metrics",
    log:
        LOGS + "markdup/{sample}.log",
    benchmark:
        BENCH + "markdup/{sample}.tsv"
    conda:
        env("picard")
    resources:
        mem_mb=scaled(config["align"]["markdup_mem_mb"]),
        runtime=scaled(240),
    params:
        # The JVM's heap, leaving room for its own overhead within the job.
        heap=lambda wildcards, resources: int(resources.mem_mb * 0.8),
        tmp=BAM_DIR + "work/{sample}.picard_tmp",
    shell:
        """
        mkdir -p {params.tmp}
        picard -Xmx{params.heap}m MarkDuplicates \
            --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT \
            --TMP_DIR {params.tmp} > {log} 2>&1
        rm -rf {params.tmp}
        """


rule final_bam:
    # A proper pair has both mates on one sequence, so leaving out the whole
    # mitochondrion (-L nuclear.bed) never splits a pair.
    input:
        bam=BAM_DIR + "work/{sample}.marked.bam",
        nuclear=NUCLEAR_BED,
    output:
        bam=BAM_DIR + "{sample}.bam",
        bai=BAM_DIR + "{sample}.bam.bai",
    log:
        LOGS + "final_bam/{sample}.log",
    conda:
        SAMTOOLS_ENV
    threads: 4
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    params:
        flags=FLAGS_BAD_DUP if config["align"]["remove_duplicates"] else FLAGS_BAD,
    shell:
        """
        samtools view -b -@ {threads} -f 2 -F {params.flags} -L {input.nuclear} \
            -o {output.bam} {input.bam} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


rule bam_stats:
    input:
        bam=BAM_DIR + "{sample}.bam",
        bai=BAM_DIR + "{sample}.bam.bai",
    output:
        flagstat=QC_DIR + "samples/{sample}/{sample}.flagstat",
        idxstats=QC_DIR + "samples/{sample}/{sample}.idxstats",
        stats=QC_DIR + "samples/{sample}/{sample}.stats",
    log:
        LOGS + "bam_stats/{sample}.log",
    conda:
        SAMTOOLS_ENV
    threads: 2
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        samtools stats -@ {threads} {input.bam} > {output.stats} 2>> {log}
        """
