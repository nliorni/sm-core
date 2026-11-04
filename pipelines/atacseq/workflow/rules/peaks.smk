###############################################################################
# Peaks, consensus peak set, counts and signal tracks
#
#   macs3           per-sample peaks on whole fragments (BAMPE), blacklisted
#                   peaks removed
#   frip            fraction of each sample's fragments in its own peaks
#   consensus_peaks all samples' peaks merged; a merged region is kept when it
#                   was called in >= peaks.min_replicates samples of at least
#                   one group (peaks.consensus_group_by)
#   annotate_peaks  nearest gene, distance to its TSS, promoter/genic/distal
#   peak_counts     fragments per consensus peak per sample (featureCounts)
#   bigwig          CPM-normalised fragment coverage per sample, and averaged
#                   per group (bigwig.group_by)
###############################################################################

ATACQC_ENV = env("atacqc")
CONSENSUS_DIR = PEAKS_DIR + "consensus/"


rule macs3:
    # BAMPE: MACS3 piles up each pair's actual fragment, so no shift or
    # extension model is needed (and --nomodel is implied). --keep-dup all:
    # duplicates were already handled in final_bam.
    input:
        bam=BAM_DIR + "{sample}.bam",
        gsize=GENOME_SIZE,
        blacklist=BLACKLIST_BED,
    output:
        peaks=PEAKS_DIR + "{sample}/{sample}_peaks.narrowPeak",
        summits=PEAKS_DIR + "{sample}/{sample}_summits.bed",
        xls=PEAKS_DIR + "{sample}/{sample}_peaks.xls",
    log:
        LOGS + "macs3/{sample}.log",
    benchmark:
        BENCH + "macs3/{sample}.tsv"
    conda:
        env("macs3")
    resources:
        mem_mb=scaled(config["peaks"]["mem_mb"]),
        runtime=scaled(240),
    params:
        outdir=PEAKS_DIR + "{sample}",
        cutoff="-p" if config["peaks"]["cutoff"] == "pvalue" else "-q",
        threshold=config["peaks"]["threshold"],
        extra=config["peaks"]["extra"],
    shell:
        r"""
        # --tempdir: MACS3 otherwise spills to /tmp, often small on compute nodes.
        macs3 callpeak -t {input.bam} -f BAMPE -n {wildcards.sample} \
            --outdir {params.outdir} --tempdir {params.outdir} -g $(cat {input.gsize}) \
            {params.cutoff} {params.threshold} --keep-dup all {params.extra} > {log} 2>&1
        raw={params.outdir}/{wildcards.sample}_peaks.narrowPeak
        if [ -s {input.blacklist} ]; then
            # Not bedtools: a narrowPeak and a blacklist are both small, and
            # awk keeps the environment to MACS3 alone.
            before=$(wc -l < $raw)
            sort -k1,1 -k2,2n $raw | awk 'BEGIN {{ OFS = "\t" }}
                NR == FNR {{ n[$1]++; s[$1, n[$1]] = $2; e[$1, n[$1]] = $3; next }}
                {{ for (i = 1; i <= n[$1]; i++) if (s[$1, i] < $3 && e[$1, i] > $2) next; print }}' \
                {input.blacklist} - > $raw.tmp
            mv $raw.tmp $raw
            echo "blacklist: kept $(wc -l < $raw) of $before peaks" >> {log}
        fi
        """


rule frip:
    input:
        bam=BAM_DIR + "{sample}.bam",
        peaks=PEAKS_DIR + "{sample}/{sample}_peaks.narrowPeak",
    output:
        saf=temp(QC_DIR + "samples/{sample}/{sample}.peaks.saf"),
        counts=temp(QC_DIR + "samples/{sample}/{sample}.frip.counts"),
        summary=QC_DIR + "samples/{sample}/{sample}.frip.counts.summary",
    log:
        LOGS + "frip/{sample}.log",
    conda:
        env("subread")
    threads: 4
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    shell:
        # Fragments (pairs) overlapping the sample's own peaks, over all of
        # its fragments: ENCODE's FRiP. -O so a fragment touching two nearby
        # peaks still counts as in a peak.
        r"""
        awk 'BEGIN {{ OFS = "\t"; print "GeneID", "Chr", "Start", "End", "Strand" }}
             {{ print $4, $1, $2 + 1, $3, "+" }}' {input.peaks} > {output.saf}
        featureCounts -T {threads} -p --countReadPairs -B -C -O -F SAF \
            -a {output.saf} -o {output.counts} {input.bam} > {log} 2>&1
        """


rule consensus_peaks:
    input:
        peaks=expand(PEAKS_DIR + "{sample}/{sample}_peaks.narrowPeak", sample=SAMPLES),
        fai=GENOME_FAI,
    output:
        bed=CONSENSUS_DIR + "consensus_peaks.bed",
        saf=CONSENSUS_DIR + "consensus_peaks.saf",
        support=CONSENSUS_DIR + "consensus_peaks.support.tsv",
    log:
        LOGS + "peaks/consensus_peaks.log",
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        samples=SAMPLES,
        groups=CONSENSUS_GROUPS,
        min_replicates=config["peaks"]["min_replicates"],
    script:
        "../scripts/consensus_peaks.py"


rule annotate_peaks:
    input:
        bed=CONSENSUS_DIR + "consensus_peaks.bed",
        genes=GENES_TSV,
    output:
        CONSENSUS_DIR + "consensus_peaks.annotated.tsv",
    log:
        LOGS + "peaks/annotate_peaks.log",
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        annotation=config["annotation"],
    script:
        "../scripts/annotate_peaks.py"


rule peak_counts:
    input:
        bams=expand(BAM_DIR + "{sample}.bam", sample=SAMPLES),
        saf=CONSENSUS_DIR + "consensus_peaks.saf",
    output:
        raw=temp(COUNTS_DIR + "featurecounts.txt"),
        summary=COUNTS_DIR + "featurecounts.txt.summary",
    log:
        LOGS + "peak_counts/featurecounts.log",
    benchmark:
        BENCH + "peak_counts/featurecounts.tsv"
    conda:
        env("subread")
    threads: 8
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(120),
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs -B -C -O -F SAF \
            -a {input.saf} -o {output.raw} {input.bams} > {log} 2>&1
        """


rule count_matrix:
    # featureCounts' table, with samples as column names and the sample
    # metadata alongside, in the shape of rnaseq's gene count matrices.
    input:
        raw=COUNTS_DIR + "featurecounts.txt",
        bed=CONSENSUS_DIR + "consensus_peaks.bed",
    output:
        counts=COUNTS_DIR + "peak_counts.tsv",
        samples=COUNTS_DIR + "samples.tsv",
    log:
        LOGS + "peak_counts/count_matrix.log",
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        samples=SAMPLES,
        bams=expand(BAM_DIR + "{sample}.bam", sample=SAMPLES),
        meta={s: SAMPLE_TABLE[s]["meta"] for s in SAMPLES},
    script:
        "../scripts/count_matrix.py"


rule bigwig:
    input:
        bam=BAM_DIR + "{sample}.bam",
        bai=BAM_DIR + "{sample}.bam.bai",
        blacklist=BLACKLIST_BED,
    output:
        BIGWIG_DIR + "{sample}.bw",
    log:
        LOGS + "bigwig/{sample}.log",
    benchmark:
        BENCH + "bigwig/{sample}.tsv"
    conda:
        env("deeptools")
    threads: config["bigwig"]["threads"]
    resources:
        mem_mb=scaled(config["bigwig"]["mem_mb"]),
        runtime=scaled(120),
    params:
        bin_size=config["bigwig"]["bin_size"],
        blacklist="--blackListFileName " + BLACKLIST_BED if HAS_BLACKLIST else "",
        tmp=BIGWIG_DIR + "tmp_{sample}",
    shell:
        # --extendReads: each pair is drawn as its whole fragment, the unit
        # the peaks were called on. CPM makes tracks of different depths
        # comparable in a browser. TMPDIR: deepTools' temporary bedGraphs go
        # next to the output instead of /tmp.
        """
        mkdir -p {params.tmp}
        export TMPDIR={params.tmp}
        bamCoverage -b {input.bam} -o {output} -p {threads} \
            --binSize {params.bin_size} --normalizeUsing CPM --extendReads \
            {params.blacklist} > {log} 2>&1
        rm -rf {params.tmp}
        """


rule bigwig_group:
    input:
        lambda wildcards: expand(BIGWIG_DIR + "{sample}.bw", sample=BIGWIG_GROUPS[wildcards.group]),
    output:
        BIGWIG_DIR + "groups/{group}.bw",
    log:
        LOGS + "bigwig/group_{group}.log",
    conda:
        env("deeptools")
    threads: config["bigwig"]["threads"]
    resources:
        mem_mb=scaled(config["bigwig"]["mem_mb"]),
        runtime=scaled(120),
    params:
        bin_size=config["bigwig"]["bin_size"],
    shell:
        """
        bigwigAverage -b {input} -o {output} -p {threads} --binSize {params.bin_size} > {log} 2>&1
        """
