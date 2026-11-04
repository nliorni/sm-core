###############################################################################
# QC
#
#   atac_qc         per sample: library complexity (NRF, PBC1, PBC2),
#                   mitochondrial fraction, fragment sizes, TSS enrichment,
#                   chrX / chrY fragments
#   somalier_*      sample identity (optional, qc.somalier_sites): genotypes
#                   at common SNPs, relatedness between every pair of samples
#                   -- including BAMs of other assays (qc.identity_bams) -- and
#                   inferred sex
#   qc_summary      every sample's metrics against the ENCODE ATAC-seq
#                   standards, as a table, figures and a MultiQC section
#   multiqc         fastp, bowtie2, Picard, samtools, MACS3, featureCounts,
#                   somalier and the tables above in one report
###############################################################################


rule atac_qc:
    input:
        bam=BAM_DIR + "work/{sample}.marked.bam",
        tss=TSS_BED,
        mito=MITO_NAME,
    output:
        metrics=QC_DIR + "samples/{sample}/{sample}.atac_qc.json",
        fragments=QC_DIR + "samples/{sample}/{sample}.fragment_sizes.tsv",
        tss=QC_DIR + "samples/{sample}/{sample}.tss_profile.tsv",
    log:
        LOGS + "atac_qc/{sample}.log",
    benchmark:
        BENCH + "atac_qc/{sample}.tsv"
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(config["qc"]["mem_mb"]),
        runtime=scaled(180),
    params:
        remove_duplicates=config["align"]["remove_duplicates"],
    script:
        "../scripts/atac_qc.py"


## ---------------------------------------------------------------- identity

IDENTITY_DIR = QC_DIR + "identity/"


rule somalier_extract:
    input:
        bam=BAM_DIR + "{sample}.bam",
        bai=BAM_DIR + "{sample}.bam.bai",
        fa=GENOME_FA,
        fai=GENOME_FAI,
        sites=config["qc"]["somalier_sites"] or [],
    output:
        IDENTITY_DIR + "extract/atac/{sample}.somalier",
    log:
        LOGS + "somalier/{sample}.log",
    conda:
        env("somalier")
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    params:
        outdir=IDENTITY_DIR + "extract/atac",
    shell:
        "somalier extract -d {params.outdir} -s {input.sites} -f {input.fa} {input.bam} > {log} 2>&1"


rule somalier_extract_other:
    # The digest is named after the BAM's read-group SM, or after the file
    # (minus .bam) when the BAM has no read group (somalier would stop). The
    # name may equal an ATAC sample's; the prefix keeps them apart, and
    # somalier relate strips it again to find pairs expected to be the same
    # cells.
    input:
        bam=lambda wildcards: IDENTITY_BAMS[int(wildcards.i)],
        fa=GENOME_FA,
        fai=GENOME_FAI,
        sites=config["qc"]["somalier_sites"] or [],
    output:
        directory(IDENTITY_DIR + "extract/other/{i}"),
    log:
        LOGS + "somalier/other_{i}.log",
    conda:
        env("somalier")
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    params:
        prefix=config["qc"]["identity_prefix"],
    shell:
        """
        mkdir -p {output}
        # awk reads the whole header: exiting early would end samtools with
        # SIGPIPE, which pipefail turns into a failure.
        name=$(samtools view -H {input.bam} | awk -F'\t' '!found && $1 == "@RG" {{
            for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) {{ print substr($i, 4); found = 1; break }} }}')
        export SOMALIER_SAMPLE_NAME=${{name:-$(basename {input.bam} .bam)}}
        echo "sample name: $SOMALIER_SAMPLE_NAME" > {log}
        somalier extract -d {output} -s {input.sites} -f {input.fa} \
            --sample-prefix {params.prefix} {input.bam} >> {log} 2>&1
        """


rule somalier_relate:
    input:
        atac=expand(IDENTITY_DIR + "extract/atac/{sample}.somalier", sample=SAMPLES),
        other=expand(IDENTITY_DIR + "extract/other/{i}", i=range(len(IDENTITY_BAMS))),
        sites=config["qc"]["somalier_sites"] or [],
    output:
        pairs=IDENTITY_DIR + "somalier.pairs.tsv",
        samples=IDENTITY_DIR + "somalier.samples.tsv",
        html=IDENTITY_DIR + "somalier.html",
    log:
        LOGS + "somalier/relate.log",
    conda:
        env("somalier")
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        prefix=IDENTITY_DIR + "somalier",
        sample_prefix=lambda wildcards: f"--sample-prefix {config['qc']['identity_prefix']}" if IDENTITY_BAMS else "",
        other=lambda wildcards, input: " ".join(f"{d}/*.somalier" for d in input.other),
        extra=config["qc"]["somalier_extra"],
    shell:
        """
        somalier relate -s {input.sites} -o {params.prefix} {params.sample_prefix} {params.extra} \
            {input.atac} {params.other} > {log} 2>&1
        """


rule identity_matches:
    # Every sample's best genetic match among the ATAC samples and among the
    # other assays' BAMs: the swap check that does not rely on names.
    input:
        pairs=IDENTITY_DIR + "somalier.pairs.tsv",
    output:
        IDENTITY_DIR + "identity_matches.tsv",
    log:
        LOGS + "somalier/identity_matches.log",
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(10),
    params:
        prefix=config["qc"]["identity_prefix"],
        has_other=bool(IDENTITY_BAMS),
    script:
        "../scripts/identity_matches.py"


## ---------------------------------------------------------------- summary

def per_sample(*patterns):
    return [path for pattern in patterns for path in expand(pattern, sample=SAMPLES)]


rule qc_summary:
    input:
        metrics=per_sample(QC_DIR + "samples/{sample}/{sample}.atac_qc.json"),
        fragments=per_sample(QC_DIR + "samples/{sample}/{sample}.fragment_sizes.tsv"),
        tss=per_sample(QC_DIR + "samples/{sample}/{sample}.tss_profile.tsv"),
        frip=per_sample(QC_DIR + "samples/{sample}/{sample}.frip.counts.summary"),
        markdup=per_sample(QC_DIR + "samples/{sample}/{sample}.markdup.metrics"),
        bowtie2=per_sample(LOGS + "bowtie2/{sample}.log"),
        peaks=per_sample(PEAKS_DIR + "{sample}/{sample}_peaks.narrowPeak"),
    output:
        table=QC_DIR + "atac_qc_summary.tsv",
        mqc_table=QC_DIR + "multiqc/atac_qc_mqc.tsv",
        mqc_fragments=QC_DIR + "multiqc/fragment_sizes_mqc.json",
        mqc_tss=QC_DIR + "multiqc/tss_enrichment_mqc.json",
        fig_metrics=QC_DIR + "qc_metrics.png",
        fig_fragments=QC_DIR + "fragment_sizes.png",
        fig_tss=QC_DIR + "tss_enrichment.png",
    log:
        LOGS + "qc/qc_summary.log",
    conda:
        ATACQC_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        smcore=SMCORE_ROOT,
        samples=SAMPLES,
        groups=CONSENSUS_GROUPS,
        remove_duplicates=config["align"]["remove_duplicates"],
    script:
        "../scripts/qc_summary.py"


def multiqc_inputs(wildcards):
    files = per_sample(RESULTS + "fastp/{sample}.fastp.json")
    files += per_sample(LOGS + "bowtie2/{sample}.log")
    files += per_sample(QC_DIR + "samples/{sample}/{sample}.markdup.metrics")
    files += per_sample(QC_DIR + "samples/{sample}/{sample}.flagstat")
    files += per_sample(QC_DIR + "samples/{sample}/{sample}.idxstats")
    files += per_sample(QC_DIR + "samples/{sample}/{sample}.stats")
    files += per_sample(PEAKS_DIR + "{sample}/{sample}_peaks.xls")
    files.append(COUNTS_DIR + "featurecounts.txt.summary")
    files += [
        QC_DIR + "multiqc/atac_qc_mqc.tsv",
        QC_DIR + "multiqc/fragment_sizes_mqc.json",
        QC_DIR + "multiqc/tss_enrichment_mqc.json",
    ]
    if IDENTITY:
        files += [IDENTITY_DIR + "somalier.pairs.tsv", IDENTITY_DIR + "somalier.samples.tsv"]
    return files


rule multiqc:
    input:
        multiqc_inputs,
    output:
        QC_DIR + "multiqc_report.html",
    log:
        LOGS + "multiqc/multiqc.log",
    conda:
        env("multiqc")
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        outdir=QC_DIR,
        config=str(PIPELINE_DIR / "workflow" / "multiqc_config.yaml"),
        extra=config["qc"]["multiqc_extra"],
        paths=lambda wildcards, input: "\n".join(sorted(input)),
    shell:
        """
        printf '%s\\n' "{params.paths}" > {params.outdir}/multiqc_files.txt
        multiqc \
            --file-list {params.outdir}/multiqc_files.txt \
            --config {params.config} \
            --title "sm-core/atacseq" \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            {params.extra} \
            > {log} 2>&1
        """
