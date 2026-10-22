###############################################################################
# QC
#
# rnaseqc     RNA-SeQC 2 on the STAR BAM -- the GTEx/ENCODE metrics: exonic,
#             intronic and intergenic rates (genomic DNA contamination,
#             pre-mRNA), rRNA rate, 3'/5' coverage bias (RNA degradation),
#             genes detected, duplicate rate, insert size.
# multiqc     one report per branch, so each branch's QC stands alone:
#               star:   fastp, STAR mapping, STAR gene counts, RNA-SeQC,
#                       strandedness, count composition
#               salmon: fastp, Salmon mapping/fragment lengths,
#                       strandedness, count composition
###############################################################################


rule rnaseqc:
    input:
        bam=STAR_DIR + "{sample}.bam",
        bai=STAR_DIR + "{sample}.bam.bai",
        gtf=COLLAPSED_GTF,
        strandedness=STAR_DIR + "{sample}.strandedness.tsv",
    output:
        metrics=RESULTS + "rnaseqc/{sample}/{sample}.metrics.tsv",
    log:
        LOGS + "rnaseqc/{sample}.log",
    benchmark:
        BENCH + "rnaseqc/{sample}.tsv"
    conda:
        env("rnaseqc")
    resources:
        mem_mb=scaled(8000),
        runtime=scaled(120),
    params:
        outdir=RESULTS + "rnaseqc/{sample}",
        extra=config["qc"]["rnaseqc_extra"],
    shell:
        # The library type comes from the STAR branch's own strandedness call
        # (column 5, `effective`): RF = read 1 antisense (dUTP), FR = read 1
        # sense. Unstranded libraries get no flag.
        r"""
        strand=$(awk -F'\t' 'NR == 2 {{ print $5 }}' {input.strandedness})
        case "$strand" in
            reverse) flag="--stranded=RF" ;;
            forward) flag="--stranded=FR" ;;
            *)       flag="" ;;
        esac
        rnaseqc {input.gtf} {input.bam} {params.outdir} \
            --sample {wildcards.sample} \
            $flag \
            {params.extra} \
            > {log} 2>&1
        """


def _multiqc_path(f):
    if f.endswith("aux_info/meta_info.json"):
        return str(Path(f).parent.parent)
    if f.endswith("lib_format_counts.json"):
        return str(Path(f).parent)
    return f


def multiqc_inputs(wildcards):
    files = expand(RESULTS + "fastp/{sample}.fastp.json", sample=SAMPLES)
    files += [QC_DIR + f"{wildcards.branch}/strandedness_mqc.tsv",
              QC_DIR + f"{wildcards.branch}/biotypes_mqc.tsv"]
    if wildcards.branch == "star":
        files += expand(STAR_DIR + "{sample}.Log.final.out", sample=SAMPLES)
        files += expand(STAR_DIR + "{sample}.ReadsPerGene.out.tab", sample=SAMPLES)
        if config["qc"]["rnaseqc"]:
            files += expand(RESULTS + "rnaseqc/{sample}/{sample}.metrics.tsv", sample=SAMPLES)
    else:
        files += expand(SALMON_DIR + "/aux_info/meta_info.json", sample=SAMPLES)
        files += expand(SALMON_DIR + "/lib_format_counts.json", sample=SAMPLES)
    return files


rule multiqc:
    input:
        multiqc_inputs,
    output:
        QC_DIR + "{branch}/multiqc_report.html",
    log:
        LOGS + "multiqc/{branch}.log",
    conda:
        env("multiqc")
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    params:
        outdir=QC_DIR + "{branch}",
        config=str(PIPELINE_DIR / "workflow" / "multiqc_config.yaml"),
        extra=config["qc"]["multiqc_extra"],
        # Salmon's module reads more than meta_info.json (fragment-length
        # distribution, library format); hand it each sample's directory.
        paths=lambda wildcards, input: "\n".join(sorted({_multiqc_path(f) for f in input})),
    shell:
        """
        printf '%s\\n' "{params.paths}" > {params.outdir}/multiqc_files.txt
        multiqc \
            --file-list {params.outdir}/multiqc_files.txt \
            --config {params.config} \
            --title "sm-core/rnaseq: {wildcards.branch} branch" \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            {params.extra} \
            > {log} 2>&1
        """
