###############################################################################
# Count matrices -- one set per branch, side by side
#
# results/counts/star/     STAR's per-gene read counts: integers, uniquely
#                          mapped fragments only, htseq-count union model.
#                          Classic raw counts for DESeq2/edgeR/limma, and the
#                          matrix to use when counts must be traceable to
#                          reads in the BAM.
# results/counts/salmon/   Salmon's bias-corrected transcript estimates,
#                          aggregated to genes with tximport. Multimapping
#                          reads are apportioned rather than discarded, and
#                          gene_counts.rds carries the average-transcript-
#                          length offsets that correct for isoform switching
#                          -- the recommended input for differential
#                          expression (Soneson et al. 2015; Love et al. 2018).
#
# Both directories also hold samples.tsv (the colData: samplesheet metadata
# plus the strandedness that branch used) and gene_counts.rds, a
# SummarizedExperiment ready for DESeq2.
###############################################################################

wildcard_constraints:
    branch="star|salmon",

TXIMPORT_OUTPUTS = [
    "gene_counts.tsv",
    "gene_counts_length_scaled.tsv",
    "gene_tpm.tsv",
    "gene_lengths.tsv",
    "transcript_counts.tsv",
    "transcript_tpm.tsv",
    "gene_counts.rds",
    "transcript_counts.rds",
]

STAR_COUNT_OUTPUTS = [
    "gene_counts.tsv",
    "gene_counts_summary.tsv",
    "gene_counts.rds",
]


def strandedness_files(wildcards):
    if wildcards.branch == "star":
        return expand(STAR_DIR + "{sample}.strandedness.tsv", sample=SAMPLES)
    return expand(SALMON_DIR + "/strandedness.tsv", sample=SAMPLES)


rule strandedness_branch:
    input:
        strandedness_files,
    output:
        tsv=QC_DIR + "{branch}/strandedness.tsv",
        mqc=QC_DIR + "{branch}/strandedness_mqc.tsv",
    params:
        fail=config["qc"]["fail_on_strandedness_mismatch"],
    run:
        rows = [Path(f).read_text().splitlines() for f in input]
        header, body = rows[0][0], [r[1] for r in rows]
        Path(output.tsv).write_text("\n".join([header] + body) + "\n")
        # MultiQC custom content: a table in the report, mismatches visible.
        Path(output.mqc).write_text(
            f"# id: 'strandedness_{wildcards.branch}'\n"
            f"# section_name: 'Strandedness ({wildcards.branch})'\n"
            "# description: 'Declared vs. inferred library strandedness. "
            "MISMATCH means the declared value is probably wrong.'\n"
            "# plot_type: 'table'\n"
            + "\n".join([header] + body) + "\n"
        )
        status = header.split("\t").index("status")
        mismatches = [r.split("\t")[0] for r in body if r.split("\t")[status] == "MISMATCH"]
        if mismatches:
            message = (
                f"{wildcards.branch}: declared strandedness disagrees with the data for "
                f"{', '.join(mismatches)}. See {output.tsv}."
            )
            if params.fail:
                Path(output.tsv).unlink()
                Path(output.mqc).unlink()
                raise WorkflowError(message)
            print(f"WARNING: {message}", file=sys.stderr)


rule strandedness_compare:
    # Two independent inferences of the same property: when STAR and Salmon
    # disagree, look at the library before trusting either matrix.
    input:
        star=QC_DIR + "star/strandedness.tsv",
        salmon=QC_DIR + "salmon/strandedness.tsv",
    output:
        QC_DIR + "strandedness.tsv",
    run:
        def read(path):
            lines = Path(path).read_text().splitlines()
            cols = lines[0].split("\t")
            return {r["sample"]: r for r in (dict(zip(cols, l.split("\t"))) for l in lines[1:])}
        star, salmon = read(input.star), read(input.salmon)
        out = ["\t".join(["sample", "declared", "star_inferred", "salmon_inferred",
                          "star_forward_fraction", "salmon_forward_fraction", "agreement"])]
        for s in SAMPLES:
            a, b = star[s], salmon[s]
            agree = "yes" if a["inferred"] == b["inferred"] else "NO"
            out.append("\t".join([s, a["declared"], a["inferred"], b["inferred"],
                                  a["forward_fraction"], b["forward_fraction"], agree]))
        Path(output[0]).write_text("\n".join(out) + "\n")


rule samples_table:
    # colData of the SummarizedExperiments: every samplesheet metadata
    # column, how many FASTQ pairs were merged, and the strandedness this
    # branch declared, inferred and used.
    input:
        QC_DIR + "{branch}/strandedness.tsv",
    output:
        COUNTS_DIR + "{branch}/samples.tsv",
    run:
        lines = Path(input[0]).read_text().splitlines()
        cols = lines[0].split("\t")
        strand = {r["sample"]: r for r in (dict(zip(cols, l.split("\t"))) for l in lines[1:])}
        meta_cols = sorted({k for s in SAMPLES for k in SAMPLE_TABLE[s]["meta"]})
        header = ["sample", *meta_cols, "fastq_pairs",
                  "strandedness_declared", "strandedness_inferred", "strandedness"]
        out = ["\t".join(header)]
        for s in SAMPLES:
            entry = SAMPLE_TABLE[s]
            out.append("\t".join([
                s,
                *[entry["meta"].get(c, "") for c in meta_cols],
                str(len(entry["units"])),
                strand[s]["declared"],
                strand[s]["inferred"],
                strand[s]["effective"],
            ]))
        Path(output[0]).write_text("\n".join(out) + "\n")


rule tximport:
    input:
        quant=expand(SALMON_DIR + "/quant.sf", sample=SAMPLES),
        tx2gene=TX2GENE,
        samples=COUNTS_DIR + "salmon/samples.tsv",
    output:
        [COUNTS_DIR + "salmon/" + name for name in TXIMPORT_OUTPUTS],
    log:
        LOGS + "counts/tximport.log",
    conda:
        env("bioconductor")
    params:
        samples=SAMPLES,
        outdir=COUNTS_DIR + "salmon",
    resources:
        mem_mb=scaled(8000),
        runtime=scaled(60),
    script:
        "../scripts/tximport.R"


rule star_counts:
    input:
        counts=expand(STAR_DIR + "{sample}.ReadsPerGene.out.tab", sample=SAMPLES),
        strandedness=QC_DIR + "star/strandedness.tsv",
        tx2gene=TX2GENE,
    output:
        counts=COUNTS_DIR + "star/gene_counts.tsv",
        summary=COUNTS_DIR + "star/gene_counts_summary.tsv",
    conda:
        PYTHON_ENV
    params:
        samples=" ".join(SAMPLES),
    shell:
        """
        python {SCRIPTS}/star_counts.py \
            --samples {params.samples} \
            --counts {input.counts} \
            --strandedness {input.strandedness} \
            --tx2gene {input.tx2gene} \
            --out-counts {output.counts} \
            --out-summary {output.summary}
        """


rule star_counts_se:
    input:
        counts=COUNTS_DIR + "star/gene_counts.tsv",
        tx2gene=TX2GENE,
        samples=COUNTS_DIR + "star/samples.tsv",
    output:
        COUNTS_DIR + "star/gene_counts.rds",
    log:
        LOGS + "counts/star_counts_se.log",
    conda:
        env("bioconductor")
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    script:
        "../scripts/star_counts_se.R"


rule biotypes:
    # Where each branch's counts go, by gene biotype and chromosome: rRNA
    # and mitochondrial shares are the first things to check in any RNA-seq
    # library (failed rRNA depletion / poly-A selection, dying cells), and
    # spike-ins show their share when present.
    input:
        counts=lambda wildcards: COUNTS_DIR + f"{wildcards.branch}/gene_counts.tsv",
        tx2gene=TX2GENE,
    output:
        QC_DIR + "{branch}/biotypes_mqc.tsv",
    conda:
        PYTHON_ENV
    shell:
        """
        python {SCRIPTS}/biotypes.py \
            --counts {input.counts} \
            --tx2gene {input.tx2gene} \
            --branch {wildcards.branch} \
            > {output}
        """
