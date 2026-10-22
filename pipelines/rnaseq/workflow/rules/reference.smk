###############################################################################
# Reference preparation -- built once, under reference.dir
#
#   genome.fa         the configured genome (+ spike-in sequences), plain FASTA
#   genes.gtf         the annotation, restricted to sequences in genome.fa
#                     (+ one gene per spike-in sequence)
#   transcripts.fa    spliced transcript sequences cut from genome.fa by
#                     genes.gtf -- so Salmon quantifies exactly the transcripts
#                     STAR's gene counts are built from, with identical IDs
#   tx2gene.tsv       transcript -> gene, with gene name, biotype, chromosome
#   genes.collapsed.gtf  one merged transcript per gene, for RNA-SeQC (star)
#   star/             STAR genome index                      (star branch)
#   salmon/           Salmon decoy-aware index               (salmon branch)
###############################################################################

PYTHON_ENV = env("python")
SAMTOOLS_ENV = env("samtools")


rule prepare_genome:
    # gzip -cdf passes an uncompressed file through unchanged, so the genome
    # and spike-in FASTA may each be plain or gzipped.
    input:
        genome=config["reference"]["genome"],
        spikein=config["reference"]["spikein_fasta"] or [],
    output:
        fa=GENOME_FA,
        fai=GENOME_FAI,
    log:
        LOGS + "reference/prepare_genome.log",
    conda:
        SAMTOOLS_ENV
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    shell:
        """
        gzip -cdf {input.genome} {input.spikein} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        dups=$(cut -f1 {output.fai} | sort | uniq -d | head -5)
        if [ -n "$dups" ]; then
            echo "Duplicate sequence names in genome + spike-ins: $dups" | tee -a {log} >&2
            exit 1
        fi
        """


rule prepare_annotation:
    input:
        gtf=config["reference"]["gtf"],
        fai=GENOME_FAI,
        spikein=config["reference"]["spikein_fasta"] or [],
    output:
        gtf=ANNOTATION_GTF,
        tx2gene=TX2GENE,
    log:
        LOGS + "reference/prepare_annotation.log",
    conda:
        PYTHON_ENV
    params:
        spikein=lambda wildcards, input: f"--spikein-fasta {input.spikein}" if input.spikein else "",
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(60),
    shell:
        """
        python {SCRIPTS}/gtf_tools.py prepare \
            --gtf {input.gtf} \
            --fai {input.fai} \
            {params.spikein} \
            --out-gtf {output.gtf} \
            --out-tx2gene {output.tx2gene} \
            2> {log}
        """


rule transcriptome:
    input:
        fa=GENOME_FA,
        fai=GENOME_FAI,
        gtf=ANNOTATION_GTF,
    output:
        TRANSCRIPTS_FA,
    log:
        LOGS + "reference/gffread.log",
    conda:
        env("gffread")
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(60),
    shell:
        "gffread -w {output} -g {input.fa} {input.gtf} 2> {log}"


rule star_index:
    input:
        fa=GENOME_FA,
        fai=GENOME_FAI,
        gtf=ANNOTATION_GTF,
    output:
        directory(REF_DIR + "star"),
    log:
        LOGS + "reference/star_index.log",
    benchmark:
        BENCH + "reference/star_index.tsv"
    conda:
        env("star")
    threads: config["reference"]["threads"]
    resources:
        mem_mb=scaled(config["reference"]["star_index_mem_mb"]),
        runtime=scaled(180),
    params:
        overhang=config["reference"]["star_sjdb_overhang"],
        tmp=REF_DIR + "star_tmp",
    shell:
        # --genomeSAindexNbases must shrink with the genome: STAR's own rule is
        # min(14, log2(genome length)/2 - 1), and the default 14 crashes or
        # hangs on small genomes (test data, spike-in-only references).
        # --limitGenomeGenerateRAM follows the job's memory, so a retry with
        # more memory (scaled()) really gets to use it.
        r"""
        sa=$(awk '{{ L += $2 }} END {{ n = int(log(L) / log(2) / 2 - 1); print (n < 14 ? n : 14) }}' {input.fai})
        # bioconda ships STAR as a launcher choosing a SIMD build (STAR-avx2,
        # STAR-avx, ...) at run time, and of these only STAR-plain is linked
        # with OpenMP -- which index generation uses to sort the suffix array
        # on all threads. Through the launcher a human index sorts on one core
        # (~6 h instead of ~40 min). Alignment threads differently and is
        # unaffected, so only this rule picks STAR-plain.
        star=$(command -v STAR-plain || command -v STAR)
        rm -rf {params.tmp}
        mkdir -p {output}
        $star --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --genomeSAindexNbases $sa \
            --limitGenomeGenerateRAM $(( {resources.mem_mb} * 1000000 )) \
            --outTmpDir {params.tmp} \
            --outFileNamePrefix {output}/ \
            > {log} 2>&1
        """


rule salmon_index:
    # Decoy-aware ("gentrome") index: the whole genome is indexed alongside
    # the transcripts as decoys, so a read that matches an unannotated or
    # intronic locus better than any transcript is discarded instead of being
    # forced onto a transcript it merely resembles (Srivastava et al. 2020).
    input:
        transcripts=TRANSCRIPTS_FA,
        fa=GENOME_FA,
        fai=GENOME_FAI,
        spikein=config["reference"]["spikein_fasta"] or [],
    output:
        directory(REF_DIR + "salmon"),
    log:
        LOGS + "reference/salmon_index.log",
    benchmark:
        BENCH + "reference/salmon_index.tsv"
    conda:
        env("salmon")
    threads: config["reference"]["threads"]
    resources:
        mem_mb=scaled(config["reference"]["salmon_index_mem_mb"]),
        runtime=scaled(180),
    params:
        k=config["reference"]["salmon_kmer"],
        gentrome=REF_DIR + "gentrome.fa",
        decoys=REF_DIR + "decoys.txt",
    shell:
        # Spike-in sequences are both a "chromosome" of genome.fa and a
        # transcript of transcripts.fa. They are left out of the decoys --
        # listing them twice would give two references one name -- and so must
        # also be left out of the genome half of the gentrome: Salmon requires
        # every decoy to come after every transcript, and the spike-ins sit at
        # the end of genome.fa, after the chromosomes.
        r"""
        cut -f1 {input.fai} > {params.decoys}
        if [ -n "{input.spikein}" ]; then
            gzip -cdf {input.spikein} | grep '^>' | sed 's/^>//; s/[[:space:]].*//' > {params.decoys}.spikein
            grep -vxFf {params.decoys}.spikein {params.decoys} > {params.decoys}.tmp || true
            mv {params.decoys}.tmp {params.decoys}
            rm -f {params.decoys}.spikein
        fi
        cp {input.transcripts} {params.gentrome}
        awk 'NR == FNR {{ decoy[$1] = 1; next }}
             /^>/ {{ name = substr($1, 2); keep = (name in decoy) }}
             keep' {params.decoys} {input.fa} >> {params.gentrome}
        salmon index \
            -t {params.gentrome} \
            -d {params.decoys} \
            -i {output} \
            -k {params.k} \
            -p {threads} \
            > {log} 2>&1
        rm -f {params.gentrome}
        """


rule collapse_annotation:
    input:
        ANNOTATION_GTF,
    output:
        COLLAPSED_GTF,
    log:
        LOGS + "reference/collapse_annotation.log",
    conda:
        PYTHON_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(60),
    shell:
        "python {SCRIPTS}/gtf_tools.py collapse --gtf {input} --out-gtf {output} 2> {log}"
