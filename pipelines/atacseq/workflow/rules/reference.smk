###############################################################################
# Reference preparation -- built once, under reference.dir
#
#   genome.fa (+ .fai)          the configured genome, plain FASTA
#   chrom.sizes                 sequence lengths (bigWigs)
#   mito_name.txt               the mitochondrial sequence's name (may be empty)
#   nuclear.bed                 every sequence but the mitochondrion
#   effective_genome_size.txt   MACS3 -g: reference.effective_genome_size, or
#                               the genome's non-N length
#   genes.tsv                   genes with name, type, strand and TSS
#   tss.bed                     transcript TSSs for the TSS enrichment score
#   blacklist.bed               the configured blacklist, sorted, restricted to
#                               the genome's sequences (empty without one)
#   bowtie2/                    bowtie2 index
###############################################################################

PYTHON_ENV = env("python")
SAMTOOLS_ENV = env("samtools")

MITO_CANDIDATES = ("chrM", "MT", "chrMT", "M")


rule prepare_genome:
    # gzip -cdf passes an uncompressed file through unchanged.
    input:
        genome=config["reference"]["genome"],
    output:
        fa=GENOME_FA,
        fai=GENOME_FAI,
        sizes=CHROM_SIZES,
        mito=MITO_NAME,
        nuclear=NUCLEAR_BED,
        gsize=GENOME_SIZE,
    log:
        LOGS + "reference/prepare_genome.log",
    conda:
        SAMTOOLS_ENV
    resources:
        mem_mb=scaled(2000),
        runtime=scaled(60),
    params:
        mito=config["reference"]["mito_name"],
        candidates=" ".join(MITO_CANDIDATES),
        gsize=config["reference"]["effective_genome_size"] or "",
    shell:
        r"""
        gzip -cdf {input.genome} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        cut -f1,2 {output.fai} > {output.sizes}
        if [ "{params.mito}" = auto ]; then
            mito=""
            for name in {params.candidates}; do
                if cut -f1 {output.fai} | grep -qxF "$name"; then mito=$name; break; fi
            done
        else
            mito="{params.mito}"
            cut -f1 {output.fai} | grep -qxF "$mito" || {{
                echo "reference.mito_name: '$mito' is not a sequence of the genome" | tee -a {log} >&2; exit 1; }}
        fi
        echo "$mito" > {output.mito}
        echo "mitochondrial sequence: ${{mito:-none found}}" >> {log}
        awk -v m="$mito" 'BEGIN {{ OFS = "\t" }} $1 != m {{ print $1, 0, $2 }}' {output.fai} > {output.nuclear}
        if [ -n "{params.gsize}" ]; then
            echo {params.gsize} > {output.gsize}
        else
            # Non-N bases: what MACS3 calls the effective (mappable) genome
            # size, slightly generous since repeats are counted too.
            grep -v '^>' {output.fa} | tr -d 'NnRYKMSWBDHVrykmswbdhv\n' | wc -c > {output.gsize}
        fi
        echo "effective genome size: $(cat {output.gsize})" >> {log}
        """


rule prepare_genes:
    # Genes and their TSSs, for the peak annotation. The script is an input so
    # that a change to it rebuilds an existing reference's genes.tsv (a shell
    # call is not tracked as code); only the cheap annotation steps follow.
    input:
        gtf=config["reference"]["gtf"],
        fai=GENOME_FAI,
        script=SCRIPTS + "/gtf_annotation.py",
    output:
        GENES_TSV,
    log:
        LOGS + "reference/prepare_genes.log",
    conda:
        PYTHON_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(60),
    shell:
        "python {input.script} --gtf {input.gtf} --fai {input.fai} --out-genes {output} 2> {log}"


rule prepare_tss:
    # The TSS set of the enrichment score. A separate rule, so that changing
    # the gene annotation does not rerun the QC (which reads the pre-filter
    # BAMs and would realign every sample).
    input:
        gtf=config["reference"]["gtf"],
        fai=GENOME_FAI,
    output:
        TSS_BED,
    log:
        LOGS + "reference/prepare_tss.log",
    conda:
        PYTHON_ENV
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(60),
    shell:
        "python {SCRIPTS}/gtf_annotation.py --gtf {input.gtf} --fai {input.fai} --out-tss {output} 2> {log}"


rule prepare_blacklist:
    # Sequence names must be the genome's: an Ensembl-named genome ("1", "MT")
    # with a UCSC-named blacklist ("chr1", "chrM") would silently filter
    # nothing, so a blacklist sharing no name with the genome stops the run.
    input:
        bed=config["reference"]["blacklist"] or [],
        fai=GENOME_FAI,
    output:
        BLACKLIST_BED,
    log:
        LOGS + "reference/prepare_blacklist.log",
    resources:
        mem_mb=scaled(1000),
        runtime=scaled(10),
    shell:
        r"""
        if [ -z "{input.bed}" ]; then
            : > {output}
            echo "no blacklist configured" > {log}
            exit 0
        fi
        gzip -cdf {input.bed} \
            | awk 'NR == FNR {{ keep[$1] = 1; next }} ($1 in keep) {{ print $1 "\t" $2 "\t" $3 }}' {input.fai} - \
            | sort -k1,1 -k2,2n > {output}
        n=$(wc -l < {output})
        total=$(gzip -cdf {input.bed} | grep -vc '^#' || true)
        echo "blacklist: $n of $total regions on genome sequences" > {log}
        if [ "$n" -eq 0 ]; then
            echo "reference.blacklist shares no sequence name with the genome (e.g. chr1 vs 1)." | tee -a {log} >&2
            exit 1
        fi
        """


rule bowtie2_index:
    input:
        fa=GENOME_FA,
    output:
        multiext(REF_DIR + "bowtie2/genome", *BOWTIE2_SUFFIXES),
    log:
        LOGS + "reference/bowtie2_index.log",
    benchmark:
        BENCH + "reference/bowtie2_index.tsv"
    conda:
        env("bowtie2")
    threads: config["reference"]["threads"]
    resources:
        mem_mb=scaled(config["reference"]["bowtie2_index_mem_mb"]),
        runtime=scaled(360),
    params:
        prefix=REF_DIR + "bowtie2/genome",
    shell:
        "bowtie2-build --threads {threads} {input.fa} {params.prefix} > {log} 2>&1"
