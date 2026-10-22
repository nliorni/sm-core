# Wrap the STAR branch's gene count matrix in a SummarizedExperiment.
#
# Run by Snakemake's `script:` directive (rule star_counts_se). The same shape
# as the Salmon branch's gene_counts.rds -- rowData gene_name, gene_type,
# chrom; colData the sample table -- but with a single integer `counts` assay
# and no length offsets: these are reads counted per gene, so
# DESeq2::DESeqDataSet(se, ~ condition) takes them as they are.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages(library(SummarizedExperiment))

counts <- read.delim(snakemake@input[["counts"]], check.names = FALSE, colClasses = c("character", "character"))
rownames(counts) <- counts$gene_id
samples <- setdiff(colnames(counts), c("gene_id", "gene_name"))
m <- as.matrix(counts[, samples, drop = FALSE])
storage.mode(m) <- "integer"

tx2gene <- read.delim(snakemake@input[["tx2gene"]], colClasses = "character")
gene_annot <- unique(tx2gene[, c("gene_id", "gene_name", "gene_type", "chrom")])
gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), ]
rownames(gene_annot) <- gene_annot$gene_id
missing <- setdiff(rownames(m), rownames(gene_annot))
if (length(missing) > 0) {
    stop(length(missing), " genes counted by STAR are not in tx2gene.tsv, e.g. ",
         paste(head(missing, 3), collapse = ", "),
         " -- was the STAR index built from a different GTF?")
}

coldata <- read.delim(snakemake@input[["samples"]], colClasses = "character", check.names = FALSE)
rownames(coldata) <- coldata$sample

se <- SummarizedExperiment(
    assays = list(counts = m),
    rowData = gene_annot[rownames(m), ],
    colData = coldata[samples, , drop = FALSE],
    metadata = list(pipeline = "sm-core/rnaseq", branch = "star",
                    counting = "STAR --quantMode GeneCounts (htseq-count union model)")
)
saveRDS(se, snakemake@output[[1]])
message("genes: ", nrow(se), "  samples: ", ncol(se))
