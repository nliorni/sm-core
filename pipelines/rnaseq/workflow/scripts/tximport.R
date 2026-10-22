# Aggregate Salmon's per-sample transcript estimates into count matrices.
#
# Run by Snakemake's `script:` directive (rule tximport), which provides the
# `snakemake` object. Writes, into params$outdir:
#
#   transcript_counts.tsv / transcript_tpm.tsv   Salmon's estimates, as is
#   gene_counts.tsv              estimated counts summed per gene
#                                (countsFromAbundance = "no")
#   gene_lengths.tsv             per-sample average transcript length per gene:
#                                the offset that makes gene_counts.tsv correct
#                                for differential expression
#   gene_counts_length_scaled.tsv  counts rescaled from TPM to the library
#                                size and average length ("lengthScaledTPM"):
#                                usable as-is by edgeR/limma-voom, no offset
#   gene_tpm.tsv                 TPM summed per gene
#   gene_counts.rds              SummarizedExperiment: assays counts,
#                                abundance (TPM), avgTxLength,
#                                counts_length_scaled; rowData gene_name,
#                                gene_type, chrom; colData the sample table.
#                                DESeq2::DESeqDataSet(se, ~ condition) uses
#                                the avgTxLength offsets automatically.
#   transcript_counts.rds        the same at transcript level
#
# Soneson, Love & Robinson 2015 (F1000Research 4:1521) for why gene-level
# counts should be built this way rather than by counting reads per gene.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
    library(tximport)
    library(SummarizedExperiment)
})

outdir <- snakemake@params[["outdir"]]
samples <- unlist(snakemake@params[["samples"]])
quant_files <- setNames(unlist(snakemake@input[["quant"]]), samples)

coldata <- read.delim(snakemake@input[["samples"]], colClasses = "character", check.names = FALSE)
rownames(coldata) <- coldata$sample
coldata <- coldata[samples, , drop = FALSE]

tx2gene <- read.delim(snakemake@input[["tx2gene"]], colClasses = "character")

# dropInfReps: Salmon is run without --numBootstraps, so there are no
# inferential replicates to import (and reading them needs jsonlite).
txi_tx <- tximport(quant_files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)

# Salmon's index collapses transcripts with identical sequence, and spliced
# sequences gffread could not build are absent too: aggregate only what was
# quantified.
missing <- setdiff(rownames(txi_tx$counts), tx2gene$transcript_id)
if (length(missing) > 0) {
    stop(length(missing), " quantified transcripts are not in tx2gene.tsv, e.g. ",
         paste(head(missing, 3), collapse = ", "))
}
tx2gene <- tx2gene[tx2gene$transcript_id %in% rownames(txi_tx$counts), ]

gene <- summarizeToGene(txi_tx, tx2gene[, c("transcript_id", "gene_id")],
                        countsFromAbundance = "no")
gene_ls <- summarizeToGene(txi_tx, tx2gene[, c("transcript_id", "gene_id")],
                           countsFromAbundance = "lengthScaledTPM")

gene_annot <- unique(tx2gene[, c("gene_id", "gene_name", "gene_type", "chrom")])
gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), ]
rownames(gene_annot) <- gene_annot$gene_id
gene_annot <- gene_annot[rownames(gene$counts), ]

tx_annot <- tx2gene[, c("transcript_id", "gene_id", "gene_name", "gene_type")]
rownames(tx_annot) <- tx_annot$transcript_id
tx_annot <- tx_annot[rownames(txi_tx$counts), ]

write_matrix <- function(m, ids, file) {
    df <- data.frame(ids, m[rownames(ids), samples, drop = FALSE], check.names = FALSE)
    write.table(df, file.path(outdir, file), sep = "\t", quote = FALSE, row.names = FALSE)
}

gene_ids <- gene_annot[, c("gene_id", "gene_name")]
tx_ids <- tx_annot[, c("transcript_id", "gene_id")]

write_matrix(gene$counts, gene_ids, "gene_counts.tsv")
write_matrix(gene_ls$counts, gene_ids, "gene_counts_length_scaled.tsv")
write_matrix(gene$abundance, gene_ids, "gene_tpm.tsv")
write_matrix(gene$length, gene_ids, "gene_lengths.tsv")
write_matrix(txi_tx$counts, tx_ids, "transcript_counts.tsv")
write_matrix(txi_tx$abundance, tx_ids, "transcript_tpm.tsv")

provenance <- list(
    pipeline = "sm-core/rnaseq",
    branch = "salmon",
    tximport = as.character(packageVersion("tximport")),
    countsFromAbundance = "no"
)

# `counts` must be the first assay: DESeqDataSet() and most tools read
# assay(se) without naming it.
gene_se <- SummarizedExperiment(
    assays = list(
        counts = gene$counts[, samples, drop = FALSE],
        abundance = gene$abundance[, samples, drop = FALSE],
        avgTxLength = gene$length[, samples, drop = FALSE],
        counts_length_scaled = gene_ls$counts[, samples, drop = FALSE]
    ),
    rowData = gene_annot,
    colData = coldata,
    metadata = provenance
)
saveRDS(gene_se, file.path(outdir, "gene_counts.rds"))

tx_se <- SummarizedExperiment(
    assays = list(
        counts = txi_tx$counts[, samples, drop = FALSE],
        abundance = txi_tx$abundance[, samples, drop = FALSE],
        length = txi_tx$length[, samples, drop = FALSE]
    ),
    rowData = tx_annot,
    colData = coldata,
    metadata = provenance
)
saveRDS(tx_se, file.path(outdir, "transcript_counts.rds"))

message("genes: ", nrow(gene_se), "  transcripts: ", nrow(tx_se), "  samples: ", ncol(gene_se))
