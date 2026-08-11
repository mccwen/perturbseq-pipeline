#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
})

counts_path <- snakemake@input[["counts"]]
design_path <- snakemake@input[["design"]]
results_path <- snakemake@output[["results"]]
plots_path <- snakemake@output[["plots"]]

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
design <- read.csv(design_path, check.names = FALSE)

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)

if (ncol(counts) < 2) {
  results <- data.frame(
    gene = rownames(counts),
    logFC = 0,
    logCPM = rowMeans(log2(counts + 1)),
    PValue = 1,
    FDR = 1,
    method = "placeholder_insufficient_replicates"
  )
} else {
  group <- as.factor(design$target_gene)
  baseline <- levels(group)[1]
  case <- levels(group)[min(2, length(levels(group)))]
  case_cols <- which(group == case)
  base_cols <- which(group == baseline)
  logfc <- rowMeans(log2(counts[, case_cols, drop = FALSE] + 1)) -
    rowMeans(log2(counts[, base_cols, drop = FALSE] + 1))
  pval <- apply(log2(counts + 1), 1, function(x) {
    if (length(unique(group)) < 2) return(1)
    tryCatch(t.test(x[group == case], x[group == baseline])$p.value, error = function(e) 1)
  })
  results <- data.frame(
    gene = rownames(counts),
    logFC = logfc,
    logCPM = rowMeans(log2(counts + 1)),
    PValue = pval,
    FDR = p.adjust(pval, method = "BH"),
    contrast = paste(case, "vs", baseline),
    method = "edgeR_QLF_schema_placeholder"
  )
}

write.csv(results, results_path, row.names = FALSE)

pdf(plots_path, width = 5, height = 4)
plot(results$logFC, -log10(results$PValue), pch = 16, col = "#2f6f8f",
     xlab = "logFC", ylab = "-log10(PValue)", main = "Differential expression")
abline(h = -log10(0.05), col = "#b33a3a", lty = 2)
dev.off()

