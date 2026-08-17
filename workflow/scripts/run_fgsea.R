#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fgsea)
})

runtime <- list(
  de = snakemake@input[["de"]],
  gmt = snakemake@input[["gmt"]],
  results = snakemake@output[["pathway"]],
  plots = snakemake@output[["plots"]],
  min_size = as.integer(snakemake@params[["min_size"]]),
  max_size = as.integer(snakemake@params[["max_size"]]),
  random_seed = as.integer(snakemake@params[["random_seed"]])
)

set.seed(runtime$random_seed)
de <- read.csv(runtime$de, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("gene", "logFC", "F", "target_gene", "reference", "stratum")
missing <- setdiff(required, colnames(de))
if (length(missing) > 0) stop("DE results are missing columns: ", paste(missing, collapse = ", "))

pathways <- gmtPathways(runtime$gmt)
if (length(pathways) == 0) stop("No pathways found in GMT file: ", runtime$gmt)

de$.contrast <- paste(de$target_gene, de$reference, de$stratum, sep = "|")
rows <- list()
dir.create(dirname(runtime$results), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(runtime$plots), recursive = TRUE, showWarnings = FALSE)
pdf(runtime$plots, width = 8, height = 5.5, onefile = TRUE)

for (contrast_name in unique(de$.contrast)) {
  contrast <- de[de$.contrast == contrast_name, , drop = FALSE]
  contrast <- contrast[order(abs(contrast$logFC), decreasing = TRUE), , drop = FALSE]
  contrast <- contrast[!duplicated(contrast$gene), , drop = FALSE]
  ranks <- sign(contrast$logFC) * sqrt(pmax(contrast$F, 0))
  names(ranks) <- contrast$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  if (length(ranks) < runtime$min_size) next

  enrichment <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = runtime$min_size,
    maxSize = runtime$max_size,
    eps = 0
  )
  if (nrow(enrichment) == 0) next

  enrichment <- as.data.frame(enrichment)
  enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, collapse = ";", character(1))
  enrichment$target_gene <- contrast$target_gene[[1]]
  enrichment$reference <- contrast$reference[[1]]
  enrichment$stratum <- contrast$stratum[[1]]
  enrichment$method <- "fgseaMultilevel"
  rows[[length(rows) + 1]] <- enrichment

  top <- enrichment[order(enrichment$padj, -abs(enrichment$NES)), , drop = FALSE]
  top <- head(top, 12)
  top <- top[order(top$NES), , drop = FALSE]
  colors <- ifelse(top$NES >= 0, "#b33a3a", "#2f6f8f")
  par(mar = c(5, 12, 4, 2))
  barplot(
    top$NES,
    names.arg = top$pathway,
    horiz = TRUE,
    las = 1,
    col = colors,
    border = NA,
    xlab = "Normalized enrichment score",
    main = paste(contrast$target_gene[[1]], "vs", contrast$reference[[1]], "-", contrast$stratum[[1]])
  )
  abline(v = 0, col = "#374151")
}
dev.off()

if (length(rows) == 0) {
  stop("No pathways passed the configured size filters for any DE contrast")
}

results <- do.call(rbind, rows)
rownames(results) <- NULL
write.csv(results, runtime$results, row.names = FALSE)
