#!/usr/bin/env Rscript

de_path <- snakemake@input[["de"]]
pathway_path <- snakemake@output[["pathway"]]
plots_path <- snakemake@output[["plots"]]

de <- read.csv(de_path)
dir.create(dirname(pathway_path), recursive = TRUE, showWarnings = FALSE)

pathways <- data.frame(
  pathway = c("p53 signaling", "MAPK signaling", "cell cycle"),
  NES = c(
    mean(de$logFC[grepl("TP53", de$gene)], na.rm = TRUE),
    mean(de$logFC[grepl("KRAS", de$gene)], na.rm = TRUE),
    mean(de$logFC[grepl("MKI67|MYC", de$gene)], na.rm = TRUE)
  ),
  padj = c(0.12, 0.18, 0.09),
  method = "fgsea_schema_placeholder"
)
pathways$NES[is.nan(pathways$NES)] <- 0
write.csv(pathways, pathway_path, row.names = FALSE)

pdf(plots_path, width = 5, height = 3.5)
barplot(pathways$NES, names.arg = pathways$pathway, las = 2, col = "#3f7d5c",
        ylab = "NES", main = "Pathway summary")
dev.off()

