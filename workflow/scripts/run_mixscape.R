#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

runtime <- list(
  counts = snakemake@input[["counts"]],
  genes = snakemake@input[["genes"]],
  barcodes = snakemake@input[["barcodes"]],
  metadata = snakemake@input[["metadata"]],
  scores = snakemake@output[["scores"]],
  plots = snakemake@output[["plots"]],
  reference = snakemake@params[["reference"]],
  split_by = snakemake@params[["split_by"]],
  neighbors = as.integer(snakemake@params[["neighbors"]]),
  min_cells = as.integer(snakemake@params[["min_cells"]]),
  min_de_genes = as.integer(snakemake@params[["min_de_genes"]]),
  probability_threshold = as.numeric(snakemake@params[["probability_threshold"]])
)

counts <- readMM(runtime$counts)
rownames(counts) <- readLines(runtime$genes, warn = FALSE)
colnames(counts) <- readLines(runtime$barcodes, warn = FALSE)
metadata <- read.csv(runtime$metadata, stringsAsFactors = FALSE, check.names = FALSE)

required <- c(
  "cell_barcode", "target_gene", "assigned_guide", "is_non_targeting",
  "include_primary_analysis"
)
missing <- setdiff(required, colnames(metadata))
if (length(missing) > 0) stop("Metadata is missing columns: ", paste(missing, collapse = ", "))
if (!(runtime$reference %in% metadata$target_gene)) {
  stop("Mixscape reference class is absent: ", runtime$reference)
}

eligible <- metadata$include_primary_analysis %in% c(TRUE, "True", "true", "1", 1)
analysis_metadata <- metadata[eligible, , drop = FALSE]
analysis_cells <- intersect(analysis_metadata$cell_barcode, colnames(counts))
analysis_metadata <- analysis_metadata[match(analysis_cells, analysis_metadata$cell_barcode), , drop = FALSE]
if (length(analysis_cells) < 20) stop("Mixscape requires at least 20 eligible guide-singlet cells")
if (sum(analysis_metadata$target_gene == runtime$reference) < runtime$neighbors) {
  stop("Mixscape requires at least num.neighbors eligible non-targeting cells")
}

rownames(analysis_metadata) <- analysis_metadata$cell_barcode
object <- CreateSeuratObject(
  counts = counts[, analysis_cells, drop = FALSE],
  meta.data = analysis_metadata,
  min.cells = 0,
  min.features = 0
)
object <- NormalizeData(object, verbose = FALSE)
variable_count <- min(2000, nrow(object))
object <- FindVariableFeatures(object, nfeatures = variable_count, verbose = FALSE)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
npcs <- min(20, ncol(object) - 1, length(VariableFeatures(object)) - 1)
if (npcs < 2) stop("Not enough cells or variable genes to run Mixscape PCA")
object <- RunPCA(object, features = VariableFeatures(object), npcs = npcs, verbose = FALSE)

split_by <- runtime$split_by
if (!nzchar(split_by) || !(split_by %in% colnames(object[[]]))) split_by <- NULL
object <- CalcPerturbSig(
  object = object,
  assay = "RNA",
  slot = "data",
  gd.class = "target_gene",
  nt.cell.class = runtime$reference,
  split.by = split_by,
  num.neighbors = runtime$neighbors,
  reduction = "pca",
  ndims = npcs,
  new.assay.name = "PRTB",
  verbose = FALSE
)
DefaultAssay(object) <- "PRTB"
object <- ScaleData(object, do.scale = FALSE, do.center = TRUE, verbose = FALSE)
object <- RunMixscape(
  object = object,
  assay = "PRTB",
  slot = "scale.data",
  labels = "target_gene",
  nt.class.name = runtime$reference,
  new.class.name = "mixscape_class",
  min.de.genes = runtime$min_de_genes,
  min.cells = runtime$min_cells,
  de.assay = "RNA",
  split.by = split_by,
  fine.mode = TRUE,
  fine.mode.labels = "assigned_guide",
  prtb.type = "KO",
  verbose = FALSE
)

mixscape_metadata <- object[[]]
global_class_column <- if ("mixscape_class.global" %in% colnames(mixscape_metadata)) {
  "mixscape_class.global"
} else {
  "mixscape_class"
}
probability_columns <- grep("p_ko$", colnames(mixscape_metadata), value = TRUE)
score_columns <- grep("perturbation.*score", colnames(mixscape_metadata), value = TRUE, ignore.case = TRUE)

eligible_results <- data.frame(
  cell_barcode = rownames(mixscape_metadata),
  target_gene = mixscape_metadata$target_gene,
  assigned_guide = mixscape_metadata$assigned_guide,
  mixscape_class = as.character(mixscape_metadata[[global_class_column]]),
  perturbation_probability = if (length(probability_columns) > 0) {
    as.numeric(mixscape_metadata[[probability_columns[[1]]]])
  } else {
    NA_real_
  },
  perturbation_score = if (length(score_columns) > 0) {
    as.numeric(mixscape_metadata[[score_columns[[1]]]])
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

is_nt <- eligible_results$target_gene == runtime$reference
has_probability <- !is.na(eligible_results$perturbation_probability)
eligible_results$mixscape_class[is_nt] <- "NT"
eligible_results$mixscape_class[!is_nt & has_probability] <- ifelse(
  eligible_results$perturbation_probability[!is_nt & has_probability] >=
    runtime$probability_threshold,
  "KO",
  "NP"
)

results <- data.frame(
  cell_barcode = metadata$cell_barcode,
  target_gene = metadata$target_gene,
  assigned_guide = metadata$assigned_guide,
  mixscape_class = "excluded_low_moi",
  perturbation_probability = NA_real_,
  perturbation_score = NA_real_,
  stringsAsFactors = FALSE
)
matched <- match(results$cell_barcode, eligible_results$cell_barcode)
has_match <- !is.na(matched)
results[has_match, colnames(eligible_results)[-1]] <- eligible_results[matched[has_match], -1]
results$escaped_perturbation <- results$mixscape_class == "NP"

dir.create(dirname(runtime$scores), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(runtime$plots), recursive = TRUE, showWarnings = FALSE)
write.csv(results, runtime$scores, row.names = FALSE)

pdf(runtime$plots, width = 8, height = 5, onefile = TRUE)
for (target in setdiff(unique(eligible_results$target_gene), runtime$reference)) {
  target_rows <- eligible_results$target_gene == target
  values <- eligible_results$perturbation_probability[target_rows]
  values <- values[is.finite(values)]
  if (length(values) == 0) next
  hist(
    values,
    breaks = 20,
    col = "#2f6f8f",
    border = "white",
    xlim = c(0, 1),
    xlab = "Posterior perturbation probability",
    main = paste("Mixscape classification:", target)
  )
  abline(v = runtime$probability_threshold, col = "#b33a3a", lty = 2, lwd = 2)
}
dev.off()
