#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

runtime <- list(
  counts = snakemake@input[["counts"]],
  design = snakemake@input[["design"]],
  results = snakemake@output[["results"]],
  plots = snakemake@output[["plots"]],
  summary = snakemake@output[["summary"]],
  reference = snakemake@params[["reference"]],
  covariates = snakemake@params[["covariates"]],
  strata = snakemake@params[["strata"]],
  min_replicates = as.integer(snakemake@params[["min_replicates"]])
)

split_csv <- function(value) {
  value <- trimws(value)
  if (!nzchar(value)) character(0) else trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

append_summary <- function(rows, stratum, target, n_target, n_reference, status, reason, formula) {
  rows[[length(rows) + 1]] <- data.frame(
    stratum = stratum,
    target_gene = target,
    reference = runtime$reference,
    n_target_samples = n_target,
    n_reference_samples = n_reference,
    status = status,
    reason = reason,
    model_formula = formula,
    stringsAsFactors = FALSE
  )
  rows
}

covariates <- split_csv(runtime$covariates)
strata <- split_csv(runtime$strata)

counts <- read.csv(runtime$counts, row.names = 1, check.names = FALSE)
sample_data <- read.csv(runtime$design, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- unique(c("pseudobulk_id", "target_gene", covariates, strata))
missing_columns <- setdiff(required_columns, colnames(sample_data))
if (length(missing_columns) > 0) {
  stop("Design matrix is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(sample_data$pseudobulk_id)) {
  stop("Design matrix contains duplicated pseudobulk_id values")
}
if (!all(colnames(counts) %in% sample_data$pseudobulk_id)) {
  stop("Every count-matrix column must have a matching pseudobulk_id")
}

sample_data <- sample_data[match(colnames(counts), sample_data$pseudobulk_id), , drop = FALSE]
if (!identical(colnames(counts), sample_data$pseudobulk_id)) {
  stop("Failed to align count columns with the design matrix")
}
if (anyNA(counts) || any(counts < 0) || any(abs(as.matrix(counts) - round(as.matrix(counts))) > 0)) {
  stop("edgeR requires non-negative integer counts without missing values")
}
counts <- as.matrix(round(counts))
storage.mode(counts) <- "integer"

if (!(runtime$reference %in% sample_data$target_gene)) {
  stop("Reference target is absent from the design matrix: ", runtime$reference)
}

if (length(strata) == 0) {
  sample_data$.stratum <- "all"
} else {
  sample_data$.stratum <- apply(sample_data[, strata, drop = FALSE], 1, paste, collapse = "|")
}

dir.create(dirname(runtime$results), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(runtime$plots), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(runtime$summary), recursive = TRUE, showWarnings = FALSE)

result_rows <- list()
summary_rows <- list()
pdf(runtime$plots, width = 6.5, height = 5.2, onefile = TRUE)

for (stratum_name in unique(sample_data$.stratum)) {
  stratum_idx <- which(sample_data$.stratum == stratum_name)
  stratum_design <- sample_data[stratum_idx, , drop = FALSE]
  targets <- setdiff(unique(stratum_design$target_gene), runtime$reference)

  for (target in targets) {
    keep_samples <- stratum_design$target_gene %in% c(runtime$reference, target)
    contrast_design <- stratum_design[keep_samples, , drop = FALSE]
    contrast_counts <- counts[, stratum_idx[keep_samples], drop = FALSE]
    n_target <- sum(contrast_design$target_gene == target)
    n_reference <- sum(contrast_design$target_gene == runtime$reference)

    if (n_target < runtime$min_replicates || n_reference < runtime$min_replicates) {
      summary_rows <- append_summary(
        summary_rows, stratum_name, target, n_target, n_reference, "skipped",
        "insufficient_replicates", ""
      )
      next
    }

    contrast_design$target_gene <- relevel(factor(contrast_design$target_gene), runtime$reference)
    usable_covariates <- covariates[
      vapply(contrast_design[, covariates, drop = FALSE], function(x) length(unique(x)) > 1, logical(1))
    ]

    formula <- reformulate(c(usable_covariates, "target_gene"))
    model <- model.matrix(formula, data = contrast_design)

    formula_text <- paste(deparse(formula), collapse = "")
    if (qr(model)$rank < ncol(model) || nrow(model) <= ncol(model)) {
      summary_rows <- append_summary(
        summary_rows, stratum_name, target, n_target, n_reference, "skipped",
        "rank_deficient_or_no_residual_df_no_covariates_dropped", formula_text
      )
      next
    }

    y <- DGEList(counts = contrast_counts)
    expressed <- filterByExpr(y, design = model)
    if (sum(expressed) < 2) {
      summary_rows <- append_summary(
        summary_rows, stratum_name, target, n_target, n_reference, "skipped",
        "no_genes_after_filtering", formula_text
      )
      next
    }

    y <- y[expressed, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, model, robust = TRUE)
    fit <- glmQLFit(y, model, robust = TRUE)
    target_coef <- grep("^target_gene", colnames(model))
    if (length(target_coef) != 1) {
      stop("Expected exactly one target coefficient for contrast: ", target)
    }
    qlf <- glmQLFTest(fit, coef = target_coef)
    table <- topTags(qlf, n = Inf, sort.by = "none")$table
    table$gene <- rownames(table)
    table$target_gene <- target
    table$reference <- runtime$reference
    table$stratum <- stratum_name
    table$n_target_samples <- n_target
    table$n_reference_samples <- n_reference
    table$method <- "edgeR_QLF"
    result_rows[[length(result_rows) + 1]] <- table[
      , c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "target_gene",
          "reference", "stratum", "n_target_samples", "n_reference_samples", "method")
    ]
    summary_rows <- append_summary(
      summary_rows, stratum_name, target, n_target, n_reference, "tested", "", formula_text
    )

    significant <- table$FDR < 0.05 & abs(table$logFC) >= 0.5
    colors <- ifelse(significant, "#b33a3a", "#6b7280")
    plot(
      table$logFC, -log10(pmax(table$PValue, .Machine$double.xmin)),
      pch = 16, cex = 0.65, col = colors,
      xlab = "log2 fold change", ylab = "-log10(P value)",
      main = paste(target, "vs", runtime$reference, "-", stratum_name)
    )
    abline(h = -log10(0.05), v = c(-0.5, 0.5), col = "#374151", lty = 2)
  }
}
dev.off()

if (length(result_rows) == 0) {
  stop("No differential-expression contrasts could be tested; inspect contrast summary settings")
}

results <- do.call(rbind, result_rows)
summary_table <- do.call(rbind, summary_rows)
rownames(results) <- NULL
rownames(summary_table) <- NULL
write.csv(results, runtime$results, row.names = FALSE)
write.csv(summary_table, runtime$summary, row.names = FALSE)
