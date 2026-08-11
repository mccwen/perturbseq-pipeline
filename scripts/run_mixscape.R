#!/usr/bin/env Rscript

# Production hook for Seurat/Mixscape. The Python fallback documents the expected
# output schema and keeps the demo workflow lightweight.
suppressPackageStartupMessages({
  library(Seurat)
})

message("Load Seurat object, run CalcPerturbSig(), RunMixscape(), and export KO/NP/NT labels.")

