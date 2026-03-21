#!/usr/bin/env Rscript
# One-time: save metacell h5ad as RDS for envs where zellkonverter/basilisk fails.
library(zellkonverter)
library(SummarizedExperiment)
library(Matrix)

H5AD <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot_clean.h5ad"
OUT  <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.rds"

cat("Reading", H5AD, "\n")
sce <- readH5AD(H5AD)
cat("Shape:", nrow(sce), "genes x", ncol(sce), "cells\n")

hvg_genes <- rownames(sce)[rowData(sce)$highly_variable == TRUE]
cat("HVG genes:", length(hvg_genes), "\n")

saveRDS(list(
  counts     = assay(sce, "X"),
  metadata   = as.data.frame(colData(sce)),
  gene_names = rownames(sce),
  hvg_genes  = hvg_genes
), file = OUT)
cat("Saved to", OUT, "\n")
