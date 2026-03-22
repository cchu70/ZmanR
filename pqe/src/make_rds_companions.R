#!/usr/bin/env Rscript
# Create RDS companion files for new adata h5ad objects.
# These cache the counts matrix + metadata + gene info for fast loading by run_*.R scripts.
library(zellkonverter); library(SummarizedExperiment); library(Matrix)

make_rds <- function(h5ad_path) {
  rds_path <- sub("\\.h5ad$", ".rds", h5ad_path)
  cat("Reading", h5ad_path, "...\n")
  sce <- readH5AD(h5ad_path)
  counts   <- assay(sce, "X")
  metadata <- as.data.frame(colData(sce))
  gene_names <- rownames(counts)
  hvg_genes  <- gene_names[rowData(sce)$highly_variable]
  cat("  Cells:", ncol(counts), "  Genes:", nrow(counts), "  HVGs:", length(hvg_genes), "\n")
  saved <- list(counts = counts, metadata = metadata, gene_names = gene_names, hvg_genes = hvg_genes)
  saveRDS(saved, rds_path)
  cat("  Saved:", rds_path, "\n")
}

RESULTS <- "/home/unix/cchu/projects/ZmanR/pqe/results/06"
make_rds(file.path(RESULTS, "zmanseq_momac_metacells_annot_clean.h5ad"))
make_rds(file.path(RESULTS, "zmanseq_momac_metacells_annot_clean_mcumi.h5ad"))
cat("Done.\n")
