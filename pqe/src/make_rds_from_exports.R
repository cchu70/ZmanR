#!/usr/bin/env Rscript
# Create RDS companion files from Python-exported components.
library(Matrix)

make_rds_from_exports <- function(base_path) {
  counts_npz   <- paste0(base_path, "_counts.npz")  # legacy, unused
  gene_names_f <- paste0(base_path, "_gene_names.txt")
  hvg_genes_f  <- paste0(base_path, "_hvg_genes.txt")
  metadata_f   <- paste0(base_path, "_metadata.csv")
  rds_path     <- paste0(base_path, ".rds")

  cat("Loading exports for", base_path, "\n")

  gene_names <- readLines(gene_names_f)
  hvg_genes  <- readLines(hvg_genes_f)
  metadata   <- read.csv(metadata_f, row.names = 1, check.names = FALSE)

  # Load sparse matrix from Matrix Market format (genes x cells)
  counts_mtx <- paste0(base_path, "_counts.mtx")
  counts <- readMM(counts_mtx)
  rownames(counts) <- gene_names
  colnames(counts) <- rownames(metadata)

  cat("  Cells:", ncol(counts), "  Genes:", nrow(counts), "  HVGs:", length(hvg_genes), "\n")
  saved <- list(counts = counts, metadata = metadata, gene_names = gene_names, hvg_genes = hvg_genes)
  saveRDS(saved, rds_path)
  cat("  Saved:", rds_path, "\n")
}

RESULTS <- "/home/unix/cchu/projects/ZmanR/pqe/results/06"
make_rds_from_exports(file.path(RESULTS, "zmanseq_momac_metacells_annot_clean"))
make_rds_from_exports(file.path(RESULTS, "zmanseq_momac_metacells_annot_clean_mcumi"))
cat("Done.\n")
