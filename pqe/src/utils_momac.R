# Gene-list loading helpers for momac pseudotime scripts.
# Source this after utils.R in each _momac / _momac_smooth / _momac_hvg script.

REPO_ROOT <- "/home/unix/cchu/projects/ZmanR/pqe"

NAIVE_IGG   <- file.path(REPO_ROOT, "results/06/igg_hvg_spearman_df.naive.csv")
NAIVE_ATREM <- file.path(REPO_ROOT, "results/06/atrem_hvg_spearman_df.naive.csv")
SMOOTH_IGG  <- file.path(REPO_ROOT, "results/06/smoothed_igg_hvg_spearman_df.csv")
SMOOTH_ATREM <- file.path(REPO_ROOT, "results/06/smoothed_atrem_hvg_spearman_df.csv")

# Naive CSVs: genes as columns, metacells as rows, values are -log10(p_value).
# Returns intersection of genes with p < pval_thresh in at least one metacell.
load_naive_genes <- function(pval_thresh = 0.05) {
  threshold <- -log10(pval_thresh)
  igg   <- read.csv(NAIVE_IGG,   row.names = 1, check.names = FALSE)
  atrem <- read.csv(NAIVE_ATREM, row.names = 1, check.names = FALSE)
  igg_sig   <- colnames(igg)[apply(igg   > threshold, 2, any)]
  atrem_sig <- colnames(atrem)[apply(atrem > threshold, 2, any)]
  genes <- intersect(igg_sig, atrem_sig)
  cat("Naive gene intersection (p <", pval_thresh, "):", length(genes), "genes\n")
  genes
}

# Smoothed CSVs: genes as rows, columns include p_value.
# Returns intersection of genes with p_value < pval_thresh.
load_smooth_genes <- function(pval_thresh = 0.05) {
  igg_s   <- read.csv(SMOOTH_IGG,   row.names = 1)
  atrem_s <- read.csv(SMOOTH_ATREM, row.names = 1)
  igg_sig   <- rownames(igg_s)[igg_s$p_value   < pval_thresh]
  atrem_sig <- rownames(atrem_s)[atrem_s$p_value < pval_thresh]
  genes <- intersect(igg_sig, atrem_sig)
  cat("Smooth gene intersection (p <", pval_thresh, "):", length(genes), "genes\n")
  genes
}

ZMAN_S3_ATREM <- "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_aTREM2_Time.csv"
ZMAN_S3_IGG   <- "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_Isotype_Control_Time.csv"

# Table S3 CSVs: Gene, SpearmanCorrelation, Pvalue columns.
# Returns intersection of genes with Pvalue < pval_thresh in both conditions.
load_zman_s3_genes <- function(pval_thresh = 0.05) {
  atrem <- read.csv(ZMAN_S3_ATREM)
  igg   <- read.csv(ZMAN_S3_IGG)
  atrem_sig <- atrem$Gene[atrem$Pvalue < pval_thresh]
  igg_sig   <- igg$Gene[igg$Pvalue   < pval_thresh]
  genes <- intersect(atrem_sig, igg_sig)
  cat("Zman S3 gene intersection (p <", pval_thresh, "):", length(genes), "genes\n")
  genes
}

# HVG genes stored in the .rds companion to the h5ad.
load_hvg_genes <- function(h5ad_path) {
  rds_path <- sub("\\.h5ad$", ".rds", h5ad_path)
  saved <- readRDS(rds_path)
  genes <- saved$hvg_genes
  cat("HVG genes from adata.var:", length(genes), "genes\n")
  genes
}
