# Gene-list loading helpers for momac pseudotime scripts.
# Source this after utils.R in each _momac / _momac_smooth / _momac_hvg script.

REPO_ROOT <- "/home/unix/cchu/projects/ZmanR/pqe"

NAIVE_IGG   <- file.path(REPO_ROOT, "results/06/igg_hvg_spearman_df.naive.csv")
NAIVE_ATREM <- file.path(REPO_ROOT, "results/06/atrem_hvg_spearman_df.naive.csv")
SMOOTH_IGG  <- file.path(REPO_ROOT, "results/06/smoothed_igg_hvg_spearman_df.csv")
SMOOTH_ATREM <- file.path(REPO_ROOT, "results/06/smoothed_atrem_hvg_spearman_df.csv")

# Naive HVG CSVs: genes as rows, columns include p_value.
# Returns intersection of genes with p_value < pval_thresh.
load_naive_genes <- function(pval_thresh = 0.05) {
  igg   <- read.csv(NAIVE_IGG,   row.names = 1)
  atrem <- read.csv(NAIVE_ATREM, row.names = 1)
  igg_sig   <- rownames(igg)[igg$p_value   < pval_thresh]
  atrem_sig <- rownames(atrem)[atrem$p_value < pval_thresh]
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

FULL_IGG     <- file.path(REPO_ROOT, "results/06/igg_full_spearman_df.naive.csv")
FULL_ATREM   <- file.path(REPO_ROOT, "results/06/atrem_full_spearman_df.naive.csv")
SMOOTH_FULL_IGG   <- file.path(REPO_ROOT, "results/06/smoothed_igg_full_spearman_df.csv")
SMOOTH_FULL_ATREM <- file.path(REPO_ROOT, "results/06/smoothed_atrem_full_spearman_df.csv")

# Full naive CSVs: genes as rows, columns include p_value.
load_full_genes <- function(pval_thresh = 0.05) {
  igg   <- read.csv(FULL_IGG,   row.names = 1)
  atrem <- read.csv(FULL_ATREM, row.names = 1)
  igg_sig   <- rownames(igg)[igg$p_value   < pval_thresh]
  atrem_sig <- rownames(atrem)[atrem$p_value < pval_thresh]
  genes <- intersect(igg_sig, atrem_sig)
  cat("Full naive gene intersection (p <", pval_thresh, "):", length(genes), "genes\n")
  genes
}

# Smoothed full CSVs: genes as rows, columns include p_value.
load_smooth_full_genes <- function(pval_thresh = 0.05) {
  igg_s   <- read.csv(SMOOTH_FULL_IGG,   row.names = 1)
  atrem_s <- read.csv(SMOOTH_FULL_ATREM, row.names = 1)
  igg_sig   <- rownames(igg_s)[igg_s$p_value   < pval_thresh]
  atrem_sig <- rownames(atrem_s)[atrem_s$p_value < pval_thresh]
  genes <- intersect(igg_sig, atrem_sig)
  cat("Smooth full gene intersection (p <", pval_thresh, "):", length(genes), "genes\n")
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

CTET_MOS_SPEARMAN_DIR <- file.path(REPO_ROOT, "results/06/ctet_mos/spearman")

# Load intersection of significant genes from ctet_mos spearman tables.
# label: "scumi" or "mcumi"
# geneset: "hvg" or "full"
# timecol: "ctet" or "smooth_ctet"
load_ctet_mos_genes <- function(label, geneset = "hvg", timecol = "ctet", pval_thresh = 0.05) {
  f_atrem <- file.path(CTET_MOS_SPEARMAN_DIR, paste0(label, "_atrem_", geneset, "_", timecol, "_spearman.csv"))
  f_igg   <- file.path(CTET_MOS_SPEARMAN_DIR, paste0(label, "_igg_",   geneset, "_", timecol, "_spearman.csv"))
  df_atrem <- read.csv(f_atrem, row.names = 1)
  df_igg   <- read.csv(f_igg,   row.names = 1)
  sig_atrem <- rownames(df_atrem)[df_atrem$p_value < pval_thresh]
  sig_igg   <- rownames(df_igg)[df_igg$p_value   < pval_thresh]
  genes <- intersect(sig_atrem, sig_igg)
  cat("ctet_mos genes (", label, geneset, timecol, "p <", pval_thresh, "):", length(genes), "genes\n")
  genes
}
