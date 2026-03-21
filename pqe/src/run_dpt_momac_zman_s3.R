#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_momac_zman_s3/dpt_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(destiny)

genes <- load_zman_s3_genes()
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Building DiffusionMap...\n"); set.seed(42)
dm  <- DiffusionMap(dat$expr, n_pcs = min(30, ncol(dat$expr) - 1), n_eigs = 10, verbose = TRUE)
cat("Computing DPT...\n")
dpt <- DPT(dm, tips = 1L)
save_pseudotime_plot(dat$sc_x, dat$sc_y, dpt$DPT1, OUT_PLOT, "DPT pseudotime (momac zman_s3)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, dpt$DPT1, sub("\\.png$", ".tsv", OUT_PLOT))
