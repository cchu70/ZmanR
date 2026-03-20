#!/usr/bin/env Rscript
# Run DPT (destiny) pseudotime on zmanseq data.
# Change ADATA_PATH to point to a different h5ad file.

ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime/dpt_pseudotime.png"

source(file.path(dirname(sys.frame(1)$ofile), "utils.R"))
library(destiny)

# Load and preprocess
dat <- load_adata(ADATA_PATH, max_cells = 5000, n_hvg = 2000)

# Build diffusion map from cell × gene expression matrix
cat("Building DiffusionMap...\n")
set.seed(42)
dm  <- DiffusionMap(dat$expr, n_pcs = 30, n_eigs = 10, verbose = TRUE)

# Compute DPT from a single random root tip
cat("Computing DPT...\n")
dpt <- DPT(dm, tips = 1L)

# dpt$DPT1 holds the pseudotime from the first tip
save_pseudotime_plot(dat$sc_x, dat$sc_y, dpt$DPT1, OUT_PLOT, "DPT (destiny) pseudotime")
