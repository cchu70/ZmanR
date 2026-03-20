#!/usr/bin/env Rscript
# Run redPATH pseudotime on zmanseq data.
# Change ADATA_PATH to point to a different h5ad file.
# NOTE: redPATH uses a Hamiltonian-path approach designed for small datasets.
#       Cells are capped at 300 to keep runtime feasible.

ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime/redpath_pseudotime.png"

source(file.path(dirname(sys.frame(1)$ofile), "utils.R"))
library(redPATH)

# Load and preprocess — hard cap at 300 cells
dat <- load_adata(ADATA_PATH, max_cells = 300, n_hvg = 2000)

# redPATH expects cells × genes; further filter by dispersion
# filter_exp returns cell × gene with high-dispersion genes retained
cat("Filtering genes by dispersion...\n")
filtered <- filter_exp(dat$expr, dispersion_threshold = 1, threads = 4)
cat("Genes after dispersion filter:", ncol(filtered), "\n")

# Run redPATH — base_path_range should bracket estimated number of cell types
cat("Running redPATH...\n")
pseudotime <- redpath(filtered, threadnum = 4, base_path_range = c(3:7))

save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "redPATH pseudotime")
