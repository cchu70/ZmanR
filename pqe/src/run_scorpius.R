#!/usr/bin/env Rscript
# Run SCORPIUS pseudotime on zmanseq data.
# Change ADATA_PATH to point to a different h5ad file.

ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime/scorpius_pseudotime.png"

.script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.script_dir, "utils.R"))
library(SCORPIUS)

# Load and preprocess — SCORPIUS uses MDS (O(n^2)), cap cells for speed
dat <- load_adata(ADATA_PATH, max_cells = 3000, n_hvg = 2000)

# Reduce to 3D space using Spearman distance, then infer linear trajectory
cat("Running SCORPIUS...\n")
space <- reduce_dimensionality(dat$expr, dist = "spearman", ndim = 3)
traj  <- infer_trajectory(space)

# traj$time is a [0, 1] pseudotime for each cell
save_pseudotime_plot(dat$sc_x, dat$sc_y, traj$time, OUT_PLOT, "SCORPIUS pseudotime")
save_pseudotime_tsv(dat$cell_names, dat$coldata, traj$time,
                    sub("\\.png$", ".tsv", OUT_PLOT))
