#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot_clean_mcumi.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_mcumi_hvg/scorpius_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(SCORPIUS)

genes <- load_hvg_genes(ADATA_PATH)
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Running SCORPIUS...\n")
space <- reduce_dimensionality(dat$expr, dist = "spearman", ndim = 3)
traj  <- infer_trajectory(space)
start_idx <- which.min(dat$coldata$cTET)
if (traj$time[start_idx] > 0.5) traj$time <- 1 - traj$time
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\n")
save_pseudotime_plot(dat$sc_x, dat$sc_y, traj$time, OUT_PLOT, "SCORPIUS pseudotime (mcumi hvg)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, traj$time, sub("[.]png$", ".tsv", OUT_PLOT))
