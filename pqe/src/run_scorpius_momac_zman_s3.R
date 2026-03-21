#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_momac_zman_s3/scorpius_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(SCORPIUS)

genes <- load_zman_s3_genes()
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Running SCORPIUS...\n")
space <- reduce_dimensionality(dat$expr, dist = "spearman", ndim = 3)
traj  <- infer_trajectory(space)
save_pseudotime_plot(dat$sc_x, dat$sc_y, traj$time, OUT_PLOT, "SCORPIUS pseudotime (momac zman_s3)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, traj$time, sub("\\.png$", ".tsv", OUT_PLOT))
