#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_momac_zman_s3/redpath_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(redPATH)

genes <- load_zman_s3_genes()
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Filtering genes by dispersion...\n")
filtered <- filter_exp(dat$expr, dispersion_threshold = 1, threads = 4)
cat("Genes after dispersion filter:", ncol(filtered), "\n")
cat("Running redPATH...\n")
pseudotime <- redpath(filtered, threadnum = 4, base_path_range = c(3:7))
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "redPATH pseudotime (momac zman_s3)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pseudotime, sub("\\.png$", ".tsv", OUT_PLOT))
