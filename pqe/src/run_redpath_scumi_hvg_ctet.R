#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot_clean.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_scumi_hvg_ctet/redpath_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(redPATH)

genes <- load_ctet_mos_genes("scumi", geneset="hvg", timecol="ctet")
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Genes:", ncol(dat$expr), "\n")
cat("Running redPATH...\n")
pseudotime <- redpath(dat$expr, threadnum = 4, base_path_range = c(2:6))
start_idx <- which.min(dat$coldata$cTET)
if (pseudotime[start_idx] > median(pseudotime, na.rm = TRUE)) pseudotime <- max(pseudotime, na.rm = TRUE) - pseudotime + min(pseudotime, na.rm = TRUE)
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\n")
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "redPATH pseudotime (scumi hvg_ctet)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pseudotime, sub("[.]png$", ".tsv", OUT_PLOT))
