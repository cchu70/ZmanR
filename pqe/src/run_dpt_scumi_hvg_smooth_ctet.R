#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot_clean.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_scumi_hvg_smooth_ctet/dpt_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(destiny)

genes <- load_ctet_mos_genes("scumi", geneset="hvg", timecol="smooth_ctet")
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Building DiffusionMap...\n"); set.seed(42)
dm  <- DiffusionMap(dat$expr, n_pcs = min(30, ncol(dat$expr) - 1), n_eigs = 10, verbose = TRUE)
start_idx  <- which.min(dat$coldata$cTET)
atrem_idxs <- which(dat$coldata$enrichment > 0)
igg_idxs   <- which(dat$coldata$enrichment < 0)
tip_atrem  <- atrem_idxs[which.max(dat$coldata$cTET[atrem_idxs])]
tip_igg    <- igg_idxs[which.max(dat$coldata$cTET[igg_idxs])]
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\n")
cat("Tip aTrem2:", dat$cell_names[tip_atrem], "(cTET =", dat$coldata$cTET[tip_atrem], ")\n")
cat("Tip IgG:  ", dat$cell_names[tip_igg],   "(cTET =", dat$coldata$cTET[tip_igg],   ")\n")
cat("Computing DPT...\n")
dpt <- DPT(dm, tips = c(start_idx, tip_atrem, tip_igg))
dpt_mat  <- as.matrix(dpt)
pt_root  <- dpt_mat[, start_idx]
pt_atrem <- dpt_mat[, tip_atrem]
pt_igg   <- dpt_mat[, tip_igg]
dpt_cols <- data.frame(DPT_root = pt_root, DPT_atrem = pt_atrem, DPT_igg = pt_igg)
save_pseudotime_plot(dat$sc_x, dat$sc_y, pt_root, OUT_PLOT, "DPT pseudotime (scumi hvg_smooth_ctet)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pt_root, sub("\.png$", ".tsv", OUT_PLOT), extra_cols = dpt_cols)
