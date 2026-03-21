#!/usr/bin/env Rscript
ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_momac/monocle2_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(monocle)

genes <- load_naive_genes()
dat   <- load_adata(ADATA_PATH, gene_list = genes)
expr_mat <- t(dat$expr)
pd <- new("AnnotatedDataFrame", data = data.frame(row.names = dat$cell_names,
          sc_x = dat$sc_x, sc_y = dat$sc_y))
fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(expr_mat),
          gene_short_name = rownames(expr_mat)))
cat("Building CellDataSet...\n")
cds <- newCellDataSet(expr_mat, phenoData = pd, featureData = fd, expressionFamily = tobit())
cds <- estimateSizeFactors(cds)
cds <- setOrderingFilter(cds, rownames(expr_mat))
cat("Running Monocle2 DDRTree...\n")
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
pseudotime <- pData(cds)$Pseudotime
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "Monocle2 pseudotime (momac)")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pseudotime, sub("\\.png$", ".tsv", OUT_PLOT))
