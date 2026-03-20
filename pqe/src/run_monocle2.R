#!/usr/bin/env Rscript
# Run Monocle2 pseudotime on zmanseq data.
# Change ADATA_PATH to point to a different h5ad file.

ADATA_PATH <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_PLOT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime/monocle2_pseudotime.png"

source(file.path(dirname(sys.frame(1)$ofile), "utils.R"))
library(monocle)

# Load and preprocess
dat <- load_adata(ADATA_PATH, max_cells = NULL, n_hvg = 2000)

# Build CellDataSet — use normalized log counts with tobit expressionFamily
expr_mat <- t(dat$expr)  # genes × cells, as monocle expects
pd <- new("AnnotatedDataFrame",
          data = data.frame(row.names = dat$cell_names,
                            sc_x = dat$sc_x, sc_y = dat$sc_y))
fd <- new("AnnotatedDataFrame",
          data = data.frame(row.names = rownames(expr_mat),
                            gene_short_name = rownames(expr_mat)))

cat("Building CellDataSet...\n")
cds <- newCellDataSet(expr_mat, phenoData = pd, featureData = fd,
                      expressionFamily = tobit())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Use high-dispersion genes for ordering
disp_table     <- dispersionTable(cds)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.1 &
                         dispersion_empirical >= 1 * dispersion_fit)$gene_id
cat("Ordering genes:", length(ordering_genes), "\n")
cds <- setOrderingFilter(cds, ordering_genes)

# Reduce dimensionality with DDRTree and order cells
cat("Running Monocle2 DDRTree...\n")
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

pseudotime <- pData(cds)$Pseudotime
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "Monocle2 pseudotime")
