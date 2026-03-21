#!/usr/bin/env Rscript
# One-time helper: read zmanseq.h5ad and save raw counts + metadata as RDS.
# Run with dpt_env (or any env with working zellkonverter).
# Output is consumed by utils.R when zellkonverter fails (e.g. monocle2_env).

library(zellkonverter)
library(SummarizedExperiment)
library(Matrix)

H5AD  <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT   <- "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.rds"

cat("Reading", H5AD, "\n")
sce <- readH5AD(H5AD)

cat("Extracting counts and metadata...\n")
counts   <- assay(sce, "X")          # genes × cells (sparse)
metadata <- as.data.frame(colData(sce))

cat("Cells:", ncol(sce), "Genes:", nrow(sce), "\n")
saveRDS(list(counts = counts, metadata = metadata, gene_names = rownames(sce)),
        file = OUT)
cat("Saved to", OUT, "\n")
