# Shared helpers for pseudotime run scripts.
# Source this file at the top of each run_*.R script.

library(zellkonverter)
library(Matrix)

# Read h5ad and return a list with preprocessed data.
# max_cells: downsample to at most this many annotated cells (NULL = no limit)
# n_hvg:     number of highly variable genes to retain
load_adata <- function(h5ad_path, max_cells = NULL, n_hvg = 2000) {
  cat("Reading", h5ad_path, "\n")
  sce <- readH5AD(h5ad_path)

  # Keep only cells with spatial coordinates (annotated cells)
  sc_x <- suppressWarnings(as.numeric(as.character(colData(sce)$sc_x)))
  keep <- !is.na(sc_x)
  sce  <- sce[, keep]
  cat("Annotated cells:", ncol(sce), "\n")

  # Optional downsampling
  if (!is.null(max_cells) && ncol(sce) > max_cells) {
    set.seed(42)
    sce <- sce[, sample(ncol(sce), max_cells)]
    cat("Downsampled to:", max_cells, "\n")
  }

  # Log-normalize: counts per 10k, then log1p
  counts    <- assay(sce, "X")               # genes × cells
  lib_size  <- colSums(counts)
  norm      <- t(t(counts) * (1e4 / lib_size))
  log_norm  <- log1p(norm)

  # Select top HVGs by variance (computed on sparse matrix)
  row_mean  <- rowMeans(log_norm)
  row_mean2 <- rowMeans(log_norm^2)
  gene_var  <- row_mean2 - row_mean^2
  n_hvg     <- min(n_hvg, sum(gene_var > 0))
  hvg       <- order(gene_var, decreasing = TRUE)[seq_len(n_hvg)]
  expr      <- as.matrix(t(log_norm[hvg, ]))   # cells × genes

  sc_x_vals <- suppressWarnings(as.numeric(as.character(colData(sce)$sc_x)))
  sc_y_vals <- suppressWarnings(as.numeric(as.character(colData(sce)$sc_y)))

  list(
    expr       = expr,
    sc_x       = sc_x_vals,
    sc_y       = sc_y_vals,
    cell_names = colnames(sce),
    coldata    = as.data.frame(colData(sce))
  )
}

# Save a pseudotime scatter plot coloured by pseudotime using sc_x/sc_y.
save_pseudotime_plot <- function(sc_x, sc_y, pseudotime, out_path, title = "Pseudotime") {
  library(ggplot2)

  df <- data.frame(
    sc_x       = sc_x,
    sc_y       = sc_y,
    pseudotime = as.numeric(pseudotime)
  )
  # Remove cells where pseudotime could not be computed (NA)
  df <- df[!is.na(df$pseudotime), ]

  p <- ggplot(df, aes(x = sc_x, y = sc_y, colour = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_viridis_c(option = "plasma") +
    labs(title = title, colour = "Pseudotime") +
    theme_minimal(base_size = 11)

  ggsave(out_path, p, width = 6, height = 5, dpi = 150)
  cat("Plot saved to", out_path, "\n")
}
