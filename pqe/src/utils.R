# Shared helpers for pseudotime run scripts.
# Source this file at the top of each run_*.R script.

library(Matrix)

# Read h5ad and return a list with preprocessed data.
# max_cells: downsample to at most this many annotated cells (NULL = no limit)
# n_hvg:     number of highly variable genes to retain
load_adata <- function(h5ad_path, max_cells = NULL, n_hvg = 2000, gene_list = NULL) {
  rds_path <- sub("\\.h5ad$", ".rds", h5ad_path)
  if (file.exists(rds_path)) {
    cat("Loading pre-saved RDS:", rds_path, "\n")
    saved    <- readRDS(rds_path)
    counts   <- saved$counts
    metadata <- saved$metadata
    rownames(counts) <- saved$gene_names
  } else {
    cat("Reading", h5ad_path, "\n")
    library(zellkonverter)
    library(SummarizedExperiment)
    sce      <- readH5AD(h5ad_path)
    counts   <- assay(sce, "X")
    metadata <- as.data.frame(colData(sce))
  }

  # Keep only cells with spatial coordinates (annotated cells)
  sc_x <- suppressWarnings(as.numeric(as.character(metadata$sc_x)))
  keep <- !is.na(sc_x)
  counts   <- counts[, keep]
  metadata <- metadata[keep, ]
  cat("Annotated cells:", sum(keep), "\n")

  # Optional downsampling
  if (!is.null(max_cells) && ncol(counts) > max_cells) {
    set.seed(42)
    idx      <- sample(ncol(counts), max_cells)
    counts   <- counts[, idx]
    metadata <- metadata[idx, ]
    cat("Downsampled to:", max_cells, "\n")
  }

  # X is already log-normalized; use directly
  log_norm <- counts

  # Gene selection: use gene_list if provided, else select top HVGs by variance
  if (!is.null(gene_list)) {
    available <- intersect(gene_list, rownames(log_norm))
    cat("Using", length(available), "of", length(gene_list), "provided genes present in data\n")
    expr <- as.matrix(t(log_norm[available, ]))
  } else {
    row_mean  <- rowMeans(log_norm)
    row_mean2 <- rowMeans(log_norm^2)
    gene_var  <- row_mean2 - row_mean^2
    n_hvg     <- min(n_hvg, sum(gene_var > 0))
    hvg       <- order(gene_var, decreasing = TRUE)[seq_len(n_hvg)]
    expr      <- as.matrix(t(log_norm[hvg, ]))   # cells × genes
  }

  sc_x_vals <- suppressWarnings(as.numeric(as.character(metadata$sc_x)))
  sc_y_vals <- suppressWarnings(as.numeric(as.character(metadata$sc_y)))

  list(
    expr       = expr,
    sc_x       = sc_x_vals,
    sc_y       = sc_y_vals,
    cell_names = colnames(counts),
    coldata    = metadata
  )
}

# Save cell metadata + pseudotime column as TSV.
# extra_cols: optional data.frame of additional columns to append (e.g. DPT2, DPT3)
save_pseudotime_tsv <- function(cell_names, coldata, pseudotime, out_path, extra_cols = NULL) {
  df <- coldata
  df$cell_id    <- cell_names
  df$pseudotime <- as.numeric(pseudotime)
  if (!is.null(extra_cols)) df <- cbind(df, extra_cols)
  write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("TSV saved to", out_path, "\n")
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
