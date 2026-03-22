#!/usr/bin/env Rscript
# Convert mData_CosMx1k_CosMx6k.RData + spatial_transcriptomics.zip → h5ad
# Also plots cell types with spatial coordinates.
# Requires: anndata, Matrix, ggplot2, dplyr

library(anndata)
library(Matrix)
library(ggplot2)
library(dplyr)

DATA_DIR <- "/mnt/thechenlab/ClaudiaC/gbm_migliozzi"
OUT_DIR  <- "/home/unix/cchu/projects/ZmanR/pqe/results/10"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

THEME <- theme_void(base_size = 10) +
  theme(strip.text      = element_text(face = "bold", size = 8),
        legend.key.size  = unit(0.4, "cm"),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.title     = element_text(face = "bold"))


# ══════════════════════════════════════════════════════════════════════════════
# 1 — Load metadata
# ══════════════════════════════════════════════════════════════════════════════

cat("Loading metadata ...\n")
load(file.path(DATA_DIR, "mData_CosMx1k_CosMx6k.RData"))
md <- mData_CosMx1k_CosMx6k
rm(mData_CosMx1k_CosMx6k)

num_cols <- c("simplicityGPM", "simplicityMTC", "simplicityNEU", "simplicityPPR", "TSPS")
md[num_cols] <- lapply(md[num_cols], as.numeric)

# Unique obs key: Sample_ID + Cell_ID
md$cell_key  <- paste(md$Sample_ID, md$Cell_ID, sep = "__")
rownames(md) <- md$cell_key

cat(sprintf("  %s cells x %d metadata columns\n", format(nrow(md), big.mark = ","), ncol(md)))

tsv_out <- file.path(OUT_DIR, "mData_CosMx1k_CosMx6k.tsv.gz")
write.table(md, gzfile(tsv_out), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved %s\n", basename(tsv_out)))


# ══════════════════════════════════════════════════════════════════════════════
# 2 — Iterate over per-sample zips: read expression + coordinates
# ══════════════════════════════════════════════════════════════════════════════

zip_path     <- file.path(DATA_DIR, "spatial_transcriptomics.zip")
zip_list_raw <- system(paste("unzip -l", shQuote(zip_path)), intern = TRUE)
zip_entries  <- zip_list_raw[grepl("\\.zip$", zip_list_raw)]

parsed <- lapply(zip_entries, function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  data.frame(Name = parts[4], Size_B = as.numeric(parts[1]), stringsAsFactors = FALSE)
})
zip_df <- do.call(rbind, parsed)
zip_df$sample <- sub("data_deposit_cosmx_paper/", "", sub("\\.zip$", "", zip_df$Name))

# Human CosMx 1k samples (UM01-UM07; 986-gene panel, metadata available)
HUMAN_1K <- c("UM01", "UM02", "UM03", "UM04", "UM05", "UM06-P", "UM06-R", "UM07-R")
zip_df <- zip_df[zip_df$sample %in% HUMAN_1K, ]
cat(sprintf("Human CosMx 1k samples to process: %s\n", paste(zip_df$sample, collapse = ", ")))

tmp_dir <- file.path(tempdir(), "cosmx_h5ad")
dir.create(tmp_dir, showWarnings = FALSE)

mat_list    <- list()
coord_list  <- list()
gene_names  <- NULL

for (i in seq_len(nrow(zip_df))) {
  nm     <- zip_df$Name[i]
  sname  <- zip_df$sample[i]
  size_b <- zip_df$Size_B[i]

  cat(sprintf("[%d/%d] %s (%.0f MB) ...\n", i, nrow(zip_df), sname, size_b / 1e6))

  # Extract inner zip from outer zip
  system(paste("unzip -o", shQuote(zip_path), shQuote(nm), "-d", shQuote(tmp_dir)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  inner_zip <- file.path(tmp_dir, nm)
  inner_dir <- file.path(tmp_dir, sname)
  dir.create(inner_dir, showWarnings = FALSE)
  system(paste("unzip -o", shQuote(inner_zip), "-d", shQuote(inner_dir)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)

  mat_f    <- list.files(inner_dir, "matrix\\.mtx$",      recursive = TRUE, full.names = TRUE)
  cell_f   <- list.files(inner_dir, "cells\\.tsv$",       recursive = TRUE, full.names = TRUE)
  gene_f   <- list.files(inner_dir, "genes\\.tsv$",       recursive = TRUE, full.names = TRUE)
  coord_f  <- list.files(inner_dir, "coordinates\\.tsv$", recursive = TRUE, full.names = TRUE)

  if (!length(mat_f) || !length(cell_f) || !length(gene_f)) {
    cat("    Skipping — missing MTX/cells/genes\n"); next
  }

  # MTX is already cells x genes
  mat   <- readMM(mat_f[1])
  cells <- readLines(cell_f[1])
  genes <- readLines(gene_f[1])
  stopifnot(nrow(mat) == length(cells), ncol(mat) == length(genes))

  if (is.null(gene_names)) gene_names <- genes

  rownames(mat) <- paste(sname, cells, sep = "__")
  colnames(mat) <- genes
  mat_list[[sname]] <- mat

  # Coordinates
  if (length(coord_f)) {
    coords <- read.table(coord_f[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    coords$cell_key <- paste(sname, coords$cell_id, sep = "__")
    # Extract FOV from cell_id (format: {sample}_{fov}_{cell})
    parts <- strsplit(coords$cell_id, "_")
    coords$fov <- sapply(parts, function(p) p[length(p) - 1])
    coord_list[[sname]] <- coords[, c("cell_key", "x", "y", "fov")]
  }

  # Free extracted files
  unlink(inner_dir, recursive = TRUE)
  unlink(inner_zip)
}


all_plot_dfs <- list()

# ══════════════════════════════════════════════════════════════════════════════
# 3 — Per-sample: align to metadata, write h5ad, plot
# (samples use different gene panels so are kept separate)
# ══════════════════════════════════════════════════════════════════════════════

for (sname in names(mat_list)) {
  cat(sprintf("\n── Processing %s ──\n", sname))

  mat        <- mat_list[[sname]]
  has_coords <- !is.null(coord_list[[sname]])

  common <- intersect(rownames(mat), md$cell_key)
  cat(sprintf("  Cells matched to metadata: %s / %s  |  coordinates: %s\n",
              format(length(common), big.mark = ","),
              format(nrow(mat), big.mark = ","),
              ifelse(has_coords, "yes", "no")))

  all_keys  <- rownames(mat)
  meta_cols <- setdiff(colnames(md), "cell_key")
  obs_aln   <- as.data.frame(md[match(all_keys, md$cell_key), meta_cols, drop = FALSE])
  rownames(obs_aln) <- all_keys
  # Store the cell ID (without sname prefix) so Python can use it for joining
  obs_aln$cell_id_raw <- sub(paste0("^", sname, "__"), "", all_keys)

  # anndata cannot write NA in character columns — replace with empty string
  chr_cols <- sapply(obs_aln, is.character)
  obs_aln[chr_cols] <- lapply(obs_aln[chr_cols], function(x) ifelse(is.na(x), "", x))

  genes_s <- colnames(mat)
  var_df  <- data.frame(row.names = genes_s)

  # ── h5ad ──────────────────────────────────────────────────────────────────
  cat("  Building AnnData ...\n")
  if (has_coords) {
    coords     <- coord_list[[sname]]
    rownames(coords) <- coords$cell_key
    coords_aln <- as.matrix(coords[all_keys, c("x", "y")])
    rownames(coords_aln) <- all_keys
    adata <- AnnData(X = mat, obs = obs_aln, var = var_df,
                     obsm = list(spatial = coords_aln))
  } else {
    coords_aln <- NULL
    adata <- AnnData(X = mat, obs = obs_aln, var = var_df)
  }

  out_h5ad <- file.path(OUT_DIR, sprintf("cosmx_%s.h5ad", tolower(sname)))
  adata$write_h5ad(out_h5ad)
  cat(sprintf("  Saved %s\n", basename(out_h5ad)))

  # ── Per-patient spatial plot (cell types) ─────────────────────────────────
  if (has_coords) {
    plot_df <- obs_aln %>%
      select(cellType, cellTypeWithMyeloidSubsets) %>%
      mutate(x = coords_aln[, "x"], y = coords_aln[, "y"])

    p_ct <- plot_df %>%
      filter(!is.na(cellType)) %>%
      ggplot(aes(x, y, colour = cellType)) +
      geom_point(size = 0.15, alpha = 0.6, stroke = 0) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = sprintf("CosMx %s — cell types", sname), colour = "Cell type") +
      THEME

    ggsave(file.path(OUT_DIR, sprintf("10_spatial_celltypes_%s.pdf", tolower(sname))),
           p_ct, width = 8, height = 7)
    cat(sprintf("  Saved 10_spatial_celltypes_%s.pdf\n", tolower(sname)))

    p_myl <- plot_df %>%
      filter(!is.na(cellTypeWithMyeloidSubsets)) %>%
      ggplot(aes(x, y, colour = cellTypeWithMyeloidSubsets)) +
      geom_point(size = 0.15, alpha = 0.6, stroke = 0) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = sprintf("CosMx %s — myeloid subsets", sname), colour = "Cell type") +
      THEME

    ggsave(file.path(OUT_DIR, sprintf("11_spatial_myeloid_%s.pdf", tolower(sname))),
           p_myl, width = 8, height = 7)
    cat(sprintf("  Saved 11_spatial_myeloid_%s.pdf\n", tolower(sname)))

    # Accumulate for combined multi-patient plot
    all_plot_dfs[[sname]] <- plot_df %>% mutate(patient = sname)
  } else {
    cat("  Skipping spatial plots (no coordinates)\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 4b — Combined multi-patient spatial plots
# ══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating combined multi-patient spatial plots ...\n")
all_df <- do.call(rbind, all_plot_dfs)

p_all_ct <- all_df %>%
  filter(!is.na(cellType)) %>%
  ggplot(aes(x, y, colour = cellType)) +
  geom_point(size = 0.1, alpha = 0.5, stroke = 0) +
  facet_wrap(~ patient, scales = "free", ncol = 4) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "CosMx 1k — cell types by patient", colour = "Cell type") +
  THEME

ggsave(file.path(OUT_DIR, "12_spatial_celltypes_all_patients.pdf"),
       p_all_ct, width = 22, height = 12)
cat("Saved 12_spatial_celltypes_all_patients.pdf\n")

p_all_myl <- all_df %>%
  filter(!is.na(cellTypeWithMyeloidSubsets)) %>%
  ggplot(aes(x, y, colour = cellTypeWithMyeloidSubsets)) +
  geom_point(size = 0.1, alpha = 0.5, stroke = 0) +
  facet_wrap(~ patient, scales = "free", ncol = 4) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "CosMx 1k — myeloid subsets by patient", colour = "Cell type") +
  THEME

ggsave(file.path(OUT_DIR, "13_spatial_myeloid_all_patients.pdf"),
       p_all_myl, width = 22, height = 12)
cat("Saved 13_spatial_myeloid_all_patients.pdf\n")

cat("\nDone. Output in", OUT_DIR, "\n")
