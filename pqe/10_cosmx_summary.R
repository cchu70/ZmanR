#!/usr/bin/env Rscript
# Summary of CosMx spatial transcriptomics data:
#   1. mData_CosMx1k_CosMx6k.RData  — cell metadata (cellType, simplicity, SWED, TSPS, etc.)
#   2. spatial_transcriptomics.zip  — per-sample MTX + coordinates

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)
library(gridExtra)

DATA_DIR <- "/mnt/thechenlab/ClaudiaC/gbm_migliozzi"
OUT_DIR  <- "/home/unix/cchu/projects/ZmanR/pqe/results/10"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

THEME <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), plot.title = element_text(face = "bold"))


# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — mData metadata table
# ══════════════════════════════════════════════════════════════════════════════

cat("Loading mData_CosMx1k_CosMx6k.RData ...\n")
load(file.path(DATA_DIR, "mData_CosMx1k_CosMx6k.RData"))
md <- mData_CosMx1k_CosMx6k
rm(mData_CosMx1k_CosMx6k)

cat(sprintf("mData: %d cells × %d columns\n", nrow(md), ncol(md)))
cat("Columns:", paste(colnames(md), collapse = ", "), "\n")
cat("Samples:", paste(sort(unique(md$Sample_ID)), collapse = ", "), "\n")

# All columns are character — convert numeric ones
num_cols <- c("simplicityGPM", "simplicityMTC", "simplicityNEU", "simplicityPPR", "SWED", "TSPS")
md[num_cols] <- lapply(md[num_cols], as.numeric)


# ── 1a. Cells per sample ──────────────────────────────────────────────────────

p_cells <- md %>%
  count(Sample_ID) %>%
  arrange(desc(n)) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = Sample_ID)) %>%
  ggplot(aes(x = Sample_ID, y = n, fill = Sample_ID)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.3, size = 3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cells per sample", x = NULL, y = "# cells") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "01_cells_per_sample.pdf"), p_cells, width = 8, height = 4)


# ── 1b. Cell type composition per sample ─────────────────────────────────────

ct_order <- md %>% count(cellType) %>% arrange(desc(n)) %>% pull(cellType)

p_ct <- md %>%
  count(Sample_ID, cellType) %>%
  group_by(Sample_ID) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  mutate(cellType = factor(cellType, levels = ct_order)) %>%
  ggplot(aes(x = Sample_ID, y = pct, fill = cellType)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(title = "Cell type composition per sample", x = NULL, y = "Fraction", fill = "Cell type") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(OUT_DIR, "02_celltype_composition.pdf"), p_ct, width = 10, height = 5)


# ── 1c. Myeloid subsets per sample ───────────────────────────────────────────

myl_order <- md %>% count(cellTypeWithMyeloidSubsets) %>% arrange(desc(n)) %>% pull(cellTypeWithMyeloidSubsets)

p_myl <- md %>%
  count(Sample_ID, cellTypeWithMyeloidSubsets) %>%
  group_by(Sample_ID) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  mutate(cellTypeWithMyeloidSubsets = factor(cellTypeWithMyeloidSubsets, levels = myl_order)) %>%
  ggplot(aes(x = Sample_ID, y = pct, fill = cellTypeWithMyeloidSubsets)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(title = "Cell type (with myeloid subsets) per sample",
       x = NULL, y = "Fraction", fill = "Cell type") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(OUT_DIR, "03_myeloid_subsets_composition.pdf"), p_myl, width = 12, height = 5)


# ── 1d. Simplicity scores distributions ──────────────────────────────────────

simp_cols <- c("simplicityGPM", "simplicityMTC", "simplicityNEU", "simplicityPPR")

p_simp <- md %>%
  select(Sample_ID, all_of(simp_cols)) %>%
  pivot_longer(all_of(simp_cols), names_to = "score", values_to = "value") %>%
  mutate(score = sub("simplicity", "", score)) %>%
  ggplot(aes(x = value, fill = score)) +
  geom_histogram(bins = 60, show.legend = FALSE) +
  facet_grid(score ~ Sample_ID, scales = "free_y") +
  labs(title = "Simplicity score distributions", x = "Score", y = "# cells") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text.x = element_text(size = 7))

ggsave(file.path(OUT_DIR, "04_simplicity_scores.pdf"), p_simp, width = 16, height = 6)


# ── 1e. SWED and TSPS distributions ──────────────────────────────────────────

p_swed <- md %>%
  ggplot(aes(x = Sample_ID, y = SWED, fill = Sample_ID)) +
  geom_violin(show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", show.legend = FALSE) +
  labs(title = "SWED per sample", x = NULL, y = "SWED") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_tsps <- md %>%
  ggplot(aes(x = Sample_ID, y = TSPS, fill = Sample_ID)) +
  geom_violin(show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", show.legend = FALSE) +
  labs(title = "TSPS per sample", x = NULL, y = "TSPS") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "05_SWED_TSPS.pdf"),
       arrangeGrob(p_swed, p_tsps, ncol = 2), width = 14, height = 5)


# ── 1f. Spatial homotypic pattern frequency ───────────────────────────────────

p_pattern <- md %>%
  count(Sample_ID, spatial_homotypic_pattern) %>%
  group_by(Sample_ID) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = Sample_ID, y = pct, fill = spatial_homotypic_pattern)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(title = "Spatial homotypic pattern per sample",
       x = NULL, y = "Fraction", fill = "Pattern") +
  THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "06_spatial_homotypic_pattern.pdf"), p_pattern, width = 10, height = 5)


# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — spatial_transcriptomics.zip: sample inventory and cell counts
# ══════════════════════════════════════════════════════════════════════════════

cat("\nSummarising spatial_transcriptomics.zip ...\n")

zip_path <- file.path(DATA_DIR, "spatial_transcriptomics.zip")
zip_list <- unzip(zip_path, list = TRUE)
cat(sprintf("Top-level entries: %d files\n", nrow(zip_list)))

# Each entry is a nested zip: extract small ones for detailed peek
nested_zips <- zip_list$Name
tmp_dir <- tempdir()

summary_rows <- list()

for (nz in nested_zips) {
  sample_name <- sub("data_deposit_cosmx_paper/", "", nz)
  sample_name <- sub("\\.zip$", "", sample_name)

  file_size_mb <- zip_list$Length[zip_list$Name == nz] / 1e6
  cat(sprintf("  %s (%.0f MB)\n", sample_name, file_size_mb))

  # Only extract small samples (<100 MB) for cell/gene counts; others describe by size
  n_cells <- NA_integer_
  n_genes <- NA_integer_

  if (file_size_mb < 100) {
    tryCatch({
      unzip(zip_path, files = nz, exdir = tmp_dir, overwrite = TRUE)
      inner_path <- file.path(tmp_dir, nz)
      inner_files <- unzip(inner_path, list = TRUE)$Name

      cells_file <- inner_files[grepl("cells\\.tsv$", inner_files)]
      genes_file <- inner_files[grepl("genes\\.tsv$", inner_files)]
      mtx_file   <- inner_files[grepl("matrix\\.mtx$", inner_files)]

      if (length(cells_file)) {
        cells_tmp <- tempfile()
        unzip(inner_path, files = cells_file, exdir = tempdir(), overwrite = TRUE)
        n_cells <- length(readLines(file.path(tempdir(), cells_file)))
      }
      if (length(genes_file)) {
        genes_tmp <- tempfile()
        unzip(inner_path, files = genes_file, exdir = tempdir(), overwrite = TRUE)
        n_genes <- length(readLines(file.path(tempdir(), genes_file)))
      }
    }, error = function(e) cat("    Warning:", conditionMessage(e), "\n"))
  }

  summary_rows[[sample_name]] <- data.frame(
    Sample    = sample_name,
    Size_MB   = round(file_size_mb, 1),
    N_cells   = n_cells,
    N_genes   = n_genes
  )
}

zip_summary <- bind_rows(summary_rows)
write.csv(zip_summary, file.path(OUT_DIR, "zip_sample_summary.csv"), row.names = FALSE)
print(zip_summary)

# Plot: sample file sizes and cell counts where available
p_zip_size <- zip_summary %>%
  mutate(Sample = factor(Sample, levels = Sample[order(Size_MB)])) %>%
  ggplot(aes(x = Sample, y = Size_MB, fill = !is.na(N_cells))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#4C9BE8", "FALSE" = "#AAAAAA"),
                    labels = c("TRUE" = "Extracted", "FALSE" = "Size only"),
                    name = NULL) +
  labs(title = "spatial_transcriptomics.zip: file sizes",
       x = NULL, y = "Compressed size (MB)") +
  THEME

p_zip_cells <- zip_summary %>%
  filter(!is.na(N_cells)) %>%
  mutate(Sample = factor(Sample, levels = Sample[order(N_cells)])) %>%
  ggplot(aes(x = Sample, y = N_cells)) +
  geom_col(fill = "#4C9BE8") +
  geom_text(aes(label = comma(N_cells)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Cell counts (small samples extracted)", x = NULL, y = "# cells") +
  THEME

ggsave(file.path(OUT_DIR, "07_zip_sample_overview.pdf"),
       arrangeGrob(p_zip_size, p_zip_cells, ncol = 2), width = 14, height = 5)


cat("\nAll figures saved to", OUT_DIR, "\n")
