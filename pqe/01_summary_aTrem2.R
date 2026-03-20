library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

outdir    <- "results/01"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

data_dir  <- "/mnt/thechenlab/ClaudiaC/zmanseq"
meta_file <- file.path(data_dir, "metadata.txt")

# ==============================================================================
# 1. LOAD METADATA
# Metadata uses space delimiters but mouse names can contain spaces (e.g. "Mouse 3")
# so we parse from both ends using fixed column positions.
# Column positions (1-indexed, field 1 = row index):
#   4  = Amp_batch_ID
#   10 = Cell_barcode
#   NF = Treatment  |  NF-1 = sc_y  |  NF-2 = sc_x
#   NF-3 = cluster_colors  |  NF-4 = time_assignment
# ==============================================================================

cat("Loading metadata...\n")
lines <- readLines(meta_file)

parse_meta <- function(lines) {
  pb <- txtProgressBar(min = 0, max = length(lines) - 1, style = 3)
  rows <- lapply(seq_along(lines[-1]), function(i) {
    setTxtProgressBar(pb, i)
    f <- strsplit(lines[i + 1], " ")[[1]]
    n <- length(f)
    list(
      Well_ID        = f[2],
      Amp_batch_ID   = f[4],
      Cell_barcode   = f[10],
      time_assignment = f[n - 4],
      celltype       = f[n - 3],
      sc_x           = suppressWarnings(as.numeric(f[n - 2])),
      sc_y           = suppressWarnings(as.numeric(f[n - 1])),
      Treatment      = f[n]
    )
  })
  close(pb)
  as.data.frame(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)),
                stringsAsFactors = FALSE)
}

meta <- parse_meta(lines)
cat("\nMetadata rows:", nrow(meta), "\n")
cat("Treatment counts:\n")
print(table(meta$Treatment))

# Filter to aTrem2 and IgG only
meta_sub <- meta[meta$Treatment %in% c("aTrem2", "IgG"), ]
cat("aTrem2 + IgG cells:", nrow(meta_sub), "\n")

# ==============================================================================
# 2. COMPUTE PER-CELL QC FROM COUNT MATRICES
# Load each plate, compute total UMI and genes detected per cell, then discard.
# ==============================================================================

cat("\nComputing per-cell QC from count matrices...\n")

plates_atrem2 <- list.files(file.path(data_dir, "aTrem2"), full.names = TRUE)
plates_igg    <- list.files(file.path(data_dir, "IgG"),    full.names = TRUE)
all_plates    <- c(plates_atrem2, plates_igg)

qc_list <- lapply(all_plates, function(f) {
  plate_id <- sub(".*_(AB\\d+)\\.txt\\.gz", "\\1", basename(f))
  treatment <- if (grepl("aTrem2", f)) "aTrem2" else "IgG"
  cat(sprintf("  %s (%s)\n", plate_id, treatment))

  mat <- read.table(gzfile(f), header = TRUE, row.names = 1,
                    check.names = FALSE, sep = "\t")
  data.frame(
    Well_ID      = colnames(mat),
    Amp_batch_ID = plate_id,
    Treatment    = treatment,
    total_UMI    = colSums(mat),
    genes_det    = colSums(mat > 0),
    stringsAsFactors = FALSE
  )
})
qc <- do.call(rbind, qc_list)
cat("QC rows:", nrow(qc), "\n")

# Join QC with metadata
merged <- merge(qc, meta_sub[, c("Well_ID", "time_assignment", "celltype", "sc_x", "sc_y")],
                by = "Well_ID", all.x = TRUE)

# Clean up
merged$Treatment    <- factor(merged$Treatment, levels = c("aTrem2", "IgG"))
merged$time_assignment <- factor(merged$time_assignment,
                                  levels = c("12H", "24H", "36H", "Negative", "Not_Assigned"))
merged$is_doublet   <- merged$celltype == "Doublet" & !is.na(merged$celltype)

# ==============================================================================
# 3. SUMMARY TABLES
# ==============================================================================

cat("\nWriting summary tables...\n")

## 3a. Per-plate summary
plate_summary <- merged %>%
  group_by(Treatment, Amp_batch_ID) %>%
  summarise(
    n_cells      = n(),
    median_UMI   = median(total_UMI),
    mean_UMI     = round(mean(total_UMI), 1),
    median_genes = median(genes_det),
    mean_genes   = round(mean(genes_det), 1),
    pct_zero_umi = round(mean(total_UMI == 0) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(Treatment, Amp_batch_ID)

write.csv(plate_summary, file.path(outdir, "plate_summary.csv"), row.names = FALSE)

## 3b. Per-treatment summary
treatment_summary <- merged %>%
  group_by(Treatment) %>%
  summarise(
    n_cells      = n(),
    n_plates     = n_distinct(Amp_batch_ID),
    median_UMI   = median(total_UMI),
    mean_UMI     = round(mean(total_UMI), 1),
    median_genes = median(genes_det),
    mean_genes   = round(mean(genes_det), 1),
    pct_zero_umi = round(mean(total_UMI == 0) * 100, 1),
    .groups = "drop"
  )

write.csv(treatment_summary, file.path(outdir, "treatment_summary.csv"), row.names = FALSE)

## 3c. Cell type composition by treatment
celltype_summary <- merged %>%
  filter(!is.na(celltype)) %>%
  group_by(Treatment, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(pct = round(n / sum(n) * 100, 2)) %>%
  arrange(Treatment, desc(n))

write.csv(celltype_summary, file.path(outdir, "celltype_by_treatment.csv"), row.names = FALSE)

## 3d. Time assignment by treatment
time_summary <- merged %>%
  filter(!is.na(time_assignment)) %>%
  group_by(Treatment, time_assignment) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

write.csv(time_summary, file.path(outdir, "time_assignment_by_treatment.csv"), row.names = FALSE)

## 3e. Doublet summary per treatment and plate
doublet_summary <- merged %>%
  group_by(Treatment, Amp_batch_ID) %>%
  summarise(
    n_cells    = n(),
    n_doublets = sum(is_doublet),
    pct_doublets = round(mean(is_doublet) * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(Treatment, Amp_batch_ID)

write.csv(doublet_summary, file.path(outdir, "doublet_summary.csv"), row.names = FALSE)

cat("Doublets by treatment:\n")
print(merged %>% group_by(Treatment) %>%
        summarise(n_doublets = sum(is_doublet),
                  pct = round(mean(is_doublet) * 100, 2), .groups = "drop"))

cat("Tables written.\n")

# ==============================================================================
# 4. FIGURES
# ==============================================================================

cat("Generating figures...\n")

pal_treat <- c("aTrem2" = "#E64B35", "IgG" = "#4DBBD5")

# Cell type color palette based on Fig 2A
pal_celltype <- c(
  # T cells (pinks/peach)
  "Treg"              = "#F4C0C8",
  "CD4"               = "#F5B48A",
  "CD8"               = "#E07878",
  "CD8_Dysfunctional" = "#C85050",
  # NK cells (oranges)
  "NK_Chemotactic"    = "#F0C060",
  "Chemotactic"       = "#F0C060",
  "NK_Dysfunctional"  = "#E89030",
  "Dysfunctional"     = "#E89030",
  "NK_intermediate"   = "#E07820",
  "Cytotoxic"         = "#C84010",
  # DC / Monocyte-derived (yellows to orange-red)
  "MigDC"             = "#D45030",
  "cDC1"              = "#E8B030",
  "cDC2"              = "#F0E080",
  "MonDC"             = "#F8D060",
  "pDC"               = "#D8A8D0",
  # Monocytes / macrophages (purples)
  "inflammatory"      = "#7050B0",
  "Monocytes"         = "#8070C0",
  "Monocyte"          = "#9080C8",
  "MoMac"             = "#A088C8",
  "MoMac1"            = "#9878B8",
  "MoMac2"            = "#B098D0",
  "transitory"        = "#C0B0D8",
  # TAMs (teals/cyans)
  "TAM"               = "#70CED8",
  "Acp5_TAM"          = "#50B8C8",
  "Arg1_TAM"          = "#80D8C8",
  "Gpnmb_TAM"         = "#90C8B8",
  # Unannotated / QC
  "Not_annotated"     = "#C8C8C8",
  "Doublet"           = "#909090"
)

## -- QC by treatment ----------------------------------------------------------

p_umi_treat <- ggplot(merged, aes(x = Treatment, y = log10(total_UMI + 1), fill = Treatment)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = pal_treat) +
  labs(title = "Total UMI per cell by treatment",
       x = NULL, y = "log10(total UMI + 1)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_genes_treat <- ggplot(merged, aes(x = Treatment, y = genes_det, fill = Treatment)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = pal_treat) +
  labs(title = "Genes detected per cell by treatment",
       x = NULL, y = "Genes detected") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

pdf(file.path(outdir, "qc_by_treatment.pdf"), width = 6, height = 5)
gridExtra::grid.arrange(p_umi_treat, p_genes_treat, ncol = 2)
dev.off()

## -- QC by plate --------------------------------------------------------------

p_umi_plate <- ggplot(merged, aes(x = reorder(Amp_batch_ID, total_UMI, median),
                                   y = log10(total_UMI + 1), fill = Treatment)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = pal_treat) +
  labs(title = "Total UMI per cell by plate",
       x = "Plate (Amp batch)", y = "log10(total UMI + 1)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

p_genes_plate <- ggplot(merged, aes(x = reorder(Amp_batch_ID, genes_det, median),
                                     y = genes_det, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = pal_treat) +
  labs(title = "Genes detected per cell by plate",
       x = "Plate (Amp batch)", y = "Genes detected") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

pdf(file.path(outdir, "qc_by_plate.pdf"), width = 14, height = 5)
gridExtra::grid.arrange(p_umi_plate, p_genes_plate, nrow = 2)
dev.off()

## -- Cell type composition by treatment ---------------------------------------

p_celltype <- ggplot(celltype_summary, aes(x = Treatment, y = pct, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_celltype, na.value = "#C8C8C8") +
  labs(title = "Cell type composition by treatment",
       x = NULL, y = "% of cells", fill = "Cell type") +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 8))

pdf(file.path(outdir, "celltype_composition.pdf"), width = 7, height = 5)
print(p_celltype)
dev.off()

## -- Cell type composition by plate ------------------------------------------

celltype_plate <- merged %>%
  filter(!is.na(celltype)) %>%
  group_by(Treatment, Amp_batch_ID, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment, Amp_batch_ID) %>%
  mutate(pct = n / sum(n) * 100)

p_celltype_plate <- ggplot(celltype_plate,
                            aes(x = Amp_batch_ID, y = pct, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_celltype, na.value = "#C8C8C8") +
  facet_wrap(~ Treatment, scales = "free_x") +
  labs(title = "Cell type composition by plate",
       x = "Plate", y = "% of cells", fill = "Cell type") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        legend.text = element_text(size = 8))

pdf(file.path(outdir, "celltype_by_plate.pdf"), width = 14, height = 5)
print(p_celltype_plate)
dev.off()

## -- Time assignment by treatment ---------------------------------------------

p_time <- ggplot(time_summary, aes(x = Treatment, y = pct, fill = time_assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2", na.value = "grey80") +
  labs(title = "Time assignment by treatment",
       x = NULL, y = "% of cells", fill = "Time bin") +
  theme_bw(base_size = 12)

pdf(file.path(outdir, "time_assignment.pdf"), width = 6, height = 5)
print(p_time)
dev.off()

## -- UMI vs genes scatter by treatment ----------------------------------------

set.seed(42)
merged_samp <- merged[sample(nrow(merged), min(5000, nrow(merged))), ]

p_scatter <- ggplot(merged_samp, aes(x = log10(total_UMI + 1), y = genes_det,
                                      color = Treatment)) +
  geom_point(alpha = 0.3, size = 0.6) +
  scale_color_manual(values = pal_treat) +
  facet_wrap(~ Treatment) +
  labs(title = "UMI vs genes detected (5k cell sample)",
       x = "log10(total UMI + 1)", y = "Genes detected") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

pdf(file.path(outdir, "umi_vs_genes.pdf"), width = 8, height = 4)
print(p_scatter)
dev.off()

## -- Doublet QC: UMI and genes vs non-doublets --------------------------------

p_doublet_umi <- ggplot(merged, aes(x = is_doublet, y = log10(total_UMI + 1), fill = is_doublet)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#FF7F00"),
                    labels = c("Singlet", "Doublet")) +
  scale_x_discrete(labels = c("FALSE" = "Singlet", "TRUE" = "Doublet")) +
  facet_wrap(~ Treatment) +
  labs(title = "UMI: doublets vs singlets", x = NULL, y = "log10(total UMI + 1)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_doublet_genes <- ggplot(merged, aes(x = is_doublet, y = genes_det, fill = is_doublet)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#FF7F00")) +
  scale_x_discrete(labels = c("FALSE" = "Singlet", "TRUE" = "Doublet")) +
  facet_wrap(~ Treatment) +
  labs(title = "Genes detected: doublets vs singlets", x = NULL, y = "Genes detected") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

pdf(file.path(outdir, "doublet_qc.pdf"), width = 8, height = 4)
gridExtra::grid.arrange(p_doublet_umi, p_doublet_genes, nrow = 2)
dev.off()

# ==============================================================================
# 5. DOUBLET MARKER EXPRESSION
# Check whether doublets co-express T cell and tumor/glioma markers
# ==============================================================================

cat("\nExtracting marker expression for doublets...\n")

# T cell markers subdivided by functional state
t_markers_general    <- c("Cd3e", "Cd3d", "Cd3g", "Trac", "Trbc1", "Trbc2")
t_markers_lineage    <- c("Cd4", "Cd8a", "Cd8b1")
t_markers_activation <- c("Cd44", "Cd69", "Il2ra", "Ifng", "Tnf")
t_markers_effector   <- c("Gzmb", "Gzma", "Prf1", "Fasl")
t_markers_exhaustion <- c("Pdcd1", "Havcr2", "Lag3", "Tigit", "Tox", "Ctla4", "Entpd1")

# GL261 GBM tumor cell markers subdivided by category
tumor_markers_neural  <- c("Nes", "Sox2", "Prom1")        # neural stem / progenitor
tumor_markers_glioma  <- c("Gfap", "Olig2", "S100b")      # glioma lineage
tumor_markers_prolif  <- c("Mki67", "Top2a")               # proliferation
tumor_markers_invasive <- c("Vim", "Fn1", "Cd44")          # mesenchymal / invasive
tumor_markers_receptor <- c("Egfr", "Pdgfra", "Met")       # receptor / signaling

t_markers     <- c(t_markers_general, t_markers_lineage, t_markers_activation,
                   t_markers_effector, t_markers_exhaustion)
tumor_markers <- c(tumor_markers_neural, tumor_markers_glioma, tumor_markers_prolif,
                   tumor_markers_invasive, tumor_markers_receptor)
all_markers   <- c(t_markers, tumor_markers)

# Annotation lookup: gene -> subtype label
marker_subtype <- c(
  setNames(rep("T cell: general",    length(t_markers_general)),    t_markers_general),
  setNames(rep("T cell: lineage",    length(t_markers_lineage)),    t_markers_lineage),
  setNames(rep("T cell: activation", length(t_markers_activation)), t_markers_activation),
  setNames(rep("T cell: effector",   length(t_markers_effector)),   t_markers_effector),
  setNames(rep("T cell: exhaustion", length(t_markers_exhaustion)), t_markers_exhaustion),
  setNames(rep("Tumor: neural stem", length(tumor_markers_neural)),   tumor_markers_neural),
  setNames(rep("Tumor: glioma",      length(tumor_markers_glioma)),   tumor_markers_glioma),
  setNames(rep("Tumor: proliferation", length(tumor_markers_prolif)), tumor_markers_prolif),
  setNames(rep("Tumor: invasive",    length(tumor_markers_invasive)), tumor_markers_invasive),
  setNames(rep("Tumor: receptor",    length(tumor_markers_receptor)), tumor_markers_receptor)
)

# For each plate, extract marker expression for all cells, tag doublets
marker_list <- lapply(all_plates, function(f) {
  plate_id  <- sub(".*_(AB\\d+)\\.txt\\.gz", "\\1", basename(f))
  treatment <- if (grepl("aTrem2", f)) "aTrem2" else "IgG"
  cat(sprintf("  %s\n", plate_id))

  mat <- read.table(gzfile(f), header = TRUE, row.names = 1,
                    check.names = FALSE, sep = "\t")

  # Keep only marker genes present in this matrix
  present <- intersect(all_markers, rownames(mat))
  sub_mat  <- as.data.frame(t(mat[present, , drop = FALSE]))
  sub_mat$Well_ID <- rownames(sub_mat)
  sub_mat
})
marker_df <- do.call(rbind, marker_list)

# Ensure all marker columns exist (fill missing with 0)
for (g in all_markers) {
  if (!g %in% colnames(marker_df)) marker_df[[g]] <- 0
}

# Join with doublet annotation
marker_merged <- merge(marker_df,
                       merged[, c("Well_ID", "Treatment", "Amp_batch_ID",
                                  "celltype", "is_doublet")],
                       by = "Well_ID")
marker_merged$cell_class <- ifelse(marker_merged$is_doublet, "Doublet", "Singlet")

# Compute mean expression per group (doublet vs singlet) per treatment
marker_long <- marker_merged %>%
  select(Well_ID, Treatment, cell_class, all_of(all_markers)) %>%
  pivot_longer(cols = all_of(all_markers), names_to = "gene", values_to = "expr") %>%
  mutate(
    marker_type    = ifelse(gene %in% t_markers, "T cell", "Tumor/Glia (GL261)"),
    marker_subtype = marker_subtype[gene]
  )

# Summary table: mean expression per gene, cell class, and treatment
marker_summary <- marker_long %>%
  group_by(Treatment, cell_class, marker_type, marker_subtype, gene) %>%
  summarise(
    mean_expr   = round(mean(expr), 4),
    pct_expr    = round(mean(expr > 0) * 100, 2),
    .groups = "drop"
  )

write.csv(marker_summary, file.path(outdir, "doublet_marker_summary.csv"), row.names = FALSE)

## -- Plot: % cells expressing each marker, doublets vs singlets ---------------

# Ordered factor for subtype faceting
subtype_levels <- c(
  "T cell: general", "T cell: lineage", "T cell: activation",
  "T cell: effector", "T cell: exhaustion",
  "Tumor: neural stem", "Tumor: glioma", "Tumor: proliferation",
  "Tumor: invasive", "Tumor: receptor"
)
marker_summary$marker_subtype <- factor(marker_summary$marker_subtype, levels = subtype_levels)

p_markers <- ggplot(marker_summary,
                    aes(x = gene, y = pct_expr, fill = cell_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Singlet" = "grey60", "Doublet" = "#FF7F00")) +
  facet_grid(marker_subtype ~ Treatment, scales = "free_x", space = "free_x") +
  labs(title = "% cells expressing T cell and GL261 tumor markers: doublets vs singlets",
       x = "Gene", y = "% cells with > 0 counts", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9))

pdf(file.path(outdir, "doublet_markers.pdf"), width = 12, height = 14)
print(p_markers)
dev.off()

## -- Plot: mean expression heatmap of markers in doublets vs singlets ---------

library(tidyr)
heat_df <- marker_summary %>%
  mutate(label = paste(Treatment, cell_class, sep = "\n")) %>%
  select(label, gene, mean_expr) %>%
  pivot_wider(names_from = label, values_from = mean_expr)

heat_mat <- as.matrix(heat_df[, -1])
rownames(heat_mat) <- heat_df$gene

# Scale per gene for visibility
heat_scaled <- t(scale(t(heat_mat)))
heat_scaled[is.nan(heat_scaled)] <- 0

pdf(file.path(outdir, "doublet_marker_heatmap.pdf"), width = 6, height = 5)
heatmap(heat_scaled, scale = "none", margins = c(10, 8),
        col = colorRampPalette(c("navy", "white", "firebrick"))(50),
        main = "Marker expression (z-scored)\ndoublets vs singlets")
dev.off()

# ==============================================================================
# 6. CO-EXPRESSION OF T CELL AND TUMOR MARKERS IN DOUBLETS VS SINGLETS
# ==============================================================================

cat("\nComputing co-expression plots...\n")

# Per-cell aggregate scores (sum of UMIs across each marker group)
marker_merged$t_score    <- rowSums(marker_merged[, t_markers])
marker_merged$tumor_score <- rowSums(marker_merged[, tumor_markers])

# Also compute subtype-level scores
marker_merged$t_general    <- rowSums(marker_merged[, t_markers_general])
marker_merged$t_lineage    <- rowSums(marker_merged[, t_markers_lineage])
marker_merged$t_activation <- rowSums(marker_merged[, t_markers_activation])
marker_merged$t_effector   <- rowSums(marker_merged[, t_markers_effector])
marker_merged$t_exhaustion <- rowSums(marker_merged[, t_markers_exhaustion])
marker_merged$tumor_neural  <- rowSums(marker_merged[, tumor_markers_neural])
marker_merged$tumor_glioma  <- rowSums(marker_merged[, tumor_markers_glioma])
marker_merged$tumor_prolif  <- rowSums(marker_merged[, tumor_markers_prolif])
marker_merged$tumor_invasive <- rowSums(marker_merged[, tumor_markers_invasive])
marker_merged$tumor_receptor <- rowSums(marker_merged[, tumor_markers_receptor])

## -- Plot 1: scatter T cell score vs tumor score per cell ---------------------

# Jitter doublets on top; sample singlets for legibility
set.seed(42)
singlets <- marker_merged[!marker_merged$is_doublet, ]
singlets  <- singlets[sample(nrow(singlets), min(3000, nrow(singlets))), ]
doublets  <- marker_merged[marker_merged$is_doublet, ]
plot_df   <- rbind(singlets, doublets)
plot_df$cell_class <- factor(plot_df$cell_class, levels = c("Singlet", "Doublet"))
plot_df   <- plot_df[order(plot_df$cell_class), ]   # doublets on top

p_scatter_coexpr <- ggplot(plot_df,
                            aes(x = log1p(t_score), y = log1p(tumor_score),
                                color = cell_class, alpha = cell_class, size = cell_class)) +
  geom_point() +
  scale_color_manual(values = c("Singlet" = "grey70", "Doublet" = "#FF7F00")) +
  scale_alpha_manual(values = c("Singlet" = 0.3, "Doublet" = 0.85)) +
  scale_size_manual(values  = c("Singlet" = 0.5, "Doublet" = 1.2)) +
  facet_wrap(~ Treatment) +
  labs(title = "T cell score vs GL261 tumor score per cell",
       subtitle = "Singlets (3k sampled) + all doublets",
       x = "log1p(T cell marker UMI sum)", y = "log1p(Tumor marker UMI sum)",
       color = NULL, alpha = NULL, size = NULL) +
  theme_bw(base_size = 12) +
  guides(alpha = "none", size = "none")

pdf(file.path(outdir, "coexpr_scatter.pdf"), width = 9, height = 4.5)
print(p_scatter_coexpr)
dev.off()

## -- Plot 2: subtype score heatmap — doublets vs singlets, by treatment -------

score_cols <- c("t_general", "t_lineage", "t_activation", "t_effector", "t_exhaustion",
                "tumor_neural", "tumor_glioma", "tumor_prolif", "tumor_invasive", "tumor_receptor")
score_labels <- c("T: general", "T: lineage", "T: activation", "T: effector", "T: exhaustion",
                  "Tumor: neural stem", "Tumor: glioma", "Tumor: proliferation",
                  "Tumor: invasive", "Tumor: receptor")

score_summary <- marker_merged %>%
  group_by(Treatment, cell_class) %>%
  summarise(across(all_of(score_cols), ~ round(mean(. > 0) * 100, 2)), .groups = "drop") %>%
  pivot_longer(cols = all_of(score_cols), names_to = "score", values_to = "pct_positive") %>%
  mutate(
    score = factor(score, levels = score_cols, labels = score_labels),
    group = paste(Treatment, cell_class, sep = "\n")
  )

group_levels <- c("aTrem2\nDoublet", "aTrem2\nSinglet", "IgG\nDoublet", "IgG\nSinglet")
score_summary$group <- factor(score_summary$group, levels = group_levels)

p_score_heat <- ggplot(score_summary, aes(x = group, y = score, fill = pct_positive)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_positive, 1)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "#FDD49E", high = "#D7301F",
                       midpoint = 15, limits = c(0, NA),
                       name = "% cells\nexpressing") +
  scale_y_discrete(limits = rev(score_labels)) +
  labs(title = "% cells with any expression in marker group",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

pdf(file.path(outdir, "coexpr_score_heatmap.pdf"), width = 7, height = 5)
print(p_score_heat)
dev.off()

## -- Plot 3: pairwise co-expression — T cell subtype × tumor subtype ----------
# For each cell: does it express ≥1 gene in T cell subtype AND ≥1 in tumor subtype?

t_subtypes    <- list("T: general"    = t_markers_general,
                      "T: lineage"    = t_markers_lineage,
                      "T: activation" = t_markers_activation,
                      "T: effector"   = t_markers_effector,
                      "T: exhaustion" = t_markers_exhaustion)

tumor_subtypes <- list("Tumor: neural stem"   = tumor_markers_neural,
                       "Tumor: glioma"        = tumor_markers_glioma,
                       "Tumor: proliferation" = tumor_markers_prolif,
                       "Tumor: invasive"      = tumor_markers_invasive,
                       "Tumor: receptor"      = tumor_markers_receptor)

coexpr_rows <- lapply(names(t_subtypes), function(tname) {
  t_pos <- rowSums(marker_merged[, t_subtypes[[tname]], drop = FALSE]) > 0
  lapply(names(tumor_subtypes), function(tuname) {
    tu_pos <- rowSums(marker_merged[, tumor_subtypes[[tuname]], drop = FALSE]) > 0
    marker_merged %>%
      mutate(coexpr = t_pos & tu_pos) %>%
      group_by(Treatment, cell_class) %>%
      summarise(pct_coexpr = round(mean(coexpr) * 100, 2), .groups = "drop") %>%
      mutate(t_subtype = tname, tumor_subtype = tuname)
  })
})
coexpr_df <- do.call(rbind, do.call(c, coexpr_rows))

coexpr_df$t_subtype     <- factor(coexpr_df$t_subtype,     levels = names(t_subtypes))
coexpr_df$tumor_subtype <- factor(coexpr_df$tumor_subtype, levels = names(tumor_subtypes))
coexpr_df$group         <- factor(paste(coexpr_df$Treatment, coexpr_df$cell_class, sep = "\n"),
                                   levels = group_levels)

p_coexpr <- ggplot(coexpr_df, aes(x = tumor_subtype, y = t_subtype, size = pct_coexpr,
                                   color = pct_coexpr)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10), name = "% co-expressing") +
  scale_color_gradient(low = "#FEE0D2", high = "#CB181D", name = "% co-expressing") +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "Co-expression: T cell subtype × GL261 tumor subtype",
       x = "Tumor marker group", y = "T cell marker group") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9))

pdf(file.path(outdir, "coexpr_bubble.pdf"), width = 13, height = 5)
print(p_coexpr)
dev.off()

write.csv(coexpr_df, file.path(outdir, "coexpr_summary.csv"), row.names = FALSE)

# ==============================================================================
# 7. METACELL CONSTRUCTION AND ANNOTATION
# ==============================================================================

library(metacell)
library(Matrix)

mc_dir <- file.path(outdir, "mc_db")
dir.create(mc_dir, recursive = TRUE, showWarnings = FALSE)
scdb_init(mc_dir, force_reinit = TRUE)

mc_id <- "aTrem2_IgG"

cat("\nBuilding sparse count matrix from plate files...\n")
mat_list <- lapply(all_plates, function(f) {
  plate_id <- sub(".*_(AB\\d+)\\.txt\\.gz", "\\1", basename(f))
  treatment <- if (grepl("aTrem2", f)) "aTrem2" else "IgG"
  cat(sprintf("  %s\n", plate_id))
  m <- read.table(gzfile(f), header = TRUE, row.names = 1,
                  check.names = FALSE, sep = "\t")
  as(as.matrix(m), "dgCMatrix")
})
# Align genes across plates
all_genes <- Reduce(union, lapply(mat_list, rownames))
mat_list  <- lapply(mat_list, function(m) {
  missing <- setdiff(all_genes, rownames(m))
  if (length(missing) > 0) {
    pad <- Matrix(0, nrow = length(missing), ncol = ncol(m),
                  dimnames = list(missing, colnames(m)), sparse = TRUE)
    m   <- rbind(m, pad)
  }
  m[all_genes, , drop = FALSE]
})
combined_mat <- do.call(cbind, mat_list)

# Keep only cells present in metadata (aTrem2 + IgG)
keep_cells   <- intersect(colnames(combined_mat), meta_sub$Well_ID)
combined_mat <- combined_mat[, keep_cells]
cat(sprintf("Sparse matrix: %d genes x %d cells\n",
            nrow(combined_mat), ncol(combined_mat)))

# Build cell metadata data.frame for metacell
idx      <- match(keep_cells, meta_sub$Well_ID)
cell_md  <- data.frame(
  row.names      = keep_cells,
  Treatment      = meta_sub$Treatment[idx],
  Mouse          = meta_sub$Amp_batch_ID[idx],   # use plate as proxy for mouse
  celltype       = meta_sub$celltype[idx],
  time_assignment = meta_sub$time_assignment[idx],
  stringsAsFactors = FALSE
)
cell_md$time_assignment[is.na(cell_md$time_assignment)] <- "Not_Assigned"

# Create scmat object and add to scdb
scmat <- scm_new_matrix(combined_mat, cell_metadata = cell_md, stat_type = "umi")
scdb_add_mat(mc_id, scmat)

# Filter: ignore mitochondrial, immunoglobulin genes and low-UMI cells
ig_genes <- grep("^mt-|^Igk|^Igh|^Igl", rownames(combined_mat), value = TRUE)
if (length(ig_genes) > 0)
  mcell_mat_ignore_genes(mc_id, mc_id, ig_genes)
mcell_mat_ignore_small_cells(mc_id, mc_id, min_umis = 200)

# Compute gene statistics, then select feature genes
cat("Computing gene statistics...\n")
mcell_add_gene_stat(mc_id, mc_id, force = TRUE)
mcell_gset_filter_varmean(mc_id, mc_id, T_vm = 0.08, force_new = TRUE)
mcell_gset_filter_cov(mc_id, mc_id, T_tot = 100, T_top3 = 2)
mcell_gset_filter_szcor(mc_id, mc_id, T_szcor = -0.1)
cat(sprintf("Feature genes selected: %d\n", length(scdb_gset(mc_id)@gene_set)))

# Build kNN graph -> co-clustering -> metacells
cat("Building kNN graph (K=100)...\n")
mcell_add_cgraph_from_mat_bknn(mc_id, mc_id, mc_id, K = 100)
cat("Co-clustering (500 resamples)...\n")
mcell_coclust_from_graph_resamp(mc_id, mc_id, min_mc_size = 20, p_resamp = 0.75, n_resamp = 500)
cat("Building metacells...\n")
mcell_mc_from_coclust_balanced(mc_id, mc_id, mc_id, K = 30, min_mc_size = 20)

# 2D projection for visualization
cat("Computing 2D projection...\n")
mcell_mc2d_force_knn(mc_id, K = 20)

# Extract mc assignments and annotate per metacell
mc   <- scdb_mc(mc_id)
mc2d <- scdb_mc2d(mc_id)

# Dominant cell type per metacell
mc_cells    <- names(mc@mc)
mc_celltype <- tapply(cell_md[mc_cells, "celltype"], mc@mc[mc_cells],
                      function(x) names(sort(table(x), decreasing = TRUE))[1])

# Per-metacell QC summary: n_cells, dominant celltype, time distribution
mc_qc <- data.frame(
  mc_id_num  = as.integer(names(mc_celltype)),
  n_cells    = as.integer(table(mc@mc)[names(mc_celltype)]),
  celltype   = as.character(mc_celltype),
  mc_x       = mc2d@mc_x[as.integer(names(mc_celltype))],
  mc_y       = mc2d@mc_y[as.integer(names(mc_celltype))],
  stringsAsFactors = FALSE
)

# Fraction of cells per metacell assigned to each time bin
time_bins <- c("12H", "24H", "36H", "Negative", "Not_Assigned")
mc_time   <- do.call(rbind, lapply(sort(unique(mc@mc)), function(m) {
  ta <- cell_md[names(mc@mc)[mc@mc == m], "time_assignment"]
  ta[is.na(ta)] <- "Not_Assigned"
  tb <- table(factor(ta, levels = time_bins))
  as.data.frame(t(as.matrix(tb / sum(tb))), stringsAsFactors = FALSE)
}))
mc_qc <- cbind(mc_qc, mc_time)

write.csv(mc_qc, file.path(outdir, "metacell_qc.csv"), row.names = FALSE)
cat("Wrote metacell QC to", file.path(outdir, "metacell_qc.csv"), "\n")
cat(sprintf("  %d metacells, median %d cells/mc\n",
            nrow(mc_qc), median(mc_qc$n_cells)))

# Add per-cell mc assignment to merged (NA for cells filtered out by metacell)
mc_assignment <- data.frame(
  Well_ID   = names(mc@mc),
  mc_id_num = as.integer(mc@mc),
  stringsAsFactors = FALSE
)
# Add mc 2D coords per cell (from mc2d@sc_x / sc_y)
mc_assignment$mc_sc_x <- mc2d@sc_x[names(mc@mc)]
mc_assignment$mc_sc_y <- mc2d@sc_y[names(mc@mc)]
# Add dominant celltype of the metacell each cell belongs to
mc_assignment$mc_celltype <- mc_celltype[as.character(mc_assignment$mc_id_num)]

merged <- merge(merged, mc_assignment, by = "Well_ID", all.x = TRUE)

# ==============================================================================
# 8. EXPORT FOR SCANPY / ANNDATA
# ==============================================================================

cat("\nExporting AnnData-compatible files for Python/scanpy...\n")

# Attach QC metrics from merged (now includes mc columns)
qc_cols   <- c("Well_ID", "total_UMI", "genes_det", "time_assignment",
               "sc_x", "sc_y", "mc_id_num", "mc_sc_x", "mc_sc_y", "mc_celltype")
export_df <- merge(marker_merged,
                   merged[, qc_cols],
                   by = "Well_ID", all.x = TRUE)

# -- Observation metadata (one row per cell) --
obs_cols <- c("Well_ID", "Amp_batch_ID", "Treatment", "celltype",
              "is_doublet", "cell_class", "time_assignment",
              "sc_x", "sc_y", "total_UMI", "genes_det",
              "mc_id_num", "mc_sc_x", "mc_sc_y", "mc_celltype")
write.csv(export_df[, obs_cols],
          file.path(outdir, "cells_obs.csv"), row.names = FALSE)

# -- Marker gene count matrix (cells x genes) --
write.csv(export_df[, c("Well_ID", all_markers)],
          file.path(outdir, "marker_counts.csv"), row.names = FALSE)

cat("Exported:\n")
cat(" ", file.path(outdir, "cells_obs.csv"), "\n")
cat(" ", file.path(outdir, "marker_counts.csv"), "\n")
cat(" ", file.path(outdir, "metacell_qc.csv"), "\n")
cat("Run 01_summary_aTrem2.py to build cells_markers.h5ad and reproduce figures.\n")

cat("Done. Outputs in", outdir, "\n")
