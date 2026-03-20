library(ZmanR)
library(metacell)
library(ggplot2)
library(ggpubr)

# Output directory for figures
dir.create("results/00", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PART 1: TIME ASSIGNMENT
# ==============================================================================

# --- GLM model evaluation ---
# Check whether the fluorescence signal fits a normal distribution based on
# unstained cells. Outputs Q-Q, fitted vs residual, etc. plots.
pdf("results/00/glm_fit.pdf", width = 10, height = 8)
FACS_model_eval(well_fcs_mc_Blood,
                fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"))
dev.off()

# --- Fluorophore distribution: unstained vs stained ---
pdf("results/00/fluo_hist.pdf", width = 10, height = 6)
plot_FACS_histogram(well_fcs_mc_Blood,
                    fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"))
dev.off()

# --- Assign time bins to cells ---
# sd_threshold: SD cutoff per fluorophore for calling a cell "stained"
# timebins: time labels corresponding to each fluorophore
well_fcs_mc_Blood <- FACS_model(well_fcs_mc_Blood,
                                sd_threshold = c(2.5, 1.5, 2.5, 2),
                                fluorophores  = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"),
                                timebins      = c("12H", "24H", "36H", "48H"))

# --- Evaluate time-bin assignments ---
# Contour shows GLM classification boundary; points colored by assigned timebin
pdf("results/00/stain_assignment.pdf", width = 12, height = 8)
plot_FACS_eval(well_fcs_mc_Blood,
               fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"),
               timebins     = c("12H", "24H", "36H", "48H"),
               y_axis       = "log_APC.A")
dev.off()

# ==============================================================================
# PART 2: TRAJECTORY ANALYSIS
# ==============================================================================

# Load the metacell database for the GBM T/NK example (Figure 3 of paper)
# GBM_T_mc_annotations: columns "mc_id" and "celltype"
# well_fcs_GBM_T_time: column "groundtruth_group" with assigned time bin
scdb_init(file.path(dirname(getwd()), "mc_example"))
new_id <- "T_clean"

# --- Compute metacell CDF and AUC time ---
# select_mcs: T/NK metacells (1:37); kept separate from myeloids
# time_for_auc: actual injection times in hours (0 = unlabeled baseline)
mc_cdf <- compute_mc_cdf(new_id,
                         well_fcs_GBM_T_time,
                         select_mcs     = 1:37,
                         mc_annotations = GBM_T_mc_annotations,
                         time_points    = c("12H", "24H", "36H"),
                         time_for_auc   = c(0, 12, 24, 36),
                         Gate           = "CD45_high")

# --- Plot raw CDF curves colored by AUC, faceted by cell type ---
# Red dashed line = linear reference (intercept/slope from GLM fit)
p_auc <- ggplot(mc_cdf$mc_auc_time,
                aes(x = time, y = cums, group = mc, color = auc)) +
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~ factor(celltype,
                      levels = c("Treg", "CD4", "CD8",
                                 "chemotactic", "cytotoxic",
                                 "intermediate", "dysfunctional")),
             ncol = 7) +
  theme(text = element_text(size = 20), legend.position = "none") +
  ylab("% of cells") +
  xlab("time (hours)") +
  geom_abline(intercept = 0.01412, slope = 0.02810,
              color = "red", linetype = "dashed", alpha = 0.5)

ggsave("results/00/celltype_auc.pdf", p_auc, width = 21, height = 3.5)

# --- Normalize expression for NK/cytotoxic cell types ---
NK_mc_exprs <- normalize_mc_exprs(new_id,
                                  mc_annotations  = GBM_T_mc_annotations,
                                  mc_cdf$mc2d_auc_time,
                                  select_celltypes = c("chemotactic", "cytotoxic",
                                                       "intermediate", "dysfunctional"))

# Subset to GO genes and filter by dispersion
select_GO_mc  <- get_GO_exp(NK_mc_exprs$select_exprs,
                             gene_type = "gene_names",
                             organism  = "mouse",
                             takelog   = FALSE)
filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold = 0.05, threads = 1)

# --- Smooth Zman trajectory ---
# ref_k = number of distinct cell type annotations in the data
NK_smoothed_res <- smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs, ref_k = 4)

# --- Gene-time correlations (Spearman) ---
nk_corr <- calculate_corr_genes(new_id, NK_smoothed_res, "spearman")

# Print correlations and p-values for key genes
key_genes <- c("Tigit", "Xcl1", "Pmepa1", "Igflr1", "Gzmc")
cat("Spearman correlations:\n")
print(signif(nk_corr$correlation[match(key_genes, names(nk_corr$correlation))], 3))
cat("P-values:\n")
print(signif(nk_corr$pvalues[match(key_genes, names(nk_corr$correlation))], 3))

# --- Predict and plot smoothed expression along time ---
NK_predicted_res <- predict_expression_along_time(new_id,
                                                  NK_smoothed_res,
                                                  out_len     = 36,
                                                  loess_degree = 1,
                                                  loess_span  = 0.9)

traj_plot <- plot_smoothed_trajectory(NK_smoothed_res, mc_cdf)
gene_plot <- plot_zman_genes_heatmap(NK_predicted_res,
                                     NK_smoothed_res,
                                     up_regulated_genes   = c("Tigit", "Xcl1", "Pmepa1", "Igflr1", "Gzmc"),
                                     down_regulated_genes = c("Ccl3", "Prf1", "Gzma", "Gzmb", "Nkg7"),
                                     k = 2)

pdf("results/00/traj_plot.pdf", width = 8, height = 5)
print(traj_plot)
dev.off()

pdf("results/00/gene_heatmap.pdf", width = 6, height = 5)
ComplexHeatmap::draw(gene_plot)
dev.off()

pdf("results/00/traj_plots.pdf", width = 14, height = 5)
p_combined <- ggarrange(traj_plot,
                        ggplotify::as.ggplot(gene_plot),
                        widths = c(8, 6))
print(p_combined)
dev.off()

cat("Done. Figures saved to figures/\n")
