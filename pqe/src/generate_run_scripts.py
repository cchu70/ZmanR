#!/usr/bin/env python3
"""Generate all pseudotime run scripts for scumi and mcumi adata × 6 gene sets × 5 tools."""
import os

RESULTS = "/home/unix/cchu/projects/ZmanR/pqe/results"
REPO    = os.path.join(RESULTS, "06")
SRC     = os.path.dirname(os.path.abspath(__file__))
SPEARMAN_DIR = os.path.join(REPO, "ctet_mos/spearman")
ZMAN_S3_ATREM = "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_aTREM2_Time.csv"
ZMAN_S3_IGG   = "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_Isotype_Control_Time.csv"

# Original adata paths (R scripts use these + RDS cache)
ADATA_PATHS_R = {
    "scumi": os.path.join(REPO, "zmanseq_momac_metacells_annot_clean.h5ad"),
    "mcumi": os.path.join(REPO, "zmanseq_momac_metacells_annot_clean_mcumi.h5ad"),
}
# Fixed adata paths for Python (Palantir): uns/log1p/base removed, highly_variable added
ADATA_PATHS_PY = {
    "scumi": os.path.join(REPO, "zmanseq_momac_metacells_annot_clean_fixed.h5ad"),
    "mcumi": os.path.join(REPO, "zmanseq_momac_metacells_annot_clean_mcumi_fixed.h5ad"),
}
ADATA_PATHS = ADATA_PATHS_R  # used by R scripts

# Gene set definitions: (short_label, human_label, r_gene_loader, py_gene_loader)
# r_gene_loader: R code snippet that assigns `genes`
# py_gene_loader: Python code snippet that assigns `genes`

def r_ctet_mos(label, geneset, timecol):
    return f'genes <- load_ctet_mos_genes("{label}", geneset="{geneset}", timecol="{timecol}")'

def py_ctet_mos(label, geneset, timecol):
    return (
        f'igg_df   = pd.read_csv(f"{SPEARMAN_DIR}/{label}_igg_{geneset}_{timecol}_spearman.csv",   index_col=0)\n'
        f'atrem_df = pd.read_csv(f"{SPEARMAN_DIR}/{label}_atrem_{geneset}_{timecol}_spearman.csv", index_col=0)\n'
        f'genes = list(set(igg_df.index[igg_df["p_value"] < 0.05]) & set(atrem_df.index[atrem_df["p_value"] < 0.05]))\n'
        f'print(f"{label} {geneset} {timecol} gene intersection: {{len(genes)}} genes")'
    )

GENE_SETS = {
    "hvg_ctet":        ("hvg_ctet",        lambda lbl: r_ctet_mos(lbl, "hvg",  "ctet"),        lambda lbl: py_ctet_mos(lbl, "hvg",  "ctet")),
    "hvg_smooth_ctet": ("hvg_smooth_ctet", lambda lbl: r_ctet_mos(lbl, "hvg",  "smooth_ctet"), lambda lbl: py_ctet_mos(lbl, "hvg",  "smooth_ctet")),
    "full_ctet":       ("full_ctet",       lambda lbl: r_ctet_mos(lbl, "full", "ctet"),        lambda lbl: py_ctet_mos(lbl, "full", "ctet")),
    "full_smooth_ctet":("full_smooth_ctet",lambda lbl: r_ctet_mos(lbl, "full", "smooth_ctet"), lambda lbl: py_ctet_mos(lbl, "full", "smooth_ctet")),
    "hvg":             ("hvg",
                        lambda lbl: 'genes <- load_hvg_genes(ADATA_PATH)',
                        lambda lbl: 'import scanpy as _sc2; _sc2.pp.highly_variable_genes(adata_full, n_top_genes=1000)\ngenes = list(adata_full.var_names[adata_full.var["highly_variable"]])\nprint(f"HVG genes: {len(genes)}")'),
    "zman_s3":         ("zman_s3",
                        lambda lbl: 'genes <- load_zman_s3_genes()',
                        lambda lbl: (f'atrem_s3 = pd.read_csv("{ZMAN_S3_ATREM}")\n'
                                     f'igg_s3   = pd.read_csv("{ZMAN_S3_IGG}")\n'
                                     f'genes = list(set(atrem_s3.loc[atrem_s3["Pvalue"] < 0.05, "Gene"]) & '
                                     f'set(igg_s3.loc[igg_s3["Pvalue"] < 0.05, "Gene"]))\n'
                                     f'print(f"Zman S3 gene intersection: {{len(genes)}} genes")')),
}

# ── R script templates ────────────────────────────────────────────────────────

def write_dpt_r(label, geneset_key, adata_path, out_dir, gene_loader_r):
    fname = os.path.join(SRC, f"run_dpt_{label}_{geneset_key}.R")
    content = f"""#!/usr/bin/env Rscript
ADATA_PATH <- "{adata_path}"
OUT_PLOT   <- "{RESULTS}/{out_dir}/dpt_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(destiny)

{gene_loader_r}
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Building DiffusionMap...\\n"); set.seed(42)
dm  <- DiffusionMap(dat$expr, n_pcs = min(30, ncol(dat$expr) - 1), n_eigs = 10, verbose = TRUE)
start_idx  <- which.min(dat$coldata$cTET)
atrem_idxs <- which(dat$coldata$enrichment > 0)
igg_idxs   <- which(dat$coldata$enrichment < 0)
tip_atrem  <- atrem_idxs[which.max(dat$coldata$cTET[atrem_idxs])]
tip_igg    <- igg_idxs[which.max(dat$coldata$cTET[igg_idxs])]
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\\n")
cat("Tip aTrem2:", dat$cell_names[tip_atrem], "(cTET =", dat$coldata$cTET[tip_atrem], ")\\n")
cat("Tip IgG:  ", dat$cell_names[tip_igg],   "(cTET =", dat$coldata$cTET[tip_igg],   ")\\n")
cat("Computing DPT...\\n")
dpt <- DPT(dm, tips = c(start_idx, tip_atrem, tip_igg))
dpt_mat  <- as.matrix(dpt)
pt_root  <- dpt_mat[, start_idx]
pt_atrem <- dpt_mat[, tip_atrem]
pt_igg   <- dpt_mat[, tip_igg]
dpt_cols <- data.frame(DPT_root = pt_root, DPT_atrem = pt_atrem, DPT_igg = pt_igg)
save_pseudotime_plot(dat$sc_x, dat$sc_y, pt_root, OUT_PLOT, "DPT pseudotime ({label} {geneset_key})")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pt_root, sub("[.]png$", ".tsv", OUT_PLOT), extra_cols = dpt_cols)
"""
    with open(fname, "w") as f:
        f.write(content)
    print(f"  Wrote {fname}")

def write_monocle2_r(label, geneset_key, adata_path, out_dir, gene_loader_r):
    fname = os.path.join(SRC, f"run_monocle2_{label}_{geneset_key}.R")
    content = f"""#!/usr/bin/env Rscript
ADATA_PATH <- "{adata_path}"
OUT_PLOT   <- "{RESULTS}/{out_dir}/monocle2_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(monocle)

{gene_loader_r}
dat   <- load_adata(ADATA_PATH, gene_list = genes)
expr_mat <- t(dat$expr)
pd <- new("AnnotatedDataFrame", data = data.frame(row.names = dat$cell_names,
          sc_x = dat$sc_x, sc_y = dat$sc_y))
fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(expr_mat),
          gene_short_name = rownames(expr_mat)))
cat("Building CellDataSet...\\n")
cds <- newCellDataSet(expr_mat, phenoData = pd, featureData = fd, expressionFamily = tobit())
cds <- estimateSizeFactors(cds)
cds <- setOrderingFilter(cds, rownames(expr_mat))
cat("Running Monocle2 DDRTree...\\n")
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
start_cell <- dat$cell_names[which.min(dat$coldata$cTET)]
root_state <- pData(cds)[start_cell, "State"]
cat("Start cell:", start_cell, "root state:", as.character(root_state), "\\n")
cds <- orderCells(cds, root_state = root_state)
pseudotime <- pData(cds)$Pseudotime
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "Monocle2 pseudotime ({label} {geneset_key})")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pseudotime, sub("[.]png$", ".tsv", OUT_PLOT))
"""
    with open(fname, "w") as f:
        f.write(content)
    print(f"  Wrote {fname}")

def write_redpath_r(label, geneset_key, adata_path, out_dir, gene_loader_r):
    fname = os.path.join(SRC, f"run_redpath_{label}_{geneset_key}.R")
    content = f"""#!/usr/bin/env Rscript
ADATA_PATH <- "{adata_path}"
OUT_PLOT   <- "{RESULTS}/{out_dir}/redpath_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(redPATH)

{gene_loader_r}
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Genes:", ncol(dat$expr), "\\n")
cat("Running redPATH...\\n")
pseudotime <- redpath(dat$expr, threadnum = 4, base_path_range = c(2:6))
start_idx <- which.min(dat$coldata$cTET)
if (pseudotime[start_idx] > median(pseudotime, na.rm = TRUE)) pseudotime <- max(pseudotime, na.rm = TRUE) - pseudotime + min(pseudotime, na.rm = TRUE)
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\\n")
save_pseudotime_plot(dat$sc_x, dat$sc_y, pseudotime, OUT_PLOT, "redPATH pseudotime ({label} {geneset_key})")
save_pseudotime_tsv(dat$cell_names, dat$coldata, pseudotime, sub("[.]png$", ".tsv", OUT_PLOT))
"""
    with open(fname, "w") as f:
        f.write(content)
    print(f"  Wrote {fname}")

def write_scorpius_r(label, geneset_key, adata_path, out_dir, gene_loader_r):
    fname = os.path.join(SRC, f"run_scorpius_{label}_{geneset_key}.R")
    content = f"""#!/usr/bin/env Rscript
ADATA_PATH <- "{adata_path}"
OUT_PLOT   <- "{RESULTS}/{out_dir}/scorpius_pseudotime.png"

.sd <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(.sd, "utils.R")); source(file.path(.sd, "utils_momac.R"))
library(SCORPIUS)

{gene_loader_r}
dat   <- load_adata(ADATA_PATH, gene_list = genes)
cat("Running SCORPIUS...\\n")
space <- reduce_dimensionality(dat$expr, dist = "spearman", ndim = 3)
traj  <- infer_trajectory(space)
start_idx <- which.min(dat$coldata$cTET)
if (traj$time[start_idx] > 0.5) traj$time <- 1 - traj$time
cat("Start cell:", dat$cell_names[start_idx], "(cTET =", dat$coldata$cTET[start_idx], ")\\n")
save_pseudotime_plot(dat$sc_x, dat$sc_y, traj$time, OUT_PLOT, "SCORPIUS pseudotime ({label} {geneset_key})")
save_pseudotime_tsv(dat$cell_names, dat$coldata, traj$time, sub("[.]png$", ".tsv", OUT_PLOT))
"""
    with open(fname, "w") as f:
        f.write(content)
    print(f"  Wrote {fname}")

def write_palantir_py(label, geneset_key, adata_path, adata_path_py, out_dir, gene_loader_py):
    base = adata_path.replace(".h5ad", "")
    fname = os.path.join(SRC, f"run_palantir_{label}_{geneset_key}.py")
    content = f"""\"\"\"Run Palantir pseudotime on {label} metacells using {geneset_key} gene set.\"\"\"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, scanpy as sc, palantir

# Load adata from numpy/CSV exports (avoids anndata version compatibility issues)
_X   = np.load("{base}_X.npy")
_obs = pd.read_csv("{base}_obs.csv", index_col=0)
_var = pd.read_csv("{base}_var.csv", index_col=0)
adata_full = sc.AnnData(X=_X, obs=_obs, var=_var)
OUT_PLOT   = "{RESULTS}/{out_dir}/palantir_pseudotime.png"

{gene_loader_py}

genes_in = [g for g in genes if g in adata_full.var_names]
print(f"Genes present in adata: {{len(genes_in)}}")
adata = adata_full[:, genes_in].copy()

# Mark all genes as highly_variable so run_pca uses the full gene set
adata.var["highly_variable"] = True

n_pcs = min(10, len(genes_in) - 1)
palantir.utils.run_pca(adata, n_components=n_pcs)
palantir.utils.run_diffusion_maps(adata, n_components=min(5, n_pcs))
palantir.utils.determine_multiscale_space(adata)

start_cell = adata.obs["cTET"].idxmin()
print(f"Start cell: {{start_cell}} (cTET = {{adata.obs.loc[start_cell, 'cTET']:.4f}})")

pr_res = palantir.core.run_palantir(adata, start_cell, num_waypoints=40, knn=10, use_early_cell_as_start=True)
adata.obs["palantir_pseudotime"] = pr_res.pseudotime

sc_x, sc_y, pt = adata.obs["sc_x"].astype(float), adata.obs["sc_y"].astype(float), adata.obs["palantir_pseudotime"].astype(float)
fig, ax = plt.subplots(figsize=(6, 5))
sc_plot = ax.scatter(sc_x, sc_y, c=pt, cmap="plasma", s=8, alpha=0.8, linewidths=0)
plt.colorbar(sc_plot, ax=ax, label="Pseudotime")
ax.set_title("Palantir pseudotime ({label} {geneset_key})"); ax.set_xlabel("sc_x"); ax.set_ylabel("sc_y")
plt.tight_layout(); plt.savefig(OUT_PLOT, dpi=150); print(f"Plot saved to {{OUT_PLOT}}")

out_tsv = OUT_PLOT.replace(".png", ".tsv")
adata.obs.to_csv(out_tsv, sep="\\t"); print(f"TSV saved to {{out_tsv}}")
"""
    with open(fname, "w") as f:
        f.write(content)
    print(f"  Wrote {fname}")

# ── Generate all scripts ──────────────────────────────────────────────────────

for label in ADATA_PATHS_R:
    adata_path_r  = ADATA_PATHS_R[label]
    adata_path_py = ADATA_PATHS_PY[label]
    for geneset_key, (_, r_loader_fn, py_loader_fn) in GENE_SETS.items():
        out_dir = f"pseudotime_{label}_{geneset_key}"
        os.makedirs(os.path.join(RESULTS, out_dir), exist_ok=True)
        r_code  = r_loader_fn(label)
        py_code = py_loader_fn(label)
        write_dpt_r(label, geneset_key, adata_path_r, out_dir, r_code)
        write_monocle2_r(label, geneset_key, adata_path_r, out_dir, r_code)
        write_redpath_r(label, geneset_key, adata_path_r, out_dir, r_code)
        write_scorpius_r(label, geneset_key, adata_path_r, out_dir, r_code)
        write_palantir_py(label, geneset_key, adata_path_r, adata_path_py, out_dir, py_code)

print("\nDone. Created 60 run scripts (12 variants × 5 tools).")
