#!/usr/bin/env python
"""
06_ctet_mos.py
Clean annotation + cTET pipeline for myeloid metacells.

Produces 4 h5ad objects in results/06/:
  1. zmanseq_myeloid_metacells_annot_clean.h5ad
        All 72 myeloid metacells; downsampled to median single-cell UMI;
        cTET normalized across all myeloid metacells.

  2. zmanseq_myeloid_metacells_annot_clean_mcumi.h5ad
        All 72 myeloid metacells; downsampled to mean metacell UMI;
        cTET normalized across all myeloid metacells.

  3. zmanseq_momac_metacells_annot_clean.h5ad          (replaces existing)
        Mo-type metacells only; downsampled to median single-cell UMI;
        cTET renormalized within mo-type metacells + smooth_cTET.

  4. zmanseq_momac_metacells_annot_clean_mcumi.h5ad
        Mo-type metacells only; downsampled to mean metacell UMI;
        cTET renormalized within mo-type metacells + smooth_cTET.

AUC is calculated once on the full madata; norm_AUC / cTET and smooth_cTET
are subset-dependent and differ between the all-myeloid and mo-only objects.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import spearmanr
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage
from sklearn import metrics
from sklearn.cluster import AgglomerativeClustering
from matplotlib.colors import Normalize
from pathlib import Path

DATA_DIR = Path("results/06")
OUT_DIR  = DATA_DIR

MYELOID_CELL_TYPE_COLORS = {
    'Arg1_TAM':   "mediumorchid",
    'Acp5_TAM':   "slateblue",
    'MoMac1':     "tab:blue",
    'MoMac2':     "tab:red",
    'Monocytes':  "tab:orange",
    'cDC2':       "goldenrod",
    'MonDC':      "gold",
    'Gpnmb_TAM':  "lightblue",
    'cDC1':       "darkkhaki",
    'MigDC':      "tab:brown",
    'transitory': "teal",
    'pDC':        "tab:green",
    'nan':        "gray",
}

MO_TYPES    = ["Monocytes", "MoMac1", "MoMac2", "Gpnmb_TAM", "Arg1_TAM", "Acp5_TAM"]
ATREM_TYPES = ["Monocytes", "MoMac1", "Acp5_TAM"]
IGG_TYPES   = ["Monocytes", "MoMac2", "Gpnmb_TAM", "Arg1_TAM"]
TIME_BINS   = ["12H", "24H", "36H", "48H"]
CUM_BINS    = ["0H", "12H", "24H", "36H"]


# ══════════════════════════════════════════════════════════════════════════════
# 1 — Load data
# ══════════════════════════════════════════════════════════════════════════════

print("Loading data ...")
adata  = sc.read_h5ad(DATA_DIR / "zmanseq_myeloid_cells.h5ad")
madata = sc.read_h5ad(DATA_DIR / "zmanseq_myeloid_metacells.h5ad")
print(f"  Single cells: {adata.shape}")
print(f"  Metacells:    {madata.shape}")

# Downsampling targets
sc_umi_target = int(adata.obs['total_counts'].median())
mc_umi_target = int(madata.obs['total_umis'].mean())
print(f"  Downsample target (median sc UMI):  {sc_umi_target:,}")
print(f"  Downsample target (mean MC UMI):    {mc_umi_target:,}")


# ══════════════════════════════════════════════════════════════════════════════
# 2 — Annotation (computed once on full madata)
# ══════════════════════════════════════════════════════════════════════════════

print("\nAnnotating metacells ...")

# Drop outlier cells (metacell_name == "Outliers") before any groupby
cells_obs = adata.obs[adata.obs["metacell_name"].isin(madata.obs_names)].copy()
# Reset categorical levels so groupby doesn't produce phantom "Outliers" groups
if hasattr(cells_obs["metacell_name"], "cat"):
    cells_obs["metacell_name"] = cells_obs["metacell_name"].cat.remove_unused_categories()

# ── 2a. Cell type (majority vote from single cells) ───────────────────────────
mc_ct_counts = cells_obs.groupby(["metacell_name", "cluster_colors"]).size().unstack(fill_value=0)
mc_ct_pct    = mc_ct_counts.div(mc_ct_counts.sum(axis=1), axis=0)
madata.obs["metacell_cell_type"]     = mc_ct_pct.idxmax(axis=1).loc[madata.obs_names]
madata.obs["metacell_cell_type_pct"] = mc_ct_pct.max(axis=1).loc[madata.obs_names]

# ── 2b. Spatial position (median single-cell coordinates) ─────────────────────
madata.obs["sc_x"] = cells_obs.groupby("metacell_name")["sc_x"].median()
madata.obs["sc_y"] = cells_obs.groupby("metacell_name")["sc_y"].median()

# ── 2c. Treatment enrichment ──────────────────────────────────────────────────
trt_counts = cells_obs.groupby(["metacell_name", "Treatment"]).size().unstack(fill_value=0)
trt_counts["enrichment"] = 0.5 - (trt_counts["aTrem2"] / trt_counts.sum(axis=1))
madata.obs[trt_counts.columns] = trt_counts.loc[madata.obs_names].values

# ── 2d. Cumulative time-bin proportions ───────────────────────────────────────
tb_counts  = cells_obs.groupby(["metacell_name", "time_assignment"]).size().unstack(fill_value=0)
tb_pct     = tb_counts[TIME_BINS].div(tb_counts[TIME_BINS].sum(axis=1), axis=0)
cum_pct_df = tb_pct.cumsum(axis=1).rename(columns={"12H": "12H", "24H": "24H", "36H": "36H"})
cum_pct_df["0H"] = 0.0
cum_pct_df = cum_pct_df[CUM_BINS]
cum_pct_df["metacell_cell_type"] = madata.obs.loc[cum_pct_df.index, "metacell_cell_type"]


# ══════════════════════════════════════════════════════════════════════════════
# 3 — AUC (once, on all metacells)
# ══════════════════════════════════════════════════════════════════════════════

print("Computing AUC ...")
cum_pct_df["AUC"] = cum_pct_df.apply(
    lambda r: metrics.auc(np.arange(len(CUM_BINS)), r[CUM_BINS].values.astype(float)),
    axis=1,
)


# ══════════════════════════════════════════════════════════════════════════════
# 4 — norm_AUC / cTET (subset-dependent)
# ══════════════════════════════════════════════════════════════════════════════

def add_ctet(df):
    """Normalise AUC within df rows and compute cTET = 1 - norm_AUC."""
    auc = df["AUC"].astype(float)
    df = df.copy()
    df["norm_AUC"] = (auc - auc.min()) / (auc.max() - auc.min())
    df["cTET"]     = 1 - df["norm_AUC"]
    return df

# All myeloid
all_cum_pct_df = add_ctet(cum_pct_df)

# Mo types only — renormalize within subset
mo_mask        = cum_pct_df["metacell_cell_type"].isin(MO_TYPES)
mo_cum_pct_df  = add_ctet(cum_pct_df[mo_mask].copy())


# ══════════════════════════════════════════════════════════════════════════════
# 5 — Raw summed counts layer (same for all 4 objects)
# ══════════════════════════════════════════════════════════════════════════════

print("Summing raw counts per metacell ...")
mc_raw = np.zeros((madata.n_obs, madata.n_vars))
for i, mc in enumerate(madata.obs_names):
    cells_in_mc = cells_obs.query("metacell_name == @mc").index
    mc_raw[i, :] = np.asarray(adata[cells_in_mc].X.sum(axis=0)).ravel()

madata.layers["freq"]       = madata.X          # normalized frequencies
madata.layers["total_umis"] = sp.csr_matrix(mc_raw)
madata.layers["zeros"]      = madata.layers.get("zeros", madata.layers["total_umis"])


# ══════════════════════════════════════════════════════════════════════════════
# 6 — smooth_cTET (mo types only, per treatment arm)
# ══════════════════════════════════════════════════════════════════════════════

def get_distance_matrix(mc_data, gene_mask):
    # following redPath 
    rho = spearmanr(mc_data[:, gene_mask].X.T.toarray())
    d_rho = (rho.statistic + 1) / 2
    row_mean   = d_rho.mean(axis=1, keepdims=True)
    col_mean   = d_rho.mean(axis=0, keepdims=True)
    grand_mean = d_rho.mean()
    return d_rho - row_mean - col_mean + grand_mean

def smooth_auc(mc_data, gene_mask, a=0.5, b=0.5, time_col="cTET",
               cell_type_col="metacell_cell_type"):
    N     = mc_data.shape[0]
    D     = get_distance_matrix(mc_data, gene_mask)
    ref_K = mc_data.obs[cell_type_col].nunique()
    ks    = range(ref_K + 1, N - 1)
    k_cols = [f"k={k}" for k in ks]
    df = pd.DataFrame(index=mc_data.obs_names, columns=[f"ref_k={ref_K}"] + k_cols)
    df[cell_type_col] = mc_data.obs[cell_type_col].values
    df[time_col]      = mc_data.obs[time_col].values
    # cell-type mean as anchor
    ct_avg = df.groupby(cell_type_col)[time_col].mean()
    df[f"ref_k={ref_K}"] = df[cell_type_col].map(ct_avg)
    # kNN smoothing
    for k in ks:
        labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit(D).labels_
        df[f"k={k}"] = pd.Series(labels, index=mc_data.obs_names)
        k_avg = df.groupby(f"k={k}")[time_col].mean()
        df[f"k={k}"] = df[f"k={k}"].map(k_avg)
    df[f"smooth_{time_col}"] = (
        b * df[k_cols + [time_col]].astype(float).mean(axis=1) +
        a * df[f"ref_k={ref_K}"].astype(float)
    )
    return df


def compute_smooth_ctet(madata_proc, mo_cum_pct_df):
    """
    Compute smooth_cTET on a processed (log-normalised) madata.
    Requires highly_variable_genes to have been run on madata_proc.
    Returns a Series indexed by obs_names with smooth_cTET values (NaN for non-mo cells).
    """
    sc.pp.highly_variable_genes(madata_proc, n_top_genes=1000)
    hvg_mask = madata_proc.var["highly_variable"].values

    smooth = pd.Series(np.nan, index=madata_proc.obs_names)
    for arm_types, arm_name in [(ATREM_TYPES, "aTrem2"), (IGG_TYPES, "IgG")]:
        arm_idx  = mo_cum_pct_df.index[mo_cum_pct_df["metacell_cell_type"].isin(arm_types)]
        arm_idx  = arm_idx.intersection(madata_proc.obs_names)
        arm_mc   = madata_proc[arm_idx].copy()
        # temporarily attach cTET
        arm_mc.obs["cTET"] = mo_cum_pct_df.loc[arm_idx, "cTET"].astype(float)
        arm_mc.obs["metacell_cell_type"] = mo_cum_pct_df.loc[arm_idx, "metacell_cell_type"]
        sdf = smooth_auc(arm_mc, hvg_mask)
        smooth.loc[arm_idx] = sdf["smooth_cTET"].values
    return smooth


# ══════════════════════════════════════════════════════════════════════════════
# 7 — Build + save each h5ad
# ══════════════════════════════════════════════════════════════════════════════

def build_anndata(madata_base, cum_df, downsample_target, label):
    """
    Build a processed AnnData from raw counts.
      madata_base    : base AnnData with layers['total_umis'] and annotations
      cum_df         : cum_pct_df with AUC, norm_AUC, cTET columns for this subset
      downsample_target : int, counts_per_cell for sc.pp.downsample_counts
      label          : string used in print messages
    """
    print(f"\n  Building {label} (n={madata_base.n_obs}, target={downsample_target:,}) ...")
    ad = madata_base.copy()
    ad.X = ad.layers["total_umis"].copy()

    sc.pp.downsample_counts(ad, counts_per_cell=downsample_target, copy=False)
    sc.pp.calculate_qc_metrics(ad, inplace=True, log1p=True)
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)

    # Attach cTET columns
    for col in CUM_BINS + ["AUC", "norm_AUC", "cTET"]:
        ad.obs[col] = cum_df.reindex(ad.obs_names)[col].astype(float).values

    return ad


print("\nBuilding h5ad objects ...")

# ── 7a. All myeloid — annotate madata ─────────────────────────────────────────
for col in CUM_BINS + ["AUC", "norm_AUC", "cTET", "metacell_cell_type"]:
    madata.obs[col] = all_cum_pct_df.reindex(madata.obs_names)[col].values

# Object 1: all myeloid, median sc UMI
ad_all_scumi = build_anndata(madata, all_cum_pct_df, sc_umi_target, "all-myeloid / sc UMI")
out1 = OUT_DIR / "zmanseq_myeloid_metacells_annot_clean.h5ad"
ad_all_scumi.write_h5ad(out1)
print(f"  Saved {out1.name}")

# Object 2: all myeloid, mean MC UMI
ad_all_mcumi = build_anndata(madata, all_cum_pct_df, mc_umi_target, "all-myeloid / MC UMI")
out2 = OUT_DIR / "zmanseq_myeloid_metacells_annot_clean_mcumi.h5ad"
ad_all_mcumi.write_h5ad(out2)
print(f"  Saved {out2.name}")

# ── 7b. Mo types only ─────────────────────────────────────────────────────────
mo_madata = madata[madata.obs["metacell_cell_type"].isin(MO_TYPES)].copy()
# Replace cTET with mo-renormalized values
for col in ["norm_AUC", "cTET"]:
    mo_madata.obs[col] = mo_cum_pct_df.reindex(mo_madata.obs_names)[col].astype(float).values

# Object 3: mo only, median sc UMI + smooth_cTET
ad_mo_scumi = build_anndata(mo_madata, mo_cum_pct_df, sc_umi_target, "mo-only / sc UMI")
print("  Computing smooth_cTET (sc UMI) ...")
ad_mo_scumi.obs["smooth_cTET"] = compute_smooth_ctet(ad_mo_scumi.copy(), mo_cum_pct_df)
ad_mo_scumi.obs["smooth_cTET"] = ad_mo_scumi.obs["smooth_cTET"].astype(float)
out3 = OUT_DIR / "zmanseq_momac_metacells_annot_clean.h5ad"
ad_mo_scumi.write_h5ad(out3)
print(f"  Saved {out3.name}")

# Object 4: mo only, mean MC UMI + smooth_cTET
ad_mo_mcumi = build_anndata(mo_madata, mo_cum_pct_df, mc_umi_target, "mo-only / MC UMI")
print("  Computing smooth_cTET (MC UMI) ...")
ad_mo_mcumi.obs["smooth_cTET"] = compute_smooth_ctet(ad_mo_mcumi.copy(), mo_cum_pct_df)
ad_mo_mcumi.obs["smooth_cTET"] = ad_mo_mcumi.obs["smooth_cTET"].astype(float)
out4 = OUT_DIR / "zmanseq_momac_metacells_annot_clean_mcumi.h5ad"
ad_mo_mcumi.write_h5ad(out4)
print(f"  Saved {out4.name}")


# ══════════════════════════════════════════════════════════════════════════════
# 8 — Summary plots (cTET barplots for all 4 objects)
# ══════════════════════════════════════════════════════════════════════════════

print("\nSaving summary plots ...")

def plot_ctet_bar(ad, title, ax):
    df = ad.obs[["cTET", "metacell_cell_type"]].copy()
    df = df.sort_values("cTET").reset_index()
    df["metacell_name"] = df["index"].astype(str)
    df["metacell_cell_type"] = df["metacell_cell_type"].astype(str)
    pal = {k: v for k, v in MYELOID_CELL_TYPE_COLORS.items() if k in df["metacell_cell_type"].values}
    sns.barplot(x=df["cTET"], y=df["metacell_name"], hue=df["metacell_cell_type"],
                width=1, linewidth=0, palette=pal, order=df["metacell_name"].tolist(), ax=ax)
    ax.set_title(title)
    ax.set_xlabel("cTET")
    ax.legend(loc="lower right", fontsize=6)

fig, axes = plt.subplots(1, 4, figsize=(22, 12), sharey=False)
for ax, (ad, title) in zip(axes, [
    (ad_all_scumi, "All myeloid\nsc UMI"),
    (ad_all_mcumi, "All myeloid\nMC UMI"),
    (ad_mo_scumi,  "Mo only\nsc UMI"),
    (ad_mo_mcumi,  "Mo only\nMC UMI"),
]):
    plot_ctet_bar(ad, title, ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "06_ctet_barplots.pdf", bbox_inches="tight")
plt.close()
print("  Saved 06_ctet_barplots.pdf")

print("\nDone.")
print(f"  {out1.name}: {ad_all_scumi.shape}")
print(f"  {out2.name}: {ad_all_mcumi.shape}")
print(f"  {out3.name}: {ad_mo_scumi.shape}")
print(f"  {out4.name}: {ad_mo_mcumi.shape}")
