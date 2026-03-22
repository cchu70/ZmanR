"""
Run metacell construction on human CosMx 1k MoMac data (UM01).

Input : results/10/momac_um01_adata.h5ad  (3,042 cells × 914 genes)
Output: results/12/
  momac_um01_cells.h5ad        — cells adata with metacell assignments
  momac_um01_metacells.h5ad    — metacell adata with annotations
  12_*.pdf                     — summary plots
"""

import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import metacells as mc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

ADATA_PATH = Path("/home/unix/cchu/projects/ZmanR/pqe/results/10/momac_um01_adata.h5ad")
OUT_DIR    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/12")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. Load ────────────────────────────────────────────────────────────────────
print("Loading adata ...")
adata = sc.read_h5ad(ADATA_PATH)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Ensure integer counts in X
adata.X = adata.X.astype(np.float32)

# ── 2. Gene filtering (human panel) ───────────────────────────────────────────
pattern = re.compile(
    r"^RPL|^RPS|^MT-|^IGH|^IGK|^IGL|^HIST|^SNOR|^MIR",
    re.IGNORECASE,
)
genes_keep = [g for g in adata.var_names if not pattern.search(g)]
adata = adata[:, genes_keep].copy()
print(f"  After gene filter: {adata.n_vars} genes")

# ── 3. Set up metacells ────────────────────────────────────────────────────────
mc.ut.set_name(adata, "momac_um01")
mc.pl.mark_lateral_genes(adata, lateral_gene_names=[])
mc.pl.mark_noisy_genes(adata, noisy_gene_names=[])

# ── 4. Run metacell construction ───────────────────────────────────────────────
# ~3k cells; target=30 → ~100 metacells
print("Running divide-and-conquer metacells ...")
mc.pl.compute_divide_and_conquer_metacells(
    adata,
    target_metacell_size=30,
    select_min_gene_relative_variance=0.1,
    select_min_gene_total=50,
    select_downsample_min_samples=500,
    random_seed=42,
)

# ── 5. Collect metacells ───────────────────────────────────────────────────────
mcdata = mc.pl.collect_metacells(adata, name="momac_um01.metacells", random_seed=42)
print(f"  Metacells: {mcdata.n_obs} × {mcdata.n_vars} genes")

# ── 6. Annotate metacell object ────────────────────────────────────────────────
cells = adata.obs.copy()
cells["metacell"] = adata.obs["metacell"]   # integer index from mc

ann_cols = ["cellType", "cellTypeWithMyeloidSubsets", "FOV_ID", "SWED", "TSPS"]

rows = []
for mc_idx in range(mcdata.n_obs):
    mask = cells["metacell"] == mc_idx
    grp  = cells[mask]
    row  = {"metacell_idx": mc_idx, "n_cells": mask.sum()}

    # Proportion of each cellType
    for col in ["cellType", "cellTypeWithMyeloidSubsets"]:
        if col in grp.columns:
            vc = grp[col].value_counts(normalize=True)
            for cat, prop in vc.items():
                row[f"{col}__{cat}"] = prop

    # Dominant / distribution of FOV_ID and SWED
    for col in ["FOV_ID", "SWED"]:
        if col in grp.columns:
            vc = grp[col].value_counts(normalize=True)
            row[f"{col}__dominant"] = vc.index[0] if len(vc) > 0 else ""
            for cat, prop in vc.items():
                row[f"{col}__{cat}"] = prop

    # Numeric columns: mean
    for col in ["TSPS"]:
        if col in grp.columns:
            row[f"{col}_mean"] = grp[col].astype(float).mean()

    rows.append(row)

mc_ann = pd.DataFrame(rows).set_index("metacell_idx")
mc_ann = mc_ann.fillna(0)

# Add to mcdata.obs
for col in mc_ann.columns:
    mcdata.obs[col] = mc_ann[col].values

# Add adjacency matrix via kNN on log-normalized expression
_tmp = mcdata.copy()
sc.pp.normalize_total(_tmp)
sc.pp.log1p(_tmp)
n_comps = min(30, _tmp.n_obs - 2, _tmp.n_vars - 1)
sc.pp.pca(_tmp, n_comps=n_comps)
n_neighbors = min(15, _tmp.n_obs - 1)
sc.pp.neighbors(_tmp, n_neighbors=n_neighbors)
sc.tl.umap(_tmp)
mcdata.obsp["connectivities"] = _tmp.obsp["connectivities"]
mcdata.obsp["distances"]      = _tmp.obsp["distances"]
mcdata.uns["neighbors"]       = _tmp.uns["neighbors"]
mcdata.obsm["X_umap"]         = _tmp.obsm["X_umap"]

# ── 7. Save h5ad ──────────────────────────────────────────────────────────────
adata.write_h5ad(OUT_DIR / "momac_um01_cells.h5ad")
mcdata.write_h5ad(OUT_DIR / "momac_um01_metacells.h5ad")
print("  Saved h5ad files")

# ── 8. Summary plots ───────────────────────────────────────────────────────────
PALETTE = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# 8a. Cells per metacell
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(mcdata.obs["n_cells"], bins=20, edgecolor="white")
ax.set_xlabel("Cells per metacell")
ax.set_ylabel("Count")
ax.set_title("Cells per metacell distribution")
fig.tight_layout()
fig.savefig(OUT_DIR / "12_cells_per_metacell.pdf")
plt.close(fig)

# 8b. cellType composition stacked bar (sorted by Macrophage proportion)
ct_cols = [c for c in mcdata.obs.columns if c.startswith("cellType__") and "WithMyeloid" not in c]
if ct_cols:
    ct_df = mcdata.obs[ct_cols].copy()
    ct_df.columns = [c.replace("cellType__", "") for c in ct_cols]
    ct_df = ct_df.sort_values(ct_df.columns[0], ascending=False)
    fig, ax = plt.subplots(figsize=(max(8, len(ct_df) * 0.12), 4))
    bottom = np.zeros(len(ct_df))
    for i, col in enumerate(ct_df.columns):
        ax.bar(range(len(ct_df)), ct_df[col].values, bottom=bottom,
               label=col, color=PALETTE[i % len(PALETTE)], width=1.0, linewidth=0)
        bottom += ct_df[col].values
    ax.set_xlabel("Metacell (sorted by first cell type)")
    ax.set_ylabel("Proportion")
    ax.set_title("Cell type composition per metacell")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_xlim(-0.5, len(ct_df) - 0.5)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_celltype_composition.pdf")
    plt.close(fig)

# 8c. cellTypeWithMyeloidSubsets composition stacked bar
sub_cols = [c for c in mcdata.obs.columns if c.startswith("cellTypeWithMyeloidSubsets__") and "dominant" not in c]
if sub_cols:
    sub_df = mcdata.obs[sub_cols].copy()
    sub_df.columns = [c.replace("cellTypeWithMyeloidSubsets__", "") for c in sub_cols]
    sub_df = sub_df.sort_values(sub_df.columns[0], ascending=False)
    fig, ax = plt.subplots(figsize=(max(8, len(sub_df) * 0.12), 4))
    bottom = np.zeros(len(sub_df))
    for i, col in enumerate(sub_df.columns):
        ax.bar(range(len(sub_df)), sub_df[col].values, bottom=bottom,
               label=col, color=PALETTE[i % len(PALETTE)], width=1.0, linewidth=0)
        bottom += sub_df[col].values
    ax.set_xlabel("Metacell (sorted by first subset)")
    ax.set_ylabel("Proportion")
    ax.set_title("Myeloid subset composition per metacell")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    ax.set_xlim(-0.5, len(sub_df) - 0.5)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_myeloid_subset_composition.pdf")
    plt.close(fig)

# 8d. SWED distribution stacked bar per metacell
swed_cols = [c for c in mcdata.obs.columns if c.startswith("SWED__") and "dominant" not in c]
if swed_cols:
    swed_df = mcdata.obs[swed_cols].copy()
    swed_df.columns = [c.replace("SWED__", "") for c in swed_cols]
    # Sort metacells by dominant SWED
    swed_df = swed_df.loc[swed_df.idxmax(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(max(8, len(swed_df) * 0.12), 4))
    bottom = np.zeros(len(swed_df))
    for i, col in enumerate(sorted(swed_df.columns)):
        ax.bar(range(len(swed_df)), swed_df[col].values, bottom=bottom,
               label=col, color=PALETTE[i % len(PALETTE)], width=1.0, linewidth=0)
        bottom += swed_df[col].values
    ax.set_xlabel("Metacell (sorted by dominant SWED)")
    ax.set_ylabel("Proportion")
    ax.set_title("SWED distribution per metacell")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_xlim(-0.5, len(swed_df) - 0.5)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_swed_distribution.pdf")
    plt.close(fig)

# 8e. FOV distribution stacked bar per metacell
fov_cols = [c for c in mcdata.obs.columns if c.startswith("FOV_ID__") and "dominant" not in c]
if fov_cols:
    fov_df = mcdata.obs[fov_cols].copy()
    fov_df.columns = [c.replace("FOV_ID__", "") for c in fov_cols]
    fov_df = fov_df.loc[fov_df.idxmax(axis=1).sort_values().index]
    n_fovs = len(fov_df.columns)
    cmap = plt.cm.get_cmap("tab20", n_fovs)
    fig, ax = plt.subplots(figsize=(max(8, len(fov_df) * 0.12), 4))
    bottom = np.zeros(len(fov_df))
    for i, col in enumerate(sorted(fov_df.columns, key=lambda x: int(x) if x.isdigit() else 999)):
        ax.bar(range(len(fov_df)), fov_df[col].values, bottom=bottom,
               label=col, color=cmap(i), width=1.0, linewidth=0)
        bottom += fov_df[col].values
    ax.set_xlabel("Metacell (sorted by dominant FOV)")
    ax.set_ylabel("Proportion")
    ax.set_title("FOV distribution per metacell")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    ax.set_xlim(-0.5, len(fov_df) - 0.5)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_fov_distribution.pdf")
    plt.close(fig)

# 8f. TSPS violin per cellType (dominant cell type label per metacell)
if "TSPS_mean" in mcdata.obs.columns and ct_cols:
    ct_df2 = mcdata.obs[[c for c in ct_cols]].copy()
    ct_df2.columns = [c.replace("cellType__", "") for c in ct_cols]
    mcdata.obs["dominant_cellType"] = ct_df2.idxmax(axis=1)
    fig, ax = plt.subplots(figsize=(5, 4))
    groups = [mcdata.obs.loc[mcdata.obs["dominant_cellType"] == ct, "TSPS_mean"].dropna().values
              for ct in ct_df2.columns]
    ax.violinplot(groups, positions=range(len(ct_df2.columns)), showmedians=True)
    ax.set_xticks(range(len(ct_df2.columns)))
    ax.set_xticklabels(ct_df2.columns, rotation=30, ha="right")
    ax.set_ylabel("Mean TSPS")
    ax.set_title("TSPS per dominant cell type")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_tsps_by_celltype.pdf")
    plt.close(fig)

# 8g. UMAP colored by dominant cell type
if "X_umap" in mcdata.obsm and ct_cols:
    ct_df3 = mcdata.obs[[c for c in ct_cols]].copy()
    ct_df3.columns = [c.replace("cellType__", "") for c in ct_cols]
    mcdata.obs["dominant_cellType"] = ct_df3.idxmax(axis=1)
    fig, ax = plt.subplots(figsize=(5, 5))
    cats = mcdata.obs["dominant_cellType"].unique()
    for i, cat in enumerate(cats):
        mask = mcdata.obs["dominant_cellType"] == cat
        xy = mcdata.obsm["X_umap"][mask]
        ax.scatter(xy[:, 0], xy[:, 1], s=30, label=cat,
                   color=PALETTE[i % len(PALETTE)], alpha=0.8)
    ax.set_title("Metacell UMAP — dominant cell type")
    ax.legend(fontsize=8)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_umap_celltype.pdf")
    plt.close(fig)

print(f"\nAll outputs saved to {OUT_DIR}")
print("Done.")
