"""
Run metacell construction on human CosMx 1k MoMac data (UM01).
target_metacell_size=15 (finer resolution than s30 run; s10 was below the
algorithm's viable minimum for this dataset size, giving ~56 metacells).

Input : results/10/momac_um01_adata.h5ad  (3,042 cells × 914 genes)
Output: results/12/s10/
  momac_um01_cells_s10.h5ad      — cells adata with metacell assignments
  momac_um01_metacells_s10.h5ad  — metacell adata with annotations
  12_s10_*.pdf                   — summary plots
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
from pathlib import Path

ADATA_PATH = Path("/home/unix/cchu/projects/ZmanR/pqe/results/10/momac_um01_adata.h5ad")
OUT_DIR    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/12/s10")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TAG = "s10"   # suffix for filenames

# ── 1. Load ────────────────────────────────────────────────────────────────────
print("Loading adata ...")
adata = sc.read_h5ad(ADATA_PATH)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

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
mc.ut.set_name(adata, f"momac_um01_{TAG}")  # TAG=s10 kept for file naming continuity
mc.pl.mark_lateral_genes(adata, lateral_gene_names=[])
mc.pl.mark_noisy_genes(adata, noisy_gene_names=[])

# ── 4. Run metacell construction ───────────────────────────────────────────────
# ~3k cells; target=15 → ~56 metacells (target=10 is below algorithm minimum)
print("Running divide-and-conquer metacells (target_metacell_size=15) ...")
mc.pl.compute_divide_and_conquer_metacells(
    adata,
    target_metacell_size=15,
    select_min_gene_relative_variance=0.1,
    select_min_gene_total=50,
    select_downsample_min_samples=500,
    random_seed=42,
)

# ── 5. Collect metacells ───────────────────────────────────────────────────────
mcdata = mc.pl.collect_metacells(adata, name=f"momac_um01_{TAG}.metacells", random_seed=42)
print(f"  Metacells: {mcdata.n_obs} × {mcdata.n_vars} genes")

# ── 6. Annotate metacell object ────────────────────────────────────────────────
cells = adata.obs.copy()

rows = []
for mc_idx in range(mcdata.n_obs):
    mask = cells["metacell"] == mc_idx
    grp  = cells[mask]
    row  = {"metacell_idx": mc_idx, "n_cells": mask.sum()}

    for col in ["cellType", "cellTypeWithMyeloidSubsets"]:
        if col in grp.columns:
            vc = grp[col].value_counts(normalize=True)
            for cat, prop in vc.items():
                row[f"{col}__{cat}"] = prop

    for col in ["FOV_ID", "SWED"]:
        if col in grp.columns:
            vc = grp[col].value_counts(normalize=True)
            row[f"{col}__dominant"] = vc.index[0] if len(vc) > 0 else ""
            for cat, prop in vc.items():
                row[f"{col}__{cat}"] = prop

    for col in ["TSPS"]:
        if col in grp.columns:
            row[f"{col}_mean"] = grp[col].astype(float).mean()

    rows.append(row)

mc_ann = pd.DataFrame(rows).set_index("metacell_idx").fillna(0)
for col in mc_ann.columns:
    mcdata.obs[col] = mc_ann[col].values

# kNN adjacency + UMAP on log-normalized expression
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
adata.write_h5ad(OUT_DIR / f"momac_um01_cells_{TAG}.h5ad")
mcdata.write_h5ad(OUT_DIR / f"momac_um01_metacells_{TAG}.h5ad")
print("  Saved h5ad files")

# ── 8. Summary plots ───────────────────────────────────────────────────────────
PALETTE = plt.rcParams["axes.prop_cycle"].by_key()["color"]

def stacked_bar(df, title, xlabel, ylabel, out_path, fontsize=7, ncol=1, cmap=None):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.08), 4))
    bottom = np.zeros(len(df))
    for i, col in enumerate(df.columns):
        color = cmap(i) if cmap else PALETTE[i % len(PALETTE)]
        ax.bar(range(len(df)), df[col].values, bottom=bottom,
               label=col, color=color, width=1.0, linewidth=0)
        bottom += df[col].values
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=fontsize, ncol=ncol)
    ax.set_xlim(-0.5, len(df) - 0.5)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)

# 8a. Cells per metacell
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(mcdata.obs["n_cells"], bins=30, edgecolor="white")
ax.set_xlabel("Cells per metacell")
ax.set_ylabel("Count")
ax.set_title(f"Cells per metacell distribution (target={TAG})")
fig.tight_layout()
fig.savefig(OUT_DIR / f"12_{TAG}_cells_per_metacell.pdf")
plt.close(fig)

# 8b. cellType composition
ct_cols = [c for c in mcdata.obs.columns if c.startswith("cellType__") and "WithMyeloid" not in c]
if ct_cols:
    ct_df = mcdata.obs[ct_cols].copy()
    ct_df.columns = [c.replace("cellType__", "") for c in ct_cols]
    ct_df = ct_df.sort_values(ct_df.columns[0], ascending=False)
    stacked_bar(ct_df, "Cell type composition per metacell",
                "Metacell (sorted by first cell type)", "Proportion",
                OUT_DIR / f"12_{TAG}_celltype_composition.pdf")

# 8c. Myeloid subset composition
sub_cols = [c for c in mcdata.obs.columns if c.startswith("cellTypeWithMyeloidSubsets__") and "dominant" not in c]
if sub_cols:
    sub_df = mcdata.obs[sub_cols].copy()
    sub_df.columns = [c.replace("cellTypeWithMyeloidSubsets__", "") for c in sub_cols]
    sub_df = sub_df.sort_values(sub_df.columns[0], ascending=False)
    stacked_bar(sub_df, "Myeloid subset composition per metacell",
                "Metacell (sorted by first subset)", "Proportion",
                OUT_DIR / f"12_{TAG}_myeloid_subset_composition.pdf")

# 8d. SWED distribution
swed_cols = [c for c in mcdata.obs.columns if c.startswith("SWED__") and "dominant" not in c]
if swed_cols:
    swed_df = mcdata.obs[swed_cols].copy()
    swed_df.columns = [c.replace("SWED__", "") for c in swed_cols]
    swed_df = swed_df.loc[swed_df.idxmax(axis=1).sort_values().index]
    stacked_bar(swed_df, "SWED distribution per metacell",
                "Metacell (sorted by dominant SWED)", "Proportion",
                OUT_DIR / f"12_{TAG}_swed_distribution.pdf")

# 8e. FOV distribution
fov_cols = [c for c in mcdata.obs.columns if c.startswith("FOV_ID__") and "dominant" not in c]
if fov_cols:
    fov_df = mcdata.obs[fov_cols].copy()
    fov_df.columns = [c.replace("FOV_ID__", "") for c in fov_cols]
    fov_df = fov_df.loc[fov_df.idxmax(axis=1).sort_values().index]
    n_fovs = len(fov_df.columns)
    cmap = matplotlib.colormaps.get_cmap("tab20").resampled(n_fovs)
    cols_sorted = sorted(fov_df.columns, key=lambda x: int(x) if x.isdigit() else 999)
    fov_df = fov_df[cols_sorted]
    stacked_bar(fov_df, "FOV distribution per metacell",
                "Metacell (sorted by dominant FOV)", "Proportion",
                OUT_DIR / f"12_{TAG}_fov_distribution.pdf",
                fontsize=6, ncol=2, cmap=cmap)

# 8f. TSPS violin per dominant cell type
if "TSPS_mean" in mcdata.obs.columns and ct_cols:
    ct_df2 = mcdata.obs[[c for c in ct_cols]].copy()
    ct_df2.columns = [c.replace("cellType__", "") for c in ct_cols]
    mcdata.obs["dominant_cellType"] = ct_df2.idxmax(axis=1)
    fig, ax = plt.subplots(figsize=(5, 4))
    groups = [mcdata.obs.loc[mcdata.obs["dominant_cellType"] == ct, "TSPS_mean"].dropna().values
              for ct in ct_df2.columns]
    vp = ax.violinplot(groups, positions=range(len(ct_df2.columns)), showmedians=True)
    ax.set_xticks(range(len(ct_df2.columns)))
    ax.set_xticklabels(ct_df2.columns, rotation=30, ha="right")
    ax.set_ylabel("Mean TSPS")
    ax.set_title(f"TSPS per dominant cell type ({TAG})")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"12_{TAG}_tsps_by_celltype.pdf")
    plt.close(fig)

# 8g. UMAP colored by dominant cell type
if "X_umap" in mcdata.obsm and ct_cols:
    ct_df3 = mcdata.obs[[c for c in ct_cols]].copy()
    ct_df3.columns = [c.replace("cellType__", "") for c in ct_cols]
    mcdata.obs["dominant_cellType"] = ct_df3.idxmax(axis=1)
    fig, ax = plt.subplots(figsize=(5, 5))
    for i, cat in enumerate(ct_df3.columns):
        mask = mcdata.obs["dominant_cellType"] == cat
        xy = mcdata.obsm["X_umap"][mask]
        ax.scatter(xy[:, 0], xy[:, 1], s=15, label=cat,
                   color=PALETTE[i % len(PALETTE)], alpha=0.8)
    ax.set_title(f"Metacell UMAP — dominant cell type ({TAG})")
    ax.legend(fontsize=8)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"12_{TAG}_umap_celltype.pdf")
    plt.close(fig)

print(f"\nAll outputs saved to {OUT_DIR}")
print("Done.")
