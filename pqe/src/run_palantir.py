"""
Run Palantir pseudotime on zmanseq data.
Change ADATA_PATH to point to a different h5ad file.
"""

ADATA_PATH = "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_PLOT   = "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime/palantir_pseudotime.png"

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import palantir

# ── 1. Load and preprocess ───────────────────────────────────────────────────
print("Reading h5ad...")
adata = sc.read_h5ad(ADATA_PATH)

# Keep only annotated cells (have sc_x/sc_y coordinates)
has_coords = ~adata.obs["sc_x"].isna()
adata = adata[has_coords].copy()
print(f"Annotated cells: {adata.n_obs}")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()

# ── 2. Run Palantir preprocessing ────────────────────────────────────────────
print("Running PCA and diffusion maps...")
palantir.preprocess.run_pca(adata, n_components=30)
palantir.utils.run_diffusion_maps(adata, n_components=10)

# ── 3. Pick start cell — use a cell from the "Negative" (uninfected) group ──
neg_cells = adata.obs_names[adata.obs["time_assignment"] == "Negative"]
if len(neg_cells) == 0:
    neg_cells = adata.obs_names  # fallback: use all cells
# Among Negative cells, pick the one with median sc_x (approximate "center")
neg_sc_x = adata.obs.loc[neg_cells, "sc_x"].astype(float)
start_cell = neg_cells[np.argmin(np.abs(neg_sc_x - neg_sc_x.median()))]
print(f"Start cell: {start_cell}")

# ── 4. Run Palantir ──────────────────────────────────────────────────────────
print("Running Palantir...")
pr_res = palantir.core.run_palantir(
    adata,
    start_cell,
    num_waypoints=500,
    use_early_cell_as_start=True,
)

adata.obs["palantir_pseudotime"] = pr_res.pseudotime

# ── 5. Plot using sc_x / sc_y ────────────────────────────────────────────────
sc_x = adata.obs["sc_x"].astype(float)
sc_y = adata.obs["sc_y"].astype(float)
pt   = adata.obs["palantir_pseudotime"].astype(float)

fig, ax = plt.subplots(figsize=(6, 5))
sc_plot = ax.scatter(sc_x, sc_y, c=pt, cmap="plasma", s=1.5, alpha=0.7, linewidths=0)
plt.colorbar(sc_plot, ax=ax, label="Pseudotime")
ax.set_title("Palantir pseudotime")
ax.set_xlabel("sc_x")
ax.set_ylabel("sc_y")
plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=150)
print(f"Plot saved to {OUT_PLOT}")
