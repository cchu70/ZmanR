"""Run Palantir pseudotime on mcumi metacells using hvg gene set."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, scanpy as sc, palantir

ADATA_PATH = "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot_clean_mcumi.h5ad"
OUT_PLOT   = "/home/unix/cchu/projects/ZmanR/pqe/results/pseudotime_mcumi_hvg/palantir_pseudotime.png"

adata_full = sc.read_h5ad(ADATA_PATH)
genes = list(adata_full.var_names[adata_full.var["highly_variable"]])
print(f"HVG genes: {len(genes)}")
del adata_full

adata = sc.read_h5ad(ADATA_PATH)
genes = [g for g in genes if g in adata.var_names]
print(f"Genes present in adata: {len(genes)}")
adata = adata[:, genes].copy()

# Mark all genes as highly_variable so run_pca uses the full gene set
adata.var["highly_variable"] = True

n_pcs = min(10, len(genes) - 1)
palantir.utils.run_pca(adata, n_components=n_pcs)
palantir.utils.run_diffusion_maps(adata, n_components=min(5, n_pcs))
palantir.utils.determine_multiscale_space(adata)

start_cell = adata.obs["cTET"].idxmin()
print(f"Start cell: {start_cell} (cTET = {adata.obs.loc[start_cell, 'cTET']:.4f})")

pr_res = palantir.core.run_palantir(adata, start_cell, num_waypoints=40, knn=10, use_early_cell_as_start=True)
adata.obs["palantir_pseudotime"] = pr_res.pseudotime

sc_x, sc_y, pt = adata.obs["sc_x"].astype(float), adata.obs["sc_y"].astype(float), adata.obs["palantir_pseudotime"].astype(float)
fig, ax = plt.subplots(figsize=(6, 5))
sc_plot = ax.scatter(sc_x, sc_y, c=pt, cmap="plasma", s=8, alpha=0.8, linewidths=0)
plt.colorbar(sc_plot, ax=ax, label="Pseudotime")
ax.set_title("Palantir pseudotime (mcumi hvg)"); ax.set_xlabel("sc_x"); ax.set_ylabel("sc_y")
plt.tight_layout(); plt.savefig(OUT_PLOT, dpi=150); print(f"Plot saved to {OUT_PLOT}")

out_tsv = OUT_PLOT.replace(".png", ".tsv")
adata.obs.to_csv(out_tsv, sep="\t"); print(f"TSV saved to {out_tsv}")
