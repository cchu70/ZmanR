#!/usr/bin/env python
"""
Build annotated mouse CosMx h5ad files from scratch:
  - Expression matrix from spatial_transcriptomics.zip
  - Spatial coordinates from coordinates.tsv (if present)
  - Cell-type annotations from mmc7 Table S6A

Output: cosmx_mhp.h5ad, cosmx_pdx.h5ad in results/10/
Annotation columns in obs:
  mouse_model_type, malignant_flag, tumor_cell_state, TSPS, SWED,
  GPM_prob, MTC_prob, NEU_prob, PPR_prob, entropy, std_entropy
"""

import zipfile, subprocess, tempfile
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.io import mmread
from pathlib import Path

DOCS_DIR   = Path("/home/unix/cchu/projects/ZmanR/pqe/docs")
RESULTS    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/10")
ZIP_PATH   = DOCS_DIR / "1-s2.0-S1535610825003666-mmc7.S6A.csv.zip"
CSV_NAME   = "1-s2.0-S1535610825003666-mmc7.S6A.csv"
EXPR_ZIP   = Path("/mnt/thechenlab/ClaudiaC/gbm_migliozzi/spatial_transcriptomics.zip")

# zip folder name → annotation model label
SAMPLE_MAP = {"mHP": "mhp", "PDX": "pdx"}

# ── 1. Load annotation table ──────────────────────────────────────────────────
print("Loading annotation CSV ...")
with zipfile.ZipFile(ZIP_PATH) as zf:
    with zf.open(CSV_NAME) as f:
        ann = pd.read_csv(f, skiprows=1, dtype=str)

# Standardise column names
ann.columns = [c.strip().replace(" ", "_") for c in ann.columns]
print(f"  {len(ann):,} rows, columns: {list(ann.columns)}")

# Convert numeric columns
num_cols = ["TSPS", "GPM_prob", "MTC_prob", "NEU_prob", "PPR_prob",
            "entropy", "std_entropy"]
for c in num_cols:
    if c in ann.columns:
        ann[c] = pd.to_numeric(ann[c], errors="coerce")

# annotation columns to keep
ann_cols = ["mouse_model_type", "malignant_flag", "tumor_cell_state", "SWED"] + num_cols

# ── 2. Build each mouse h5ad from scratch ─────────────────────────────────────
for sname, model_label in SAMPLE_MAP.items():
    print(f"\n── {sname} ──")
    h5ad_path = RESULTS / f"cosmx_{sname.lower()}.h5ad"

    ann_sub = ann[ann["mouse_model_type"] == model_label].set_index("Cell_ID")
    ann_sub = ann_sub[[c for c in ann_cols if c in ann_sub.columns]]
    print(f"  Annotation rows: {len(ann_sub):,}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # Extract inner zip
        subprocess.run(["unzip", "-o", str(EXPR_ZIP),
                        f"data_deposit_cosmx_paper/{sname}.zip", "-d", str(tmp)],
                       check=True, capture_output=True)
        inner_zip = tmp / "data_deposit_cosmx_paper" / f"{sname}.zip"
        subprocess.run(["unzip", "-o", str(inner_zip), "-d", str(tmp)],
                       check=True, capture_output=True)

        mat_f    = next(tmp.rglob("matrix.mtx"))
        cells_f  = next(tmp.rglob("cells.tsv"))
        genes_f  = next(tmp.rglob("genes.tsv"))
        coord_f  = list(tmp.rglob("coordinates.tsv"))

        cells = pd.read_csv(cells_f, header=None, names=["cell_id"])["cell_id"].tolist()
        genes = pd.read_csv(genes_f, header=None, names=["gene"])["gene"].tolist()
        X     = sp.csr_matrix(mmread(mat_f))  # cells x genes
        print(f"  Expression: {X.shape[0]:,} cells x {X.shape[1]:,} genes")

        coords_df = None
        if coord_f:
            coords_df = pd.read_csv(coord_f[0], sep="\t").set_index("cell_id")

    # obs: annotations aligned to cell order
    obs = ann_sub.reindex(cells).reset_index(drop=False).rename(columns={"Cell_ID": "cell_id"})
    obs.index = cells
    # fill NaN strings
    for col in obs.select_dtypes(include="object").columns:
        obs[col] = obs[col].fillna("")

    var = pd.DataFrame(index=genes)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = cells

    if coords_df is not None:
        xy = coords_df.reindex(cells)[["x", "y"]].values.astype(float)
        adata.obsm["spatial"] = xy

    adata.write_h5ad(h5ad_path)
    print(f"  Saved {h5ad_path.name}")
    print(f"  tumor_cell_state value counts:\n{adata.obs['tumor_cell_state'].value_counts().to_string()}")

print("\nDone.")
