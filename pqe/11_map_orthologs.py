#!/usr/bin/env python
"""
Convert ZmanSeq myeloid adata objects (mouse) to human gene names using mousipy,
then subset to genes present in the CosMx 1k panel (cosmx_um01.h5ad).

For objects with layers['total_umis']:
  - X is set to raw UMI counts before mousipy merges counts
  - After translation, raw counts stored in layers['raw']
  - X is then log1p normalized

Inputs:
  results/06/zmanseq_myeloid_cells.h5ad
  results/06/zmanseq_momac_metacells_annot_clean.h5ad
  results/06/zmanseq_momac_metacells_annot_clean_mcumi.h5ad

Outputs (results/11/):
  zmanseq_myeloid_cells_human.h5ad
  zmanseq_momac_metacells_annot_clean_human.h5ad
  zmanseq_momac_metacells_annot_clean_mcumi_human.h5ad
  human_to_mouse_orthologs.tsv
"""

import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import mousipy
from pathlib import Path

RESULTS_06 = Path("/home/unix/cchu/projects/ZmanR/pqe/results/06")
RESULTS_10 = Path("/home/unix/cchu/projects/ZmanR/pqe/results/10")
OUT_DIR    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/11")
OUT_DIR.mkdir(parents=True, exist_ok=True)

INPUTS = [
    ("zmanseq_myeloid_cells.h5ad",                      "zmanseq_myeloid_cells_human.h5ad"),
    ("zmanseq_momac_metacells_annot_clean.h5ad",         "zmanseq_momac_metacells_annot_clean_human.h5ad"),
    ("zmanseq_momac_metacells_annot_clean_mcumi.h5ad",   "zmanseq_momac_metacells_annot_clean_mcumi_human.h5ad"),
]

# ── Load CosMx 1k gene panel ───────────────────────────────────────────────────
print("Loading CosMx 1k human adata ...")
cosmx = ad.read_h5ad(RESULTS_10 / "cosmx_um01.h5ad", backed="r")
cosmx_genes = set(cosmx.var_names)
print(f"  CosMx 1k genes: {len(cosmx_genes):,}\n")


def clean_var(adata):
    """Cast var columns to writable types for h5ad."""
    for col in adata.var.columns:
        if adata.var[col].dtype == object:
            adata.var[col] = adata.var[col].astype(str).replace("nan", "")
        else:
            adata.var[col] = pd.to_numeric(adata.var[col], errors="coerce")
    return adata


mapping_saved = False

for in_name, out_name in INPUTS:
    print(f"── {in_name} ──")
    adata = ad.read_h5ad(RESULTS_06 / in_name)
    print(f"  Input: {adata.shape[0]:,} cells x {adata.n_vars:,} genes")

    # Swap X to raw UMI counts before mousipy merges
    has_raw = "total_umis" in adata.layers
    if has_raw:
        print("  Setting X = layers['total_umis'] (raw UMI counts) ...")
        adata.X = adata.layers["total_umis"].copy()

    print("  Humanizing ...")
    h_adata = mousipy.translate(adata)
    print(f"  Humanized: {h_adata.shape[0]:,} cells x {h_adata.n_vars:,} genes")

    # Save mapping table once (same gene space for all objects)
    if not mapping_saved:
        mapping_df = h_adata.var.copy().reset_index()
        mapping_df.columns = ["human_gene"] + list(mapping_df.columns[1:])
        mapping_df["in_cosmx_1k"] = mapping_df["human_gene"].isin(cosmx_genes)
        mapping_df.to_csv(OUT_DIR / "human_to_mouse_orthologs.tsv", sep="\t", index=False)
        print(f"  Saved mapping table ({len(mapping_df):,} genes)")
        mapping_saved = True

    overlap = [g for g in h_adata.var_names if g in cosmx_genes]
    print(f"  CosMx 1k overlap: {len(overlap):,} / {h_adata.n_vars:,} genes")

    sub = h_adata[:, overlap].copy()

    if has_raw:
        # Store merged raw counts in layers['raw'], then log1p normalize X
        import scipy.sparse as sp
        raw = sub.X.copy()
        sub.layers["raw"] = sp.csr_matrix(raw) if not sp.issparse(raw) else raw
        print("  Log1p normalizing X ...")
        sc.pp.normalize_total(sub)
        sc.pp.log1p(sub)

    sub = clean_var(sub)
    sub.write_h5ad(OUT_DIR / out_name)
    print(f"  Saved {out_name}\n")

print("Done.")
