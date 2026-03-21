"""
Run metacell construction on Zman-seq data, restricted to myeloid cell types
from IgG and aTrem2 treatments.

Cell filtering (in addition to QC):
  - Treatment in ["IgG", "aTrem2"]
  - cluster_colors in myeloid cell type list

Gene filtering (applied before calling metacells):
  1. Remove mitochondrial, ribosomal, Ig, poorly-supported model genes (regex, from notebook)
  2. Remove unwanted RNA biotypes via Ensembl biomart (from notebook)

Additional gene selection handled internally by metacells:
  - select_min_gene_relative_variance = 0.1  (T_vm threshold)
  - select_min_gene_total              = 100  (minimum total UMI per gene)

Metacell construction parameters:
  - target_metacell_size         = 100
  - select_downsample_min_samples = 750
"""

import re
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import metacells as mc
from pathlib import Path

ADATA_PATH = "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_DIR    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/06")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TREATMENTS = ["IgG", "aTrem2"]

MYELOID_CELLTYPES = [
    "Monocyte", "MoMac1", "Monocytes", "Gpnmb_TAM", "Arg1_TAM",
    "MigDC", "Acp5_TAM", "cDC1", "cDC2", "MoMac2", "MoMac", "pDC", "MonDC",
]


# ── 1. Load raw data ──────────────────────────────────────────────────────────

print("Loading h5ad ...")
adata = sc.read_h5ad(ADATA_PATH)
print(f"Raw: {adata.n_obs} cells × {adata.n_vars} genes")


# ── 2. Gene filtering: biotypes (from notebook) ───────────────────────────────

unwanted_set = set()
try:
    from pybiomart import Server
    server = Server(host="http://www.ensembl.org")
    mart   = server["ENSEMBL_MART_ENSEMBL"]["mmusculus_gene_ensembl"]
    unwanted_biotypes = ["rRNA", "rRNA_pseudogene", "Mt_rRNA", "snRNA", "snoRNA", "misc_RNA"]
    unwanted_genes    = mart.query(
        attributes=["mgi_symbol", "gene_biotype"],
        filters={"biotype": unwanted_biotypes},
    )
    unwanted_set = set(unwanted_genes["MGI symbol"].dropna())
    print(f"Biomart: {len(unwanted_set)} unwanted-biotype genes")
except Exception as e:
    print(f"Biomart query failed ({e}), skipping biotype filter")


# ── 3. Gene filtering: regex (from notebook) ──────────────────────────────────

pattern = re.compile(
    r"^Rpl|^Rps|^mt-|^Igh|^Igk|^Igl|Rik$|^Gm\d|^AK\d|^AY\d|^BC\d|^BG\d|^BM\d|^LOC|^EG\d|^ENSMUSG",
    re.IGNORECASE,
)
genes_keep = [
    g for g in adata.var_names
    if not pattern.search(g) and g not in unwanted_set
]
adata = adata[:, genes_keep].copy()
print(f"After regex/biotype filter: {adata.n_vars} genes")


# ── 4. Cell filtering: total UMI > 300 (from notebook) ───────────────────────

sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=False)
adata = adata[adata.obs["total_counts"] > 300].copy()
print(f"After cell UMI filter (>300): {adata.n_obs} cells")


# ── 5. Cell filtering: treatment and cell type ────────────────────────────────

treatment_mask = adata.obs["Treatment"].isin(TREATMENTS)
celltype_mask  = adata.obs["cluster_colors"].isin(MYELOID_CELLTYPES)
adata = adata[treatment_mask & celltype_mask].copy()
print(f"After treatment ({TREATMENTS}) + cell type filter: {adata.n_obs} cells")
print("Treatment counts:\n", adata.obs["Treatment"].value_counts().to_string())
print("Cell type counts:\n", adata.obs["cluster_colors"].value_counts().to_string())


# ── 6. Set up metacells ───────────────────────────────────────────────────────

mc.ut.set_name(adata, "zmanseq_myeloid")

# Initialize required boolean masks
mc.pl.mark_lateral_genes(adata, lateral_gene_names=[])
mc.pl.mark_noisy_genes(adata, noisy_gene_names=[])


# ── 7. Run metacell construction ──────────────────────────────────────────────

print("Running divide-and-conquer metacells ...")
mc.pl.compute_divide_and_conquer_metacells(
    adata,
    target_metacell_size=100,
    select_min_gene_relative_variance=0.1,
    select_min_gene_total=100,
    select_downsample_min_samples=750,
    random_seed=42,
)


# ── 8. Collect and save ───────────────────────────────────────────────────────

mcdata = mc.pl.collect_metacells(adata, name="zmanseq_myeloid.metacells", random_seed=42)
print(f"Metacells: {mcdata.n_obs} metacells × {mcdata.n_vars} genes")

cells_out = OUT_DIR / "zmanseq_myeloid_cells.h5ad"
mc_out    = OUT_DIR / "zmanseq_myeloid_metacells.h5ad"

adata.write_h5ad(cells_out)


# ── 9. Metacell adjacency matrix via kNN on log-transformed expression ────────

print("Computing metacell adjacency matrix ...")
_tmp = mcdata.copy()
sc.pp.log1p(_tmp)
n_comps = min(50, _tmp.n_obs - 2, _tmp.n_vars - 1)
sc.pp.pca(_tmp, n_comps=n_comps)
n_neighbors = min(15, _tmp.n_obs - 1)
sc.pp.neighbors(_tmp, n_neighbors=n_neighbors)

mcdata.obsp["connectivities"] = _tmp.obsp["connectivities"]
mcdata.obsp["distances"]      = _tmp.obsp["distances"]
mcdata.uns["neighbors"]       = _tmp.uns["neighbors"]

adj_out = OUT_DIR / "zmanseq_myeloid_metacell_adjacency.npz"
sp.save_npz(adj_out, _tmp.obsp["connectivities"])

mcdata.write_h5ad(mc_out)

print(f"Saved cells      → {cells_out}")
print(f"Saved metacells  → {mc_out}")
print(f"Saved adjacency  → {adj_out}")
print("Done.")
