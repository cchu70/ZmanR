"""
Run metacell construction on Zman-seq data using tanaylab/metacells 0.9.x.

Gene filtering (applied before calling metacells):
  1. Remove mitochondrial, ribosomal, Ig, poorly-supported model genes (regex, from notebook)
  2. Remove unwanted RNA biotypes via Ensembl biomart (from notebook)

Additional gene selection handled internally by metacells:
  - select_min_gene_relative_variance = 0.1  (T_vm threshold)
  - select_min_gene_total              = 100  (minimum total UMI per gene)

Metacell construction parameters:
  - knn_k               = 100   (K nearest neighbors)
  - target_metacell_size = 100  (target cells per metacell)
  - select_downsample_min_samples = 750  (bootstrap downsampling for gene selection)
"""

import re
import numpy as np
import scanpy as sc
import metacells as mc
from pathlib import Path

ADATA_PATH = "/home/unix/cchu/projects/ZmanR/pqe/results/04/zmanseq.h5ad"
OUT_DIR    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/05")
OUT_DIR.mkdir(parents=True, exist_ok=True)


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


# ── 5. Set up metacells ───────────────────────────────────────────────────────

mc.ut.set_name(adata, "zmanseq")

# Initialize required boolean masks (empty — no lateral or noisy genes beyond
# what the regex/biotype filtering already removed)
mc.pl.mark_lateral_genes(adata, lateral_gene_names=[])
mc.pl.mark_noisy_genes(adata, noisy_gene_names=[])


# ── 6. Run metacell construction ──────────────────────────────────────────────
# Parameters:
#   target_metacell_size         = 100  (K: target cells per metacell;
#                                        knn_k is auto-computed as n_obs/target_metacell_size
#                                        per pile — passing knn_k explicitly triggers a bug
#                                        in metacells 0.9.5)
#   select_min_gene_relative_variance = 0.1  (T_vm threshold)
#   select_min_gene_total        = 100  (minimum total UMI per gene)
#   select_downsample_min_samples = 750 (bootstrap downsampling for gene selection)

print("Running divide-and-conquer metacells ...")
mc.pl.compute_divide_and_conquer_metacells(
    adata,
    target_metacell_size=100,
    select_min_gene_relative_variance=0.1,
    select_min_gene_total=100,
    select_downsample_min_samples=750,
    random_seed=42,
)


# ── 7. Collect and save ───────────────────────────────────────────────────────

mcdata = mc.pl.collect_metacells(adata, name="zmanseq.metacells", random_seed=42)
print(f"Metacells: {mcdata.n_obs} metacells × {mcdata.n_vars} genes")

cells_out = OUT_DIR / "zmanseq_cells.h5ad"
mc_out    = OUT_DIR / "zmanseq_metacells.h5ad"

adata.write_h5ad(cells_out)
mcdata.write_h5ad(mc_out)

print(f"Saved cells     → {cells_out}")
print(f"Saved metacells → {mc_out}")
print("Done.")
