"""
Save metacell h5ad data as RDS via Python+rpy2, bypassing basilisk.
Writes: zmanseq_momac_metacells_annot.rds (used by R pseudotime scripts)
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
numpy2ri.activate(); pandas2ri.activate()

H5AD = "/home/unix/cchu/projects/ZmanR/pqe/results/06/zmanseq_momac_metacells_annot.h5ad"
OUT  = H5AD.replace(".h5ad", ".rds")

print("Reading", H5AD)
adata = sc.read_h5ad(H5AD)
print(f"Shape: {adata.shape}")

hvg_genes = list(adata.var_names[adata.var["highly_variable"]])
print(f"HVG genes: {len(hvg_genes)}")

# Convert sparse X to dense for rpy2 (counts are small: 54 x 27k)
X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
# X is cells x genes; R expects genes x cells
X_t = X.T

base    = importr("base")
Matrix  = importr("Matrix")

r_mat = ro.r.matrix(ro.FloatVector(X_t.flatten(order='F')),
                    nrow=X_t.shape[0], ncol=X_t.shape[1])
r_mat.rownames = ro.StrVector(list(adata.var_names))
r_mat.colnames = ro.StrVector(list(adata.obs_names))

r_sparse = Matrix.Matrix(r_mat, sparse=True)
r_meta   = pandas2ri.py2rpy(adata.obs)
r_genes  = ro.StrVector(list(adata.var_names))
r_hvg    = ro.StrVector(hvg_genes)

r_list = ro.ListVector({
    "counts":     r_sparse,
    "metadata":   r_meta,
    "gene_names": r_genes,
    "hvg_genes":  r_hvg,
})
base.saveRDS(r_list, file=OUT)
print(f"Saved to {OUT}")
