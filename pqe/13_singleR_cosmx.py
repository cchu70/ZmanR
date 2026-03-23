"""
Assign CosMx MoMac cells (UM01) to ZmanSeq momac metacells using SingleR.

Test  : results/10/momac_um01_adata.h5ad    (3,042 cells × 914 genes, raw counts)
Ref   : results/11/zmanseq_myeloid_cells_human.h5ad
        (single cells, human orthologs; filtered to Treatment=="IgG" and
         metacell_name present in zmanseq_momac_metacells_annot_clean_human.h5ad)

Output: results/13/singleR_cosmx_um01.tsv
  Columns: cell_id, predicted_label, pruned_label, delta_next,
           score_<metacell_name> × N, plus any extra SingleR metadata.
"""

import os
import subprocess
import tempfile
import textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.io
from pathlib import Path

RSCRIPT     = "Rscript"   # system R 4.5.2 — has SingleR 2.12.0
TEST_H5AD   = Path("/home/unix/cchu/projects/ZmanR/pqe/results/10/momac_um01_adata.h5ad")
REF_H5AD    = Path("/home/unix/cchu/projects/ZmanR/pqe/results/11/zmanseq_momac_metacells_annot_clean_human.h5ad")
REF_H5AD_SC = Path("/home/unix/cchu/projects/ZmanR/pqe/results/11/zmanseq_myeloid_cells_human.h5ad")
OUT_DIR     = Path("/home/unix/cchu/projects/ZmanR/pqe/results/13")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV     = OUT_DIR / "singleR_cosmx_um01.tsv"

# ── 1. Load test cells ─────────────────────────────────────────────────────────
print("Loading test adata ...")
test = sc.read_h5ad(TEST_H5AD)
print(f"  Test: {test.shape[0]:,} cells × {test.n_vars:,} genes")

# Normalize + log1p raw counts for SingleR
sc.pp.normalize_total(test)
sc.pp.log1p(test)

# ── 2. Load metacell reference (to get valid metacell names) ───────────────────
print("Loading metacell reference (to get valid labels) ...")
ref_mc = sc.read_h5ad(REF_H5AD)
valid_metacells = set(ref_mc.obs_names[ref_mc.obs["enrichment"] > 0].tolist())
print(f"  Valid metacell labels (enrichment > 0): {len(valid_metacells)}")

# ── 3. Load single-cell reference and filter ───────────────────────────────────
print("Loading single-cell reference ...")
ref_sc = sc.read_h5ad(REF_H5AD_SC)
print(f"  Ref SC before filter: {ref_sc.shape[0]:,} cells × {ref_sc.n_vars:,} genes")
print(f"  Obs columns: {list(ref_sc.obs.columns)}")

# Filter: Treatment == "IgG" and metacell_name in valid metacells
mask_igg = ref_sc.obs["Treatment"] == "IgG"
mask_mc  = ref_sc.obs["metacell_name"].isin(valid_metacells)
ref_sc   = ref_sc[mask_igg & mask_mc].copy()
print(f"  Ref SC after filter (IgG + valid metacell): {ref_sc.shape[0]:,} cells")
print(f"  Unique metacell labels in ref: {ref_sc.obs['metacell_name'].nunique()}")

# Normalize + log1p (ref_sc X is raw counts from mousipy translation)
sc.pp.normalize_total(ref_sc)
sc.pp.log1p(ref_sc)

labels = ref_sc.obs["metacell_name"].tolist()

# ── 4. Align genes ─────────────────────────────────────────────────────────────
shared = list(ref_sc.var_names)   # already subset to CosMx genes in script 11
test   = test[:, shared].copy()
print(f"  Shared genes: {len(shared)}")

# ── 5. Write sparse matrices + metadata to temp files ─────────────────────────
print("Writing sparse matrices to temp files ...")
with tempfile.TemporaryDirectory() as tmpdir:
    tmpdir = Path(tmpdir)

    def to_csc(X):
        return sp.csc_matrix(X) if not isinstance(X, sp.csc_matrix) else X

    # Write test (genes × cells) as sparse MTX
    test_mtx   = str(tmpdir / "test_expr.mtx")
    test_genes = str(tmpdir / "test_genes.txt")
    test_cells = str(tmpdir / "test_cells.txt")
    scipy.io.mmwrite(test_mtx, to_csc(test.X.T))
    Path(test_genes).write_text("\n".join(test.var_names.tolist()))
    Path(test_cells).write_text("\n".join(test.obs_names.tolist()))

    # Write ref (genes × cells) as sparse MTX
    ref_mtx    = str(tmpdir / "ref_expr.mtx")
    ref_genes  = str(tmpdir / "ref_genes.txt")
    ref_cells  = str(tmpdir / "ref_cells.txt")
    labels_txt = str(tmpdir / "labels.txt")
    out_csv    = str(tmpdir / "singler_out.csv")
    scipy.io.mmwrite(ref_mtx, to_csc(ref_sc.X.T))
    Path(ref_genes).write_text("\n".join(ref_sc.var_names.tolist()))
    Path(ref_cells).write_text("\n".join(ref_sc.obs_names.tolist()))
    Path(labels_txt).write_text("\n".join(labels))

    # ── 6. Run SingleR via R subprocess ───────────────────────────────────────
    r_script = textwrap.dedent(f"""
        suppressPackageStartupMessages({{
            library(Matrix)
            library(SingleR)
            library(BiocParallel)
        }})
        cat("Reading sparse matrices ...\\n")
        test_mat <- readMM("{test_mtx}")
        ref_mat  <- readMM("{ref_mtx}")

        test_genes <- readLines("{test_genes}")
        test_cells <- readLines("{test_cells}")
        ref_genes  <- readLines("{ref_genes}")
        ref_cells  <- readLines("{ref_cells}")
        labels     <- readLines("{labels_txt}")

        rownames(test_mat) <- test_genes
        colnames(test_mat) <- test_cells
        rownames(ref_mat)  <- ref_genes
        colnames(ref_mat)  <- ref_cells

        cat(sprintf("  Test: %d genes x %d cells\\n", nrow(test_mat), ncol(test_mat)))
        cat(sprintf("  Ref:  %d genes x %d cells\\n",  nrow(ref_mat),  ncol(ref_mat)))
        cat(sprintf("  Labels: %d (unique: %d)\\n", length(labels), length(unique(labels))))

        cat("Running SingleR ...\\n")
        pred <- SingleR(
            test    = test_mat,
            ref     = ref_mat,
            labels  = labels,
            BPPARAM = SerialParam()
        )

        cat("Building output table ...\\n")
        result <- data.frame(
            cell_id         = colnames(test_mat),
            predicted_label = pred$labels,
            pruned_label    = as.character(pred$pruned.labels),
            delta_next      = pred$delta.next,
            row.names       = NULL,
            stringsAsFactors = FALSE
        )

        # Score columns (one per unique metacell label)
        scores_df <- as.data.frame(pred$scores)
        colnames(scores_df) <- paste0("score__", colnames(scores_df))
        result <- cbind(result, scores_df)

        write.table(result, "{out_csv}", sep=",", quote=FALSE, row.names=FALSE)
        cat(sprintf("Wrote %d rows x %d cols\\n", nrow(result), ncol(result)))
    """)

    r_script_path = str(tmpdir / "run_singleR.R")
    Path(r_script_path).write_text(r_script)

    print("Running SingleR in R ...")
    proc = subprocess.run(
        [RSCRIPT, "--vanilla", r_script_path],
        capture_output=True, text=True
    )
    print(proc.stdout)
    if proc.returncode != 0:
        print("STDERR:", proc.stderr[-3000:])
        raise RuntimeError(f"R script failed (exit {proc.returncode})")

    # ── 7. Read result and save ────────────────────────────────────────────────
    result_df = pd.read_csv(out_csv)

print(f"  Result shape: {result_df.shape}")
result_df.to_csv(OUT_TSV, sep="\t", index=False)
print(f"Saved {OUT_TSV}")
print("Done.")
