"""
Combine Zman-seq count matrices for cells listed in metadata.txt and save as AnnData h5ad.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


# ── 1. Parse metadata ────────────────────────────────────────────────────────

def preprocess_metadata(path):
    """
    Parse metadata.txt, handling spaces in Mouse and cluster_colors columns.
    Uses regex anchors for Mouse (Mouse_3 / Mouse 3) and time_assignment (12H / Negative / Not_Assigned).
    """
    mouse_pat = re.compile(r'(Mouse[\s_]\d+)')
    time_pat  = re.compile(r'\b(\d+H|Negative|Not_Assigned)\b')

    with open(path) as f:
        lines = f.readlines()

    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        raw = line.strip()
        mouse_m = mouse_pat.search(raw)
        time_m  = time_pat.search(raw)

        mouse   = mouse_m.group(1).replace(" ", "_") if mouse_m else None
        time_as = time_m.group(1) if time_m else None

        left  = raw[:mouse_m.start()].strip().split() if mouse_m else []
        after = raw[time_m.end():].strip().split()    if time_m  else []

        # after time_assignment: [cluster_colors..., sc_x, sc_y, Treatment]
        cluster   = " ".join(after[:-3]) if len(after) > 3 else (after[0] if after else None)
        sc_x, sc_y, treatment = (after[-3], after[-2], after[-1]) if len(after) >= 3 else (None, None, None)

        rows.append(left + [mouse, time_as, cluster, sc_x, sc_y, treatment])

    cols = ["row_index"] + header
    meta = pd.DataFrame(rows, columns=cols).set_index("row_index")
    meta["sc_x"] = pd.to_numeric(meta["sc_x"], errors="coerce")
    meta["sc_y"] = pd.to_numeric(meta["sc_y"], errors="coerce")
    return meta


# ── 2. Load count data ────────────────────────────────────────────────────────

def load_counts_for_metadata(count_dir, meta):
    """
    For each Amp_batch_ID in meta, find the matching .txt.gz file and load counts
    for cells present in meta. Returns a combined (cells × genes) sparse matrix
    along with aligned cell and gene lists.
    """
    count_dir = Path(count_dir)
    batches   = meta["Amp_batch_ID"].unique()
    well_set  = set(meta["Well_ID"])

    matrices = []  # list of (cells, sparse matrix)
    gene_index = None

    for batch in batches:
        matches = list(count_dir.glob(f"*_{batch}.txt.gz"))
        if not matches:
            print(f"Warning: no count file found for batch {batch}, skipping")
            continue
        path = matches[0]
        print(f"Loading {path.name} ...")

        df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

        # keep only wells that are in the metadata
        wells_in_file = [w for w in df.columns if w in well_set]
        df = df[wells_in_file]

        # ensure consistent gene order across batches
        if gene_index is None:
            gene_index = df.index
        else:
            df = df.loc[gene_index]

        # genes × cells → transpose to cells × genes, store as sparse
        mat = sp.csr_matrix(df.values.T.astype(np.float32))
        matrices.append((wells_in_file, mat))

    all_cells = [w for wells, _ in matrices for w in wells]
    combined  = sp.vstack([mat for _, mat in matrices])
    return all_cells, gene_index.tolist(), combined


# ── 3. Build and save AnnData ─────────────────────────────────────────────────

def main():
    data_dir  = Path("/mnt/thechenlab/ClaudiaC/zmanseq")
    meta_path = data_dir / "metadata.txt"
    out_path  = Path("pqe/results/04/zmanseq.h5ad")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Parsing metadata ...")
    meta = preprocess_metadata(meta_path)

    cells, genes, matrix = load_counts_for_metadata(data_dir, meta)

    print(f"Building AnnData: {len(cells)} cells × {len(genes)} genes")
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes

    # attach metadata, aligned to the cell order in adata
    meta_indexed = meta.set_index("Well_ID")
    adata.obs = meta_indexed.loc[cells].copy()

    print(f"Saving to {out_path} ...")
    adata.write_h5ad(out_path)
    print("Done.")


if __name__ == "__main__":
    main()
