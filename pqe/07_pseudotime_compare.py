"""
Combine all momac pseudotime outputs into a single multi-indexed TSV.

Rows:    54 metacells (indexed by cell name)
Columns: two levels —
  level 0: "meta"                → shared obs metadata (sc_x, sc_y, cTET, etc.)
           "{tool}_{gene_set}"   → tool + gene-set combination
  level 1: column name within that group;
           the harmonized pseudotime is always labelled "pseudotime"
           Palantir also carries "entropy"

Outputs:
  results/07/pseudotime_momac_all.tsv      — combined pseudotime table
  results/07/gene_set_sizes.tsv            — genes used per variant after adata filtering
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc

RESULTS  = "/home/unix/cchu/projects/ZmanR/pqe/results"
REPO     = os.path.join(RESULTS, "06")
ADATA    = os.path.join(REPO, "zmanseq_momac_metacells_annot_clean.h5ad")
ADATA_MCUMI = os.path.join(REPO, "zmanseq_momac_metacells_annot_clean_mcumi.h5ad")
SPEARMAN_DIR = os.path.join(REPO, "ctet_mos/spearman")
OUT_DIR  = os.path.join(RESULTS, "07")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Tool definitions ────────────────────────────────────────────────────────
# Each entry: (file_stem, cell_id_col, pseudotime_col, extra_cols)
# cell_id_col=None means the row index is the cell name (Palantir TSVs)
TOOLS = {
    "scorpius":  ("scorpius_pseudotime.tsv",  "cell_id", "pseudotime",          []),
    "dpt":       ("dpt_pseudotime.tsv",        "cell_id", "pseudotime",          []),
    "redpath":   ("redpath_pseudotime.tsv",    "cell_id", "pseudotime",          []),
    "monocle2":  ("monocle2_pseudotime.tsv",   "cell_id", "pseudotime",          []),
    "palantir":  ("palantir_pseudotime.tsv",   None,      "palantir_pseudotime", ["palantir_entropy"]),
}

# ── Gene-set variants and their result directories ──────────────────────────
VARIANTS = {
    # Original momac variants (zmanseq_momac_metacells_annot.h5ad)
    "momac":            "pseudotime_momac",
    "momac_smooth":     "pseudotime_momac_smooth",
    "momac_hvg":        "pseudotime_momac_hvg",
    "momac_zman_s3":    "pseudotime_momac_zman_s3",
    "momac_full":       "pseudotime_momac_full",
    "momac_smooth_full":"pseudotime_momac_smooth_full",
    # scumi variants (zmanseq_momac_metacells_annot_clean.h5ad)
    "scumi_hvg_ctet":        "pseudotime_scumi_hvg_ctet",
    "scumi_hvg_smooth_ctet": "pseudotime_scumi_hvg_smooth_ctet",
    "scumi_full_ctet":       "pseudotime_scumi_full_ctet",
    "scumi_full_smooth_ctet":"pseudotime_scumi_full_smooth_ctet",
    "scumi_hvg":             "pseudotime_scumi_hvg",
    "scumi_zman_s3":         "pseudotime_scumi_zman_s3",
    # mcumi variants (zmanseq_momac_metacells_annot_clean_mcumi.h5ad)
    "mcumi_hvg_ctet":        "pseudotime_mcumi_hvg_ctet",
    "mcumi_hvg_smooth_ctet": "pseudotime_mcumi_hvg_smooth_ctet",
    "mcumi_full_ctet":       "pseudotime_mcumi_full_ctet",
    "mcumi_full_smooth_ctet":"pseudotime_mcumi_full_smooth_ctet",
    "mcumi_hvg":             "pseudotime_mcumi_hvg",
    "mcumi_zman_s3":         "pseudotime_mcumi_zman_s3",
}

# ── Gene-set size table ──────────────────────────────────────────────────────
# Load adata once to get var_names and highly_variable flags.
print("Loading adata for gene-set size computation...")
_ad       = sc.read_h5ad(ADATA)
adata_var = set(_ad.var_names)
# Compute HVGs the same way as spearman script
import scanpy as _sc; _ad_tmp = _ad.copy(); _sc.pp.highly_variable_genes(_ad_tmp, n_top_genes=1000)
hvg_genes_scumi = set(_ad_tmp.var_names[_ad_tmp.var["highly_variable"]]); del _ad_tmp
del _ad

_ad_mc = sc.read_h5ad(ADATA_MCUMI)
_ad_tmp2 = _ad_mc.copy(); _sc.pp.highly_variable_genes(_ad_tmp2, n_top_genes=1000)
hvg_genes_mcumi = set(_ad_tmp2.var_names[_ad_tmp2.var["highly_variable"]]); del _ad_tmp2, _ad_mc

def _sig(df, col="p_value", thresh=0.05):
    return set(df.index[df[col] < thresh])

def _sig_named(df, gene_col, pval_col, thresh=0.05):
    return set(df.loc[df[pval_col] < thresh, gene_col])

def _ctet_mos(label, geneset, timecol):
    igg   = pd.read_csv(f"{SPEARMAN_DIR}/{label}_igg_{geneset}_{timecol}_spearman.csv",   index_col=0)
    atrem = pd.read_csv(f"{SPEARMAN_DIR}/{label}_atrem_{geneset}_{timecol}_spearman.csv", index_col=0)
    return _sig(igg) & _sig(atrem)

_zman_s3 = (
    _sig_named(pd.read_csv("/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_aTREM2_Time.csv"),
               "Gene", "Pvalue") &
    _sig_named(pd.read_csv("/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_Isotype_Control_Time.csv"),
               "Gene", "Pvalue")
)

GENE_SET_SOURCES = {
    # Original momac variants
    "momac": (
        _sig(pd.read_csv(f"{REPO}/igg_hvg_spearman_df.naive.csv",   index_col=0)) &
        _sig(pd.read_csv(f"{REPO}/atrem_hvg_spearman_df.naive.csv", index_col=0))
    ),
    "momac_smooth": (
        _sig(pd.read_csv(f"{REPO}/smoothed_igg_hvg_spearman_df.csv",   index_col=0)) &
        _sig(pd.read_csv(f"{REPO}/smoothed_atrem_hvg_spearman_df.csv", index_col=0))
    ),
    "momac_hvg": hvg_genes_scumi,
    "momac_zman_s3": _zman_s3,
    "momac_full": (
        _sig(pd.read_csv(f"{REPO}/igg_full_spearman_df.naive.csv",   index_col=0)) &
        _sig(pd.read_csv(f"{REPO}/atrem_full_spearman_df.naive.csv", index_col=0))
    ),
    "momac_smooth_full": (
        _sig(pd.read_csv(f"{REPO}/smoothed_igg_full_spearman_df.csv",   index_col=0)) &
        _sig(pd.read_csv(f"{REPO}/smoothed_atrem_full_spearman_df.csv", index_col=0))
    ),
    # scumi variants
    "scumi_hvg_ctet":         _ctet_mos("scumi", "hvg",  "ctet"),
    "scumi_hvg_smooth_ctet":  _ctet_mos("scumi", "hvg",  "smooth_ctet"),
    "scumi_full_ctet":        _ctet_mos("scumi", "full", "ctet"),
    "scumi_full_smooth_ctet": _ctet_mos("scumi", "full", "smooth_ctet"),
    "scumi_hvg":              hvg_genes_scumi,
    "scumi_zman_s3":          _zman_s3,
    # mcumi variants
    "mcumi_hvg_ctet":         _ctet_mos("mcumi", "hvg",  "ctet"),
    "mcumi_hvg_smooth_ctet":  _ctet_mos("mcumi", "hvg",  "smooth_ctet"),
    "mcumi_full_ctet":        _ctet_mos("mcumi", "full", "ctet"),
    "mcumi_full_smooth_ctet": _ctet_mos("mcumi", "full", "smooth_ctet"),
    "mcumi_hvg":              hvg_genes_mcumi,
    "mcumi_zman_s3":          _zman_s3,
}

gene_size_rows = []
for variant_name, raw_genes in GENE_SET_SOURCES.items():
    filtered_genes = raw_genes & adata_var
    gene_size_rows.append({
        "gene_set":         variant_name,
        "n_genes_raw":      len(raw_genes),
        "n_genes_in_adata": len(filtered_genes),
    })

gene_sizes = pd.DataFrame(gene_size_rows).set_index("gene_set")
gene_sizes_path = os.path.join(OUT_DIR, "gene_set_sizes.tsv")
gene_sizes.to_csv(gene_sizes_path, sep="\t")
print("Gene set sizes:")
print(gene_sizes.to_string())
print(f"Saved to {gene_sizes_path}\n")

# ── Shared obs columns (everything that is NOT tool-specific) ───────────────
# Determined from the R-tool TSVs: all columns except cell_id and pseudotime
OBS_EXCLUDE = {"cell_id", "pseudotime",
               "palantir_pseudotime", "palantir_entropy",
               "DPT_root", "DPT_atrem", "DPT_igg"}

def load_r_tool(path, cell_id_col, pt_col, extra_cols):
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(cell_id_col)
    return df, pt_col, extra_cols

def load_palantir(path, pt_col, extra_cols):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df, pt_col, extra_cols

# ── Collect shared metadata from one reference file ─────────────────────────
ref_path = os.path.join(RESULTS, "pseudotime_momac", "scorpius_pseudotime.tsv")
ref = pd.read_csv(ref_path, sep="\t").set_index("cell_id")
meta_cols = [c for c in ref.columns if c not in OBS_EXCLUDE]
meta_df = ref[meta_cols].copy()
# Rename columns that R prefixed with X due to leading digits/underscores
meta_df.columns = [c.lstrip("X") if c.startswith("X__") or
                   (c.startswith("X") and c[1:2].isdigit()) else c
                   for c in meta_df.columns]

# Build multi-indexed metadata block
meta_mi = pd.concat({"meta": meta_df}, axis=1)

# ── Build tool blocks ────────────────────────────────────────────────────────
tool_frames = []

for variant_name, variant_dir in VARIANTS.items():
    variant_path = os.path.join(RESULTS, variant_dir)
    for tool_name, (fname, cell_id_col, pt_col, extra_cols) in TOOLS.items():
        fpath = os.path.join(variant_path, fname)
        if not os.path.exists(fpath):
            print(f"  MISSING: {fpath}")
            continue

        if cell_id_col is not None:
            df, pt_col_used, extras = load_r_tool(fpath, cell_id_col, pt_col, extra_cols)
            # For DPT, dynamically collect all DPT* columns (DPT1, DPT2, DPT3, ...)
            if tool_name == "dpt":
                extras = [c for c in df.columns if c.startswith("DPT_")]
        else:
            df, pt_col_used, extras = load_palantir(fpath, pt_col, extra_cols)

        # Harmonize pseudotime column name
        tool_df = df[[pt_col_used] + extras].rename(columns={pt_col_used: "pseudotime"})
        # Rename extras to strip tool-name prefix for cleanliness (e.g. palantir_entropy → entropy)
        tool_df.columns = [c.replace(f"{tool_name}_", "") for c in tool_df.columns]

        label = f"{tool_name}_{variant_name}"
        tool_mi = pd.concat({label: tool_df}, axis=1)
        tool_frames.append(tool_mi)
        print(f"  Loaded {label}: {list(tool_df.columns)}")

# ── Combine all ──────────────────────────────────────────────────────────────
combined = pd.concat([meta_mi] + tool_frames, axis=1)
combined.index.name = "cell_id"

out_path = os.path.join(OUT_DIR, "pseudotime_momac_all.tsv")
combined.to_csv(out_path, sep="\t")
print(f"\nSaved {combined.shape[0]} cells × {combined.shape[1]} columns to {out_path}")
print(f"Column groups: {combined.columns.get_level_values(0).unique().tolist()}")
