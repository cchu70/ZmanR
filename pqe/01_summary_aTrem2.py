"""
01_summary_aTrem2.py

Python reproduction of 01_summary_aTrem2.R.
Loads cells_obs.csv + marker_counts.csv (written by the R script), builds an
AnnData object saved as cells_markers.h5ad, then reproduces all summary tables
and figures in results/01/.

Requirements:  scanpy pandas numpy matplotlib seaborn scipy
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from pathlib import Path

# ── output directory ────────────────────────────────────────────────────────────
outdir = Path("results/01/py_version")
outdir.mkdir(parents=True, exist_ok=True)

# ── color palettes (matching R script) ─────────────────────────────────────────
pal_treat = {"aTrem2": "#E64B35", "IgG": "#4DBBD5"}

pal_celltype = {
    "Treg": "#F4C0C8", "CD4": "#F5B48A", "CD8": "#E07878",
    "CD8_Dysfunctional": "#C85050",
    "NK_Chemotactic": "#F0C060", "Chemotactic": "#F0C060",
    "NK_Dysfunctional": "#E89030", "Dysfunctional": "#E89030",
    "NK_intermediate": "#E07820", "Cytotoxic": "#C84010",
    "MigDC": "#D45030", "cDC1": "#E8B030", "cDC2": "#F0E080",
    "MonDC": "#F8D060", "pDC": "#D8A8D0",
    "inflammatory": "#7050B0", "Monocytes": "#8070C0",
    "Monocyte": "#9080C8", "MoMac": "#A088C8",
    "MoMac1": "#9878B8", "MoMac2": "#B098D0", "transitory": "#C0B0D8",
    "TAM": "#70CED8", "Acp5_TAM": "#50B8C8", "Arg1_TAM": "#80D8C8",
    "Gpnmb_TAM": "#90C8B8",
    "Not_annotated": "#C8C8C8", "Doublet": "#909090",
}

# ── marker gene lists ───────────────────────────────────────────────────────────
t_markers_general    = ["Cd3e", "Cd3d", "Cd3g", "Trac", "Trbc1", "Trbc2"]
t_markers_lineage    = ["Cd4", "Cd8a", "Cd8b1"]
t_markers_activation = ["Cd44", "Cd69", "Il2ra", "Ifng", "Tnf"]
t_markers_effector   = ["Gzmb", "Gzma", "Prf1", "Fasl"]
t_markers_exhaustion = ["Pdcd1", "Havcr2", "Lag3", "Tigit", "Tox", "Ctla4", "Entpd1"]

tumor_markers_neural   = ["Nes", "Sox2", "Prom1"]
tumor_markers_glioma   = ["Gfap", "Olig2", "S100b"]
tumor_markers_prolif   = ["Mki67", "Top2a"]
tumor_markers_invasive = ["Vim", "Fn1", "Cd44"]
tumor_markers_receptor = ["Egfr", "Pdgfra", "Met"]

t_markers     = (t_markers_general + t_markers_lineage + t_markers_activation
                 + t_markers_effector + t_markers_exhaustion)
tumor_markers = (tumor_markers_neural + tumor_markers_glioma + tumor_markers_prolif
                 + tumor_markers_invasive + tumor_markers_receptor)
all_markers   = t_markers + tumor_markers

marker_subtype = (
    {g: "T cell: general"       for g in t_markers_general}
    | {g: "T cell: lineage"     for g in t_markers_lineage}
    | {g: "T cell: activation"  for g in t_markers_activation}
    | {g: "T cell: effector"    for g in t_markers_effector}
    | {g: "T cell: exhaustion"  for g in t_markers_exhaustion}
    | {g: "Tumor: neural stem"  for g in tumor_markers_neural}
    | {g: "Tumor: glioma"       for g in tumor_markers_glioma}
    | {g: "Tumor: proliferation" for g in tumor_markers_prolif}
    | {g: "Tumor: invasive"     for g in tumor_markers_invasive}
    | {g: "Tumor: receptor"     for g in tumor_markers_receptor}
)
marker_type_map = (
    {g: "T cell" for g in t_markers}
    | {g: "Tumor/Glia (GL261)" for g in tumor_markers}
)

subtype_levels = [
    "T cell: general", "T cell: lineage", "T cell: activation",
    "T cell: effector", "T cell: exhaustion",
    "Tumor: neural stem", "Tumor: glioma", "Tumor: proliferation",
    "Tumor: invasive", "Tumor: receptor",
]

# ==============================================================================
# LOAD / BUILD ANNDATA
# ==============================================================================

r_outdir  = Path("results/01")   # where R wrote cells_obs.csv / marker_counts.csv
h5ad_path = outdir / "cells_markers.h5ad"

if h5ad_path.exists():
    print("Loading existing AnnData from", h5ad_path)
    adata = sc.read_h5ad(h5ad_path)
else:
    print("Building AnnData from CSV exports...")
    obs  = pd.read_csv(r_outdir / "cells_obs.csv").set_index("Well_ID")
    cnts = pd.read_csv(r_outdir / "marker_counts.csv").set_index("Well_ID")

    # keep genes present in count matrix; deduplicate (Cd44 appears in both T cell and tumor lists)
    genes = list(dict.fromkeys(g for g in all_markers if g in cnts.columns))
    # align rows
    common = obs.index.intersection(cnts.index)
    obs  = obs.loc[common]
    cnts = cnts.loc[common, genes]

    var = pd.DataFrame({
        "marker_type":    [marker_type_map.get(g, "unknown") for g in genes],
        "marker_subtype": [marker_subtype.get(g, "unknown")  for g in genes],
    }, index=genes)

    adata = sc.AnnData(X=cnts.values.astype(np.float32), obs=obs, var=var)

    # tidy dtypes
    adata.obs["is_doublet"] = adata.obs["is_doublet"].astype(bool)
    for col in ["Treatment", "celltype", "time_assignment", "cell_class",
                "Amp_batch_ID", "mc_celltype"]:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")
    for col in ["mc_id_num"]:
        if col in adata.obs:
            adata.obs[col] = pd.to_numeric(adata.obs[col], errors="coerce")

    adata.write_h5ad(h5ad_path)
    print("Saved", h5ad_path)

print(adata)

# Convenience references
merged   = adata.obs.copy()
merged["total_UMI"] = pd.to_numeric(merged["total_UMI"], errors="coerce")
merged["genes_det"] = pd.to_numeric(merged["genes_det"],  errors="coerce")

genes_present = list(adata.var_names)
expr_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=genes_present)

# ==============================================================================
# 3. SUMMARY TABLES
# ==============================================================================

print("\nWriting summary tables...")

## 3a. Per-plate summary
plate_summary = (
    merged.groupby(["Treatment", "Amp_batch_ID"], observed=True)
    .agg(
        n_cells      = ("total_UMI", "count"),
        median_UMI   = ("total_UMI", "median"),
        mean_UMI     = ("total_UMI", "mean"),
        median_genes = ("genes_det", "median"),
        mean_genes   = ("genes_det", "mean"),
        pct_zero_umi = ("total_UMI", lambda x: round((x == 0).mean() * 100, 1)),
    )
    .reset_index()
    .sort_values(["Treatment", "Amp_batch_ID"])
)
plate_summary["mean_UMI"]   = plate_summary["mean_UMI"].round(1)
plate_summary["mean_genes"] = plate_summary["mean_genes"].round(1)
plate_summary.to_csv(outdir / "plate_summary.csv", index=False)

## 3b. Per-treatment summary
treatment_summary = (
    merged.groupby("Treatment", observed=True)
    .agg(
        n_cells      = ("total_UMI", "count"),
        n_plates     = ("Amp_batch_ID", "nunique"),
        median_UMI   = ("total_UMI", "median"),
        mean_UMI     = ("total_UMI", "mean"),
        median_genes = ("genes_det", "median"),
        mean_genes   = ("genes_det", "mean"),
        pct_zero_umi = ("total_UMI", lambda x: round((x == 0).mean() * 100, 1)),
    )
    .reset_index()
)
treatment_summary["mean_UMI"]   = treatment_summary["mean_UMI"].round(1)
treatment_summary["mean_genes"] = treatment_summary["mean_genes"].round(1)
treatment_summary.to_csv(outdir / "treatment_summary.csv", index=False)

## 3c. Cell type composition by treatment
celltype_summary = (
    merged.dropna(subset=["celltype"])
    .groupby(["Treatment", "celltype"], observed=True)
    .size().reset_index(name="n")
)
celltype_summary["pct"] = (
    celltype_summary.groupby("Treatment", observed=True)["n"]
    .transform(lambda x: (x / x.sum() * 100).round(2))
)
celltype_summary = celltype_summary.sort_values(["Treatment", "n"], ascending=[True, False])
celltype_summary.to_csv(outdir / "celltype_by_treatment.csv", index=False)

## 3d. Time assignment by treatment
time_summary = (
    merged.dropna(subset=["time_assignment"])
    .groupby(["Treatment", "time_assignment"], observed=True)
    .size().reset_index(name="n")
)
time_summary["pct"] = (
    time_summary.groupby("Treatment", observed=True)["n"]
    .transform(lambda x: (x / x.sum() * 100).round(2))
)
time_summary.to_csv(outdir / "time_assignment_by_treatment.csv", index=False)

## 3e. Doublet summary per plate
doublet_summary = (
    merged.groupby(["Treatment", "Amp_batch_ID"], observed=True)
    .agg(
        n_cells      = ("is_doublet", "count"),
        n_doublets   = ("is_doublet", "sum"),
        pct_doublets = ("is_doublet", lambda x: round(x.mean() * 100, 2)),
    )
    .reset_index()
    .sort_values(["Treatment", "Amp_batch_ID"])
)
doublet_summary.to_csv(outdir / "doublet_summary.csv", index=False)

by_treat = (
    merged.groupby("Treatment", observed=True)["is_doublet"]
    .agg(n_doublets="sum", pct=lambda x: round(x.mean() * 100, 2))
    .reset_index()
)
print("Doublets by treatment:\n", by_treat.to_string(index=False))
print("Tables written.")

# ==============================================================================
# 4. FIGURES
# ==============================================================================

print("\nGenerating figures (section 4)...")

sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 11})

TREAT_ORDER = ["aTrem2", "IgG"]


def savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def violin_box(ax, data, x, y, palette, log=False, xlabel="", ylabel="", title=""):
    """Violin + overlaid thin boxplot."""
    order = [t for t in TREAT_ORDER if t in data[x].values]
    tmp   = data.copy()
    tmp["_y"] = np.log10(tmp[y] + 1) if log else tmp[y]
    sns.violinplot(data=tmp, x=x, y="_y", palette=palette, order=order,
                   inner=None, alpha=0.7, scale="width", ax=ax, linewidth=0.8)
    sns.boxplot(data=tmp, x=x, y="_y", order=order, width=0.12,
                flierprops=dict(markersize=1), color="white",
                boxprops=dict(alpha=0.9), linewidth=0.8, ax=ax)
    ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.legend([], frameon=False)


## -- QC by treatment ----------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(6, 5))
violin_box(axes[0], merged, "Treatment", "total_UMI", pal_treat,
           log=True, xlabel="", ylabel="log10(total UMI + 1)",
           title="Total UMI per cell by treatment")
violin_box(axes[1], merged, "Treatment", "genes_det", pal_treat,
           log=False, xlabel="", ylabel="Genes detected",
           title="Genes detected per cell by treatment")
fig.tight_layout()
savefig(fig, outdir / "qc_by_treatment.pdf")

## -- QC by plate --------------------------------------------------------------
plate_treat_map = merged.groupby("Amp_batch_ID", observed=True)["Treatment"].first()
plate_ord_umi   = (merged.groupby("Amp_batch_ID", observed=True)["total_UMI"]
                   .median().sort_values().index.tolist())
plate_ord_genes = (merged.groupby("Amp_batch_ID", observed=True)["genes_det"]
                   .median().sort_values().index.tolist())

fig, axes = plt.subplots(2, 1, figsize=(14, 8))
for ax, col, log, order, ylabel, title in [
    (axes[0], "total_UMI", True,  plate_ord_umi,
     "log10(total UMI + 1)", "Total UMI per cell by plate"),
    (axes[1], "genes_det",  False, plate_ord_genes,
     "Genes detected",       "Genes detected per cell by plate"),
]:
    tmp = merged.copy()
    tmp["_y"] = np.log10(tmp[col] + 1) if log else tmp[col]
    pal = {p: pal_treat[plate_treat_map[p]] for p in order if p in plate_treat_map}
    sns.boxplot(data=tmp, x="Amp_batch_ID", y="_y", order=order,
                palette=pal, flierprops=dict(markersize=1),
                linewidth=0.7, ax=ax)
    ax.set_title(title); ax.set_xlabel("Plate (Amp batch)"); ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    handles = [mpatches.Patch(color=c, label=t) for t, c in pal_treat.items()]
    ax.legend(handles=handles, title="Treatment", fontsize=8)
fig.tight_layout()
savefig(fig, outdir / "qc_by_plate.pdf")

## -- Cell type composition by treatment ---------------------------------------
ct_pivot = (
    celltype_summary
    .pivot(index="celltype", columns="Treatment", values="pct")
    .fillna(0)
    .reindex(columns=[t for t in TREAT_ORDER if t in celltype_summary["Treatment"].values])
)
colors = [pal_celltype.get(str(ct), "#C8C8C8") for ct in ct_pivot.index]

fig, ax = plt.subplots(figsize=(7, 5))
ct_pivot.T.plot.bar(stacked=True, color=colors, ax=ax, legend=False, width=0.5)
ax.set_xlabel(None); ax.set_ylabel("% of cells")
ax.set_title("Cell type composition by treatment")
ax.tick_params(axis="x", rotation=0)
handles = [mpatches.Patch(color=pal_celltype.get(str(ct), "#C8C8C8"), label=str(ct))
           for ct in ct_pivot.index]
ax.legend(handles=handles, title="Cell type", bbox_to_anchor=(1.01, 1),
          loc="upper left", fontsize=7)
fig.tight_layout()
savefig(fig, outdir / "celltype_composition.pdf")

## -- Cell type by plate -------------------------------------------------------
ct_plate = (
    merged.dropna(subset=["celltype"])
    .groupby(["Treatment", "Amp_batch_ID", "celltype"], observed=True)
    .size().reset_index(name="n")
)
ct_plate["pct"] = (
    ct_plate.groupby(["Treatment", "Amp_batch_ID"], observed=True)["n"]
    .transform(lambda x: x / x.sum() * 100)
)

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for ax, treat in zip(axes, TREAT_ORDER):
    sub    = ct_plate[ct_plate["Treatment"] == treat]
    plates = sorted(sub["Amp_batch_ID"].unique())
    mat    = (sub.pivot(index="Amp_batch_ID", columns="celltype", values="pct")
              .fillna(0).reindex(plates))
    colors_p = [pal_celltype.get(str(ct), "#C8C8C8") for ct in mat.columns]
    mat.plot.bar(stacked=True, color=colors_p, ax=ax, legend=False, width=0.8)
    ax.set_title(f"Cell type by plate — {treat}")
    ax.set_xlabel("Plate"); ax.set_ylabel("% of cells")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
fig.tight_layout()
savefig(fig, outdir / "celltype_by_plate.pdf")

## -- Time assignment by treatment ---------------------------------------------
ta_order = ["12H", "24H", "36H", "Negative", "Not_Assigned"]
time_pivot = (
    time_summary
    .pivot(index="time_assignment", columns="Treatment", values="pct")
    .fillna(0)
    .reindex(index=[t for t in ta_order if t in time_summary["time_assignment"].values])
    .reindex(columns=[t for t in TREAT_ORDER if t in time_summary["Treatment"].values])
)
ta_colors = sns.color_palette("Set2", n_colors=len(time_pivot))

fig, ax = plt.subplots(figsize=(6, 5))
time_pivot.T.plot.bar(stacked=True, color=ta_colors, ax=ax, width=0.5)
ax.set_xlabel(None); ax.set_ylabel("% of cells")
ax.set_title("Time assignment by treatment")
ax.tick_params(axis="x", rotation=0)
ax.legend(title="Time bin", bbox_to_anchor=(1.01, 1), loc="upper left")
fig.tight_layout()
savefig(fig, outdir / "time_assignment.pdf")

## -- UMI vs genes scatter -----------------------------------------------------
np.random.seed(42)
samp = merged.sample(n=min(5000, len(merged)))

fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
for ax, treat in zip(axes, TREAT_ORDER):
    sub = samp[samp["Treatment"] == treat]
    ax.scatter(np.log10(sub["total_UMI"] + 1), sub["genes_det"],
               alpha=0.3, s=2, color=pal_treat[treat])
    ax.set_title(treat)
    ax.set_xlabel("log10(total UMI + 1)")
    ax.set_ylabel("Genes detected" if treat == TREAT_ORDER[0] else "")
fig.suptitle("UMI vs genes detected (5k cell sample)")
fig.tight_layout()
savefig(fig, outdir / "umi_vs_genes.pdf")

## -- Doublet QC ---------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
for row, (col, log, ylabel, title) in enumerate([
    ("total_UMI", True,  "log10(total UMI + 1)", "UMI: doublets vs singlets"),
    ("genes_det", False, "Genes detected",        "Genes detected: doublets vs singlets"),
]):
    for col_i, treat in enumerate(TREAT_ORDER):
        ax  = axes[row][col_i]
        sub = merged[merged["Treatment"] == treat].copy()
        sub["_y"]  = np.log10(sub[col] + 1) if log else sub[col]
        sub["_cl"] = sub["is_doublet"].map({True: "Doublet", False: "Singlet"})
        sns.violinplot(data=sub, x="_cl", y="_y", order=["Singlet", "Doublet"],
                       palette={"Singlet": "grey", "Doublet": "#FF7F00"},
                       inner=None, alpha=0.7, scale="width", ax=ax, linewidth=0.8)
        sns.boxplot(data=sub, x="_cl", y="_y", order=["Singlet", "Doublet"],
                    width=0.12, color="white", flierprops=dict(markersize=1),
                    linewidth=0.8, ax=ax)
        ax.set_title(f"{title}\n({treat})", fontsize=9)
        ax.set_xlabel(""); ax.set_ylabel(ylabel if col_i == 0 else "")
        ax.legend([], frameon=False)
fig.tight_layout()
savefig(fig, outdir / "doublet_qc.pdf")

print("Figures written (section 4).")

# ==============================================================================
# 5. DOUBLET MARKER EXPRESSION
# ==============================================================================

print("\nComputing doublet marker expression...")

# Build a per-cell marker expression dataframe with metadata
marker_merged = expr_df.copy()
marker_merged = marker_merged.join(merged[["Treatment", "Amp_batch_ID",
                                           "celltype", "is_doublet", "cell_class"]])

# Long-format
marker_long = (
    marker_merged
    .reset_index()
    .rename(columns={"index": "Well_ID"})
    .melt(id_vars=["Well_ID", "Treatment", "cell_class"],
          value_vars=list(dict.fromkeys(g for g in all_markers if g in marker_merged.columns)),
          var_name="gene", value_name="expr")
)
marker_long["marker_type"]    = marker_long["gene"].map(marker_type_map)
marker_long["marker_subtype"] = marker_long["gene"].map(marker_subtype)

# Summary per gene × treatment × cell_class
marker_summary = (
    marker_long
    .groupby(["Treatment", "cell_class", "marker_type", "marker_subtype", "gene"])
    .agg(
        mean_expr = ("expr", lambda x: round(x.mean(), 4)),
        pct_expr  = ("expr", lambda x: round((x > 0).mean() * 100, 2)),
    )
    .reset_index()
)
marker_summary["marker_subtype"] = pd.Categorical(
    marker_summary["marker_subtype"], categories=subtype_levels, ordered=True)
marker_summary = marker_summary.sort_values(["marker_subtype", "gene"])
marker_summary.to_csv(outdir / "doublet_marker_summary.csv", index=False)

## -- Bar chart: % expressing each marker, doublets vs singlets ----------------
fig, axes = plt.subplots(len(subtype_levels), 2,
                         figsize=(12, len(subtype_levels) * 1.6),
                         squeeze=False)
for row, subtype in enumerate(subtype_levels):
    for col_i, treat in enumerate(TREAT_ORDER):
        ax  = axes[row][col_i]
        sub = marker_summary[
            (marker_summary["marker_subtype"] == subtype)
            & (marker_summary["Treatment"] == treat)
        ]
        genes = sub["gene"].unique()
        x     = np.arange(len(genes))
        w     = 0.35
        for i, (cl, color) in enumerate([("Singlet", "grey"), ("Doublet", "#FF7F00")]):
            vals = (sub[sub["cell_class"] == cl]
                    .groupby("gene")["pct_expr"].first()
                    .reindex(genes).fillna(0))
            ax.bar(x + i * w, vals.values, width=w, color=color,
                   label=cl if (row == 0 and col_i == 0) else "_nolegend_")
        ax.set_xticks(x + w / 2)
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=7)
        ax.set_title(f"{subtype} | {treat}", fontsize=7)
        ax.set_ylabel("% cells > 0" if col_i == 0 else "", fontsize=7)

fig.legend(["Singlet", "Doublet"], loc="upper right", fontsize=9)
fig.suptitle("% cells expressing T cell / GL261 tumor markers: doublets vs singlets",
             y=1.002)
fig.tight_layout()
savefig(fig, outdir / "doublet_markers.pdf")

## -- Heatmap: z-scored mean expression ----------------------------------------
heat_df = (
    marker_summary
    .assign(label=lambda d: d["Treatment"].astype(str) + "\n" + d["cell_class"].astype(str))
    .pivot_table(index="gene", columns="label", values="mean_expr")
    .fillna(0)
)
row_std  = heat_df.std(axis=1).replace(0, 1)
heat_scaled = heat_df.sub(heat_df.mean(axis=1), axis=0).div(row_std, axis=0).fillna(0)

fig, ax = plt.subplots(figsize=(6, max(5, len(heat_scaled) * 0.22)))
sns.heatmap(heat_scaled, cmap="RdBu_r", center=0, ax=ax,
            xticklabels=True, yticklabels=True,
            cbar_kws={"label": "z-score"}, linewidths=0.3)
ax.set_title("Marker expression (z-scored)\ndoublets vs singlets")
ax.tick_params(axis="y", labelsize=6)
fig.tight_layout()
savefig(fig, outdir / "doublet_marker_heatmap.pdf")

print("Figures written (section 5).")

# ==============================================================================
# 6. CO-EXPRESSION OF T CELL AND TUMOR MARKERS
# ==============================================================================

print("\nComputing co-expression plots...")

# Per-cell aggregate scores
def score(genes):
    cols = [g for g in genes if g in expr_df.columns]
    return expr_df[cols].sum(axis=1) if cols else pd.Series(0, index=expr_df.index)

merged["t_score"]        = score(t_markers)
merged["tumor_score"]    = score(tumor_markers)
merged["t_general"]      = score(t_markers_general)
merged["t_lineage"]      = score(t_markers_lineage)
merged["t_activation"]   = score(t_markers_activation)
merged["t_effector"]     = score(t_markers_effector)
merged["t_exhaustion"]   = score(t_markers_exhaustion)
merged["tumor_neural"]   = score(tumor_markers_neural)
merged["tumor_glioma"]   = score(tumor_markers_glioma)
merged["tumor_prolif"]   = score(tumor_markers_prolif)
merged["tumor_invasive"] = score(tumor_markers_invasive)
merged["tumor_receptor"] = score(tumor_markers_receptor)

## -- Plot 1: scatter T score vs tumor score -----------------------------------
np.random.seed(42)
singlets = merged[~merged["is_doublet"]].sample(n=min(3000, (~merged["is_doublet"]).sum()))
doublets = merged[merged["is_doublet"]]
plot_df  = pd.concat([singlets, doublets])

fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
for ax, treat in zip(axes, TREAT_ORDER):
    sub = plot_df[plot_df["Treatment"] == treat]
    for cl, (color, alpha, s) in {
        "Singlet": ("grey",      0.3,  3),
        "Doublet": ("#FF7F00",  0.85,  8),
    }.items():
        m = sub["cell_class"] == cl
        ax.scatter(np.log1p(sub.loc[m, "t_score"]),
                   np.log1p(sub.loc[m, "tumor_score"]),
                   c=color, alpha=alpha, s=s, label=cl)
    ax.set_title(treat)
    ax.set_xlabel("log1p(T cell marker UMI sum)")
    ax.set_ylabel("log1p(Tumor marker UMI sum)" if treat == TREAT_ORDER[0] else "")
    ax.legend(fontsize=8)
fig.suptitle("T cell score vs GL261 tumor score per cell\n"
             "(Singlets: 3k sampled + all doublets)")
fig.tight_layout()
savefig(fig, outdir / "coexpr_scatter.pdf")

## -- Plot 2: subtype score heatmap --------------------------------------------
score_cols   = ["t_general", "t_lineage", "t_activation", "t_effector", "t_exhaustion",
                "tumor_neural", "tumor_glioma", "tumor_prolif", "tumor_invasive", "tumor_receptor"]
score_labels = ["T: general", "T: lineage", "T: activation", "T: effector", "T: exhaustion",
                "Tumor: neural stem", "Tumor: glioma", "Tumor: proliferation",
                "Tumor: invasive", "Tumor: receptor"]

group_levels = ["aTrem2\nDoublet", "aTrem2\nSinglet", "IgG\nDoublet", "IgG\nSinglet"]
score_rows = []
for (treat, cl), grp in merged.groupby(["Treatment", "cell_class"], observed=True):
    treat, cl = str(treat), str(cl)
    row = {"Treatment": treat, "cell_class": cl,
           "group": f"{treat}\n{cl}"}
    for col, lbl in zip(score_cols, score_labels):
        row[lbl] = round((grp[col] > 0).mean() * 100, 2)
    score_rows.append(row)

score_summary = pd.DataFrame(score_rows)
score_summary["group"] = pd.Categorical(score_summary["group"],
                                        categories=group_levels)
heat_mat = (score_summary.set_index("group")[score_labels]
            .reindex(group_levels).T
            .reindex(score_labels[::-1]))

fig, ax = plt.subplots(figsize=(7, 5))
im = ax.imshow(heat_mat.values.astype(float), aspect="auto",
               cmap=plt.cm.YlOrRd, vmin=0)
plt.colorbar(im, ax=ax, label="% cells expressing")
ax.set_xticks(range(len(group_levels)))
ax.set_xticklabels(group_levels, fontsize=9)
ax.set_yticks(range(len(heat_mat.index)))
ax.set_yticklabels(heat_mat.index, fontsize=9)
for r in range(heat_mat.shape[0]):
    for c in range(heat_mat.shape[1]):
        val = heat_mat.values[r, c]
        ax.text(c, r, f"{float(val):.1f}", ha="center", va="center",
                fontsize=8, color="white" if float(val) > 40 else "black")
ax.set_title("% cells with any expression in marker group")
fig.tight_layout()
savefig(fig, outdir / "coexpr_score_heatmap.pdf")

## -- Plot 3: pairwise co-expression bubble plot -------------------------------
t_subtypes = {
    "T: general":    t_markers_general,
    "T: lineage":    t_markers_lineage,
    "T: activation": t_markers_activation,
    "T: effector":   t_markers_effector,
    "T: exhaustion": t_markers_exhaustion,
}
tumor_subtypes = {
    "Tumor: neural stem":   tumor_markers_neural,
    "Tumor: glioma":        tumor_markers_glioma,
    "Tumor: proliferation": tumor_markers_prolif,
    "Tumor: invasive":      tumor_markers_invasive,
    "Tumor: receptor":      tumor_markers_receptor,
}

coexpr_rows = []
for tname, tgenes in t_subtypes.items():
    t_pos = score(tgenes) > 0
    for tuname, tugenes in tumor_subtypes.items():
        tu_pos = score(tugenes) > 0
        coexpr = t_pos & tu_pos
        for (treat, cl), grp in merged.groupby(["Treatment", "cell_class"], observed=True):
            pct = round(coexpr.loc[grp.index].mean() * 100, 2)
            coexpr_rows.append({"t_subtype": tname, "tumor_subtype": tuname,
                                 "Treatment": treat, "cell_class": cl,
                                 "pct_coexpr": pct})

coexpr_df = pd.DataFrame(coexpr_rows)
coexpr_df["group"] = pd.Categorical(
    coexpr_df["Treatment"].astype(str) + "\n" + coexpr_df["cell_class"].astype(str),
    categories=group_levels)
coexpr_df.to_csv(outdir / "coexpr_summary.csv", index=False)

t_names  = list(t_subtypes.keys())
tu_names = list(tumor_subtypes.keys())

fig, axes = plt.subplots(1, 4, figsize=(13, 5), sharey=True, sharex=True)
for ax, grp in zip(axes, group_levels):
    sub = coexpr_df[coexpr_df["group"] == grp]
    mat = (sub.pivot(index="t_subtype", columns="tumor_subtype", values="pct_coexpr")
           .reindex(index=t_names, columns=tu_names).fillna(0))
    max_val = coexpr_df["pct_coexpr"].max() or 1
    for ri, t in enumerate(t_names):
        for ci, tu in enumerate(tu_names):
            val  = float(mat.loc[t, tu])
            size = (val / max_val) * 400
            rgba = plt.cm.Reds(min(val / max(max_val, 1), 1))
            ax.scatter([ci], [ri], s=max(size, 4), c=[rgba], vmin=0, vmax=max_val,
                       edgecolors="grey", linewidths=0.3)
            if val > 0:
                ax.annotate(f"{val:.1f}", (ci, ri), ha="center", va="bottom",
                            fontsize=6)
    ax.set_xticks(range(len(tu_names)))
    ax.set_xticklabels(tu_names, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(t_names)))
    ax.set_yticklabels(t_names, fontsize=8)
    ax.set_title(grp, fontsize=9)
    ax.set_xlim(-0.6, len(tu_names) - 0.4)
    ax.set_ylim(-0.6, len(t_names)  - 0.4)
fig.suptitle("Co-expression: T cell subtype × GL261 tumor subtype")
fig.tight_layout()
savefig(fig, outdir / "coexpr_bubble.pdf")

# ==============================================================================
# 7. METACELL QC
# ==============================================================================

print("\nComputing metacell QC metrics...")

mc_qc_path = r_outdir / "metacell_qc.csv"
if not mc_qc_path.exists():
    print(f"  {mc_qc_path} not found — skipping metacell QC (run R script first).")
else:
    mc_qc = pd.read_csv(mc_qc_path)

    # -- Add mc_qc as uns to AnnData and re-save ---------------------------------
    adata.uns["metacell_qc"] = mc_qc.to_dict(orient="list")
    adata.write_h5ad(h5ad_path)
    print(f"  Saved metacell_qc to adata.uns in {h5ad_path}")
    print(f"  {len(mc_qc)} metacells, median {mc_qc['n_cells'].median():.0f} cells/mc")

    # -- Table: cells per metacell summary by treatment --------------------------
    mc_merged = merged.dropna(subset=["mc_id_num"]).copy()
    mc_merged["mc_id_num"] = mc_merged["mc_id_num"].astype(int)

    mc_cell_summary = (
        mc_merged.groupby(["Treatment", "mc_id_num"], observed=True)
        .agg(n_cells=("mc_id_num", "count"))
        .reset_index()
        .merge(mc_qc[["mc_id_num", "celltype", "mc_x", "mc_y"]], on="mc_id_num", how="left")
    )
    mc_cell_summary.to_csv(outdir / "metacell_cell_summary.csv", index=False)

    # -- Table: metacells per cell type ------------------------------------------
    mc_per_ct = (
        mc_qc.groupby("celltype")
        .agg(n_metacells=("mc_id_num", "count"),
             total_cells=("n_cells", "sum"),
             median_cells_per_mc=("n_cells", "median"))
        .reset_index()
        .sort_values("n_metacells", ascending=False)
    )
    mc_per_ct.to_csv(outdir / "metacells_per_celltype.csv", index=False)

    # -- Figure 1: cells per metacell histogram by treatment ---------------------
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
    for ax, treat in zip(axes, TREAT_ORDER):
        sub = mc_cell_summary[mc_cell_summary["Treatment"] == treat]
        ax.hist(sub["n_cells"], bins=30, color=pal_treat[treat], edgecolor="white",
                linewidth=0.4)
        ax.axvline(sub["n_cells"].median(), color="black", linestyle="--",
                   linewidth=1, label=f"median={sub['n_cells'].median():.0f}")
        ax.set_title(treat)
        ax.set_xlabel("Cells per metacell")
        ax.set_ylabel("# metacells" if treat == TREAT_ORDER[0] else "")
        ax.legend(fontsize=8)
    fig.suptitle("Cells per metacell by treatment")
    fig.tight_layout()
    savefig(fig, outdir / "mc_cells_per_mc.pdf")

    # -- Figure 2: metacells per cell type (bar chart) ---------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = [pal_celltype.get(str(ct), "#C8C8C8") for ct in mc_per_ct["celltype"]]
    ax.barh(mc_per_ct["celltype"], mc_per_ct["n_metacells"], color=colors)
    ax.set_xlabel("# metacells")
    ax.set_title("Metacells per cell type")
    ax.invert_yaxis()
    fig.tight_layout()
    savefig(fig, outdir / "mc_per_celltype.pdf")

    # -- Figure 3: time assignment per metacell (stacked bar, sorted by 12H frac) -
    ta_cols = ["12H", "24H", "36H", "Negative", "Not_Assigned"]
    ta_present = [c for c in ta_cols if c in mc_qc.columns]
    if ta_present:
        mc_time = mc_qc[["mc_id_num", "celltype"] + ta_present].copy()
        mc_time = mc_time.sort_values("12H" if "12H" in ta_present else ta_present[0],
                                      ascending=False)
        ta_colors = sns.color_palette("Set2", n_colors=len(ta_present))

        fig, ax = plt.subplots(figsize=(10, 4))
        bottom = np.zeros(len(mc_time))
        for ta, color in zip(ta_present, ta_colors):
            ax.bar(range(len(mc_time)), mc_time[ta].values, bottom=bottom,
                   color=color, label=ta, width=1.0)
            bottom += mc_time[ta].values
        ax.set_xlabel("Metacell (sorted by 12H fraction)")
        ax.set_ylabel("Fraction of cells")
        ax.set_title("Time assignment distribution per metacell")
        ax.legend(title="Time bin", bbox_to_anchor=(1.01, 1), loc="upper left")
        fig.tight_layout()
        savefig(fig, outdir / "mc_time_assignment.pdf")

    # -- Figure 4: 2D metacell map colored by dominant cell type -----------------
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, col, title, cmap_col in [
        (axes[0], "celltype",  "Dominant cell type",    None),
        (axes[1], "n_cells",   "Cells per metacell",    "viridis"),
    ]:
        if col == "celltype":
            for ct in mc_qc["celltype"].dropna().unique():
                sub = mc_qc[mc_qc["celltype"] == ct]
                ax.scatter(sub["mc_x"], sub["mc_y"], s=sub["n_cells"] * 2,
                           c=pal_celltype.get(str(ct), "#C8C8C8"),
                           label=str(ct), alpha=0.8, edgecolors="none")
            ax.legend(fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left",
                      markerscale=0.8)
        else:
            sc = ax.scatter(mc_qc["mc_x"], mc_qc["mc_y"], s=mc_qc["n_cells"] * 2,
                            c=mc_qc["n_cells"], cmap=cmap_col,
                            alpha=0.8, edgecolors="none")
            plt.colorbar(sc, ax=ax, label="Cells per metacell")
        ax.set_title(title)
        ax.set_xlabel("MC x"); ax.set_ylabel("MC y")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Metacell 2D map")
    fig.tight_layout()
    savefig(fig, outdir / "mc_2d_map.pdf")

    print(f"  Metacell QC figures written to {outdir}/")

print(f"Done. All outputs in {outdir}/")
