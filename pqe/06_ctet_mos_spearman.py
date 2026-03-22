#!/usr/bin/env python
"""
06_ctet_mos_spearman.py

Spearman correlation of each gene with cTET / smooth_cTET for the two
mo-only madata objects (sc-UMI and MC-UMI downsampling), split by treatment arm.

Inputs (from 06_ctet_mos.py):
  results/06/zmanseq_momac_metacells_annot_clean.h5ad      (sc UMI)
  results/06/zmanseq_momac_metacells_annot_clean_mcumi.h5ad (MC UMI)

Outputs — all under results/06/ctet_mos/:
  spearman/
    {label}_{arm}_{geneset}_{timecol}_spearman.csv
        label   : scumi | mcumi
        arm     : atrem | igg
        geneset : hvg   | full
        timecol : ctet  | smooth_ctet

  plots/
    {label}_volcano_{timecol}.pdf       — 4-panel volcano per adata × time_col
    {label}_{arm}_{geneset}_{timecol}_heatmap.pdf
    {label}_{arm}_{geneset}_{timecol}_lineplot.pdf
    {label}_venn.pdf                    — smooth vs naive, arm vs arm, arm vs S3
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from scipy.stats import spearmanr
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from matplotlib_venn import venn2, venn3

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR = Path("results/06")
OUT_DIR  = DATA_DIR / "ctet_mos"
SPEAR_DIR = OUT_DIR / "spearman"
PLOT_DIR  = OUT_DIR / "plots"
for d in [SPEAR_DIR, PLOT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

S3_ATREM_FN = "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_aTREM2_Time.csv"
S3_IGG_FN   = "/mnt/thechenlab/ClaudiaC/zmanseq/mmc2.Table_S3_Isotype_Control_Time.csv"

# ── Constants ────────────────────────────────────────────────────────────────
ATREM_TYPES = ["Monocytes", "MoMac1", "Acp5_TAM"]
IGG_TYPES   = ["Monocytes", "MoMac2", "Gpnmb_TAM", "Arg1_TAM"]

MYELOID_CELL_TYPE_COLORS = {
    'Arg1_TAM':   "mediumorchid",
    'Acp5_TAM':   "slateblue",
    'MoMac1':     "tab:blue",
    'MoMac2':     "tab:red",
    'Monocytes':  "tab:orange",
    'cDC2':       "goldenrod",
    'MonDC':      "gold",
    'Gpnmb_TAM':  "lightblue",
    'cDC1':       "darkkhaki",
    'MigDC':      "tab:brown",
    'transitory': "teal",
    'pDC':        "tab:green",
    'nan':        "gray",
}

ADATAS = [
    ("scumi", DATA_DIR / "zmanseq_momac_metacells_annot_clean.h5ad"),
    ("mcumi", DATA_DIR / "zmanseq_momac_metacells_annot_clean_mcumi.h5ad"),
]

P_THRESH = 0.05

# ── Core functions ────────────────────────────────────────────────────────────

def get_spearman_df(freq_sort_df, madata, time_col="cTET"):
    """Spearman r + p-value between each gene and time_col for cells in freq_sort_df."""
    result = pd.DataFrame(index=freq_sort_df.columns, columns=["spearman_r", "p_value"])
    cells = freq_sort_df.index
    time_vals = madata.obs.loc[cells, time_col].astype(float).values
    for gene in freq_sort_df.columns:
        s = spearmanr(time_vals, freq_sort_df.loc[cells, gene].values)
        result.loc[gene, "spearman_r"] = s.statistic
        result.loc[gene, "p_value"]    = s.pvalue
    result = result.dropna()
    result["spearman_r"] = result["spearman_r"].astype(float)
    result["p_value"]    = result["p_value"].astype(float)
    result["-log10(p_value)"] = -np.log10(result["p_value"])
    return result


def build_freq_sort_df(madata, cells, genes=None):
    """Expression matrix (cells × genes) for given cells; genes=None → all genes."""
    if genes is None:
        genes = madata.var_names
    return pd.DataFrame(
        madata[cells, genes].X.toarray(),
        index=cells,
        columns=genes,
    )


def plot_hvg_spearman_heatmap(freq_sort_df, spearman_df, madata, time_col="cTET", ax_title=""):
    """Clustermap of significant genes sorted by time_col, with cell-type legend and cTET colorbar."""
    if spearman_df.empty:
        print(f"  [skip heatmap] no significant genes for {ax_title}")
        return None
    norm = Normalize(
        vmin=madata.obs.loc[freq_sort_df.index, time_col].min(),
        vmax=madata.obs.loc[freq_sort_df.index, time_col].max(),
    )
    sorted_mcs = madata.obs.loc[freq_sort_df.index].sort_values(time_col).index
    sig_genes  = [g for g in spearman_df.index if g in freq_sort_df.columns]
    g = sns.clustermap(
        freq_sort_df.loc[sorted_mcs, sig_genes].apply(zscore).apply(gaussian_filter1d, sigma=2).T,
        col_cluster=False,
        col_colors=[
            madata.obs.loc[sorted_mcs, "metacell_cell_type"].map(MYELOID_CELL_TYPE_COLORS),
            madata.obs.loc[sorted_mcs, time_col].apply(lambda x: plt.cm.viridis_r(norm(x))),
        ],
        cmap="rocket",
    )
    g.fig.suptitle(ax_title, y=1.02)

    # ── Cell type legend in the col_dendrogram space (top-left) ──────────────
    ct_present = madata.obs.loc[freq_sort_df.index, "metacell_cell_type"].unique()
    patches = [mpatches.Patch(color=MYELOID_CELL_TYPE_COLORS.get(ct, "gray"), label=ct)
               for ct in sorted(ct_present)]
    g.ax_col_dendrogram.legend(
        handles=patches, loc="upper left", ncol=min(len(patches), 3),
        bbox_to_anchor=(0.0, 1.0), frameon=True, fontsize=7, title="Cell type",
        title_fontsize=7,
    )

    # ── cTET colorbar: placed below the existing rocket cbar ─────────────────
    existing_pos = g.ax_cbar.get_position()
    ctet_cbar_ax = g.fig.add_axes([
        existing_pos.x0,
        max(0.02, existing_pos.y0 - existing_pos.height - 0.06),
        existing_pos.width,
        existing_pos.height * 0.8,
    ])
    sm_ctet = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
    sm_ctet.set_array([])
    g.fig.colorbar(sm_ctet, cax=ctet_cbar_ax, label=time_col)

    return g


def plot_hvg_spearman_lineplot(freq_sort_df, spearman_df, madata, time_col="cTET", ax_title=""):
    """Line plot of z-scored, smoothed expression vs time_col."""
    if spearman_df.empty:
        print(f"  [skip lineplot] no significant genes for {ax_title}")
        return None
    sorted_mcs  = madata.obs.loc[freq_sort_df.index].sort_values(time_col).index
    up_genes    = [g for g in spearman_df.query("spearman_r > 0").index if g in freq_sort_df.columns]
    down_genes  = [g for g in spearman_df.query("spearman_r < 0").index if g in freq_sort_df.columns]
    norm        = Normalize(vmin=-1, vmax=1)
    fig, ax     = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)
    fig.suptitle(ax_title)
    for i, (genes, title) in enumerate([(up_genes, "Up-regulated"), (down_genes, "Down-regulated")]):
        for g in genes:
            sns.lineplot(
                x=madata.obs.loc[sorted_mcs, time_col],
                y=gaussian_filter1d(zscore(freq_sort_df.loc[sorted_mcs, g]), sigma=2),
                color=plt.cm.RdBu_r(norm(spearman_df.loc[g, "spearman_r"])),
                legend=False, ax=ax[i], alpha=0.5,
            )
        ax[i].set_xlabel(time_col)
        ax[i].set_ylabel("Z-scored expression")
        ax[i].set_title(f"{title} (n={len(genes)})")
    # Dedicated colorbar axis — avoids shrinking subplots
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.87, 0.15, 0.025, 0.7])
    sm  = plt.cm.ScalarMappable(cmap=plt.cm.RdBu, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label="Spearman r")
    return fig


def plot_volcano(spearman_df, title="", ax=None):
    """Volcano: spearman_r vs -log10(p_value)."""
    own_ax = ax is None
    if own_ax:
        fig, ax = plt.subplots(figsize=(5, 4))
    n_sig = (spearman_df["p_value"] < P_THRESH).sum()
    ax.scatter(
        spearman_df["spearman_r"], spearman_df["-log10(p_value)"],
        s=8, alpha=0.6, c=spearman_df["spearman_r"],
        cmap="RdBu_r", vmin=-1, vmax=1,
    )
    ax.axhline(-np.log10(P_THRESH), color="red", linestyle="--", lw=0.8)
    ax.set_xlabel("Spearman r")
    ax.set_ylabel("-log10(p_value)")
    ax.set_title(f"{title}\n(n_sig={n_sig})")
    if own_ax:
        plt.tight_layout()
        return fig


# ── Load S3 reference ────────────────────────────────────────────────────────
print("Loading S3 reference tables ...")
s3_atrem = pd.read_csv(S3_ATREM_FN)
s3_igg   = pd.read_csv(S3_IGG_FN)
s3_atrem["-log10(p_value)"] = -np.log10(s3_atrem["Pvalue"])
s3_igg["-log10(p_value)"]   = -np.log10(s3_igg["Pvalue"])
s3_atrem_sig = set(s3_atrem.query("Pvalue < @P_THRESH")["Gene"])
s3_igg_sig   = set(s3_igg.query("Pvalue < @P_THRESH")["Gene"])

# ── Main loop ────────────────────────────────────────────────────────────────
all_results = {}   # {(label, arm, geneset, timecol): spearman_df}

for label, adata_path in ADATAS:
    print(f"\n{'='*60}")
    print(f"Processing {label}: {adata_path.name}")
    madata = sc.read_h5ad(adata_path)
    print(f"  Shape: {madata.shape}")

    # ── HVG selection ────────────────────────────────────────────────────────
    madata_tmp = madata.copy()
    sc.pp.highly_variable_genes(madata_tmp, n_top_genes=1000)
    hvg_genes = madata_tmp.var_names[madata_tmp.var["highly_variable"]].tolist()
    print(f"  HVGs: {len(hvg_genes)}")
    del madata_tmp

    # ── Per-arm cell lists ────────────────────────────────────────────────────
    atrem_cells = madata.obs.query("metacell_cell_type in @ATREM_TYPES").index.tolist()
    igg_cells   = madata.obs.query("metacell_cell_type in @IGG_TYPES").index.tolist()
    print(f"  aTrem2 cells: {len(atrem_cells)},  IgG cells: {len(igg_cells)}")

    # ── Build expression matrices ─────────────────────────────────────────────
    atrem_hvg_df  = build_freq_sort_df(madata, atrem_cells, hvg_genes)
    igg_hvg_df    = build_freq_sort_df(madata, igg_cells,   hvg_genes)
    atrem_full_df = build_freq_sort_df(madata, atrem_cells)
    igg_full_df   = build_freq_sort_df(madata, igg_cells)

    freq_dfs = {
        ("atrem", "hvg"):  atrem_hvg_df,
        ("atrem", "full"): atrem_full_df,
        ("igg",   "hvg"):  igg_hvg_df,
        ("igg",   "full"): igg_full_df,
    }

    # ── Spearman correlations ─────────────────────────────────────────────────
    for (arm, geneset), freq_df in freq_dfs.items():
        for timecol, col in [("ctet", "cTET"), ("smooth_ctet", "smooth_cTET")]:
            key = (label, arm, geneset, timecol)
            print(f"  Spearman {label} {arm} {geneset} {timecol} ...", end=" ")
            sp_df = get_spearman_df(freq_df, madata, time_col=col)
            n_sig = (sp_df["p_value"] < P_THRESH).sum()
            print(f"n_sig={n_sig}")
            all_results[key] = sp_df
            out_fn = SPEAR_DIR / f"{label}_{arm}_{geneset}_{timecol}_spearman.csv"
            sp_df.to_csv(out_fn)

    # ── Volcano plots ─────────────────────────────────────────────────────────
    for timecol, col in [("ctet", "cTET"), ("smooth_ctet", "smooth_cTET")]:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"{label} — volcano plots ({col})", fontsize=13)
        panel_order = [
            ("atrem", "hvg",  axes[0, 0]),
            ("igg",   "hvg",  axes[0, 1]),
            ("atrem", "full", axes[1, 0]),
            ("igg",   "full", axes[1, 1]),
        ]
        for arm, geneset, ax in panel_order:
            sp_df = all_results[(label, arm, geneset, timecol)]
            plot_volcano(sp_df, title=f"{arm} {geneset}", ax=ax)
        plt.tight_layout()
        out_fn = PLOT_DIR / f"{label}_volcano_{timecol}.pdf"
        fig.savefig(out_fn, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {out_fn.name}")

    # ── Heatmaps + line plots ─────────────────────────────────────────────────
    for (arm, geneset), freq_df in freq_dfs.items():
        for timecol, col in [("ctet", "cTET"), ("smooth_ctet", "smooth_cTET")]:
            sp_df  = all_results[(label, arm, geneset, timecol)]
            sp_sig = sp_df.query("p_value < @P_THRESH")
            tag    = f"{label}_{arm}_{geneset}_{timecol}"

            # Heatmap
            g = plot_hvg_spearman_heatmap(freq_df, sp_sig, madata, time_col=col,
                                           ax_title=f"{label} {arm} {geneset} {col}")
            if g is not None:
                out_fn = PLOT_DIR / f"{tag}_heatmap.pdf"
                g.fig.savefig(out_fn, bbox_inches="tight")
                plt.close(g.fig)
                print(f"  Saved {out_fn.name}")

            # Line plot
            fig = plot_hvg_spearman_lineplot(freq_df, sp_sig, madata, time_col=col,
                                              ax_title=f"{label} {arm} {geneset} {col}")
            if fig is not None:
                out_fn = PLOT_DIR / f"{tag}_lineplot.pdf"
                fig.savefig(out_fn, bbox_inches="tight")
                plt.close(fig)
                print(f"  Saved {out_fn.name}")

    # ── Venn diagrams ─────────────────────────────────────────────────────────
    with PdfPages(PLOT_DIR / f"{label}_venn.pdf") as pdf:

        # 1. cTET vs smooth_cTET, per arm × gene_set
        for arm in ["atrem", "igg"]:
            for geneset in ["hvg", "full"]:
                naive_sig   = set(all_results[(label, arm, geneset, "ctet")].query("p_value < @P_THRESH").index)
                smooth_sig  = set(all_results[(label, arm, geneset, "smooth_ctet")].query("p_value < @P_THRESH").index)
                fig, ax = plt.subplots(figsize=(5, 4))
                venn2([naive_sig, smooth_sig], set_labels=("cTET", "smooth_cTET"), ax=ax)
                ax.set_title(f"{label} {arm} {geneset} — cTET vs smooth_cTET")
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 2. aTrem2 vs IgG, per gene_set × time_col
        for geneset in ["hvg", "full"]:
            for timecol in ["ctet", "smooth_ctet"]:
                atrem_sig = set(all_results[(label, "atrem", geneset, timecol)].query("p_value < @P_THRESH").index)
                igg_sig   = set(all_results[(label, "igg",   geneset, timecol)].query("p_value < @P_THRESH").index)
                fig, ax = plt.subplots(figsize=(5, 4))
                venn2([atrem_sig, igg_sig], set_labels=("aTrem2", "IgG"), ax=ax)
                ax.set_title(f"{label} {geneset} {timecol} — aTrem2 vs IgG")
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 3. Full spearman vs S3 reference, per arm × time_col
        for arm, s3_sig in [("atrem", s3_atrem_sig), ("igg", s3_igg_sig)]:
            for timecol in ["ctet", "smooth_ctet"]:
                our_sig = set(all_results[(label, arm, "full", timecol)].query("p_value < @P_THRESH").index)
                fig, ax = plt.subplots(figsize=(5, 4))
                venn2([our_sig, s3_sig], set_labels=(f"ours ({timecol})", "S3"), ax=ax)
                ax.set_title(f"{label} {arm} full — vs S3")
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 4. aTrem2 vs IgG vs S3 (full, smooth_ctet) — 3-way
        for timecol in ["ctet", "smooth_ctet"]:
            atrem_sig = set(all_results[(label, "atrem", "full", timecol)].query("p_value < @P_THRESH").index)
            igg_sig   = set(all_results[(label, "igg",   "full", timecol)].query("p_value < @P_THRESH").index)
            # aTrem2 vs IgG arms vs S3 aTrem2
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            venn3([atrem_sig, igg_sig, s3_atrem_sig],
                  set_labels=("aTrem2", "IgG", "S3 aTrem2"), ax=axes[0])
            axes[0].set_title(f"{label} full {timecol} — arms vs S3 aTrem2")
            venn3([atrem_sig, igg_sig, s3_igg_sig],
                  set_labels=("aTrem2", "IgG", "S3 IgG"), ax=axes[1])
            axes[1].set_title(f"{label} full {timecol} — arms vs S3 IgG")
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"  Saved {label}_venn.pdf")

print("\nDone.")
print(f"Spearman tables: {SPEAR_DIR}")
print(f"Plots:           {PLOT_DIR}")
