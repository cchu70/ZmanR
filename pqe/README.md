# pqe/ — Zman-seq Analysis

This directory contains the end-to-end analysis pipeline for the Zman-seq experiment, which profiles monocyte-to-macrophage (MoMac) differentiation trajectories in a GBM tumor model under aTrem2 antibody treatment vs. IgG isotype control.

---

## General Approach

The analysis follows four main phases:

1. **Data assembly and QC** — Parse the Zman-seq metadata and raw UMI count matrices, assemble an AnnData object, apply gene and cell quality filters.

2. **Metacell construction** — Compress single cells into metacells (tanaylab/metacells 0.9.x) to reduce noise and enable robust trajectory analysis. Done first on all cell types, then restricted to myeloid cells.

3. **cTET annotation and pseudotime benchmarking** — Compute the cumulative Time-labeling Enrichment Trajectory (cTET) metric from Zman-seq time-bin proportions, then benchmark five pseudotime tools (DPT, Monocle2, SCORPIUS, redPath, Palantir) across six gene sets and two UMI normalization schemes. Spearman correlations between gene expression and cTET are used to define informative gene sets.

4. **CosMx spatial integration** — Transfer the Zman-seq MoMac trajectory onto a human CosMx 1k spatial transcriptomics dataset via ortholog mapping, metacell construction on the spatial data, and SingleR-based cell-to-metacell assignment.

---

## Scripts and Notebooks

### Top-level numbered scripts

| File | Description | Key outputs |
|------|-------------|-------------|
| `00_test_example.R` | Tests the ZmanR R package on example data: GLM model evaluation, time assignment, FACS model QC plots. | `results/00/glm_fit.pdf` |
| `01_summary_aTrem2.R` | Loads `metadata.txt`, parses cell metadata, counts marker gene expression per cell, and generates summary plots of cell-type composition and marker expression by treatment (aTrem2 vs. IgG). | `results/01/` — CSVs and plots |
| `01_summary_aTrem2.py` | Python reproduction of `01_summary_aTrem2.R`. Loads `cells_obs.csv` + `marker_counts.csv` written by the R script, builds an AnnData (`cells_markers.h5ad`), and reproduces all summary tables and figures. | `results/01/py_version/` |
| `02_simulator.py` | Interactive Dash web app for visualizing Zman-seq labeling trajectories. Simulates cells along up to two parallel A→B→C state progressions, with configurable label times and population sizes. | Interactive app (no file output) |
| `03_Metadata_explore.ipynb` | Exploratory notebook: parses `metadata.txt`, visualizes cell spatial coordinates (`sc_x`, `sc_y`) by treatment, time assignment, and cell type; also parses a second metadata file (`metadata_q27.txt`) for colon/lung experiments. | Inline plots |
| `04_gen_adata.py` | Parses `metadata.txt` and loads per-batch raw UMI count matrices (`.txt.gz`) to build a combined AnnData. | `results/04/zmanseq.h5ad` — 43,392 cells × 52,634 genes |
| `04_qc_adata.ipynb` | QC notebook: queries Ensembl BioMart for unwanted RNA biotypes, applies a regex filter (mt, ribo, Ig, poorly-annotated genes), removes cells with < 300 UMIs, identifies HVGs. Documents the gene-filtering criteria used in downstream scripts. | Inline QC plots; defines filtering logic |
| `05_metacell.py` | Runs tanaylab/metacells divide-and-conquer metacell construction on all cells (target size 100). Computes a kNN adjacency matrix on log-normalized metacell expression. | `results/05/zmanseq_cells.h5ad`, `zmanseq_metacells.h5ad`, `zmanseq_metacell_adjacency.npz` |
| `05_metacell.ipynb` | Notebook companion to `05_metacell.py`; interactive exploration of metacell results. | Inline plots |
| `06_metacell_myeloid.py` | Same metacell construction as `05_metacell.py` but restricted to myeloid cell types (Monocyte, MoMac1/2, TAMs, DCs, etc.) from IgG and aTrem2 conditions only. | `results/06/zmanseq_myeloid_cells.h5ad`, `zmanseq_myeloid_metacells.h5ad`, `zmanseq_myeloid_metacell_adjacency.npz` |
| `06_metacell_myeloid.ipynb` | Notebook companion to `06_metacell_myeloid.py`. | Inline plots |
| `06_ctet_mos.py` | Annotates myeloid metacells (cell type by majority vote, spatial coordinates, treatment enrichment) and computes cTET from cumulative time-bin proportions. Builds four processed h5ad objects covering all-myeloid vs. Mo-type-only subsets and two UMI downsampling schemes (median single-cell UMI vs. mean metacell UMI). Adds `smooth_cTET` (kNN + cell-type anchor smoothing) to Mo-type objects. | `results/06/zmanseq_myeloid_metacells_annot_clean.h5ad` (all myeloid, sc-UMI), `..._mcumi.h5ad` (all myeloid, MC-UMI), `zmanseq_momac_metacells_annot_clean.h5ad` (Mo-only, sc-UMI), `..._mcumi.h5ad` (Mo-only, MC-UMI), `06_ctet_barplots.pdf` |
| `06_ctet_mos_spearman.py` | Computes Spearman correlation of each gene with cTET and smooth_cTET for the two Mo-only adata objects, stratified by treatment arm (aTrem2 vs. IgG) and gene set (HVG vs. full). Generates volcano plots, ordered heatmaps, line plots, and Venn diagrams comparing significant gene sets. | `results/06/ctet_mos/spearman/*.csv`, `results/06/ctet_mos/plots/*.pdf` |
| `07_pseudotime_compare.py` | Collects pseudotime outputs from all tool × gene-set × normalization combinations into a single multi-indexed TSV. Also reports gene-set sizes after adata filtering. | `results/07/pseudotime_momac_all.tsv`, `results/07/gene_set_sizes.tsv` |
| `07_pseudotime_compare.ipynb` | Notebook companion to `07_pseudotime_compare.py`; comparison plots and correlation analyses across pseudotime tools. | Inline plots |
| `10_cosmx_summary.R` | Summarizes the CosMx 1k/6k spatial transcriptomics dataset: cells per sample, cell-type composition, TSPS/simplicity score distributions, and spatial scatter plots. | `results/10/` — PDF summary plots |
| `10_cosmx_adata.R` | Converts `mData_CosMx1k_CosMx6k.RData` + `spatial_transcriptomics.zip` to per-sample AnnData h5ad files. Builds expression + coordinate matrices for human CosMx 1k samples (UM01–UM07). | `results/10/*.h5ad`, `mData_CosMx1k_CosMx6k.tsv.gz` |
| `10_cosmx_annotate_mouse.py` | Builds annotated mouse CosMx h5ad files (mHP and PDX tumor models) from the spatial transcriptomics zip, merging cell-type annotations from Table S6A (tumor cell state, TSPS, GPM/MTC/NEU/PPR probabilities). | `results/10/cosmx_mhp.h5ad`, `results/10/cosmx_pdx.h5ad` |
| `10_cosmx_explore.ipynb` | Exploratory notebook for the CosMx dataset; spatial and expression visualizations. | Inline plots |
| `11_map_orthologs.py` | Translates mouse ZmanSeq myeloid adata (single cells and MoMac metacells) to human gene symbols using mousipy, then subsets to genes in the CosMx 1k panel. | `results/11/zmanseq_myeloid_cells_human.h5ad`, `zmanseq_momac_metacells_annot_clean_human.h5ad`, `zmanseq_momac_metacells_annot_clean_mcumi_human.h5ad`, `human_to_mouse_orthologs.tsv` |
| `12_metacell_cosmx.py` | Metacell construction (target size 30, ~100 metacells) on human CosMx 1k MoMac cells (UM01 sample, ~3,042 cells). Applies human-specific gene filtering, annotates metacells, and saves summary plots. | `results/12/momac_um01_cells.h5ad`, `momac_um01_metacells.h5ad`, `12_*.pdf` |
| `12_metacell_cosmx_s10.py` | Same as `12_metacell_cosmx.py` but with finer resolution (target size 15, ~56 metacells). | `results/12/s10/momac_um01_cells_s10.h5ad`, `momac_um01_metacells_s10.h5ad` |
| `13_singleR_cosmx.py` | Assigns CosMx MoMac cells (UM01) to ZmanSeq MoMac metacell labels using SingleR, with IgG single cells as reference. Runs an embedded R script for SingleR classification and returns per-cell predicted labels and scores. | `results/13/singleR_cosmx_um01.tsv` |

---

### src/ — Pseudotime runner infrastructure

| File(s) | Description |
|---------|-------------|
| `gen_rds.R` | One-time helper: reads `zmanseq.h5ad` via zellkonverter and saves raw counts + metadata as an RDS. Used by R environments that cannot read h5ad directly. |
| `gen_momac_rds.R`, `gen_momac_rds.py` | Same as above for the MoMac metacell objects. |
| `make_rds_companions.R`, `make_rds_from_exports.R` | Helper scripts to create RDS companion files for exported h5ad objects. |
| `utils.R` | Shared R utilities for pseudotime scripts: `load_adata()` (h5ad or RDS), `load_hvg_genes()`, `load_zman_s3_genes()`, harmonized TSV output. |
| `utils_momac.R` | Extension of `utils.R` with `load_ctet_mos_genes()` — loads the Spearman-significant gene intersection from `06_ctet_mos_spearman.py` outputs. |
| `generate_run_scripts.py` | Meta-script that programmatically generates all 120 pseudotime runner scripts (5 tools × 6 gene sets × 2 normalizations: scumi and mcumi). |
| `run_dpt_{variant}.R` | Run DPT (diffusion pseudotime via destiny) on the MoMac metacell adata for a given gene set and normalization. Outputs `dpt_pseudotime.tsv`. |
| `run_monocle2_{variant}.R` | Run Monocle 2 DDRTree trajectory on the MoMac metacell adata. Outputs `monocle2_pseudotime.tsv`. |
| `run_scorpius_{variant}.R` | Run SCORPIUS principal curve on the MoMac metacell adata. Outputs `scorpius_pseudotime.tsv`. |
| `run_redpath_{variant}.R` | Run redPath (Spearman-distance minimum spanning tree) on the MoMac metacell adata. Outputs `redpath_pseudotime.tsv`. |
| `run_palantir_{variant}.py` | Run Palantir diffusion pseudotime on the MoMac metacell adata. Outputs `palantir_pseudotime.tsv` (includes entropy). |

Each `{variant}` encodes the normalization (`scumi` or `mcumi`) and gene set (`hvg`, `hvg_ctet`, `hvg_smooth_ctet`, `full_ctet`, `full_smooth_ctet`, `zman_s3`). All runner scripts write their output under `results/06/pseudotime_{variant}/`.

**Gene sets used for pseudotime:**

| Gene set | Source |
|----------|--------|
| `hvg` | Top 1,000 HVGs from the adata |
| `hvg_ctet` | HVG genes with Spearman p < 0.05 for cTET (intersection of aTrem2 and IgG arms) |
| `hvg_smooth_ctet` | Same but correlated with smooth_cTET |
| `full_ctet` | All genes with Spearman p < 0.05 for cTET |
| `full_smooth_ctet` | All genes with Spearman p < 0.05 for smooth_cTET |
| `zman_s3` | Published gene list from Zman-seq Table S3 (p < 0.05 in both arms) |

---

## Data Flow Summary

```
metadata.txt + UMI count matrices
         |
    04_gen_adata.py
         |
    zmanseq.h5ad  (43k cells)
         |
    05_metacell.py          06_metacell_myeloid.py
         |                           |
  all-cell metacells         myeloid metacells (72)
                                     |
                              06_ctet_mos.py
                                     |
                    4x processed h5ad (scumi/mcumi x all/mo-only)
                    + cTET, smooth_cTET annotations
                                     |
                    06_ctet_mos_spearman.py
                    (Spearman gene rankings)
                                     |
                    generate_run_scripts.py
                    (5 tools x 6 gene sets x 2 norms)
                                     |
                    07_pseudotime_compare.py
                    (combined pseudotime_momac_all.tsv)
                                     |
         +-----------CosMx integration-----------+
         |                                       |
  10_cosmx_adata.R             11_map_orthologs.py
  10_cosmx_annotate_mouse.py        |
         |                  human-translated adata
         |                          |
  human CosMx 1k           12_metacell_cosmx.py
  MoMac cells                       |
         |                  CosMx metacells
         +----------13_singleR_cosmx.py----------+
                    (cell → ZmanSeq metacell assignment)
```
