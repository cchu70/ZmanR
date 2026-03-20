#!/bin/bash
# Creates the monocle2_env conda environment and installs Monocle2 from submodule.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PREFIX="/data/ClaudiaC/envs/monocle2_env"

echo "=== Creating monocle2_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge -c bioconda \
    r-base=4.1 r-biocmanager r-devtools \
    r-ggplot2 r-matrix r-rcpp r-irlba r-mass r-matrixstats \
    r-vgam r-ddrtree r-fastica r-combinat r-qlcmatrix \
    r-pheatmap r-slam r-viridis r-proxy r-rann \
    r-reshape2 r-dplyr r-tibble r-stringr \
    bioconductor-zellkonverter bioconductor-biocgenerics \
    bioconductor-biobase bioconductor-limma

echo "=== Installing HSMMSingleCell via BiocManager ==="
conda run --prefix "$PREFIX" Rscript -e "
  BiocManager::install('HSMMSingleCell', ask = FALSE, update = FALSE)
"

echo "=== Installing Monocle2 from submodule ==="
conda run --prefix "$PREFIX" Rscript -e "
  devtools::install_local('$REPO_ROOT/pqe/src/monocle2', upgrade = 'never')
"

echo "=== Done: monocle2_env ==="
