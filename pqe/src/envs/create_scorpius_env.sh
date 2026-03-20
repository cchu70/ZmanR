#!/bin/bash
# Creates the scorpius_env conda environment and installs SCORPIUS from submodule.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PREFIX="/data/ClaudiaC/envs/scorpius_env"

echo "=== Creating scorpius_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge -c bioconda \
    r-base=4.3 r-biocmanager r-devtools \
    r-ggplot2 r-dplyr r-tidyr r-mass r-matrix r-mclust \
    r-rann r-reshape2 r-pheatmap r-ranger r-rcpp \
    r-pbapply r-purrr r-rcolorbrewer r-tsp r-princurve \
    r-igraph \
    bioconductor-zellkonverter

echo "=== Installing dynverse packages ==="
conda run --prefix "$PREFIX" Rscript -e "
  options(Ncpus = 4)
  remotes::install_github('dynverse/dynutils', upgrade = 'never')
  remotes::install_github('dynverse/lmds',    upgrade = 'never')
  remotes::install_github('dynverse/dynwrap',  upgrade = 'never')
"

echo "=== Installing SCORPIUS from submodule ==="
conda run --prefix "$PREFIX" Rscript -e "
  devtools::install_local('$REPO_ROOT/pqe/src/SCORPIUS', upgrade = 'never')
"

echo "=== Done: scorpius_env ==="
