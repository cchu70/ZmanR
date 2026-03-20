#!/bin/bash
# Creates the dpt_env conda environment and installs destiny from submodule.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PREFIX="/data/ClaudiaC/envs/dpt_env"

echo "=== Creating dpt_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge -c bioconda \
    r-base=4.3 r-biocmanager r-devtools \
    r-ggplot2 r-matrix r-rcpp r-rcppeigen r-irlba \
    r-rspectra r-vim r-proxy r-scales r-ggthemes r-tidyr r-tidyselect \
    r-scatterplot3d r-smoother r-rcpphnsw \
    "r-ggplot.multistats" "r-knn.covertree" \
    bioconductor-zellkonverter bioconductor-biocgenerics \
    bioconductor-biobase bioconductor-summarizedexperiment \
    bioconductor-singlecellexperiment bioconductor-pcamethods

echo "=== Installing destiny from submodule ==="
conda run --prefix "$PREFIX" Rscript -e "
  devtools::install_local('$REPO_ROOT/pqe/src/destiny', upgrade = 'never')
"

echo "=== Done: dpt_env ==="
