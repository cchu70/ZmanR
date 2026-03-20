#!/bin/bash
# Creates the redpath_env conda environment and installs redPATH from submodule.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PREFIX="/data/ClaudiaC/envs/redpath_env"

echo "=== Creating redpath_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge -c bioconda \
    r-base=4.3 r-biocmanager r-devtools \
    r-ggplot2 r-ggalt r-mass r-mclust r-rcpp r-rcpparmadillo r-rcppalgos \
    r-rcppparallel \
    r-car r-doparallel r-dplyr r-energy r-foreach r-combinat r-gmp \
    bioconductor-zellkonverter bioconductor-scater \
    bioconductor-gosummaries

echo "=== Installing redPATH from submodule ==="
conda run --prefix "$PREFIX" Rscript -e "
  devtools::install_local('$REPO_ROOT/pqe/src/redPATH', upgrade = 'never')
"

echo "=== Done: redpath_env ==="
