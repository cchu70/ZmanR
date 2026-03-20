#!/bin/bash
# Creates the metacells_env conda environment and installs metacells via pip.
set -e

PREFIX="/data/ClaudiaC/envs/metacells_env"

echo "=== Creating metacells_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge \
    python=3.10 pip \
    cmake make cxx-compiler

echo "=== Installing metacells ==="
"$PREFIX/bin/pip" install metacells

echo "=== Installing analysis dependencies ==="
"$PREFIX/bin/pip" install scanpy anndata pybiomart matplotlib seaborn

echo "=== Done: metacells_env at $PREFIX ==="
