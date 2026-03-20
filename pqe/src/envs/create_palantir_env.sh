#!/bin/bash
# Creates the palantir_env conda environment and installs Palantir from submodule.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PREFIX="/data/ClaudiaC/envs/palantir_env"

echo "=== Creating palantir_env ==="
conda create --prefix "$PREFIX" -y -c conda-forge python=3.10 pip

echo "=== Installing Palantir and dependencies from submodule ==="
"$PREFIX/bin/pip" install -e "$REPO_ROOT/pqe/src/palantir"

echo "=== Installing scanpy and matplotlib ==="
"$PREFIX/bin/pip" install scanpy matplotlib

echo "=== Done: palantir_env ==="
