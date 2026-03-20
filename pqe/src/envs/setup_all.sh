#!/bin/bash
# Run all environment setup scripts sequentially.
set -e

DIR="$(cd "$(dirname "$0")" && pwd)"

bash "$DIR/create_palantir_env.sh"   # Python — fastest, run first
bash "$DIR/create_dpt_env.sh"
bash "$DIR/create_scorpius_env.sh"
bash "$DIR/create_monocle2_env.sh"
bash "$DIR/create_redpath_env.sh"

echo "=== All environments created ==="
