#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p .jupyter/config .jupyter/data .jupyter/runtime .jupyter/ipython .jupyter/cache

export JUPYTER_CONFIG_DIR="$PWD/.jupyter/config"
export JUPYTER_DATA_DIR="$PWD/.jupyter/data"
export JUPYTER_RUNTIME_DIR="$PWD/.jupyter/runtime"
export IPYTHONDIR="$PWD/.jupyter/ipython"
export XDG_CACHE_HOME="$PWD/.jupyter/cache"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-comparison-cache}"

NOTEBOOK="${1:-compare_intensity_profiles.ipynb}"

exec ../.venv/bin/jupyter lab --no-browser --ip=127.0.0.1 --port=8888 "$NOTEBOOK"
