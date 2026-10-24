#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2 MPLBACKEND=Agg
python scripts/download_data.py
mkdir -p results
/usr/bin/time -v -o results/runtime.txt python -m pbmc_annotation.pipeline
python scripts/build_report.py
python scripts/verify_results.py
python -m pytest -q
