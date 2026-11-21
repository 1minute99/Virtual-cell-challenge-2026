#!/usr/bin/env bash
# Task 11 offline runner: one rq2_source_curve.py subcommand in a transient scope capped at MemoryHigh 7G / MemoryMax 8G,
# nice 10, OMP/MKL/OpenBLAS 4 threads (resource rules 2026-11-05). Usage: rq2_source_curve_run.sh SUBCMD [ARGS...]
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 PYTHONUNBUFFERED=1
cd /home/wonmin-kim/AIVC/output/data-analysis
exec systemd-run --user --scope --quiet --slice=aivc-jobs.slice -p MemoryHigh=7G -p MemoryMax=8G -p MemorySwapMax=1G \
  /usr/bin/time -v nice -n 10 /home/wonmin-kim/miniconda3/envs/vcc2026/bin/python -P rq2_source_curve.py "$@"
