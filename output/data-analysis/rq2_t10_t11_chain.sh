#!/usr/bin/env bash
# 2026-11-06 resume chain: Task 10 queue (rq2_arch_replicate.sh, resumes where it stopped) THEN Task 11 spot check
# (rq2_source_curve_spot.sh). Serial, so one queue at a time; both scripts take results/.celleval.lock per build/scorer step.
# T11 runs even if T10 fails (independent tasks); each keeps writing to its own log.
set -uo pipefail
D=/home/wonmin-kim/AIVC/output/data-analysis
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 PYTHONUNBUFFERED=1
echo "[$(date -Is)] chain start: T10"
/usr/bin/bash "$D/rq2_arch_replicate.sh" >> "$D/results/rq2_arch/logs/paper-t10-queue.log" 2>&1; r10=$?
echo "[$(date -Is)] T10 exit $r10; starting T11 spot"
/usr/bin/bash "$D/rq2_source_curve_spot.sh" >> "$D/results/logs/paper-t11-spot.log" 2>&1; r11=$?
echo "[$(date -Is)] T11 spot exit $r11; CHAIN_DONE"
exit $(( r10 != 0 || r11 != 0 ))
