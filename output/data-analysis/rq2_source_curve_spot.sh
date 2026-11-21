#!/usr/bin/env bash
# Task 11 exactness spot check (real scorer, H1, k = 1, ceil(K/2), K). Waits until the T18, T10 and T17-rest queues and the
# H1 curve unit have exited; package preparation/build and every scorer run hold results/.celleval.lock (one at a time).
set -euo pipefail
LOCK=/home/wonmin-kim/AIVC/output/data-analysis/results/.celleval.lock
D=/home/wonmin-kim/AIVC/output/data-analysis
PY=/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 PYTHONUNBUFFERED=1
busy() { local s; s=$(systemctl --user show -p SubState --value "$1" 2>/dev/null || true); [[ "$s" == running || "$s" == start* ]]; }
for u in aivc-paper-t18-queue aivc-paper-t10-queue aivc-aivc-paper-t17-rest aivc-paper-t11-curve-h1; do
  while busy "$u.service"; do sleep 120; done
done
echo "[$(date -Is)] queues exited; spot prep + build under lock"
flock "$LOCK" "$D/rq2_source_curve_run.sh" spot_prep
flock "$LOCK" "$D/rq2_source_curve_run.sh" spot_build
for i in 1 2 3; do
  echo "[$(date -Is)] waiting for lock: scorer run $i"
  cd "$D"; flock "$LOCK" nice -n 10 "$PY" -P "$D/rq2_source_curve.py" spot_real
  "$D/rq2_source_curve_run.sh" collect || echo "collect failed (non-fatal)"
done
echo "[$(date -Is)] SPOT_QUEUE_DONE"
