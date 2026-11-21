#!/usr/bin/env bash
# Task 18 real-scorer queue. Waits for unit aivc-paper-t17-real to exit, builds the 7 packages under the lock, then
# runs each H1 arm under the shared cell-eval2 lock (one scorer at a time, AGENTS.md), verdict-relevant arms first.
set -euo pipefail
LOCK=/home/wonmin-kim/AIVC/output/data-analysis/results/.celleval.lock
PY=/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python
S=/home/wonmin-kim/AIVC/output/data-analysis/rq4_sameline.py
busy() { local s; s=$(systemctl --user show -p ActiveState --value "$1" 2>/dev/null || true); [[ "$s" == active || "$s" == activating || "$s" == deactivating ]] && [[ "$(systemctl --user show -p SubState --value "$1")" != exited ]]; }
while busy aivc-paper-t17-real.service; do sleep 120; done
# packages are built HERE, under the scorer lock, so the build never overlaps a scorer (oomd killed a concurrent build)
echo "[$(date -Is)] build (under lock)"
flock "$LOCK" "$PY" -P "$S" build
for a in w0 w0.25 w0.5 w0.75 w1 noM2 noM3; do [[ -f /home/wonmin-kim/AIVC/output/data-analysis/results/rq4/pkg/$a/meta.json ]] || { echo "package $a missing (build failed?)"; exit 1; }; done
for a in w0 w1 noM2 noM3 w0.5 w0.25 w0.75; do
  echo "[$(date -Is)] waiting for lock: $a"
  flock "$LOCK" "$PY" -P "$S" real "$a"
  "$PY" -P "$S" verdict
done
echo "[$(date -Is)] QUEUE_DONE"
