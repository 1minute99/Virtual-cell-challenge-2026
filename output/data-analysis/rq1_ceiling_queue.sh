#!/usr/bin/env bash
# Task 13 real-scorer queue (H1, PDS): waits for the T18, T10 and T17-rest queues to exit, then builds the two
# packages and scores them one after the other under the shared cell-eval2 lock, then collects into rq1_h1.json.
set -uo pipefail
LOCK=/home/wonmin-kim/AIVC/output/data-analysis/results/.celleval.lock
PY=/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python
S=/home/wonmin-kim/AIVC/output/data-analysis/rq1_ceiling.py
cd /home/wonmin-kim/AIVC/output/data-analysis
for u in aivc-paper-t18-queue aivc-paper-t10-queue aivc-aivc-paper-t17-rest; do
  while s=$(systemctl --user show -p SubState --value "$u.service"); [ "$s" = running ] || [ "$s" = start ] || [ "$s" = start-pre ]; do sleep 120; done
  echo "[$(date -Is)] $u exited (SubState=$s)"
done
echo "[$(date -Is)] build (under lock)"
/usr/bin/flock "$LOCK" /usr/bin/nice -n 10 "$PY" -P "$S" build || { echo "build failed"; exit 1; }
for a in oracle pipeline; do
  echo "[$(date -Is)] waiting for lock: $a"
  /usr/bin/flock "$LOCK" /usr/bin/nice -n 10 "$PY" -P "$S" real "$a" || echo "real $a FAILED"
  "$PY" -P "$S" collect
done
echo "[$(date -Is)] QUEUE_DONE"
