#!/usr/bin/env bash
# Task 10 (plan input/ideas/zero-shot-limits-context-complexity-plan.md): seed 1 + 2 replication of six arch-screen
# architectures (additive, moe, bilinear, full_hyper, token_transformer, cross_attn; selection recorded in
# results/rq2_arch/selection.json) on K562 + HepG2 through the seed-0 protocol (arch_screen.py training + harness,
# ext_preset_score.py official cell-eval2 --preset vcc2026 --set de.backend=pdex).
# Resource rules: waits until BOTH aivc-paper-t17-real and aivc-paper-t18-queue have exited; every build and scorer
# step holds results/.celleval.lock (one at a time); non-scorer steps run in a transient scope capped at
# MemoryHigh 8G / MemoryMax 10G (tree: python + harness worker), nice 10, OMP/MKL 4 threads. Idempotent: rerun resumes.
# Launch: cd /home/wonmin-kim/AIVC && scripts/run_detached.sh paper-t10-queue output/data-analysis/results/rq2_arch/logs \
#         /home/wonmin-kim/AIVC/output/data-analysis/rq2_arch_replicate.sh
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONUNBUFFERED=1
# 2026-11-06 execution fix: arch_screen.py -> b19_optimization_data.Data() loads the base table F.TABLE, whose y.npy /
# controls.npy were deleted by the 2026-10-11 cleanup (first run: FileNotFoundError 3 s into seed-0 training). The
# documented NaN stub (AGENTS.md; reports/final_2026-10-27/README.md) is used: the arch screen reads only
# features / feature_map / targets from the base table (the training/refit data come from ext_contexts table_v1).
export TV2_BASE_TABLE=/home/wonmin-kim/AIVC/vcc2026_data/reports/final_2026-10-27/base_table_stub
LOCK=/home/wonmin-kim/AIVC/output/data-analysis/results/.celleval.lock
D=/home/wonmin-kim/AIVC/output/data-analysis
A=$D/results/rq2_arch
CT=/home/wonmin-kim/miniconda3/envs/ct_transfer/bin/python      # torch env (as ext_chain.sh used for ext_film)
VP=/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python          # scorer env
SIX=additive,moe,bilinear,full_hyper,token_transformer,cross_attn
mkdir -p "$A/logs"
log() { echo "[$(date -Is)] $*"; }
busy() { local s; s=$(systemctl --user show -p SubState --value "$1" 2>/dev/null || true)
         [[ -n "$s" && "$s" != exited && "$s" != dead && "$s" != failed ]]; }
cap() { systemd-run --user --scope --quiet --slice=aivc-jobs.slice -p MemoryHigh=8G -p MemoryMax=10G \
        -p MemorySwapMax=1G /usr/bin/time -v nice -n 10 "$@"; }
# 2026-11-06 main-session exception (option a): the TRAINING step only gets MemoryHigh 12G / MemoryMax 14G. The unchanged
# seed-0 refit inside arch_screen (ext_film.fold -> LeanSourceBaseline on table_v1) needs ~9.1 GB anon and crawled at the
# 8G cap. Everything else stays at the 8 GB rule. Guard: start only with >= 20 GB available and no other AIVC job running
# a scorer or holding > 4 GB (checked under the lock, re-checked every 5 min).
capT() { systemd-run --user --scope --quiet --slice=aivc-jobs.slice -p MemoryHigh=18G -p MemoryMax=20G \
        -p MemorySwapMax=1G /usr/bin/time -v nice -n 10 "$@"; }
guard() { local av big sc u m
  while :; do
    av=$(free -g | awk '/^Mem:/{print $7}'); big=""
    for u in $(systemctl --user list-units --type=service --type=scope --state=running --no-legend --plain | awk '{print $1}'); do
      [[ "$(systemctl --user show -p Slice --value "$u")" == aivc-jobs.slice ]] || continue
      m=$(systemctl --user show -p MemoryCurrent --value "$u"); [[ "$m" =~ ^[0-9]+$ ]] && (( m > 4294967296 )) && big="$big $u"
    done
    sc=$(pgrep -f '[c]ell-eval2 run|[c]ell_eval2' | tr '\n' ' ' || true)
    if (( av >= 20 )) && [[ -z "$big" && -z "${sc// /}" ]]; then echo "[$(date -Is)] guard ok: avail=${av}G"; return 0; fi
    echo "[$(date -Is)] guard wait: avail=${av}G big=[$big] scorer=[$sc]"; sleep 300
  done; }

while busy aivc-paper-t17-real.service || busy aivc-paper-t18-queue.service; do sleep 120; done
log "T17 + T18 queue exited; starting Task 10"

train() {  # seed line dir archs
  log "train seed=$1 line=$2 -> $3 (waiting for lock)"
  mkdir -p "$3"
  flock "$LOCK" bash -c "$(declare -f capT guard); guard >&2; capT '$CT' '$D/rq2_arch_seed_screen.py' --seed $1 --line $2 --out '$3' --archs $4" >> "$3/screen.log" 2>&1
  grep -q '"event": "done"' "$3/screen.log" || { log "TRAIN FAILED $3"; exit 1; }
}
preset() {  # seed line dir arms   (scorer: lock, no extra cap -- the scorer is exempt from the 8 GB rule)
  log "preset seed=$1 line=$2 arms=$4 (waiting for lock)"
  flock "$LOCK" nice -n 10 "$VP" "$D/rq2_arch_preset.py" score --src "$3" --line $2 --seed $1 --arms "$4" >> "$3/preset.log" 2>&1
  log "preset done: receipts=$(find "$3/preset/runs" -name receipt.json | wc -l) scored=$(find "$3/preset/runs" -name scored.csv | wc -l)"
}
collect() { nice -n 10 "$VP" "$D/rq2_arch_seeds.py" >> "$A/logs/collect.log" 2>&1 || log "collector failed (see collect.log)"; }
arms() { local s="base_refit_v2"; for a in ${1//,/ }; do s="$s,${a}_correct"; done; echo "$s"; }

# 1. regenerate the frozen real panels (K562 file is missing on disk; HepG2 = determinism control), under the lock
for L in HepG2 K562; do
  log "realpanel $L (waiting for lock)"
  flock "$LOCK" bash -c "$(declare -f cap); cap '$VP' '$D/rq2_arch_preset.py' realpanel --line $L" >> "$A/logs/realpanel.log" 2>&1
done
grep '"event": "realpanel' "$A/logs/realpanel.log" | tail -2

# 2. seed-0 reproduction gate on K562 (base re-score must equal the stored 0.3480608966819578 to 1e-6)
R0=$A/s0repro_K562
train 0 K562 "$R0" additive,moe
preset 0 K562 "$R0" "$(arms additive,moe)"
collect
"$VP" - <<'PY' || { log "SEED-0 REPRODUCTION GATE FAILED -- queue stopped"; exit 1; }
import json, sys
s = json.load(open('/home/wonmin-kim/AIVC/output/data-analysis/results/rq2_arch_seeds_summary.json'))['seed0_reproduction']
print(json.dumps(s)); sys.exit(0 if s.get('base_refit_v2_seed0_composite', {}).get('pass_1e6') else 1)
PY

# 3. the 24 architecture runs + 4 paired baseline runs
for S in 1 2; do
  for L in K562 HepG2; do
    O=$A/s${S}_${L}
    train $S $L "$O" $SIX
    preset $S $L "$O" "$(arms $SIX)"
    collect
  done
done
cp -p "$D/results/rq2_arch_seeds_summary.json" "$A/rq2_arch_seeds_summary.json"   # copy at the path named in the brief
log "QUEUE_DONE"
