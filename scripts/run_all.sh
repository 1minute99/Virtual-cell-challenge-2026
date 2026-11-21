#!/usr/bin/env bash
# End-to-end reproduction. Env: ~/miniconda3/envs/AIVC (python 3.12).
set -euo pipefail
PY=~/miniconda3/envs/AIVC/bin/python
cd "$(dirname "$0")/.."

# Task 1 - coverage counts (gates everything)
$PY scripts/task1_coverage.py

# One-off cache: replogle_concat.h5ad -> per-cell-line .npy
[ -f work/cellline_cache/X_k562.npy ] || $PY scripts/build_cache.py

# Task 2 - VCC predictions (submission shape; NOT locally scorable)
$PY scripts/task2_vcc_predict.py

# Task 2 - surrogate scored with cell-eval (ground truth exists here)
for t in rpe1 jurkat hepg2; do
  $PY -u scripts/task2_baselines.py --source K562_gwps --target "$t" --seed 0
done
for t in jurkat hepg2 k562; do
  $PY -u scripts/task2_baselines.py --source rpe1 --target "$t" --seed 0
done
# target-gene inclusion effect
$PY -u scripts/task2_baselines.py --source K562_gwps --target hepg2 --seed 0 --drop-target-gene
$PY scripts/task2_aggregate.py

# Task 3 - SKIPPED: STATE ST pert vocabulary ∩ VCC targets = 0 (see results/task1_coverage.csv)

# Task 4 - representations vs transferability
$PY scripts/task4_transfer.py
$PY scripts/task4_embed.py --rep pca
$PY scripts/task4_embed.py --rep sclong
$PY scripts/task4_embed.py --rep state_se
$PY scripts/task4_correlate.py
