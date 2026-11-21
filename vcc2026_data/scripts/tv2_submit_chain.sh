#!/usr/bin/env bash
# transfer_v2 official build chain: build -> vcc prep -> sha. Stops BEFORE upload (READY marker).
# Usage (via run_detached): tv2_submit_chain.sh SPEC_JSON TAG
set -u
R=/home/wonmin-kim/AIVC/vcc2026_data
O=$R/reports/submission_2026-09-23
SPEC=$1; TAG=$2
VCCPY=$HOME/miniconda3/envs/vcc2026/bin/python
VCC=$HOME/miniconda3/envs/vcc2026/bin/vcc
OFF=$("$VCCPY" -c "import json,sys; print(json.load(open(sys.argv[1])).get('controls_dir','/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls'))" "$SPEC")
mkdir -p "$O"; cd "$R" || exit 1
step () { echo "=== $1 $(date -u +%FT%TZ) ==="; }
step build
"$VCCPY" -P scripts/tv2_submit_build.py "$SPEC" "$O/$TAG.pred.h5ad" > "$O/$TAG.build.log" 2>&1 || { echo "BUILD FAILED"; tail -20 "$O/$TAG.build.log"; exit 1; }
tail -4 "$O/$TAG.build.log"
step prep
"$VCC" prep "$O/$TAG.pred.h5ad" -g "$OFF/gene_names.csv" --perts "$OFF/pert_counts.csv" -o "$O/$TAG.vcc" > "$O/$TAG.prep.log" 2>&1 || { echo "PREP FAILED"; tail -20 "$O/$TAG.prep.log"; exit 1; }
cat "$O/$TAG.prep.log" | tail -8
step sha
sha256sum "$O/$TAG.vcc" | tee "$O/$TAG.sha256"
ls -la "$O/$TAG.vcc"
echo "READY_FOR_UPLOAD $(date -u +%FT%TZ)"
