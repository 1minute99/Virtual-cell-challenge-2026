#!/usr/bin/env bash
# Usage: scripts/run_detached.sh NAME LOG_DIRECTORY /absolute/program [arguments...]
# Runs from the caller's working directory; survives editor/terminal closure.
set -euo pipefail
if (( $# < 3 )); then
  echo 'Usage: run_detached.sh NAME LOG_DIRECTORY /absolute/program [arguments...]' >&2
  exit 2
fi
name=$1
log_directory=$2
shift 2
[[ "$name" =~ ^[a-zA-Z0-9][a-zA-Z0-9_-]*$ ]] || { echo 'Invalid job name' >&2; exit 2; }
[[ "$1" = /* && -x "$1" ]] || { echo 'Use an absolute executable path (including the Python environment).' >&2; exit 2; }
unit="aivc-${name}.service"
if [[ $(systemctl --user show "$unit" -p LoadState --value) != not-found ]]; then
  echo "Job $unit already exists; inspect it before choosing a new name." >&2
  exit 2
fi
mkdir -p "$log_directory"
log_directory=$(realpath "$log_directory")
systemd-run --user --unit="$unit" --slice=aivc-jobs.slice \
  --working-directory="$PWD" \
  --property=Type=exec --property=RemainAfterExit=yes --property=Restart=no \
  --property=MemoryAccounting=yes --property=MemoryHigh=40G --property=MemoryMax=48G \
  --property=MemorySwapMax=2G --property=CPUQuota=800% --property=TasksMax=512 \
  --property=Nice=10 --property=StandardInput=null \
  --property="StandardOutput=append:$log_directory/$name.log" \
  --property="StandardError=append:$log_directory/$name.log" \
  --setenv=PYTHONUNBUFFERED=1 --setenv=OMP_NUM_THREADS=4 \
  --setenv=OPENBLAS_NUM_THREADS=4 --setenv=MKL_NUM_THREADS=4 \
  -- "$@"
printf 'Status: systemctl --user status %s\nLog: %s/%s.log\n' "$unit" "$log_directory" "$name"
