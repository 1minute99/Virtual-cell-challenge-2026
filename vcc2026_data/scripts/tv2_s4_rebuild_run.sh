#!/usr/bin/env bash
# Usage: tv2_s4_rebuild_run.sh OUT_DIR CONTROLS_DIR CTX_LIST
set -u
export TV2_S4_OUT=$1 VCC_OFFICIAL_DIR=$2 TV2_CTX=$3
cd /home/wonmin-kim/AIVC/vcc2026_data
$HOME/miniconda3/envs/ct_transfer/bin/python scripts/tv2_s4_rebuild.py --full
echo S4_REBUILD_DONE
