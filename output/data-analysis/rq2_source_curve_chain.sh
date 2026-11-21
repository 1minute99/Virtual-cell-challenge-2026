#!/usr/bin/env bash
# Task 11: after the H1 curve unit exits, run the remaining six panels serially (one capped process at a time), then collect.
set -euo pipefail
while [[ "$(systemctl --user show -p SubState --value aivc-paper-t11-curve-h1.service)" == running ]]; do sleep 60; done
D=/home/wonmin-kim/AIVC/output/data-analysis
"$D/rq2_source_curve_run.sh" curve OR_HCT116X KOLF RN_K562 RN_HepG2 RN_Jurkat RN_RPE1
"$D/rq2_source_curve_run.sh" collect
echo "CHAIN_DONE $(date -Is)"
