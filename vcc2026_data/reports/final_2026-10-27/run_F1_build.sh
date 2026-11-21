#!/usr/bin/bash
# Build F1 (final D/E/F) -> vcc prep -> sha; delete the 17 GB pred.h5ad once the .vcc exists (disk).
set -u
cd /home/wonmin-kim/AIVC/vcc2026_data
O=reports/submission_2026-09-23
bash scripts/tv2_submit_chain.sh reports/final_2026-10-27/F1/F1_spec.json F1 || { echo F1_BUILD_FAILED; exit 1; }
[ -s $O/F1.vcc ] && rm -f $O/F1.pred.h5ad && echo "F1_BUILD_OK pred.h5ad removed"
