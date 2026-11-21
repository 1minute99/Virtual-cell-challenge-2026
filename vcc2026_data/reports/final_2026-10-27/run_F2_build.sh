#!/usr/bin/bash
# Build F2 (F1 + ew 50,0.25 kept: ew_adaptive off) (final D/E/F) -> vcc prep -> sha; delete the 17 GB pred.h5ad once the .vcc exists (disk).
set -u
cd /home/wonmin-kim/AIVC/vcc2026_data
O=reports/submission_2026-09-23
bash scripts/tv2_submit_chain.sh reports/final_2026-10-27/F1/F2_spec.json F2 || { echo F2_BUILD_FAILED; exit 1; }
[ -s $O/F2.vcc ] && rm -f $O/F2.pred.h5ad && echo "F2_BUILD_OK pred.h5ad removed"
