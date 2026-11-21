#!/usr/bin/bash
# user go-ahead 2026-10-11 ("업로드해봐"): S17d first (budget-down measurement), then S18 (so the displayed entry is S18).
cd /home/wonmin-kim/AIVC/vcc2026_data
bash scripts/tv2_upload.sh S17d "S16 with the sum budget moved DOWN (margin 0.6, floor 0): official MSE budget gradient"
bash scripts/tv2_upload.sh S18 "S17 + control-cell PC deflation of the weaker half of targets (k8, eta 0.75)"
echo UPLOADS_DONE
