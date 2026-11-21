#!/usr/bin/bash
# 2026-10-27 final phase D/E/F: AGENTS 10/22 command + S17 families (RN x3 + DLD1) + S18 control-PC deflation (8,0.75,0.5). Spec only.
set -u
cd /home/wonmin-kim/AIVC/vcc2026_data
export TV2_BASE_TABLE=/home/wonmin-kim/AIVC/vcc2026_data/reports/final_2026-10-27/base_table_stub   # y/controls deleted 10-11 (stub)
PY=$HOME/miniconda3/envs/vcc2026/bin/python
$PY -P scripts/tv2_final.py --controls-dir data/final_2026-10-22/controls --contexts D,E,F --out reports/final_2026-10-27/F1 --tag F1 \
  --rn-sources --fam-weights DLD1=1 --s4 --blend 0.3 --dual 1.0,3000,0.5,ag,cpc,ci5,lsb --kink 1.17,0.8 \
  --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g --sum-policy auto --sum-ens 5,10,15,30,50,100,200 --de-follow-sum \
  --sum-resid 0.25 --s4-excl KOLF,FENG --de-denoise intended --sum-tamp 0.5 --same-line 0.5 \
  --sum-cens 5,10,20,100 --de-drc 0.8,1.1 --sum-defl 8,0.75,0.5 && echo F1_SPEC_DONE
