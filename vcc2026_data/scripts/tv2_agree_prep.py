"""Source-family sign agreement for a target panel (fh_dual 'ag' option): for every (t, g), the fraction of source
families whose centred LFC (c=50) has the same sign as the family-mean LFC, rescaled clip(2f - 1, 0, 1).
The builder multiplies it by |delta| to order the dual re-allocation strength. Uses the harness library (all families).
Usage: tv2_agree_prep.py TARGETS_CSV(first column) OUT_NPZ [LIB_DIR]
"""
import sys
import numpy as np, pandas as pd
from cell_eval2.moments import jackknife_correction  # noqa: F401
_argv = list(sys.argv); sys.argv = ["x"]
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H
if len(_argv) > 3:
    import pathlib
    H.TM.SRC = pathlib.Path(_argv[3]); H.TM._cache.clear()
p = _argv[1]
targets = (pd.read_csv(p).iloc[:, 0] if p.endswith(".csv") else pd.Series(open(p).read().split())).astype(str).tolist()
fl = H.fam_lfc(targets, set())
S = np.stack([np.nan_to_num(v) for v in fl.values()]); M = np.stack([np.isfinite(v) for v in fl.values()])
mean = np.where(M.sum(0) > 0, S.sum(0) / np.maximum(M.sum(0), 1), 0.0)
agree = (M & (np.sign(S) == np.sign(mean)[None])).sum(0)
frac = np.where(M.sum(0) > 0, agree / np.maximum(M.sum(0), 1), 0.5)
np.savez_compressed(_argv[2], targets=np.array(targets), agree=np.clip(2 * frac - 1, 0, 1).astype(np.float32))
print("targets", len(targets), "families", len(fl), flush=True)
