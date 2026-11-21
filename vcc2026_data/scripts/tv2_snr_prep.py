"""Source-family SNR for a target panel (sum-side reliability weight, sch_emit 'swsnr'): for every (t, g) with >= 2
measuring families, snr = mu^2 / (mu^2 + sd^2 / n_fam) of the families' centred LFCs (c=50); 0.5 where fewer than two
families measure. The builder multiplies the SUM-side delta by it (the DE side keeps the unweighted delta).
Harness 2026-09-27 (six panels, offline exact evaluator): PDS +0.001..+0.038, with centring 0.5 PDS up on 6/6.
Usage: tv2_snr_prep.py TARGETS_CSV(first column)|TARGETS_TXT OUT_NPZ [LIB_DIR]"""
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
F = np.stack(list(fl.values())); meas = np.isfinite(F); nme = meas.sum(0)
mu = np.nanmean(np.where(meas, F, np.nan), 0); sd = np.nanstd(np.where(meas, F, np.nan), 0)
snr = np.where(nme >= 2, mu ** 2 / (mu ** 2 + np.nan_to_num(sd) ** 2 / np.maximum(nme, 1) + 1e-6), 0.5)
np.savez_compressed(_argv[2], targets=np.array(targets), snr=np.nan_to_num(snr, nan=0.5).astype(np.float32))
print("targets", len(targets), "families", len(fl), "mean snr", float(np.nanmean(snr)), flush=True)
