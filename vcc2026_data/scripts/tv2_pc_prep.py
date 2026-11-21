"""Top principal axes of the source-library responses on the official gene axis (sum-side PC removal, harness option pr<k>x<l>).
Same construction as sch_emit's `pcr`: the recipe delta of every library target (fh_harness lib_all, nothing excluded; each
target's own knockdown column zeroed), covered rows only, centred over targets, randomized SVD (random_state 0).
Usage: tv2_pc_prep.py OUT_NPZ [K=10]"""
import sys
import numpy as np
from cell_eval2.moments import jackknife_correction  # noqa: F401
argv = list(sys.argv); sys.argv = ["x"]
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H
from sklearn.utils.extmath import randomized_svd
out = argv[1]; K = int(argv[2]) if len(argv) > 2 else 10
import os as _os
_HOLD = set(x for x in _os.environ.get("TV2_HOLDOUT_FAMS", "").split(",") if x)   # 2026-10-03 holdout rehearsals only
dlib, _ = H.base(H.LIB, _HOLD, "lib_all" + ("_ex" + "-".join(sorted(_HOLD)) if _HOLD else ""))
for i, t in enumerate(H.LIB):
    j = H.GIC.get(H.canon(t))
    if j is not None:
        dlib[i, j] = 0.0
Dl = dlib[np.abs(dlib).sum(1) > 0]
_, S, Vt = randomized_svd(Dl - Dl.mean(0), n_components=K, random_state=0)
np.savez(out, Vt=Vt.astype(np.float64), S=S, n_rows=Dl.shape[0])
tot = float(((Dl - Dl.mean(0)) ** 2).sum())
print("library rows", Dl.shape[0], "energy share of top-k:", [round(float((S[:k] ** 2).sum() / tot), 4) for k in (1, 3, 5, 10) if k <= K], flush=True)
