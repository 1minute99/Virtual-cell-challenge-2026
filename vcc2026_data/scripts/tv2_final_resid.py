"""Final phase: the S10k/S11/S12 no-embedding learned residual on the sums (builder key sum.resid) for a NEW panel.

  tv2_final_resid.py TARGETS_TXT CONTROLS_DIR CTX[:EXCL+EXCL],... SRC_DIR OUT_NPZ TAG

Per context (as ml_dr_official_prep.py, which built resid_S10k.npz for A/B/C): genes G = genes >= 10 CPM in the context's
controls with the panel targets removed; input X = W_NEW-weighted mean of the per-family source LFCs (fh_harness.fam_lfc
on SRC_DIR, the final pipeline's knn/sources -- same builder as the harness library, whole-screen centres) with the
context's excluded families left out (somatic destinations: KOLF) -> disc_residual/panel_<TAG>_<ctx>.npz; then
ml_dr_train.py noesm_strict at seeds 0/1/2 (the panel's targets are removed from every pseudo-panel; pseudo-panels come
from the fixed harness library) and the seed-mean residual -> OUT_NPZ (targets, r_<ctx>, gG_<ctx>, xnorm_<ctx>)."""
import sys, pathlib, subprocess
import numpy as np
from cell_eval2.moments import jackknife_correction  # noqa: F401
argv = list(sys.argv); sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H
import anndata as ad
from scipy import sparse
DR = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/disc_residual")
TRAIN_PY = pathlib.Path.home() / "miniconda3/envs/ct_transfer/bin/python"
tg_file, D, ctx_spec, src_dir, out, tag = pathlib.Path(argv[1]), pathlib.Path(argv[2]), argv[3], pathlib.Path(argv[4]), pathlib.Path(argv[5]), argv[6]
targets = [t for t in tg_file.read_text().split() if t]; tset = set(targets)
ctxs = []
for tok in ctx_spec.split(","):
    c, _, ex = tok.partition(":"); ctxs.append((c, set(x for x in ex.split("+") if x)))
H.TM.SRC = src_dir; H.TM._cache.clear()
seeds = ("", "@1", "@2")
res = dict(targets=np.array(targets))
for c, excl in ctxs:
    name = f"{tag}_{c}"
    a = ad.read_h5ad(D / f"context_{c}.h5ad"); assert list(map(str, a.var_names)) == list(H.TM.GENES)
    tot = np.asarray(sparse.csr_matrix(a.X).sum(0)).ravel().astype(float); cpm = tot / tot.sum() * 1e6
    gG = np.flatnonzero((cpm >= 10) & np.array([g not in tset for g in a.var_names.astype(str)]))
    fl = H.fam_lfc(targets, excl)
    num = np.zeros((len(targets), len(gG))); den = np.zeros_like(num)
    for f, M in fl.items():
        v = M[:, gG]; m = np.isfinite(v); w = H.W_NEW.get(f, 1.0); num[m] += w * v[m]; den[m] += w
    X = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
    np.savez(DR / f"panel_{name}.npz", gP=gG, gG=gG, X=X.astype(np.float32), D=np.zeros((1, 1), np.float32), Od=np.zeros((1, 1), np.float32),
             targets=np.array(targets), excl=np.array(sorted(excl)))
    print(f"context {c}: excluded {sorted(excl)}, families {sorted(fl)}, genes >= 10 CPM {len(gG)}, targets with input "
          f"{int((np.abs(X).sum(1) > 0).sum())}/{len(targets)}", flush=True)
    subprocess.run([str(TRAIN_PY), str(pathlib.Path(__file__).with_name("ml_dr_train.py")), name] + [f"noesm_strict{s}" for s in seeds], check=True)
    Rs = [np.load(DR / f"res_{name}_noesm_strict{s}.npz")["r"] for s in seeds]
    res[f"r_{c}"] = np.mean(Rs, 0).astype(np.float32); res[f"gG_{c}"] = gG; res[f"xnorm_{c}"] = np.linalg.norm(X, axis=1).astype(np.float32)
    res[f"excl_{c}"] = np.array(sorted(excl), dtype=str)          # checked on reuse by tv2_final.py (Codex R16)
    print(c, "residual |r|/|X| median", round(float(np.median(np.linalg.norm(res[f"r_{c}"], axis=1) / np.maximum(res[f"xnorm_{c}"], 1e-9))), 4), flush=True)
np.savez(out, **res); print("wrote", out, flush=True)
