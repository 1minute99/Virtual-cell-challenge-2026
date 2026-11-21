"""WHAT IS SAME-LINEAGE DATA WORTH? GSE295214 (Weissman lab, H1, 10x 5', genome-scale) vs Arc's VCC2025 H1 truth (Flex).
DIAGNOSTIC ONLY -- H1GW is never a source in the H1 harness for decisions. Scorer space as in ml_ceiling_v2.py.
(1) Per covered panel target: cosine with the truth of H1GW alone, our S11 sum side, the best per-target NNLS mix of our
    families + S4 (half genes), and that mix + H1GW -- by strength stratum.
(2) Offline evaluator: our S11 sum side with H1GW blended into the covered targets' delta (energy-matched per target)
    at w = 0.25 / 0.5 / 1.0 -- dPDS / dMSE / J on the full H1 panel, i.e. what "having same-line data" buys in score."""
import sys, pathlib, json
import numpy as np
from cell_eval2.moments import jackknife_correction  # noqa
sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H, sch_sumeval as SE, sch_emit as EMT, tv2_emit as EM, tv2_model as TM
from sklearn.utils.extmath import randomized_svd
from scipy.optimize import nnls
CACHE0 = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/src_ablation/cache_hgv")
H.CACHE = CACHE0; CACHE0.mkdir(parents=True, exist_ok=True)
MC = (5, 10, 15, 30, 50, 100, 200); p = "H1"
z = np.load(TM.SRC_EXTRA / "H1GW.npz"); pos = {TM.canon(str(t)): i for i, t in enumerate(z["targets"])}; cen = np.nan_to_num(z["center_c50"])
r = H.real(p); ok = r["gmap"] >= 0; n = H.PANELS[p]["n_cells"]; ctrl = r["ctrl"]; tg = list(r["targets"]); genes = r["genes"]
Rm = SE.real_moments(p); rp = list(Rm["perts"]); ir = np.array([rp.index(t) for t in tg])
mr = Rm["means"][ir]; cr = Rm["jk"][ir]; m0 = Rm["means"][rp.index("non-targeting")]; c0 = float(Rm["jk"][rp.index("non-targeting")])
fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]
gi = {g: k for k, g in enumerate(genes)}; pc = np.array(sorted({gi[t] for t in tg if t in gi})); keep = np.ones(len(genes), bool); keep[pc] = False
def sdelta(L):
    mp = SE.bulk_lognorm_means(np.round(fm[None, :] * np.exp(np.nan_to_num(L)) * n), SE.TS); return (mp - m0[None, :])[:, keep]
R = (mr - m0[None, :])[:, keep]; tot = (R ** 2).sum(1); sig = tot - cr - c0; o = np.argsort(-sig)
Lg = np.full((len(tg), len(genes)), np.nan); cov = np.zeros(len(tg), bool)
for k, t in enumerate(tg):
    i = pos.get(TM.canon(t))
    if i is None or z["ncell"][i] < 20: continue
    v = np.log((z["fp"][i].astype(np.float64) * 1e6 + 50) / (z["fc"][i].astype(np.float64) * 1e6 + 50)) - cen
    Lg[k, ok] = v[r["gmap"][ok]]; cov[k] = True
print(f"H1 panel targets covered by H1GW (>= 20 cells): {cov.sum()} / {len(tg)}", flush=True)
cand = {"H1GW alone": sdelta(Lg)}
fl = H.fam_lfc(tg, H.EXCL[p])
for f, Lf in fl.items():
    B = np.full((len(tg), len(genes)), np.nan); B[:, ok] = Lf[:, r["gmap"][ok]]; cand["fam:" + f] = sdelta(B)
d, ov = EMT.delta_s4(p, "Yorig", 0.3, 0.5, mc_list=MC)
L0 = H.on_panel(d, p); Od = np.full(L0.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
dlib = H.lib_delta(p, set(H.EXCL[p])); Dl = np.zeros((dlib.shape[0], len(genes))); Dl[:, ok] = dlib[:, r["gmap"][ok]]
Dl = Dl[np.abs(Dl).sum(1) > 0]; _, _, V3 = randomized_svd(Dl - Dl.mean(0), n_components=3, random_state=0); del dlib, Dl
gk, _ = EMT.kink_gamma(ctrl, L0, Od, n, factor=1.17, margin=0.8); Bt = EM.budget_B(ctrl, L0, Od, n)(max(1.0, gk))
def ops(L):
    mt = np.nanmean(np.where(np.isfinite(Od), np.nan, L), 0); X = L - 0.5 * np.nan_to_num(mt)[None, :]
    P_ = X @ V3.T; f_ = float((P_ ** 2).sum() / (X ** 2).sum()); X = X - 0.75 * float(np.clip((0.25 - f_) / 0.10, 0, 1)) * P_ @ V3
    return X * EM.budget_match(ctrl, X, Od, n, Bt)
X0 = ops(L0); cand["OURS(S11 sum side)"] = sdelta(np.where(np.isfinite(Od), Od, X0))
fS = EMT.S4_NEW / "H1" / "film_orig.npy"; cand["comp:S4"] = sdelta(H.on_panel(3.0 * np.load(fS).astype(np.float64), p))
cosv = lambda A, B_: np.sum(A * B_, 1) / (np.linalg.norm(A, axis=1) * np.linalg.norm(B_, axis=1) + 1e-12)
cos = {k: cosv(v, R) for k, v in cand.items()}
rng = np.random.default_rng(0); half = rng.random(keep.sum()) < 0.5
def nnls_half(keys):
    out = np.full(len(tg), np.nan)
    for t in range(len(tg)):
        A = np.stack([cand[k][t] for k in keys]); A = A[np.linalg.norm(A, axis=1) > 0]
        if not len(A): continue
        w, _ = nnls(A[:, half].T, R[t, half]); y = w @ A[:, ~half]; out[t] = y @ R[t, ~half] / (np.linalg.norm(y) * np.linalg.norm(R[t, ~half]) + 1e-12)
    return out
fams = [k for k in cand if k.startswith("fam:")] + ["comp:S4"]
cos["NNLS ours-sources (half)"] = nnls_half(fams); cos["NNLS ours-sources + H1GW (half)"] = nnls_half(fams + ["H1GW alone"])
ncl = np.sqrt(np.clip(sig, 0, None) / np.maximum(tot, 1e-12))
oc = [t for t in o if cov[t]]
strata = (("top-10c", oc[:10]), ("top-30c", oc[:30]), ("top-60c", oc[:60]), ("rest-c", oc[60:]))
print("  (strata over COVERED targets, by strength)" + "".join(f"{lab:>10}" for lab, _ in strata), flush=True)
print(f"  {'NOISE CEILING':<42}" + "".join(f"{np.nanmean(ncl[s]):10.3f}" for _, s in strata), flush=True)
for k in ["H1GW alone", "NNLS ours-sources + H1GW (half)", "NNLS ours-sources (half)", "comp:S4", "OURS(S11 sum side)"] + [k for k in cand if k.startswith("fam:")]:
    print(f"  {k:<42}" + "".join(f"{np.nanmean(cos[k][s]):10.3f}" for _, s in strata), flush=True)
base = SE.evaluate(p, tg, genes, ctrl, X0, Od, n); print(f"\n  offline S11 sum side: PDS {base['pds']:.4f} MSE {base['mse']:.4f}", flush=True)
res = dict(cos={k: {lab: float(np.nanmean(v[s])) for lab, s in strata} for k, v in cos.items()}, base=base)
for w in (0.25, 0.5, 1.0):
    L = L0.copy()
    for k in np.flatnonzero(cov):
        g = np.nan_to_num(Lg[k]); nb = np.linalg.norm(np.nan_to_num(L0[k])); ng = np.linalg.norm(g)
        if ng > 0: L[k] = (1 - w) * L0[k] + w * g * (nb / ng)
    o2 = SE.evaluate(p, tg, genes, ctrl, ops(L), Od, n); J = 2.23 * (o2["pds"] - base["pds"]) - 1.05 * (o2["mse"] - base["mse"])
    res[f"blend{w}"] = dict(o2, J=J)
    print(f"  H1GW blended into covered targets w={w}: dPDS {o2['pds']-base['pds']:+.4f} dMSE {o2['mse']-base['mse']:+.4f} J {J:+.4f} (~{J/6:+.4f} composite)", flush=True)
json.dump(res, open("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/err_decomp/h1gw_value.json", "w"), indent=1, default=float)
import shutil; shutil.rmtree(CACHE0, ignore_errors=True)
print("H1GWV_DONE", flush=True)
