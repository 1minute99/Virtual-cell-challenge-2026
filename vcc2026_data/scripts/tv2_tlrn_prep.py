"""Learned per-target amplitude signal for the official build (harness token `ta<beta>_tlrn`, FINDINGS 23.12):
1. train the GBM on ALL harness panels (features rank-normalised per panel, label = rank of the per-target true cosine; same
   hyper-parameters as ml_amp_learn.py's leave-one-panel-out runs);
2. compute the same 8 test-time features for the official contexts: cross-family agreement, #families measuring >= 500 genes,
   |sum override before resid| (builder dump *_dbase_<ctx>.npy), S4-part vs fusion-part cosine, log1p(total family cells),
   max family knockdown, log1p(target CPM in the context's controls), #promoter neighbours (0 in training -> 0 here);
3. write predictions in the `sum.tamp` file format (keys agree_<ctx> hold the predicted score; the builder ranks them).
Usage: tv2_tlrn_prep.py OUT.npz DBASE_PREFIX CTX:kind ...   (DBASE_PREFIX + f"_dbase_{ctx}.npy"; kind somatic|stem)"""
import sys, json, pathlib, pickle, numpy as np, pandas as pd
from cell_eval2.moments import jackknife_correction  # noqa
_a = list(sys.argv); sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_model as TM, tv2_submit_build as B, ml_h1fuse as F
from sklearn.ensemble import GradientBoostingRegressor
R = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data"); HF = R / "reports/sc_hybrid_2026-09-25/h1fuse"
PANELS = ["H1", "KOLF", "RN_Jurkat", "RN_HepG2", "RN_RPE1", "RN_K562", "OR_HCT116X", "OR_HCT116"]
FAMS = ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1", "HEPG2", "RNJURKAT", "RPE1"]
W_OFF = {"K562": 8.0, "HCT116": 1.0, "HEK293T": 0.5, "CD4": 1.0, "FENG": 0.25, "H1": 2.0}      # official fusion part (spec_S14fk2)
rank = lambda v: (np.argsort(np.argsort(v, kind="stable"), kind="stable") + 1) / len(v)
out, dprefix, ctx_args = pathlib.Path(_a[1]), _a[2], _a[3:]
data = {p: np.load(HF / f"{p}_ampfeat.npz") for p in PANELS}
X = np.concatenate([np.apply_along_axis(rank, 0, data[p]["F"]) for p in PANELS]); y = np.concatenate([rank(data[p]["y"]) for p in PANELS])
gbm = GradientBoostingRegressor(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=0).fit(X, y)
pickle.dump(gbm, open(out.with_suffix(".gbm.pkl"), "wb"))
names = [str(n) for n in data["H1"]["names"]]; print("trained on", len(y), "targets x", len(PANELS), "panels; importances", dict(zip(names, np.round(gbm.feature_importances_, 3))), flush=True)
targets = pd.read_csv("/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls/pert_counts.csv").iloc[:, 0].astype(str).tolist()
res = {"targets": np.array(targets)}; feats = {}
S4DIR = R / "reports/submission_2026-09-30"
for a in ctx_args:
    ctx, _, kind = a.partition(":"); kind = kind or "somatic"
    fin = [f for f in FAMS if not (kind == "somatic" and f == "KOLF")]
    Fm = F.fam_means(targets, fin)                                         # min_cells 20, as the harness agreement
    fams = [f for f in fin if f in Fm]
    ag = np.full(len(targets), np.nan); nfam = np.zeros(len(targets))
    for k in range(len(targets)):
        cols = [f for f in fams if np.isfinite(Fm[f][k]).sum() > 500]; nfam[k] = len(cols); cs = []
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                x_, y_ = Fm[cols[i]][k], Fm[cols[j]][k]; mm = np.isfinite(x_) & np.isfinite(y_)
                cs.append(float(x_[mm] @ y_[mm] / (np.linalg.norm(x_[mm]) * np.linalg.norm(y_[mm]) + 1e-12)))
        if cs: ag[k] = float(np.mean(cs))
    ag = np.nan_to_num(ag, nan=float(np.nanmedian(ag)))
    Db = np.load(f"{dprefix}_dbase_{ctx}.npy").astype(np.float64); nrm = np.linalg.norm(Db, axis=1)
    S3 = 3.0 * np.load(S4DIR / f"S4noKF_ctx{ctx}.delta.npy").astype(np.float64) if kind == "somatic" else None
    d1 = np.zeros((len(targets), TM.NG))
    for mc in (20, 5, 10, 15, 30, 50, 100, 200):                          # fusion part: fixed official weights, min_cells ensemble
        Fmc = F.fam_means(targets, [f for f in W_OFF if f in fin], mc); num = 0; den = 0
        for f, w in W_OFF.items():
            if f in Fmc: m_ = np.isfinite(Fmc[f]); num = num + w * np.where(m_, Fmc[f], 0.0); den = den + w * m_
        d1 += np.clip(np.where(den > 0, num / np.maximum(den, 1e-300), 0.0), -3.0, 3.0) / 8   # Cfg.clip 3.0
    for k, t in enumerate(targets):
        j = TM.GIC.get(TM.canon(t))
        if j is not None: d1[k, j] = 0.0
    s4f = np.array([float(S3[k] @ d1[k] / (np.linalg.norm(S3[k]) * np.linalg.norm(d1[k]) + 1e-12)) for k in range(len(targets))]) if S3 is not None else np.zeros(len(targets))
    cells, kd, _ = F.fam_stats(targets, fams)
    cl = np.log1p(np.sum([cells[f] for f in fams], 0)); kdmax = np.max([kd[f] for f in fams], 0)
    M = B.load_controls(ctx); tot = np.asarray(M.sum(0)).ravel(); cpm = tot / tot.sum() * 1e6
    tex = np.array([np.log1p(cpm[TM.GIC[TM.canon(t)]]) if TM.GIC.get(TM.canon(t)) is not None else 0.0 for t in targets])
    Fx = np.stack([ag, nfam, nrm, s4f, cl, kdmax, tex, np.zeros(len(targets))], 1)
    pred = gbm.predict(np.apply_along_axis(rank, 0, Fx))
    res[f"agree_{ctx}"] = pred; res[f"feat_{ctx}"] = Fx; feats[ctx] = Fx
    print(ctx, kind, "families", fams, "| feature medians", dict(zip(names, np.round(np.median(Fx, 0), 3))),
          "| Spearman(pred, agreement)", round(float(pd.Series(pred).corr(pd.Series(ag), method="spearman")), 3), flush=True)
res["feature_names"] = np.array(names); np.savez(out, **res); print("TLRN_PREP_DONE", out)
