"""Final-phase harness: six held-out cell lines x ML arms on top of the final-phase recipe, official cell-eval2 preset.

Panels (real cells): H1 (VCC2025 Flex), KOLF, RN K562 / HepG2 / Jurkat / RPE1. For a held-out panel the sources of
its own family are excluded everywhere (base prediction AND the training data of every learned arm); learned arms are
trained on the other five panels' real pseudobulk only.

Base = final-phase recipe without the S4 part (tv2_final default): 0.5 x equal-weight uncentred fold changes (c=100)
x3 + 0.5 x centred fusion (c=50, K562 8 ...) x1.6, promoter prior + knockdown override, shared_cellshift_gated.
Arms (each a transform of the base target x gene delta D):
  base            the recipe
  amp             per-target amplitude m_t (rescaled to the base's total energy) from a regressor of the per-target least-squares slope (real vs D) on
                  target features (|D_t|, coverage, source knockdown, DepMap, destination expression of t)
  amp_shrink      same, m_t capped at 1 (only weakens targets)
  expr            per-gene multiplier m_g (rescaled to the base's total energy) from a regressor of the per-gene slope on destination / source expression
  lr64            D projected on 64 response programs (SVD of single-source target rows, held-out family excluded)
  lr128_half      0.5 D + 0.5 projection on 128 programs
  unseen_knn      panel targets treated as absent from every source: D_t = weighted mean of DepMap co-essentiality
                  neighbours' D (panel targets excluded as neighbours); promoter/knockdown override kept
  unseen_mean     same setting, D_t = mean response of all covered non-panel targets
  unseen_zero     same setting, D_t = 0 (override only)
  blend20         seen targets: 0.8 D + 0.2 DepMap-neighbour estimate (neighbours: any library target but t)
  blend20_e       blend20 rescaled to the base's total energy
  sim_soft/strong source-family weights x exp(beta * z(similarity of control expression)), beta 0.5 / 1.5
  gbm / gbm_half  gradient-boosted stacking: predicts the real LFC per (t, g) from per-family LFCs, base delta, family
                  agreement and control expression (trained on the other five panels); energy-matched; half = 0.5 base
  rf_half         same with a random forest, 0.5 base + 0.5 prediction
  blend_adapt     same as blend20 but neighbour weight 0.5 / 0.35 / 0.1 when t is covered by <=1 / 2 / >=3 source families
Usage: fh_harness.py [ARM[@sSEED] ...]   (default: all, seed 0). Idempotent: scores/<arm>/<panel>/result.json.
"""
import json, os, pathlib, subprocess, sys, time, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2.moments import jackknife_correction  # noqa: F401  (before scripts/ joins sys.path)
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_model as TM, tv2_emit as EM, tv2_score as S
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.utils.extmath import randomized_svd

R = TM.R; D = R / "reports/final_harness_2026-09-24"
TM.SRC = D / "sources"; TM._cache.clear()
TM.FAMILY.update({"RN_HepG2": "HEPG2", "RN_Jurkat": "RNJURKAT", "RN_RPE1": "RPE1"})
PANELS = json.loads((D / "panels.json").read_text())
NB = json.loads((D / "neighbours.json").read_text())
EXCL = {"H1": {"H1"}, "KOLF": {"KOLF"}, "RN_K562": {"K562", "VIP"}, "RN_HepG2": {"HEPG2"},
        "RN_Jurkat": {"RNJURKAT", "JURKAT"}, "RN_RPE1": {"RPE1"},
        "OR_HCT116": {"HCT116"}, "OR_HCT116X": {"HCT116", "HEK293T"}}          # Orion HCT116 somatic panel (2026-09-27)
FAMS = ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1", "HEPG2", "RNJURKAT", "RPE1"]
W_OLD = {f: 1.0 for f in FAMS}
W_NEW = {"K562": 8.0, "HCT116": 1.0, "HEK293T": 0.5, "CD4": 1.0, "KOLF": 0.25, "FENG": 0.25, "H1": 2.0,
         "HEPG2": 1.0, "RNJURKAT": 1.0, "RPE1": 1.0}
PARTS = [(0.5, dict(c=100.0, center=False, weights=W_OLD, promoter=True), 3.0),
         (0.5, dict(c=50.0, center=True, center_frac=1.0, weights=W_NEW, agree_shrink=0.0), 1.6)]
CACHE = D / "cache"; CACHE.mkdir(exist_ok=True)
NG, GIC, canon = TM.NG, TM.GIC, TM.canon
LIB = sorted({str(t) for s in TM.FAMILY if (TM.SRC / f"{s}.npz").exists() for t in np.load(TM.SRC / f"{s}.npz")["targets"]})


GX = {}   # harness-wide cfg override merged into EVERY part (arm token xw<p>, 2026-09-30); non-empty -> cache-key suffix


def base(targets, excl, key, wmul=None, wset=None):
    """Combined recipe delta (n x NG) and multiplicative override for a target list, sources in excl left out.
    wmul: optional {family: multiplier} applied to both parts' source weights (cache key must then differ)."""
    if GX:
        key = f"{key}_gx" + "_".join(f"{k}{v:g}" for k, v in sorted(GX.items()))
    f = CACHE / f"{key}__{'-'.join(sorted(excl))}.npz"
    if f.exists():
        z = np.load(f); return z["d"].astype(np.float64), z["ov"].astype(np.float64)
    d = 0.0; ov = None
    for w, c, a in PARTS:
        c = dict(c, **GX)
        if wset:
            c["weights"] = dict(wset)
        if wmul:
            c["weights"] = {k: v * wmul.get(k, 1.0) for k, v in c["weights"].items()}
        cfg = S.build_cfg(c)
        L, prom, tidx = TM.predict_lfc(list(targets), cfg, exclude=excl)
        dd, oo = TM.split(L, prom, tidx, a, cfg)
        d = d + w * dd; ov = oo if ov is None else np.where(np.isfinite(ov), ov, oo)
    np.savez(f, d=d.astype(np.float32), ov=ov.astype(np.float32))
    return d, ov


_real = {}
def real(panel):
    """Panel controls (CSR), gene map to the official axis, control CPM, real fold change (c=50 CPM) per target."""
    if panel in _real:
        return _real[panel]
    P = PANELS[panel]; a = ad.read_h5ad(P["path"])
    lab = a.obs["target"].astype(str).values; X = sparse.csr_matrix(a.X)
    genes = a.var_names.astype(str).tolist(); gmap = np.array([GIC.get(canon(g), -1) for g in genes])
    ctrl = X[lab == "non-targeting"]; cs = np.asarray(ctrl.sum(0)).ravel().astype(float); cpm = cs / cs.sum() * 1e6
    targets = P["targets"]; Rl = np.zeros((len(targets), len(genes)))
    for k, t in enumerate(targets):
        s = np.asarray(X[lab == t].sum(0)).ravel().astype(float)
        Rl[k] = np.log((s / s.sum() * 1e6 + 50) / (cpm + 50))
    _real[panel] = dict(ctrl=ctrl, genes=genes, gmap=gmap, cpm=cpm, R=Rl, targets=targets)
    return _real[panel]


def own_mask(panel):
    r = real(panel); gi = {g: j for j, g in enumerate(r["genes"])}
    m = np.zeros((len(r["targets"]), len(r["genes"])), bool)
    for k, t in enumerate(r["targets"]):
        j = gi.get(t)
        if j is not None:
            m[k, j] = True
    return m


def on_panel(d, panel):
    r = real(panel); ok = r["gmap"] >= 0
    out = np.zeros((d.shape[0], len(r["genes"]))); out[:, ok] = d[:, r["gmap"][ok]]
    return out


def tfeat(panel, excl, d):
    """Target-level features for the amplitude model."""
    r = real(panel); tg = r["targets"]; dp = on_panel(d, panel)
    kd = np.full(len(tg), np.nan); nf = np.zeros(len(tg))
    for s, fam in TM.FAMILY.items():
        if fam in excl or not (TM.SRC / f"{s}.npz").exists():
            continue
        z = TM.load(s)
        for k, t in enumerate(tg):
            i = z["pos"].get(canon(t)); j = GIC.get(canon(t))
            if i is None or j is None:
                continue
            nf[k] += 1
            fp, fc = z["fp"][i][j], z["fc"][i][j]
            if np.isfinite(fp) and np.isfinite(fc) and fc * 1e6 >= 5:
                v = np.log((fp * 1e6 + 1) / (fc * 1e6 + 1)); kd[k] = v if np.isnan(kd[k]) else min(kd[k], v)
    dm = DEPMAP
    tcpm = np.array([np.log10(r["cpm"][r["genes"].index(t)] + 0.01) if t in r["genes"] else np.nan for t in tg])
    srt = -np.sort(-np.abs(dp), 1)
    return np.column_stack([np.log(np.linalg.norm(dp, axis=1) + 1e-6), np.log(srt[:, :100].mean(1) + 1e-6),
                            np.log(srt[:, :10].mean(1) + 1e-6), nf, kd, [dm.get(canon(t), np.nan) for t in tg], tcpm])


def depmap_mean():
    z = np.load(R / "processed/depmap_23q4/gene_effects.npz", allow_pickle=True)
    return {canon(str(g)): float(m) for g, m in zip(z["genes"], np.nanmean(z["values"], 1))}
DEPMAP = depmap_mean()


def slopes_t(panel, d):
    dp = on_panel(d, panel); dp[own_mask(panel)] = 0.0; Rl = real(panel)["R"]
    return np.clip((dp * Rl).sum(1) / np.maximum((dp * dp).sum(1), 1e-9), -0.5, 3.0)


def gfeat(panel, excl, d):
    r = real(panel); ok = r["gmap"] >= 0; g = np.where(ok, r["gmap"], 0)
    se = []; nm = np.zeros(len(r["genes"]))
    for s, fam in TM.FAMILY.items():
        if fam in excl or not (TM.SRC / f"{s}.npz").exists():
            continue
        z = TM.load(s); fc = np.nanmedian(z["fc"][:, g], 0) * 1e6
        se.append(np.where(ok, np.log10(fc + 0.01), np.nan)); nm += np.isfinite(fc) & ok
    se = np.nanmean(np.stack(se), 0); ec = np.log10(r["cpm"] + 0.01)
    dp = on_panel(d, panel)
    return np.column_stack([ec, se, ec - se, np.log((dp ** 2).mean(0) + 1e-9), nm])


def learned_multiplier(test, kind):
    """Train on the other five panels (their own and the test family excluded from features), predict for test."""
    X, y, w = [], [], []
    for p in PANELS:
        if p == test:
            continue
        ex = EXCL[p] | EXCL[test]
        d, _ = base(PANELS[p]["targets"], ex, p)
        if kind == "amp":
            X.append(tfeat(p, ex, d)); y.append(slopes_t(p, d)); w.append(np.ones(len(y[-1])))
        else:
            dp = on_panel(d, p); dp[own_mask(p)] = 0.0; Rl = real(p)["R"]; den = (dp * dp).sum(0)
            sl = np.clip((dp * Rl).sum(0) / np.maximum(den, 1e-9), -1.0, 4.0); keep = den > 0
            X.append(gfeat(p, ex, d)[keep]); y.append(sl[keep]); w.append(den[keep] / den[keep].sum())
    X = np.vstack(X); y = np.concatenate(y); w = np.concatenate(w)
    reg = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_leaf_nodes=15,
                                        min_samples_leaf=20 if kind == "amp" else 200, random_state=0)
    reg.fit(X, y, sample_weight=w * len(w))
    d, _ = base(PANELS[test]["targets"], EXCL[test], test)
    if kind == "amp":
        yh = reg.predict(tfeat(test, EXCL[test], d)); m = np.clip(yh, 0.05, None)
        return np.clip(m / m.mean(), 0.2, 3.0), dict(train_n=int(len(y)), pred_mean=float(yh.mean()))
    F = gfeat(test, EXCL[test], d); yh = reg.predict(F)
    dp = on_panel(d, test); wt = (dp * dp).sum(0); m = np.clip(yh, 0.0, None)
    m = m / max((m * wt).sum() / max(wt.sum(), 1e-12), 1e-6)
    return np.clip(m, 0.0, 3.0), dict(train_n=int(len(y)), pred_mean=float(np.average(yh, weights=wt + 1e-12)))


def fam_similarity(panel, excl):
    """Spearman correlation between the panel's control expression and each source family's control expression."""
    from scipy.stats import spearmanr
    r = real(panel); ok = r["gmap"] >= 0; g = r["gmap"][ok]; ec = np.log10(r["cpm"][ok] + 0.01)
    sims = {}
    for fam in sorted({f for f in TM.FAMILY.values()} - set(excl)):
        vals = []
        for s_, f_ in TM.FAMILY.items():
            if f_ != fam or not (TM.SRC / f"{s_}.npz").exists():
                continue
            fc = np.nanmedian(TM.load(s_)["fc"][:, g], 0) * 1e6
            m = np.isfinite(fc)
            if m.sum() > 1000:
                vals.append(spearmanr(ec[m], np.log10(fc[m] + 0.01))[0])
        if vals:
            sims[fam] = float(np.mean(vals))
    return sims


def fam_lfc(targets, excl):
    """Per-family mean source LFC (c=50, centred) for a target list; NaN where the family does not measure (t, g)."""
    out = {}
    for fam in sorted({f for f in TM.FAMILY.values()} - set(excl)):
        num = np.zeros((len(targets), NG), np.float32); den = np.zeros((len(targets), NG), np.float32)
        for s_, f_ in TM.FAMILY.items():
            if f_ != fam or not (TM.SRC / f"{s_}.npz").exists():
                continue
            z = TM.load(s_); cen = np.nan_to_num(z["center_c50"])
            for k, t in enumerate(targets):
                i = z["pos"].get(canon(t))
                if i is None or z["ncell"][i] < 20:
                    continue
                v = np.log((z["fp"][i].astype(np.float64) * 1e6 + 50) / (z["fc"][i].astype(np.float64) * 1e6 + 50)) - cen
                m = np.isfinite(v); num[k, m] += v[m]; den[k, m] += 1
        if den.any():
            out[fam] = np.where(den > 0, num / np.maximum(den, 1), np.nan).astype(np.float32)
    return out


STACK_FAMS = ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1", "HEPG2", "RNJURKAT", "RPE1"]


def stack_features(panel, excl, d):
    """(nt*ng_ok, F) features on the panel's mapped genes: per-family LFC (NaN if excluded/missing), base delta,
    family mean / sd / count, destination and source control expression."""
    r = real(panel); ok = r["gmap"] >= 0; g = r["gmap"][ok]; tg = PANELS[panel]["targets"]
    fl = fam_lfc(tg, excl); nt = len(tg); ng = int(ok.sum())
    cols = [fl[f][:, g] if f in fl else np.full((nt, ng), np.nan, np.float32) for f in STACK_FAMS]
    F = np.stack(cols, -1)
    with np.errstate(all="ignore"):
        mu = np.nanmean(F, -1); sd = np.nanstd(F, -1); n = np.isfinite(F).sum(-1)
    ec = np.broadcast_to(np.log10(r["cpm"][ok] + 0.01), (nt, ng))
    se = []
    for s_, f_ in TM.FAMILY.items():
        if f_ in excl or not (TM.SRC / f"{s_}.npz").exists():
            continue
        se.append(np.log10(np.nanmedian(TM.load(s_)["fc"][:, g], 0) * 1e6 + 0.01))
    se = np.broadcast_to(np.nanmean(np.stack(se), 0), (nt, ng))
    X = np.concatenate([F, np.stack([d[:, g], mu, sd, n, ec, se, ec - se], -1)], -1).astype(np.float32)
    return X.reshape(nt * ng, -1)


def stacked_delta(test, kind):
    rng = np.random.default_rng(0); Xs, ys = [], []
    for p in PANELS:
        if p == test:
            continue
        ex = EXCL[p] | EXCL[test]
        d, _ = base(PANELS[p]["targets"], ex, p)
        X = stack_features(p, ex, d); r = real(p); ok = r["gmap"] >= 0
        y = np.clip(r["R"][:, ok], -3, 3); own = own_mask(p)[:, ok]
        y = y.reshape(-1); keep = np.flatnonzero(~own.reshape(-1))
        keep = rng.choice(keep, size=min(len(keep), 400000), replace=False)
        Xs.append(X[keep]); ys.append(y[keep])
    X = np.vstack(Xs); y = np.concatenate(ys)
    keepc = np.isfinite(X).mean(0) > 0.01                # drop feature columns that are (almost) empty in training
    X = X[:, keepc]
    if kind == "rf":
        from sklearn.ensemble import RandomForestRegressor
        mdl = RandomForestRegressor(n_estimators=100, max_depth=14, min_samples_leaf=50, max_features=0.5,
                                    n_jobs=8, random_state=0)
        mdl.fit(np.nan_to_num(X, nan=-9.0), y)
    else:
        mdl = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=63, min_samples_leaf=200,
                                            l2_regularization=1.0, random_state=0)
        mdl.fit(X, y)
    d, _ = base(PANELS[test]["targets"], EXCL[test], test)
    Xt = stack_features(test, EXCL[test], d)[:, keepc]
    yh = mdl.predict(np.nan_to_num(Xt, nan=-9.0) if kind == "rf" else Xt)
    r = real(test); ok = r["gmap"] >= 0
    out = np.zeros_like(d); out[:, r["gmap"][ok]] = yh.reshape(d.shape[0], -1)
    return out, dict(train_rows=int(len(y)))



def ematch(dn, dref, panel):
    """Rescale dn to the total energy of dref on the panel's genes."""
    return dn * np.sqrt((on_panel(dref, panel) ** 2).sum() / max((on_panel(dn, panel) ** 2).sum(), 1e-12))


def lib_delta(panel, ex):
    dlib, _ = base(LIB, ex, f"lib_{panel}")
    for i, t in enumerate(LIB):
        j = GIC.get(canon(t))
        if j is not None:
            dlib[i, j] = 0.0
    return dlib


def knn_estimate(panel, ex, tg):
    dlib = lib_delta(panel, ex)
    cov = np.abs(dlib).sum(1) > 0; pos = {canon(t): i for i, t in enumerate(LIB)}
    knn = np.zeros((len(tg), NG)); has = np.zeros(len(tg), bool)
    for k_, t in enumerate(tg):
        cand = [(pos[x], c) for x, c in NB.get(canon(t), []) if x in pos and x != canon(t) and cov[pos[x]] and c > 0][:10]
        if cand:
            wv = np.array([c for _, c in cand]) ** 2
            knn[k_] = (wv[:, None] * dlib[[i for i, _ in cand]]).sum(0) / wv.sum(); has[k_] = True
    return knn, has


def train_rows(test, fn, per_panel=400000, seed=0):
    """Collect (features, label, extra) from the five training panels; fn(p, ex, d) -> (X, y, keep_mask, extra)."""
    rng = np.random.default_rng(seed); out = []
    for p in PANELS:
        if p == test:
            continue
        ex = EXCL[p] | EXCL[test]
        d, _ = base(PANELS[p]["targets"], ex, p)
        X, y, m, extra = fn(p, ex, d)
        idx = np.flatnonzero(m); idx = rng.choice(idx, size=min(len(idx), per_panel), replace=False)
        out.append((X[idx], y[idx], None if extra is None else extra[idx]))
    X = np.vstack([o[0] for o in out]); y = np.concatenate([o[1] for o in out])
    ex_ = None if out[0][2] is None else np.concatenate([o[2] for o in out])
    return X, y, ex_


def learn_family_weights(test):
    """Low-capacity stacking: non-negative family weights w_f for a weighted mean of family LFCs (c=50, centred),
    fitted (with a free global scale) to the real LFC of the five training panels."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(0); blocks = []
    for p in PANELS:
        if p == test:
            continue
        ex = EXCL[p] | EXCL[test]; tg = PANELS[p]["targets"]; r = real(p); ok = r["gmap"] >= 0
        fl = fam_lfc(tg, ex); own = own_mask(p)[:, ok]
        F = np.stack([fl[f][:, r["gmap"][ok]] if f in fl else np.full((len(tg), int(ok.sum())), np.nan, np.float32)
                      for f in STACK_FAMS], -1).reshape(-1, len(STACK_FAMS))
        y = np.clip(r["R"][:, ok], -3, 3).reshape(-1); m = ~own.reshape(-1) & np.isfinite(F).any(1)
        idx = rng.choice(np.flatnonzero(m), size=min(int(m.sum()), 300000), replace=False)
        blocks.append((F[idx], y[idx]))
    F = np.vstack([b[0] for b in blocks]); y = np.concatenate([b[1] for b in blocks])
    M = np.isfinite(F); F0 = np.nan_to_num(F)
    avail = M.any(0)
    def loss(theta):
        w = np.log1p(np.exp(theta[:-1])); a = theta[-1]
        num = (F0 * w).sum(1); den = (M * w).sum(1)
        pred = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
        return float(((a * pred - y) ** 2).mean())
    th0 = np.zeros(len(STACK_FAMS) + 1); th0[-1] = 0.5
    res = minimize(loss, th0, method="L-BFGS-B", options=dict(maxiter=200))
    w = np.log1p(np.exp(res.x[:-1])); w = w / w[avail].max()
    return {f: (float(x) if avail[i] else 0.0) for i, (f, x) in enumerate(zip(STACK_FAMS, w))}, float(res.x[-1])


def gene_graph(panel, ex, k=15, ntop=6000):
    """Gene-gene kNN graph from co-response across library targets (held-out family excluded)."""
    f = CACHE / f"genegraph_{k}_{ntop}__{'-'.join(sorted(ex))}.npz"
    if f.exists():
        z = np.load(f); return z["genes"], z["nbr"], z["w"]
    dlib = lib_delta(panel, ex); cov = np.abs(dlib).sum(1) > 0; D = dlib[cov]
    v = D.var(0); genes = np.argsort(-v)[:ntop]
    Z = D[:, genes]; Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    C = (Z.T @ Z) / Z.shape[0]; np.fill_diagonal(C, -1)
    nbr = np.argsort(-C, 1)[:, :k]; w = np.take_along_axis(C, nbr, 1).clip(0)
    np.savez(f, genes=genes, nbr=nbr, w=w); return genes, nbr, w


def resid_pred(test, kind="gbm"):
    """Target-centred residual prediction (real LFC scale) on the test panel's mapped genes, cached."""
    f = CACHE / f"resid_{kind}_{test}.npz"
    if f.exists():
        z = np.load(f); return z["R"].astype(np.float64), float(z["a"])
    tg = PANELS[test]["targets"]; ex = EXCL[test]
    d, _ = base(tg, ex, test)
    def fn(p, ex_, dd):
        X = stack_features(p, ex_, dd); r = real(p); ok = r["gmap"] >= 0
        dp = dd[:, r["gmap"][ok]]; y = r["R"][:, ok]; own = own_mask(p)[:, ok]
        a = float((dp[~own] * y[~own]).sum() / max((dp[~own] ** 2).sum(), 1e-9))
        return X, np.clip(y - a * dp, -3, 3).reshape(-1), ~own.reshape(-1), np.full(y.size, a)
    X, y, av = train_rows(test, fn)
    keepc = np.isfinite(X).mean(0) > 0.01
    mdl = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=500,
                                        l2_regularization=1.0, random_state=0).fit(X[:, keepc], y)
    R_ = mdl.predict(stack_features(test, ex, d)[:, keepc]).reshape(len(tg), -1)
    R_ = R_ - R_.mean(0, keepdims=True)
    np.savez(f, R=R_.astype(np.float32), a=float(av.mean()))
    return R_, float(av.mean())


def programs(excl, k):
    f = CACHE / f"programs_{k}__{'-'.join(sorted(excl))}.npy"
    if f.exists():
        return np.load(f)
    rng = np.random.default_rng(0); rows = []
    for s, fam in TM.FAMILY.items():
        if fam in excl or not (TM.SRC / f"{s}.npz").exists():
            continue
        z = TM.load(s); cen = np.nan_to_num(z["center_c50"])
        idx = np.flatnonzero(z["ncell"] >= 20)
        idx = rng.choice(idx, size=min(len(idx), 1500), replace=False)
        v = np.log((z["fp"][idx].astype(np.float64) * 1e6 + 50) / (z["fc"][idx].astype(np.float64) * 1e6 + 50)) - cen
        rows.append(np.nan_to_num(v).astype(np.float32))
    M = np.vstack(rows)
    _, _, Vt = randomized_svd(M, n_components=k, random_state=0)
    np.save(f, Vt.astype(np.float32)); return Vt


def arm_delta(arm, panel):
    tg = PANELS[panel]["targets"]; ex = EXCL[panel]; info = {}
    d, ov = base(tg, ex, panel)
    if arm == "base":
        pass
    elif arm in ("amp", "amp_shrink"):
        m, info = learned_multiplier(panel, "amp")
        if arm == "amp_shrink":
            m = np.minimum(m, 1.0)                       # only weakens targets
        else:                                            # redistribute at the base's total energy
            e = (on_panel(d, panel) ** 2).sum(1); m = m / np.sqrt((m ** 2 * e).sum() / e.sum())
        d = d * m[:, None]; info["m_quantiles"] = np.percentile(m, [5, 25, 50, 75, 95]).round(3).tolist()
    elif arm == "expr":
        m, info = learned_multiplier(panel, "expr")
        r = real(panel); ok = r["gmap"] >= 0
        e = (on_panel(d, panel) ** 2).sum(0); m = m / np.sqrt((m ** 2 * e).sum() / e.sum())
        mg = np.ones(NG); mg[r["gmap"][ok]] = m[ok]
        d = d * mg[None, :]; info["m_quantiles"] = np.percentile(m[ok], [5, 25, 50, 75, 95]).round(3).tolist()
    elif arm in ("lr64", "lr128_half"):
        k, lam = (64, 1.0) if arm == "lr64" else (128, 0.5)
        V = programs(ex, k).astype(np.float64)
        d = lam * (d @ V.T) @ V + (1 - lam) * d
    elif arm.startswith("unseen_"):
        pan = {canon(t) for t in tg}
        if arm == "unseen_zero":
            d = np.zeros_like(d)
        else:
            pool = [t for t in LIB if canon(t) not in pan]
            dpool, _ = base(pool, ex, f"pool_{panel}")
            for i, t in enumerate(pool):                 # a neighbour's own knockdown is not part of t's response
                j = GIC.get(canon(t))
                if j is not None:
                    dpool[i, j] = 0.0
            cov = np.abs(dpool).sum(1) > 0; pos = {canon(t): i for i, t in enumerate(pool)}
            if arm == "unseen_mean":
                d = np.repeat(dpool[cov].mean(0, keepdims=True), len(tg), 0)
            else:
                out = np.zeros_like(d); nn = []
                for k_, t in enumerate(tg):
                    cand = [(pos[x], c) for x, c in NB.get(canon(t), []) if x in pos and cov[pos[x]] and c > 0][:10]
                    nn.append(len(cand))
                    if cand:
                        wv = np.array([c for _, c in cand]) ** 2
                        out[k_] = (wv[:, None] * dpool[[i for i, _ in cand]]).sum(0) / wv.sum()
                    else:
                        out[k_] = dpool[cov].mean(0)
                d = out; info["neighbours_median"] = float(np.median(nn)); info["no_neighbour"] = int((np.array(nn) == 0).sum())
    elif arm in ("blend20", "blend20_e", "blend_adapt"):
        # seen targets: mix the measured-source delta with the DepMap-neighbour estimate; neighbours may be any library
        # target except t itself (other panel targets included: they have source predictions at test time too)
        dlib, _ = base(LIB, ex, f"lib_{panel}")
        for i, t in enumerate(LIB):
            j = GIC.get(canon(t))
            if j is not None:
                dlib[i, j] = 0.0
        cov = np.abs(dlib).sum(1) > 0; pos = {canon(t): i for i, t in enumerate(LIB)}
        knn = np.zeros_like(d); has = np.zeros(len(tg), bool)
        for k_, t in enumerate(tg):
            cand = [(pos[x], c) for x, c in NB.get(canon(t), []) if x in pos and x != canon(t) and cov[pos[x]] and c > 0][:10]
            if cand:
                wv = np.array([c for _, c in cand]) ** 2
                knn[k_] = (wv[:, None] * dlib[[i for i, _ in cand]]).sum(0) / wv.sum(); has[k_] = True
        if arm in ("blend20", "blend20_e"):
            w = np.where(has, 0.2, 0.0)
        else:
            nf = np.zeros(len(tg))
            for fam in {f for f in TM.FAMILY.values() if f not in ex}:
                srcs = [s_ for s_, f_ in TM.FAMILY.items() if f_ == fam and (TM.SRC / f"{s_}.npz").exists()]
                cf = np.zeros(len(tg), bool)
                for s_ in srcs:
                    z = TM.load(s_)
                    cf |= np.array([(z["pos"].get(canon(t)) is not None) and z["ncell"][z["pos"][canon(t)]] >= 20 for t in tg])
                nf += cf
            w = np.where(nf <= 1, 0.5, np.where(nf == 2, 0.35, 0.1)) * has
            info["n_fam_hist"] = {str(k): int((nf == k).sum()) for k in range(0, 8) if (nf == k).any()}
        e0 = (on_panel(d, panel) ** 2).sum()
        d = (1 - w)[:, None] * d + w[:, None] * knn
        if arm == "blend20_e":                           # same total energy as the base
            d = d * np.sqrt(e0 / max((on_panel(d, panel) ** 2).sum(), 1e-12))
        info["with_neighbours"] = int(has.sum()); info["mean_w"] = float(w.mean())
    elif arm in ("sim_soft", "sim_strong"):
        sims = fam_similarity(panel, ex); v = np.array(list(sims.values()))
        beta = 0.5 if arm == "sim_soft" else 1.5
        wmul = {f: float(np.exp(beta * (x - v.mean()) / (v.std() + 1e-9))) for f, x in sims.items()}
        d, ov = base(tg, ex, f"{panel}_{arm}", wmul=wmul)
        info = dict(similarity={k: round(x, 3) for k, x in sims.items()}, wmul={k: round(x, 2) for k, x in wmul.items()})
    elif arm in ("gbm", "gbm_half", "rf_half"):
        e0 = (on_panel(d, panel) ** 2).sum()
        yh, info = stacked_delta(panel, "rf" if arm == "rf_half" else "gbm")
        yh = yh * np.sqrt(e0 / max((on_panel(yh, panel) ** 2).sum(), 1e-12))
        r = real(panel); ok = r["gmap"] >= 0
        info["corr_with_base"] = float(np.mean([np.corrcoef(on_panel(d, panel)[i], on_panel(yh, panel)[i])[0, 1]
                                                for i in range(d.shape[0]) if on_panel(d, panel)[i].std() > 0]))
        d = yh if arm == "gbm" else 0.5 * d + 0.5 * yh
        d = d * np.sqrt(e0 / max((on_panel(d, panel) ** 2).sum(), 1e-12))
    elif arm.startswith("blend") and arm.endswith("_e") and arm != "blend20_e":      # blend30_e, blend40_e ...
        w0 = int(arm[5:-2]) / 100.0
        knn, has = knn_estimate(panel, ex, tg); w = w0 * has
        d = ematch((1 - w)[:, None] * d + w[:, None] * knn, d, panel); info = dict(w=w0)
    elif arm == "combo":                                  # sim_soft source weights + blend20_e
        sims = fam_similarity(panel, ex); v = np.array(list(sims.values()))
        wmul = {f: float(np.exp(0.5 * (x - v.mean()) / (v.std() + 1e-9))) for f, x in sims.items()}
        ds, ov = base(tg, ex, f"{panel}_sim_soft", wmul=wmul)
        knn, has = knn_estimate(panel, ex, tg); w = 0.2 * has
        d = ematch((1 - w)[:, None] * ds + w[:, None] * knn, ds, panel)
    elif arm == "gain":                                   # per-(t,g) multiplier g(features): E[r] ~ d * g
        def fn(p, ex_, dd):
            X = stack_features(p, ex_, dd); r = real(p); ok = r["gmap"] >= 0
            dp = dd[:, r["gmap"][ok]].reshape(-1); y = r["R"][:, ok].reshape(-1)
            m = (np.abs(dp) > 1e-3) & ~own_mask(p)[:, ok].reshape(-1)
            ratio = np.clip(np.where(np.abs(dp) > 1e-3, y / np.where(dp == 0, 1, dp), 0), -2, 4)
            return X, ratio, m, dp ** 2
        X, y, wt = train_rows(panel, fn)
        keepc = np.isfinite(X).mean(0) > 0.01
        reg = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=500,
                                            l2_regularization=1.0, random_state=0)
        reg.fit(X[:, keepc], y, sample_weight=wt / wt.mean())
        r = real(panel); ok = r["gmap"] >= 0
        g = np.clip(reg.predict(stack_features(panel, ex, d)[:, keepc]), 0, 3).reshape(len(tg), -1)
        dn = d.copy(); dn[:, r["gmap"][ok]] = d[:, r["gmap"][ok]] * g
        d = ematch(dn, d, panel); info = dict(g_quantiles=np.percentile(g, [5, 25, 50, 75, 95]).round(3).tolist())
    elif arm == "nnls_fam":                               # learned family weights (10 numbers), both parts
        wl, a = learn_family_weights(panel)
        dn, ov = base(tg, ex, f"{panel}_nnls", wset={f: max(x, 1e-6) for f, x in wl.items()})
        d = ematch(dn, d, panel); info = dict(weights={k: round(x, 3) for k, x in wl.items()}, scale=round(a, 3))
    elif arm == "gene_smooth":                            # gene-axis neighbours: 0.8 d + 0.2 d over co-responding genes
        genes, nbr, w = gene_graph(panel, ex)
        ds = d.copy(); ws = w.sum(1)
        sm = (d[:, genes[nbr]] * w[None]).sum(-1) / np.maximum(ws, 1e-9)[None]
        ds[:, genes] = np.where(ws[None] > 0, sm, d[:, genes])
        d = ematch(0.8 * d + 0.2 * ds, d, panel)
    elif arm == "resp_knn":                               # response-space target neighbours (library), 20 % blend
        dlib = lib_delta(panel, ex); cov = np.abs(dlib).sum(1) > 0
        L = dlib[cov]; names = [t for t, c in zip(LIB, cov) if c]
        Ln = L / (np.linalg.norm(L, axis=1, keepdims=True) + 1e-9)
        dn_ = d / (np.linalg.norm(d, axis=1, keepdims=True) + 1e-9)
        S_ = dn_ @ Ln.T
        for k_, t in enumerate(tg):
            S_[k_, [i for i, x in enumerate(names) if canon(x) == canon(t)]] = -1
        nb = np.argsort(-S_, 1)[:, :10]; sw = np.take_along_axis(S_, nb, 1).clip(0) ** 2
        est = (sw[:, :, None] * L[nb]).sum(1) / np.maximum(sw.sum(1), 1e-9)[:, None]
        d = ematch(0.8 * d + 0.2 * est, d, panel); info = dict(mean_top_cos=float(np.take_along_axis(S_, nb, 1).mean()))
    elif arm.startswith("gbm_resid_tc_l") or arm.startswith("b30r"):
        # gbm_resid_tc_lXX: base + (XX/100) x target-centred GBM residual; b30rXX: blend30_e + the same residual
        lam = int(arm.split("_l")[-1] if arm.startswith("gbm") else arm[4:]) / 100.0
        R_, a_mean = resid_pred(panel, "gbm")
        dref = d
        if arm.startswith("b30r"):
            knn, has = knn_estimate(panel, ex, tg); w = 0.3 * has
            d = ematch((1 - w)[:, None] * d + w[:, None] * knn, d, panel)
        r = real(panel); ok = r["gmap"] >= 0
        dn = d.copy(); dn[:, r["gmap"][ok]] += lam * R_ / max(a_mean, 1e-3)
        d = ematch(dn, dref, panel); info = dict(lam=lam, a_mean=a_mean)
    elif arm in ("gbm_resid_tc", "mlp_resid_tc"):         # residual on the calibrated base, target-centred, lambda 0.5
        def fn(p, ex_, dd):
            X = stack_features(p, ex_, dd); r = real(p); ok = r["gmap"] >= 0
            dp = dd[:, r["gmap"][ok]]; y = r["R"][:, ok]; own = own_mask(p)[:, ok]
            a = float((dp[~own] * y[~own]).sum() / max((dp[~own] ** 2).sum(), 1e-9))
            res = np.clip(y - a * dp, -3, 3).reshape(-1)
            return X, res, ~own.reshape(-1), np.full(res.shape, a)
        X, y, av = train_rows(panel, fn)
        keepc = np.isfinite(X).mean(0) > 0.01; X = X[:, keepc]
        Xt = stack_features(panel, ex, d)[:, keepc]
        if arm == "gbm_resid_tc":
            mdl = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=500,
                                                l2_regularization=1.0, random_state=0).fit(X, y)
            yh = mdl.predict(Xt)
        else:
            from sklearn.neural_network import MLPRegressor
            from sklearn.preprocessing import StandardScaler
            Xm = np.concatenate([np.nan_to_num(X), np.isnan(X)], 1); Xtm = np.concatenate([np.nan_to_num(Xt), np.isnan(Xt)], 1)
            sc = StandardScaler().fit(Xm)
            mdl = MLPRegressor(hidden_layer_sizes=(64, 32), alpha=1e-2, batch_size=4096, learning_rate_init=1e-3,
                               max_iter=40, early_stopping=True, validation_fraction=0.1, n_iter_no_change=4, random_state=0)
            mdl.fit(sc.transform(Xm)[:1000000], y[:1000000]); yh = mdl.predict(sc.transform(Xtm))
        r = real(panel); ok = r["gmap"] >= 0
        R_ = yh.reshape(len(tg), -1); R_ = R_ - R_.mean(0, keepdims=True)          # target-specific part only
        dn = d.copy(); dn[:, r["gmap"][ok]] += 0.5 * R_ / max(float(av.mean()), 1e-3)
        d = ematch(dn, d, panel); info = dict(a_mean=float(av.mean()), resid_rms=float(np.sqrt((R_ ** 2).mean())))
    else:
        raise ValueError(arm)
    return d, ov, info


def score(arm, panel):
    out = D / "scores" / arm / panel
    if (out / "result.json").exists():
        return json.loads((out / "result.json").read_text())
    out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    name, _, sd = arm.partition("@s")                   # "arm@sN": same arm, emission seed N
    seed = int(sd) if sd else 0
    d, ov, info = arm_delta(name, panel)
    r = real(panel); P = PANELS[panel]
    Ld = on_panel(d, panel); ok = r["gmap"] >= 0
    Od = np.full(Ld.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
    X = EM.emit(r["ctrl"], Ld, P["n_cells"], seed, "shared_cellshift_gated", override=Od)
    tg = r["targets"]
    pred = ad.AnnData(X=sparse.vstack([X, r["ctrl"]], format="csr"),
                      obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(tg, P["n_cells"]), *["non-targeting"] * r["ctrl"].shape[0]])},
                                       index=[f"p{i}" for i in range(X.shape[0] + r["ctrl"].shape[0])]),
                      var=pd.DataFrame(index=r["genes"]))
    pp = out / "pred.h5ad"; pred.write_h5ad(pp); del pred, X
    cache = (R / "reports/transfer_v2_2026-09-23/cache_real_H1") if panel == "H1" else (D / f"cache_real_{panel}")
    cmd = [str(S.CE), "run", "-ap", str(pp), "-ar", P["path"], "--preset", "vcc2026", "--pert-col", "target",
           "--control", "non-targeting", "--input-type", "counts", "-o", str(out / "cell_eval2"), "--cache-real", str(cache),
           "--set", "de.backend=pdex", "--set", "num_threads=2", "--set", "gather_threads=2"]
    env = dict(os.environ, NUMBA_NUM_THREADS="2", OMP_NUM_THREADS="2")
    with (out / "scorer.log").open("w") as g:
        subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
    agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
    raw = {k: float(agg.loc["mean", k]) for k in S.CONST}
    comp, scaled = S.composite(raw)
    res = dict(arm=arm, panel=panel, raw=raw, scaled=scaled, composite=comp, info=info, seconds=round(time.time() - t0, 1))
    (out / "result.json").write_text(json.dumps(res, indent=1)); pp.unlink()
    print(json.dumps(dict(arm=arm, panel=panel, composite=round(comp, 4), pds=round(raw["pds_cosine"], 4),
                          mse=round(raw["expr_mse_unbiased_capped_norm"], 3), info=info, secs=res["seconds"])), flush=True)
    return res


ARMS = ["base", "amp", "expr", "lr64", "amp_shrink", "lr128_half", "unseen_knn", "unseen_mean", "unseen_zero"]
if __name__ == "__main__":
    arms = sys.argv[1:] or ARMS
    for arm in arms:
        for panel in PANELS:
            try:
                score(arm, panel)
            except Exception as ex:
                import traceback; traceback.print_exc()
                print(json.dumps(dict(arm=arm, panel=panel, error=repr(ex))), flush=True)
    print("ALL_DONE", flush=True)
