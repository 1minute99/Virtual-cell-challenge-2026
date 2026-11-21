"""Learned fusion trained on the H1-2025 TRUTH -- step 1 data assembly and exact OOF evaluation
(design: docs/superpowers/specs/2026-10-01-h1-truth-fusion-design.md; Codex R25 corrections applied).

  assemble PANEL          donor family means for the panel targets (8 min_cells members), donor side statistics,
                          destination truth on the official axis -> reports/sc_hybrid_2026-09-25/h1fuse/PANEL_data.npz
  eval PANEL OOF.npz TAG  per-FOLD full-panel construction (fold k's model predicts all targets; the complete sum
                          construction -- S4 part, neighbour blend, energy match, ct, guarded pcr, budget fixed from the
                          base -- runs on all of them; only fold k's targets are kept), exact offline PDS / MSE per target
                          aggregated over folds (MSE = sum num / sum dist), base (fixed weights, parity vs delta_s4) and
                          per-target NNLS oracle; fusion-only diagnostic. -> h1fuse/PANEL_eval_TAG.json
The sum construction is ml_gate_pilot.py's path, unchanged (its base reproduces delta_s4 to 1e-5)."""
import sys, json, pathlib, time
import numpy as np
from scipy.optimize import nnls
from cell_eval2.moments import jackknife_correction  # noqa (import order rule)
_argv = list(sys.argv); sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H, sch_sumeval as SE, sch_emit as EMT, tv2_emit as EM, tv2_model as TM
from sklearn.utils.extmath import randomized_svd

OUTER = {"H1": {"H1", "H1GW"}, "KOLF": {"KOLF"}, "OR_HCT116": {"HCT116"}, "OR_HCT116X": {"HCT116", "HEK293T"},
         "RN_K562": {"K562", "VIP"}, "RN_HepG2": {"HEPG2"}, "RN_Jurkat": {"RNJURKAT", "JURKAT"}, "RN_RPE1": {"RPE1"}}
S4KEY = {"H1": "Yorig", "KOLF": "F", "OR_HCT116": "YorRnoKOLFFENG", "OR_HCT116X": "YorXnoKOLFFENG", "RN_K562": "YnoKOLFFENG",
         "RN_HepG2": "YnoKOLFFENG", "RN_Jurkat": "YnoKOLFFENG", "RN_RPE1": "YnoKOLFFENG"}
R = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data"); OUTD = R / "reports/sc_hybrid_2026-09-25/h1fuse"; OUTD.mkdir(parents=True, exist_ok=True)
W_NEW = dict(H.W_NEW); FAMS = [f for f in W_NEW if W_NEW[f] > 0]
MC = [None] + list(EMT.MC7)
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)


def fam_means(targets, fams, min_cells=20):
    out = {}
    for fam in fams:
        num = np.zeros((len(targets), H.NG)); den = np.zeros((len(targets), H.NG))
        for s_, f_ in TM.FAMILY.items():
            if f_ != fam or not TM._src_path(s_).exists():
                continue
            z = TM.load(s_); cen = np.where(np.isfinite(z["center_c50"]), z["center_c50"], 0.0)
            for k, t in enumerate(targets):
                i = z["pos"].get(H.canon(t))
                if i is None or z["ncell"][i] < min_cells or H.canon(t) in TM.DROP_TARGETS:
                    continue
                v = np.log((z["fp"][i].astype(np.float64) * 1e6 + 50) / (z["fc"][i].astype(np.float64) * 1e6 + 50)) - cen
                m = np.isfinite(v); num[k, m] += v[m]; den[k, m] += 1.0
        if den.any():
            out[fam] = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    return out


def fam_stats(targets, fams):
    cells = {f: np.zeros(len(targets)) for f in fams}; kd = {f: np.zeros(len(targets)) for f in fams}; cpm = {}
    for fam in fams:
        acc = np.zeros(H.NG); cnt = np.zeros(H.NG)
        for s_, f_ in TM.FAMILY.items():
            if f_ != fam or not TM._src_path(s_).exists():
                continue
            z = TM.load(s_); fc = z["fc"].astype(np.float64) * 1e6
            fin = np.isfinite(fc); acc += np.where(fin, fc, 0).sum(0); cnt += fin.sum(0)
            for k, t in enumerate(targets):
                i = z["pos"].get(H.canon(t))
                if i is None or z["ncell"][i] < 20:
                    continue
                cells[fam][k] += float(z["ncell"][i]); j = H.GIC.get(H.canon(t))
                if j is not None and np.isfinite(z["fc"][i][j]) and z["fc"][i][j] * 1e6 >= 5:
                    kd[fam][k] = max(kd[fam][k], 1.0 - (float(z["fp"][i][j]) * 1e6 + 1) / (float(z["fc"][i][j]) * 1e6 + 1))
        cpm[fam] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return cells, kd, cpm


def pieces(p):
    """Everything the sum construction and the scorer need for panel p (ml_gate_pilot.panel_eval step 2, verbatim)."""
    outer = OUTER[p]; somatic = p not in ("H1", "KOLF")
    wbase = {f: (0.0 if (somatic and f == "KOLF") else W_NEW[f]) for f in FAMS}
    for f in outer: wbase[f] = 0.0
    P0_, C0_ = H.PARTS, H.CACHE
    try:
        if somatic:
            H.PARTS = [(ww, dict(cc, weights={**cc["weights"], "KOLF": 0.0}), a) for ww, cc, a in P0_]
            H.CACHE = C0_ / "dk0"; H.CACHE.mkdir(parents=True, exist_ok=True)
        ex = H.EXCL[p]; tg = list(H.PANELS[p]["targets"]); c1 = dict(H.PARTS[1][1]); a1 = H.PARTS[1][2]
        cfg = TM.Cfg()
        for k_, v_ in c1.items(): setattr(cfg, k_, v_)
        L_ref, prom, tidx = TM.predict_lfc(tg, cfg, exclude=ex)
        ref_total, ref_ov = EMT.delta_s4(p, S4KEY[p], 0.3, 0.5, mc_list=EMT.MC7)
        lib = [None, None]
        for k, part in enumerate(H.PARTS):
            H.PARTS = [(1.0, part[1], part[2])]; lib[k] = H.base(H.LIB, ex, f"sch_lib_p{k}_{p}")[0]
        H.PARTS = P0_ if not somatic else [(ww, dict(cc, weights={**cc["weights"], "KOLF": 0.0}), a) for ww, cc, a in P0_]
        dlp = H.lib_delta(p, set(ex))
    finally:
        H.PARTS, H.CACHE = P0_, C0_
    line = EMT.S4_LINE[p]; key = S4KEY[p]
    fS = {"F": EMT.S4_EXT / line / "film_v2_a0.25.npy"}.get(key, EMT.S4_NEW / line / f"film_{key[1:]}.npy")
    S3 = 3.0 * np.load(fS).astype(np.float64)
    dlib = 0.5 * lib[0] + 0.5 * lib[1]
    for i, t in enumerate(H.LIB):
        j = H.GIC.get(H.canon(t))
        if j is not None: dlib[i, j] = 0.0
    cov = np.abs(dlib).sum(1) > 0; pos = {H.canon(t): i for i, t in enumerate(H.LIB)}
    Kn = np.zeros((len(tg), H.NG)); has = np.zeros(len(tg), bool)
    for k_, t in enumerate(tg):
        cand = [(pos[x], c) for x, c in H.NB.get(H.canon(t), []) if x in pos and x != H.canon(t) and cov[pos[x]] and c > 0][:EMT.NB_TOPK]
        if cand:
            wv = np.array([c for _, c in cand]) ** 2; Kn[k_] = (wv[:, None] * dlib[[i for i, _ in cand]]).sum(0) / wv.sum(); has[k_] = True
    fin = [f for f in FAMS if f not in outer and wbase[f] > 0]
    r = H.real(p); ok = r["gmap"] >= 0; gpan = r["gmap"][ok]
    cpm_dst = np.full(H.NG, np.nan); cpm_dst[gpan] = r["cpm"][ok]
    genes = r["genes"]; ctrl = r["ctrl"]; n = H.PANELS[p]["n_cells"]
    Dl = np.zeros((dlp.shape[0], len(genes))); Dl[:, ok] = dlp[:, gpan]; Dl = Dl[np.abs(Dl).sum(1) > 0]
    _, _, V3 = randomized_svd(Dl - Dl.mean(0), n_components=3, random_state=0); del Dl, dlp
    return dict(p=p, tg=tg, ex=ex, cfg=cfg, prom=prom, tidx=tidx, a1=a1, ref_total=ref_total, ref_ov=ref_ov, S3=S3, Kn=Kn, has=has,
                fin=fin, wb=np.array([wbase[f] for f in fin], float), r=r, ok=ok, gpan=gpan, cpm_dst=cpm_dst, genes=genes, ctrl=ctrl, n=n, V3=V3)


def build(pc, yh_members):
    """yh_members: list over members of fused LFC (T x NG, 0 where no donor). Exact ml_gate_pilot.build path."""
    cfg = pc["cfg"]; d1s = []; ov = None
    for yh in yh_members:
        L = np.clip(yh, -cfg.clip, cfg.clip)
        dd, oo = TM.split(L, pc["prom"], pc["tidx"], pc["a1"], cfg); d1s.append(dd); ov = oo if ov is None else ov
    d1 = np.mean(d1s, 0); total = 0.5 * pc["S3"] + 0.5 * d1; e0 = float((total ** 2).sum()); w = 0.3 * pc["has"]
    total = (1 - w)[:, None] * total + w[:, None] * pc["Kn"]
    return total * np.sqrt(e0 / max(float((total ** 2).sum()), 1e-12)), ov, d1


def fixed_fusion(x, m, wb):
    """x (M, T, G, F) donor LFC (0 where missing), m mask -> the incumbent fixed-weight fusion per member."""
    out = []
    for k in range(x.shape[0]):
        a = wb * m[k]; A = a.sum(-1); out.append(np.where(A > 0, (a * x[k]).sum(-1) / np.maximum(A, 1e-300), 0.0))
    return out


def evaluate_pt(panel, targets, genes, ctrl_counts, L, OV, n):
    """sch_sumeval.evaluate with per-target outputs (num, dist, pds_t); identical arithmetic."""
    Rm = SE.real_moments(panel); rp = list(Rm["perts"]); c = "non-targeting"
    ir = np.array([rp.index(t) for t in targets]); mr = Rm["means"][ir]; cr = Rm["jk"][ir]
    m0 = Rm["means"][rp.index(c)]; c0 = float(Rm["jk"][rp.index(c)])
    fm = np.asarray(ctrl_counts.sum(0)).ravel() / ctrl_counts.shape[0]
    e = np.exp(L); j = np.isfinite(OV); e[j] = np.exp(OV[j])
    mp = SE.bulk_lognorm_means(np.round(fm[None, :] * e * n), SE.TS)
    gi = {g: k for k, g in enumerate(genes)}
    own = np.zeros(mp.shape, bool)
    for k, t in enumerate(targets):
        if t in gi: own[k, gi[t]] = True
    sq = np.where(own, 0.0, (mp - mr) ** 2).sum(1)
    claim = np.minimum(cr, cr)
    mpc = np.where(own, np.nan, mp); B = float(np.nansum((mpc - np.nanmean(mpc, 0)) ** 2))
    r = min(1.0, B / claim.sum())
    num = sq - r * claim - cr
    dist = np.where(own, 0.0, (mr - m0[None, :]) ** 2).sum(1) - cr - c0
    panel_cols = np.array(sorted({gi[t] for t in targets if t in gi}))
    keep = np.ones(len(genes), bool); keep[panel_cols] = False
    P_ = (mp - m0[None, :])[:, keep]; R_ = (mr - m0[None, :])[:, keep]
    Pn = P_ / (np.linalg.norm(P_, axis=1, keepdims=True) + 1e-12); Rn = R_ / (np.linalg.norm(R_, axis=1, keepdims=True) + 1e-12)
    Dm = 1.0 - Pn @ Rn.T; nT = len(targets); ranks = np.empty(nT)
    for k in range(nT):
        row = Dm[k]; x = row[k]; less = (row < x).sum(); eq = (row == x).sum(); ranks[k] = less + (eq - 1) / 2.0
    pds_t = 1.0 - ranks / nT
    return dict(num=num, dist=dist, pds_t=pds_t, mse=float(num.sum() / dist.sum()), pds=float(pds_t.mean()), r=r)


class Scorer:
    """ct 0.5 + guarded pcr + budget match, budget FIXED from the base total (ml_gate_pilot score), then evaluate_pt."""
    def __init__(self, pc, base_total, base_ov):
        self.pc = pc; p = pc["p"]; r = pc["r"]; ok = pc["ok"]
        L0 = H.on_panel(base_total, p); Od = np.full(L0.shape, np.nan); Od[:, ok] = base_ov[:, pc["gpan"]]
        gk, _ = EMT.kink_gamma(pc["ctrl"], L0, Od, pc["n"], factor=1.17, margin=0.8)
        self.Bt = EM.budget_B(pc["ctrl"], L0, Od, pc["n"])(max(1.0, gk)); self.Od = Od
    def transform(self, total):
        pc = self.pc; L = H.on_panel(total, pc["p"]); mt = np.nanmean(np.where(np.isfinite(self.Od), np.nan, L), 0)
        X = L - 0.5 * np.nan_to_num(mt)[None, :]; P_ = X @ pc["V3"].T; f_ = float((P_ ** 2).sum() / (X ** 2).sum())
        X = X - 0.75 * float(np.clip((0.25 - f_) / 0.10, 0, 1)) * P_ @ pc["V3"]
        return X * EM.budget_match(pc["ctrl"], X, self.Od, pc["n"], self.Bt)
    def score(self, total):
        pc = self.pc
        return evaluate_pt(pc["p"], [str(t) for t in pc["tg"]], pc["genes"], pc["ctrl"], self.transform(total), self.Od, pc["n"])


def assemble(p):
    pc = pieces(p); tg = pc["tg"]; fin = pc["fin"]; log(p, "pieces ready; donors", fin, "weights", pc["wb"].tolist())
    LIB = list(H.LIB); FL = fam_means(LIB, fin); LC, LKD, CPM = fam_stats(LIB, fin)
    SD = {}
    for f in fin:
        v = FL[f]; med = np.nanmedian(v, 0); SD[f] = np.maximum(np.nan_to_num(np.nanmedian(np.abs(v - med[None, :]), 0) * 1.4826, nan=1.0), 1e-3)
    pi = np.array([LIB.index(t) if t in LIB else -1 for t in tg])
    cells_t = np.stack([np.where(pi >= 0, LC[f][np.maximum(pi, 0)], 0.0) for f in fin]); kd_t = np.stack([np.where(pi >= 0, LKD[f][np.maximum(pi, 0)], 0.0) for f in fin])
    log(p, "library stats ready")
    X = np.zeros((len(MC), len(tg), H.NG, len(fin)), np.float32)   # float32: float16 broke the 1e-4 parity (4e-4)
    Mk = np.zeros((len(MC), len(tg), H.NG, len(fin)), bool)
    for mi, mc in enumerate(MC):
        Xm = fam_means(tg, fin, 20 if mc is None else mc)
        for fi, f in enumerate(fin):
            v = Xm.get(f, np.full((len(tg), H.NG), np.nan)); Mk[mi, :, :, fi] = np.isfinite(v); X[mi, :, :, fi] = np.nan_to_num(v).astype(np.float32)
        log(p, f"member {mi} (min_cells {20 if mc is None else mc}) ready; coverage per donor", Mk[mi].any(1).mean(0).round(2).tolist())
    # truth on the official axis; own gene + promoter neighbours masked; destination gate cpm >= 1
    r = pc["r"]; y = np.full((len(tg), H.NG), np.nan); y[:, pc["gpan"]] = r["R"][:, pc["ok"]]
    gate = np.isfinite(y) & (np.nan_to_num(pc["cpm_dst"]) >= 1.0)[None, :]
    for k, t in enumerate(tg):
        j = H.GIC.get(H.canon(t))
        if j is not None: gate[k, j] = False
        for jj in (pc["prom"][k] if isinstance(pc["prom"], (list, tuple)) and k < len(pc["prom"]) else []):
            try: gate[k, int(jj)] = False
            except Exception: pass
    # base parity
    base_members = fixed_fusion(X.astype(np.float32), Mk, pc["wb"]); tot0, ov0, _ = build(pc, base_members)
    par = float(np.abs(tot0 - pc["ref_total"]).max()); log(p, f"parity fixed-weight fusion vs delta_s4: {par:.2e}"); assert par < 1e-4, par
    np.savez(OUTD / f"{p}_data.npz", x=X, m=Mk, y=np.nan_to_num(y).astype(np.float32), gate=gate, wb=pc["wb"], fin=np.array(fin),
             cpm_fam=np.stack([np.nan_to_num(CPM[f]) for f in fin]).astype(np.float32), cpm_dst=np.nan_to_num(pc["cpm_dst"]).astype(np.float32),
             cells_t=cells_t.astype(np.float32), kd_t=kd_t.astype(np.float32), sd=np.stack([SD[f] for f in fin]).astype(np.float32),
             tg=np.array([str(t) for t in tg]), mc=np.array([20 if mc is None else mc for mc in MC]))
    log(p, "DATA_SAVED", OUTD / f"{p}_data.npz")


def evaluate(p, oof_file, tag):
    pc = pieces(p); tg = pc["tg"]; d = np.load(OUTD / f"{p}_data.npz"); X = d["x"].astype(np.float32); Mk = d["m"]; wb = d["wb"]
    base_members = fixed_fusion(X, Mk, wb); tot0, ov0, d1_0 = build(pc, base_members)
    par = float(np.abs(tot0 - pc["ref_total"]).max()); log(p, f"parity base vs delta_s4: {par:.2e}"); assert par < 1e-4
    sc = Scorer(pc, tot0, ov0); base = sc.score(tot0); log(p, f"BASE offline pds {base['pds']:.4f} mse {base['mse']:.4f} (harness real-scorer base: see result.json)")
    o = np.load(oof_file); yhat = o["yhat"]; folds = o["fold"]; K = yhat.shape[0]
    num = np.zeros(len(tg)); pds_t = np.zeros(len(tg)); num_f = np.zeros(len(tg)); pds_f = np.zeros(len(tg)); dist = base["dist"]
    for k in range(K):
        mem = [yhat[k, mi].astype(np.float32) for mi in range(yhat.shape[1])]
        tot, ov, d1 = build(pc, mem); s = sc.score(tot); rows = np.flatnonzero(folds == k)
        num[rows] = s["num"][rows]; pds_t[rows] = s["pds_t"][rows]
        # fusion-only diagnostic: part 1 alone (energy matched to the base's part 1), same transforms / budget
        tf = d1 * np.sqrt(float((d1_0 ** 2).sum()) / max(float((d1 ** 2).sum()), 1e-12)); sf = sc.score(tf)
        num_f[rows] = sf["num"][rows]; pds_f[rows] = sf["pds_t"][rows]
        log(p, f"fold {k}: full-panel pds {s['pds']:.4f} mse {s['mse']:.4f} r {s['r']:.3f}; kept {len(rows)} targets")
    base_f = sc.score(d1_0)
    oof = dict(pds=float(pds_t.mean()), mse=float(num.sum() / dist.sum())); oof_f = dict(pds=float(pds_f.mean()), mse=float(num_f.sum() / dist.sum()))
    J = lambda a, b: 2.23 * (a["pds"] - b["pds"]) - 1.05 * (a["mse"] - b["mse"])
    # oracle: per-target NNLS family weights on the panel truth (panel axis, panel-target genes excluded), member 0 inputs
    r = pc["r"]; ok = pc["ok"]; gpan = pc["gpan"]; genes = pc["genes"]; tset = set(str(t) for t in tg); keepg = np.array([g not in tset for g in genes])
    Wt = np.zeros((len(tg), len(pc["fin"])))
    for t in range(len(tg)):
        A = X[0, t][gpan][:, :][keepg[ok]]; yv = np.nan_to_num(r["R"][t][ok])[keepg[ok]]; mm = Mk[0, t].sum(0) > 100
        if mm.sum() and np.linalg.norm(yv) > 0:
            wv, _ = nnls(A[:, mm], yv); Wt[t, np.flatnonzero(mm)] = wv
        if Wt[t].sum() <= 0: Wt[t] = wb * mm
    orc_members = []
    for mi in range(X.shape[0]):
        a = Wt[:, None, :] * Mk[mi]; A_ = a.sum(-1); orc_members.append(np.where(A_ > 0, (a * X[mi]).sum(-1) / np.maximum(A_, 1e-300), 0.0))
    orc = sc.score(build(pc, orc_members)[0])
    rec = dict(panel=p, tag=tag, base=dict(pds=base["pds"], mse=base["mse"]), oof=oof, J=J(oof, base), dpds=oof["pds"] - base["pds"], dmse=oof["mse"] - base["mse"],
               fusion_only=dict(base=dict(pds=base_f["pds"], mse=base_f["mse"]), oof=oof_f, J=J(oof_f, base_f)),
               oracle=dict(pds=orc["pds"], mse=orc["mse"], J=J(orc, base)), folds=int(K), pds_t_oof=pds_t.tolist(), pds_t_base=base["pds_t"].tolist())
    verdict = "KILL" if rec["J"] < 0.03 else ("INCONCLUSIVE" if rec["J"] < 0.06 or rec["dpds"] <= 0 else "PROCEED")
    rec["verdict"] = verdict
    (OUTD / f"{p}_eval_{tag}.json").write_text(json.dumps(rec, default=float))
    log(p, f"OOF full recipe: pds {oof['pds']:.4f} (base {base['pds']:.4f}, d {rec['dpds']:+.4f}) mse {oof['mse']:.4f} (base {base['mse']:.4f}, d {rec['dmse']:+.4f}) J {rec['J']:+.4f} -> {verdict}")
    log(p, f"fusion-only: base pds {base_f['pds']:.4f} mse {base_f['mse']:.4f} -> oof pds {oof_f['pds']:.4f} mse {oof_f['mse']:.4f} J {rec['fusion_only']['J']:+.4f}")
    log(p, f"oracle (per-target NNLS, reference): pds {orc['pds']:.4f} mse {orc['mse']:.4f} J {rec['oracle']['J']:+.4f}")
    log("EVAL_DONE")


if __name__ == "__main__":
    mode = _argv[1]
    if mode == "assemble": assemble(_argv[2])
    elif mode == "eval": evaluate(_argv[2], _argv[3], _argv[4] if len(_argv) > 4 else "v1")
