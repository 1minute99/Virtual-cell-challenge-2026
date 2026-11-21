"""Emission optimisation (2026-09-25): separate the two effects of dual-moment emission.
Dual (fh_dual) = pinned multiplicative group totals (-> MSE credit, PDS) + top-K within-column re-allocation (-> DE).
The harness default shared_cellshift_gated = log1p(CP10K) additive shift on one shared control block: its amplitude
on 5-50 CPM genes is what gives NMAE/reach, but its group sums carry the shared block's sampling noise (MSE 2.4).

Arm name  em_<target>_<shape>_k<K>[_strat][_du<a>][_ag][_g<gamma>][@s<seed>]
  target  mix  : log T = (1-b) log T_mult + b log T_cs   (b from _b<beta>; expression-dependent amplitude)
          mult : T_g = n * control mean_g * exp(gamma d_g)                     (fh_dual)
          cs   : T_g = n * E_{all controls}[ shifted count ]  -- the cellshift emission's expected group total
                 computed over ALL controls (exact, chunked), i.e. its amplitude without the block noise
  shape   share: allocation proportional to the shared block's counts          (fh_dual)
          cs   : allocation proportional to the cellshift-shifted shared block
  K       top-K re-allocation of fh_dual (0 = none); du/ag as in fh_dual
  strat   the shared block drawn phase-stratified (G1/S/G2M in control proportions, sch_arms pools)
Low genes (< 5 CPM) and override entries are multiplicative in every arm (as in both parents).
Parity: em_mult_share_k3000 == fh_dual dual_g1_k3000 (same RNG stream).
Usage: sch_emit.py ARM ... [panel=P ...]  -> reports/sc_hybrid_2026-09-25/emit/<arm>/<panel>/result.json"""
import json, os, pathlib, subprocess, sys, time, shutil, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2.moments import jackknife_correction  # noqa: F401
_argv = list(sys.argv); sys.argv = ["x"]
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H
import fh_dual as FD

OUT = H.R / "reports/sc_hybrid_2026-09-25/emit"; PH = ["G1", "S", "G2M"]


DELTAS = {"Dunseen": "unseen_knn", "Dbase": "base", "Dblend20e": "blend20_e", "Dblend40e": "blend40_e", "Dblend50e": "blend50_e",
          "Dsimsoft": "sim_soft", "Dcombo": "combo"}


_SNB = {}
def string_nb():
    if not _SNB:
        _SNB.update(json.loads((OUT.parent / "string_neighbours.json").read_text()))
    return _SNB


S4_LINE = {"KOLF": "KOLF", "RN_K562": "K562", "RN_HepG2": "HepG2", "RN_Jurkat": "Jurkat", "RN_RPE1": "RPE1", "H1": "H1", "OR_HCT116": "HCT116", "OR_HCT116X": "HCT116"}
S4_EXT = H.R / "reports/ext_contexts_2026-09-21"; S4_NEW = H.R / "reports/s4_improve_2026-09-27"


def delta_s4(panel, key, bw=0.3, share=0.5, mc_list=None, c_list=None, wmul=None, fwins=0.0, crob=0.0, wadd=None):
    """Official-builder mirror with a LEARNED part 0: total = 0.5 x 3 x S4 + 0.5 x (1.6 x fusion); DepMap-neighbour blend
    bw (library = the harness parts, as knn_official), GLOBAL energy match (no per-target matching); override from the
    fusion part (promoter floor + knockdown). S4 = the saved leave-one-biology-out prediction of the panel's line:
    key F = film_v2_a0.25 (the S4 recipe), R = refit_v2 (refit only), X<tag> / Y<tag> = reports/s4_improve_2026-09-27/<line>/
    film_<tag>.npy with (X) or without (Y) per-target energy matching to film_v2_a0.25."""
    line = S4_LINE[panel]; tg = H.PANELS[panel]["targets"]; ex = H.EXCL[panel]; P0 = H.PARTS
    wsfx = ("_wm" + "_".join(f"{x}{y:g}" for x, y in sorted(wmul.items()))) if wmul else ""
    if fwins: wsfx += f"_fw{fwins:g}"          # cache key must separate winsorized builds from plain ones
    if crob:  wsfx += f"_cr{crob:g}"
    if wadd: wsfx += "_wa" + "_".join(f"{x}{y:g}" for x, y in sorted(wadd.items()))
    _ov = {}
    if fwins: _ov["fwins"] = fwins
    if crob:  _ov["center_robust"] = crob
    if wadd:  P0 = [(w_, dict(c_, weights={**c_["weights"], **wadd}), a_) for w_, c_, a_ in P0]
    P0 = [(w_, dict(c_, **_ov), a_) for w_, c_, a_ in P0] if _ov else P0
    f = {"F": S4_EXT / line / "film_v2_a0.25.npy", "R": S4_EXT / line / "refit_v2.npy"}.get(key, S4_NEW / line / f"film_{key[1:]}.npy")
    S = np.load(f).astype(np.float64); assert S.shape == (len(tg), H.NG), (S.shape, len(tg))
    if key.startswith("X"):                        # X<tag>: per-target energy matched to the original S4 recipe (direction only)
        S0 = np.load(S4_EXT / line / "film_v2_a0.25.npy").astype(np.float64)
        n0 = np.linalg.norm(S0, axis=1); n1 = np.linalg.norm(S, axis=1); S = S * np.where(n1 > 0, n0 / np.maximum(n1, 1e-12), 1.0)[:, None]
    try:
        if mc_list or c_list:                      # fusion part ensembled over min_cells (ml_hp_ensemble6.py arm A)
            d1s, ov = [], None                     # and/or over the LFC PSEUDOCOUNT c (new 2026-10-03, token `cens`):
            mem = [(None, None)] + [("min_cells", v) for v in (mc_list or ())] + [("c", v) for v in (c_list or ())]
            for fld, v in mem:                     # union, not cross-product -- mirrors ENS = _MC5 + _W6 + _C3
                cfg1 = dict(P0[1][1]) if fld is None else dict(P0[1][1], **{fld: v})
                H.PARTS = [(1.0, cfg1, P0[1][2])]
                d_, ov_ = H.base(tg, ex, f"sch_p1_{panel}" + wsfx + ("" if fld is None else f"_cx{fld}{v:g}"), wmul=wmul)[:2]
                d1s.append(d_)
                if ov is None: ov = ov_
            d1 = np.mean(d1s, axis=0)
        else:
            H.PARTS = [(1.0, P0[1][1], P0[1][2])]; d1, ov = H.base(tg, ex, f"sch_p1_{panel}" + wsfx, wmul=wmul)
        lib = [None, None]
        for k, part in enumerate(P0):
            H.PARTS = [(1.0, part[1], part[2])]; lib[k] = H.base(H.LIB, ex, f"sch_lib_p{k}_{panel}" + wsfx, wmul=wmul)[0]
    finally:
        H.PARTS = P0
    total = share * 3.0 * S + (1 - share) * d1        # official S10c/S10i: share 0.5 (S10j sum side: 0.75)
    dlib = 0.5 * lib[0] + 0.5 * lib[1]
    for i, t in enumerate(H.LIB):
        j = H.GIC.get(H.canon(t))
        if j is not None:
            dlib[i, j] = 0.0
    cov = np.abs(dlib).sum(1) > 0; pos = {H.canon(t): i for i, t in enumerate(H.LIB)}
    K = np.zeros_like(total); has = np.zeros(len(tg), bool)
    KL_ = NB_TOPK_ENS or [NB_TOPK]            # averaging over the cut is a NEW diversity axis, like the min_cells ensemble
    for k_, t in enumerate(tg):
        allc_ = [(pos[x], c) for x, c in H.NB.get(H.canon(t), []) if x in pos and x != H.canon(t) and cov[pos[x]] and c > 0]
        acc_ = []
        for kk_ in KL_:
            cand = allc_[:kk_]
            if cand:
                wv = np.array([c for _, c in cand]) ** 2; acc_.append((wv[:, None] * dlib[[i for i, _ in cand]]).sum(0) / wv.sum())
        if acc_:
            K[k_] = np.mean(acc_, 0); has[k_] = True
    e0 = float((total ** 2).sum()); w = bw * has
    total = (1 - w)[:, None] * total + w[:, None] * K; total *= np.sqrt(e0 / max(float((total ** 2).sum()), 1e-12))
    return total, ov


def _kv(b):
    """cache-key text for a cfg override value (numbers as %g, anything else as a short stable hash)"""
    if isinstance(b, (int, float)):
        return f"{b:g}"
    import hashlib
    return hashlib.md5(json.dumps(b, sort_keys=True).encode()).hexdigest()[:8]


# HYPERPARAMETER-ENSEMBLE member pools (ml_hp_ensemble*.py). Each member is (fs, bw, nb, wmul, cfgx); the ensemble
# delta is their plain MEAN, which -- unlike any re-weighting of ONE prediction -- improves PDS and MSE together
# because the members carry OPPOSITE biases (measured: every `w` member is individually J < 0 on H1, their mean is J > 0).
_MC = lambda v: {0: {"min_cells": v}, 1: {"min_cells": v}}
_W6 = [(0.5, 0.3, "D", {f: m}, None) for f, m in
       (("K562", 0.25), ("K562", 4.0), ("HCT116", 4.0), ("HEK293T", 4.0), ("CD4", 4.0), ("FENG", 4.0))]
_MC5 = [(0.5, 0.3, "D", None, None)] + [(0.5, 0.3, "D", None, _MC(v)) for v in (5, 10, 50, 100)]
_MC7 = [(0.5, 0.3, "D", None, None)] + [(0.5, 0.3, "D", None, _MC(v)) for v in (5, 10, 15, 30, 50, 100, 200)]
_C3 = [(0.5, 0.3, "D", None, {1: {"c": c}}) for c in (10.0, 20.0, 100.0)]
MC7 = (5, 10, 15, 30, 50, 100, 200)
NB_TOPK = 10          # neighbour-library cut: another discontinuous SELECTION boundary, ensemble-able like min_cells
CENS_OK = (1.0, 5.0, 10.0, 20.0, 50.0, 100.0)   # tv2_sources.CS: the pseudocounts with precomputed `center_c<c>`
NB_TOPK_ENS = []      # when non-empty, the neighbour estimate is AVERAGED over these cuts (token `nbe5.10.20.40`)      # the min_cells pool of the winning ensemble
ENS = {"mc": _MC5, "mc7": _MC7, "mcw": _MC5 + _W6, "mcwc": _MC5 + _W6 + _C3}


def delta_ens(panel, name):
    """Mean of the ENS[name] members' delta_fs outputs; the override comes from the first (base) member."""
    ds, ov0 = [], None
    for fs, bw, nb, wm, cx in ENS[name]:
        d, ov = delta_fs(panel, fs, bw, nb=nb, wmul=wm, cfgx=cx, cfgx_lib=cx)
        ds.append(d)
        if ov0 is None: ov0 = ov
    return np.mean(ds, axis=0), ov0


def delta_fs(panel, fs, bw, pt=False, nb="D", wmul=None, cfgx=None, cfgx_lib=None):
    """blend<bw>_e recipe with the two PARTS mixed (1-fs) : fs, energy-matched to the default recipe (fs 0.5, same bw).
    NOTE (2026-09-27): the loop variable used to shadow `pt`, so the per-target energy match below ALWAYS ran; every
    delta_fs arm scored so far is per-target matched. Kept that way (pt forced True) so old and new arms compare.
    wmul: {family: multiplier} on both parts' source weights (targets and neighbour library alike)."""
    tg = H.PANELS[panel]["targets"]; ex = H.EXCL[panel]; P0 = H.PARTS; pt = True
    sfx = ("_wm" + "_".join(f"{a}{b:g}" for a, b in sorted(wmul.items()))) if wmul else ""
    parts, lib = [], []
    try:
        for k, part in enumerate(P0):
            cx = (cfgx or {}).get(k) or {}
            sk = sfx + ("_cx" + "_".join(f"{a}{_kv(b)}" for a, b in sorted(cx.items())) if cx else "")
            H.PARTS = [(1.0, dict(part[1], **cx), part[2])]
            parts.append(H.base(tg, ex, f"sch_p{k}_{panel}{sk}", wmul=wmul))
            cl = cx if cfgx_lib is None else ((cfgx_lib or {}).get(k) or {})   # neighbour library may use its own cfg
            skl = sfx + ("_cx" + "_".join(f"{a}{_kv(b)}" for a, b in sorted(cl.items())) if cl else "")
            H.PARTS = [(1.0, dict(part[1], **cl), part[2])]
            lib.append(H.base(H.LIB, ex, f"sch_lib_p{k}_{panel}{skl}", wmul=wmul)[0])
    finally:
        H.PARTS = P0
    def mix(a, b): return (1 - fs) * a + fs * b
    d = mix(parts[0][0], parts[1][0]); ov = np.where(np.isfinite(parts[0][1]), parts[0][1], parts[1][1])
    dlib = mix(lib[0], lib[1])
    for i, t in enumerate(H.LIB):
        j = H.GIC.get(H.canon(t))
        if j is not None:
            dlib[i, j] = 0.0
    cov = np.abs(dlib).sum(1) > 0; pos = {H.canon(t): i for i, t in enumerate(H.LIB)}
    def knn_of(NBD):
        kn = np.zeros((len(tg), H.NG)); hs = np.zeros(len(tg), bool)
        for k_, t in enumerate(tg):
            cand = [(pos[x], c) for x, c in NBD.get(H.canon(t), []) if x in pos and x != H.canon(t) and cov[pos[x]] and c > 0][:NB_TOPK]
            if cand:
                wv = np.array([c for _, c in cand]) ** 2
                kn[k_] = (wv[:, None] * dlib[[i for i, _ in cand]]).sum(0) / wv.sum(); hs[k_] = True
        return kn, hs
    if nb == "D":
        knn, has = knn_of(H.NB)
    elif nb.startswith("P"):
        knn, has = knn_of(json.loads((OUT.parent / f"prior_neighbours_{nb[1:]}.json").read_text()))
    elif nb.startswith("C"):
        knn, has = knn_of(json.loads((OUT.parent / f"concat_neighbours_a{float(nb[1:]):g}.json").read_text()))
    elif nb == "S":
        knn, has = knn_of(string_nb())
    else:                                          # DS: average of the two estimates where both exist
        k1, h1 = knn_of(H.NB); k2, h2 = knn_of(string_nb())
        both = h1 & h2
        knn = np.where(both[:, None], 0.5 * (k1 + k2), np.where(h1[:, None], k1, k2)); has = h1 | h2
    w = bw * has; dn = H.ematch((1 - w)[:, None] * d + w[:, None] * knn, d, panel)
    ref, _, _ = H.arm_delta(f"blend{int(round(bw * 100))}_e", panel) if abs(bw - 0.2) > 1e-9 else H.arm_delta("blend20_e", panel)
    if pt:                                         # PDS is per-target scale-free, MSE is not: keep the default per-target energy
        r = H.real(panel); ok = r["gmap"] >= 0; g = r["gmap"][ok]
        nr = np.linalg.norm(ref[:, g], axis=1); nd = np.linalg.norm(dn[:, g], axis=1)
        return dn * np.where(nd > 0, nr / np.maximum(nd, 1e-12), 1.0)[:, None], ov
    return H.ematch(dn, ref, panel), ov


def parse(arm):
    name, _, sd = arm.partition("@s"); f = name.split("_")
    o = dict(target=f[1], shape=f[2], K=int(f[3][1:]), strat=False, du=0.0, ag=False, gamma=1.0, seed=int(sd) if sd else 0)
    for x in f[4:]:
        if x.startswith("ens"): o["ens"] = x[3:]        # hyperparameter ensemble, see ENS above
        elif x == "lsb": o["lsb"] = True                        # depth-balanced emitted block, see pick_block
        elif x == "bmb": o["bmb"] = True                        # budget-match to the NON-ensembled recipe (see below)
        elif x.startswith("crob"):                            # robust (trimmed-mean) per-source center: crob0.5
            o["crob"] = float(x[4:])
        elif x.startswith("fw") and x[2:].replace(".", "").isdigit():   # family-level winsorize: fw0.1
            o["fwins"] = float(x[2:])
        elif x.startswith("cens"):                            # ENSEMBLE the LFC pseudocount: cens5.10.20.100
            o["cens"] = [float(v) for v in x[4:].split(".") if v]
            bad_ = [v for v in o["cens"] if v not in CENS_OK]  # only these have precomputed centring vectors (tv2_sources.CS)
            assert not bad_, f"cens: c must be one of {CENS_OK}, got {bad_}"
        elif x.startswith("nbe"):                             # ENSEMBLE the neighbour-library cut: nbe5.10.20.40
            o["nbe"] = [int(v) for v in x[3:].split(".") if v]
        elif x.startswith("tzb"):                             # SHAPE-only variant: applied BEFORE the budget match, so
            if "x" in x[3:]:                                  # `bm` restores B_pred -> budget-INVARIANT, offline price valid
                g6_, d6_ = x[3:].split("x"); o["tzb"] = (float(g6_), float(d6_))   # ...x<a>: DE side follows at w**a
            else:
                o["tzb"] = (float(x[3:]), 0.0)
        elif x.startswith("tz") and "x" in x:                 # 2nd per-target amplitude AFTER the budget match: tz<gamma>x<gscale>
            a5_, b5_ = x[2:].split("x"); o["tz"] = (float(a5_), float(b5_))
        elif x.startswith("dis") and "x" in x:                 # SUM-side per-target HUBNESS reduction dis<k>x<lam>
            k4_, l4_ = x[3:].split("x"); o["dis"] = (int(k4_), float(l4_))
        elif x == "strat": o["strat"] = True
        elif x == "ag": o["ag"] = True
        elif x.startswith("du"): o["du"] = float(x[2:])
        elif x == "pt": o["pt"] = True
        elif x == "kp": pass
        elif x.startswith("sbw"): o["sbw"] = float(x[3:]) / 100               # SUM-side-only neighbour blend weight (DE side: 0.3)
        elif x == "tlrn": o["tlrn"] = True       # 2026-10-02: ta rank = LOPO-learned true-cosine prediction (ml_amp_learn.py)
        elif x.startswith("drp"): o["drp"] = float(x[3:])   # 2026-10-03: continuous DE rank calibration, (rank/N)^gamma, median-normalised
        elif x.startswith("drc") and "x" in x[3:]:          # 2026-10-03: DE-side RANK-CLASS amplitude calibration, sums untouched
            a_, b_ = x[3:].split("x"); o["drc"] = (float(a_), float(b_))   # top 5 % of genes by |DE lfc| x a, bottom 70 % x b
        elif x == "tcomp": o["tcomp"] = True     # 2026-10-02: ta rank = rank(rank(agree) + rank(#families) + rank(|delta|))
        elif x.startswith("tad") and x[3:].replace(".", "", 1).isdigit(): o["tad"] = float(x[3:])   # 2026-10-02: same rank^beta on the DE side
        elif x.startswith("ta") and x[2:].replace(".", "", 1).isdigit(): o["ta"] = float(x[2:])   # 2026-10-01: per-target SUM amplitude x rank(cross-family agreement)^beta
        elif x == "bm": o["bm"] = True                                           # budget-matched sums (with gk), see tv2_emit.budget_match
        elif x.startswith("dpj") and "x" in x:                              # DE-side projection onto the top-K library programs
            k_, r_ = x[3:].split("x"); o["dpj"] = [int(k_), float(r_)]
        elif x.startswith("xw") and x[2:].replace(".", "").isdigit(): o["xw"] = float(x[2:])   # Cfg.xw on the recipe (sums AND DE side),
                                                                 # not on the budget reference (brs / bmb); tv2_model.Cfg.xw
        elif x.startswith("dsrc"): o["dsrc"] = x[4:]
        elif x.startswith("dsl") and not x.startswith("dsw"): o["dslib"] = x[3:]   # dsrc arm: neighbour LIBRARY from this other variant dir
                                                                 # (mirrors the shipped S14fk DE side: K562-only fusion, full library)                  # DE side from a VARIANT source library (denoise_src/<tag>)
        elif x == "agc": o["agc"] = True                        # DE ranking = |DE delta| x agreement with the FAMILY-MEAN sign (builder)
        elif x.startswith("sls") and "x" in x[3:]:                         # sls<SRC>x<w>: same-line replacement on the SUMS only
            s_, w_ = x[3:].rsplit("x", 1); o["sl"] = [s_, float(w_)]; o["sl_sum"] = True
        elif x.startswith("sl") and "x" in x[2:]:                          # sl<SRC>x<w>: same-line replacement (builder same_line)
            s_, w_ = x[2:].rsplit("x", 1); o["sl"] = [s_, float(w_)]
        elif x == "dk0": o["dk0"] = True                        # base delta WITHOUT KOLF in the fusion + neighbour library (S10k sum side)
        elif x.startswith("deo"): o["deo"] = x[3:]              # DE side + budget ref + agreement from the OLD main recipe (S4 key, KOLF on, no ensemble)
        elif x.startswith("bro"): o["bro"] = x[3:]              # budget ref + agreement from the OLD main recipe; DE side = the base delta
        elif x.startswith("brs"): o["brs"] = x[3:]              # budget ref from the OLD main recipe only; agreement AND DE side from the base delta
        elif x.startswith("s4w"): o["s4w"] = float(x[3:])                        # S4 share of the parts mix (default 0.5)
        elif x.startswith("s4") and len(x) > 2: o["s4"] = x[2:]                  # part 0 = S4-type LOBO prediction (see delta_s4)
        elif x.startswith("dpr") and "x" in x:                                   # DE-side removal of the top library axes dpr<k>x<l>[g]
            g2_ = x.endswith("g"); x2_ = x[3:-1] if g2_ else x[3:]; k2_, l2_ = x2_.split("x")
            o["dpcr"] = (int(k2_), float(l2_)); o["dpcr_guard"] = (0.15, 0.25) if g2_ else None
        elif x.startswith("zs") and x[2:].replace(".", "").isdigit(): o["cfgx"] = {0: {"zshrink": float(x[2:])}, 1: {"zshrink": float(x[2:])}}
        elif x.startswith("zf") and x[2:].replace(".", "").isdigit(): o["cfgx"] = {1: {"zshrink": float(x[2:])}}   # fusion part only
        elif x.startswith("zlp") and x[3:].replace(".", "").isdigit(): o["cfgx_lib"] = {0: {"zshrink": float(x[3:])}}  # library part 0
        elif x.startswith("zl") and x[2:].replace(".", "").isdigit(): o["cfgx_lib"] = {0: {"zshrink": float(x[2:])}, 1: {"zshrink": float(x[2:])}}
        elif x.startswith("wa") and "x" in x:    # ADD a family that carries no weight yet: wa<FAM>x<weight>
            fa_, wa_ = x[2:].split("x"); o.setdefault("wadd", {})[fa_] = float(wa_)
        elif x.startswith("wm") and "x" in x:                                     # family weight multiplier wm<FAM>x<mult>
            fam_, m_ = x[2:].split("x"); o.setdefault("wmul", {})[fam_] = float(m_)
        elif x.startswith("sw"): o["sw"] = x[2:]
        elif x.startswith("ct"): o["ct"] = float(x[2:])            # SUM-side centring: minus ct x across-target mean
        elif x.startswith("ew"):                                    # SUM-side expression down-weighting ew<c0>k<kappa>
            c0_, k_ = x[2:].split("k"); o["ew"] = (float(c0_), float(k_))
        elif x.startswith("pr") and "x" in x:                        # SUM-side removal of top source-library axes pr<k>x<lam>
            g_ = x.endswith("g"); x = x[:-1] if g_ else x                            # pr<k>x<l>g: guarded (see tv2_submit_build)
            k_, l_ = x[2:].split("x"); o["pcr"] = (int(k_), float(l_)); o["pcr_guard"] = (0.15, 0.25) if g_ else None                   # SUM-side reliability weight: snr / snr05 / agree / agree05
        elif x.startswith("ci"): o["ci"] = int(x[2:])              # cpc fixed-point iterations on the emitted cells' libraries
        elif x.startswith("gk"): o["gk"] = float(x[2:]) if len(x) > 2 else 1.17   # sums at the #348 budget kink, floor 1
        elif x.startswith("mg"): o["mg"] = float(x[2:])                              # kink safety margin (target B/claim)
        elif x.startswith("gdf"): o["gdf"] = float(x[3:])                            # DE side follows gamma_sum^gdf
        elif x.startswith("dct"): o["dct"] = float(x[3:])                            # DE-side centring (cpc target, top-K)
        elif x.startswith("lcp"): o["lcp"] = float(x[3:])                            # re-allocation strength x (median cpm / cpm)^eta
        elif x.startswith("dsw"): o["dsw"] = float(x[3:])                            # DE-side weight snr^dsw (cpc target, top-K)
        elif x in ("so", "dso", "do"): o["tscale"] = x               # per-target ORACLE slope: sums only / both / DE side only                                     # no-op tag (a fresh result dir, e.g. to keep DE tables)
        elif x in ("nbS", "nbDS"): o["nb"] = x[2:]
        elif x.startswith("nbC"): o["nb"] = x[2:]
        elif x.startswith("nbP"): o["nb"] = x[2:]                # alternative prior neighbours (sch_prior_nb.py)                # concatenated-prior neighbours, STRING weight a               # neighbour source: STRING, or DepMap+STRING average
        elif x.startswith("gt"): o["gt"] = float(x[2:])          # expression gate (CPM) for re-allocation / cpc (5)
        elif x.startswith("pc"): o["pc"] = float(x[2:])          # pseudocount of the cellshift DE target (log1p: 1)                         # per-target energy matched to the default recipe
        elif x.startswith("fs"): o["fs"] = float(x[2:]) / 100     # fusion share of the two recipe parts (0.5)
        elif x.startswith("bw"): o["bw"] = float(x[2:]) / 100     # DepMap-neighbour blend weight (0.3)
        elif x.startswith("nt"): o["nt"] = float(x[2:])          # fraction of cells down-genes are concentrated into (0.1)
        elif x.startswith("le"): o["le"] = float(x[2:])          # re-allocation strength at the last rank (0.3)
        elif x.startswith("gd"): o["gd"] = float(x[2:])          # DE-side (cpc target) amplitude, absolute on the delta
        elif x.startswith("D"): o["delta"] = DELTAS[x]            # delta recipe (fh_harness arm) instead of blend30_e
        elif x.startswith("g"): o["gamma"] = float(x[1:])
        elif x.startswith("cb"): o["cb"] = float(x[2:])
        elif x == "cp": o["cp"] = True
        elif x in ("cpm", "cpc"): o["cp"] = True; o["cpt"] = x; o["cpall"] = True     # DE-functional pin on ALL gated genes
        elif x.startswith("b"): o["beta"] = float(x[1:])
        else: raise ValueError(x)
    return o


def pick_block(ctrl, n, rng, strat, panel, bal=None):
    """bal: choose the emitted block to MATCH the full control pool instead of sampling it at random.

    The scorer's DE log2FC is a per-cell normalize-then-average against the REAL control pool, so whatever the 400-cell
    block's own mean composition is, its difference from the pool's enters every gene's predicted LFC as a bias. Simple
    random sampling leaves that difference at ~1/sqrt(400) of the per-cell spread. `lsb` removes the part of it that is
    explained by sequencing depth (systematic sampling on the library-size order, which is unbiased and much
    lower-variance for anything correlated with depth -- and #286/§1.2 name depth-composition covariance as exactly what
    decides whether the emitted and scored functionals agree); `lsbg` then greedily swaps cells to reduce the residual
    mean-CPM discrepancy on the gated genes."""
    if bal:
        L = np.asarray(ctrl.sum(1)).ravel()
        order = np.argsort(L, kind="stable")
        step = len(order) / n
        off = rng.random()
        idx = order[np.clip(((np.arange(n) + off) * step).astype(int), 0, len(order) - 1)]
        return idx
    if not strat:
        return rng.choice(ctrl.shape[0], size=n, replace=False)
    c = pd.read_parquet(H.R / "reports/sc_hybrid_2026-09-25/phase" / f"{panel}.cells.parquet")
    ph = c[c.target == "non-targeting"].phase.values; assert len(ph) == ctrl.shape[0]
    pc = np.array([(ph == p).mean() for p in PH]); w = pc * n; k = np.floor(w).astype(int)
    k[np.argsort(-(w - k))[:n - k.sum()]] += 1
    return np.concatenate([rng.permutation(np.flatnonzero(ph == p))[:kk] for p, kk in zip(PH, k)])


def cs_expected(ctrl, Ld, low, chunk=2000):
    """Mean over ALL controls of the cellshift count ((x/L*1e4+1)e^d - 1) L/1e4 clipped at 0 (low genes: x e^d)."""
    ctrl = ctrl.tocsr(); N = ctrl.shape[0]; out = np.zeros(Ld.shape); E = np.exp(Ld)
    for a in range(0, N, chunk):
        X = np.asarray(ctrl[a:a + chunk].todense(), dtype=np.float32); L = X.sum(1, keepdims=True)
        base = X / L * 1e4 + 1.0
        for t in range(Ld.shape[0]):
            e = E[t].astype(np.float32)[None, :]
            w = (base * e - 1.0) * L / 1e4; np.clip(w, 0, None, out=w)
            w[:, low] = X[:, low] * e[:, low]
            out[t] += w.sum(0)
    return out / N


def cs_cpm_expected(ctrl, Ld, low, chunk=2000, pc=1.0):
    """Mean over ALL controls of the per-cell CPM (fraction) of the cellshift emission, per target (n_t x G)."""
    ctrl = ctrl.tocsr(); N = ctrl.shape[0]; out = np.zeros(Ld.shape); E = np.exp(Ld)
    for a in range(0, N, chunk):
        X = np.asarray(ctrl[a:a + chunk].todense(), dtype=np.float32); L = X.sum(1, keepdims=True); base = X / L * 1e4 + pc
        for t in range(Ld.shape[0]):
            e = E[t].astype(np.float32)[None, :]
            w = (base * e - pc) * L / 1e4; np.clip(w, 0, None, out=w); w[:, low] = X[:, low] * e[:, low]
            out[t] += (w / w.sum(1, keepdims=True)).sum(0)
    return out / N


def kink_gamma(ctrl, Ld, Od, n, factor=1.17, draws=4, seed=0, margin=1.0):
    """gamma at which the prediction's across-perturbation spread B_pred (bulk_lognorm means of the pinned sums, TS 5e4)
    equals the capped correction claim ~ n_perts * factor * C_real(n); C_real(n) = delete-1 jackknife trace of n controls."""
    from cell_eval2.moments import jackknife_correction as JK
    from cell_eval2.prep import bulk_lognorm_means as BLM
    rng = np.random.default_rng(seed); ctrl = ctrl.tocsr()
    cr = float(np.mean([JK(ctrl[rng.choice(ctrl.shape[0], size=n, replace=False)], np.zeros(n, np.int64), 1, 5e4)[0] for _ in range(draws)]))
    claim = Ld.shape[0] * factor * margin * cr; fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; j = np.isfinite(Od)
    def B(g):
        e = np.exp(g * Ld); e[j] = np.exp(Od[j]); mu = BLM(fm[None, :] * e * n, 5e4); return float(((mu - mu.mean(0)) ** 2).sum())
    lo, hi = 0.3, 6.0
    if B(hi) < claim:
        return hi, cr
    for _ in range(40):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(mid) < claim else (lo, mid)
    return 0.5 * (lo + hi), cr


def emit(ctrl, Ld, Od, n, o, gate, low, rs, panel, LdDE=None):
    rng = np.random.default_rng(o["seed"]); ctrl = ctrl.tocsr()
    pick = pick_block(ctrl, n, rng, o["strat"], panel, bal=o.get("lsb"))
    blk = np.asarray(ctrl[pick].todense(), dtype=np.float64); Lc = blk.sum(1)
    full_mean = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]
    base = blk / Lc[:, None] * 1e4 + 1.0
    colsum = blk.sum(0); share = np.where(colsum > 0, blk / np.maximum(colsum, 1e-12), 1.0 / n)
    K = o["K"]; lam_rank = np.linspace(1.0, o.get("le", 0.3), max(K, 1)); n_top = max(1, int(round(o.get("nt", 0.1) * n)))
    order_cells = np.argsort(-blk, axis=0, kind="stable")
    Tcs = cs_expected(ctrl, Ld, low) if o["target"] in ("cs", "mix") else None
    if o.get("cp"):                                # control per-cell CPM mean (the scorer's DE functional)
        Lall = np.asarray(ctrl.sum(1)).ravel(); cpm_mean = np.asarray((sparse.diags(1.0 / Lall) @ ctrl).mean(0)).ravel()
        logL = np.log(Lc)
        CPC = (cs_cpm_expected(ctrl, (Ld if LdDE is None else LdDE) * (o.get("gd", o["gamma"]) / o["gamma"]), low, pc=o.get("pc", 1.0))
               if o.get("cpt") == "cpc" else None)
    for t in range(Ld.shape[0]):
        j = np.isfinite(Od[t]); lt = Ld[t].copy()
        if o.get("cb", 1.0) != 1.0 and K:          # amplitude boost on the top-K called genes only (fh_dual boost, no dd)
            cand0 = np.flatnonzero(gate & ~j); lt[cand0[np.argsort(-np.abs(Ld[t, cand0]), kind="stable")[:K]]] *= o["cb"]
        e = np.exp(lt); e[j] = np.exp(Od[t, j])
        if o["target"] == "mult":
            T = np.round(full_mean * e * n)
        elif o["target"] == "mix":                 # geometric blend: log T = (1-b) log T_mult + b log T_cs
            b = o["beta"]; Tm = full_mean * e; Tc = np.maximum(Tcs[t], 0)
            T = np.where((Tm > 0) & (Tc > 0), Tm ** (1 - b) * Tc ** b, Tm); T[j] = Tm[j]; T = np.round(T * n)
        else:
            T = Tcs[t] * n; T[j] = full_mean[j] * e[j] * n; T = np.round(T)
        if o["shape"] == "share":
            W = share * T[None, :]
        else:
            w = (base * e[None, :] - 1.0) * Lc[:, None] / 1e4; np.clip(w, 0, None, out=w)
            m = low | j; w[:, m] = blk[:, m] * e[None, m]
            cw = w.sum(0); W = np.where(cw > 0, w / np.maximum(cw, 1e-12), 1.0 / n) * T[None, :]
        Lr = Ld if LdDE is None else LdDE                  # the DE side (top-K membership, up/down, cpc) follows LdDE
        if K:
            cand = np.flatnonzero(gate & ~j); top = cand[np.argsort(-np.abs(Lr[t, cand]), kind="stable")[:K]]
            lam = lam_rank[:len(top)]
            if rs is not None:
                lam = np.empty(len(top)); lam[np.argsort(-rs[t, top], kind="stable")] = lam_rank[:len(top)]
            if o.get("lcp") and len(top):                # weaken high-expression genes so the p-value order follows |d| (reach)
                cm = full_mean[top] / max(full_mean.sum(), 1e-12) * 1e6
                lam = lam * np.clip((np.median(cm) / np.maximum(cm, 1e-9)) ** o["lcp"], 0.2, 1.0)
            up = top[Lr[t, top] > 0]; lu = lam[Lr[t, top] > 0]; dn = top[Lr[t, top] < 0]; ld = lam[Lr[t, top] < 0]
            if up.size:
                if o["du"] > 0:
                    u = Lc ** (-o["du"]); u = u / u.sum()
                    W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (u[:, None] * T[up][None, :])
                else:
                    W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (T[up] / n)[None, :]
            if dn.size:
                conc = np.zeros((n, dn.size))
                np.put_along_axis(conc, order_cells[:n_top, dn], (T[dn] / n_top)[None, :].repeat(n_top, 0), axis=0)
                W[:, dn] = (1 - ld)[None, :] * W[:, dn] + ld[None, :] * conc
        if o.get("cp") and (K or o.get("cpall")):
            # depth tilt x_c ~ w_c L_c^s per gene so that sum_c x_c = T AND mean_c x_c/L_c = target (bisection on s)
            gs = np.flatnonzero(gate & ~j) if o.get("cpall") else top
            g = gs[T[gs] > 0]
            if g.size:
                W0g = W[:, g].copy()                             # pre-tilt allocation (same column sums)
                w = W[:, g] / np.maximum(W[:, g].sum(0), 1e-12)
                tgt = CPC[t, g] if CPC is not None else cpm_mean[g] * e[g]
                want = n * tgt / T[g]                           # target sum_c (x_c/T)/L_c
                Lcur = Lc
                for it in range(1 + o.get("ci", 0)):             # ci > 0: re-solve against the emitted cells' own libraries
                    lo = np.full(g.size, -4.0); hi = np.full(g.size, 4.0)
                    for _ in range(30):
                        mid = 0.5 * (lo + hi); v = w * np.exp(mid[None, :] * logL[:, None]); v /= v.sum(0, keepdims=True)
                        f = (v / Lcur[:, None]).sum(0)           # decreasing in s
                        hi = np.where(f < want, mid, hi); lo = np.where(f < want, lo, mid)
                    s_ = 0.5 * (lo + hi); v = w * np.exp(s_[None, :] * logL[:, None]); v /= v.sum(0, keepdims=True)
                    W[:, g] = v * T[g][None, :]
                    Lcur = np.maximum(W.sum(1), 1.0)
                lim = o.get("maxcell", 8e5); Wt = W[:, g].copy(); a_ = 1.0
                while W.sum(1).max() > lim and a_ > 1e-3:        # scorer limit 1e6/cell: relax the tilt, sums stay exact
                    a_ *= 0.5; W[:, g] = W0g + a_ * (Wt - W0g)
        fl = np.floor(W); rem = (T - fl.sum(0)).astype(int); frac = W - fl; cols = np.flatnonzero(rem > 0)
        if cols.size:
            oo = np.argsort(-frac[:, cols], axis=0); mask = np.arange(n)[:, None] < rem[cols][None, :]
            add = np.zeros((n, cols.size)); np.put_along_axis(add, oo, mask.astype(float), axis=0); fl[:, cols] += add
        b = sparse.csr_matrix(fl.astype(np.float32)); b.eliminate_zeros(); yield b


def score(arm, panel):
    out = OUT / arm / panel
    if (out / "result.json").exists():
        return json.loads((out / "result.json").read_text())
    out.mkdir(parents=True, exist_ok=True); t0 = time.time(); o = parse(arm); info = {}
    global NB_TOPK_ENS
    NB_TOPK_ENS = o.get("nbe") or []          # process-local; empty => the shipped single cut, bit-identical
    if o.get("delta") == "unseen_knn":             # panel targets treated as unmeasured: neighbour estimate only
        def unseen(NBD):
            keep = H.NB
            try:
                H.NB = NBD; return H.arm_delta("unseen_knn", panel)[:2]
            finally:
                H.NB = keep
        nbk = o.get("nb", "D")
        if nbk == "D":
            d, ov = unseen(H.NB)
        elif nbk == "S":
            d, ov = unseen(string_nb())
        elif nbk.startswith("C"):
            d, ov = unseen(json.loads((OUT.parent / f"concat_neighbours_a{float(nbk[1:]):g}.json").read_text()))
        else:
            d1, ov = unseen(H.NB); d2, _ = unseen(string_nb()); d = 0.5 * (d1 + d2)
    elif "s4" in o:
        P0_, C0_ = H.PARTS, H.CACHE
        try:
            if o.get("dk0"):                           # S10k/S11 sum-side sources: KOLF zeroed in the fusion AND the library
                H.PARTS = [(ww, dict(c, weights={**c["weights"], "KOLF": 0.0}), a) for ww, c, a in P0_]
                H.CACHE = C0_ / "dk0"; H.CACHE.mkdir(parents=True, exist_ok=True)   # delta_s4 cache keys ignore weights
            if o.get("xw"): H.GX = {"xw": o["xw"]}
            d, ov = delta_s4(panel, o["s4"], o.get("bw", 0.3), o.get("s4w", 0.5),
                             mc_list=MC7 if o.get("ens") == "P1" else None, c_list=o.get("cens"), wmul=o.get("wmul"),
                             fwins=float(o.get("fwins", 0.0)), crob=float(o.get("crob", 0.0)), wadd=o.get("wadd"))
            H.GX = {}
            if o.get("bmb") and o.get("ens") == "P1":      # averaging shrinks the across-target spread, so the ensemble's OWN
                d_bmb, _ = delta_s4(panel, o["s4"], o.get("bw", 0.3), o.get("s4w", 0.5))   # budget is smaller than the base's
                o["_bmb_delta"] = d_bmb                    # (kink_base 0.901 -> 0.939 on H1) and matching to it costs the MSE
        finally:
            H.PARTS, H.CACHE = P0_, C0_; H.GX = {}
        if o.get("deo") or o.get("bro") or o.get("brs"):   # the SHIPPED main recipe (official S11: DE side + budget reference +
            o["_old_delta"] = delta_s4(panel, o.get("deo") or o.get("bro") or o.get("brs"), o.get("bw", 0.3), o.get("s4w", 0.5))[0]   # agree file)

    elif "ens" in o:
        d, ov = delta_ens(panel, o["ens"])
    elif "fs" in o or "bw" in o or "nb" in o or "wmul" in o or "cfgx" in o or "cfgx_lib" in o:
        d, ov = delta_fs(panel, o.get("fs", 0.5), o.get("bw", 0.3), o.get("pt", False), o.get("nb", "D"), o.get("wmul"), o.get("cfgx"),
                         o.get("cfgx_lib"))
    else:
        d, ov, _ = H.arm_delta(o.get("delta", "blend30_e"), panel)
    if o.get("dsrc"):                             # DE side = the SAME recipe built on a variant source library (A-1 denoising):
        vsrc_ = H.R / "reports/sc_hybrid_2026-09-25/denoise_src" / o["dsrc"]   # sums, budget reference, overrides stay original;
        vcache_ = H.R / "reports/sc_hybrid_2026-09-25/denoise_src" / f"cache_{o['dsrc']}_de"; vcache_.mkdir(parents=True, exist_ok=True)
        P0_, C0_, S0_ = H.PARTS, H.CACHE, H.TM.SRC   # the DE delta and its ranking score follow the variant (via _de_pre)
        try:
            H.TM.SRC = vsrc_; H.TM._cache.clear(); H.CACHE = vcache_
            if o.get("dk0"):
                H.PARTS = [(ww, dict(c, weights={**c["weights"], "KOLF": 0.0}), a) for ww, c, a in P0_]
                H.CACHE = vcache_ / "dk0"; H.CACHE.mkdir(parents=True, exist_ok=True)
            if o.get("xw"): H.GX = {"xw": o["xw"]}
            base0_ = H.base
            if o.get("dslib"):
                lsrc_ = H.R / "reports/sc_hybrid_2026-09-25/denoise_src" / o["dslib"]
                def base_lib_(targets, excl, key, **kw):
                    if targets is H.LIB:
                        S1_, C1_ = H.TM.SRC, H.CACHE
                        try:
                            H.TM.SRC = lsrc_; H.TM._cache.clear(); H.CACHE = C1_ / f"lib_{o['dslib']}"; H.CACHE.mkdir(parents=True, exist_ok=True)
                            return base0_(targets, excl, key, **kw)
                        finally:
                            H.TM.SRC, H.CACHE = S1_, C1_; H.TM._cache.clear()
                    return base0_(targets, excl, key, **kw)
                H.base = base_lib_
            d_de_v, _ = delta_s4(panel, o["s4"], o.get("bw", 0.3), o.get("s4w", 0.5), mc_list=MC7 if o.get("ens") == "P1" else None, c_list=o.get("cens"), wmul=o.get("wmul"),
                             fwins=float(o.get("fwins", 0.0)), crob=float(o.get("crob", 0.0)), wadd=o.get("wadd"))
        finally:
            H.PARTS, H.CACHE, H.TM.SRC = P0_, C0_, S0_; H.TM._cache.clear(); H.GX = {}; H.base = base0_
        info["dsrc"] = dict(tag=o["dsrc"], cos_to_sum_delta=round(float((d_de_v * d).sum() / (np.linalg.norm(d_de_v) * np.linalg.norm(d) + 1e-12)), 5))
        print("dsrc", info["dsrc"], flush=True)
        assert info["dsrc"]["cos_to_sum_delta"] < 0.99999, "variant DE delta equals the base delta (stale cache?)"
        o["_de_pre"] = d_de_v
    if o.get("sl"):                               # SAME-LINE replacement, exactly the builder's same_line block (tv2_submit_build):
        src_, w_sl = o["sl"][0], float(min(max(o["sl"][1], 0.0), 1.0))   # covered targets (>= 20 cells) <- (1-w) d + w S |d|/|S|
        zs_ = H.TM.load(src_); cen_ = np.nan_to_num(zs_["center_c50"]); d_pre_ = d; d = d.copy(); ncov_ = 0   # on the genes the source
        for k_, t_ in enumerate(H.PANELS[panel]["targets"]):   # measures; override / own gene excluded from both; applied to the
            i_ = zs_["pos"].get(H.TM.canon(str(t_)))           # recipe delta, i.e. sums AND DE side; budget refs (_bmb/_old) untouched
            if i_ is None or zs_["ncell"][i_] < 20:
                continue
            v_ = np.log((zs_["fp"][i_].astype(np.float64) * 1e6 + 50) / (zs_["fc"][i_].astype(np.float64) * 1e6 + 50)) - cen_
            m_ = np.isfinite(v_) & ~np.isfinite(ov[k_]); jt_ = H.GIC.get(H.TM.canon(str(t_)))
            if jt_ is not None:
                m_[jt_] = False
            v_ = np.where(m_, v_, 0.0); nv_ = np.linalg.norm(v_)
            if nv_ <= 0:
                continue
            na_ = np.linalg.norm(np.where(m_, d[k_], 0.0))
            d[k_] = np.where(m_, (1 - w_sl) * d[k_] + w_sl * v_ * (na_ / nv_), d[k_]); ncov_ += 1
        info["same_line"] = dict(source=src_, w=w_sl, covered=ncov_); print("same-line", info["same_line"], flush=True)
        if o.get("sl_sum"):                       # DE side keeps the pre-replacement (untransformed) recipe delta
            o["_de_pre"] = d_pre_
    r = H.real(panel); P = H.PANELS[panel]; ok = r["gmap"] >= 0
    Ld = H.on_panel(o["gamma"] * d, panel); Od = np.full(Ld.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
    gate = r["cpm"] >= o.get("gt", 5.0); low = ~gate
    rs = FD.agreement_score(panel, (o["_de_pre"] if o.get("_de_pre") is not None else d)) if o["ag"] else None   # sls: DE-side delta
    if o["ag"] and o.get("agc"):                  # builder-exact ranking (tv2_agree_prep: agreement with the family-mean sign,
        d_de_ = o["_de_pre"] if o.get("_de_pre") is not None else d   # frozen, NOT with the candidate's sign; own family out)
        tg_ = H.PANELS[panel]["targets"]; g_ = r["gmap"][ok]; fl_ = H.fam_lfc(tg_, H.EXCL[panel])
        S_ = np.stack([np.nan_to_num(v[:, g_]) for v in fl_.values()]); M_ = np.stack([np.isfinite(v[:, g_]) for v in fl_.values()])
        mean_ = np.where(M_.sum(0) > 0, S_.sum(0) / np.maximum(M_.sum(0), 1), 0.0)
        fr_ = np.where(M_.sum(0) > 0, (M_ & (np.sign(S_) == np.sign(mean_)[None])).sum(0) / np.maximum(M_.sum(0), 1), 0.5)
        rs = np.zeros((len(tg_), len(r["genes"]))); rs[:, ok] = np.abs(d_de_[:, g_]) * np.clip(2 * fr_ - 1, 0, 1)
        del S_, M_; info["agree"] = "consensus (builder)"
    if o["target"] == "csraw":                    # the harness default emission (blend30_e as scored by fh_harness)
        from tv2_emit import emit as _cs
        X = _cs(r["ctrl"], Ld, P["n_cells"], o["seed"], "shared_cellshift_gated", override=Od)
    else:
        LdDE = None
    if o.get("tscale"):                            # oracle per-target least-squares slope of the real fold change on Ld
        Rr = r["R"]; msk = gate[None, :] & ~np.isfinite(Od)          # own gene / promoter neighbours (override) excluded
        num = np.where(msk, Ld * Rr, 0).sum(1); den = np.where(msk, Ld ** 2, 0).sum(1)
        st = np.clip(np.where(den > 0, num / np.maximum(den, 1e-12), 1.0), 0.0, 3.0)
        info["slope_median"] = float(np.median(st)); info["slope_p10_p90"] = [float(np.quantile(st, .1)), float(np.quantile(st, .9))]
        if o["tscale"] in ("so", "dso"):
            LdDE = Ld.copy() if o["tscale"] == "so" else None
            Ld = Ld * st[:, None]
            if o["tscale"] == "dso":
                LdDE = None                                         # DE side follows the scaled Ld as well
        else:                                                      # "do": DE side only
            LdDE = Ld * st[:, None]
    Lb0_ = (Ld.copy() if LdDE is None else LdDE.copy())   # untransformed recipe (budget-matching base; before any sum/DE transform)
    if o.get("_bmb_delta") is not None:
        Lb0_ = H.on_panel(o["gamma"] * o.pop("_bmb_delta"), panel); info["bmb"] = True   # pop: `o` is serialised
    if o.get("_old_delta") is not None:            # official S11 mirror: the main (DE) recipe is also the budget reference
        d_old_ = o.pop("_old_delta"); Lold_ = H.on_panel(o["gamma"] * d_old_, panel)
        Lb0_ = Lold_.copy(); info["budget_ref"] = "old"
        if not o.get("brs"):                       # brs: the ranking score follows the NEW DE delta, as the builder would
            rs = FD.agreement_score(panel, d_old_) if o["ag"] else None
        if o.get("deo"):
            LdDE = Lold_.copy(); info["de_side"] = "old"
    if o.get("sbw"):                              # SUM side from a different neighbour-blend weight; DE side keeps the default recipe
        dsb, _ = delta_fs(panel, 0.5, o["sbw"]); LdDE = Ld.copy() if LdDE is None else LdDE; Ld = H.on_panel(dsb, panel)
    if o.get("sw"):                               # SUM-side reliability weighting from source-family agreement (DE side unchanged)
        fl = H.fam_lfc(r["targets"], H.EXCL[panel]); gm = r["gmap"]; F = []
        for f_, Lf in fl.items():
            M_ = np.full(Ld.shape, np.nan); M_[:, ok] = Lf[:, gm[ok]]; F.append(M_)
        F = np.stack(F); meas = np.isfinite(F); nme = meas.sum(0)
        if o["sw"].startswith("snr"):
            mu_ = np.nanmean(np.where(meas, F, np.nan), 0); sd_ = np.nanstd(np.where(meas, F, np.nan), 0)
            wgt = np.where(nme >= 2, mu_ ** 2 / (mu_ ** 2 + np.nan_to_num(sd_) ** 2 / np.maximum(nme, 1) + 1e-6), 0.5)
        else:
            ag_ = np.where(meas, np.sign(np.nan_to_num(F)) == np.sign(Ld)[None], False).sum(0) / np.maximum(nme, 1)
            wgt = np.where(nme >= 2, np.clip(2 * ag_ - 1, 0, 1), 0.5)
        if o["sw"].endswith("05"):
            wgt = np.sqrt(wgt)
        wgt = np.nan_to_num(wgt, nan=0.5); info["sw_mean_weight"] = round(float(wgt.mean()), 3)
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld * wgt
    if o.get("ew"):                               # SUM-side weight (1 + cpm/c0)^-kappa from the DESTINATION controls' CPM
        c0_, k_ = o["ew"]; w_ = (1 + r["cpm"] / c0_) ** (-k_); w_ = w_ / np.average(w_, weights=np.maximum(r["cpm"], 1e-3))
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld * w_[None, :]
    if o.get("ct"):                               # SUM-side partial centring (override entries excluded from the mean)
        mt = np.nanmean(np.where(np.isfinite(Od), np.nan, Ld), 0)
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld - o["ct"] * np.nan_to_num(mt)[None, :]
        info["ct_shared_share"] = round(float(len(Ld) * (np.nan_to_num(mt) ** 2).sum() / max((LdDE ** 2).sum(), 1e-12)), 4)
    if o.get("pcr"):                              # SUM-side partial removal of the top principal axes of the source library responses
        from sklearn.utils.extmath import randomized_svd
        k_, l_ = o["pcr"]; dlib = H.lib_delta(panel, H.EXCL[panel]); Dl = np.zeros((dlib.shape[0], len(r["genes"]))); Dl[:, ok] = dlib[:, r["gmap"][ok]]
        Dl = Dl[np.abs(Dl).sum(1) > 0]; _, _, Vt_ = randomized_svd(Dl - Dl.mean(0), n_components=k_, random_state=0)
        P_ = Ld @ Vt_.T
        if o.get("pcr_guard"):
            lo_, hi_ = o["pcr_guard"]; f_ = float((P_ ** 2).sum() / max(float((Ld ** 2).sum()), 1e-12))
            l_ = l_ * float(np.clip((hi_ - f_) / (hi_ - lo_), 0.0, 1.0)); info.update(pcr_energy=round(f_, 4), pcr_lambda=round(l_, 3))
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld - l_ * P_ @ Vt_
    if o.get("dis"):                              # SUM-side HUBNESS reduction: remove part of each target's projection on
        k4_, l4_ = o["dis"]                       # the centroid of its own k most similar OTHER predictions (arXiv 2511.16954:
        U_ = Ld / np.maximum(np.linalg.norm(Ld, axis=1, keepdims=True), 1e-12)   # PDS is a nearest-neighbour rank and its
        G_ = U_ @ U_.T; np.fill_diagonal(G_, -np.inf)                            # hard part is hubness, which a GLOBAL
        idx_ = np.argsort(-G_, axis=1)[:, :k4_]                                  # operator like `pcr` cannot address)
        C_ = np.stack([Ld[idx_[i]].mean(0) for i in range(len(Ld))])
        Cn_ = C_ / np.maximum(np.linalg.norm(C_, axis=1, keepdims=True), 1e-12)
        LdDE = Ld.copy() if LdDE is None else LdDE
        Ld = Ld - l4_ * (np.einsum("ij,ij->i", Ld, Cn_))[:, None] * Cn_
        info["dis"] = [k4_, l4_]
    if o.get("dpcr"):                             # DE-side partial removal of the top source-library axes (sums unchanged)
        from sklearn.utils.extmath import randomized_svd
        k3_, l3_ = o["dpcr"]; dlib = H.lib_delta(panel, H.EXCL[panel]); Dl = np.zeros((dlib.shape[0], len(r["genes"]))); Dl[:, ok] = dlib[:, r["gmap"][ok]]
        Dl = Dl[np.abs(Dl).sum(1) > 0]; _, _, Vd_ = randomized_svd(Dl - Dl.mean(0), n_components=k3_, random_state=0)
        base_de = Ld.copy() if LdDE is None else LdDE; P3_ = base_de @ Vd_.T
        if o.get("dpcr_guard"):
            lo_, hi_ = o["dpcr_guard"]; f3_ = float((P3_ ** 2).sum() / max(float((base_de ** 2).sum()), 1e-12))
            l3_ = l3_ * float(np.clip((hi_ - f3_) / (hi_ - lo_), 0.0, 1.0)); info.update(dpcr_energy=round(f3_, 4), dpcr_lambda=round(l3_, 3))
        LdDE = base_de - l3_ * P3_ @ Vd_
    if o.get("dsw"):                              # DE-side source-family SNR weighting (sums unchanged)
        fl2 = H.fam_lfc(r["targets"], H.EXCL[panel]); F2 = []
        for f_, Lf in fl2.items():
            M_ = np.full(Ld.shape, np.nan); M_[:, ok] = Lf[:, r["gmap"][ok]]; F2.append(M_)
        F2 = np.stack(F2); me2 = np.isfinite(F2); nm2 = me2.sum(0)
        mu2 = np.nanmean(np.where(me2, F2, np.nan), 0); sd2 = np.nanstd(np.where(me2, F2, np.nan), 0)
        snr2 = np.nan_to_num(np.where(nm2 >= 2, mu2 ** 2 / (mu2 ** 2 + np.nan_to_num(sd2) ** 2 / np.maximum(nm2, 1) + 1e-6), 0.5), nan=0.5)
        LdDE = (Ld.copy() if LdDE is None else LdDE) * snr2 ** o["dsw"]
    if o.get("dct"):                              # DE-side partial centring over targets (override entries excluded)
        base_de = Ld.copy() if LdDE is None else LdDE
        mt_ = np.nanmean(np.where(np.isfinite(Od), np.nan, base_de), 0)
        LdDE = base_de - o["dct"] * np.nan_to_num(mt_)[None, :]
    if o.get("_de_pre") is not None:              # sls: the DE side = the untransformed PRE-replacement delta (sums keep it)
        LdDE = H.on_panel(o["gamma"] * o.pop("_de_pre"), panel); info["de_side"] = "pre_same_line"
    if o.get("dpj") and LdDE is not None:          # A-3 (2026-09-30): DE delta <- P_K(DE) + rho (DE - P_K(DE)), P_K = projection onto
        from sklearn.utils.extmath import randomized_svd   # the top-K right singular vectors of the library recipe deltas
        K_, rho_ = o["dpj"]; dlib_ = H.lib_delta(panel, H.EXCL[panel])   # (own family excluded); override genes untouched;
        Dl_ = np.zeros((dlib_.shape[0], len(r["genes"]))); Dl_[:, ok] = dlib_[:, r["gmap"][ok]]; Dl_ = Dl_[np.abs(Dl_).sum(1) > 0]
        _, _, Vk_ = randomized_svd(Dl_, n_components=int(K_), n_iter=5, random_state=0); del Dl_, dlib_
        msk_ = np.isfinite(Od); X_ = np.where(msk_, 0.0, LdDE); Pj_ = (X_ @ Vk_.T) @ Vk_
        LdDE = np.where(msk_, LdDE, Pj_ + rho_ * (X_ - Pj_))
        if o["ag"]:                                 # the ranking follows the projected DE delta (as |Dde| x agree in the builder)
            dtmp_ = np.zeros((LdDE.shape[0], H.NG)); dtmp_[:, r["gmap"][ok]] = np.nan_to_num(LdDE[:, ok]) / o["gamma"]
            rs = FD.agreement_score(panel, dtmp_)
        info["dpj"] = dict(K=int(K_), rho=rho_, energy_kept=round(float((Pj_ ** 2).sum() / max((X_ ** 2).sum(), 1e-12)), 4))
        print("dpj", info["dpj"], flush=True)
    if o.get("ta"):                               # per-target SUM amplitude (PDS-invariant): a_k = rank(cross-family agreement_k)^beta,
        z_ = np.load(H.R / "reports/sc_hybrid_2026-09-25/h1fuse" / f"{panel}_data.npz")   # family means at min_cells 20 with the
        assert [str(t) for t in z_["tg"]] == [str(t) for t in H.PANELS[panel]["targets"]], "ta: target order"   # panel's exclusions
        xa_, ma_ = z_["x"][0], z_["m"][0]; ag_ = np.full(xa_.shape[0], np.nan)
        for k_ in range(xa_.shape[0]):                # mean pairwise cosine among families measuring >= 500 genes of target k
            cols_ = [f_ for f_ in range(xa_.shape[-1]) if ma_[k_, :, f_].sum() > 500]; cs_ = []
            for i_ in range(len(cols_)):
                for j_ in range(i_ + 1, len(cols_)):
                    a_, b_ = xa_[k_, :, cols_[i_]], xa_[k_, :, cols_[j_]]; mm_ = ma_[k_, :, cols_[i_]] & ma_[k_, :, cols_[j_]]
                    cs_.append(float(a_[mm_] @ b_[mm_] / (np.linalg.norm(a_[mm_]) * np.linalg.norm(b_[mm_]) + 1e-12)))
            if cs_: ag_[k_] = float(np.mean(cs_))
        ag_ = np.nan_to_num(ag_, nan=float(np.nanmedian(ag_))); rk_ = (np.argsort(np.argsort(ag_, kind="stable"), kind="stable") + 1) / len(ag_)
        if o.get("tlrn"):                         # leave-one-panel-out GBM prediction of the per-target true cosine
            pr_ = np.load(H.R / "reports/sc_hybrid_2026-09-25/h1fuse" / f"{panel}_amppred.npz")["pred"]
            rk_ = (np.argsort(np.argsort(pr_, kind="stable"), kind="stable") + 1) / len(pr_)
        if o.get("tcomp"):                        # composite: + #families measuring >= 500 genes, + |recipe delta| (official axis)
            rnk_ = lambda v: (np.argsort(np.argsort(v, kind="stable"), kind="stable") + 1) / len(v)
            nf_ = (ma_.sum(1) > 500).sum(-1).astype(float); nr_ = np.linalg.norm(d, axis=1)
            rk_ = rnk_(rnk_(ag_) + rnk_(nf_) + rnk_(nr_))
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld * (rk_ ** o["ta"])[:, None]   # sums only; budget match below restores B_pred
        info["ta"] = dict(beta=o["ta"], agree_median=round(float(np.median(ag_)), 4), a_median=round(float(np.median(rk_ ** o["ta"])), 3))
        if o.get("tad"):                          # DE-side amplitude (top-K membership order within a target is unchanged; the cpc
            LdDE = LdDE * (rk_ ** o["tad"])[:, None]; info["tad"] = o["tad"]   # composition target and up/down LFCs shrink)
    if o.get("tzb"):                              # rho-hat MAGNITUDE reshaping of the per-target sums, placed BEFORE the
        gb6_, db6_ = o["tzb"]                     # budget match so `bm` restores B_pred (budget-INVARIANT: offline price valid)
        pb_ = np.load(H.R / "reports/sc_hybrid_2026-09-25/h1fuse" / f"{panel}_amppred.npz")["pred"]
        wb_ = np.maximum(pb_ - pb_.min(), 1e-6) ** gb6_; wb_ = wb_ / np.median(wb_)
        LdDE = Ld.copy() if LdDE is None else LdDE
        Ld = Ld * wb_[:, None]
        if db6_:                                  # the DE side follows the shrink at w**a, so the cpc composition pin stays
            LdDE = LdDE * (wb_ ** db6_)[:, None]  # reachable for the shrunk targets (OR_X NMAE failure mode, FINDINGS 25.6)
        info["tzb"] = [gb6_, db6_]; print("tzb", info["tzb"], flush=True)
    if o.get("gk") and o.get("bm") and LdDE is not None:   # budget matching: transformed sums keep the base recipe's B_pred
        import tv2_emit as EMO
        gkb_, cr_ = kink_gamma(r["ctrl"], Lb0_, Od, P["n_cells"], factor=o["gk"], margin=o.get("mg", 1.0)); gb_ = max(1.0, gkb_)
        Bt_ = EMO.budget_B(r["ctrl"], Lb0_, Od, P["n_cells"])(gb_); gs_ = EMO.budget_match(r["ctrl"], Ld, Od, P["n_cells"], Bt_)
        info.update(kink_base=round(gkb_, 3), gamma_base=round(gb_, 3), gamma_sum=round(gs_, 3), c_real_ctrl=round(cr_, 2), budget_match=True)
        Ld = Ld * gs_
    elif o.get("gk"):                             # test-time kink calibration of the SUM amplitude (DE side unchanged)
        gk_, cr_ = kink_gamma(r["ctrl"], Ld, Od, P["n_cells"], factor=o["gk"], margin=o.get("mg", 1.0))
        gs_ = max(1.0, gk_); info.update(kink=round(gk_, 3), gamma_sum=round(gs_, 3), c_real_ctrl=round(cr_, 2))
        LdDE = Ld.copy() if LdDE is None else LdDE; Ld = Ld * gs_
        if o.get("gdf"):
            LdDE = LdDE * gs_ ** o["gdf"]; info["gamma_de"] = round(gs_ ** o["gdf"], 3)
    if o.get("tz"):                               # SECOND per-target amplitude, applied AFTER the budget match: uses the
        gz_, sz_ = o["tz"]                        # rho-hat MAGNITUDE (not its rank, as `ta` does) and deliberately lets
        pz_ = np.load(H.R / "reports/sc_hybrid_2026-09-25/h1fuse" / f"{panel}_amppred.npz")["pred"]   # B_pred land near the
        wz_ = np.maximum(pz_ - pz_.min(), 1e-6) ** gz_; wz_ = wz_ / np.median(wz_)                    # rebate cap (B/claim ~1)
        LdDE = Ld.copy() if LdDE is None else LdDE                                   # DE side keeps the pre-tz sums
        Ld = Ld * (wz_ * sz_)[:, None]
        info["tz"] = [gz_, sz_, round(float(np.median(wz_ * sz_)), 3)]; print("tz", info["tz"], flush=True)
    if o.get("drp"):                              # CONTINUOUS version of drc: within each target multiply the DE lfc by
        g_ = o["drp"]                             # (descending |lfc| rank / N)^gamma, normalised to median 1 -> a smooth
        LdDE = Ld.copy() if LdDE is None else LdDE.copy()   # monotone magnitude calibration; the BULK target is untouched,
        for t_ in range(LdDE.shape[0]):           # so PDS / MSE stay invariant exactly as for drc
            idx_ = np.flatnonzero(~np.isfinite(Od[t_]))
            o_ = idx_[np.argsort(-np.abs(LdDE[t_, idx_]), kind="stable")]
            w_ = ((np.arange(len(o_)) + 1) / len(o_)) ** g_
            LdDE[t_, o_] *= w_ / np.median(w_)
        info["drp"] = g_; print("drp", g_, flush=True)
    if o.get("drc"):                              # public description of leaderboard entry v17a (team William, 0.2423): within each
        a_, b_ = o["drc"]                         # target the top 5 % of genes by |DE lfc| x a and the bottom 70 % x b, the BULK
        LdDE = Ld.copy() if LdDE is None else LdDE.copy()   # (sum) target unchanged -> PDS / MSE invariant, NMAE / reach move
        for t_ in range(LdDE.shape[0]):
            idx_ = np.flatnonzero(~np.isfinite(Od[t_]))
            o_ = idx_[np.argsort(-np.abs(LdDE[t_, idx_]), kind="stable")]
            nt_ = int(round(0.05 * len(o_))); nb_ = int(round(0.70 * len(o_)))
            LdDE[t_, o_[:nt_]] *= a_
            if nb_: LdDE[t_, o_[-nb_:]] *= b_
        info["drc"] = [a_, b_]; print("drc", info["drc"], flush=True)
    X = sparse.vstack(list(emit(r["ctrl"], Ld, Od, P["n_cells"], o, gate, low, rs, panel, LdDE)), format="csr")
    tg = r["targets"]
    pred = ad.AnnData(X=sparse.vstack([X, r["ctrl"]], format="csr"),
                      obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(tg, P["n_cells"]), *["non-targeting"] * r["ctrl"].shape[0]])},
                                       index=[f"p{i}" for i in range(X.shape[0] + r["ctrl"].shape[0])]),
                      var=pd.DataFrame(index=r["genes"]))
    pp = out / "pred.h5ad"; pred.write_h5ad(pp); del pred, X
    cache = (H.R / "reports/transfer_v2_2026-09-23/cache_real_H1") if panel == "H1" else (H.D / f"cache_real_{panel}")
    cmd = [str(H.S.CE), "run", "-ap", str(pp), "-ar", P["path"], "--preset", "vcc2026", "--pert-col", "target",
           "--control", "non-targeting", "--input-type", "counts", "-o", str(out / "cell_eval2"), "--cache-real", str(cache),
           "--set", "de.backend=pdex", "--set", "num_threads=2", "--set", "gather_threads=2"]
    if os.environ.get("SCH_KEEP"):
        cmd += ["--cache-pred", str(out / "cache_pred")]
    env = dict(os.environ, NUMBA_NUM_THREADS="2", OMP_NUM_THREADS="2")
    with (out / "scorer.log").open("w") as g:
        subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
    agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
    raw = {k: float(agg.loc["mean", k]) for k in H.S.CONST}; comp, scaled = H.S.composite(raw)
    res = dict(arm=arm, panel=panel, opts=o, info=info, raw=raw, scaled=scaled, composite=comp, seconds=round(time.time() - t0, 1))
    (out / "result.json").write_text(json.dumps(res, indent=1)); pp.unlink()
    if not os.environ.get("SCH_KEEP"): shutil.rmtree(out / "cell_eval2", ignore_errors=True)
    print(json.dumps(dict(arm=arm, panel=panel, composite=round(comp, 4), **{k[:12]: round(v, 3) for k, v in raw.items()},
                          secs=res["seconds"])), flush=True)
    return res


if __name__ == "__main__":
    a = _argv[1:]; only = [x[6:] for x in a if x.startswith("panel=")]
    for arm in [x for x in a if not x.startswith("panel=")]:
        for panel in (only or list(H.PANELS)):
            try:
                score(arm, panel)
            except Exception as ex:
                import traceback; traceback.print_exc(); print(json.dumps(dict(arm=arm, panel=panel, error=repr(ex))), flush=True)
    print("ALL_DONE", flush=True)
