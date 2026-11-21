"""Task 18 (plan `input/ideas/zero-shot-limits-context-complexity-plan.md`, RQ4 / H4): same-line data value vs modelling
gains on the H1 harness panel, real scorer (`cell-eval2 run --preset vcc2026 --set de.backend=pdex`), FULL composite.

Pre-registered H4: on the H1 panel with the real scorer, adding H1GW (w = 1) gains more composite than the largest single
modelling gain (M2, M3, or the best of the 17 architectures relative to its transfer baseline).
Rejected if composite(base + H1GW, w = 1) - composite(base) <= the largest single modelling gain.

Subcommands (vcc2026 env; nice -n 10; OMP/MKL/OpenBLAS 4 threads):
  build          DE-side deltas (harness mirror of F2's denoised DE side), the 7 arm packages
                 (results/rq4/pkg/<arm>/{pred.h5ad, inputs.npz, meta.json}), offline sum-side sanity (ERA pack evaluator)
  real ARM ...   score packages with the real scorer, one after the other (the queue wraps EACH call in
                 `flock results/.celleval.lock`; the queue starts only after unit aivc-paper-t17-real has exited)
  verdict        collect -> results/rq4_sameline.parquet + results/rq4_h4.json (provisional until the 4 verdict runs exist)

ARMS (H1 panel, 300 targets x 300 cells, 12,000 controls):
  w0      base = frozen recipe as Task 14 (harness M1 + F2 sum side incl. M2 and M3) + F2 DE side, no same-line data
  w0.25 / w0.5 / w0.75 / w1   base + H1GW (GSE295214) same-line replacement at weight w (sums AND DE side)
  noM2    base - M2 (no `--sum-tamp` amplitude; M3 kept, its weak half still from the agreement ranking)
  noM3    base - M3 (= Task 14 pre-M3 sums); DE side = base
ASSUMPTIONS (stated; checked in code; knowledge/synthesis/rq4_h4-assumptions.md)
  A1  Sums = Task 14 recipe verbatim (helpers copied from rq3_build_preds.py; parity asserted against the stored
      results/preds/H1_{m3,prem3}.npz lfc). Budget target Bt = the ERA pack Bt (Task 14 recomputed it, rel < 1e-5) and is
      NOT changed by the same-line replacement (builder: the budget reference D and the override OV are untouched).
  A2  DE side "as in the frozen F2 build" (tv2_submit_build dual_mult): Dde = the sum-side recipe built on the
      DENOISED source library (F2 `--de-denoise intended` = K562 sources program-denoised, rank 100, r 0.25; harness mirror
      `denoise_src/dnK100r0.25`, sch_emit token `dsrcdnK100r0.25`), the same recipe as the sum side's M1 (S4 Yorig part
      + 8-member min_cells fusion ensemble + 0.3 DepMap-neighbour blend from the denoised library, energy matched,
      own family H1 excluded), at gamma 1 (g0 * Dde; the budget scale acts on the sums only), then the same-line
      replacement (the builder applies it to Dover as well), then `--de-drc 0.8,1.1`. Emission = builder dual
      (K 3000, depth_up 0.5, agree ranking |lfc_de| x consensus agreement (builder formula, own family out = harness
      `agc`), cpc pin, ci 5, lsb depth-balanced block), through the harness emitter `sch_emit.emit`
      (arm `em_mult_share_k3000_du0.5_ag_cpc_ci5_lsb`, seed 0), i.e. the emitter behind every earlier harness
      composite (FINDINGS 20.20). F2 M1 items that Task 14 did not reproduce stay out on the DE side too: the `cens`
      pseudocount members, DLD1 family, promoter-off fusion, knn_noK (H1 is a stem line: knn with KOLF anyway),
      `sum.resid` (sums only in F2 anyway), S4noFK.
  A3  Same-line replacement = tv2_submit_build `same_line` verbatim on the official axis: for each panel target the
      source covers with >= 20 cells, on the genes the source measures, minus override entries and the target's own gene:
      d <- (1 - w) d + w v |d| / |v|, v = log((fp 1e6 + 50)/(fc 1e6 + 50)) - center_c50. Source H1GW =
      transfer_v2_2026-09-23/sources_new/H1GW.npz (sha256 417c861c...; frozen input manifest).
  A4  The M2 rank and the M3 weak half (rank(agreement)^0.75 <= median) are source-side and unchanged by w (the
      builder's agreement file is precomputed from the source library).
  A5  composite = tv2_score.composite (the vcc2026 preset scaling constants), one run per arm (emission seed 0).
      Seed-to-seed emission noise is NOT re-measured here; the reference spread of the earlier harness same-line arm
      (@s0/@s1/@s2) is reported for context only.
  A6  Best-architecture gain = arch_screen_2026-09-22/final_ranking.json `<line>_c` (official composite delta vs the
      refit_v2 transfer base, correct context, alpha 1, seed 0). PRIMARY (pre-specified before any Task 18 real result):
      the largest line-MEAN over K562 and HepG2 (the screen's ranking statistic); SENSITIVITY: the largest single-line
      value (max over 34 cells; winner's-curse-prone, the most favourable reading for modelling). Different panel
      (K562 / HepG2, not H1), different base (refit_v2, not the fusion pipeline), different emission; Task 10 re-checks.
  A7  Per-context value (Task 11 LOFO, converted to composite on H1) is not available yet -> null.
Seeds: emission rng 0 (sch_emit.emit), ctrl-PC draw rng 0 / randomized_svd random_state 0 (Task 14), kink jackknife rng 0."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import json, pathlib, shutil, subprocess, sys, time
import numpy as np

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); V = ROOT / "vcc2026_data"; RP = V / "reports"; SCR = V / "scripts"
OUTD = ROOT / "output/data-analysis"; RES = OUTD / "results"; PRED = RES / "preds"; W4 = RES / "rq4"; PKG = W4 / "pkg"
PACKS = RP / "era_2026-10-07/packs"; CACHE_REAL = RES / "rq3_cache/cache_real_H1"
CE = pathlib.Path("/home/wonmin-kim/miniconda3/envs/vcc2026/bin/cell-eval2")
P = "H1"; SL_SRC = "H1GW"; DSRC = "dnK100r0.25"; DRC = (0.8, 1.1)
EMIT_ARM = "em_mult_share_k3000_du0.5_ag_cpc_ci5_lsb"
WS = (0.0, 0.25, 0.5, 0.75, 1.0)
ARMS = [f"w{w:g}" for w in WS] + ["noM2", "noM3"]
VERDICT_ARMS = ("w0", "w1", "noM2", "noM3")
ARCH_JSON = RP / "arch_screen_2026-09-22/final_ranking.json"
REF_SL = "em_mult_share_k3000_du0.5_ag_cpc_ct0.5_pr3x0.75g_gk_mg0.8_bm_s4Yorig_ensP1_bmb_ci5_lsb_brsYorig_slH1GWx1"
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def jdump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)))


def imports():
    from cell_eval2.moments import jackknife_correction  # noqa: F401  (import order: BEFORE scripts/ joins sys.path)
    from cell_eval2.prep import bulk_lognorm_means  # noqa: F401
    sys.argv = ["x"]; sys.path.append(str(SCR)); sys.path.insert(0, str(OUTD))
    sys.path.insert(0, str(ROOT / "external/era-aivc/aivc"))
    import sch_emit as SE
    import fh_harness as H
    import common as CM
    from tasks.vcc_sum import evaluator as EV
    return SE, H, CM, EV


# ----------------------------------------------------------------------------- Task 14 sum-side helpers (verbatim)
def Bfun(fm, L, Od, n):
    from cell_eval2.prep import bulk_lognorm_means as BLM
    j = np.isfinite(Od)

    def B(g):
        e = np.exp(g * L); e[j] = np.exp(Od[j]); mu = BLM(fm[None, :] * e * n, 5e4); return float(((mu - mu.mean(0)) ** 2).sum())
    return B


def bmatch(fm, X, Od, n, Bt):
    B = Bfun(fm, X, Od, n); lo, hi = 0.05, 8.0
    for _ in range(40):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(mid) < Bt else (lo, mid)
    s = 0.5 * (lo + hi); assert 0.06 < s < 7.9, f"budget scale at the bisection edge: {s}"
    return s


def sum_side(L0, Od, V3, cpm, rk, beta=0.5, ew=(50.0, 0.25), Vd=None, wdf=None):
    ovj = np.isfinite(Od); X = L0.copy()
    if ew:
        wew = (1 + cpm / ew[0]) ** (-ew[1]); X = X * (wew / np.average(wew, weights=np.maximum(cpm, 1e-3)))[None, :]
    mt = np.nanmean(np.where(ovj, np.nan, X), 0); X = X - 0.5 * np.nan_to_num(mt)[None, :]
    P_ = X @ V3.T; f_ = float((P_ ** 2).sum() / (X ** 2).sum()); lam = 0.75 * float(np.clip((0.25 - f_) / 0.10, 0, 1))
    X = X - lam * P_ @ V3
    if rk is not None:
        X = X * (rk ** beta)[:, None]
    if Vd is not None:
        X = X - wdf[:, None] * ((X @ Vd.T) @ Vd)
    return X, dict(pcr_energy=round(f_, 4), pcr_lambda=round(lam, 4))


def cpm_of(ctrl, ok):
    tot = np.asarray(ctrl.sum(0)).ravel().astype(np.float64) * ok
    return tot / tot.sum() * 1e6


def on_panel(d, gmap):
    ok = gmap >= 0; out = np.zeros((d.shape[0], len(gmap))); out[:, ok] = d[:, gmap[ok]]; return out


# ----------------------------------------------------------------------------- builder same_line / drc (verbatim)
def same_line(arrs, targets, OV, w, TM):
    """tv2_submit_build same_line block; arrs (official axis) modified in place copies; returns (arrs, n_covered)."""
    zs = TM.load(SL_SRC); cen_ = np.nan_to_num(zs["center_c50"]); w_ = float(np.clip(w, 0.0, 1.0))
    arrs = [a.copy() for a in arrs]; ncov = 0
    if w_ == 0.0:
        return arrs, 0
    for k_, t_ in enumerate(targets):
        i_ = zs["pos"].get(TM.canon(t_))
        if i_ is None or zs["ncell"][i_] < 20:
            continue
        v_ = np.log((zs["fp"][i_].astype(np.float64) * 1e6 + 50) / (zs["fc"][i_].astype(np.float64) * 1e6 + 50)) - cen_
        m_ = np.isfinite(v_) & ~np.isfinite(OV[k_])
        jt_ = TM.GIC.get(TM.canon(t_))
        if jt_ is not None: m_[jt_] = False
        v_ = np.where(m_, v_, 0.0); nv_ = np.linalg.norm(v_)
        if nv_ <= 0:
            continue
        for arr_ in arrs:
            na_ = np.linalg.norm(np.where(m_, arr_[k_], 0.0))
            arr_[k_] = np.where(m_, (1 - w_) * arr_[k_] + w_ * v_ * (na_ / nv_), arr_[k_])
        ncov += 1
    return arrs, ncov


def drc(Dde, OV, a_, b_):
    Dde = Dde.copy()
    for t_ in range(Dde.shape[0]):
        idx_ = np.flatnonzero(~np.isfinite(OV[t_]))
        o_ = idx_[np.argsort(-np.abs(Dde[t_, idx_]), kind="stable")]
        nt_ = int(round(0.05 * len(o_))); nb_ = int(round(0.70 * len(o_)))
        Dde[t_, o_[:nt_]] *= a_
        if nb_: Dde[t_, o_[-nb_:]] *= b_
    return Dde


# ----------------------------------------------------------------------------- build
def cmd_build(only=None):
    import anndata as ad, pandas as pd
    from scipy import sparse
    SE, H, CM, EV = imports(); TM = H.TM
    W4.mkdir(parents=True, exist_ok=True); rep = dict(recipe_note=__doc__.split("ASSUMPTIONS")[0].strip())
    zm = np.load(PRED / "H1_m3.npz", allow_pickle=True); zp = np.load(PRED / "H1_prem3.npz", allow_pickle=True)
    sp = np.load(PACKS / "score_H1.npz"); zi = np.load(PACKS / "inputs_H1.npz", allow_pickle=False)
    tg = [str(t) for t in zi["targets"]]; assert tg == [str(t) for t in zm["targets"]] == [str(t) for t in H.PANELS[P]["targets"]]
    n = int(sp["n"]); gmap = sp["gmap"]; ok = gmap >= 0; Od = sp["Od"].astype(np.float64); V3 = sp["V3"].astype(np.float64)
    Bt = float(sp["Bt"]); NG = H.NG
    r = H.real(P); ctrl = r["ctrl"]; genes = np.asarray(r["genes"]).astype(str)   # ONE control load (memory; oomd 2026-11-05)
    r.pop("R", None)
    assert list(genes) == [str(g) for g in zm["genes"]] and np.array_equal(r["gmap"], gmap)
    assert str(CM.truth_path(P)).endswith(str(H.PANELS[P]["path"]).split("vcc2026_data/")[-1]), "truth path != Task 14"
    fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; assert np.allclose(fm, sp["fm"], rtol=1e-6, atol=1e-9)
    cpm = cpm_of(ctrl, ok); rk = zm["m2_rank"].astype(np.float64); wdf = np.where(zm["m3_weak"], 0.75, 0.0)
    Vd = zm["ctrl_V"].astype(np.float64)
    base_total = zi["base_total"].astype(np.float64); assert base_total.shape == (len(tg), NG)
    L0 = on_panel(base_total, gmap)

    # --- A1 parity: Task 14 sums
    Xm, tinfo = sum_side(L0, Od, V3, cpm, rk, Vd=Vd, wdf=wdf); Xp, _ = sum_side(L0, Od, V3, cpm, rk)
    Lm = Xm * bmatch(fm, Xm, Od, n, Bt); Lp = Xp * bmatch(fm, Xp, Od, n, Bt)
    par = dict(m3_maxrel=float(np.abs(Lm - zm["lfc"]).max() / np.abs(zm["lfc"]).max()),
               prem3_maxrel=float(np.abs(Lp - zp["lfc"]).max() / np.abs(zp["lfc"]).max()))
    log("A1 parity vs Task 14", par); assert par["m3_maxrel"] < 1e-5 and par["prem3_maxrel"] < 1e-5, par
    rep["parity_task14"] = par; rep["sum_side_info"] = tinfo

    # --- A2 DE side: harness delta_s4 (S4 Yorig + 8-member ensemble + 0.3 knn), original (parity) and denoised library
    for kk in ["sch_p1_H1"] + [f"sch_p1_H1_cxmin_cells{v}" for v in SE.MC7] + ["sch_lib_p0_H1", "sch_lib_p1_H1"]:
        assert (H.CACHE / f"{kk}__H1.npz").exists(), f"harness cache key missing (would be WRITTEN): {kk}"
    d_orig, ov = SE.delta_s4(P, "Yorig", 0.3, 0.5, mc_list=SE.MC7)
    par2 = dict(base_total_maxabs=float(np.abs(d_orig - base_total).max()),
                base_total_cos=float((d_orig * base_total).sum() / np.linalg.norm(d_orig) / np.linalg.norm(base_total)))
    Odp = np.full((len(tg), len(genes)), np.nan); Odp[:, ok] = ov[:, gmap[ok]]
    par2["Od_equal"] = bool(np.array_equal(np.isfinite(Odp), np.isfinite(Od)) and np.allclose(np.nan_to_num(Odp), np.nan_to_num(Od), atol=1e-6))
    log("A2 harness s4Yorig_ensP1 vs ERA base_total", par2); rep["parity_harness_recipe_vs_pack"] = par2
    assert par2["base_total_cos"] > 0.999 and par2["Od_equal"], par2
    OV = ov.astype(np.float64)
    vsrc = RP / "sc_hybrid_2026-09-25/denoise_src" / DSRC; vcache = W4 / f"cache_{DSRC}_de"; vcache.mkdir(parents=True, exist_ok=True)
    assert not [f.name for f in vsrc.iterdir() if f.is_symlink() and not f.exists()], "broken symlinks in the denoised library"
    dde_f = W4 / "de_delta_H1.npz"
    if dde_f.exists():
        d_de = np.load(dde_f)["d_de"].astype(np.float64)
    else:
        P0_, C0_, S0_ = H.PARTS, H.CACHE, TM.SRC
        try:
            TM.SRC = vsrc; TM._cache.clear(); H.CACHE = vcache
            d_de, ov_de = SE.delta_s4(P, "Yorig", 0.3, 0.5, mc_list=SE.MC7)
        finally:
            H.PARTS, H.CACHE, TM.SRC = P0_, C0_, S0_; TM._cache.clear()
        # builder: overrides OV (and the budget reference) stay the MAIN recipe's; only the DE delta follows the denoised
        # library (tv2_submit_build: "budget reference + overrides stay the main recipe") -> record, do not assert
        np.savez_compressed(dde_f, d_de=d_de.astype(np.float32), ov=ov.astype(np.float32), ov_de=ov_de.astype(np.float32),
                            targets=np.array(tg))
    zd = np.load(dde_f); ov_de = zd["ov_de"].astype(np.float64)
    fo, fd = np.isfinite(ov), np.isfinite(ov_de)
    rep["de_override_vs_main"] = dict(main_finite=int(fo.sum()), de_finite=int(fd.sum()), pattern_mismatch=int((fo != fd).sum()),
                                      maxabs_on_common=float(np.abs(np.where(fo & fd, ov - ov_de, 0)).max()),
                                      used="main OV (builder convention)")
    log("DE-side override vs main", rep["de_override_vs_main"])
    rep["de_side"] = dict(source_library=str(vsrc), cos_de_vs_sum_recipe=float((d_de * d_orig).sum() / np.linalg.norm(d_de) / np.linalg.norm(d_orig)),
                          reference_harness_dsrc_cos=0.9723)
    log("DE delta", rep["de_side"])

    # --- agc consensus agreement (harness, own family out) on the panel axis
    fl = H.fam_lfc(tg, H.EXCL[P]); g_ = gmap[ok]
    S_ = np.stack([np.nan_to_num(v[:, g_]) for v in fl.values()]); M_ = np.stack([np.isfinite(v[:, g_]) for v in fl.values()])
    mean_ = np.where(M_.sum(0) > 0, S_.sum(0) / np.maximum(M_.sum(0), 1), 0.0)
    fr_ = np.where(M_.sum(0) > 0, (M_ & (np.sign(S_) == np.sign(mean_)[None])).sum(0) / np.maximum(M_.sum(0), 1), 0.5)
    AG = np.zeros((len(tg), len(genes))); AG[:, ok] = np.clip(2 * fr_ - 1, 0, 1); del S_, M_, fl

    o = SE.parse(EMIT_ARM); gate = r["cpm"] >= o.get("gt", 5.0); low = ~gate
    ev = EV.Pack(PACKS / "score_H1.npz")
    zs = TM.load(SL_SRC); rep["same_line_source"] = dict(name=SL_SRC, path=str(TM._src_path(SL_SRC)), n_targets=int(len(zs["targets"])))
    import gc
    del zi; gc.collect()

    def make_arm(a):                                      # lazily, one arm in memory at a time
        if a.startswith("w"):
            w = float(a[1:]); (bt_w, dde_w), ncov = same_line([base_total, d_de], tg, OV, w, TM)
            Xw, _ = sum_side(on_panel(bt_w, gmap), Od, V3, cpm, rk, Vd=Vd, wdf=wdf); s_ = bmatch(fm, Xw, Od, n, Bt)
            if w == 0.0: assert np.abs(Xw * s_ - Lm).max() == 0.0
            return dict(w=w, m2=True, m3=True, Lsum=Xw * s_, Dde=dde_w, scale=s_, covered=ncov)
        if a == "noM2":
            Xn2, _ = sum_side(L0, Od, V3, cpm, None, Vd=Vd, wdf=wdf); s2 = bmatch(fm, Xn2, Od, n, Bt)
            return dict(w=0.0, m2=False, m3=True, Lsum=Xn2 * s2, Dde=d_de, scale=s2, covered=0)
        assert a == "noM3"
        return dict(w=0.0, m2=True, m3=False, Lsum=Lp, Dde=d_de, scale=float(zp["budget_scale"]), covered=0)
    for a in (only or ARMS):
        out = PKG / a; out.mkdir(parents=True, exist_ok=True)
        if (out / "pred.h5ad").exists() and (out / "meta.json").exists():
            log(a, "package exists"); continue
        ta = time.time(); A = make_arm(a); gc.collect()
        LdDE = drc(on_panel(A["Dde"], gmap), Od, *DRC)
        rs = np.abs(LdDE) * AG
        e0 = ev.evaluate(A["Lsum"])
        if os.environ.get("RQ4_SMOKE"):                  # smoke: 3 targets through the emitter + pinned-sum check, nothing written
            sl_ = slice(0, 3); bl_ = list(SE.emit(r["ctrl"], A["Lsum"][sl_], Od[sl_], n, o, gate, low, rs[sl_], P, LdDE[sl_]))
            for t in range(3):
                e = np.exp(A["Lsum"][t]); j = np.isfinite(Od[t]); e[j] = np.exp(Od[t, j])
                assert np.array_equal(np.asarray(bl_[t].sum(0)).ravel(), np.round(fm * e * n)), "smoke: pinned sums"
            log("SMOKE_OK", a, dict(pds_offline=e0["pds"], covered=A["covered"], scale=A["scale"], nnz=[int(b.nnz) for b in bl_],
                                    secs=round(time.time() - ta, 1)))
            shutil.rmtree(out, ignore_errors=True); return
        blocks = list(SE.emit(r["ctrl"], A["Lsum"], Od, n, o, gate, low, rs, P, LdDE))
        # pinned sums check (PDS / MSE side): emitted column totals == round(fm exp(L) n), override entries exp(Od)
        for t in (0, len(tg) // 2, len(tg) - 1):
            e = np.exp(A["Lsum"][t]); j = np.isfinite(Od[t]); e[j] = np.exp(Od[t, j])
            assert np.array_equal(np.asarray(blocks[t].sum(0)).ravel(), np.round(fm * e * n)), f"{a}: pinned sums broken at target {t}"
        nrow = len(tg) * n + r["ctrl"].shape[0]
        X = sparse.vstack(blocks + [r["ctrl"].astype(np.float32)], format="csr"); del blocks; gc.collect()   # ONE copy (memory)
        pred = ad.AnnData(X=X, obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(tg, n), *["non-targeting"] * r["ctrl"].shape[0]])},
                                                index=[f"p{i}" for i in range(nrow)]),
                          var=pd.DataFrame(index=list(genes)))
        pred.write_h5ad(out / "pred.h5ad"); del pred, X; gc.collect()
        np.savez_compressed(out / "inputs.npz", Lsum=A["Lsum"].astype(np.float32), LdDE=LdDE.astype(np.float32), rs=rs.astype(np.float32),
                            Od=Od.astype(np.float32), targets=np.array(tg), genes=genes)
        meta = dict(arm=a, w=A["w"], m2=A["m2"], m3=A["m3"], same_line_covered=A["covered"], budget_scale=A["scale"],
                    pds_offline=e0["pds"], mse_offline=e0["mse"], emit_arm=EMIT_ARM, emit_seed=o["seed"], de_drc=DRC,
                    de_source=DSRC, seconds_emit=round(time.time() - ta, 1),
                    sha_note="pred.h5ad deterministic given inputs.npz + controls (rng 0)")
        jdump(out / "meta.json", meta); log(a, {k: v for k, v in meta.items() if k not in ("sha_note",)})
        rep.setdefault("arms", {})[a] = meta; del A, LdDE, rs; gc.collect()
    old = json.loads((W4 / "build_report.json").read_text()) if (W4 / "build_report.json").exists() else {}
    old.update({k: v for k, v in rep.items() if k != "arms"}); old.setdefault("arms", {}).update(rep.get("arms", {}))
    jdump(W4 / "build_report.json", old); log("BUILD_DONE")


# ----------------------------------------------------------------------------- real scorer
def cmd_real(arms):
    import pandas as pd
    from cell_eval2.moments import jackknife_correction  # noqa: F401  (import order: BEFORE scripts/ joins sys.path)
    sys.path.append(str(SCR)); import tv2_score as S
    truth = RP / "transfer_v2_2026-09-23/h1/real_thin.h5ad"
    for a in arms:
        out = PKG / a; res_f = out / "result.json"
        if res_f.exists():
            log(a, "already scored"); continue
        pp = out / "pred.h5ad"; assert pp.exists(), f"package missing: {pp}"
        ta = time.time()
        cmd = [str(CE), "run", "-ap", str(pp), "-ar", str(truth), "--preset", "vcc2026", "--pert-col", "target",
               "--control", "non-targeting", "--input-type", "counts", "-o", str(out / "cell_eval2"), "--cache-real", str(CACHE_REAL),
               "--set", "de.backend=pdex", "--set", "num_threads=4", "--set", "gather_threads=4"]
        env = dict(os.environ, NUMBA_NUM_THREADS="4", OMP_NUM_THREADS="4")
        with (out / "scorer.log").open("w") as g:
            subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
        agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
        raw = {k: float(agg.loc["mean", k]) for k in S.CONST}; comp, scaled = S.composite(raw)
        shutil.copy(out / "cell_eval2/agg_results.csv", out / "agg_results.csv")
        meta = json.loads((out / "meta.json").read_text())
        res = dict(arm=a, raw=raw, scaled=scaled, composite=comp, seconds=round(time.time() - ta, 1), meta=meta,
                   scorer=" ".join(cmd[1:]))
        jdump(res_f, res); shutil.rmtree(out / "cell_eval2", ignore_errors=True)
        log(a, "REAL composite", round(comp, 5), {k[:14]: round(v, 4) for k, v in raw.items()}, "secs", res["seconds"])
    log("REAL_DONE")


# ----------------------------------------------------------------------------- verdict
def cmd_verdict():
    import pandas as pd
    rows = []
    for a in ARMS:
        m = json.loads((PKG / a / "meta.json").read_text()) if (PKG / a / "meta.json").exists() else {"arm": a}
        rf = PKG / a / "result.json"; r = json.loads(rf.read_text()) if rf.exists() else None
        row = dict(arm=a, w=m.get("w"), m2=m.get("m2"), m3=m.get("m3"), same_line_covered=m.get("same_line_covered"),
                   pds_offline=m.get("pds_offline"), mse_offline=m.get("mse_offline"), status="scored" if r else "queued",
                   composite=(r or {}).get("composite"))
        for k, v in ((r or {}).get("raw") or {}).items(): row[f"raw_{k}"] = v
        for k, v in ((r or {}).get("scaled") or {}).items(): row[f"scaled_{k}"] = v
        rows.append(row)
    df = pd.DataFrame(rows); df.to_parquet(RES / "rq4_sameline.parquet", index=False)
    comp = {r["arm"]: r["composite"] for r in rows}
    # A6: best architecture
    fr = json.loads(ARCH_JSON.read_text())["rows"]
    mean_ = {x["arch"]: 0.5 * (x["K562_c"] + x["HepG2_c"]) for x in fr}
    cells = {(x["arch"], ln): x[f"{ln}_c"] for x in fr for ln in ("K562", "HepG2")}
    best_mean = max(mean_, key=mean_.get); best_cell = max(cells, key=cells.get)
    arch = dict(n_architectures=len(fr), primary=dict(statistic="line-mean (K562, HepG2) official composite delta vs refit_v2, alpha 1, correct context, seed 0",
                                                        arch=best_mean, gain=mean_[best_mean]),
                sensitivity=dict(statistic="largest single-line value (max over 34 line x architecture cells)",
                                 arch=best_cell[0], line=best_cell[1], gain=cells[best_cell]),
                setup_differences=["panels K562 / HepG2 (not H1)", "base = refit_v2 extended transfer (not the fusion pipeline)",
                                   "emission / budget = dcp preset path (not the F2 dual emission)", "single seed 0; Task 10 re-checks with seeds",
                                   "same scorer preset (cell-eval2 vcc2026)"])
    have = all(comp.get(a) is not None for a in VERDICT_ARMS)
    g = lambda x, y: (comp[x] - comp[y]) if comp.get(x) is not None and comp.get(y) is not None else None
    gains = dict(M2=g("w0", "noM2"), M3=g("w0", "noM3"),
                 **{f"H1GW_w{w:g}": g(f"w{w:g}", "w0") for w in WS[1:]})
    ref = {}
    for sd in ("", "@s1", "@s2"):
        f = RP / "sc_hybrid_2026-09-25/emit" / (REF_SL + sd) / "H1/result.json"
        if f.exists(): ref[sd or "@s0"] = json.loads(f.read_text())["composite"]
    out = dict(plan="input/ideas/zero-shot-limits-context-complexity-plan.md", task=18, hypothesis="H4",
               rule="Rejected if composite(base + H1GW, w = 1) - composite(base) <= the largest single modelling gain (M2, M3, best architecture), real scorer, H1 panel",
               composites=comp, gains=gains, best_architecture=arch,
               per_context_value=None, per_context_value_note="Task 11 (LOFO per-context value, converted to composite on H1) not yet run; fill from Task 11",
               H1GW_per_context_ratio=None,
               reference_emission_seed_spread=dict(arm=REF_SL, composites=ref,
                                                    range=(max(ref.values()) - min(ref.values())) if len(ref) > 1 else None,
                                                    note="earlier harness same-line w1 arm (different sum side); context for single-run noise only"),
               related_line_cases=dict(status="skipped", reason="lowest priority (cut first); the scorer queue is occupied by Task 17 and "
                                                              "then the 7 H1 runs (~28 min each); not started in this dispatch"))
    if have:
        mod = dict(M2=gains["M2"], M3=gains["M3"], arch_primary=arch["primary"]["gain"])
        lg = max(mod, key=mod.get); verdict = "rejected" if gains["H1GW_w1"] <= mod[lg] else "confirmed"
        mod_s = dict(mod, arch_primary=arch["sensitivity"]["gain"]); lg_s = max(mod_s, key=mod_s.get)
        out.update(status="final", largest_modelling_gain=dict(source=lg, value=mod[lg]), verdict=verdict,
                   sensitivity=dict(largest_modelling_gain=dict(source=lg_s, value=mod_s[lg_s]),
                                    verdict="rejected" if gains["H1GW_w1"] <= mod_s[lg_s] else "confirmed"),
                   margin=gains["H1GW_w1"] - mod[lg])
    else:
        out.update(status="provisional", verdict=None, missing=[a for a in VERDICT_ARMS if comp.get(a) is None],
                   note="verdict requires the real-scorer runs w0, w1, noM2, noM3")
    jdump(RES / "rq4_h4.json", out); log("VERDICT", json.dumps({k: out[k] for k in ("status", "verdict", "gains")}, default=str))


if __name__ == "__main__":
    cmd = sys.argv[1]; args = sys.argv[2:]
    if cmd == "build":
        cmd_build(args or None)
    elif cmd == "real":
        cmd_real(args or ARMS)
    elif cmd == "verdict":
        cmd_verdict()
    else:
        raise SystemExit(f"unknown subcommand {cmd}")
