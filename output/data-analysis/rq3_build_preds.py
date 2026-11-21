"""Task 14 (plan `input/ideas/zero-shot-limits-context-complexity-plan.md`, RQ3): frozen pre-M3 / M3 predictions on the
7 harness panels, leave-destination-out, plus the KOLF leak fix (Amendment A3) and the control-dependence check.

Subcommands (run in this order; each writes under output/data-analysis/results/preds/):
  feng     [vcc2026 env]     kolf-free FENG_T / FENG_GWF / FENG_GWN source stores, rebuilt from the cell-level Feng
                             pseudobulk sums per (target x iPSC line) (pseudobulk/feng_*_target_line, sums of raw UMI over
                             cells; dropping the kolf_2 / kolf_3 rows == dropping those cells). Parity: rebuilding with ALL
                             lines reproduces the harness stores (final_harness_2026-09-24/sources/FENG_*.npz).
  s4kolf   [ct_transfer env] S4 for the KOLF destination without kolf cells: refit_v2 exactly as ext_baselines.refit, with the
                             POOLED Feng rows removed (table sources xl/feng, fengpooled/fitness and the pooled Feng
                             empirical stores contain kolf_2 / kolf_3 cells); the per-line non-kolf Feng rows are kept.
                             S4 = refit + 0.25 x the saved FiLM residual (film_ext_a1 - refit_v1), the s4i_refit_excl.py
                             convention. Parity: the unmodified refit reproduces the saved refit_v2.npy.
  build    [vcc2026 env]     per panel: pre-M3 and M3 sum-side predictions, M2 ranks, sch_sumeval-type sanity check,
                             control-dependence test, three-way control splits (seeds 0-4) and C_a u C_b pre-M3 variants.

ASSUMPTIONS (stated; the checks are in the code and in results/preds/rq3_build_preds_summary.json):
  A1  "Frozen recipe" on a harness panel = the harness M1 sum side (ERA pack `inputs_<P>.npz:base_total`: fixed-weight
      family fusion over the 8 min_cells members, 0.5 x S4, 0.3 DepMap-neighbour blend, energy-matched; destination family
      left out, panels.json `excluded_sources`) + F2's sum-side settings (tv2_submit_build.py sum_side order):
      ew (1 + cpm/50)^-0.25 [F2 ew_adaptive false -> always on] -> centring 0.5 -> guarded top-3 library-PC removal
      (lambda 0.75, guard 0.15/0.25) -> M2 = F2 `--sum-tamp 0.5`: rank(cross-family agreement)^0.5 -> [M3] -> budget
      match to the base recipe's budget (kink factor 1.17, margin 0.8, floor 1). NOT mirrored on the harness (F2 M1 items
      that exist only for the official contexts): `--sum-cens` pseudocount members, the DLD1 family, promoter-off
      fusion, knn_noK, the learned `--sum-resid 0.25` residual, same-line replacement, S4noFK (somatic-only branch).
  A2  M2 ranks: agreement_t = mean pairwise cosine of the donor-family mean LFCs (c 50, centred, min_cells 20, families
      measuring > 500 genes of t) over the panel's fusion families (W_NEW > 0, minus the destination's family, minus KOLF
      for somatic panels) -- tv2_tamp_prep.agreement verbatim; rank = (stable rank + 1)/T. The historical harness M2
      (S16/S18: learned tlrn rho_hat, beta 0.75) is NOT used; it is reported only as a reference sanity arm.
  A3  M3 = builder `sum.ctrl_defl` {k 8, eta 0.75, weak 0.5, beta 0.75, ncells 6000, seed 0}: V = top-8 right singular
      vectors (randomized_svd n_iter 5, random_state 0) of 6,000 random (rng 0) destination control cells in log1p CP10K
      on the OFFICIAL axis (panel genes mapped with the harness map; library size over mapped genes), gene-centred;
      X <- X - 0.75 (X V^T) V on targets whose agreement amplitude is <= its median.
  A4  Pseudobulk predicted delta (deliverable) = log1p(CP10K(fm * exp(L))) - log1p(CP10K(fm)), fm = mean control counts per
      panel gene, override entries (own gene / promoter floor) multiplicative as in the scorer, NO rounding (the scorer
      rounds expected counts; the delta here is the noise-free expectation). Gene axis = the truth file's var_names.
  A5  The offline scorer is EXACT for M3 vs pre-M3 because both are budget-matched to the same budget (AGENTS.md
      2026-10-03 rule); it reads the ERA pack's real moments (`score_<P>.npz`), numerically = sch_sumeval.evaluate.
  A6  Three-way control split: common.split_halves(n_ctrl, seed, 3) -> (C_pc, C_a, C_b), seeds 0-4; indices are positions
      among the truth file's control cells (file order) plus the absolute obs rows.
Seeds: everything deterministic (rng 0 for ctrl-PC cell draw and kink jackknife draws, as in the frozen builder;
split seeds 0-4)."""
from __future__ import annotations
import json, os, pathlib, sys, time

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); V = ROOT / "vcc2026_data"; RP = V / "reports"; SCR = V / "scripts"
OUT = ROOT / "output/data-analysis/results/preds"; FIX = OUT / "kolf_fix"
LOG = ROOT / "output/data-analysis/results/logs"
PACKS = RP / "era_2026-10-07/packs"; HSRC = RP / "final_harness_2026-09-24/sources"
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
KOLF_LINES = ("kolf_2", "kolf_3")
FENG_ARMS = [("FENG_T", "feng_targetedscreen"), ("FENG_GWF", "feng_genomewidescreen_fitnessgenes"),
             ("FENG_GWN", "feng_genomewidescreen_nonfitnessgenes")]
POOLED_FENG_TABLE_SOURCES = ("xl/feng", "fengpooled/fitness")
RECIPE = ("harness M1 = ERA pack base_total (S14fk-structure: fixed-weight fusion x 8 min_cells members, 0.5 S4, 0.3 knn "
          "blend, energy-matched; leave-destination-out) | F2 sum side: ew 50,0.25 (always on) -> center 0.5 -> pcr 3x0.75g "
          "(guard 0.15/0.25) -> tamp beta 0.5 on rank(cross-family agreement) [M2, F2 --sum-tamp 0.5] -> "
          "[M3: ctrl_defl k8 eta0.75 weak0.5 beta0.75 ncells6000 seed0 = --sum-defl 8,0.75,0.5] -> budget match to the "
          "base budget (kink 1.17, margin 0.8, floor 1)")
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def jdump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)))


# ======================================================================================= feng (vcc2026 env)
def cmd_feng():
    import numpy as np, pandas as pd
    from scipy import sparse
    sys.path.append(str(SCR)); import tv2_sources as TS
    src = FIX / "sources_kolffree"; src.mkdir(parents=True, exist_ok=True)
    kolf_targets = json.loads((RP / "final_harness_2026-09-24/panels.json").read_text())["KOLF"]["targets"]
    kt = {TS.canon(t) for t in kolf_targets}
    report = {}

    def rebuild(store, drop):
        o, M, o2s = TS.load_store(store + "_target_line")
        is_c = o["is_control"].astype(str).str.lower().isin(["true", "1"]).to_numpy()
        keep = ~o["donor"].astype(str).isin(drop).to_numpy()
        tgt = np.array([TS.canon(x) for x in o["target_symbol"].astype(str)])
        Tot = o["total_counts"].to_numpy(np.float64); nc = o["n_cells"].to_numpy(np.int64)
        rows_c = np.flatnonzero(keep & is_c); rows_p = np.flatnonzero(keep & ~is_c)
        keys = sorted(set(tgt[rows_p])); kix = {k: i for i, k in enumerate(keys)}
        G = sparse.csr_matrix((np.ones(len(rows_p)), ([kix[t] for t in tgt[rows_p]], rows_p)), shape=(len(keys), M.shape[0]))
        S = (G @ M.astype(np.float64)).toarray()                          # exact integer sums (float64)
        Sc = np.asarray(M[rows_c].astype(np.float64).sum(0)).ravel()
        okm = o2s >= 0

        def off(v):
            w = np.full(TS.NG, np.nan); w[okm] = v[o2s[okm]]; return w
        fc = off(Sc) / Tot[rows_c].sum()
        Tt = np.asarray(G @ Tot).ravel(); Nt = np.asarray(G @ nc.astype(np.float64)).ravel().astype(np.int64)
        return keys, S, Tt, Nt, fc, off, dict(ctrl_cells=int(nc[rows_c].sum()), pert_cells=int(nc[rows_p].sum()),
                                              ctrl_rows=int(len(rows_c)), lines=sorted(set(o["donor"].astype(str)[keep])))

    for name, store in FENG_ARMS:
        old = np.load(HSRC / f"{name}.npz", allow_pickle=True); need = set(str(t) for t in old["targets"])
        res = {}
        for tag, drop in (("all_lines", ()), ("kolf_free", KOLF_LINES)):
            keys, S, Tt, Nt, fc, off, meta = rebuild(store, drop)
            acc = TS.Acc(need)
            for i, k in enumerate(keys):
                if Tt[i] > 0:
                    acc.add(k, off(S[i]) / Tt[i], fc, int(Nt[i]))
            ks = sorted(acc.fp)
            arr = dict(targets=np.array(ks), fp=np.vstack([acc.fp[k] for k in ks]).astype(np.float32),
                       fc=np.vstack([acc.fc[k] for k in ks]).astype(np.float32), ncell=np.array([acc.n[k] for k in ks]),
                       **{f"center_c{c:g}": np.where(acc.cc[c] > 0, acc.cs[c] / np.maximum(acc.cc[c], 1), np.nan) for c in TS.CS},
                       **{f"ncenter_c{c:g}": acc.cc[c] for c in TS.CS})
            on_kolf = int(sum(acc.n[k] for k in ks if k in kt and acc.n[k] >= 20))
            meta.update(n_targets_stored=len(ks), kolf_panel_target_cells_ge20=on_kolf,
                        kolf_panel_targets_ge20=int(sum(1 for k in ks if k in kt and acc.n[k] >= 20)))
            res[tag] = meta
            if tag == "all_lines":                                      # parity vs the harness store
                assert list(arr["targets"]) == [str(t) for t in old["targets"]], name
                assert np.array_equal(arr["ncell"], old["ncell"]), name
                for k in ("fp", "fc", "center_c50"):
                    a, b = arr[k].astype(np.float64), old[k].astype(np.float64); m = np.isfinite(b)
                    assert np.array_equal(np.isfinite(a), m), (name, k)
                    rel = float(np.abs(a[m] - b[m]).max() / max(np.abs(b[m]).max(), 1e-30))
                    meta[f"parity_maxrel_{k}"] = rel; assert rel < 1e-5, (name, k, rel)
            else:
                np.savez(src / f"{name}.npz", **arr)
                (src / f"{name}.json").write_text(json.dumps(dict(store=store + "_target_line", dropped_lines=list(drop), **meta), indent=1))
            log(name, tag, {k: v for k, v in meta.items() if k != "lines"})
        res["removed_pert_cells"] = res["all_lines"]["pert_cells"] - res["kolf_free"]["pert_cells"]
        res["removed_ctrl_cells"] = res["all_lines"]["ctrl_cells"] - res["kolf_free"]["ctrl_cells"]
        res["removed_pert_frac"] = round(res["removed_pert_cells"] / res["all_lines"]["pert_cells"], 4)
        report[name] = res
    for f in HSRC.glob("*.npz"):                                        # every other harness source unchanged (symlink)
        if not f.name.startswith("FENG_") and not (src / f.name).exists():
            (src / f.name).symlink_to(f)
    jdump(FIX / "feng_rebuild.json", report)
    log("FENG_DONE", FIX / "feng_rebuild.json")


# ======================================================================================= s4kolf (ct_transfer env)
def cmd_s4kolf():
    import gc
    import numpy as np
    sys.path.insert(0, str(SCR)); import ext_baselines as XB
    line = "KOLF"; EXT = XB.ROOT
    base = XB.A.Data(); targets = base.targets[line]
    hp = json.loads((RP / "final_harness_2026-09-24/panels.json").read_text())["KOLF"]["targets"]
    assert [str(t) for t in targets] == [str(t) for t in hp], "S4 target order != harness KOLF panel"
    canon = [XB.D.canonical_target(t) for t in targets]; feats = base.features(targets)
    indices = {g: i for i, g in enumerate(XB.V.official_features())}
    table = XB.D.Table.load(EXT / "table_v2"); columns = np.array([indices[g] for g in table.genes])
    srcs = np.asarray(table.sources).astype(str)

    def refit(drop_sources):
        ids = XB.F.permitted_ids(table, [line])
        conf = XB.CONFLICT_GENE.get(XB.D.canonical_bio(line), [])
        if conf:
            ids = ids[~np.isin(np.asarray(table.bios)[ids].astype(str), conf)]
        n0 = len(ids)
        if drop_sources:
            ids = ids[~np.isin(srcs[ids], list(drop_sources))]
        bios_ = np.asarray(table.bios)[ids].astype(str)
        held = {XB.D.canonical_bio(x) for x in bios_}
        assert XB.D.canonical_bio(line) not in held
        assert not np.isin(srcs[ids], ["fengline/kolf_2", "fengline/kolf_3"]).any()
        y = np.asarray(table.y[ids], np.float32)
        meta = {k: np.asarray(getattr(table, k))[ids] for k in ('targets', 'bios', 'modality', 'features', 'weights', 'sources')}
        sub = XB.D.Table(y=y, controls=np.zeros((len(ids), 1), np.float32), genes=list(table.genes), **meta)
        fitted = XB.LB.LeanSourceBaseline(sub, np.arange(len(ids)), rank=32)
        common, seen = fitted.predict(canon, feats)
        nfeng = dict(zip(*np.unique(srcs[ids][bios_ == "FENG"], return_counts=True)))
        del sub, fitted, y, meta; gc.collect()
        return common.astype(np.float32), dict(rows=int(len(ids)), rows_before_drop=int(n0), bios=len(held),
                                               feng_rows_by_source={str(k): int(v) for k, v in nfeng.items()})

    rep = {}
    # parity: the unmodified path reproduces the saved refit_v2 / film_v2_a0.25
    emp0, nv0 = XB.empirical(targets, [line], "v2"); c0, m0 = refit(())
    pred0 = emp0.copy(); pred0[:, columns] = c0
    ref = np.load(EXT / line / "refit_v2.npy").astype(np.float64)
    rep["parity_refit_v2_maxabs"] = float(np.abs(pred0 - ref).max()); rep["parity_refit_v2_cos"] = float((pred0 * ref).sum() / np.linalg.norm(pred0) / np.linalg.norm(ref))
    res = np.load(EXT / line / "film_ext_a1.npy").astype(np.float64) - np.load(EXT / line / "refit_v1.npy").astype(np.float64)
    f_ref = np.load(EXT / line / "film_v2_a0.25.npy").astype(np.float64)
    rep["parity_film_v2_from_saved_maxabs"] = float(np.abs(ref + 0.25 * res - f_ref).max())
    rep["refit_original"] = m0; log("parity", rep)
    # fixed: pooled Feng (contains kolf_2/kolf_3 cells) out of the refit rows AND the empirical means; per-line non-kolf kept
    B = XB.B; meta0 = B.metadata; nm0 = XB.new_meta

    def meta_(store):
        p = meta0(store).copy()
        if store == "feng_genomewidescreen_nonfitnessgenes_target":
            p["bio"] = "FENGPOOL"
        return p

    def new_meta_(store):
        p = nm0(store)
        if store == "feng_fit":
            p["bio"] = "FENGPOOL"
        return p
    B.metadata, XB.new_meta = meta_, new_meta_
    try:
        emp1, nv1 = XB.empirical(targets, [line, "FENGPOOL"], "v2")
    finally:
        B.metadata, XB.new_meta = meta0, nm0
    c1, m1 = refit(POOLED_FENG_TABLE_SOURCES)
    pred1 = emp1.copy(); pred1[:, columns] = c1
    film1 = pred1.astype(np.float64) + 0.25 * res
    np.save(FIX / "s4_KOLF_refit_nopooledFENG.npy", pred1.astype(np.float32))
    np.save(FIX / "s4_KOLF_film_nopooledFENG.npy", film1.astype(np.float32))
    rep.update(refit_fixed=m1, empirical_voters_original=int(nv0), empirical_voters_fixed=int(nv1),
               cos_film_fixed_vs_original=float((film1 * f_ref).sum() / np.linalg.norm(film1) / np.linalg.norm(f_ref)),
               film_residual_note="saved FiLM residual (film_ext_a1 - refit_v1) kept, weight 0.25 (s4i_refit_excl convention); "
                                  "it was trained with the kolf HipSci lines held out as contexts but with pooled Feng rows "
                                  "present -> a residual kolf exposure through 0.25 x residual x 0.5 S4 share remains")
    jdump(FIX / "s4kolf.json", rep); log("S4KOLF_DONE", rep)


# ======================================================================================= build (vcc2026 env)
def _imports():
    from cell_eval2.moments import jackknife_correction  # noqa: F401  (import order rule: BEFORE scripts/ joins sys.path)
    from cell_eval2.prep import bulk_lognorm_means  # noqa: F401
    sys.argv = ["x"]; sys.path.append(str(SCR)); sys.path.insert(0, str(ROOT / "output/data-analysis"))
    sys.path.insert(0, str(ROOT / "external/era-aivc/aivc"))
    import ml_h1fuse as F, fh_harness as H, sch_emit as EMT, tv2_emit as EM, tv2_model as TM
    import common as CM
    from tasks.vcc_sum import evaluator as EV
    return F, H, EMT, EM, TM, CM, EV


def agreement(Fm, fams, T):
    """tv2_tamp_prep.agreement, verbatim."""
    import numpy as np
    ag = np.full(T, np.nan)
    for k in range(T):
        cols = [f for f in fams if f in Fm and np.isfinite(Fm[f][k]).sum() > 500]; cs = []
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                a, b = Fm[cols[i]][k], Fm[cols[j]][k]; mm = np.isfinite(a) & np.isfinite(b)
                cs.append(float(a[mm] @ b[mm] / (np.linalg.norm(a[mm]) * np.linalg.norm(b[mm]) + 1e-12)))
        if cs: ag[k] = float(np.mean(cs))
    return ag


def rank01(ag):
    import numpy as np
    ag = np.nan_to_num(np.asarray(ag, np.float64), nan=float(np.nanmedian(ag)))
    return (np.argsort(np.argsort(ag, kind="stable"), kind="stable") + 1) / len(ag)


def build_main(panels, splits=range(5), nofilm=False):
    import numpy as np
    from scipy import sparse
    from sklearn.utils.extmath import randomized_svd
    F, H, EMT, EM, TM, CM, EV = _imports()
    from cell_eval2.prep import bulk_lognorm_means as BLM
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads((OUT / "rq3_build_preds_summary.json").read_text()) if (OUT / "rq3_build_preds_summary.json").exists() else {}

    # ------------------------------------------------------------------ helpers
    def ctrl_pcs(ctrl, gmap, ok, k=8, ncells=6000, seed=0):
        """builder sum.ctrl_defl V on the OFFICIAL axis (panel genes mapped), returned on the panel axis."""
        okc = np.flatnonzero(ok); g = gmap[okc]; first = np.unique(g, return_index=True)[1]; okc, g = okc[first], g[first]
        M = sparse.csr_matrix(ctrl, dtype=np.float64)[:, okc].tocoo()
        Mo = sparse.csr_matrix((M.data, (M.row, g[M.col])), shape=(ctrl.shape[0], H.NG))
        L_ = np.asarray(Mo.sum(1)).ravel(); Xc = Mo.multiply(1e4 / np.maximum(L_, 1)[:, None]).tocsr(); Xc.data = np.log1p(Xc.data)
        ix = np.random.default_rng(seed).choice(Xc.shape[0], min(ncells, Xc.shape[0]), replace=False)
        Dc = np.asarray(Xc[ix].todense()); Dc -= Dc.mean(0)
        Vo = randomized_svd(Dc, n_components=k, n_iter=5, random_state=0)[2].astype(np.float64); del Dc
        Vp = np.zeros((k, len(gmap))); Vp[:, ok] = Vo[:, gmap[ok]]
        return Vp, int(len(g)), int(ok.sum() - len(g))

    def Bfun(fm, L, Od, n):
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

    def budget(ctrl, L0, Od, n):
        gk, cr = EMT.kink_gamma(ctrl, L0, Od, n, factor=1.17, margin=0.8)
        fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]
        return Bfun(fm, L0, Od, n)(max(1.0, gk)), float(gk), float(cr)

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

    def delta_lcp(L, Od, fm):
        j = np.isfinite(Od); e = np.exp(L); e[j] = np.exp(Od[j]); S = fm[None, :] * e
        return (np.log1p(S / S.sum(1, keepdims=True) * 1e4) - np.log1p(fm / fm.sum() * 1e4)[None, :]).astype(np.float32)

    def cpm_of(ctrl, ok):
        tot = np.asarray(ctrl.sum(0)).ravel().astype(np.float64) * ok
        return tot / tot.sum() * 1e6

    # ------------------------------------------------------------------ KOLF M1 rebuild (A3)
    def kolf_m1(src_dir, cache_dir, s4_path, tg, ctrl, genes, gmap):
        """ml_h1fuse.pieces('KOLF') + assemble's member means + build, with a chosen source library / S4 / cache.
        Returns base_total, ov, Od (panel), V3, Bt and the agreement inputs. Nothing is written outside cache_dir."""
        p = "KOLF"; P0_, C0_, S0_, L0_ = H.PARTS, H.CACHE, TM.SRC, H.LIB
        try:
            TM.SRC = pathlib.Path(src_dir); TM._cache.clear()
            H.LIB = sorted({str(t) for s in TM.FAMILY if (TM.SRC / f"{s}.npz").exists() for t in np.load(TM.SRC / f"{s}.npz")["targets"]})
            H.CACHE = pathlib.Path(cache_dir)
            if H.CACHE == C0_:                      # original cache is READ-ONLY here
                for kk in ("lib_KOLF__KOLF", "sch_lib_p0_KOLF__KOLF", "sch_lib_p1_KOLF__KOLF"):
                    assert (C0_ / f"{kk}.npz").exists(), kk
            else:
                H.CACHE.mkdir(parents=True, exist_ok=True)
            outer = F.OUTER[p]; wbase = {f: F.W_NEW[f] for f in F.FAMS}
            for f in outer: wbase[f] = 0.0
            ex = H.EXCL[p]; c1 = dict(H.PARTS[1][1]); a1 = H.PARTS[1][2]
            cfg = TM.Cfg()
            for k_, v_ in c1.items(): setattr(cfg, k_, v_)
            _, prom, tidx = TM.predict_lfc(tg, cfg, exclude=ex)
            lib = [None, None]
            for k, part in enumerate(P0_):
                H.PARTS = [(1.0, part[1], part[2])]; lib[k] = H.base(H.LIB, ex, f"sch_lib_p{k}_{p}")[0]
            H.PARTS = P0_
            dlp = H.lib_delta(p, set(ex))
            S3 = 3.0 * np.load(s4_path).astype(np.float64)
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
            del dlib, lib
            fin = [f for f in F.FAMS if f not in outer and wbase[f] > 0]
            ok = gmap >= 0; gpan = gmap[ok]
            Dl = np.zeros((dlp.shape[0], len(genes))); Dl[:, ok] = dlp[:, gpan]; Dl = Dl[np.abs(Dl).sum(1) > 0]
            V3 = randomized_svd(Dl - Dl.mean(0), n_components=3, random_state=0)[2]; del Dl, dlp
            pc = dict(p=p, tg=tg, cfg=cfg, prom=prom, tidx=tidx, a1=a1, S3=S3, Kn=Kn, has=has)
            members = []; Fm20 = None
            for mi, mc in enumerate(F.MC):
                Xm = F.fam_means(tg, fin, 20 if mc is None else mc)
                if mi == 0: Fm20 = Xm
                xs = np.stack([np.nan_to_num(Xm.get(f, np.full((len(tg), H.NG), np.nan))) for f in fin], -1).astype(np.float32)
                ms = np.stack([np.isfinite(Xm.get(f, np.full((len(tg), H.NG), np.nan))) for f in fin], -1)
                a = np.array([wbase[f] for f in fin]) * ms; A = a.sum(-1)
                members.append(np.where(A > 0, (a * xs).sum(-1) / np.maximum(A, 1e-300), 0.0)); del xs, ms, a
                if mi > 0: del Xm
            tot, ov, _ = F.build(pc, members)
            L0 = np.zeros((len(tg), len(genes))); L0[:, ok] = tot[:, gpan]
            Od = np.full(L0.shape, np.nan); Od[:, ok] = ov[:, gpan]
            Bt, gk, cr = budget(ctrl, L0, Od, CM_n["KOLF"])
            ag = agreement(Fm20, fin, len(tg))
            fam_cells = {f: int(np.isfinite(Fm20[f]).any(1).sum()) for f in fin if f in Fm20}
            return dict(base_total=tot, ov=ov, L0=L0, Od=Od, V3=V3, Bt=Bt, kink=gk, c_real=cr, agree=ag, fin=fin,
                        wb=[wbase[f] for f in fin], lib_n=len(H.LIB), targets_with_family=fam_cells)
        finally:
            H.PARTS, H.CACHE, TM.SRC, H.LIB = P0_, C0_, S0_, L0_; TM._cache.clear()

    CM_n = {}
    for p in panels:
        tp = time.time(); log("=====", p)
        sp = np.load(PACKS / f"score_{p}.npz"); n = int(sp["n"]); CM_n[p] = n
        zi = np.load(PACKS / f"inputs_{p}.npz", allow_pickle=False)
        tg = [str(t) for t in zi["targets"]]; fin = [str(f) for f in zi["families"]]
        ctrl, genes = CM.controls(p); genes = np.asarray(genes)
        gmap = CM.to_official(genes); ok = gmap >= 0
        assert np.array_equal(gmap, sp["gmap"]), "gene map != pack"
        fm_full = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]
        assert np.allclose(fm_full, sp["fm"], rtol=1e-6, atol=1e-9), "control mean != pack fm (different control cells?)"
        rec = dict(panel=p, n_targets=len(tg), n_genes_panel=len(genes), n_ctrl=int(ctrl.shape[0]), cells_per_target=n,
                   families=fin, recipe=RECIPE)
        pk = EV.Pack(PACKS / f"score_{p}.npz")
        if p == "KOLF":
            src_fix = FIX / "sources_kolffree"
            s4_fix = FIX / ("s4_KOLF_refit_nopooledFENG.npy" if nofilm else "s4_KOLF_film_nopooledFENG.npy")
            assert (src_fix / "FENG_T.npz").exists() and s4_fix.exists(), "run `feng` and `s4kolf` first"
            if not nofilm:
                orig = kolf_m1(HSRC, H.CACHE, EMT.S4_EXT / "KOLF/film_v2_a0.25.npy", tg, ctrl, genes, gmap)
                par = dict(base_total_maxabs=float(np.abs(orig["base_total"] - zi["base_total"]).max()),
                           V3_absdot_min=float(np.abs((orig["V3"] * sp["V3"]).sum(1)).min()),
                           Bt_rel=float(abs(orig["Bt"] - float(sp["Bt"])) / float(sp["Bt"])),
                           Od_equal=bool(np.array_equal(np.isfinite(orig["Od"]), np.isfinite(sp["Od"])) and
                                         np.allclose(np.nan_to_num(orig["Od"]), np.nan_to_num(sp["Od"]), atol=1e-6)))
                log("KOLF rebuild parity (original sources/S4 vs ERA pack)", par)
                assert par["base_total_maxabs"] < 1e-4 and par["V3_absdot_min"] > 0.999 and par["Bt_rel"] < 1e-5 and par["Od_equal"], par
                fx = kolf_m1(src_fix, FIX / "cache_kolffree", s4_fix, tg, ctrl, genes, gmap)
                fr = json.loads((FIX / "feng_rebuild.json").read_text()); s4r = json.loads((FIX / "s4kolf.json").read_text())
                cosb = float((fx["base_total"] * orig["base_total"]).sum() / np.linalg.norm(fx["base_total"]) / np.linalg.norm(orig["base_total"]))
                rec["kolf_fix"] = dict(
                    method="rebuild (A3 primary option): kolf_2/kolf_3 cells removed from every FENG arm -- fusion family column "
                           "(8 min_cells members), equal-weight first part + DepMap-neighbour library + library PCs (all via the "
                           "kolf-free FENG_T/GWF/GWN stores rebuilt from per-line Feng pseudobulk sums), M2 agreement, and S4 "
                           "(refit without the pooled Feng rows; per-line non-kolf Feng rows kept; saved FiLM residual kept)",
                    feng_cells_before_after=fr, s4=s4r, parity_original=par, cos_base_total_fixed_vs_original=cosb,
                    Bt_original=orig["Bt"], Bt_fixed=fx["Bt"], lib_targets=dict(original=orig["lib_n"], fixed=fx["lib_n"]),
                    agreement_spearman_fixed_vs_original=float(np.corrcoef(rank01(fx["agree"]), rank01(orig["agree"]))[0, 1]))
                log("KOLF fixed M1: cos(base_total fixed, original)", round(cosb, 5))
                del orig
            else:
                fx = kolf_m1(src_fix, FIX / "cache_kolffree", s4_fix, tg, ctrl, genes, gmap)
            L0, Od, V3, Bt, ag = fx["L0"], fx["Od"], fx["V3"], fx["Bt"], fx["agree"]
            pk.Od, pk.ovj, pk.V3, pk.Bt = Od, np.isfinite(Od), V3, Bt     # scorer: pack real moments + the fixed build's Od/V3/Bt
            del fx
        else:
            L0 = pk.on_panel(zi["base_total"].astype(np.float64)); Od, V3, Bt = sp["Od"], sp["V3"].astype(np.float64), float(sp["Bt"])
            TM._cache.clear()
            Fm = F.fam_means(tg, fin, 20); ag = agreement(Fm, fin, len(tg)); del Fm; TM._cache.clear()
            Bt_re, gk, cr = budget(ctrl, L0, Od, n)
            rec["Bt_recomputed_rel"] = float(abs(Bt_re - Bt) / Bt); assert rec["Bt_recomputed_rel"] < 1e-5, rec["Bt_recomputed_rel"]
        rk = rank01(ag); amp_m3 = rk ** 0.75; wdf = np.where(amp_m3 <= np.quantile(amp_m3, 0.5), 0.75, 0.0)
        cpm = cpm_of(ctrl, ok)
        Vd, n_off, n_dup = ctrl_pcs(ctrl, gmap, ok)
        rec["ctrl_pc_genes_official"] = n_off; rec["ctrl_pc_dup_panel_genes_dropped"] = n_dup
        Xp, tinfo = sum_side(L0, Od, V3, cpm, rk); Xm, _ = sum_side(L0, Od, V3, cpm, rk, Vd=Vd, wdf=wdf)
        sp_, sm_ = bmatch(fm_full, Xp, Od, n, Bt), bmatch(fm_full, Xm, Od, n, Bt)
        Lp, Lm = Xp * sp_, Xm * sm_
        ep, em = pk.evaluate(Lp), pk.evaluate(Lm)
        # reference arm = the historical harness/ERA setting (S16/S18 M2 = rho_hat^0.75 post amplitude, no ew, ERA ctrlpc V[:8])
        rho = np.load(PACKS / f"aux_{p}.npz")["rho_hat"]; rr = rank01(rho)
        Vera = np.load(PACKS / f"ctrlpc_{p}.npz")["V"][:8].astype(np.float64)
        wref = np.where(rr ** 0.75 <= np.quantile(rr ** 0.75, 0.5), 0.75, 0.0)
        Xr, _ = sum_side(L0, Od, V3, cpm, rr, beta=0.75, ew=None); Xr3 = Xr - wref[:, None] * ((Xr @ Vera.T) @ Vera)
        er = pk.evaluate(Xr * bmatch(fm_full, Xr, Od, n, Bt)); er3 = pk.evaluate(Xr3 * bmatch(fm_full, Xr3, Od, n, Bt))
        # F2 M2 + ERA panel-axis PCs (isolates the PC-estimation difference)
        Xe, _ = sum_side(L0, Od, V3, cpm, rk, Vd=Vera, wdf=wdf); ee = pk.evaluate(Xe * bmatch(fm_full, Xe, Od, n, Bt))
        rec["sanity"] = dict(prem3=dict(pds=ep["pds"], mse=ep["mse"]), m3=dict(pds=em["pds"], mse=em["mse"]),
                             dPDS_m3_minus_prem3=em["pds"] - ep["pds"], dMSE_m3_minus_prem3=em["mse"] - ep["mse"],
                             ref_S18_M2_ERA_PCs=dict(prem3_pds=er["pds"], m3_pds=er3["pds"], dPDS=er3["pds"] - er["pds"]),
                             F2_M2_ERA_PCs=dict(m3_pds=ee["pds"], dPDS=ee["pds"] - ep["pds"]),
                             budget_scale=dict(prem3=sp_, m3=sm_), **tinfo,
                             spearman_M2_F2_vs_rho_hat=float(np.corrcoef(rk, rr)[0, 1]))
        if nofilm:                                             # sensitivity arm: S4 = refit only (FiLM residual weight 0)
            res = dict(arm="KOLF S4 FiLM-residual weight 0 (S4 = refit without pooled Feng rows); everything else = main KOLF build")
            for nm, Lx, ex_ in (("prem3", Lp, ep), ("m3", Lm, em)):
                zm = np.load(OUT / f"KOLF_{nm}.npz"); Lmain = zm["lfc"].astype(np.float64); dm = zm["delta"].astype(np.float64)
                dx = delta_lcp(Lx, Od, fm_full).astype(np.float64)
                res[nm] = dict(pds_noFiLM=ex_["pds"], mse_noFiLM=ex_["mse"], pds_main=float(zm["pds_offline"]), mse_main=float(zm["mse_offline"]),
                               cos_lfc=float((Lx * Lmain).sum() / np.linalg.norm(Lx) / np.linalg.norm(Lmain)),
                               relL2_lfc=float(np.linalg.norm(Lx - Lmain) / np.linalg.norm(Lmain)),
                               cos_delta=float((dx * dm).sum() / np.linalg.norm(dx) / np.linalg.norm(dm)),
                               relL2_delta=float(np.linalg.norm(dx - dm) / np.linalg.norm(dm)))
                meta = dict(json.loads(str(zm["meta"])), variant=res["arm"])
                np.savez_compressed(OUT / f"KOLF_{nm}_noFiLM.npz", delta=dx.astype(np.float32), lfc=Lx.astype(np.float32),
                                    budget_scale=(sp_ if nm == "prem3" else sm_), pds_offline=ex_["pds"], mse_offline=ex_["mse"],
                                    targets=np.array(tg), genes=genes, gmap=gmap, fm=fm_full, n_cells=n, Od=Od.astype(np.float32),
                                    agreement=ag, m2_rank=rk, m2_amp=rk ** 0.5, m3_weak=wdf > 0, meta=np.array(json.dumps(meta, default=float)))
            res["dPDS_m3_minus_prem3_noFiLM"] = em["pds"] - ep["pds"]; res["seconds"] = round(time.time() - tp, 1)
            jdump(OUT / "KOLF_noFiLM.json", res); log("NOFILM_DONE", res); return
        if p in ("RN_HepG2",):                                 # evaluator parity with sch_sumeval.evaluate (real-moments cache)
            try:
                import sch_sumeval as SE
                o = SE.evaluate(p, tg, list(genes), ctrl, Lm, Od, n)
                rec["sanity"]["sch_sumeval_parity"] = dict(pds=o["pds"], mse=o["mse"], pack_pds=round(em["pds"], 4), pack_mse=round(em["mse"], 4))
            except Exception as e:  # noqa
                rec["sanity"]["sch_sumeval_parity"] = f"not run: {e}"
        log(p, "sanity", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec["sanity"].items() if not isinstance(v, dict)},
            "dPDS", round(em["pds"] - ep["pds"], 5), "ref", round(er3["pds"] - er["pds"], 5))
        meta = dict(recipe=RECIPE, panel=p, gene_axis="truth var_names (panel axis); gmap = official index (-1 none)",
                    truth_path=str(CM.truth_path(p)), sources_used=fin, sources_excluded=json.dumps(
                        json.loads((ROOT / "output/data-analysis/results/panels.json").read_text())["panels"][p]["excluded_sources"]),
                    M2="F2 --sum-tamp 0.5: rank(cross-family agreement)^0.5 (tv2_tamp_prep.agreement over the panel's fusion families); "
                       "historical harness M2 (S16/S18 learned tlrn rho_hat beta 0.75) not used",
                    kolf_fix=json.dumps(rec.get("kolf_fix", "n/a (not the KOLF destination)"), default=float))
        common_arrays = dict(targets=np.array(tg), genes=genes, gmap=gmap, fm=fm_full, n_cells=n, Od=Od.astype(np.float32),
                             agreement=ag, m2_rank=rk, m2_amp=rk ** 0.5, m3_weak=wdf > 0,
                             meta=np.array(json.dumps(meta, default=float)))
        np.savez_compressed(OUT / f"{p}_prem3.npz", delta=delta_lcp(Lp, Od, fm_full), lfc=Lp.astype(np.float32),
                            budget_scale=sp_, pds_offline=ep["pds"], mse_offline=ep["mse"], arm=np.array("pre-M3"), **common_arrays)
        np.savez_compressed(OUT / f"{p}_m3.npz", delta=delta_lcp(Lm, Od, fm_full), lfc=Lm.astype(np.float32), ctrl_V=Vd.astype(np.float32),
                            budget_scale=sm_, pds_offline=em["pds"], mse_offline=em["mse"], arm=np.array("M3"), **common_arrays)

        # ------------------------------------------------------------ control dependence (item 4) + C_a u C_b variants
        nC = ctrl.shape[0]; dep = {}
        for k in splits:
            pc_i, a_i, b_i = CM.split_halves(nC, k, 3)
            lab = CM.labels(p) if k == splits[0] else lab
            cobs = np.flatnonzero(lab == "non-targeting"); assert len(cobs) == nC
            np.savez(OUT / f"{p}_ctrl_split_seed{k}.npz", C_pc=pc_i, C_a=a_i, C_b=b_i, obs_C_pc=cobs[pc_i], obs_C_a=cobs[a_i],
                     obs_C_b=cobs[b_i], n_ctrl=nC, seed=k, truth_path=str(CM.truth_path(p)),
                     rule="common.split_halves(n_ctrl, seed, 3) -> (C_pc, C_a, C_b); indices = positions among control cells in file order")
            ab = np.sort(np.concatenate([a_i, b_i])); cab = ctrl[ab]
            fm_ab = np.asarray(cab.sum(0)).ravel() / cab.shape[0]; cpm_ab = cpm_of(cab, ok)
            Bt_ab, gk_ab, cr_ab = budget(cab, L0, Od, n)
            Xab, _ = sum_side(L0, Od, V3, cpm_ab, rk); s_ab = bmatch(fm_ab, Xab, Od, n, Bt_ab); Lab = Xab * s_ab
            if k == splits[0]:                                   # decomposition: which control quantity moves the prediction
                s_prof = bmatch(fm_ab, sum_side(L0, Od, V3, cpm_ab, rk)[0], Od, n, Bt)       # profile from C_ab, budget Bt full
                Xpf, _ = sum_side(L0, Od, V3, cpm, rk); s_bud = bmatch(fm_full, Xpf, Od, n, Bt_ab)   # profile full, budget C_ab
                nrm = np.linalg.norm(Lp)
                dep = dict(
                    uses=["fm = control mean counts (baseline profile): pseudobulk baseline + budget match",
                          "cpm of controls (baseline profile): ew expression weights",
                          "C_real(n) = delete-1 jackknife trace of n random control cells (rng 0, 4 draws) inside kink_gamma -> the "
                          "budget Bt -> one GLOBAL scale of the prediction: a cell-level (variance) quantity, i.e. BEYOND the profile",
                          "M2 ranks, fusion, S4, neighbour blend, library PCs: source data only (S4's FiLM context = the line's control "
                          "profile from the source-side line data, a profile quantity)"],
                    seed0_rel_change_all=float(np.linalg.norm(Lab - Lp) / nrm),
                    seed0_rel_change_profile_only=float(np.linalg.norm(sum_side(L0, Od, V3, cpm_ab, rk)[0] * s_prof - Lp) / nrm),
                    seed0_rel_change_budget_only=float(np.linalg.norm(Xpf * s_bud - Lp) / nrm),
                    seed0_budget_scale=dict(full=sp_, ctrlAB=s_ab, budget_only=s_bud), seed0_c_real=dict(ctrlAB=cr_ab),
                    seed0_Bt_rel=float(Bt_ab / Bt - 1))
            np.savez_compressed(OUT / f"{p}_prem3_ctrlAB_seed{k}.npz", delta=delta_lcp(Lab, Od, fm_ab), lfc=Lab.astype(np.float32),
                                fm=fm_ab, budget_scale=s_ab, Bt=Bt_ab, kink=gk_ab, c_real=cr_ab, ctrl_rows=ab, split_seed=k,
                                targets=np.array(tg), genes=genes, gmap=gmap, m2_rank=rk, Od=Od.astype(np.float32),
                                meta=np.array(json.dumps(dict(meta, variant="pre-M3 with every control-derived quantity (fm, ew cpm, "
                                                              "kink C_real / budget Bt, budget match) from C_a u C_b only"), default=float)))
        rec["control_dependence"] = dep; rec["seconds"] = round(time.time() - tp, 1)
        summary[p] = rec; jdump(OUT / "rq3_build_preds_summary.json", summary)
        log(p, "dependence", {k: v for k, v in dep.items() if k != "uses"}, "secs", rec["seconds"])
        del zi, ctrl, L0
    allp = [q for q in PANELS7 if q in summary]
    d = {q: summary[q]["sanity"]["dPDS_m3_minus_prem3"] for q in allp}
    summary["_aggregate"] = dict(dPDS=d, n_positive=int(sum(v > 0 for v in d.values())), n=len(d),
                                 ref_S18_dPDS={q: summary[q]["sanity"]["ref_S18_M2_ERA_PCs"]["dPDS"] for q in allp})
    jdump(OUT / "rq3_build_preds_summary.json", summary); log("BUILD_DONE", summary["_aggregate"])


if __name__ == "__main__":
    cmd = sys.argv[1]; args = sys.argv[2:]
    if cmd == "feng":
        cmd_feng()
    elif cmd == "s4kolf":
        cmd_s4kolf()
    elif cmd == "kolf_nofilm":
        build_main(["KOLF"], nofilm=True)
    elif cmd == "build":
        build_main(args or PANELS7)
    else:
        raise SystemExit(f"unknown subcommand {cmd}")
