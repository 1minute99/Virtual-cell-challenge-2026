"""Task 17 (plan `input/ideas/zero-shot-limits-context-complexity-plan.md`, RQ3, E1b / H3b): M3 (destination control-PC
deflation) vs five alternative corrections, inserted at M3's position in the frozen sum-side chain (after M2, before the
budget match), on the 7 harness panels; offline per-target PDS + within-panel target bootstrap; H1 real-scorer runs;
leave-one-panel-out (LOPO) choice of (k, eta) for M3. Inputs: Task 14 `results/preds/<panel>_{prem3,m3}.npz`.

Subcommands (vcc2026 env; nice -n 10; OMP/MKL/OpenBLAS 4 threads):
  offline            all arms x 7 panels -> results/rq3_corr/offline_*.{npz,json}, per-target table, H1 arm lfcs
  real ARM [ARM..]   H1 real scorer (`cell-eval2 --preset vcc2026 --set de.backend=pdex`), one arm after the other
  verdict            bootstrap, budget-invariance check, LOPO, H3b verdict -> results/rq3_h3b.json, rq3_corrections.parquet

ARMS (all act only on the weaker half of targets by M2 rank -- Task 14's `m3_weak`, i.e. rank(agreement) <= median --
  with k = 8, eta = 0.75, at the M3 position; every arm is then budget-matched to the same base budget Bt):
  none        X (the pre-M3 sum-side matrix after M2) unchanged
  M3_dest     X_w <- X_w - eta (X_w V^T) V, V = Task 14's frozen destination control PCs (`ctrl_V`, 8 x G)
  M3_other_j  same with V = top-8 PCs of other harness line j's controls (6 lines; the builder recipe -- 6,000 random
              control cells rng 0, log1p CP10K on the official axis, gene-centred, randomized_svd n_iter 5 rs 0 -- then
              mapped to the destination's gene axis and re-orthonormalised); primary = mean of the 6 lines' per-target PDS
  M3_perm_i   same with V = top-8 PCs of the destination's 6,000-cell control matrix after an independent per-gene
              permutation across cells (i = 0..19, rng default_rng([i, 11])); primary = mean over the 20 permutations
  predPC      same with V = top-8 right singular vectors of X itself (override entries zeroed, centred across targets)
  sys_eta     Systema-style: X_w <- X_w - eta * mean_t X (mean over all targets, override entries excluded, as the
              chain's own centring), eta in {0.25, 0.5, 0.75, 1.0}; the best eta per panel (panel-mean PDS) is kept --
              this selection on the scored data favours the alternative, as the plan specifies
  LOPO grid   M3_dest with k in {8, 24} (k 24: the same recipe with 24 components), eta in {0.5, 0.75, 1.0}

ASSUMPTIONS (stated; checked in the code; see knowledge/synthesis/rq3_h3b-assumptions.md)
  A1  X = lfc_prem3 / budget_scale_prem3 (Task 14 stored the budget-matched pre-M3 lfc and its scale; float32 storage)
      -- parity: re-matching `none` reproduces Task 14's pre-M3 offline PDS and M3_dest reproduces its M3 PDS (|d| < 2e-4).
  A2  Offline scorer = the ERA-pack evaluator (`tasks/vcc_sum/evaluator.Pack`: real moments from the pack, pinned
      sums, midrank PDS) = `sch_sumeval.evaluate` numerically (Task 14 parity on RN_HepG2). Exact for budget-INVARIANT arms
      (AGENTS.md 2026-10-03 rule); every arm here is budget-matched to the same Bt, and step 1 checks it on H1 against
      the real scorer (threshold |PDS_real - PDS_offline| <= 0.001). KOLF: Od / Bt from the kolf-fixed Task 14 build.
  A3  Per-target PDS_t (1 - midrank/T) is the bootstrap unit; targets are exchangeable within a panel; panels are
      weighted equally; B = 2,000 percentile CIs; the same resampled target sets are used for every alternative
      (paired). The CI is conditional on the fixed prediction and truth (it does not include truth-cell resampling).
  A4  Real-scorer emission (H1): per target, the integer pinned sums T_g = round(fm_g * exp(L_tg) * n) (override entries
      exp(Od)) -- the very sums the offline evaluator uses -- are distributed over n = 300 cells drawn without
      replacement from the destination controls (rng 0) in proportion to each cell's own counts (largest-remainder
      rounding; exact column sums). PDS / MSE depend on the sums only (AGENTS.md: totals fix the pseudobulk). The DE side
      of this simple emission is not the frozen F2 emission, so DE metrics and the composite are recorded but are NOT
      comparable with the official pipeline; only PDS is used.
  A5  M3_perm and M3_other enter the real-scorer runs as a single matrix: permutation 0 (pre-specified, non-selective)
      and the M3_other line with the highest offline H1 panel PDS (the plan's "best M3_other line").
  A6  The H1 real-moments cache is missing (Task 4); the first real-scorer run writes a fresh one under
      results/rq3_cache/cache_real_H1 (protected report dirs are not touched); `sch_sumeval.evaluate` is then also run on
      that cache to confirm the pack moments.
Seeds: control-cell draw rng 0 and randomized_svd random_state 0 (builder), permutations default_rng([i, 11]) i = 0..19,
emission rng 0, bootstrap default_rng(20261105)."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import json, pathlib, shutil, subprocess, sys, time
import numpy as np

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); OUTD = ROOT / "output/data-analysis"; RES = OUTD / "results"
PRED = RES / "preds"; CW = RES / "rq3_corr"; CACHE = RES / "rq3_cache"; PACKS = ROOT / "vcc2026_data/reports/era_2026-10-07/packs"
CE = pathlib.Path("/home/wonmin-kim/miniconda3/envs/vcc2026/bin/cell-eval2")
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
K, ETA = 8, 0.75
SYS_ETAS = (0.25, 0.5, 0.75, 1.0)
N_PERM = 20; B_BOOT = 2000; BOOT_SEED = 20261105
LOPO_K = (8, 24); LOPO_ETA = (0.5, 0.75, 1.0)
REAL_TOL = 0.001
ALTS = ["none", "M3_other", "M3_perm", "predPC", "sys_best"]
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def jdump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)))


def imports():
    from cell_eval2.moments import jackknife_correction  # noqa: F401  (import order rule: before scripts/ joins sys.path)
    sys.path.insert(0, str(OUTD)); sys.path.insert(0, str(ROOT / "external/era-aivc/aivc"))
    import common as CM
    from tasks.vcc_sum import evaluator as EV
    return CM, EV


# ------------------------------------------------------------------------------------------------ subspaces
def ctrl_official(ctrl, gmap, NG, ncells=6000, seed=0):
    """Builder sum.ctrl_defl control matrix: log1p CP10K on the official axis, gene-centred (rq3_build_preds.ctrl_pcs)."""
    from scipy import sparse
    ok = gmap >= 0; okc = np.flatnonzero(ok); g = gmap[okc]; first = np.unique(g, return_index=True)[1]; okc, g = okc[first], g[first]
    M = sparse.csr_matrix(ctrl, dtype=np.float64)[:, okc].tocoo()
    Mo = sparse.csr_matrix((M.data, (M.row, g[M.col])), shape=(ctrl.shape[0], NG))
    L_ = np.asarray(Mo.sum(1)).ravel(); Xc = Mo.multiply(1e4 / np.maximum(L_, 1)[:, None]).tocsr(); Xc.data = np.log1p(Xc.data)
    ix = np.random.default_rng(seed).choice(Xc.shape[0], min(ncells, Xc.shape[0]), replace=False)
    Dc = np.asarray(Xc[ix].todense()); Dc -= Dc.mean(0)
    return Dc


def top_v(D, k):
    from sklearn.utils.extmath import randomized_svd
    return randomized_svd(D, n_components=k, n_iter=5, random_state=0)[2].astype(np.float64)


def to_panel(Vo, gmap):
    ok = gmap >= 0; Vp = np.zeros((Vo.shape[0], len(gmap))); Vp[:, ok] = Vo[:, gmap[ok]]; return Vp


def orthonormal_rows(V):
    return np.linalg.svd(V, full_matrices=False)[2]


def official_subspaces(CM):
    """Per panel (cached): Vo8, Vo24 and the 20 permuted Vo8 on the official axis."""
    NG = len(CM.official_genes()); out = {}
    CACHE.mkdir(parents=True, exist_ok=True)
    for p in PANELS7:
        f = CACHE / f"rq3corr_Vofficial_{p}.npz"
        if f.exists():
            z = np.load(f); out[p] = dict(V8=z["V8"].astype(np.float64), V24=z["V24"].astype(np.float64),
                                          Vperm=z["Vperm"].astype(np.float64)); continue
        tp = time.time(); ctrl, genes = CM.controls(p); gmap = CM.to_official(np.asarray(genes))
        Dc = ctrl_official(ctrl, gmap, NG); del ctrl
        V8 = top_v(Dc, 8); V24 = top_v(Dc, 24)
        Vperm = np.stack([top_v(np.random.default_rng([i, 11]).permuted(Dc, axis=0), 8) for i in range(N_PERM)])
        del Dc
        np.savez(f, V8=V8.astype(np.float32), V24=V24.astype(np.float32), Vperm=Vperm.astype(np.float32))
        out[p] = dict(V8=V8, V24=V24, Vperm=Vperm); log(p, "official subspaces", round(time.time() - tp, 1), "s")
    return out


def defl(X, w, V):
    return X - w[:, None] * ((X @ V.T) @ V)


# ------------------------------------------------------------------------------------------------ offline
def cmd_offline():
    CM, EV = imports()
    from sklearn.utils.extmath import randomized_svd
    summ = json.loads((PRED / "rq3_build_preds_summary.json").read_text())
    VO = official_subspaces(CM)
    CW.mkdir(parents=True, exist_ok=True)
    rows = []; panel_tab = {}; checks = {}; h1_real = {}
    for p in PANELS7:
        tp = time.time()
        z = np.load(PRED / f"{p}_prem3.npz", allow_pickle=True); zm = np.load(PRED / f"{p}_m3.npz", allow_pickle=True)
        assert list(z["targets"]) == list(zm["targets"])
        tg = [str(t) for t in z["targets"]]; gmap = z["gmap"]; T = len(tg)
        X0 = z["lfc"].astype(np.float64) / float(z["budget_scale"])
        Od = z["Od"].astype(np.float64); ovj = np.isfinite(Od)
        weak = zm["m3_weak"].astype(bool); w = np.where(weak, ETA, 0.0)
        pk = EV.Pack(PACKS / f"score_{p}.npz")
        if p == "KOLF":
            pk.Od, pk.ovj, pk.Bt = Od, ovj, float(summ["KOLF"]["kolf_fix"]["Bt_fixed"])
        else:
            assert np.array_equal(np.isfinite(pk.Od), ovj)
        assert np.allclose(pk.fm, z["fm"], rtol=1e-5)

        def run(X):
            s = pk.budget_match(X); e = pk.evaluate(X * s)
            return dict(pds=e["pds"], mse=e["mse"], pds_t=e["pds_t"], scale=s)

        res = {}
        res["none"] = run(X0)
        Vd = zm["ctrl_V"].astype(np.float64)
        res["M3_dest"] = run(defl(X0, w, Vd))
        # parity / subspace checks
        Vd_re = to_panel(VO[p]["V8"], gmap)
        ck = dict(none_pds_vs_T14=res["none"]["pds"] - float(z["pds_offline"]), m3_pds_vs_T14=res["M3_dest"]["pds"] - float(zm["pds_offline"]),
                  ctrlV_recompute_absdot_min=float(np.min(np.abs((Vd_re * Vd).sum(1)) / np.linalg.norm(Vd_re, axis=1) / np.linalg.norm(Vd, axis=1))),
                  Vd_gram_maxdev=float(np.abs(Vd @ Vd.T - np.eye(8)).max()), n_weak=int(weak.sum()), T=T)
        assert abs(ck["none_pds_vs_T14"]) < 2e-4 and abs(ck["m3_pds_vs_T14"]) < 2e-4, ck
        assert ck["ctrlV_recompute_absdot_min"] > 0.99, ck
        # M3_other (6 lines)
        cov = {}
        for j in PANELS7:
            if j == p:
                continue
            Vj = to_panel(VO[j]["V8"], gmap); cov[j] = float((Vj ** 2).sum() / 8)
            res[f"M3_other_{j}"] = run(defl(X0, w, orthonormal_rows(Vj)))
        ck["M3_other_energy_on_dest_axis"] = cov
        # M3_perm (20)
        for i in range(N_PERM):
            res[f"M3_perm_{i}"] = run(defl(X0, w, to_panel(VO[p]["Vperm"][i], gmap)))
        # predPC
        Xz = np.where(ovj, 0.0, X0); Xz = Xz - Xz.mean(0)
        Vpred = randomized_svd(Xz, n_components=8, n_iter=5, random_state=0)[2].astype(np.float64); del Xz
        res["predPC"] = run(defl(X0, w, Vpred))
        # Systema-style
        mt = np.nan_to_num(np.nanmean(np.where(ovj, np.nan, X0), 0))
        for e_ in SYS_ETAS:
            res[f"sys_{e_:g}"] = run(X0 - np.where(weak, e_, 0.0)[:, None] * mt[None, :])
        # LOPO grid (M3_dest, k in {8, 24})
        V24 = to_panel(VO[p]["V24"], gmap)
        for k_ in LOPO_K:
            for e_ in LOPO_ETA:
                Vk = Vd if k_ == 8 else V24
                key = f"lopo_k{k_}_e{e_:g}"
                res[key] = res["M3_dest"] if (k_ == 8 and e_ == ETA) else run(defl(X0, np.where(weak, e_, 0.0), Vk))
        # composites
        other_keys = [k for k in res if k.startswith("M3_other_")]; perm_keys = [k for k in res if k.startswith("M3_perm_")]
        sys_keys = [f"sys_{e_:g}" for e_ in SYS_ETAS]; best_sys = max(sys_keys, key=lambda k: res[k]["pds"])
        comp = {"M3_other": np.mean([res[k]["pds_t"] for k in other_keys], 0),
                "M3_perm": np.mean([res[k]["pds_t"] for k in perm_keys], 0),
                "sys_best": res[best_sys]["pds_t"]}
        panel_tab[p] = dict(T=T, sys_best_eta=best_sys, **{k: dict(pds=v["pds"], mse=v["mse"], scale=v["scale"]) for k, v in res.items()},
                            **{f"{k}_pds": float(v.mean()) for k, v in comp.items()})
        checks[p] = ck
        for k, v in res.items():
            grp = ("M3_other" if k.startswith("M3_other_") else "M3_perm" if k.startswith("M3_perm_") else
                   "systema" if k.startswith("sys_") else "lopo" if k.startswith("lopo_") else k)
            for t in range(T):
                rows.append((p, k, grp, tg[t], bool(weak[t]), float(v["pds_t"][t]), v["pds"], v["mse"], v["scale"]))
        for k, v in comp.items():
            for t in range(T):
                rows.append((p, k, "composite", tg[t], bool(weak[t]), float(v[t]), float(v.mean()), np.nan, np.nan))
        if p == "H1":                                  # matrices for the real-scorer runs
            best_other = max(other_keys, key=lambda k: res[k]["pds"])
            sel = dict(none=X0, M3_dest=defl(X0, w, Vd),
                       M3_other_best=defl(X0, w, orthonormal_rows(to_panel(VO[best_other.split("M3_other_")[1]]["V8"], gmap))),
                       M3_perm0=defl(X0, w, to_panel(VO[p]["Vperm"][0], gmap)), predPC=defl(X0, w, Vpred),
                       sys_best=X0 - np.where(weak, float(best_sys.split("_")[1]), 0.0)[:, None] * mt[None, :])
            src = dict(none="none", M3_dest="M3_dest", M3_other_best=best_other, M3_perm0="M3_perm_0", predPC="predPC", sys_best=best_sys)
            arr = {}
            for a, X in sel.items():
                rr = run(X); assert abs(rr["pds"] - res[src[a]]["pds"]) < 1e-9, a
                arr[f"lfc_{a}"] = (X * rr["scale"]).astype(np.float64); h1_real[a] = dict(source_arm=src[a], pds_offline=rr["pds"],
                                                                                        mse_offline=rr["mse"], scale=rr["scale"])
            np.savez(CW / "H1_real_arms.npz", Od=Od, targets=np.array(tg), genes=z["genes"].astype(str), fm=z["fm"],
                     n_cells=int(z["n_cells"]), **arr, **{f"pdst_{a}": res[src[a]]["pds_t"] for a in sel})
        log(p, "done", round(time.time() - tp, 1), "s | none", round(res["none"]["pds"], 4), "M3", round(res["M3_dest"]["pds"], 4),
            "other", round(float(comp["M3_other"].mean()), 4), "perm", round(float(comp["M3_perm"].mean()), 4),
            "predPC", round(res["predPC"]["pds"], 4), best_sys, round(res[best_sys]["pds"], 4))
    import pandas as pd
    df = pd.DataFrame(rows, columns=["panel", "arm", "arm_group", "target", "weak_half", "pds_t", "panel_pds", "panel_mse", "budget_scale"])
    df.to_parquet(CW / "offline_per_target.parquet", index=False)
    jdump(CW / "offline_summary.json", dict(panels=panel_tab, checks=checks, h1_real_arms=h1_real))
    log("OFFLINE_DONE", len(df), "rows")


# ------------------------------------------------------------------------------------------------ real scorer (H1)
def largest_remainder(W, T):
    fl = np.floor(W); rem = (T - fl.sum(0)).astype(int); frac = W - fl; cols = np.flatnonzero(rem > 0)
    if cols.size:
        oo = np.argsort(-frac[:, cols], axis=0); mask = np.arange(W.shape[0])[:, None] < rem[cols][None, :]
        add = np.zeros((W.shape[0], cols.size)); np.put_along_axis(add, oo, mask.astype(float), axis=0); fl[:, cols] += add
    return fl


def emit_share(ctrl, L, Od, n, seed=0):
    from scipy import sparse
    rng = np.random.default_rng(seed); ctrl = ctrl.tocsr(); N = ctrl.shape[0]
    fm = np.asarray(ctrl.sum(0)).ravel().astype(np.float64) / N; out = []; Ts = []
    for t in range(L.shape[0]):
        blk = ctrl[np.sort(rng.choice(N, n, replace=False))].toarray().astype(np.float64)
        e = np.exp(L[t]); j = np.isfinite(Od[t]); e[j] = np.exp(Od[t, j]); T = np.round(fm * e * n)
        cs = blk.sum(0); share = np.where(cs[None, :] > 0, blk / np.maximum(cs, 1e-12)[None, :], 1.0 / n)
        fl = largest_remainder(share * T[None, :], T); assert np.array_equal(fl.sum(0), T)
        out.append(sparse.csr_matrix(fl.astype(np.float32))); Ts.append(T)
    return out, np.array(Ts)


def cmd_real(arms):
    CM, EV = imports()
    import anndata as ad, pandas as pd
    from scipy import sparse
    z = np.load(CW / "H1_real_arms.npz", allow_pickle=True)
    tg = [str(t) for t in z["targets"]]; n = int(z["n_cells"]); Od = z["Od"]
    ctrl, genes = CM.controls("H1"); genes = np.asarray(genes).astype(str)
    assert list(genes) == list(z["genes"]); truth = CM.truth_path("H1")
    lab = set(CM.labels("H1")) - {"non-targeting"}; assert set(tg) == lab, "target set != truth labels"
    fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; assert np.allclose(fm, z["fm"], rtol=1e-5)
    cache = CACHE / "cache_real_H1"; cache.mkdir(parents=True, exist_ok=True)
    for a in arms:
        out = CW / "real" / a
        if (out / "result.json").exists():
            log(a, "cached"); continue
        out.mkdir(parents=True, exist_ok=True); ta = time.time()
        blocks, Ts = emit_share(ctrl, z[f"lfc_{a}"], Od, n)
        X = sparse.vstack(blocks, format="csr"); del blocks
        pred = ad.AnnData(X=sparse.vstack([X, ctrl], format="csr"),
                          obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(tg, n), *["non-targeting"] * ctrl.shape[0]])},
                                           index=[f"p{i}" for i in range(X.shape[0] + ctrl.shape[0])]),
                          var=pd.DataFrame(index=genes))
        pp = out / "pred.h5ad"; pred.write_h5ad(pp); del pred, X
        log(a, "emitted", round(time.time() - ta, 1), "s; scoring")
        cmd = [str(CE), "run", "-ap", str(pp), "-ar", str(truth), "--preset", "vcc2026", "--pert-col", "target",
               "--control", "non-targeting", "--input-type", "counts", "-o", str(out / "cell_eval2"), "--cache-real", str(cache),
               "--set", "de.backend=pdex", "--set", "num_threads=4", "--set", "gather_threads=4"]
        env = dict(os.environ, NUMBA_NUM_THREADS="4", OMP_NUM_THREADS="4")
        with (out / "scorer.log").open("w") as g:
            subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
        agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
        raw = {k: float(agg.loc["mean", k]) for k in agg.columns if k in agg.columns and np.isfinite(pd.to_numeric(agg.loc["mean", k], errors="coerce"))}
        rl = pd.read_csv(out / "cell_eval2/results.csv")
        pt = rl[rl["metric"] == "pds_cosine"].set_index("perturbation")["value"]
        pt.to_csv(out / "pds_per_target.csv")
        shutil.copy(out / "cell_eval2/agg_results.csv", out / "agg_results.csv")
        res = dict(arm=a, raw=raw, pds_real=raw.get("pds_cosine"), n_targets_scored=int(len(pt)), seconds=round(time.time() - ta, 1),
                   emission="share of n=300 random controls (rng 0), pinned sums round(fm exp(L) n); DE side NOT the F2 emission")
        jdump(out / "result.json", res); pp.unlink(); shutil.rmtree(out / "cell_eval2", ignore_errors=True)
        log(a, "REAL pds", res["pds_real"], "secs", res["seconds"])
    log("REAL_DONE")


# ------------------------------------------------------------------------------------------------ verdict
def boot_delta(df, base, alt, panels):
    """7-panel mean of (base - alt) per-target PDS; per-panel paired target bootstrap."""
    rng = np.random.default_rng(BOOT_SEED); pm = []; bm = []
    for p in panels:
        a = df[(df.panel == p) & (df.arm == base)].set_index("target")["pds_t"]
        b = df[(df.panel == p) & (df.arm == alt)].set_index("target")["pds_t"].reindex(a.index)
        d = (a - b).to_numpy(); T = len(d); idx = rng.integers(0, T, (B_BOOT, T))
        pm.append(d.mean()); bm.append(d[idx].mean(1))
    pm = np.array(pm); bm = np.mean(bm, 0)
    return dict(mean=float(pm.mean()), ci95=[float(np.quantile(bm, 0.025)), float(np.quantile(bm, 0.975))],
                per_panel={p: float(v) for p, v in zip(panels, pm)}, n_panels_positive=int((pm > 0).sum()))


def cmd_verdict():
    import pandas as pd
    df = pd.read_parquet(CW / "offline_per_target.parquet"); off = json.loads((CW / "offline_summary.json").read_text())
    P = off["panels"]
    # sys_best rows -> arm "sys_best" exists as a composite; the plan's five alternatives
    alt_arm = dict(none="none", M3_other="M3_other", M3_perm="M3_perm", predPC="predPC", sys_best="sys_best")
    deltas = {k: boot_delta(df, "M3_dest", v, PANELS7) for k, v in alt_arm.items()}
    # each individual other line, aggregated as "the line's PCs" is destination-dependent -> report per (dest, line) panel deltas
    per_line = {p: {k.split("M3_other_")[1]: P[p]["M3_dest"]["pds"] - v["pds"] for k, v in P[p].items()
                    if isinstance(v, dict) and k.startswith("M3_other_")} for p in PANELS7}
    worst_line = {p: min(v, key=v.get) for p, v in per_line.items()}
    # LOPO
    grid = [f"lopo_k{k}_e{e:g}" for k in LOPO_K for e in LOPO_ETA]
    lopo = {}
    for p in PANELS7:
        tr = [q for q in PANELS7 if q != p]
        sc = {g: float(np.mean([P[q][g]["pds"] - P[q]["none"]["pds"] for q in tr])) for g in grid}
        ch = max(sc, key=sc.get)
        lopo[p] = dict(chosen=ch, train_gain=sc[ch], heldout_gain_vs_none=P[p][ch]["pds"] - P[p]["none"]["pds"],
                       heldout_gain_frozen_vs_none=P[p]["M3_dest"]["pds"] - P[p]["none"]["pds"],
                       heldout_chosen_minus_frozen=P[p][ch]["pds"] - P[p]["M3_dest"]["pds"])
    lopo_sum = dict(mean_heldout_gain_vs_none=float(np.mean([v["heldout_gain_vs_none"] for v in lopo.values()])),
                    mean_frozen_gain_vs_none=float(np.mean([v["heldout_gain_frozen_vs_none"] for v in lopo.values()])),
                    n_heldout_positive=int(sum(v["heldout_gain_vs_none"] > 0 for v in lopo.values())),
                    choices={p: v["chosen"] for p, v in lopo.items()})
    # real scorer (H1)
    real = {}; h1r = off["h1_real_arms"]
    for a in h1r:
        f = CW / "real" / a / "result.json"
        if f.exists():
            r = json.loads(f.read_text()); real[a] = dict(pds_real=r["pds_real"], pds_offline=h1r[a]["pds_offline"],
                                                          diff=r["pds_real"] - h1r[a]["pds_offline"], source_arm=h1r[a]["source_arm"])
            pt = pd.read_csv(CW / "real" / a / "pds_per_target.csv", index_col=0).iloc[:, 0]
            zz = np.load(CW / "H1_real_arms.npz", allow_pickle=True); tg = [str(t) for t in zz["targets"]]
            off_t = pd.Series(zz[f"pdst_{a}"], index=tg)
            real[a]["per_target_maxabs_real_minus_offline"] = float((pt.reindex(tg) - off_t).abs().max())
    inv = {a: dict(diff=real[a]["diff"], invariant=abs(real[a]["diff"]) <= REAL_TOL) for a in ("M3_dest", "sys_best") if a in real}
    real_alt = dict(none="none", M3_other="M3_other_best", M3_perm="M3_perm0", predPC="predPC", sys_best="sys_best")
    real_sign = {}
    if "M3_dest" in real:
        for k, a in real_alt.items():
            if a in real:
                real_sign[k] = dict(dPDS_real=real["M3_dest"]["pds_real"] - real[a]["pds_real"],
                                    dPDS_offline=real["M3_dest"]["pds_offline"] - real[a]["pds_offline"])
    # sch_sumeval on the freshly written H1 cache (A6)
    sse = {}
    try:
        from cell_eval2.moments import jackknife_correction  # noqa: F401
        CM, EV = imports(); sys.path.append(str(ROOT / "vcc2026_data/scripts")); import sch_sumeval as SE, glob
        cf = glob.glob(str(CACHE / "cache_real_H1/pseudobulk_moments_bulk_lognorm-*.moments.npz"))
        if cf:
            SE._REAL["H1"] = np.load(cf[0]); zz = np.load(CW / "H1_real_arms.npz", allow_pickle=True)
            ctrl, genes = CM.controls("H1")
            for a in h1r:
                o = SE.evaluate("H1", [str(t) for t in zz["targets"]], list(np.asarray(genes).astype(str)), ctrl, zz[f"lfc_{a}"], zz["Od"], int(zz["n_cells"]))
                sse[a] = dict(sch_sumeval_pds=o["pds"], pack_pds=round(h1r[a]["pds_offline"], 4))
    except Exception as e:  # noqa
        sse = dict(error=repr(e))
    # verdict
    neg = [k for k, v in deltas.items() if v["mean"] <= 0]
    all_ci = all(v["ci95"][0] > 0 for v in deltas.values())
    real_ok = (len(real_sign) == 5 and all(v["dPDS_real"] > 0 for v in real_sign.values()))
    inv_ok = all(v["invariant"] for v in inv.values()) and len(inv) == 2
    if neg:
        verdict = "rejected"; why = f"7-panel mean dPDS <= 0 against: {neg}"
    elif all_ci and real_ok:
        verdict = "confirmed"; why = "all five 7-panel mean dPDS > 0 with 95% CIs excluding 0 and H1 real-scorer dPDS > 0 against all five"
    else:
        verdict = "inconclusive"
        why = ("all five means > 0 but " + ", ".join(x for x, c in (("not every CI excludes 0", not all_ci),
                                                                   ("H1 real-scorer sign not positive against all five", not real_ok)) if c))
    out = dict(task="Task 17 / E1b / H3b", plan="input/ideas/zero-shot-limits-context-complexity-plan.md",
               rule="rejected if any 7-panel mean dPDS(M3 - alt) <= 0; confirmed if all five > 0 with 95% CIs excluding 0 AND the H1 "
                    "real-scorer run has the same sign against all five; otherwise inconclusive",
               verdict=verdict, verdict_reason=why, dPDS_7panel=deltas, budget_invariance_H1=inv, budget_invariance_ok=inv_ok,
               real_scorer_H1=real, real_scorer_H1_dPDS=real_sign, sch_sumeval_on_new_H1_cache=sse,
               panel_pds={p: dict(none=P[p]["none"]["pds"], M3_dest=P[p]["M3_dest"]["pds"], M3_other=P[p]["M3_other_pds"],
                                  M3_perm=P[p]["M3_perm_pds"], predPC=P[p]["predPC"]["pds"], sys_best=P[p]["sys_best_pds"],
                                  sys_best_eta=P[p]["sys_best_eta"]) for p in PANELS7},
               M3_other_per_line_dPDS=per_line, M3_other_worst_line=worst_line, lopo=lopo, lopo_summary=lopo_sum,
               checks=off["checks"], settings=dict(k=K, eta=ETA, sys_etas=SYS_ETAS, n_perm=N_PERM, B=B_BOOT, boot_seed=BOOT_SEED,
                                                   lopo_k=LOPO_K, lopo_eta=LOPO_ETA, real_tol=REAL_TOL))
    jdump(RES / "rq3_h3b.json", out)
    df.to_parquet(RES / "rq3_corrections.parquet", index=False)
    log("VERDICT", verdict, why)
    for k, v in deltas.items():
        log(f"  dPDS M3 - {k}: {v['mean']:+.5f} [{v['ci95'][0]:+.5f}, {v['ci95'][1]:+.5f}] pos {v['n_panels_positive']}/7")
    log("  real", real_sign); log("  inv", inv); log("  lopo", lopo_sum); log("  sse", sse)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "offline":
        cmd_offline()
    elif cmd == "real":
        cmd_real(sys.argv[2:])
    elif cmd == "verdict":
        cmd_verdict()
    else:
        raise SystemExit(f"unknown subcommand {cmd}")
