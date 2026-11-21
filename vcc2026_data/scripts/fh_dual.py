"""Dual-moment emission on the harness (LOCAL MEASUREMENT; assumes the organisers allow it — see
reports/submission_2026-09-23/ORGANIZER_INQUIRY_DRAFT.md).

Group sums (drive PDS and the pseudobulk MSE): for every target and gene the n-cell group total is pinned to
  round(full control mean x exp(gamma * d) x n), with the knockdown / promoter override applied multiplicatively,
  where d is the blend30_e delta (true log fold change scale of the recipe) and gamma a small amplitude.
Cell distribution (drives the Wilcoxon DE members): starting from the shared control block allocation of each
  pinned total, the top-K genes of the target by |d| among genes the context expresses (>= 5 CPM) are re-allocated
  WITHIN the column (sum unchanged): up-genes are flattened towards the uniform allocation, down-genes are
  concentrated into the top 10 % of cells; the strength lambda decays with the confidence rank (1 -> 0.3) so the
  p-value order follows |d| (reach). Every other gene keeps the proportional allocation.
Usage: fh_dual.py PANEL gamma:K[:boost[:dd]] ...   -> scores/dual_g{gamma}_k{K}[_b{boost}][_dd]/PANEL/result.json
  boost: extra amplitude multiplier on the K called genes only. Flags (comma-separated): dd = concentrate down-genes
  into the deepest cells; du<a> = flatten up-genes towards x ~ L^-a; ag = order re-allocation strength by source sign
  agreement x |d| (membership still top-K by |d|).
"""
import json, os, subprocess, sys, time
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2.moments import jackknife_correction  # noqa: F401
_argv = list(sys.argv); sys.argv = ["x"]
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H


def emit_dual(ctrl, Ld, Od, n, K, seed, gate, call_boost=1.0, depth_down=False, depth_up=0.0, rank_score=None):
    rng = np.random.default_rng(seed)
    ctrl = ctrl.tocsr()
    pick = rng.choice(ctrl.shape[0], size=n, replace=False)
    blk = np.asarray(ctrl[pick].todense(), dtype=np.float64)
    full_mean = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]
    colsum = blk.sum(0)
    share = np.where(colsum > 0, blk / np.maximum(colsum, 1e-12), 1.0 / n)
    lam_rank = np.linspace(1.0, 0.3, max(K, 1))
    n_top = max(1, int(round(0.1 * n)))
    order_cells = np.argsort(-blk, axis=0, kind="stable")          # per gene, cells by control count (for concentration)
    Lc = blk.sum(1)
    if depth_down:                                                  # concentrate down-genes into the deepest cells
        order_cells = np.repeat(np.argsort(-Lc, kind="stable")[:, None], blk.shape[1], 1)
    out = []
    for t in range(Ld.shape[0]):
        j = np.isfinite(Od[t])
        cand = np.flatnonzero(gate & ~j)
        top = cand[np.argsort(-np.abs(Ld[t, cand]), kind="stable")[:K]]
        lt = Ld[t].copy(); lt[top] *= call_boost                    # larger amplitude only on the called genes
        e = np.exp(lt); e[j] = np.exp(Od[t, j])
        T = np.round(full_mean * e * n)                             # pinned integer group totals
        W = share * T[None, :]
        lam = lam_rank[:len(top)]
        if rank_score is not None:                                  # strength order by confidence, membership by |d|
            lam = np.empty(len(top)); lam[np.argsort(-rank_score[t, top], kind="stable")] = lam_rank[:len(top)]
        up = top[Ld[t, top] > 0]; lu = lam[Ld[t, top] > 0]
        dn = top[Ld[t, top] < 0]; ld = lam[Ld[t, top] < 0]
        if up.size:
            if depth_up > 0:                                        # flatten towards x_i ~ L_i^-a (shallow cells get more)
                u = Lc ** (-depth_up); u = u / u.sum()
                W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (u[:, None] * T[up][None, :])
            else:
                W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (T[up] / n)[None, :]
        if dn.size:
            conc = np.zeros((n, dn.size))
            rows = order_cells[:n_top, dn]                          # top-10 % cells of each gene
            np.put_along_axis(conc, rows, (T[dn] / n_top)[None, :].repeat(n_top, 0), axis=0)
            W[:, dn] = (1 - ld)[None, :] * W[:, dn] + ld[None, :] * conc
        fl = np.floor(W); rem = (T - fl.sum(0)).astype(int); frac = W - fl
        cols = np.flatnonzero(rem > 0)
        if cols.size:
            o = np.argsort(-frac[:, cols], axis=0)
            mask = np.arange(n)[:, None] < rem[cols][None, :]
            add = np.zeros((n, cols.size)); np.put_along_axis(add, o, mask.astype(float), axis=0)
            fl[:, cols] += add
        b = sparse.csr_matrix(fl.astype(np.float32)); b.eliminate_zeros(); out.append(b)
    return sparse.vstack(out, format="csr")


def agreement_score(panel, d):
    """|d| x sign agreement of the source families with d (fraction agreeing, rescaled 0..1), on panel genes."""
    tg = H.PANELS[panel]["targets"]; r = H.real(panel); ok = r["gmap"] >= 0
    fl = H.fam_lfc(tg, H.EXCL[panel]); g = r["gmap"][ok]
    sd = np.sign(d[:, g]); agree = np.zeros(sd.shape); n = np.zeros(sd.shape)
    for f, L in fl.items():
        v = L[:, g]; m = np.isfinite(v); agree += m * (np.sign(np.nan_to_num(v)) == sd); n += m
    frac = np.where(n > 0, agree / np.maximum(n, 1), 0.5)
    out = np.zeros((len(tg), len(r["genes"]))); out[:, ok] = np.abs(d[:, g]) * np.clip(2 * frac - 1, 0, 1)
    return out


def score_dual(panel, gamma, K, seed=0, boost=1.0, depth_down=False, depth_up=0.0, agree=False):
    tag = (f"dual_g{gamma:g}_k{K}" + (f"_b{boost:g}" if boost != 1.0 else "") + ("_dd" if depth_down else "")
           + (f"_du{depth_up:g}" if depth_up else "") + ("_ag" if agree else "") + (f"@s{seed}" if seed else ""))
    out = H.D / "scores" / tag / panel
    if (out / "result.json").exists():
        return json.loads((out / "result.json").read_text())
    out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    d, ov, _ = H.arm_delta("blend30_e", panel)
    r = H.real(panel); P = H.PANELS[panel]; ok = r["gmap"] >= 0
    Ld = H.on_panel(gamma * d, panel)
    Od = np.full(Ld.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
    gate = r["cpm"] >= 5.0
    rs = agreement_score(panel, d) if agree else None
    X = emit_dual(r["ctrl"], Ld, Od, P["n_cells"], K, seed, gate, call_boost=boost, depth_down=depth_down,
                  depth_up=depth_up, rank_score=rs)
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
    env = dict(os.environ, NUMBA_NUM_THREADS="2", OMP_NUM_THREADS="2")
    with (out / "scorer.log").open("w") as g:
        subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
    agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
    raw = {k: float(agg.loc["mean", k]) for k in H.S.CONST}
    comp, scaled = H.S.composite(raw)
    res = dict(arm=tag, panel=panel, gamma=gamma, K=K, boost=boost, depth_down=depth_down, raw=raw, scaled=scaled, composite=comp, seconds=round(time.time() - t0, 1))
    (out / "result.json").write_text(json.dumps(res, indent=1)); pp.unlink()
    print(json.dumps(dict(arm=tag, panel=panel, composite=round(comp, 4), **{k: round(v, 4) for k, v in raw.items()},
                          secs=res["seconds"])), flush=True)
    return res


if __name__ == "__main__":
    panel = _argv[1]; a = _argv[2:]
    for spec in a:                                                  # gamma:K[:boost[:dd]][@sSEED]
        spec_, _, sd = spec.partition("@s"); seed = int(sd) if sd else 0
        f = spec_.split(":"); g, k = float(f[0]), int(f[1]); bst = float(f[2]) if len(f) > 2 else 1.0
        fl_ = f[3] if len(f) > 3 else ""                            # flags: dd, du<a>, ag  (e.g. "dd,du0.5,ag")
        opts = fl_.split(",") if fl_ else []
        dd = "dd" in opts; ag = "ag" in opts
        du = next((float(o[2:]) for o in opts if o.startswith("du")), 0.0)
        try:
            score_dual(panel, g, k, seed=seed, boost=bst, depth_down=dd, depth_up=du, agree=ag)
        except Exception as ex:
            import traceback; traceback.print_exc(); print(json.dumps(dict(spec=spec, error=repr(ex))), flush=True)
    print("ALL_DONE", flush=True)
