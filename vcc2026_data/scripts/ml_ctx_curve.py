"""HOW MUCH IS ONE MORE SOURCE CONTEXT WORTH? The curve PDS(k) over the number of source families.

Section 12 established that we already beat PERFECT same-lab transfer from THREE lines (RN_Jurkat: three real source
lines through our pipeline give PDS 0.8974, we get 0.9078), and concluded that the remaining gap to the top is "most
consistent with more source CONTEXTS". That conclusion was never quantified: nobody measured what our pipeline gains per
added family, so "acquire more contexts" and "estimate better from the ones we have" have never been priced against each
other. This measures it directly.

For each k, sample R random subsets of the panel's AVAILABLE families (its own family is always excluded, as everywhere
in the harness), rebuild both recipe parts and the neighbour library restricted to that subset, and score with the exact
offline evaluator. PDS is per-target scale-free, so it is the honest readout when the family count changes the amplitude.

Read the SLOPE at the top of the curve: it is the expected PDS from one more context, and therefore the price cap on
acquiring one. Usage: ml_ctx_curve.py PANEL..."""
import sys, pathlib, json, itertools
import numpy as np
from cell_eval2.moments import jackknife_correction  # noqa
argv = list(sys.argv); sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import fh_harness as H, sch_sumeval as SE, sch_emit as EMT, tv2_emit as EM
from sklearn.utils.extmath import randomized_svd
H.CACHE = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/src_ablation/cache_ctx"); H.CACHE.mkdir(exist_ok=True)
EMT.S4_LINE["OR_HCT116"] = "HCT116"
FAMS = ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1"]
R = 4
out = {}
for p in argv[1:]:
    r = H.real(p); ok = r["gmap"] >= 0; n = H.PANELS[p]["n_cells"]; ctrl = r["ctrl"]; tg = r["targets"]; genes = r["genes"]
    avail = [f for f in FAMS if f not in set(H.EXCL[p])]
    print(f"\n=== {p}: available families {avail}", flush=True)
    dlib = H.lib_delta(p, set(H.EXCL[p])); Dl = np.zeros((dlib.shape[0], len(genes))); Dl[:, ok] = dlib[:, r["gmap"][ok]]
    Dl = Dl[np.abs(Dl).sum(1) > 0]; _, _, V3 = randomized_svd(Dl - Dl.mean(0), n_components=3, random_state=0); del dlib, Dl

    def run(sub):
        wm = {f: 0.0 for f in avail if f not in sub}
        d, ov = EMT.delta_fs(p, 0.5, 0.3, wmul=wm or None)
        L = H.on_panel(d, p); Od = np.full(L.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
        mt = np.nanmean(np.where(np.isfinite(Od), np.nan, L), 0); X = L - 0.5 * np.nan_to_num(mt)[None, :]
        P_ = X @ V3.T; f_ = float((P_ ** 2).sum() / (X ** 2).sum())
        X = X - 0.75 * float(np.clip((0.25 - f_) / 0.10, 0, 1)) * P_ @ V3
        gk, _ = EMT.kink_gamma(ctrl, L, Od, n, factor=1.17, margin=0.8)
        Bt = EM.budget_B(ctrl, L, Od, n)(max(1.0, gk))
        s = EM.budget_match(ctrl, X, Od, n, Bt)
        return SE.evaluate(p, tg, genes, ctrl, X * s, Od, n)

    rng = np.random.default_rng(0); rec = {}
    for k in range(1, len(avail) + 1):
        subs = list(itertools.combinations(avail, k))
        if len(subs) > R:
            subs = [tuple(sorted(rng.choice(avail, size=k, replace=False))) for _ in range(R)]
            subs = list(dict.fromkeys(subs))
        vals = []
        for sub in subs:
            try:
                o = run(set(sub)); vals.append((o["pds"], o["mse"], "+".join(sub)))
            except Exception as e:
                print(f"    k={k} {sub} FAILED {e!r}", flush=True)
        if not vals: continue
        pds = np.array([v[0] for v in vals]); mse = np.array([v[1] for v in vals])
        rec[k] = dict(pds_mean=float(pds.mean()), pds_sd=float(pds.std()), pds_max=float(pds.max()),
                      mse_mean=float(mse.mean()), n=len(vals), best="+".join(sorted(vals[int(np.argmax(pds))][2].split("+"))))
        print(f"    k={k}  PDS {pds.mean():.4f} +- {pds.std():.4f}  (max {pds.max():.4f}: {rec[k]['best']})  "
              f"MSE {mse.mean():.4f}  n={len(vals)}", flush=True)
    ks = sorted(rec)
    if len(ks) >= 2:
        d1 = rec[ks[-1]]["pds_mean"] - rec[ks[-2]]["pds_mean"]
        print(f"    marginal PDS of the LAST family added: {d1:+.4f}  -> composite {d1*2.234/6:+.5f}", flush=True)
    out[p] = rec
json.dump(out, open("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/ctx_curve.json", "w"), indent=1)
