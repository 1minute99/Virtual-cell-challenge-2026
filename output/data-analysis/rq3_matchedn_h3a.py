"""EXPLORATORY, MATCHED-N sensitivity for H3a (plan input/ideas/zero-shot-limits-context-complexity-plan.md: H3a,
operational definition 7, Task 16). NOT part of the registered analysis; the registered verdict (results/rq3_h3a.json)
is untouched. Trigger: the T15+T16 review (knowledge/_meta/log.md, 2026-11-05): U_j is built from ALL control cells of
panel j (800-12,000 cells) while U_dest is built from one third C_pc (267-4,000 cells) -- a sample-size mismatch.

The estimator, residuals, splits, target halves, bootstrap draws (W) and two-level control resampling are those of
rq3_energy.py, imported as a library and not modified. Only the SUBSPACES change:
  reg      registered design (U_dest = top-k PCs of C_pc; U_j = cached top-k PCs of all controls of j). Reproduced
           here and asserted against results/rq3_energy.parquet (validation of this script, not a result).
  mn       matched-n: for each comparison j, n' = min(n_d, n_j) (n_d = |C_pc| with L > 0, n_j = #controls of j with
           L > 0); U_j from a random n'-subsample of j's controls; U_dest from C_pc if n' = n_d, else from a random
           n'-subsample of C_pc (recorded). 5 subsample reps per (dest, split seed, j); pair (dest, j, r) differences
           F(U_dest^(j,r)) - F(U_j^(r)) averaged over j and r  -> D_other.
  c267     common n = 267 (the smallest C_pc): every subspace (dest and the 6 others) from exactly 267 cells; 5 reps per
           (dest, seed); D_other = mean_r [F(U_dest,r) - mean_j F(U_j,r)].
  c800     common n = 800: as c267 with 800 cells. NOT attainable for the RN destinations (C_pc = 267 cells < 800):
           reported as NA there (no mixing of sizes).
Assumptions (stated; checked where testable)
  X1 PC-estimation noise is the only quantity the mismatch changes; the residual estimator, the target halves and the
     baselines (C_a, C_b) are identical to the registered run (asserted via the 'reg' reproduction).
  X2 Subsamples are drawn without replacement from cells with library size L > 0 (the registered PC inputs drop L = 0);
     subsample indices are sorted, so a subsample equal to the whole pool reproduces the registered PCs exactly
     (top_pcs random_state 0) and is reused instead of being recomputed.
  X3 Bootstrap: identical W (default_rng(20261105), B = 2,000) and two-level control resampling
     (default_rng([20261105, seed, 17])) as rq3_energy; the CIs are CONDITIONAL on the drawn subsamples (subsampling
     variability is reported separately as the SD of the per-rep point estimates, not folded into the CI).
  X4 Exchangeability of reps: the 5 reps are i.i.d. subsamples; their mean is the matched-n estimate.
Seeds: subsamples default_rng([20261105, 31|32|33|34, dest index, split seed, j index, rep, n]) (see code); PCs
  random_state 0; bootstrap as above.
CLI:  python rq3_matchedn_h3a.py --panels H1,OR_HCT116X     (per-dest runs -> results/rq3_matchedn/h3a_<dest>.parquet)
      python rq3_matchedn_h3a.py --collect                  (-> results/rq3_h3a_matchedn.{json,parquet,png})
Environment: ~/miniconda3/envs/ct_transfer; nice -n 10; OMP/MKL 4 threads.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import argparse, json, pathlib, sys, time
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402
import rq3_energy as E  # noqa: E402

R = common.RESULTS
OUT = R / "rq3_matchedn"
PANELS = common.PANELS7
SEEDS = E.SEEDS; KS = E.KS; KMAX = E.KMAX; B = E.B_BOOT
NREP = 5
COMMON_N = (267, 800)
LABEL = "exploratory, matched-n"
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def rng_(*parts):
    return np.random.default_rng([20261105, *[int(p) for p in parts]])


def sub(pool, n, rng):
    return pool if n >= len(pool) else np.sort(rng.choice(pool, n, replace=False))


def pcs(X, L, idx):
    return E.top_pcs(E.cell_log(X[idx], L[idx]), KMAX, seed=0)[0]


def W_matrix(T):
    rng = np.random.default_rng(E.BOOT_SEED)
    W = np.zeros((B + 1, T)); W[0] = 1.0
    draws = rng.integers(0, T, size=(B, T))
    for b in range(B):
        W[b + 1] = np.bincount(draws[b], minlength=T)
    return W


def F_seed(st, W, two_level):
    """F for every subspace of one seed's boot_store: (B+1, U, nk). Same algebra as rq3_energy.summarize2."""
    a, b = st["a"], st["b"]; T, U, K = a.shape; n = T
    num = (W @ (a * b).reshape(T, -1)).reshape(B + 1, U, K)
    den = W @ st["den_t"]
    if two_level:
        Sa = (W @ a.reshape(T, -1)).reshape(B + 1, U, K); Sb = (W @ b.reshape(T, -1)).reshape(B + 1, U, K)
        x, y = st["x"], st["y"]
        num = num - x * Sb - y * Sa + n * x * y
        den = den - (W * st["PA"]).sum(1) - (W * st["PB"]).sum(1) + n * st["dd"]
    num = np.cumsum(num, axis=2)[:, :, [k - 1 for k in KS]]
    return num / den[:, None, None], den


def run_dest(d):
    di = PANELS.index(d)
    # ---- control cells of every panel on G-cap (counts, L) ----
    ctrl = {}
    for j in PANELS:
        P = E.PanelData(j)
        try:
            X, L = P.read(P.ctrl_rows)
        finally:
            P.close()
        ctrl[j] = (X, L, np.flatnonzero(L > 0))
    U_all = {j: E.other_panel_pcs(j) for j in PANELS if j != d}
    others = [j for j in PANELS if j != d]
    # ---- residuals exactly as rq3_energy.run(mode='real', pred='prem3', ctrl_ab=True) ----
    P = E.PanelData(d)
    try:
        cX, cL = ctrl[d][0], ctrl[d][1]
        counts = {t: int((P.lab == t).sum()) for t in P.targets()}
        ptg, _, _ = E.load_pred(d, "prem3", P.gcap)
        tg = [t for t in ptg if counts.get(t, 0) >= 2]
        sel = [list(ptg).index(t) for t in tg]
        dpred = {}
        for s in SEEDS:
            t2, d2, _ = E.load_pred(d, f"prem3_ctrlAB_seed{s}", P.gcap)
            assert list(t2) == list(ptg)
            dpred[s] = d2[sel]
        halves = {s: [] for s in SEEDS}
        for t in tg:
            X, L = P.read(np.flatnonzero(P.lab == t))
            for s in SEEDS:
                halves[s].append(E.target_half_sums(X, L, s, t))
        log(f"[{d}] {len(tg)} targets read")
        splits = {s: E.ctrl_split(d, s, len(P.ctrl_rows), "truth", False, P.ctrl_rows)[0] for s in SEEDS}
    finally:
        P.close()
    T = len(tg); W = W_matrix(T)
    meta = dict(panel=d, label=LABEL, n_targets=T, n_ctrl_pool={j: int(len(ctrl[j][2])) for j in PANELS}, seeds={})
    pairs = {}                    # variant -> list of (dest_name, oth_name, j, r) -- same layout in every seed
    Fs = {"target": [], "two_level": []}; dens = {"target": [], "two_level": []}
    names_ref = None
    for s in SEEDS:
        ts = time.time()
        split = splits[s]
        pc = split[0][cL[split[0]] > 0]; n_d = len(pc)
        subs = {"reg:dest": pcs(cX, cL, pc)}
        for j in others:
            subs[f"reg:oth:{j}"] = U_all[j]
        pr = {"reg": [("reg:dest", f"reg:oth:{j}", j, -1) for j in others]}
        sm = dict(n_d=n_d, mn_n={}, mn_dest_subsampled={})
        # matched-n per comparison
        pr["mn"] = []
        for j in others:
            Xj, Lj, poolj = ctrl[j]; ji = PANELS.index(j)
            n_ = min(n_d, len(poolj)); sm["mn_n"][j] = n_; sm["mn_dest_subsampled"][j] = bool(n_ < n_d)
            for r in range(NREP):
                o = sub(poolj, n_, rng_(31, di, s, ji, r, n_))
                on = f"mn:{j}:{r}:oth"
                subs[on] = U_all[j] if len(o) == len(poolj) else pcs(Xj, Lj, o)
                if n_ < n_d:
                    dn = f"mn:{j}:{r}:dest"
                    subs[dn] = pcs(cX, cL, sub(pc, n_, rng_(32, di, s, ji, r, n_)))
                else:
                    dn = "reg:dest"
                pr["mn"].append((dn, on, j, r))
        # common n
        for n_ in COMMON_N:
            v = f"c{n_}"
            if n_d < n_:
                sm[f"{v}_NA"] = f"C_pc has {n_d} cells < {n_}"; continue
            pr[v] = []
            for r in range(NREP):
                dsel = sub(pc, n_, rng_(33, di, s, r, n_))
                dn = "reg:dest" if len(dsel) == n_d else f"{v}:{r}:dest"
                if dn != "reg:dest":
                    subs[dn] = pcs(cX, cL, dsel)
                for j in others:
                    Xj, Lj, poolj = ctrl[j]; ji = PANELS.index(j)
                    o = sub(poolj, n_, rng_(34, di, s, ji, r, n_))
                    on = f"{v}:{r}:oth:{j}"
                    subs[on] = U_all[j] if len(o) == len(poolj) else pcs(Xj, Lj, o)
                    pr[v].append((dn, on, j, r))
        if names_ref is None:
            names_ref = list(subs); pairs = pr
        else:
            assert list(subs) == names_ref and pr == pairs, "subspace layout differs across seeds"
        rA, rB = E.residuals(halves[s], cX, cL, split, dpred[s])
        rng = np.random.default_rng([E.BOOT_SEED, s, 17])
        dA = E.ctrl_boot_deltas(cX, cL, split[1], B, rng); dB = E.ctrl_boot_deltas(cX, cL, split[2], B, rng)
        st = E.boot_store(rA, rB, subs, dA, dB, KMAX); del dA, dB
        for bt in ("target", "two_level"):
            F, den = F_seed(st, W, bt == "two_level"); Fs[bt].append(F); dens[bt].append(den)
        del st
        meta["seeds"][s] = sm
        log(f"[{d}] seed {s}: n_d {n_d}, {len(subs)} subspaces, {time.time() - ts:.0f}s")
    idx = {nm: i for i, nm in enumerate(names_ref)}
    rows = []
    for bt in ("target", "two_level"):
        F = np.stack(Fs[bt], 1)                                   # (B+1, S, U, nk)
        den = np.stack(dens[bt], 1)
        for v, pl in pairs.items():
            id_ = np.array([idx[p[0]] for p in pl]); io = np.array([idx[p[1]] for p in pl])
            Fd = F[:, :, id_]; Fo = F[:, :, io]                    # (B+1, S, P, nk)
            q = {"F_dest": Fd.mean(2).mean(1), "F_other_mean": Fo.mean(2).mean(1), "D_other": (Fd - Fo).mean(2).mean(1)}
            for j in others:
                m = np.array([p[2] == j for p in pl])
                q[f"F_other:{j}"] = Fo[:, :, m].mean(2).mean(1)
                q[f"F_dest_vs:{j}"] = Fd[:, :, m].mean(2).mean(1)
                q[f"D_other:{j}"] = (Fd[:, :, m] - Fo[:, :, m]).mean(2).mean(1)
            reps = sorted({p[3] for p in pl})
            rep_pts = {r: (Fd[0, :, np.array([p[3] == r for p in pl])] - Fo[0, :, np.array([p[3] == r for p in pl])])
                       for r in reps}
            for name, val in q.items():
                for c, k in enumerate(KS):
                    bs = val[1:, c]
                    row = dict(panel=d, variant=v, quantity=name, k=k, boot=bt, point=float(val[0, c]),
                               se=float(bs.std(ddof=1)), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)),
                               frac_le0=float((bs <= 0).mean()), n_targets=T, n_pairs=len(pl), B=B,
                               den_sign_ok=bool((den[0] > 0).all()), label=LABEL)
                    if name == "D_other" and len(reps) > 1:
                        pts = [float(rep_pts[r].mean(axis=(0, 1))[c]) for r in reps]   # mean over seeds and j at rep r
                        row.update(rep_sd=float(np.std(pts, ddof=1)), rep_min=min(pts), rep_max=max(pts))
                    rows.append(row)
    import pandas as pd
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_parquet(OUT / f"h3a_{d}.parquet", index=False)
    meta["runtime_s"] = time.time() - t0
    (OUT / f"h3a_{d}.json").write_text(json.dumps(meta, indent=1, default=str))
    # ---- validation: 'reg' reproduces the registered rq3_energy.parquet ----
    Eg = pd.read_parquet(R / "rq3_energy.parquet")
    Eg = Eg[(Eg.pred == "prem3_ctrlAB") & (Eg.space == "gcap") & (Eg.stratum == "all") & (Eg.panel == d)]
    worst = 0.0
    for _, rr in df[df.variant == "reg"].iterrows():
        if rr.quantity not in ("F_dest", "F_other_mean", "D_other"):
            continue
        g = Eg[(Eg.quantity == rr.quantity) & (Eg.k == rr.k) & (Eg.boot == rr.boot)].iloc[0]
        worst = max(worst, abs(g.point - rr.point), abs(g.lo - rr.lo), abs(g.hi - rr.hi))
    log(f"[{d}] reg vs registered max |diff| (point/lo/hi) = {worst:.2e}")
    assert worst < 1e-8, f"{d}: reg does not reproduce the registered run ({worst})"
    meta["reg_reproduction_maxabs"] = worst
    (OUT / f"h3a_{d}.json").write_text(json.dumps(meta, indent=1, default=str))
    log(f"[{d}] done")


def collect():
    import pandas as pd
    fs = [OUT / f"h3a_{p}.parquet" for p in PANELS if (OUT / f"h3a_{p}.parquet").exists()]
    df = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    df.to_parquet(R / "rq3_h3a_matchedn.parquet", index=False)
    panels = [p for p in PANELS if p in set(df.panel)]
    g = lambda p, v, q, k, bt: df[(df.panel == p) & (df.variant == v) & (df.quantity == q) & (df.k == k) & (df.boot == bt)]
    out = dict(label=LABEL, status="exploratory (not pre-registered); the registered H3a verdict in results/rq3_h3a.json is unchanged",
               plan="input/ideas/zero-shot-limits-context-complexity-plan.md", script="output/data-analysis/rq3_matchedn_h3a.py",
               variants=dict(reg="registered design reproduced (U_dest from C_pc, U_j from all controls of j)",
                             mn="matched n per comparison: n' = min(|C_pc|, n_j); 5 reps",
                             c267="every subspace from 267 cells; 5 reps", c800="every subspace from 800 cells; 5 reps; NA for RN destinations"),
               table={}, meta={})
    for p in panels:
        m = json.loads((OUT / f"h3a_{p}.json").read_text()); out["meta"][p] = m
        out["table"][p] = {}
        for k in KS:
            rec = {}
            for v in ("reg", "mn", "c267", "c800"):
                t = g(p, v, "D_other", k, "target"); w = g(p, v, "D_other", k, "two_level")
                if len(t) == 0:
                    rec[v] = None; continue
                t = t.iloc[0]; w = w.iloc[0]
                rec[v] = dict(D_other=t.point, target_ci=[t.lo, t.hi], two_level_ci=[w.lo, w.hi],
                              F_dest=float(g(p, v, "F_dest", k, "target").iloc[0].point),
                              F_other_mean=float(g(p, v, "F_other_mean", k, "target").iloc[0].point),
                              rep_sd=None if not np.isfinite(t.get("rep_sd", np.nan)) else float(t.rep_sd),
                              per_line_D={j: float(g(p, v, f"D_other:{j}", k, "target").iloc[0].point)
                                          for j in PANELS if j != p})
            out["table"][p][f"k{k}"] = rec
    jd = lambda o: json.dumps(o, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    (R / "rq3_h3a_matchedn.json").write_text(jd(out))
    # ---- plot: D_other at k = 8 per panel and variant ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.2))
    cols = dict(reg="#4c4c4c", mn="#2a6fdb", c267="#d9822b", c800="#3a9d5d")
    off = dict(reg=-0.27, mn=-0.09, c267=0.09, c800=0.27)
    for v in cols:
        xs, ys, lo, hi = [], [], [], []
        for i, p in enumerate(panels):
            r = out["table"][p]["k8"][v]
            if r is None:
                continue
            xs.append(i + off[v]); ys.append(r["D_other"]); lo.append(r["D_other"] - r["target_ci"][0]); hi.append(r["target_ci"][1] - r["D_other"])
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=4, color=cols[v], capsize=2, label=v)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(range(len(panels))); ax.set_xticklabels(panels, rotation=20)
    ax.set_ylabel("D_other at k = 8 (target-bootstrap 95% CI)")
    ax.set_title("H3a D_other: registered vs matched-n (exploratory, matched-n)")
    ax.legend(frameon=False, ncol=4, fontsize=8)
    fig.tight_layout(); fig.savefig(R / "rq3_h3a_matchedn.png", dpi=150)
    for p in panels:
        r = out["table"][p]["k8"]
        print(p, " | ".join(f"{v} {r[v]['D_other']:+.4f} [{r[v]['target_ci'][0]:+.4f},{r[v]['target_ci'][1]:+.4f}]" if r[v] else f"{v} NA"
                            for v in ("reg", "mn", "c267", "c800")))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels"); ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    if a.collect:
        collect()
    else:
        for p in a.panels.split(","):
            run_dest(p)
