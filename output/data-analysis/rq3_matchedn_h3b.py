"""EXPLORATORY, MATCHED-N sensitivity for H3b (plan input/ideas/zero-shot-limits-context-complexity-plan.md: H3b,
Task 17). NOT part of the registered analysis; results/rq3_h3b.json (registered verdict) is untouched. Offline pack
evaluator only (tasks/vcc_sum/evaluator.Pack = sch_sumeval numerically); no cell-eval2 run.

Trigger: T15+T16 review (knowledge/_meta/log.md, 2026-11-05). In the registered Task 17 run every control-PC basis
follows the builder recipe "min(6,000, n_ctrl) random control cells", i.e. 6,000 cells for H1 / OR_HCT116X / KOLF
(12,000 / 12,000 / 11,600 controls) but all 800 cells for the four RN lines -- so M3_dest and M3_other_j are built from
different numbers of cells whenever one of the pair is an RN line.

Design (budget-invariant M3-type arms only: M3_dest and M3_other_j, each followed by the same budget match):
  every control-PC basis (destination and the 6 other lines) is estimated from exactly n cells, n in
  {800 (= the smallest control set among the 7 panels, so "smallest" and "800" coincide), 267 (extra: the H3a
  common-n, the smallest C_pc)}; 5 subsample reps r per line (seed default_rng([20261105, 41, line index, n, r]));
  rep r of the destination is paired with rep r of each other line. Per panel: M3_dest_mn = mean over reps of the
  per-target PDS; M3_other_mn = mean over the 6 lines x 5 reps. Delta = M3_dest_mn - M3_other_mn, 7-panel mean
  (panels weighted equally), paired per-panel target bootstrap (B = 2,000, default_rng(20261105), the registered
  boot_delta algorithm). For the RN lines n = 800 means all cells -> the 5 reps are identical there.
Assumptions
  Y1 Chain, weak half, k = 8, eta = 0.75, budget match and evaluator exactly as rq3_corrections.cmd_offline
     (imported, not modified); only the cell count behind V changes.
  Y2 The destination basis is recomputed with the builder recipe (ctrl_official + top_v) and re-orthonormalised on the
     destination axis like M3_other (the registered M3_dest used Task 14's frozen ctrl_V; the 'recipe6000' arm checks
     that this re-implementation reproduces the registered M3_dest PDS).
  Y3 The bootstrap CI is conditional on the drawn subsamples and on the fixed truth; subsampling variability is
     reported as the spread of the per-rep 7-panel mean Delta.
Validation: the 'registered' arms recomputed here (frozen ctrl_V; cached 6,000-cell V of the other lines) must equal
  rq3_h3b.json panel PDS (|d| < 1e-9).
Environment: ~/miniconda3/envs/vcc2026 (python -P); nice -n 10; OMP/MKL/OpenBLAS 4 threads.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import json, pathlib, sys, time
import numpy as np

HERE = pathlib.Path("/home/wonmin-kim/AIVC/output/data-analysis")
sys.path.insert(0, str(HERE))
import rq3_corrections as RC  # noqa: E402

PANELS = RC.PANELS7; NS = (800, 267); NREP = 5; LABEL = "exploratory, matched-n"
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def boot_delta(A, Bm):
    """A, Bm: panel -> per-target PDS arrays (same target order). Registered algorithm (rq3_corrections.boot_delta)."""
    rng = np.random.default_rng(RC.BOOT_SEED); pm = []; bm = []
    for p in PANELS:
        d = A[p] - Bm[p]; T = len(d); idx = rng.integers(0, T, (RC.B_BOOT, T))
        pm.append(d.mean()); bm.append(d[idx].mean(1))
    pm = np.array(pm); bm = np.mean(bm, 0)
    return dict(mean=float(pm.mean()), ci95=[float(np.quantile(bm, 0.025)), float(np.quantile(bm, 0.975))],
                per_panel={p: float(v) for p, v in zip(PANELS, pm)}, n_panels_positive=int((pm > 0).sum()))


def main():
    CM, EV = RC.imports()
    NG = len(CM.official_genes())
    summ = json.loads((RC.PRED / "rq3_build_preds_summary.json").read_text())
    reg = json.loads((RC.RES / "rq3_h3b.json").read_text())["panel_pds"]
    VO = RC.official_subspaces(CM)                                       # cached registered bases (6,000-cell recipe)
    # ---- matched-n bases: V[j][n][r] on the official axis ----
    V = {}; ncell = {}
    for ji, j in enumerate(PANELS):
        tj = time.time(); ctrl, genes = CM.controls(j); gmap = CM.to_official(np.asarray(genes)); ncell[j] = int(ctrl.shape[0])
        V[j] = {n: [RC.top_v(RC.ctrl_official(ctrl, gmap, NG, ncells=n, seed=[20261105, 41, ji, n, r]), 8) for r in range(NREP)]
                for n in NS}
        del ctrl
        log(j, "bases", ncell[j], "controls", round(time.time() - tj, 1), "s")
    pdst = {}; panel = {}; checks = {}
    for p in PANELS:
        tp = time.time()
        z = np.load(RC.PRED / f"{p}_prem3.npz", allow_pickle=True); zm = np.load(RC.PRED / f"{p}_m3.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]; gmap = z["gmap"]
        X0 = z["lfc"].astype(np.float64) / float(z["budget_scale"])
        Od = z["Od"].astype(np.float64); ovj = np.isfinite(Od)
        weak = zm["m3_weak"].astype(bool); w = np.where(weak, RC.ETA, 0.0)
        pk = EV.Pack(RC.PACKS / f"score_{p}.npz")
        if p == "KOLF":
            pk.Od, pk.ovj, pk.Bt = Od, ovj, float(summ["KOLF"]["kolf_fix"]["Bt_fixed"])

        def run(X):
            s = pk.budget_match(X); e = pk.evaluate(X * s)
            return np.asarray(e["pds_t"], float), float(e["pds"])

        basis = lambda Vo: RC.orthonormal_rows(RC.to_panel(Vo, gmap))
        res = {"none": run(X0), "reg_M3_dest": run(RC.defl(X0, w, zm["ctrl_V"].astype(np.float64))),
               "recipe6000_M3_dest": run(RC.defl(X0, w, basis(VO[p]["V8"])))}
        oth = [j for j in PANELS if j != p]
        for j in oth:
            res[f"reg_M3_other_{j}"] = run(RC.defl(X0, w, basis(VO[j]["V8"])))
        reg_other = np.mean([res[f"reg_M3_other_{j}"][0] for j in oth], 0)
        ck = dict(reg_M3_dest_minus_registered=res["reg_M3_dest"][1] - reg[p]["M3_dest"],
                  reg_M3_other_minus_registered=float(reg_other.mean()) - reg[p]["M3_other"],
                  none_minus_registered=res["none"][1] - reg[p]["none"],
                  recipe6000_minus_frozen_M3_dest=res["recipe6000_M3_dest"][1] - res["reg_M3_dest"][1])
        assert abs(ck["reg_M3_dest_minus_registered"]) < 1e-9 and abs(ck["reg_M3_other_minus_registered"]) < 1e-9, ck
        checks[p] = ck
        pdst[p] = {"none": res["none"][0], "reg_M3_dest": res["reg_M3_dest"][0], "reg_M3_other": reg_other}
        panel[p] = {"T": len(tg), "none": res["none"][1], "reg_M3_dest": res["reg_M3_dest"][1], "reg_M3_other": float(reg_other.mean()),
                    "recipe6000_M3_dest": res["recipe6000_M3_dest"][1]}
        for n in NS:
            dest_r = [run(RC.defl(X0, w, basis(V[p][n][r])))[0] for r in range(NREP)]
            oth_r = [[run(RC.defl(X0, w, basis(V[j][n][r])))[0] for j in oth] for r in range(NREP)]
            pdst[p][f"n{n}_M3_dest"] = np.mean(dest_r, 0)
            pdst[p][f"n{n}_M3_other"] = np.mean([np.mean(o, 0) for o in oth_r], 0)
            pdst[p][f"n{n}_rep_dest"] = dest_r; pdst[p][f"n{n}_rep_other"] = [np.mean(o, 0) for o in oth_r]
            panel[p][f"n{n}"] = dict(M3_dest=float(np.mean(dest_r)), M3_other=float(np.mean(oth_r)),
                                     M3_dest_rep=[float(x.mean()) for x in dest_r],
                                     M3_other_rep=[float(np.mean(o)) for o in oth_r],
                                     M3_other_per_line={j: float(np.mean([oth_r[r][jj].mean() for r in range(NREP)])) for jj, j in enumerate(oth)})
        log(p, "done", round(time.time() - tp, 1), "s |", {k: round(v, 4) if isinstance(v, float) else v for k, v in panel[p].items() if not isinstance(v, dict)},
            {n: (round(panel[p][f"n{n}"]["M3_dest"], 4), round(panel[p][f"n{n}"]["M3_other"], 4)) for n in NS})
    sel = lambda key: {p: pdst[p][key] for p in PANELS}
    out = dict(label=LABEL, status="exploratory (not pre-registered); registered H3b verdict in results/rq3_h3b.json unchanged",
               plan="input/ideas/zero-shot-limits-context-complexity-plan.md", script="output/data-analysis/rq3_matchedn_h3b.py",
               evaluator="offline ERA pack evaluator (sch_sumeval-equivalent); no cell-eval2",
               n_controls=ncell, note_n="n = 800 is both the smallest control set (RN lines) and the requested n = 800; n = 267 added",
               registered=dict(M3_dest_minus_M3_other=boot_delta(sel("reg_M3_dest"), sel("reg_M3_other")),
                               M3_dest_minus_none=boot_delta(sel("reg_M3_dest"), sel("none"))),
               matched={}, panel_pds=panel, checks=checks,
               settings=dict(k=RC.K, eta=RC.ETA, n=NS, reps=NREP, B=RC.B_BOOT, boot_seed=RC.BOOT_SEED))
    for n in NS:
        rep_means = []
        for r in range(NREP):
            rep_means.append(float(np.mean([pdst[p][f"n{n}_rep_dest"][r].mean() - pdst[p][f"n{n}_rep_other"][r].mean() for p in PANELS])))
        out["matched"][f"n{n}"] = dict(
            M3_dest_minus_M3_other=boot_delta(sel(f"n{n}_M3_dest"), sel(f"n{n}_M3_other")),
            M3_dest_minus_none=boot_delta(sel(f"n{n}_M3_dest"), sel("none")),
            matchedDest_minus_regOther=boot_delta(sel(f"n{n}_M3_dest"), sel("reg_M3_other")),
            regDest_minus_matchedOther=boot_delta(sel("reg_M3_dest"), sel(f"n{n}_M3_other")),
            per_rep_7panel_mean_delta=rep_means, per_rep_sd=float(np.std(rep_means, ddof=1)))
    RC.jdump(RC.RES / "rq3_h3b_matchedn.json", out)
    r0 = out["registered"]["M3_dest_minus_M3_other"]
    log(f"registered dPDS M3_dest - M3_other {r0['mean']:+.5f} [{r0['ci95'][0]:+.5f}, {r0['ci95'][1]:+.5f}] pos {r0['n_panels_positive']}/7")
    for n in NS:
        m = out["matched"][f"n{n}"]["M3_dest_minus_M3_other"]
        log(f"n={n} dPDS M3_dest - M3_other {m['mean']:+.5f} [{m['ci95'][0]:+.5f}, {m['ci95'][1]:+.5f}] pos {m['n_panels_positive']}/7 "
            f"rep sd {out['matched'][f'n{n}']['per_rep_sd']:.5f} per panel {m['per_panel']}")
    log("DONE")


if __name__ == "__main__":
    main()
