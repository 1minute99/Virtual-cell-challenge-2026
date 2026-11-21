"""Task 16 synthetic calibration of the target vs two-level bootstrap (true D = 0, isotropic residual).
Usage: python rq3_t16_calib.py N_CTRL N_REPLICATES N_SPLIT_SEEDS. Seeds: worlds 500+rep, residuals 600+rep; bootstrap as rq3_energy.
Plan: input/ideas/zero-shot-limits-context-complexity-plan.md (Task 16). Reported in knowledge/synthesis/rq3_h3a-assumptions.md."""
import sys, numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent)); sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent / "tests"))
import test_rq3_energy as T, rq3_energy as E
nctrl = int(sys.argv[1]); R = int(sys.argv[2]); S = int(sys.argv[3])
qs = ("D_perm", "D_other", "N_perm", "N_other")
res = {q: dict(pt=[], ct=[], c2=[], wt=[], w2=[]) for q in qs}
for rep in range(R):
    w = T.World(seed=500 + rep, n_ctrl=nctrl, T=100)
    rho = np.random.default_rng(600 + rep).normal(0, 1 / np.sqrt(T.G), (len(w.names), T.G))
    st = T._stores(w, rho, seeds=tuple(range(S)), n_perm=3, B=300)
    a, b = E.summarize2(st, T.KS, B=300), E.summarize2(st, T.KS, B=300, two_level=True)
    for q in qs:
        ra, rb = T._row(a, q), T._row(b, q)
        d = res[q]; d["pt"].append(ra["point"])
        d["ct"].append(ra["lo"] <= 0 <= ra["hi"]); d["c2"].append(rb["lo"] <= 0 <= rb["hi"])
        d["wt"].append((ra["hi"] - ra["lo"]) / 3.92); d["w2"].append((rb["hi"] - rb["lo"]) / 3.92)
for q, d in res.items():
    print(f"{q:8s} n_ctrl {nctrl} seeds {S}: sd(pt) {np.std(d['pt'], ddof=1):.3g} mean {np.mean(d['pt']):+.3g} | "
          f"CIw/3.92 target {np.median(d['wt']):.3g} two-level {np.median(d['w2']):.3g} | cover(0) target {np.mean(d['ct']):.2f} two-level {np.mean(d['c2']):.2f}")
