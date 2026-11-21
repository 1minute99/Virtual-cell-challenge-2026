"""Score a transfer_v2 arm on a proxy panel with the OFFICIAL cell-eval2 vcc2026 preset.

Usage: tv2_score.py ARMS_JSON   (list of arms; each arm: tag, panel H1|KOLF, cfg overrides, alpha, mode, seed)
Composite = mean of the six members scaled with the official leaderboard constants (base/anchor fitted on
the 2026-09-22 API: pds .500/.9515, mse .984/.043 clamped [0,1], nmae 1.0008/.3921, fid .5126/.8065,
reach .0790/.9653, jac .0303/.4010).
"""
import json, os, pathlib, subprocess, sys, time
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2.moments import jackknife_correction  # noqa: F401  (before scripts/ joins sys.path)
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_model as TM, tv2_emit as EM

R = TM.R
OUT = R / "reports/transfer_v2_2026-09-23/scores"
CE = pathlib.Path(sys.executable).with_name("cell-eval2")
PANELS = {"H1": (R / "reports/transfer_v2_2026-09-23/h1/real_thin.h5ad", 300, "H1"),
          "KOLF": (R / "reports/b19_sweep2_2026-09-18/full_metrics/KOLF/real_disjoint.h5ad", 200, "KOLF")}
CONST = {"pds_cosine": (0.5, 0.9515), "expr_mse_unbiased_capped_norm": (0.984, 0.043), "de_wilcoxon_lfc_nmae": (1.0008, 0.3921),
         "de_wilcoxon_direction_fidelity_yield_raw": (0.5126, 0.8065), "de_wilcoxon_direction_reach_raw": (0.0790, 0.9653),
         "de_wilcoxon_sig_jaccard": (0.0303, 0.4010)}


def composite(raw):
    s = {}
    for k, (b, a) in CONST.items():
        v = (raw[k] - b) / (a - b)
        if k.startswith("expr_mse"):
            v = min(max(v, 0.0), 1.0)
        s[k] = v
    return float(np.mean(list(s.values()))), s


def build_cfg(over):
    cfg = TM.Cfg()
    for k, v in over.items():
        setattr(cfg, k, v)
    return cfg


def score_arm(arm):
    tag = arm["tag"]; out = OUT / tag
    if (out / "result.json").exists():
        return json.loads((out / "result.json").read_text())
    out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    real_path, n_cells, family = PANELS[arm["panel"]]
    real = ad.read_h5ad(real_path)
    lab = real.obs["target"].astype(str).values
    targets = sorted(set(lab) - {"non-targeting"})
    ctrl = real.X[lab == "non-targeting"].tocsr()
    genes = real.var_names.astype(str).tolist()
    del real
    excl = {family}
    parts = arm.get("mix") or [{"w": 1.0, "cfg": arm.get("cfg", {}), "alpha": arm.get("alpha", 1.0)}]
    Loff = 0.0; Lov = None
    for part in parts:   # optional linear blend of several fusion configs (each with its own amplitude)
        cfg = build_cfg(part.get("cfg", {}))
        L, prom, tidx = TM.predict_lfc(targets, cfg, exclude=excl)
        d, ov = TM.split(L, prom, tidx, part["alpha"], cfg)
        Loff = Loff + part["w"] * d
        Lov = ov if Lov is None else np.where(np.isfinite(Lov), Lov, ov)
    gmap = np.array([TM.GIC.get(TM.canon(g), -1) for g in genes])
    Ld = np.zeros((len(targets), len(genes))); ok = gmap >= 0
    Ld[:, ok] = Loff[:, gmap[ok]]
    Od = np.full((len(targets), len(genes)), np.nan); Od[:, ok] = Lov[:, gmap[ok]]
    X = EM.emit(ctrl, Ld, n_cells, arm.get("seed", 0), arm.get("mode", "iid_mult"), override=Od)
    pred = ad.AnnData(X=sparse.vstack([X, ctrl], format="csr"),
                      obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(targets, n_cells), *["non-targeting"] * ctrl.shape[0]])},
                                       index=[f"p{i}" for i in range(X.shape[0] + ctrl.shape[0])]),
                      var=pd.DataFrame(index=genes))
    pp = out / "pred.h5ad"; pred.write_h5ad(pp); del pred, X
    cache = R / f"reports/transfer_v2_2026-09-23/cache_real_{arm['panel']}"
    cmd = [str(CE), "run", "-ap", str(pp), "-ar", str(real_path), "--preset", "vcc2026", "--pert-col", "target",
           "--control", "non-targeting", "--input-type", "counts", "-o", str(out / "cell_eval2"), "--cache-real", str(cache),
           "--set", "de.backend=pdex", "--set", "num_threads=2", "--set", "gather_threads=2"]
    env = dict(os.environ, NUMBA_NUM_THREADS="2", OMP_NUM_THREADS="2")
    with (out / "scorer.log").open("w") as g:
        subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
    agg = pd.read_csv(out / "cell_eval2/agg_results.csv", index_col=0)
    raw = {k: float(agg.loc["mean", k]) for k in CONST}
    comp, scaled = composite(raw)
    res = dict(arm=arm, raw=raw, scaled=scaled, composite=comp, seconds=round(time.time() - t0, 1))
    (out / "result.json").write_text(json.dumps(res, indent=1))
    pp.unlink()
    print(json.dumps(dict(tag=tag, composite=round(comp, 4), **{k.split("_")[-1] if not k.startswith("expr") else "mse": round(v, 4) for k, v in raw.items()}, secs=res["seconds"])), flush=True)
    return res


if __name__ == "__main__":
    arms = json.loads(pathlib.Path(sys.argv[1]).read_text())
    for arm in arms:
        try:
            score_arm(arm)
        except Exception as ex:
            print(json.dumps(dict(tag=arm["tag"], error=repr(ex))), flush=True)
    print("ALL_DONE", flush=True)
