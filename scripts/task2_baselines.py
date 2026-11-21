#!/usr/bin/env python
"""Task 2 - Source-effect-copy baseline + floor baselines, scored with Arc cell-eval.

Predictors (all training-free - no gradient is ever taken):
  source_copy_<SRC>  : pred_cell = basal_control_cell + delta_source[target_gene]
  return_control     : pred_cell = basal_control_cell                      (floor)
  train_mean_delta   : pred_cell = basal_control_cell + mean_p delta[p]    (floor)

WHY A SURROGATE: Arc withholds the VCC 2026 validation ground truth (it *is* the
leaderboard) and there is no 2026 challenge training set, so VCC predictions
cannot be scored locally by anyone. We therefore score the identical predictors
on unseen Replogle target cell lines, which do have ground truth.
  source K562_gwps = Replogle K562 genome-scale bulk (the same source used for
  the real VCC predictions), plus each Replogle cell line as an alternate source.
  target = rpe1 / jurkat / hepg2 (unseen contexts); k562 = same-context reference.

Predictions are made per cell on sampled control cells, so the pseudobulk mean
equals control_pseudobulk + delta while retaining realistic cell-level variance
(cell-eval computes differential expression and needs non-degenerate variance).

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task2_baselines.py \
        --source K562_gwps --target rpe1 --seed 0
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genemap import GeneMap
from source_delta import load_source_deltas

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "work/cellline_cache")
OUT = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)
CONTROL = "non-targeting"


def load_cached(cell_line):
    X = np.load(os.path.join(CACHE, f"X_{cell_line}.npy"))
    p = np.load(os.path.join(CACHE, f"pert_{cell_line}.npy"), allow_pickle=True)
    g = np.load(os.path.join(CACHE, "genes.npy"), allow_pickle=True)
    return X, np.asarray(p, dtype=str), [str(x) for x in g]


def deltas_from_cellline(cell_line, gm):
    """Pseudobulk delta per perturbation computed inside one Replogle cell line."""
    X, p, genes = load_cached(cell_line)
    ctrl = X[p == CONTROL].mean(0)
    out = {}
    for q in np.unique(p):
        if q == CONTROL:
            continue
        out[gm.canon(q)] = (X[p == q].mean(0) - ctrl).astype(np.float32)
    del X
    return out, genes


def subsample(X, p, max_per_pert, max_control, rng):
    keep = []
    for q in np.unique(p):
        idx = np.where(p == q)[0]
        cap = max_control if q == CONTROL else max_per_pert
        if len(idx) > cap:
            idx = rng.choice(idx, cap, replace=False)
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    return X[keep], p[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="K562_gwps",
                    help="K562_gwps (genome-scale bulk) or a cell line name")
    ap.add_argument("--target", default="rpe1")
    ap.add_argument("--max-cells-per-pert", type=int, default=64)
    ap.add_argument("--max-control", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--drop-target-gene", action="store_true",
                    help="DEPRECATED: zeroing the target gene in both real and pred "
                         "injects a matched artificial shift and manufactures DE "
                         "agreement. Use --exclude-target-genespace instead.")
    ap.add_argument("--exclude-target-genespace", action="store_true",
                    help="Remove the union of all perturbed genes from the gene "
                         "space before scoring, so no target gene can contribute "
                         "to any metric. This is the correct way to measure how "
                         "much the trivially-predictable target gene inflates scores.")
    a = ap.parse_args()

    from cell_eval import MetricsEvaluator
    t0 = time.time(); gm = GeneMap()

    X, perts, genes = load_cached(a.target)
    rng = np.random.default_rng(a.seed)
    X, perts = subsample(X, perts, a.max_cells_per_pert, a.max_control, rng)

    if a.source == "K562_gwps":
        deltas, info, _ = load_source_deltas(genes, gm)
    else:
        deltas, _ = deltas_from_cellline(a.source, gm)
        info = dict(source_cell_line=a.source, n_perturbations=len(deltas))
    print(f"source={a.source} {info}")

    perts = np.array([gm.canon(q) if q != CONTROL else q for q in perts])
    covered = [q for q in np.unique(perts) if q != CONTROL and q in deltas]
    n_all = len(np.unique(perts)) - 1
    print(f"target={a.target} cells={X.shape[0]} perts={n_all} covered={len(covered)}")
    keep = np.isin(perts, np.array(covered + [CONTROL]))
    X, perts = X[keep], perts[keep]
    mean_delta = np.mean([deltas[q] for q in covered], axis=0)

    if a.exclude_target_genespace:
        tgt_set = {q for q in np.unique(perts) if q != CONTROL}
        keep_g = np.array([i for i, g in enumerate(genes) if g not in tgt_set])
        print(f"excluding {len(genes)-len(keep_g)} perturbed genes from the gene space "
              f"-> {len(keep_g)} genes remain")
        X = X[:, keep_g]
        genes = [genes[i] for i in keep_g]
        deltas = {k: v[keep_g] for k, v in deltas.items()}
        mean_delta = mean_delta[keep_g]

    gpos = {g: i for i, g in enumerate(genes)}
    ctrl_idx = np.where(perts == CONTROL)[0]
    tag = f"{a.source}_to_{a.target}_seed{a.seed}" + ("_notgt" if a.drop_target_gene else "") \
          + ("_noTGspace" if a.exclude_target_genespace else "")

    allres = []
    de_real_path = None   # real DE is identical across predictors in a run -> compute once
    for mode in [f"source_copy_{a.source}", "return_control", "train_mean_delta"]:
        r2 = np.random.default_rng(a.seed)
        P = np.empty_like(X); Xr = X.copy()
        for q in np.unique(perts):
            idx = np.where(perts == q)[0]
            if q == CONTROL:
                P[idx] = X[idx]; continue
            basal = X[r2.choice(ctrl_idx, len(idx), replace=True)]
            d = 0.0 if mode == "return_control" else (
                mean_delta if mode == "train_mean_delta" else deltas[q])
            P[idx] = np.maximum(basal + d, 0.0)
            if a.drop_target_gene and q in gpos:
                P[idx, gpos[q]] = 0.0
                Xr[idx, gpos[q]] = 0.0
        real = ad.AnnData(Xr, obs=pd.DataFrame({"target": perts}), var=pd.DataFrame(index=genes))
        pred = ad.AnnData(P,  obs=pd.DataFrame({"target": perts}), var=pd.DataFrame(index=genes))
        od = os.path.join(OUT, "celleval", f"{tag}_{mode}")
        ev = MetricsEvaluator(adata_pred=pred, adata_real=real, control_pert=CONTROL,
                              pert_col="target", num_threads=-1, prefix=mode,
                              de_real=de_real_path, outdir=od)
        res, agg = ev.compute(profile="full")
        df = res.to_pandas() if hasattr(res, "to_pandas") else pd.DataFrame(res)
        df["predictor"] = mode; df["source"] = a.source; df["target_cell_line"] = a.target
        df["seed"] = a.seed
        df["target_gene_dropped"] = a.drop_target_gene or a.exclude_target_genespace
        allres.append(df)
        if de_real_path is None:
            cand = os.path.join(od, f"{mode}_real_de.csv")
            if os.path.exists(cand):
                de_real_path = cand
        print(f"  {mode}: {len(df)} perts scored ({time.time()-t0:.0f}s)", flush=True)
        del P, Xr, real, pred

    out = pd.concat(allres, ignore_index=True)
    f = os.path.join(OUT, f"task2_perpert_{tag}.csv")
    out.to_csv(f, index=False)
    num = [c for c in out.select_dtypes("number").columns if c != "seed"]
    print(out.groupby("predictor")[num].mean().T.to_string())
    print("wrote", f)


if __name__ == "__main__":
    main()
