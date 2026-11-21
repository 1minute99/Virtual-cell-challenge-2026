#!/usr/bin/env python
"""Task 4c - does context distance under scLong / STATE SE / PCA predict transfer?

For each representation we compute pairwise context distances (cosine and
euclidean) and correlate them against per-perturbation transfer success from
task4_transfer.py. A representation "wins" if more-similar contexts really do
transfer better, i.e. a strongly NEGATIVE distance-vs-success correlation.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task4_correlate.py
"""
import itertools, json, os
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")
REPS = ["sclong", "state_se", "pca"]


def load_rep(rep):
    f = os.path.join(OUT, f"task4_emb_{rep}.npy")
    j = os.path.join(OUT, f"task4_emb_{rep}.json")
    if not (os.path.exists(f) and os.path.exists(j)):
        return None, None
    Z = np.load(f); meta = json.load(open(j))
    return Z, meta["contexts"]


def dists(Z, ctxs):
    Zn = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12)
    cos = 1 - Zn @ Zn.T
    eu = np.linalg.norm(Z[:, None] - Z[None], axis=2)
    return {"cosine": pd.DataFrame(cos, index=ctxs, columns=ctxs),
            "euclidean": pd.DataFrame(eu, index=ctxs, columns=ctxs)}


def main():
    pairs = pd.read_csv(os.path.join(OUT, "task4_transfer_pairs.csv"))
    cross = pairs[~pairs.same_context].copy()
    cross["key"] = cross.apply(lambda r: tuple(sorted((r.source, r.target))), axis=1)
    cross = cross.groupby("key", as_index=False).agg(
        mean_pearson_delta=("mean_pearson_delta", "mean"),
        n_perturbations=("n_perturbations", "mean"))
    print(f"{len(cross)} unique cross-context pairs")

    rows, avail = [], []
    for rep in REPS:
        Z, ctxs = load_rep(rep)
        if Z is None:
            print(f"[missing] {rep} - skipped"); continue
        avail.append(rep)
        D = dists(Z, ctxs)
        pd.DataFrame(Z, index=ctxs).to_csv(os.path.join(OUT, f"task4_emb_{rep}_vectors.csv"))
        for metric, M in D.items():
            M.to_csv(os.path.join(OUT, f"task4_dist_{rep}_{metric}.csv"))
            d, s = [], []
            for _, r in cross.iterrows():
                a, b = r.key
                ka, kb = f"replogle:{a}", f"replogle:{b}"
                if ka in M.index and kb in M.index:
                    d.append(M.loc[ka, kb]); s.append(r.mean_pearson_delta)
            if len(d) >= 3:
                pr, pp = pearsonr(d, s); sr, sp = spearmanr(d, s)
                rows.append(dict(representation=rep, distance=metric, n_pairs=len(d),
                                 pearson_r=pr, pearson_p=pp, spearman_r=sr, spearman_p=sp))
    res = pd.DataFrame(rows)
    if len(res):
        res.to_csv(os.path.join(OUT, "task4_correlations.csv"), index=False)
        print(res.to_string(index=False))

    # VCC contexts: nearest Replogle source under each representation
    near = []
    for rep in avail:
        Z, ctxs = load_rep(rep)
        M = dists(Z, ctxs)["cosine"]
        for c in [x for x in ctxs if x.startswith("vcc:")]:
            r = M.loc[c, [x for x in ctxs if x.startswith("replogle:")]]
            near.append(dict(representation=rep, vcc_context=c,
                             nearest=r.idxmin(), cosine_distance=float(r.min()),
                             **{f"d_{k.split(':')[1]}": float(v) for k, v in r.items()}))
    if near:
        nd = pd.DataFrame(near)
        nd.to_csv(os.path.join(OUT, "task4_vcc_nearest_source.csv"), index=False)
        print("\nNearest Replogle context to each VCC context:")
        print(nd.to_string(index=False))


if __name__ == "__main__":
    main()
