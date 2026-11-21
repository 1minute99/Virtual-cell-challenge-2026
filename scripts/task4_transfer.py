#!/usr/bin/env python
"""Task 4a - per-perturbation transfer success for every (source, target) context pair.

transfer_success(s, t, p) = Pearson( delta_s[p], delta_t[p] ) over genes,
i.e. how well the pseudobulk effect of perturbation p measured in context s
reproduces its effect measured in context t. This is the pseudobulk analogue of
cell-eval's pearson_delta and is what Task 4 correlates context distance against.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task4_transfer.py
"""
import os, sys, itertools, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genemap import GeneMap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "work/cellline_cache")
OUT = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)
CONTROL = "non-targeting"
LINES = ["k562", "rpe1", "jurkat", "hepg2"]


def pseudobulk_deltas(cell_line, gm):
    X = np.load(os.path.join(CACHE, f"X_{cell_line}.npy"))
    p = np.asarray(np.load(os.path.join(CACHE, f"pert_{cell_line}.npy"),
                           allow_pickle=True), dtype=str)
    ctrl = X[p == CONTROL].mean(0)
    D, keys = [], []
    for q in np.unique(p):
        if q == CONTROL:
            continue
        D.append(X[p == q].mean(0) - ctrl); keys.append(gm.canon(q))
    del X
    return dict(zip(keys, np.asarray(D, dtype=np.float32))), ctrl


def rowwise_pearson(A, B):
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    num = (A * B).sum(1)
    den = np.sqrt((A ** 2).sum(1) * (B ** 2).sum(1))
    den[den == 0] = np.nan
    return num / den


def main():
    gm = GeneMap()
    genes = [str(x) for x in np.load(os.path.join(CACHE, "genes.npy"), allow_pickle=True)]
    D, C = {}, {}
    for L in LINES:
        D[L], C[L] = pseudobulk_deltas(L, gm)
        print(f"{L}: {len(D[L])} perturbation deltas")

    rows, pairs = [], []
    for s, t in itertools.product(LINES, LINES):
        common = sorted(set(D[s]) & set(D[t]))
        if not common:
            continue
        A = np.stack([D[s][p] for p in common])
        B = np.stack([D[t][p] for p in common])
        r = rowwise_pearson(A, B)
        for p, v in zip(common, r):
            rows.append(dict(source=s, target=t, perturbation=p, pearson_delta=v))
        pairs.append(dict(source=s, target=t, n_perturbations=len(common),
                          mean_pearson_delta=float(np.nanmean(r)),
                          median_pearson_delta=float(np.nanmedian(r)),
                          frac_above_0p7=float(np.nanmean(r > 0.7)),
                          same_context=(s == t)))
    per = pd.DataFrame(rows); per.to_csv(os.path.join(OUT, "task4_transfer_perpert.csv"), index=False)
    pr = pd.DataFrame(pairs); pr.to_csv(os.path.join(OUT, "task4_transfer_pairs.csv"), index=False)
    print(pr.to_string(index=False))
    np.save(os.path.join(OUT, "control_pseudobulk_replogle.npy"),
            np.stack([C[L] for L in LINES]))
    json.dump(LINES, open(os.path.join(OUT, "control_pseudobulk_lines.json"), "w"))


if __name__ == "__main__":
    main()
