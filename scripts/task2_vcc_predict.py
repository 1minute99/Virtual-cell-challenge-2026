#!/usr/bin/env python
"""Task 2 (VCC side) - produce the actual VCC 2026 validation predictions.

NOTE: Arc withholds the VCC validation ground truth (it IS the leaderboard), and
there is no 2026 challenge training set. These predictions therefore CANNOT be
scored locally by cell-eval; they are emitted in submission shape only. The
scored evidence for these predictors lives in scripts/task2_baselines.py, which
runs the identical predictors on Replogle target cell lines that do have ground
truth.

Emits, per context and per predictor, a pseudobulk prediction matrix
(300 targets x 18,533 VCC genes) plus the per-gene delta actually applied.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task2_vcc_predict.py
"""
import json, os, sys
import numpy as np, pandas as pd, h5py
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genemap import GeneMap
from source_delta import load_source_deltas

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VCC = os.path.join(ROOT, "data/2026_AIVC _validation_controls")
OUT = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)


def control_pseudobulk(ctx):
    """log1p(CP10K) pseudobulk of the VCC non-targeting control cells."""
    with h5py.File(os.path.join(VCC, f"context_{ctx}.h5ad"), "r") as h:
        g = h["X"]; shape = tuple(g.attrs["shape"])
        M = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
    tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
    L = M.astype(np.float64)
    L.data = np.log1p(L.data / np.repeat(tot, np.diff(L.indptr)) * 1e4)
    return np.asarray(L.mean(0)).ravel(), shape[0]


def main():
    gm = GeneMap()
    genes = [l.strip() for l in open(os.path.join(VCC, "gene_names.csv"))][1:]
    genes = [g for g in genes if g]
    targets = [l.strip() for l in open(os.path.join(VCC, "pert_counts.csv"))][1:]
    targets = [t for t in targets if t]

    deltas, info, idx_mapped = load_source_deltas(genes, gm)
    print("source deltas:", info)
    mapped = np.zeros(len(genes), bool); mapped[idx_mapped] = True

    covered = [t for t in targets if gm.canon(t) in deltas]
    missing = [t for t in targets if gm.canon(t) not in deltas]
    print(f"targets covered by K562 genome-scale source: {len(covered)}/{len(targets)}")
    print(f"missing ({len(missing)}): {missing}")

    mean_delta = np.mean([deltas[gm.canon(t)] for t in covered], axis=0)
    summary = []
    for ctx in ["A", "B", "C"]:
        cpb, ncells = control_pseudobulk(ctx)
        for mode in ["source_copy_K562", "return_control", "train_mean_delta"]:
            M = np.zeros((len(targets), len(genes)), dtype=np.float32)
            for i, t in enumerate(targets):
                if mode == "return_control":
                    d = 0.0
                elif mode == "train_mean_delta":
                    d = mean_delta
                else:
                    d = deltas.get(gm.canon(t), 0.0)
                M[i] = np.maximum(cpb + d, 0.0)
            np.save(os.path.join(OUT, f"vcc_pred_{ctx}_{mode}.npy"), M)
            summary.append(dict(context=ctx, predictor=mode, n_targets=len(targets),
                                n_genes=len(genes), n_control_cells=ncells,
                                n_genes_with_source_signal=int(mapped.sum()),
                                n_genes_passthrough=int((~mapped).sum()),
                                mean_abs_delta=float(np.abs(M - cpb).mean())))
        print(f"context {ctx}: control pseudobulk from {ncells} cells done")

    df = pd.DataFrame(summary)
    df.to_csv(os.path.join(OUT, "task2_vcc_prediction_summary.csv"), index=False)
    json.dump(dict(covered=covered, missing=missing, source_info=info),
              open(os.path.join(OUT, "task2_vcc_coverage.json"), "w"), indent=1)
    np.save(os.path.join(OUT, "vcc_gene_has_source_signal.npy"), mapped)
    print(df.to_string(index=False))
    print("\nNOTE: not locally scorable - VCC validation ground truth is withheld by Arc.")


if __name__ == "__main__":
    main()
