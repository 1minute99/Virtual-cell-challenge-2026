#!/usr/bin/env python3
"""Build a prediction for a held-out context, using exactly the submission's maths.

This is `build_submission.py` pointed at the harness instead of the leaderboard.
It transfers a pseudobulk log-CP10K delta onto that context's own control cells
and reconstructs integer counts the same way, so a strategy measured here is the
strategy that would be submitted -- not an approximation of it.

The aggregation flag is the point of the whole harness. `mean` is what scored
+0.0771; `unanimous` is what scored -0.0224 and must come out WORSE here or the
harness is not measuring the right thing. `weighted` is the untested proposal:
cross-source agreement multiplies a gene's magnitude instead of gating it, so an
unreliable gene sinks down the |log2FC| ranking that `reach` walks rather than
leaving the call set that `fid`'s denominator counts.

Control cells are carried through bit-exact. Under control_source="real" the
prediction's DE is computed against the REAL control pool, so these rows change
no metric -- but `io.validate_pair` requires both sides to carry the same
perturbation set, and the control label is part of it.
"""
import argparse
import json
import pathlib
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import build_baseline_delta as BBD
import vcclib as V

CONTROL = "non-targeting"
WMIN = 0.05          # a demoted gene must still be CALLED; see aggregate()
STORES = ["orion_target_cellline_all", "cd4_target_condition_all", "kolf_target_batch",
          "replogle_k562_gwps", "feng_genomewidescreen_nonfitnessgenes_target",
          "vcc2025_target_batch"]


def aggregate(stack, how, tilt=0.0, prec=None):
    """stack is (n_sources, n_targets, n_genes), NaN where a source did not measure.

    `tilt` redistributes magnitude toward genes the sources agree on, WITHOUT changing
    how much magnitude there is in total. That separation matters: the earlier
    `weighted` mode multiplied every gene by a factor <= 1, so it conflated "make the
    delta non-uniform" with "make the delta smaller" -- and shrinking has lost every
    time it was tried, which could hide whatever the reshaping was worth. Here the
    result is rescaled so sum|delta| is identical to the plain mean at tilt=0, so the
    only thing that varies is the SHAPE.

    Amplifying confident genes is also the operation the metrics ask for: `reach` walks
    a prefix ranked by predicted |log2FC|, so magnitude is a confidence ranking, and
    `fid`'s max(n_pred, n_conf) denominator punishes dropping a gene but not demoting
    one. Nothing is ever zeroed.
    """
    if how == "precision":
        # Within a source, rows are already inverse-variance weighted (build_baseline_delta.py:43-46
        # uses w = 1/(1/n_pert + 1/n_ctrl)). ACROSS sources the mean is uniform, so a source with 23
        # cells per target (feng) counts as much as one with 518 (cd4) -- a 22x difference in noise
        # variance. `prec` is each source's accumulated weight, already computed and previously
        # discarded, so this is the inverse-variance combination the within-source step already does.
        ok = np.isfinite(stack) & (prec > 0)
        w = np.where(ok, prec, 0.0)
        num = (np.where(ok, stack, 0.0) * w).sum(0)
        den = w.sum(0)
        return np.nan_to_num(np.where(den > 0, num / np.maximum(den, 1e-12), 0.0))
    ok = np.isfinite(stack)
    n = ok.sum(0)
    vals = np.where(ok, stack, 0.0)
    with np.errstate(invalid="ignore"):
        mean = np.where(n > 0, vals.sum(0) / np.maximum(n, 1), 0.0)
    mean = np.nan_to_num(mean)
    pos_, neg_ = (vals > 0).sum(0), (vals < 0).sum(0)
    if how == "mean":
        if tilt <= 0:
            return mean
        agree = np.maximum(pos_, neg_) / np.maximum(n, 1)
        w = np.where(n >= 2, (1.0 - tilt) + tilt * np.clip(2.0 * agree - 1.0, 0.0, 1.0), 1.0)
        w = np.maximum(w, WMIN)
        tot = np.abs(mean).sum()
        out = mean * w
        got = np.abs(out).sum()
        return out * (tot / got) if got > 0 else out
    pos, neg = pos_, neg_
    if how == "unanimous":
        keep = (n == 1) | ((n >= 2) & ((pos == n) | (neg == n)))
        return np.where(keep, mean, 0.0)
    if how == "weighted":
        # agreement fraction in [0.5, 1]: 1.0 when every source that measured the gene
        # took the same side, 0.5 on an even split. Rescaled to [0, 1] so a coin-flip
        # gene keeps none of its magnitude and a unanimous one keeps all of it. Genes
        # with a single source are left alone -- there is no disagreement to measure.
        agree = np.maximum(pos, neg) / np.maximum(n, 1)
        w = np.where(n >= 2, np.clip(2.0 * agree - 1.0, 0.0, 1.0), 1.0)
        return mean * w
    raise ValueError(how)


def reconstructed_blocks(M, D0, cells=400, seed=0, rounding="stochastic"):
    """Yield original reconstruction blocks for full-cell and exact bulk checks.

    This is the original main() arithmetic and RNG order, extracted unchanged.
    Each tuple is (target index, changed-gene counts, unchanged COO, S, NS).
    Consumers may sum blocks without storing all predicted cells on disk.
    """
    L = np.asarray(M.sum(1)).ravel().astype(np.float64)
    if L.min() <= 0:
        raise ValueError("a control cell has zero counts")
    S = np.nonzero(np.abs(D0).max(0) > 0)[0]
    NS = np.setdiff1d(np.arange(M.shape[1]), S, assume_unique=True)
    expd = np.exp(D0[:, S].astype(np.float64))
    rng = np.random.default_rng(seed)
    pick = rng.choice(M.shape[0], size=cells, replace=False)
    blk = M[pick]
    Lb = L[pick][:, None]
    dense = np.asarray(blk[:, S].todense(), dtype=np.float64)
    aq = blk[:, NS].tocoo()
    for ti in range(len(D0)):
        want = ((dense / Lb * 1e4 + 1.0) * expd[ti] - 1.0) * Lb / 1e4
        np.clip(want, 0, None, out=want)
        if rounding == "stochastic":
            floor = np.floor(want)
            newS = floor + (rng.random(want.shape) < (want - floor))
        elif rounding == "nearest":
            newS = np.rint(want)
        else:
            raise ValueError(rounding)
        yield ti, newS, aq, S, NS


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--real", required=True, help="the harness real h5ad")
    ap.add_argument("--out", required=True)
    ap.add_argument("--agg", choices=["mean", "unanimous", "weighted", "precision"], default="mean")
    ap.add_argument("--rounding", choices=["stochastic", "nearest"], default="stochastic")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--tilt", type=float, default=0.0,
                    help="0 = every gene scaled equally (what has been submitted so far); "
                         "1 = magnitude fully redistributed by cross-source sign agreement, "
                         "total magnitude held fixed. Only valid with --agg mean.")
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--stores", default=",".join(STORES))
    ap.add_argument("--delta-file", default="",
                    help="skip the cross-source aggregation and use this (n_targets, 18533) "
                         ".npy in OFFICIAL gene order, rows in the harness's target order "
                         "(sorted categories minus the control). This is how a LEARNED delta "
                         "is scored on exactly the same count reconstruction as the "
                         "mean-transfer baseline, so the comparison isolates the delta.")
    ap.add_argument("--source-scale", default="",
                    help="comma-separated STORE=FACTOR. Sources are averaged with EQUAL "
                         "weight, but they are not on a common scale: each normalises "
                         "CP10K over the official features IT measures, and Replogle "
                         "measures 7,681 where the others measure ~18,000. Measured on the "
                         "same Replogle data, narrowing the denominator from 8,248 genes to "
                         "7,683 moves the own-gene delta from -0.2414 to -0.3318, a factor "
                         "of 1.37 -- so its deltas enter the average inflated. Dividing by "
                         "that measured factor is a correction, not a tuned weight.")
    ap.add_argument("--drop-store", default="",
                    help="comma-separated stores to EXCLUDE. Required when the held-out "
                         "context is itself one of the six sources: leaving KOLF in while "
                         "predicting KOLF reads the answer off the reference and the "
                         "harness stops measuring transfer at all.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--save-delta", default="",
                    help="also write the aggregated delta BEFORE --scale as (n_targets, 18533) "
                         ".npy in official gene order (zeros for genes the harness drops), so "
                         "the mean-transfer aggregate can be re-used or combined via --delta-file")
    ap.add_argument("--save-delta-only", action="store_true",
                    help="stop after --save-delta (skip building cells)")
    a = ap.parse_args()

    A = ad.read_h5ad(a.real)
    genes = list(A.var_names)
    targets = [t for t in A.obs["target"].cat.categories if t != CONTROL]
    ctl = A[A.obs["target"] == CONTROL]
    M = sparse.csr_matrix(ctl.X)
    print(f"real: {A.shape}, {len(targets)} targets, {M.shape[0]} control cells", flush=True)

    F = list(V.official_features())
    col = {g: i for i, g in enumerate(F)}
    take = np.array([col[g] for g in genes])          # official order -> harness order
    gm = V.GeneMap()
    drop = {x for x in a.drop_store.split(",") if x}
    use = [k for k in a.stores.split(",") if k not in drop]
    if drop:
        missing = drop - set(a.stores.split(","))
        if missing:
            raise SystemExit(f"--drop-store names stores that are not in --stores: {missing}")
        print(f"held out of the prediction: {sorted(drop)}", flush=True)
    sscale = {}
    for item in (x for x in a.source_scale.split(",") if x):
        k, v = item.split("=")
        if k not in use:
            raise SystemExit(f"--source-scale names {k!r}, which is not in use: {use}")
        sscale[k] = float(v)
    P, PREC = [], []
    for k in ([] if a.delta_file else use):
        prof, den, *_ = BBD.source_profile(k, targets, gm)
        prof = prof[:, take]; PREC.append(den[:, take])
        if k in sscale:
            prof = prof * sscale[k]
            print(f"  {k:48s} rescaled x{sscale[k]}", flush=True)
        P.append(prof)
        print(f"  {k:48s} measured {np.isfinite(P[-1]).mean()*100:5.1f}%", flush=True)  # noqa
    if a.delta_file:
        Dfull = np.load(a.delta_file).astype(np.float32)
        if Dfull.shape != (len(targets), len(F)):
            raise SystemExit(f"--delta-file has shape {Dfull.shape}, expected "
                             f"({len(targets)}, {len(F)}) in official gene order")
        if not np.all(np.isfinite(Dfull)):
            raise SystemExit("--delta-file contains non-finite values")
        D0 = Dfull[:, take] * a.scale
        print(f"delta from {a.delta_file} (aggregation skipped)", flush=True)
    else:
        D0 = aggregate(np.stack(P), a.agg, tilt=a.tilt,
                       prec=(np.stack(PREC) if a.agg == "precision" else None)) * a.scale

    if a.save_delta:
        full = np.zeros((len(targets), len(F)), np.float32)
        full[:, take] = (D0 / a.scale).astype(np.float32)
        np.save(a.save_delta, full)
        print(f"saved aggregated delta (pre-scale) -> {a.save_delta} {full.shape}", flush=True)
        if a.save_delta_only:
            return
    nz = np.abs(D0) > 0
    print(f"\ndelta: {a.agg} x{a.scale} tilt={a.tilt}  {nz.mean()*100:.1f}% of cells move, "
          f"median |d| on moved {np.median(np.abs(D0[nz])):.4f}", flush=True)

    rows, cols, vals = [], [], []
    for ti, newS, aq, S, NS in reconstructed_blocks(M, D0, a.cells, a.seed, a.rounding):
        r2, c2 = np.nonzero(newS)
        off = ti * a.cells
        rows.append(np.concatenate([aq.row, r2]) + off)
        cols.append(np.concatenate([NS[aq.col], S[c2]]))
        vals.append(np.concatenate([aq.data.astype(np.float32),
                                    newS[r2, c2].astype(np.float32)]))
        if ti % 50 == 0:
            print(f"  target {ti+1}/{len(targets)}", flush=True)

    X = sparse.csr_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(len(targets) * a.cells, len(genes)))
    X.data = np.rint(X.data).astype(np.float32)
    X = sparse.vstack([X, M], format="csr")          # controls, unmodified
    lab = np.concatenate([np.repeat(targets, a.cells),
                          np.full(M.shape[0], CONTROL)])
    out = ad.AnnData(X=X,
                     obs=pd.DataFrame({"target": pd.Categorical(lab)},
                                      index=[f"p{i}" for i in range(X.shape[0])]),
                     var=pd.DataFrame(index=pd.Index(genes, name=None)))
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(a.out, compression="gzip")
    meta = dict(agg=a.agg, rounding=a.rounding, scale=a.scale, tilt=a.tilt, cells=a.cells,
                delta_file=a.delta_file or None,
                seed=a.seed, stores=use, dropped=sorted(drop), source_scale=sscale,
                shape=list(out.shape),
                frac_cells_moved=float(nz.mean()),
                median_abs_delta_moved=float(np.median(np.abs(D0[nz]))) if nz.any() else 0.0)
    pathlib.Path(a.out).with_suffix(".meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\n-> {a.out}  {out.shape}")


if __name__ == "__main__":
    main()
