#!/usr/bin/env python3
"""Nested whole-cell-line validation of learned cross-context effect operators.

Uses integer pseudobulk counts, replicate-matched controls, a balanced target panel,
and panel-target exclusion. No response from the outer line enters fitting or model
selection. Models are supervised kernel ridge regressions on leave-one-line-out
mean effects, with an optional RBF kernel on unperturbed expression. They predict
signed vectors, rather than positive multipliers of a baseline. All hyperparameters
are selected on inner whole-line folds. The fast rank calculation is checked against
the installed cell-eval2 implementation; these are pseudobulk PDS, NOT the full
competition's normalized six-member score or cell reconstruction pipeline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import sys
import time
from dataclasses import dataclass

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

ROOT = pathlib.Path(__file__).resolve().parents[1]


def center(x):
    return x - x.mean(axis=-2, keepdims=True)


def logcp(counts):
    return np.log1p(counts / counts.sum(axis=-1, keepdims=True) * 1e4)


def derangements(n, seed=0, count=32):
    rng = np.random.default_rng(seed)
    ans = []
    while len(ans) < count:
        p = rng.permutation(n)
        if np.all(p != np.arange(n)):
            ans.append(p)
    return np.array(ans)


def ranks_from_similarity(s):
    """cell-eval2 midrank with denominator n-1 (higher similarity is closer)."""
    n = len(s)
    diag = np.diag(s)[:, None]
    ranks = (s > diag).sum(1) + ((s == diag).sum(1) - 1) / 2
    return 1 - ranks / (n - 1)


def pds(pred, truth):
    pn = np.linalg.norm(pred, axis=1)
    rn = np.linalg.norm(truth, axis=1)
    sim = np.divide(pred @ truth.T, pn[:, None] * rn[None, :],
                    out=np.zeros((len(pred), len(truth))),
                    where=(pn[:, None] * rn[None, :]) > 0)
    return ranks_from_similarity(sim)


def verify_scorer():
    # scripts/coverage.py is a data inventory script, not the optional coverage
    # package that numba imports. Avoid that name collision without changing the env.
    previous = sys.path[:]
    script_dir = pathlib.Path(__file__).resolve().parent
    sys.path[:] = [p for p in sys.path if pathlib.Path(p).resolve() != script_dir]
    try:
        from cell_eval2.metrics.discrimination import discrimination_score
    finally:
        sys.path[:] = previous
    rng = np.random.default_rng(9)
    for kind in ["random", "zero", "constant", "ties"]:
        a, b = rng.normal(size=(8, 11)), rng.normal(size=(8, 11))
        if kind == "zero":
            a[:] = 0
        if kind == "constant":
            a[:] = a[0]
        if kind == "ties":
            b[1] = b[0]
        labels = np.array(["non-targeting"] + [f"p{i}" for i in range(8)])
        # Removal has already happened when the benchmark calls pds.
        z = np.zeros((1, a.shape[1]))
        official = discrimination_score(pred_bulk=(labels, np.r_[z, a]),
                real_bulk=(labels, np.r_[z, b]), distance="cosine",
                exclude_target_gene=False, rank_denominator="n-1")
        np.testing.assert_allclose(pds(a, b), list(official.values()), atol=1e-12)


@dataclass
class Panel:
    stimulus: str
    lines: list[str]
    targets: list[str]
    genes: list[str]
    effects: np.ndarray   # line, replicate, target, gene
    controls: np.ndarray  # line, gene (unperturbed mean across replicates)
    precision: np.ndarray # line, target; effective cell n, NOT gene-wise precision
    audit: dict


def load_panel(stim, min_cells, official_only=False):
    stem = ROOT / "processed/jiang" / stim
    obs = pd.read_csv(str(stem) + ".obs.csv")
    genes = pd.read_csv(str(stem) + ".var.csv").gene.astype(str).tolist()
    counts = np.fromfile(str(stem) + ".bin", dtype=np.float64).reshape(len(obs), len(genes))
    assert obs.set_index(["cell_type", "rep", "gene"]).index.is_unique
    assert len(genes) == len(set(genes))
    assert np.isfinite(counts).all() and counts.min() >= 0
    np.testing.assert_allclose(counts.sum(1), obs.total_counts, rtol=0, atol=0)
    assert np.equal(counts, np.rint(counts)).all()
    logs = logcp(counts)
    obs["row"] = np.arange(len(obs))
    controls = obs[obs.gene == "NT"]
    reps = sorted(obs.rep.unique())
    assert len(reps) == 2
    lines = [l for l in sorted(obs.cell_type.unique())
             if set(controls[controls.cell_type == l].rep) == set(reps)]
    valid = obs[(obs.n_cells >= min_cells) & (obs.gene != "NT")]
    targets = sorted(set.intersection(*[
        set(valid[(valid.cell_type == l) & (valid.rep == r)].gene)
        for l in lines for r in reps]))
    assert len(targets) > 10
    # Remove ALL screened target genes, including those failing cell-count QC.
    # Exclude from model inputs as well as scoring to avoid on-target shortcuts.
    panel_targets = set(obs.gene) - {"NT"}
    take = np.array([g not in panel_targets for g in genes])
    if official_only:
        feature_file = ROOT.parent / "data/2026_AIVC _validation_controls/gene_names.csv"
        official = set(pd.read_csv(feature_file).gene_name.astype(str))
        take &= np.array([g in official for g in genes])
    index = obs.set_index(["cell_type", "rep", "gene"])
    d = np.zeros((len(lines), 2, len(targets), take.sum()), dtype=np.float64)
    ctrl = np.zeros((len(lines), 2, take.sum()), dtype=np.float64)
    prec = np.zeros((len(lines), 2, len(targets)))
    for li, line in enumerate(lines):
        for ri, rep in enumerate(reps):
            c = index.loc[(line, rep, "NT")]
            assert c.n_cells >= min_cells
            ctrl[li, ri] = logs[int(c.row), take]
            for ti, target in enumerate(targets):
                p = index.loc[(line, rep, target)]
                d[li, ri, ti] = logs[int(p.row), take] - ctrl[li, ri]
                prec[li, ri, ti] = 1 / (1 / p.n_cells + 1 / c.n_cells)
    audit = dict(stimulus=stim, lines=lines, n_targets=len(targets),
                 n_genes=int(take.sum()), n_screened_targets=len(panel_targets),
                 excluded_target_columns=sum(g in panel_targets for g in genes),
                 official_features_only=official_only,
                 excluded_lines=sorted(set(obs.cell_type) - set(lines)),
                 excluded_low_cell_targets=sorted(panel_targets - set(targets)),
                 min_cells=min_cells, raw_count_sums_match_metadata=True,
                 normalization="log1p(CP10K), full raw feature library before panel exclusion",
                 file_sha256=hashlib.sha256(pathlib.Path(str(stem)+".bin").read_bytes()).hexdigest())
    return Panel(stim, lines, targets, list(np.array(genes)[take]), d,
                 ctrl.mean(1), prec.sum(1), audit)


def model_names():
    names = ["mean", "effective_n", "centered_mean"]
    for family in ["ridge", "context_ridge", "residual_ridge"]:
        for alpha in [0.1, 1.0, 10.0, 100.0]:
            names.append(f"{family}:{alpha:g}")
    for rank in [8, 16, 32]:
        names.append(f"pca:{rank}")
    names += [f"replicate_ridge:{alpha:g}" for alpha in [0.1, 1.0, 10.0]]
    names += [x + ":specific" for x in names if ":" in x]
    return names


def fit_predict(panel, train, test):
    """Only accesses training perturbations; the test context supplies controls."""
    n = len(panel.targets)
    t = len(train)
    assert test not in train and t >= 2
    src = panel.effects[train].mean(1)
    means = src.mean(1)
    common = means.mean(0)
    a = (src - means[:, None, :]).reshape(t * n, -1)
    b = np.tile(np.eye(n) / t, (1, t))
    q = np.kron((np.ones((t, t)) - np.eye(t)) / (t - 1), np.eye(n))
    gram = a @ a.T
    kg = q @ gram @ q.T
    test_kernel = b @ gram @ q.T
    w = panel.precision[train]
    out = {"mean": src.mean(0), "centered_mean": b @ a,
           "effective_n": (src * w[:, :, None]).sum(0) / w.sum(0)[:, None]}
    # RBF length scale is fitted on training controls only. Test control is allowed.
    ctr = panel.controls[train]
    cd = ((ctr[:, None] - ctr[None, :]) ** 2).sum(2)
    length = max(float(np.median(cd[np.triu_indices(t, 1)])), 1e-12)
    context_kernel = np.exp(-cd / length)
    context_query = np.exp(-((ctr - panel.controls[test]) ** 2).sum(1) / length)
    for family in ["ridge", "context_ridge", "residual_ridge"]:
        k = kg.copy()
        kt = test_kernel.copy()
        if family == "context_ridge":
            k *= np.repeat(np.repeat(1 + context_kernel, n, axis=0), n, axis=1)
            kt *= np.repeat(1 + context_query, n)[None, :]
        scale = max(float(np.trace(k) / len(k)), 1e-12)
        ev, u = eigh((k + k.T) / 2, check_finite=False)
        ev = np.maximum(ev, 0)
        proj = kt @ u
        response = (np.eye(len(q)) - q) @ a if family == "residual_ridge" else a
        decoded = u.T @ response
        for alpha in [0.1, 1.0, 10.0, 100.0]:
            pred = (proj / (ev + alpha * scale)) @ decoded
            if family == "residual_ridge":
                pred += b @ a
            # No test response centroid: restore the SOURCE common response only.
            out[f"{family}:{alpha:g}"] = pred + common
    ev, u = eigh((gram + gram.T) / 2, check_finite=False)
    for rank in [8, 16, 32]:
        order = np.argsort(ev)[::-1][:rank]
        order = order[ev[order] > 1e-10]
        # Projection of query into a basis learned from source perturbations only.
        coeff = (b @ gram @ u[:, order]) / ev[order]
        out[f"pca:{rank}"] = (coeff @ u[:, order].T) @ a + common
    # Supervised denoising: source replicate 1 predicts source replicate 2 and
    # conversely, centered separately to remove their shared control offsets.
    reps = panel.effects[train]
    reps = reps - reps.mean(axis=2, keepdims=True)
    xrep = reps.reshape(t * 2 * n, -1)
    yrep = reps[:, ::-1].reshape(t * 2 * n, -1)
    krep = xrep @ xrep.T
    ev, u = eigh((krep + krep.T) / 2, check_finite=False)
    ev = np.maximum(ev, 0)
    scale = max(float(np.trace(krep) / len(krep)), 1e-12)
    proj = ((b @ a) @ xrep.T) @ u
    decoded = u.T @ yrep
    for alpha in [0.1, 1.0, 10.0]:
        out[f"replicate_ridge:{alpha:g}"] = (proj / (ev + alpha * scale)) @ decoded + common
    # Explicit ablation of the source-common response, identically for every family.
    for name, pred in list(out.items()):
        if ":" in name:
            out[name + ":specific"] = pred - common
    return out


def metrics(pred, truth_reps, perms):
    truth = truth_reps.mean(0)
    pn = np.linalg.norm(pred, axis=1)
    rn = np.linalg.norm(truth, axis=1)
    sim = np.divide(pred @ truth.T, pn[:, None] * rn[None, :],
                    out=np.zeros((len(pred), len(truth))),
                    where=(pn[:, None] * rn[None, :]) > 0)
    row = ranks_from_similarity(sim)
    null = np.array([ranks_from_similarity(sim[p]).mean() for p in perms])
    pc = center(pred)
    r1, r2 = center(truth_reps[0]), center(truth_reps[1])
    # Aggregate signed replicate cross-product; never discard negative estimates.
    signal = float(np.sum(r1 * r2))
    pred_energy = float(np.sum(pc ** 2))
    numerator = float(np.sum(pc * (r1 + r2) / 2))
    corrected = numerator / np.sqrt(pred_energy * signal) if signal > 0 and pred_energy > 0 else np.nan
    return dict(pds=float(row.mean()), shuffle=float(null.mean()),
                shuffle_seed0=float(null[0]), shuffle_sd=float(null.std()),
                pds_rep1=float(pds(pred, truth_reps[0]).mean()),
                pds_rep2=float(pds(pred, truth_reps[1]).mean()),
                centered_pds=float(pds(pc, center(truth)).mean()),
                row_cosine=float(np.diag(sim).mean()),
                centered_noise_corrected_cosine=corrected,
                centered_signal_energy=signal,
                common_energy_fraction=float(len(pred)*np.sum(pred.mean(0)**2)/max(np.sum(pred**2),1e-20))), row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stimuli", nargs="+", default=["IFNB", "INS", "TGFB", "TNFA"])
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "reports/jiang_benchmark_2026-09-14")
    ap.add_argument("--min-cells", type=int, default=25)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-nested", action="store_true", help="Exploration only; do not report selected performance")
    ap.add_argument("--official-genes", action="store_true", help="Restrict model/evaluation features to the challenge panel")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    verify_scorer()
    print("PDS parity with installed cell-eval2 passed (random, zero, constant, ties).", flush=True)
    all_rows, inner_rows, selection, target_rows, reliability = [], [], [], [], []
    start = time.monotonic()
    with threadpool_limits(args.threads):
        for stim in args.stimuli:
            panel = load_panel(stim, args.min_cells, args.official_genes)
            (args.out / f"{stim}.audit.json").write_text(json.dumps(panel.audit, indent=2))
            print(json.dumps(panel.audit), flush=True)
            perms = derangements(len(panel.targets))
            cache = {}
            # With separate stimulus models, holding out a line excludes its responses
            # from every fitting and selection path across all stimuli.
            for outer in range(len(panel.lines)):
                train = [i for i in range(len(panel.lines)) if i != outer]
                scores = {name: [] for name in model_names()}
                if not args.no_nested:
                    for inner in train:
                        subset = tuple(i for i in train if i != inner)
                        key = (subset, inner)
                        if key not in cache:
                            fitted = fit_predict(panel, list(subset), inner)
                            ref = panel.effects[inner].mean(0)
                            cache[key] = {k: float(pds(v, ref).mean()) for k, v in fitted.items()}
                        for name, val in cache[key].items():
                            scores[name].append(val)
                            inner_rows.append(dict(stimulus=stim, outer_line=panel.lines[outer],
                                  inner_line=panel.lines[inner], model=name, pds=val))
                    # Fixed ordering makes ties deterministic. Baseline is eligible.
                    chosen = max(model_names(), key=lambda k: np.mean(scores[k]))
                    learned = [x for x in model_names() if x.startswith(("ridge:", "context_ridge:", "residual_ridge:", "replicate_ridge:"))]
                    chosen_ml = max(learned, key=lambda k: np.mean(scores[k]))
                else:
                    chosen = chosen_ml = None
                fitted = fit_predict(panel, train, outer)
                baseline, _ = metrics(fitted["mean"], panel.effects[outer], perms)
                for name, pred in fitted.items():
                    m, target_metric = metrics(pred, panel.effects[outer], perms)
                    denom = baseline["pds"] - m["shuffle"]
                    m["gain_vs_mean_above_shuffle"] = (m["pds"] - m["shuffle"]) / denom if abs(denom)>1e-9 else np.nan
                    row = dict(stimulus=stim, line=panel.lines[outer], model=name,
                               n_targets=len(panel.targets), selected=(name==chosen),
                               selected_ml=(name==chosen_ml), **m)
                    all_rows.append(row)
                    for p, val in zip(panel.targets, target_metric):
                        target_rows.append(dict(stimulus=stim, line=panel.lines[outer], model=name, target=p, pds=float(val)))
                # Evaluation-only replicate ceiling; never used to choose hyperparameters.
                dc = panel.effects[outer] - panel.effects[outer].mean(1, keepdims=True)
                signal = float(np.sum(dc[0]*dc[1]))
                total = float(np.sum(dc.mean(0)**2))
                reliability.append(dict(stimulus=stim, line=panel.lines[outer],
                    signal_energy=signal, averaged_effect_energy=total,
                    averaged_reliability=signal/total,
                    replicate_pds=float((pds(panel.effects[outer,0],panel.effects[outer,1]).mean()+
                                         pds(panel.effects[outer,1],panel.effects[outer,0]).mean())/2)))
                if chosen is not None:
                    selection.append(dict(stimulus=stim, line=panel.lines[outer], model=chosen, model_ml=chosen_ml,
                                          inner_pds=float(np.mean(scores[chosen])), inner_ml_pds=float(np.mean(scores[chosen_ml]))))
                    np.savez_compressed(args.out / f"{stim}.{panel.lines[outer]}.prediction.npz",
                        targets=np.array(panel.targets), genes=np.array(panel.genes),
                        prediction=fitted[chosen].astype(np.float32),
                        prediction_ml=fitted[chosen_ml].astype(np.float32),
                        baseline=fitted["mean"].astype(np.float32), model=chosen, model_ml=chosen_ml)
                # Checkpoint completed folds, so interrupted runs retain reviewable results.
                pd.DataFrame(all_rows).to_csv(args.out / "fold_metrics.csv", index=False)
                pd.DataFrame(inner_rows).to_csv(args.out / "inner_metrics.csv", index=False)
                pd.DataFrame(selection).to_csv(args.out / "selection.csv", index=False)
                pd.DataFrame(target_rows).to_csv(args.out / "target_metrics.csv", index=False)
                pd.DataFrame(reliability).to_csv(args.out / "reliability.csv", index=False)
                print(f"{stim} {panel.lines[outer]} mean={baseline['pds']:.4f} selected={chosen} ml={chosen_ml} elapsed={time.monotonic()-start:.0f}s", flush=True)
    frame = pd.DataFrame(all_rows)
    summary = frame.groupby("model")[["pds", "shuffle", "centered_pds", "centered_noise_corrected_cosine"]].mean()
    summary.to_csv(args.out / "exploratory_model_summary.csv")
    print(summary.sort_values("pds", ascending=False).to_string(), flush=True)
    if not args.no_nested:
        selected = frame[frame.selected].copy()
        ml = frame[frame.selected_ml].copy()
        means = frame[frame.model=="mean"]
        selected["method"] = "nested_selection"
        ml["method"] = "nested_ml"
        means = means.copy(); means["method"] = "mean"
        joined = pd.concat([selected, ml, means])
        byline = joined.groupby(["line", "method"])[["pds", "shuffle", "centered_pds", "pds_rep1", "pds_rep2"]].mean()
        byline.to_csv(args.out / "nested_by_cell_line.csv")
        print("NESTED CELL-LINE RESULTS (stimuli are not independent biological folds):", flush=True)
        print(byline.to_string(), flush=True)
    (args.out / "run.json").write_text(json.dumps(dict(args={k:str(v) if isinstance(v,pathlib.Path) else v for k,v in vars(args).items()},
        elapsed_seconds=time.monotonic()-start, scorer_parity_passed=True,
        evaluation="pseudobulk PDS, not full cell-eval2 competition score",
        script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()),indent=2))


if __name__ == "__main__":
    main()
