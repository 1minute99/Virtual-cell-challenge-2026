#!/usr/bin/env python3
"""Combine effect stores into one (300, 18533) log-space delta per VCC target.

Two combination modes, which is the whole difference between the two baselines:

  single      One source. Rows for a target (its several culture conditions) are
              averaged with weight 1/(1/n_pert + 1/n_control) -- the precision
              of a two-group log ratio, not a fitted parameter.

  consensus   Several sources. Each source is first collapsed to one profile per
              target with the same precision weighting, then sources are averaged
              with EQUAL weight. Equal weighting is deliberate: a response that
              reproduces across K562 and primary CD4 T cells is better evidence
              for transfer to an unknown context than one seen only in the source
              that happens to have the most cells.

A gene with no evidence for a target gets exactly 0.0, which the submission
builder turns into bit-exact control passthrough. That is an explicit no-change
fallback, and the report says how many genes and targets it applies to. NaN
(unmeasured) is never silently treated as 0 during averaging.
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import vcclib as V

EFF = V.ROOT / "effects"
OUT = V.ROOT / "processed" / "baseline_deltas"


def source_profile(key, targets, gm):
    """Collapse one effect store to (n_targets, n_features) plus a weight sum."""
    d = EFF / key
    D = np.load(d / "delta.npy")
    P = pd.read_csv(d / "pairs.csv")
    prov = json.loads((d / "provenance.json").read_text())
    n_p = pd.to_numeric(P.n_cells_pert, errors="coerce").fillna(0).to_numpy(float)
    n_c = pd.to_numeric(P.n_cells_control, errors="coerce").fillna(0).to_numpy(float)
    with np.errstate(divide="ignore"):
        w = 1.0 / (1.0 / np.maximum(n_p, 1e-9) + 1.0 / np.maximum(n_c, 1e-9))
    sym = P.target_symbol.astype(str).to_numpy()

    # index rows by target, accepting a reviewed alias only when the exact name is absent
    by_exact, by_canon = {}, {}
    for i, s in enumerate(sym):
        by_exact.setdefault(s, []).append(i)
        by_canon.setdefault(gm.canon(s), []).append(i)

    nT, nG = len(targets), D.shape[1]
    num = np.zeros((nT, nG))
    den = np.zeros((nT, nG))
    rows_used = np.zeros(nT, dtype=np.int64)
    cells = np.zeros(nT, dtype=np.int64)
    match = ["none"] * nT
    for t_i, t in enumerate(targets):
        idx = by_exact.get(t)
        m = "exact_symbol"
        if not idx:
            idx = by_canon.get(gm.canon(t))
            m = "reviewed_alias"
        if not idx:
            continue
        match[t_i] = m
        rows_used[t_i] = len(idx)
        cells[t_i] = int(n_p[idx].sum())
        sub = D[idx]
        ok = np.isfinite(sub)
        ww = np.where(ok, w[idx][:, None], 0.0)
        num[t_i] = np.nansum(np.where(ok, sub, 0.0) * ww, axis=0)
        den[t_i] = ww.sum(axis=0)
    prof = np.full((nT, nG), np.nan)
    nz = den > 0
    prof[nz] = num[nz] / den[nz]
    return prof, den, rows_used, cells, match, prov, P


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("--stores", required=True, help="comma-separated effect store keys")
    ap.add_argument("--mode", choices=["single", "consensus", "agree"], required=True)
    ap.add_argument("--min-sources", type=int, default=2,
                    help="agree mode: sources required before a sign vote is meaningful")
    a = ap.parse_args()

    targets = V.official_targets()
    features = V.official_features()
    gm = V.GeneMap()
    keys = a.stores.split(",")
    if a.mode == "single" and len(keys) != 1:
        raise SystemExit("mode single takes exactly one store")

    profs, meta, provs = [], [], {}
    for k in keys:
        prof, den, rows, cells, match, prov, P = source_profile(k, targets, gm)
        cov = np.isfinite(prof).any(1)
        print(f"{k}: {int(cov.sum())}/300 official targets, "
              f"{int(np.isfinite(prof).any(0).sum())}/{len(features)} features with "
              f"evidence, {int(rows.sum())} effect rows used")
        profs.append(prof)
        provs[k] = prov
        meta.append(pd.DataFrame(dict(target_gene=targets, store_key=k,
                                      match=match, effect_rows=rows,
                                      cells_pert_total=cells,
                                      genes_with_evidence=np.isfinite(prof).sum(1))))

    stack = np.stack(profs)                          # (n_sources, 300, 18533)
    ok = np.isfinite(stack)
    n_src = ok.sum(0)
    with np.errstate(invalid="ignore"):
        combined = np.where(n_src > 0,
                            np.nansum(np.where(ok, stack, 0.0), axis=0) /
                            np.maximum(n_src, 1), np.nan)
    agreement = None
    if a.mode == "agree":
        # MEASURED on the submitted three-source set: where two or more sources both
        # have a value, 50.7% of (target, gene) cells DISAGREE on the sign. Averaging
        # those keeps an essentially arbitrary direction at half magnitude, and four of
        # the six official metrics score DE DIRECTION (fid, reach, jac) or DE magnitude
        # (nmae). That is the most likely reason fid is the one negative term and why it
        # got worse when sources went from one to three.
        # So: keep a cell only where every source that has it agrees on the sign.
        # A cell seen by a single source is kept as-is - there is nothing to disagree
        # with, and dropping it would throw away most of the coverage.
        vals = np.where(ok, stack, 0.0)
        pos = (vals > 0).sum(0)
        neg = (vals < 0).sum(0)
        unanimous = (pos == n_src) | (neg == n_src)
        keep = (n_src == 1) | ((n_src >= a.min_sources) & unanimous)
        multi = n_src >= a.min_sources
        agreement = dict(
            cells_with_multiple_sources=int(multi.sum()),
            sign_agree=int((multi & unanimous).sum()),
            sign_disagree=int((multi & ~unanimous).sum()),
            frac_disagree=float((multi & ~unanimous).sum() / max(1, multi.sum())),
            mass_dropped_fraction=float(
                np.abs(np.nan_to_num(combined)[multi & ~unanimous]).sum() /
                max(1e-12, np.abs(np.nan_to_num(combined)[n_src > 0]).sum())))
        print(f"  sign agreement: {agreement['sign_agree']:,} agree / "
              f"{agreement['sign_disagree']:,} disagree of "
              f"{agreement['cells_with_multiple_sources']:,} multi-source cells "
              f"({agreement['frac_disagree']*100:.1f}% disagree); dropping "
              f"{agreement['mass_dropped_fraction']*100:.1f}% of the delta mass")
        combined = np.where(keep, combined, np.nan)
    delta = np.nan_to_num(combined, nan=0.0).astype(np.float32)

    OUT.mkdir(parents=True, exist_ok=True)
    np.save(OUT / f"{a.name}_delta.npy", delta)
    M = pd.concat(meta, ignore_index=True)
    M.to_csv(OUT / f"{a.name}_target_meta.csv", index=False)

    per_target_genes = (delta != 0).sum(1)
    covered = per_target_genes > 0
    genes_movable = int((np.abs(delta).max(0) > 0).sum())
    rep = dict(
        name=a.name, mode=a.mode, stores=keys,
        n_targets=len(targets), n_features=len(features),
        targets_with_any_signal=int(covered.sum()),
        targets_with_no_signal=int((~covered).sum()),
        targets_with_no_signal_list=[t for t, c in zip(targets, covered) if not c],
        genes_that_can_move=genes_movable,
        median_genes_moved_per_target=float(np.median(per_target_genes[covered]))
        if covered.any() else 0.0,
        mean_abs_delta_on_moved=float(np.abs(delta[delta != 0]).mean()),
        max_abs_delta=float(np.abs(delta).max()),
        sources_per_covered_gene_mean=float(n_src[np.isfinite(combined)].mean())
        if np.isfinite(combined).any() else 0.0,
        sign_agreement=agreement,
        quantity="log1p(CP10K) delta vs matched control, averaged as described",
        zero_means="no evidence -> bit-exact control passthrough in the submission",
        source_provenance=provs)
    (OUT / f"{a.name}_report.json").write_text(json.dumps(rep, indent=2, default=str))
    print(f"\n{a.name}: {rep['targets_with_any_signal']}/300 targets carry a delta, "
          f"{genes_movable}/{len(features)} genes can move, "
          f"median {rep['median_genes_moved_per_target']:.0f} genes moved per target")
    print(f"  mean |delta| where nonzero {rep['mean_abs_delta_on_moved']:.4f}, "
          f"max {rep['max_abs_delta']:.3f}")
    if rep["targets_with_no_signal"]:
        print(f"  no-change fallback for: {rep['targets_with_no_signal_list']}")
    print(f"  -> {OUT/f'{a.name}_delta.npy'}")


if __name__ == "__main__":
    main()
