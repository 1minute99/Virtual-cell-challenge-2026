#!/usr/bin/env python3
"""Turn standardised pseudobulk sums into per-target control-relative effects.

The effect stored here is one specific, named quantity:

    delta[t, g] = log1p(1e4 * S_t[g] / sum_measured(S_t))
                - log1p(1e4 * S_c[g] / sum_measured(S_c))

that is, the difference of log1p CP10K profiles between a target's pseudobulk
sum and its MATCHED control's pseudobulk sum, where CP10K is computed over the
features that assay actually measured. It is not a DESeq2 log2 fold change, not
a z-score, and not a mean of log values. Unmeasured official features are NaN,
never 0, so a downstream consumer cannot mistake "not measured" for "measured
and unchanged".

Matching rule: a control must share the dataset and the biological context of
the perturbation, including stimulation and time. A resting control is never
used for a stimulated perturbation.

Outputs per store, under effects/<key>/:
    delta.npy         float32 (n_pairs, 18533), NaN where unmeasured
    pairs.csv         target, context, stimulation, n_cells, control cells, ...
    measured_mask.npy bool (18533,)
    provenance.json
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pbstore as PBS
import vcclib as V

EFF = V.ROOT / "effects"


def log_cp10k(row_sums, measured_cols):
    """log1p CP10K of one pseudobulk row, normalised over measured features only."""
    tot = row_sums[measured_cols].sum()
    if tot <= 0:
        return None
    out = np.full(row_sums.shape[0], np.nan, dtype=np.float64)
    out[measured_cols] = np.log1p(row_sums[measured_cols] / tot * 1e4)
    return out


def build(key, control_match, min_cells=10, min_control_cells=50):
    S, obs, var, prov, o2s = PBS.load_store(key)
    features = V.official_features()
    measured = o2s >= 0
    src_cols = o2s[measured]
    print(f"{key}: {S.shape[0]} groups, {int(measured.sum())}/{len(features)} official "
          f"features measured")

    obs["is_control"] = obs["is_control"].astype(str).str.lower().isin(["true", "1"])
    ctrl = obs[obs.is_control]
    pert = obs[~obs.is_control]
    if not len(ctrl):
        raise SystemExit(f"{key}: no control groups; cannot form control-relative effects")

    # index controls by their matching key
    ctrl_idx = {}
    for _, r in ctrl.iterrows():
        ctrl_idx.setdefault(tuple(str(r[c]) for c in control_match), []).append(int(r.pb_id))
    print(f"  {len(ctrl)} control groups -> {len(ctrl_idx)} matching keys on {control_match}")

    Sd = S.tocsr()
    rows, deltas, skipped = [], [], []
    for _, r in pert.iterrows():
        k = tuple(str(r[c]) for c in control_match)
        cids = ctrl_idx.get(k)
        if not cids:
            skipped.append((r.target_symbol, k, "no matched control"))
            continue
        n_p, n_c = int(r.n_cells), int(obs.loc[obs.pb_id.isin(cids), "n_cells"].sum())
        if n_p < min_cells:
            skipped.append((r.target_symbol, k, f"pert cells {n_p} < {min_cells}"))
            continue
        if n_c < min_control_cells:
            skipped.append((r.target_symbol, k, f"control cells {n_c} < {min_control_cells}"))
            continue
        sp = np.asarray(Sd[int(r.pb_id)].todense()).ravel().astype(np.float64)
        sc = np.asarray(Sd[cids].sum(0)).ravel().astype(np.float64)
        lp = log_cp10k(sp, src_cols)
        lc = log_cp10k(sc, src_cols)
        if lp is None or lc is None:
            skipped.append((r.target_symbol, k, "zero measured counts"))
            continue
        d = np.full(len(features), np.nan, dtype=np.float32)
        d[measured] = (lp[src_cols] - lc[src_cols]).astype(np.float32)
        deltas.append(d)
        rows.append(dict(target_symbol=r.target_symbol,
                         target_ensembl=r.get("target_ensembl", ""),
                         dataset_id=r.dataset_id, experiment_id=r.experiment_id,
                         cell_context=r.cell_context, stimulation=r.stimulation,
                         time=r.time, donor=r.donor, technical_batch=r.technical_batch,
                         perturbation_mode=r.perturbation_mode,
                         n_cells_pert=n_p, n_cells_control=n_c,
                         n_control_groups=len(cids),
                         n_source_groups=r.get("n_source_groups", ""),
                         qc_flags=r.get("qc_flags", ""), store_key=key))
    if not deltas:
        raise SystemExit(f"{key}: no usable target/control pairs")
    D = np.vstack(deltas)
    P = pd.DataFrame(rows)
    d = EFF / key
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "delta.npy", D)
    np.save(d / "measured_mask.npy", measured)
    P.to_csv(d / "pairs.csv", index=False)
    sk = pd.DataFrame(skipped, columns=["target_symbol", "match_key", "reason"])
    sk.to_csv(d / "skipped.csv", index=False)
    (d / "provenance.json").write_text(json.dumps(dict(
        store_key=key, source_provenance=prov,
        quantity="log1p(CP10K) difference vs matched control pseudobulk",
        not_this=["DESeq2 log2FoldChange", "z-score", "mean of log values",
                  "raw or normalised counts"],
        cp10k_normalisation="over measured features only",
        control_match_keys=control_match, min_cells=min_cells,
        min_control_cells=min_control_cells,
        n_pairs=int(D.shape[0]), n_features=int(D.shape[1]),
        official_features_measured=int(measured.sum()),
        unmeasured_encoding="NaN (never 0)",
        n_skipped=int(len(sk))), indent=2, default=str))
    print(f"  {D.shape[0]} target/control pairs, {P.target_symbol.nunique()} unique targets, "
          f"{len(sk)} skipped -> {d}")
    if len(sk):
        print(f"  skip reasons: {sk.reason.str.replace(r'\\d+', 'N', regex=True).value_counts().to_dict()}")
    return D, P


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--key", required=True)
    ap.add_argument("--match", required=True,
                    help="comma-separated obs columns a control must share")
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--min-control-cells", type=int, default=50)
    a = ap.parse_args()
    build(a.key, a.match.split(","), a.min_cells, a.min_control_cells)


if __name__ == "__main__":
    main()
