#!/usr/bin/env python3
"""The standardised pseudobulk store, and the feature mapping into the official axis.

One directory per (dataset, experiment) under pseudobulk/<key>/ holding:

  sum_counts.npz    CSR int64 - RAW UMI SUMS per group. Never a mean, never
                    normalised, never a fold change.
  obs.csv           one row per group, with the standard grouping columns and
                    n_cells. Missing information stays empty; no donor or
                    replicate is ever invented.
  var.csv           the source feature axis plus its position in the official
                    18,533-gene order (-1 when the feature is not in it).
  provenance.json   input file, hash, matrix path, count semantics, group keys,
                    and what was actually verified.

Separately from sums, `mean_raw_counts` is derivable as sum/n_cells and is
labelled as such; `log1p`, CP10K, LFC and z-scores are all different objects and
are produced only in the effects step, never stored here as if they were counts.
"""
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import vcclib as V

ROOT = V.ROOT
PB = ROOT / "pseudobulk"

OBS_COLUMNS = ["pb_id", "dataset_id", "experiment_id", "cell_context", "donor",
               "biological_replicate", "technical_batch", "stimulation", "time",
               "perturbation_mode", "target_symbol", "target_ensembl", "guide_id",
               "control_type", "is_control", "n_cells", "total_counts", "qc_flags"]


def now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def map_features(source_symbols, source_ensembl=None, gm=None, features=None):
    """Map a source feature axis onto the official 18,533-gene order.

    Resolution order per official feature: exact symbol, then Ensembl id, then a
    reviewed unique HGNC alias. A source column is claimed by at most one
    official feature, and exact matches win -- which matters because HGNC lists
    TIAF1 as an alias of MYO18A and both names are separate official features.

    Returns (var_df, official_to_source) where official_to_source[i] is the
    source column for official feature i, or -1 when that feature is unmeasured.
    """
    gm = gm or V.GeneMap()
    features = features or V.official_features()
    src_sym = [str(s).strip() for s in source_symbols]
    n_src = len(src_sym)
    by_sym = {}
    for j, s in enumerate(src_sym):
        by_sym.setdefault(s, j)
    by_ens = {}
    if source_ensembl is not None:
        for j, e in enumerate(source_ensembl):
            e = str(e).strip().split(".")[0]
            if e.startswith("ENSG"):
                by_ens.setdefault(e, j)
    by_canon = {}
    for j, s in enumerate(src_sym):
        by_canon.setdefault(gm.canon(s), j)

    o2s = np.full(len(features), -1, dtype=np.int64)
    method = ["unmatched"] * len(features)
    claimed = set()
    for pass_no in (0, 1, 2):
        for i, g in enumerate(features):
            if o2s[i] >= 0:
                continue
            j = None
            if pass_no == 0:
                j = by_sym.get(g)
                m = "exact_symbol"
            elif pass_no == 1:
                e = gm.ens(g)
                j = by_ens.get(e) if e else None
                m = "ensembl_id"
            else:
                j = by_canon.get(gm.canon(g))
                m = "reviewed_alias"
            if j is not None and j not in claimed:
                o2s[i] = j
                claimed.add(j)
                method[i] = m
    var = pd.DataFrame(dict(
        feature_index=np.arange(n_src),
        source_symbol=src_sym,
        source_ensembl=(["" if source_ensembl is None else str(e)
                         for e in (source_ensembl if source_ensembl is not None
                                   else [""] * n_src)]),
        canonical_symbol=[gm.canon(s) for s in src_sym]))
    s2o = np.full(n_src, -1, dtype=np.int64)
    for i, j in enumerate(o2s):
        if j >= 0:
            s2o[j] = i
    var["official_index"] = s2o
    return var, o2s, method


def write_store(key, sums, obs, var, provenance, official_to_source=None,
                match_method=None):
    d = PB / key
    d.mkdir(parents=True, exist_ok=True)
    if not sparse.issparse(sums):
        sums = sparse.csr_matrix(sums)
    if not np.all(np.isfinite(sums.data)):
        raise ValueError(f"{key}: non-finite values in pseudobulk sums")
    if np.any(sums.data < 0):
        raise ValueError(f"{key}: negative values in pseudobulk sums")
    nonint = np.abs(sums.data - np.rint(sums.data)) > 1e-6
    if nonint.any():
        raise ValueError(f"{key}: {int(nonint.sum())} non-integral values in what is "
                         f"recorded as a sum of raw counts")
    sums = sums.astype(np.int64)
    sparse.save_npz(d / "sum_counts.npz", sums)
    for c in OBS_COLUMNS:
        if c not in obs.columns:
            obs[c] = ""
    obs = obs[OBS_COLUMNS + [c for c in obs.columns if c not in OBS_COLUMNS]]
    obs.to_csv(d / "obs.csv", index=False)
    var.to_csv(d / "var.csv", index=False)
    if official_to_source is not None:
        o2s = np.asarray(official_to_source)
        if o2s.shape[0] != len(V.official_features()):
            raise ValueError(f"{key}: official_to_source has {o2s.shape[0]} entries, "
                             f"expected {len(V.official_features())}")
        # The mapping indexes the STORED matrix. A store whose matrix is already in
        # official gene order must carry an identity mapping, not a source-axis one.
        bad = o2s[o2s >= 0]
        if bad.size and int(bad.max()) >= sums.shape[1]:
            raise ValueError(
                f"{key}: official_to_source points at column {int(bad.max())} but the "
                f"stored matrix has only {sums.shape[1]} columns. If this matrix is "
                f"already in official gene order, pass an identity mapping "
                f"(i for measured features, -1 otherwise) instead of the source axis.")
        np.save(d / "official_to_source.npy", o2s)
        measured = o2s >= 0
        np.save(d / "official_measured_mask.npy", measured)
        provenance = dict(provenance)
        provenance["official_features_measured"] = int(measured.sum())
        provenance["official_features_unmeasured"] = int((~measured).sum())
        if match_method is not None:
            provenance["match_method_counts"] = (
                pd.Series(match_method).value_counts().to_dict())
    provenance = dict(provenance)
    provenance.update(written_utc=now(), n_groups=int(sums.shape[0]),
                      n_source_features=int(sums.shape[1]),
                      total_counts=int(sums.sum()),
                      aggregation="sum of raw RNA UMI counts",
                      mean_raw_counts="derive as sum_counts / n_cells; NOT log or LFC")
    (d / "provenance.json").write_text(json.dumps(provenance, indent=2, default=str))
    print(f"  store {key}: {sums.shape[0]} groups x {sums.shape[1]} features, "
          f"{sums.nnz} nnz, {sums.sum():,} total counts -> {d}")
    return d


def load_store(key):
    d = PB / key
    return (sparse.load_npz(d / "sum_counts.npz"),
            pd.read_csv(d / "obs.csv"),
            pd.read_csv(d / "var.csv"),
            json.loads((d / "provenance.json").read_text()),
            np.load(d / "official_to_source.npy") if (d / "official_to_source.npy").exists() else None)
