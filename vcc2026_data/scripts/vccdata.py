#!/usr/bin/env python3
"""Reader/adapter the training code uses. Four things it will not do for you.

1. It never returns a value whose meaning is ambiguous. Every accessor says in
   its name and its returned `semantics` field whether you are holding raw
   counts, a pseudobulk SUM, a pseudobulk MEAN, or a control-relative log delta.
2. It never fills an unmeasured feature with 0. `official_view` returns the
   matrix in official gene order PLUS a boolean `measured` mask, and the loss
   must be masked with it. A gene that the assay did not measure and a gene that
   was measured as zero are different facts.
3. It never invents single cells for a pseudobulk-only source. `PBSource` has no
   `cells()` method at all.
4. It never lets two views of the same biological sample land on opposite sides
   of a split. `group_key()` returns the dedup group from the dataset registry,
   which is what a split must be taken on.
"""
import json
import pathlib
import sys
from dataclasses import dataclass, field

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import vcclib as V

ROOT = V.ROOT


# ----------------------------------------------------------------- official space
class OfficialSpace:
    """The 300 targets and 18,533 output features, with the submission order fixed."""

    def __init__(self):
        self.targets = V.official_targets()
        self.features = V.official_features()
        self.target_index = {t: i for i, t in enumerate(self.targets)}
        self.feature_index = {g: i for i, g in enumerate(self.features)}
        self.manifest = json.loads((V.OFFICIAL_DIR / "manifest.json").read_text())

    def control_cells(self, context):
        """Real official control cells for one context: RAW INTEGER UMI COUNTS, CSR."""
        p = V.OFFICIAL_DIR / f"context_{context}.h5ad"
        with h5py.File(p, "r") as h:
            g = h["X"]
            M = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                                  shape=tuple(int(v) for v in g.attrs["shape"]))
            ntc = h5_col(h, "obs/ntc_id")
        return dict(X=M, semantics="raw integer UMI counts", ntc_id=ntc,
                    context=context, n_cells=M.shape[0], path=str(p))


def h5_col(h, path):
    x = h[path]
    if isinstance(x, h5py.Group):
        if "categories" in x:
            cats = np.array([c.decode() if isinstance(c, bytes) else c
                             for c in x["categories"][:]], dtype=object)
            codes = x["codes"][:]
            out = np.empty(len(codes), dtype=object)
            ok = codes >= 0
            out[ok] = cats[codes[ok]]
            return out
        v = x["values"][:]
        return np.array([t.decode() if isinstance(t, bytes) else t for t in v], dtype=object)
    v = x[:]
    return np.array([t.decode() if isinstance(t, bytes) else t for t in v], dtype=object) \
        if v.dtype.kind in "SO" else v


# ----------------------------------------------------------------- pseudobulk source
@dataclass
class PBSource:
    """One standardised pseudobulk store. Sums only; no single cells exist here."""
    key: str
    sums: sparse.csr_matrix = field(repr=False)
    obs: pd.DataFrame = field(repr=False)
    var: pd.DataFrame = field(repr=False)
    provenance: dict = field(repr=False)
    official_to_source: np.ndarray = field(repr=False)

    @property
    def semantics(self):
        return self.provenance.get("count_semantics", "unknown")

    @property
    def measured(self):
        """Boolean (18533,): which official features this assay measured."""
        return self.official_to_source >= 0

    def sum_counts(self, rows=None):
        """Pseudobulk SUM of raw counts on the SOURCE feature axis."""
        return (self.sums if rows is None else self.sums[rows]), \
            "pseudobulk sum of raw UMI counts (source feature axis)"

    def mean_raw_counts(self, rows=None):
        """SUM / n_cells. Not log, not CP10K, not a fold change."""
        S = self.sums if rows is None else self.sums[rows]
        n = self.obs.n_cells.to_numpy()[rows] if rows is not None \
            else self.obs.n_cells.to_numpy()
        n = np.maximum(n.astype(float), 1.0)
        return sparse.diags(1.0 / n) @ S, "mean raw UMI counts per cell"

    def official_view(self, rows=None, quantity="sum"):
        """Reindex onto the 18,533 official features. Returns (dense, measured, semantics).

        Unmeasured columns are 0.0 in the array and False in the mask; the mask,
        not the zeros, is the truth. Mask your loss with it.
        """
        S, sem = (self.sum_counts(rows) if quantity == "sum"
                  else self.mean_raw_counts(rows))
        A = np.zeros((S.shape[0], len(self.official_to_source)), dtype=np.float64)
        m = self.measured
        A[:, m] = np.asarray(S[:, self.official_to_source[m]].todense())
        return A, m, sem + " reindexed to the official 18,533-gene order"

    def matched_control_rows(self, row, keys=("stimulation",)):
        """Row indices of the controls that share this row's biological context."""
        r = self.obs.iloc[row]
        sel = self.obs.is_control.astype(str).str.lower().isin(["true", "1"])
        for k in keys:
            sel &= self.obs[k].astype(str) == str(r[k])
        return np.nonzero(sel.to_numpy())[0]

    def group_key(self):
        return self.provenance.get("dataset_id", self.key)


def load_pb(key):
    d = ROOT / "pseudobulk" / key
    return PBSource(
        key=key,
        sums=sparse.load_npz(d / "sum_counts.npz"),
        obs=pd.read_csv(d / "obs.csv", low_memory=False),
        var=pd.read_csv(d / "var.csv", low_memory=False),
        provenance=json.loads((d / "provenance.json").read_text()),
        official_to_source=np.load(d / "official_to_source.npy"))


def list_pb():
    p = ROOT / "pseudobulk"
    return sorted(d.name for d in p.iterdir() if (d / "sum_counts.npz").exists()) \
        if p.exists() else []


# ----------------------------------------------------------------- effect source
@dataclass
class EffectSource:
    """Control-relative log1p(CP10K) deltas. NaN means the feature was unmeasured."""
    key: str
    delta: np.ndarray = field(repr=False)
    pairs: pd.DataFrame = field(repr=False)
    measured: np.ndarray = field(repr=False)
    provenance: dict = field(repr=False)

    @property
    def semantics(self):
        return self.provenance.get("quantity", "unknown")

    def for_target(self, target):
        i = np.nonzero(self.pairs.target_symbol.astype(str).to_numpy() == target)[0]
        return self.delta[i], self.pairs.iloc[i]

    def group_key(self):
        return self.provenance.get("store_key", self.key)


def load_effect(key):
    d = ROOT / "effects" / key
    return EffectSource(
        key=key, delta=np.load(d / "delta.npy"),
        pairs=pd.read_csv(d / "pairs.csv", low_memory=False),
        measured=np.load(d / "measured_mask.npy"),
        provenance=json.loads((d / "provenance.json").read_text()))


def list_effects():
    p = ROOT / "effects"
    return sorted(d.name for d in p.iterdir() if (d / "delta.npy").exists()) \
        if p.exists() else []


# ----------------------------------------------------------------- batch iterator
def iter_pb_batches(keys, batch=64, quantity="sum", targets_only=None, seed=0):
    """Yield training batches of (official-order matrix, mask, metadata).

    Each batch carries its own measured mask, because different assays measure
    different feature subsets and a single global mask would either discard
    usable columns or silently teach the model that unmeasured means zero.
    """
    rng = np.random.default_rng(seed)
    for key in keys:
        src = load_pb(key)
        idx = np.nonzero(~src.obs.is_control.astype(str).str.lower()
                         .isin(["true", "1"]).to_numpy())[0]
        if targets_only is not None:
            keep = src.obs.target_symbol.astype(str).isin(set(targets_only)).to_numpy()
            idx = idx[keep[idx]]
        rng.shuffle(idx)
        for s in range(0, len(idx), batch):
            rows = idx[s:s + batch]
            A, m, sem = src.official_view(rows, quantity=quantity)
            yield dict(store_key=key, X=A, measured=m, semantics=sem,
                       obs=src.obs.iloc[rows].reset_index(drop=True),
                       group_key=src.group_key())
