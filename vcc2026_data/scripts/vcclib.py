#!/usr/bin/env python3
"""Shared helpers: official VCC 2026 target/feature definitions and gene mapping.

Two things this module refuses to conflate:

  * `pert_counts.csv`  -> the 300 perturbation TARGETS to predict (`target_gene`).
  * `gene_names.csv`   -> the 18,533 output RNA FEATURES, in submission order.

Loading either one as the other is the single most damaging mistake available
here, so `official_targets()` and `official_features()` each assert the row count
and the header they expect, and refuse a file that looks like the other one.
"""
import hashlib
import os
import pathlib
from collections import defaultdict

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
# Official 2026 validation controls bundle, downloaded with `vcc datasets download controls`.
OFFICIAL_DIR = pathlib.Path(
    os.environ.get("VCC_OFFICIAL_DIR",
                   "/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls"))
HGNC_TSV = pathlib.Path(os.environ.get(
    "HGNC_TSV", "/home/wonmin-kim/AIVC/vcc2026/work/hgnc.tsv"))

N_TARGETS = 300
N_FEATURES = 18533


def sha256(path, chunk=32 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while (b := f.read(chunk)):
            h.update(b)
    return h.hexdigest()


def md5(path, chunk=32 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while (b := f.read(chunk)):
            h.update(b)
    return h.hexdigest()


def official_targets(path=None):
    """The 300 perturbation targets. Rejects the 18,533-feature file."""
    p = pathlib.Path(path or OFFICIAL_DIR / "pert_counts.csv")
    df = pd.read_csv(p)
    if "gene_name" in df.columns and "target_gene" not in df.columns:
        raise ValueError(f"{p} has a 'gene_name' header: that is the OUTPUT FEATURE "
                         f"list, not the perturbation target list")
    if "target_gene" not in df.columns:
        raise ValueError(f"{p}: no 'target_gene' column; columns={list(df.columns)}")
    t = [str(x).strip() for x in df["target_gene"] if str(x).strip()]
    if len(t) != N_TARGETS or len(set(t)) != N_TARGETS:
        raise ValueError(f"{p}: expected {N_TARGETS} unique targets, "
                         f"got {len(t)} rows / {len(set(t))} unique")
    return t


def official_features(path=None):
    """The 18,533 output RNA features, in submission order. Rejects the target file."""
    p = pathlib.Path(path or OFFICIAL_DIR / "gene_names.csv")
    df = pd.read_csv(p)
    if "target_gene" in df.columns and "gene_name" not in df.columns:
        raise ValueError(f"{p} has a 'target_gene' header: that is the TARGET list, "
                         f"not the output feature list")
    col = "gene_name" if "gene_name" in df.columns else df.columns[0]
    g = [str(x).strip() for x in df[col] if str(x).strip()]
    if len(g) != N_FEATURES or len(set(g)) != N_FEATURES:
        raise ValueError(f"{p}: expected {N_FEATURES} unique features, "
                         f"got {len(g)} rows / {len(set(g))} unique")
    return g


class GeneMap:
    """Collision-safe HGNC symbol harmonisation.

    A naive alias->symbol map is unsafe here: HGNC lists MED1 as an alias of
    MBD4, so a naive map silently merges the Mediator subunit into MBD4 and
    manufactures matches. Two rules prevent that:
      * never remap a token that is itself an Approved symbol, and
      * accept an alias only when it resolves to exactly one Approved symbol.
    Ensembl gene id is preferred as the join key wherever both sides carry one.
    """

    def __init__(self, path=HGNC_TSV):
        df = pd.read_csv(path, sep="\t", low_memory=False)
        df = df[df.status == "Approved"]
        self.approved = set(df.symbol.dropna())
        self.sym2ens = {s: e for s, e in zip(df.symbol, df.ensembl_gene_id)
                        if isinstance(e, str)}
        self.ens2sym = {e: s for s, e in self.sym2ens.items()}
        cand = defaultdict(set)
        for _, r in df[["symbol", "alias_symbol", "prev_symbol"]].iterrows():
            s = r["symbol"]
            if not isinstance(s, str):
                continue
            for fld in ("alias_symbol", "prev_symbol"):
                v = r[fld]
                if isinstance(v, str):
                    for a in v.split("|"):
                        a = a.strip()
                        if a and a not in self.approved:
                            cand[a].add(s)
        self.alias = {a: next(iter(v)) for a, v in cand.items() if len(v) == 1}
        self.ambiguous = {a: sorted(v) for a, v in cand.items() if len(v) > 1}

    def canon(self, g):
        """Approved symbol for g, else a unique alias resolution, else g unchanged."""
        g = str(g).strip()
        if g in self.approved:
            return g
        return self.alias.get(g, g)

    def ens(self, g):
        """Ensembl gene id for a symbol, version suffix never assumed."""
        return self.sym2ens.get(self.canon(g))

    def from_ens(self, e):
        """Approved symbol for an Ensembl id, tolerating a version suffix."""
        e = str(e).strip()
        return self.ens2sym.get(e.split(".")[0])

    def stats(self):
        return dict(approved=len(self.approved), unique_aliases=len(self.alias),
                    ambiguous_aliases=len(self.ambiguous), with_ensembl=len(self.sym2ens))


def feature_index(features=None):
    """symbol -> position in the official submission gene order."""
    features = features or official_features()
    return {g: i for i, g in enumerate(features)}
