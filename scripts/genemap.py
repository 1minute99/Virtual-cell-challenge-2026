"""Shared, collision-safe gene symbol harmonization (HGNC complete set).

Naive alias mapping is unsafe: HGNC lists MED1 as an alias of MBD4, so a naive
alias->symbol map silently merges the Mediator subunit MED1 into MBD4 and
manufactures false gene matches. We therefore:
  * never remap a token that is itself an Approved symbol, and
  * only accept an alias when it maps to exactly one Approved symbol.
Ensembl gene id is preferred as the join key wherever both sides carry it.
"""
import os
from collections import defaultdict
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HGNC_TSV = os.path.join(ROOT, "work/hgnc.tsv")
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


class GeneMap:
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

    def canon(self, g):
        """Approved symbol for g, else unique alias resolution, else g unchanged."""
        if g in self.approved:
            return g
        return self.alias.get(g, g)

    def ens(self, g):
        return self.sym2ens.get(self.canon(g))

    def stats(self):
        return dict(approved=len(self.approved), safe_aliases=len(self.alias),
                    with_ensembl=len(self.sym2ens))
