"""Build pseudobulk perturbation deltas from Replogle K562 genome-scale bulk.

delta_source[p] = log1p(CP10K(pseudobulk_p)) - log1p(CP10K(pseudobulk_control))

Rows of K562_gwps_raw_bulk are per-perturbation mean UMI vectors; control is the
mean over the `non-targeting` rows. Deltas are returned in a requested symbol
gene space (genes absent from the source get delta 0 = pass through control).
"""
import os, sys
import numpy as np, anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genemap import GeneMap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GWPS_BULK = os.path.join(ROOT, "data/replogle_gwps/K562_gwps_raw_bulk_01.h5ad")


def _cp10k_log1p(M):
    s = M.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return np.log1p(M / s * 1e4)


def load_source_deltas(target_genes, gm=None, path=GWPS_BULK):
    """Return (delta dict {pert_symbol: vec[len(target_genes)]}, info dict)."""
    gm = gm or GeneMap()
    A = ad.read_h5ad(path)
    X = np.asarray(A.X, dtype=np.float64)
    L = _cp10k_log1p(X)

    perts = np.array([x.split("_")[1] for x in A.obs_names])
    is_ctrl = perts == "non-targeting"
    ctrl = L[is_ctrl].mean(0)

    src_sym = [gm.canon(g) for g in A.var["gene_name"].astype(str)]
    pos = {}
    for i, s in enumerate(src_sym):
        pos.setdefault(s, i)

    tgt_c = [gm.canon(g) for g in target_genes]
    idx_src, idx_tgt = [], []
    for j, s in enumerate(tgt_c):
        if s in pos:
            idx_src.append(pos[s]); idx_tgt.append(j)
    idx_src = np.asarray(idx_src); idx_tgt = np.asarray(idx_tgt)

    # average replicate rows (multiple guides/transcripts) per perturbation
    deltas = {}
    order = {}
    for i, p in enumerate(perts):
        if p == "non-targeting":
            continue
        order.setdefault(p, []).append(i)
    for p, rows in order.items():
        d = np.zeros(len(target_genes), dtype=np.float32)
        d[idx_tgt] = (L[rows].mean(0) - ctrl)[idx_src].astype(np.float32)
        deltas[gm.canon(p)] = d

    info = dict(n_source_rows=len(A), n_control_rows=int(is_ctrl.sum()),
                n_perturbations=len(deltas), n_source_genes=X.shape[1],
                n_genes_mapped=len(idx_tgt), n_genes_requested=len(target_genes),
                n_genes_dropped=len(target_genes) - len(idx_tgt))
    return deltas, info, idx_tgt
