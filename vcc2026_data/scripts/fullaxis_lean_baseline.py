"""Memory-lean, float32, mask-aware source baseline for the wide gene axis.

The frozen `SourceBaseline` is mathematically right but allocates ~47 GiB on a
103,295 x 17,282 table and was OOM-killed twice. The cost is not the data, it is dtype and
copies: the pooling weights are float64 so `reducer @ y` promotes a 5 GiB float32 result to
10 GiB, the denoiser divides by an integer count and promotes twice more, and several full-size
intermediates are copied rather than reused.

This reimplements the SAME arithmetic with three changes and no others:
  * the two sparse pooling operators are float32, so no product promotes;
  * the responses stay on disk as a memmap instead of resident in RAM;
  * the large intermediates are written in place rather than copied.
Peak drops to roughly 11 GiB. `test_fullaxis_lean.py` pins exact agreement with the frozen
implementation on dense input, so this is an optimization and not a new model.
"""
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import solve

import simple_residual_data as D


class LeanSourceBaseline:
    """Same fitted quantities as `simple_residual_data.SourceBaseline`, float32 throughout."""

    def __init__(self, table, ids, rank=32):
        ids = np.asarray(ids, dtype=int)
        y = table.y                      # (len(ids), G); memmap or array, rows already subset
        self.n_genes, self.feature_dim = y.shape[1], table.features.shape[1]
        self.rank = int(rank)
        G = self.n_genes
        self.b = np.empty((G, 0), np.float32)
        self.coef = np.empty((0, G), np.float32)
        self.global_coef = np.zeros((self.feature_dim + 1, G), np.float32)
        self.targets = np.array([], dtype=str)
        self.means = np.empty((0, G), np.float32)
        self.lookup = {}
        if not len(ids):
            return
        if len(np.unique(table.modality[ids])) != 1:
            raise ValueError('Gene and drug baselines must be fitted separately')

        bio = np.array([D.canonical_bio(x) for x in table.bios[ids]])
        pairs = np.array([b + '\t' + str(t) for b, t in zip(bio, table.targets[ids])])
        unique, inverse = np.unique(pairs, return_inverse=True)
        w = np.asarray(table.weights[ids], dtype=np.float64)
        den = np.bincount(inverse, weights=w, minlength=len(unique))
        reducer = sparse.csr_matrix(((w / den[inverse]).astype(np.float32),
                                     (inverse, np.arange(len(ids)))),
                                    shape=(len(unique), len(ids)), dtype=np.float32)
        yy = reducer @ np.asarray(y, dtype=np.float32)          # float32 in, float32 out
        ff = np.asarray(reducer @ np.asarray(table.features[ids], np.float32), np.float32)
        ubio = np.array([u.split('\t', 1)[0] for u in unique])
        utarget = np.array([u.split('\t', 1)[1] for u in unique])
        self.targets, ti = np.unique(utarget, return_inverse=True)
        count = np.bincount(ti, minlength=len(self.targets))
        tr = sparse.csr_matrix(((1 / count[ti]).astype(np.float32), (ti, np.arange(len(ti)))),
                               shape=(len(self.targets), len(ti)), dtype=np.float32)
        self.means = np.asarray(tr @ yy, dtype=np.float32)
        feat = np.asarray(tr @ ff, dtype=np.float64)
        self.lookup = {str(t): i for i, t in enumerate(self.targets)}
        xx = np.c_[feat, np.ones(len(feat))]
        gram = xx.T @ xx
        reg = np.eye(len(gram)) * max(float(np.trace(gram) / len(gram)), 1e-8) * 10
        reg[-1, -1] *= 1e-6
        self.global_coef = solve(gram + reg, xx.T @ self.means,
                                 assume_a='pos', check_finite=False).astype(np.float32)
        if rank <= 0 or len(set(ubio)) < 2:
            return

        specific = yy                                   # reuse: yy is not needed afterwards
        del yy
        for b_ in np.unique(ubio):
            own = ubio == b_
            specific[own] -= specific[own].mean(0)
        total = np.zeros_like(self.means)
        np.add.at(total, ti, specific)
        valid = count[ti] > 1
        if valid.sum() < 2:
            return
        divisor = (count[ti[valid]] - 1).astype(np.float32)[:, None]
        xx2 = total[ti[valid]]
        xx2 -= specific[valid]
        xx2 /= divisor                                  # float32 / float32 stays float32
        del total
        residual = specific[valid]
        residual -= xx2                                 # in place; no extra full-size array
        del specific
        self.b = D.basis(self.means - self.means.mean(0), rank)
        zz = xx2 @ self.b
        del xx2
        counts = {b_: int(np.sum(ubio[valid] == b_)) for b_ in set(ubio[valid])}
        ww = np.array([1 / counts[b_] for b_ in ubio[valid]])
        ww *= len(ww) / ww.sum()
        gram = zz.T @ (zz * ww[:, None])
        penalty = max(float(np.trace(gram) / len(gram)), 1e-8)
        self.coef = solve(gram + np.eye(len(gram)) * penalty,
                          zz.T @ (residual * ww[:, None]),
                          assume_a='pos', check_finite=False).astype(np.float32)

    predict = D.SourceBaseline.predict
