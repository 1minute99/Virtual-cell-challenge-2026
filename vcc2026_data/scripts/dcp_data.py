"""Fold-aware assembly: main effects, context residuals, and the training-only bases.

Everything in here takes an explicit list of TRAINING contexts. Nothing is ever fitted on a
held-out context: not the context PCA, not the perturbation main effect, not the response basis.
The held-out context contributes exactly one thing -- its untreated control profile -- which is
available at inference time in the real task as well.
"""
from pathlib import Path
import json
import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC/vcc2026_data')
PROC = ROOT / 'data/processed/drug_context_pretrain'


def pca(X, dim, weights=None):
    """Top-`dim` right singular vectors of the (optionally weighted) centred matrix.

    Randomized SVD: the residual matrices here are ~8.5k x 18k and only rank 16-32 is ever
    wanted, so a full LAPACK SVD would cost minutes per fold for no extra information.
    """
    from sklearn.utils.extmath import randomized_svd
    X = np.asarray(X, np.float32)
    if weights is None:
        Xc = X - X.mean(0)
    else:
        w = (weights / max(weights.sum(), 1e-9)).astype(np.float32)
        Xc = (X - (w[:, None] * X).sum(0)) * np.sqrt(w)[:, None]
    dim = int(min(dim, min(Xc.shape) - 1))
    if dim < 1 or not np.any(Xc):
        return np.eye(X.shape[1], max(dim, 1), dtype=np.float32)
    return randomized_svd(Xc, n_components=dim, n_iter=4, random_state=0)[2].T.copy()


class DrugData:
    """One drug dataset condensed to a (context x perturbation x gene) response cube."""

    def __init__(self, name):
        z = np.load(PROC / f'{name}.npz', allow_pickle=True)
        self.name = name
        self.ctx_names = [str(x) for x in z['ctx_names']]
        self.ctx_keys = [f'{name}:{c}' for c in self.ctx_names]
        self.pert_names = [str(x) for x in z['pert_names']]
        self.genes = [str(x) for x in z['genes']]
        C, P, G = len(self.ctx_names), len(self.pert_names), len(self.genes)
        resp, ctx, pert, w = z['resp'], z['ctx'], z['pert'], z['weight']
        rep = z['is_repeat']
        self.R = np.zeros((C, P, G), np.float32)
        self.W = np.zeros((C, P), np.float64)
        # split halves for reliability: `rep` marks Tahoe's plate14 / GSE's replicate 2
        self.RA = np.zeros((C, P, G), np.float32); self.WA = np.zeros((C, P))
        self.RB = np.zeros((C, P, G), np.float32); self.WB = np.zeros((C, P))
        acc = np.zeros((C, P, G), np.float64)
        for r in range(len(resp)):
            c, p, ww = ctx[r], pert[r], w[r]
            acc[c, p] += resp[r] * ww
            self.W[c, p] += ww
            if rep[r]:
                self.RB[c, p] += resp[r] * ww; self.WB[c, p] += ww
            else:
                self.RA[c, p] += resp[r] * ww; self.WA[c, p] += ww
        m = self.W > 0
        self.R[m] = (acc[m] / self.W[m][:, None]).astype(np.float32)
        for Rx, Wx in ((self.RA, self.WA), (self.RB, self.WB)):
            mm = Wx > 0
            Rx[mm] = (Rx[mm] / Wx[mm][:, None]).astype(np.float32)
        del acc
        self.M = m
        self.ctx_profile = None       # filled by attach_contexts

    def coverage(self):
        return dict(contexts=len(self.ctx_names), perturbations=len(self.pert_names),
                    filled=int(self.M.sum()), density=float(self.M.mean()),
                    both_halves=int(((self.WA > 0) & (self.WB > 0)).sum()))

    def reliability(self):
        """Split-half response reproducibility, pooled over filled cells with both halves."""
        m = (self.WA > 0) & (self.WB > 0)
        if m.sum() == 0:
            return dict(n=0)
        a, b = self.RA[m], self.RB[m]
        num = (a * b).sum(1)
        den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
        cs = num / np.maximum(den, 1e-9)
        return dict(n=int(m.sum()), median_split_half_cosine=float(np.median(cs)),
                    mean_split_half_cosine=float(cs.mean()),
                    frac_above_0p3=float((cs > 0.3).mean()))

    # ---------------- fold-conditional quantities ----------------

    def main_effect(self, train_ctx):
        """mean response of each perturbation across TRAINING contexts only."""
        tr = np.asarray(train_ctx, int)
        w = self.W[tr] * self.M[tr]
        s = w.sum(0)
        main = np.zeros((len(self.pert_names), self.R.shape[2]), np.float32)
        ok = s > 0
        main[ok] = (np.einsum('cp,cpg->pg', w, self.R[tr].astype(np.float64))[ok]
                    / s[ok][:, None]).astype(np.float32)
        n_ctx = self.M[tr].sum(0)
        return main, ok, n_ctx

    def fold(self, train_ctx, rank, min_ctx=3):
        """Everything a model needs for one leave-context-out split."""
        tr = np.asarray(train_ctx, int)
        main, ok, n_ctx = self.main_effect(tr)
        usable = ok & (n_ctx >= min_ctx)
        pert_basis = pca(main[usable], rank)
        zp = (main - main[usable].mean(0)) @ pert_basis
        rows = [(c, p) for c in tr for p in np.flatnonzero(usable & self.M[c])]
        resid = np.stack([self.R[c, p] - main[p] for c, p in rows]).astype(np.float32)
        resp_basis = pca(resid, rank)
        y = resid @ resp_basis
        return dict(main=main, usable=usable, pert_basis=pert_basis, zp=zp,
                    resp_basis=resp_basis, rows=np.array(rows),
                    y=y, w=np.array([self.W[c, p] for c, p in rows]),
                    resid_rms=float(np.sqrt((resid ** 2).mean())),
                    basis_var_kept=float((y ** 2).sum() / max((resid ** 2).sum(), 1e-9)))

    def eval_rows(self, ctxs, usable):
        return np.array([(c, p) for c in np.atleast_1d(ctxs)
                         for p in np.flatnonzero(usable & self.M[c])])


def context_space(train_names, dim=64):
    """z_context for EVERY context, from a PCA fitted on the training contexts only.

    `train_names` are DATASET-QUALIFIED keys ("tahoe:A549", "gse:A549", "gene:HEK293T").

    Uses the per-dataset-centred control profiles, because the raw ones separate almost
    perfectly by assay (see preprocessing_diagnostics.json).
    """
    z = np.load(PROC / 'contexts.npz', allow_pickle=True)
    names = [str(x) for x in z['names']]
    ds = np.array([str(x) for x in z['dataset']])
    Xp = np.asarray(z['profile_adj'], np.float64)
    idx = {n: i for i, n in enumerate(names)}
    tr = np.array([idx[n] for n in train_names])
    cm = Xp[tr].mean(0)
    B = pca(Xp[tr] - cm, dim)
    Z = ((Xp - cm) @ B).astype(np.float32)
    s = float(np.sqrt(np.mean(Z[tr] ** 2)))
    Z = Z / max(s, 1e-9)
    var = float((((Xp[tr] - cm) @ B) ** 2).sum() / max(((Xp[tr] - cm) ** 2).sum(), 1e-9))
    return dict(Z=Z, names=names, dataset=ds, index=idx, var_kept=var, dim=int(B.shape[1]))


def gene_context_names():
    z = np.load(PROC / 'gene_contexts.npz', allow_pickle=True)
    return [str(x) for x in z['ctx_names']]


def shared_genes():
    return json.load(open(PROC / 'context_axis.json'))['genes']
