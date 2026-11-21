"""Read-only data and source-only anchors for the small residual experiment.

The baseline is a refitted b16-style covariance denoiser, not the archived b16
checkpoint. Every instance owns only the permitted response rows. In particular,
training anchors refit the ENTIRE baseline after excluding their own biology.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import sparse
from scipy.linalg import solve
from sklearn.utils.extmath import randomized_svd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reports/simple_residual_2026-09-16'


def canonical_bio(name):
    key = str(name).upper().replace('-', '').replace('_', '').replace(' ', '')
    return {'CVCL0023': 'A549', 'CVCL0320': 'HT29', 'HEPG2/C3A': 'HEPG2',
            'HEPG2C3A': 'HEPG2', 'C3A': 'HEPG2', 'K562ESSENTIALSCREEN': 'K562',
            'FENGIPSC': 'FENG', 'KOLF2.1J': 'KOLF', 'KOLF21J': 'KOLF',
            'VCC2025': 'H1'}.get(key, key)


@lru_cache(None)
def gene_map():
    import vcclib as V
    return V.GeneMap()


def canonical_target(name):
    return gene_map().canon(str(name))


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


@dataclass
class Table:
    y: np.ndarray
    controls: np.ndarray
    targets: np.ndarray
    bios: np.ndarray
    modality: np.ndarray
    features: np.ndarray
    weights: np.ndarray
    sources: np.ndarray
    genes: list

    def subset(self, ids):
        ids = np.asarray(ids, dtype=int)
        return Table(*(np.asarray(getattr(self, k)[ids]).copy() for k in
            ('y', 'controls', 'targets', 'bios', 'modality', 'features', 'weights', 'sources')),
            list(self.genes))

    def save(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=False)
        for key in ('y', 'controls', 'targets', 'bios', 'modality', 'features', 'weights', 'sources'):
            np.save(folder / (key+'.npy'), getattr(self, key))
        write_json(folder/'genes.json', list(self.genes))

    @classmethod
    def load(cls, folder, mmap='r'):
        folder = Path(folder)
        return cls(*(np.load(folder/(k+'.npy'), mmap_mode=mmap) for k in
            ('y', 'controls', 'targets', 'bios', 'modality', 'features', 'weights', 'sources')),
            json.loads((folder/'genes.json').read_text()))


def scope(table, excluded_bios=(), held_targets=(), modality=None):
    banned = {canonical_bio(x) for x in excluded_bios}
    bio = np.array([canonical_bio(x) for x in table.bios])
    keep = ~np.isin(bio, list(banned))
    # A held-out gene is removed in every genetic source and every stimulus.
    if held_targets:
        held = {canonical_target(t) for t in held_targets}
        canonical = np.array([canonical_target(t) if m==0 else str(t)
                              for t,m in zip(table.targets,table.modality)])
        keep &= ~((table.modality == 0) & np.isin(canonical, list(held)))
    if modality is not None:
        keep &= table.modality == modality
    return np.flatnonzero(keep)


def basis(x, rank):
    x = np.asarray(x, dtype=np.float32)
    k = min(int(rank), min(x.shape)-1)
    if k < 1 or not np.any(x):
        return np.eye(x.shape[1], min(max(k, 1), x.shape[1]), dtype=np.float32)
    return randomized_svd(x, n_components=k, n_iter=2, random_state=0)[2].T.copy()


class SourceBaseline:
    """Equal-biology empirical mean, source-only embedding ridge, covariance correction.

    Duplicated studies/states are pooled inside a biology using the supplied
    effective-count weights, then biologies receive equal weight. The denoiser
    uses centered source discrepancies, source-only PCA, scaled ridge=1 and a
    prediction-only amplitude match, following the archived b16 family.
    """
    def __init__(self, table, ids, rank=32):
        ids = np.asarray(ids, dtype=int)
        self.n_genes, self.feature_dim = table.y.shape[1], table.features.shape[1]
        self.rank = int(rank)
        self.b = np.empty((self.n_genes, 0), np.float32)
        self.coef = np.empty((0, self.n_genes), np.float32)
        self.global_coef = np.zeros((self.feature_dim+1, self.n_genes), np.float32)
        self.targets = np.array([], dtype=str)
        self.means = np.empty((0, self.n_genes), np.float32)
        self.lookup = {}
        if not len(ids):
            return
        if len(np.unique(table.modality[ids])) != 1:
            raise ValueError('Gene and drug baselines must be fitted separately')
        if not np.isfinite(table.y[ids]).all():
            raise ValueError('Baseline requires the explicit common measured gene axis')
        bio = np.array([canonical_bio(x) for x in table.bios[ids]])
        pairs = np.array([b+'\t'+str(t) for b, t in zip(bio, table.targets[ids])])
        unique, inverse = np.unique(pairs, return_inverse=True)
        w = np.asarray(table.weights[ids], dtype=np.float64)
        den = np.bincount(inverse, weights=w, minlength=len(unique))
        reducer = sparse.csr_matrix((w/den[inverse], (inverse, np.arange(len(ids)))),
                                    shape=(len(unique), len(ids)))
        yy = np.asarray(reducer @ table.y[ids], dtype=np.float32)
        ff = np.asarray(reducer @ table.features[ids], dtype=np.float32)
        ubio = np.array([u.split('\t', 1)[0] for u in unique])
        utarget = np.array([u.split('\t', 1)[1] for u in unique])
        self.targets, ti = np.unique(utarget, return_inverse=True)
        count = np.bincount(ti, minlength=len(self.targets))
        tr = sparse.csr_matrix((1/count[ti], (ti, np.arange(len(ti)))),
                              shape=(len(self.targets), len(ti)))
        self.means = np.asarray(tr @ yy, dtype=np.float32)
        feat = np.asarray(tr @ ff, dtype=np.float64)
        self.lookup = {str(t): i for i, t in enumerate(self.targets)}
        xx = np.c_[feat, np.ones(len(feat))]
        gram = xx.T @ xx
        reg = np.eye(len(gram)) * max(float(np.trace(gram)/len(gram)), 1e-8) * 10
        reg[-1, -1] *= 1e-6
        self.global_coef = solve(gram+reg, xx.T @ self.means,
            assume_a='pos', check_finite=False).astype(np.float32)
        if rank <= 0 or len(set(ubio)) < 2:
            return
        specific = yy.copy()
        for b in np.unique(ubio):
            own = ubio == b
            specific[own] -= specific[own].mean(0)
        total = np.zeros_like(self.means)
        np.add.at(total, ti, specific)
        valid = count[ti] > 1
        if valid.sum() < 2:
            return
        xx = (total[ti[valid]] - specific[valid]) / (count[ti[valid], None]-1)
        residual = specific[valid] - xx
        self.b = basis(self.means-self.means.mean(0), rank)
        zz = xx @ self.b
        # Equal biological contribution to denoiser fitting as well.
        counts = {b: int(np.sum(ubio[valid] == b)) for b in set(ubio[valid])}
        ww = np.array([1/counts[b] for b in ubio[valid]])
        ww *= len(ww)/ww.sum()
        gram = zz.T @ (zz * ww[:, None])
        penalty = max(float(np.trace(gram)/len(gram)), 1e-8)
        self.coef = solve(gram + np.eye(len(gram))*penalty,
            zz.T @ (residual * ww[:, None]), assume_a='pos', check_finite=False).astype(np.float32)

    def predict(self, targets, features, kind='denoiser'):
        f = np.asarray(features, dtype=np.float32)
        if f.shape != (len(targets), self.feature_dim):
            raise ValueError('Query feature shape does not match the gene/drug adapter')
        out = (np.c_[f, np.ones(len(f))] @ self.global_coef).astype(np.float32)
        seen = np.array([str(t) in self.lookup for t in targets])
        for j, t in enumerate(targets):
            if seen[j]:
                out[j] = self.means[self.lookup[str(t)]]
        if kind == 'mean':
            return out, seen
        if kind != 'denoiser':
            raise ValueError(kind)
        raw = out.copy()
        out -= out.mean(0)
        out += (out @ self.b) @ self.coef
        norm = np.abs(out.astype(np.float64)).sum()
        if norm > 0:
            out *= np.abs(raw.astype(np.float64)).sum()/norm
        # No empirical response exists for these targets. Their declared anchor
        # remains the source-trained global predictor, even for a singleton query.
        out[~seen] = raw[~seen]
        return out, seen


def training_anchors(table, ids, rank=32, kind='denoiser'):
    """Refit all mean/global/PCA/denoising parameters without the destination."""
    ids = np.asarray(ids, dtype=int)
    out = np.zeros((len(ids), table.y.shape[1]), dtype=np.float32)
    bios = np.array([canonical_bio(x) for x in table.bios[ids]])
    for mod in np.unique(table.modality[ids]):
        use_mod = table.modality[ids] == mod
        for b in np.unique(bios[use_mod]):
            query = use_mod & (bios == b)
            allowed = ids[use_mod & (bios != b)]
            fitted = SourceBaseline(table, allowed, rank=rank if mod == 0 else 0)
            out[query] = fitted.predict(table.targets[ids[query]],
                table.features[ids[query]], kind=kind if mod == 0 else 'mean')[0]
    return out


def balanced_weights(table, ids):
    ids = np.asarray(ids, dtype=int)
    groups = np.array([str(m)+'|'+canonical_bio(b) for m,b in
                       zip(table.modality[ids], table.bios[ids])])
    out = np.ones(len(ids), dtype=np.float64)
    for g in np.unique(groups):
        own = groups == g
        out[own] = 1/own.sum()
    return out/out.sum()
