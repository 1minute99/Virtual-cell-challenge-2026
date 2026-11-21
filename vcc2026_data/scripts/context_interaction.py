#!/usr/bin/env python3
"""Small signed response/context interaction model with whole-source exclusions.

Every Engine contains only the caller's allowed training indices. Source-only
denoising inputs for its training examples additionally exclude that example's
destination source. PCA and regression use training data; queries supply controls.
There is no context-only intercept, so the new correction cannot be identical for
every target. The output remains a signed vector, rather than a positive gate.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh, solve


def basis(x, rank):
    """Deterministic right singular vectors using the smaller row Gram matrix."""
    ev, u = eigh(x @ x.T, check_finite=False)
    order = np.argsort(ev)[::-1]
    order = order[ev[order] > max(float(ev.max()), 1e-30) * 1e-10][:rank]
    return (x.T @ u[:, order]) / np.sqrt(ev[order])[None, :]


def centered(x):
    return x - x.mean(axis=-2, keepdims=True)


def match_amplitude(pred, original):
    numerator, denominator = np.abs(original).sum(), np.abs(pred).sum()
    return pred * (numerator / denominator) if denominator > 0 else pred.copy()


def derangement(n, rng):
    while True:
        order = rng.permutation(n)
        if np.all(order != np.arange(n)):
            return order


@dataclass
class Denoising:
    original: np.ndarray
    mean: np.ndarray
    correction: np.ndarray

    def predict(self, strength=1.0):
        return match_amplitude(self.mean + strength * self.correction, self.original)


class Engine:
    def __init__(self, effects, controls, train, operator_rank=128, operator_alpha=1.0):
        # Copy only allowed data: no method below can read a held-out response.
        self.train = tuple(train)
        self.responses = np.asarray(effects[list(train)], dtype=np.float64).mean(axis=1)
        self.controls = np.asarray(controls[list(train)], dtype=np.float64).copy()
        self.specific = centered(self.responses)
        self.operator_rank, self.operator_alpha = operator_rank, operator_alpha
        self.cache = {}
        self.encoding_cache = {}
        self.full = self.denoiser(tuple(range(len(train))))

    def denoiser(self, indices):
        indices = tuple(indices)
        if indices in self.cache:
            return self.cache[indices]
        src = self.specific[list(indices)]
        raw = self.responses[list(indices)].mean(0)
        mean = src.mean(0)
        c = len(src)
        if c < 2:
            correction = np.zeros_like(mean)
        else:
            b = basis(mean, self.operator_rank)
            x = (src.sum(0)[None, :, :] - src) / (c - 1)
            z = (x @ b).reshape(-1, b.shape[1])
            y = (src - x).reshape(-1, mean.shape[1])
            gram = z.T @ z
            scale = max(float(np.trace(gram) / len(gram)), 1e-12)
            coef = solve(gram + np.eye(len(gram)) * self.operator_alpha * scale,
                         z.T @ y, assume_a="pos", check_finite=False)
            correction = (mean @ b) @ coef
        out = Denoising(raw, mean, correction)
        self.cache[indices] = out
        return out

    def encode(self, rank, control_rank=2):
        key = (rank, control_rank)
        if key in self.encoding_cache:
            return self.encoding_cache[key]
        t, n, g = self.responses.shape
        bx = basis(self.specific.mean(0), rank)
        by = basis(self.specific.reshape(t*n, g), rank)
        leave_mean, residual = [], []
        for j in range(t):
            d = self.denoiser(tuple(i for i in range(t) if i != j))
            leave_mean.append(d.mean)
            residual.append(self.specific[j] - d.predict())
        z = np.stack(leave_mean) @ bx
        zscale = np.maximum(np.sqrt(np.mean(z*z, axis=(0, 1))), 1e-10)
        z /= zscale
        control_mean = self.controls.mean(0)
        bc = basis(self.controls - control_mean, min(control_rank, t-1))
        h = (self.controls - control_mean) @ bc
        hscale = np.maximum(np.sqrt(np.mean(h*h, axis=0)), 1e-10)
        h /= hscale
        out = dict(bx=bx, by=by, z=z, zscale=zscale, h=h, hscale=hscale,
                   bc=bc, control_mean=control_mean,
                   residual=np.stack(residual) @ by,
                   query_z=(self.full.mean @ bx)/zscale)
        self.encoding_cache[key] = out
        return out

    def fit(self, rank, alpha, use_context=True, control_rank=2,
            shuffle_residual_labels=False, shuffle_training_contexts=False, seed=903):
        a = self.encode(rank, control_rank)
        t, n, r = a["z"].shape
        rng = np.random.default_rng(seed)
        h = a["h"]
        if shuffle_training_contexts:
            h = h[derangement(t, rng)]
        phi = make_features(a["z"], h, use_context).reshape(t*n, -1)
        pscale = np.maximum(np.sqrt(np.mean(phi*phi, axis=0)), 1e-10)
        phi /= pscale
        y = a["residual"].copy()
        if shuffle_residual_labels:
            # Only the new correction's labels are shuffled. The denoiser stays
            # fixed, making this a control for its incremental supervised gain.
            for i in range(t):
                y[i] = y[i, derangement(n, rng)]
        gram = phi.T @ phi
        scale = max(float(np.trace(gram) / len(gram)), 1e-12)
        coef = solve(gram + np.eye(len(gram))*alpha*scale,
                     phi.T @ y.reshape(t*n, -1), assume_a="pos", check_finite=False)
        return Model(a, coef, pscale, self.full.predict(), use_context,
                     dict(rank=rank, actual_response_rank=r,
                          actual_output_rank=a["by"].shape[1],
                          actual_control_rank=a["h"].shape[1], alpha=alpha,
                          use_context=use_context, n_features=phi.shape[1],
                          n_coefficients=int(coef.size), training_sources=t,
                          training_examples=t*n,
                          shuffle_residual_labels=shuffle_residual_labels,
                          shuffle_training_contexts=shuffle_training_contexts))


def make_features(z, h, use_context):
    if not use_context or h.shape[-1] == 0:
        return z.copy()
    cross = z[..., :, :, None] * h[..., None, None, :]
    return np.concatenate([z, cross.reshape(*z.shape[:-1], -1)], axis=-1)


@dataclass
class Model:
    encoding: dict
    coef: np.ndarray
    feature_scale: np.ndarray
    baseline: np.ndarray
    use_context: bool
    metadata: dict

    def predict(self, control, zero_context=False):
        a = self.encoding
        h = ((control-a["control_mean"]) @ a["bc"])/a["hscale"]
        if zero_context:
            h = np.zeros_like(h)
        phi = make_features(a["query_z"], h, self.use_context)/self.feature_scale
        correction = (phi @ self.coef) @ a["by"].T
        return self.baseline + correction

    def save(self, path):
        a = self.encoding
        np.savez_compressed(path, coef=self.coef, feature_scale=self.feature_scale,
                            baseline=self.baseline, query_z=a["query_z"],
                            response_basis=a["bx"], response_scale=a["zscale"],
                            output_basis=a["by"], control_mean=a["control_mean"],
                            control_basis=a["bc"], control_scale=a["hscale"],
                            use_context=self.use_context)
