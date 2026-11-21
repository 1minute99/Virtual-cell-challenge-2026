#!/usr/bin/env python3
"""New context-free denoising comparators; historical sources stay read-only."""
from __future__ import annotations

import hashlib
import json
import pathlib
import time

import numpy as np
import torch
from scipy.linalg import solve
from torch import nn

import context_architectures as A
import jiang_benchmark as J
from run_context_architectures import extra_metrics

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/context_signal_2026-09-15/denoising"
RANKS = (8, 16)
STRENGTHS = (0.0, 0.25, 0.5, 1.0)
SEEDS = (0, 1, 2)
FAMILIES = ("moe", "linear", "ridge")


def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=True) + "\n")
    temp.replace(path)


class LinearResidual(A.ResidualModule):
    """One matrix, initialized to the identically seeded MoE's effective matrix."""
    def __init__(self, dims, use_context=False, **settings):
        super().__init__(dims, False, **settings)
        reference = A.MixtureOfExperts(dims, False, **settings)
        with torch.no_grad():
            weights = reference.gate.bias.softmax(-1)
            initial = torch.einsum("k,kpr->pr", weights, reference.experts)
        self.matrix = nn.Parameter(initial.clone())

    def forward(self, batch):
        return dict(mean=batch["p"] @ self.matrix)


# This process-local registration does not edit the archived implementation.
A.MODULES["linear"] = LinearResidual


class RidgeResidual:
    """Analytic common MSE + residual energy + scaled coefficient ridge objective."""
    family = "ridge"

    def __init__(self, prepared, shuffle_labels=False):
        start = time.monotonic()
        batch, reps = prepared.training_inputs(shuffle_labels, False)
        x = batch["p"].numpy().astype(np.float64)
        x -= x.mean(1, keepdims=True)
        y = reps.numpy().astype(np.float64).mean(1)
        x = x.reshape(-1, x.shape[-1])
        y = y.reshape(-1, y.shape[-1])
        gram = x.T @ x
        coefficient_alpha = 0.001
        energy_penalty = A.DEFAULTS["residual_penalty"]
        scale = max(np.trace(gram) / len(gram), 1e-12)
        self.matrix = solve((1 + energy_penalty) * gram +
                            np.eye(len(gram)) * coefficient_alpha * scale,
                            x.T @ y, assume_a="pos", check_finite=False)
        self.params = self.matrix.size
        self.seconds = time.monotonic() - start

    def predict(self, batch):
        out = batch["p"].numpy().astype(np.float64) @ self.matrix
        return dict(mean=out - out.mean(1, keepdims=True))

    def save(self, path):
        np.savez_compressed(path, matrix=self.matrix, seconds=self.seconds,
                            params=self.params)


def fit(prepared, family, seed=0, shuffle_labels=False):
    if family == "ridge":
        return RidgeResidual(prepared, shuffle_labels)
    return A.fit(prepared, family, use_context=False, seed=seed,
                 shuffle_labels=shuffle_labels)


def collapse_moe(model, batch):
    if model.family != "moe":
        raise ValueError("Only a context-free MoE has this collapse check")
    assert not model.module.use_context
    with torch.no_grad():
        weights = model.module.gate.bias.softmax(-1)
        matrix = torch.einsum("k,kpr->pr", weights, model.module.experts)
        explicit = A.center_correction(batch["p"] @ matrix).numpy()
    observed = model.predict(batch)["mean"]
    error = float(np.max(np.abs(explicit - observed)))
    np.testing.assert_allclose(explicit, observed, atol=3e-7, rtol=3e-6)
    return matrix.numpy(), error


def all_metrics(pred, reps, baseline, perms):
    values, target = J.metrics(pred, reps, perms)
    values.update(extra_metrics(pred, reps, baseline))
    base = float(J.pds(baseline, reps.mean(0)).mean())
    den = base - values["shuffle"]
    values["shuffle_normalized_ratio_vs_b16"] = (
        (values["pds"] - values["shuffle"]) / den if abs(den) > 1e-12 else np.nan)
    return values, target


def choose(rows):
    groups = {}
    for row in rows:
        key = (int(row["rank"]), float(row["strength"]))
        groups.setdefault(key, []).append(float(row["pds"]))
    average = {k: float(np.mean(v)) for k, v in groups.items()}
    key = max(average, key=lambda k: (average[k], -k[1], -k[0]))
    return dict(rank=key[0], strength=key[1], inner_pds=average[key])


def freeze_spec(path, spec, sources):
    path = pathlib.Path(path)
    payload = {**spec, "source_hashes": {p: sha(ROOT / p) for p in sources}}
    if path.exists():
        if json.loads(path.read_text()) != payload:
            raise RuntimeError("Existing spec differs; preserve it and use a new run folder")
    else:
        write_json(path, payload)
        snap = path.parent / "sources"
        snap.mkdir(exist_ok=True)
        for source in sources:
            (snap / pathlib.Path(source).name).write_bytes((ROOT / source).read_bytes())
    return payload
