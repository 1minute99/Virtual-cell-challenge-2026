#!/usr/bin/env python3
"""Eight compact residual modules sharing frozen source-only Jiang features.

No module receives a Panel or held-out responses. Each sees the same Batch:
source responses, source contexts, query context, perturbation representation,
and a fixed permitted-source denoiser prediction. Targets are a batch axis,
never an embedding/identifier. All corrections are centered across targets.
"""
from __future__ import annotations

from dataclasses import dataclass
import copy
import hashlib
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from context_interaction import Engine, centered, derangement


FAMILIES = ("pairwise", "hypernetwork", "deepsets", "moe", "tensor", "cnp",
            "metric", "gp")
DEFAULTS = dict(steps=160, lr=0.01, weight_decay=0.001, residual_penalty=0.5,
                hidden=16, interaction_rank=4, experts=3, context_rank=2)


def as_tensor(x):
    return torch.as_tensor(np.asarray(x, dtype=np.float32).copy())


def masked_mean(x, mask):
    # Responses are episode x target x source x component.
    return (x * mask[:, None, :, None]).sum(2) / mask.sum(1)[:, None, None]


def context_values(batch, use_context):
    if use_context:
        return batch["source_c"], batch["query_c"]
    return torch.zeros_like(batch["source_c"]), torch.zeros_like(batch["query_c"])


def center_correction(x):
    return x - x.mean(dim=1, keepdim=True)


@dataclass
class Prepared:
    train: tuple
    rank: int
    basis_x: np.ndarray
    basis_y: np.ndarray
    scale_x: np.ndarray
    scale_y: float
    basis_c: np.ndarray
    mean_c: np.ndarray
    scale_c: np.ndarray
    baseline: np.ndarray
    original_mean: np.ndarray
    training: dict
    residual_reps: torch.Tensor
    destinations: np.ndarray
    query_template: dict
    training_h: np.ndarray
    feature_seconds: float

    @property
    def dims(self):
        return (self.basis_x.shape[1], self.basis_y.shape[1], self.basis_c.shape[1])

    def query(self, control, zero_context=False):
        h = ((np.asarray(control) - self.mean_c) @ self.basis_c) / self.scale_c
        if zero_context:
            h = np.zeros_like(h)
        return {**self.query_template, "query_c": as_tensor(h[None])}

    def training_inputs(self, shuffle_labels=False, shuffle_contexts=False):
        batch = {k: v.clone() for k, v in self.training.items()}
        labels = self.residual_reps.clone()
        rng = np.random.default_rng(903)
        if shuffle_labels:
            n = labels.shape[2]
            perms = [derangement(n, rng) for _ in self.train]
            for e, dst in enumerate(self.destinations):
                labels[e] = labels[e, :, perms[dst], :]
        if shuffle_contexts:
            order = derangement(len(self.train), rng)
            hs = as_tensor(self.training_h[order])
            batch["query_c"] = hs[as_tensor(self.destinations).long()]
            source_ids = batch.pop("source_ids").long()
            batch["source_c"] = hs[source_ids] * batch["mask"][:, :, None]
        else:
            batch.pop("source_ids")
        return batch, labels

    def decode(self, correction, strength=1.0):
        delta = np.asarray(correction, dtype=np.float64) * self.scale_y
        return self.baseline + float(strength) * (delta @ self.basis_y.T)

    def save_features(self, path, targets, genes):
        arrays = dict(basis_x=self.basis_x, basis_y=self.basis_y,
                      scale_x=self.scale_x, scale_y=self.scale_y,
                      basis_c=self.basis_c, mean_c=self.mean_c, scale_c=self.scale_c,
                      baseline=self.baseline, original_mean=self.original_mean,
                      train=np.asarray(self.train), targets=np.asarray(targets),
                      genes=np.asarray(genes))
        arrays.update({"query_"+k: v.numpy() for k, v in self.query_template.items()})
        np.savez_compressed(path, **arrays)


def prepare(effects, controls, train, rank=8, context_rank=2):
    start = time.monotonic()
    train = tuple(sorted(train))
    assert len(train) >= 3, "Nested Jiang folds must retain at least three sources"
    engine = Engine(effects, controls, train, operator_rank=128, operator_alpha=1.0)
    enc = engine.encode(rank, context_rank)
    bx, by = enc["bx"], enc["by"]
    scale_y = max(float(np.sqrt(np.mean((engine.specific @ by)**2))), 1e-8)
    source_x = (engine.specific @ bx) / enc["zscale"]
    source_y = (engine.specific @ by) / scale_y
    # Preserve both replicate labels. The feature input uses their mean, as before.
    allowed_reps = centered(np.asarray(effects[list(train)], dtype=np.float64))
    truth_reps = (allowed_reps @ by) / scale_y
    h = enc["h"]
    c, n, _ = source_x.shape

    def episode(indices, query_h):
        d = engine.denoiser(indices)
        baseline = d.predict()
        srcx = np.zeros((n, c, bx.shape[1]), dtype=np.float32)
        srcy = np.zeros((n, c, by.shape[1]), dtype=np.float32)
        srcc = np.zeros((c, h.shape[1]), dtype=np.float32)
        mask, ids = np.zeros(c, np.float32), np.zeros(c, np.int64)
        for slot, i in enumerate(indices):
            srcx[:, slot] = source_x[i]
            srcy[:, slot] = source_y[i]
            srcc[slot] = h[i]
            mask[slot], ids[slot] = 1, i
        return dict(p=source_x[list(indices)].mean(0), source_x=srcx,
                    source_y=srcy, source_c=srcc, query_c=query_h, mask=mask,
                    baseline=(baseline @ by)/scale_y, source_ids=ids), baseline

    episodes, labels, destinations = [], [], []
    for dst in range(c):
        full = tuple(i for i in range(c) if i != dst)
        # Identical variable-set episodes for every architecture. There is no
        # destination response in either the full or reduced source set.
        drop = full[dst % len(full)]
        smaller = tuple(i for i in full if i != drop)
        for sources in (full, smaller):
            b, _ = episode(sources, h[dst])
            episodes.append(b)
            labels.append(truth_reps[dst] - b["baseline"][None])
            destinations.append(dst)
    training = {k: as_tensor(np.stack([e[k] for e in episodes])) for k in episodes[0]}
    query, baseline = episode(tuple(range(c)), np.zeros(h.shape[1]))
    query.pop("query_c")
    query.pop("source_ids")
    query = {k: as_tensor(v[None]) for k, v in query.items()}
    baseline.setflags(write=False)
    return Prepared(train, rank, bx, by, enc["zscale"], scale_y,
                    enc["bc"], enc["control_mean"], enc["hscale"], baseline,
                    engine.full.original.copy(), training, as_tensor(np.stack(labels)),
                    np.asarray(destinations), query, h.copy(), time.monotonic()-start)


def mlp(in_dim, hidden, out_dim, zero_output=True):
    net = nn.Sequential(nn.Linear(in_dim, hidden), nn.Tanh(), nn.Linear(hidden, out_dim))
    if zero_output:
        nn.init.zeros_(net[-1].weight)
        nn.init.zeros_(net[-1].bias)
    return net


class ResidualModule(nn.Module):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__()
        self.pdim, self.rdim, self.cdim = dims
        self.use_context = use_context
        self.settings = {**DEFAULTS, **settings}


class PairwiseDeltaNet(ResidualModule):
    """Ordered-source-to-destination changes, averaged as an innovation on b16."""
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        self.net = mlp(2*self.pdim + self.cdim, self.settings["hidden"], self.rdim)

    def forward(self, b):
        cs, cq = context_values(b, self.use_context)
        e, n, s, _ = b["source_x"].shape
        p = b["p"][:, :, None].expand(-1, -1, s, -1)
        dc = (cq[:, None]-cs)[:, None].expand(-1, n, -1, -1)
        pair_delta = self.net(torch.cat([b["source_x"], p, dc], -1))
        pair_delta = pair_delta - pair_delta.mean(1, keepdim=True)
        return dict(mean=masked_mean(pair_delta, b["mask"]), pair_delta=pair_delta)


class LowRankHypernetwork(ResidualModule):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        k = min(self.settings["interaction_rank"], self.pdim, self.rdim)
        self.v = nn.Linear(self.pdim, k, bias=False)
        self.u = nn.Linear(k, self.rdim, bias=False)
        self.context_net = mlp(self.cdim, 8, k, zero_output=False)
        nn.init.zeros_(self.u.weight)

    def forward(self, b):
        _, cq = context_values(b, self.use_context)
        modulation = 1 + torch.tanh(self.context_net(cq))
        return dict(mean=self.u(self.v(b["p"]) * modulation[:, None]))


class ConditionalDeepSets(ResidualModule):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        h = self.settings["hidden"]
        self.encoder = mlp(self.pdim+2*self.cdim, h, h, zero_output=False)
        self.decoder = mlp(h+self.pdim, h, self.rdim)

    def forward(self, b):
        cs, cq = context_values(b, self.use_context)
        n = b["p"].shape[1]
        source_c = cs[:, None].expand(-1, n, -1, -1)
        dc = (cq[:, None]-cs)[:, None].expand(-1, n, -1, -1)
        encoded = self.encoder(torch.cat([b["source_x"], source_c, dc], -1))
        pooled = masked_mean(encoded, b["mask"])
        return dict(mean=self.decoder(torch.cat([pooled, b["p"]], -1)))


class MixtureOfExperts(ResidualModule):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        k = self.settings["experts"]
        self.experts = nn.Parameter(torch.empty(k, self.pdim, self.rdim))
        nn.init.normal_(self.experts, std=0.015)
        self.gate = nn.Linear(self.cdim, k)

    def forward(self, b):
        _, cq = context_values(b, self.use_context)
        weights = self.gate(cq).softmax(-1)
        predictions = torch.einsum("enp,kpr->enkr", b["p"], self.experts)
        mean = (predictions * weights[:, None, :, None]).sum(2)
        k = weights.shape[-1]
        balance = 0.02*k*((weights.mean(0)-1/k)**2).sum()
        entropy_penalty = 0.002*(weights*(weights.clamp_min(1e-8)*k).log()).sum(-1).mean()
        return dict(mean=mean, regularizer=balance+entropy_penalty, expert_weights=weights)


class LowRankTensor(ResidualModule):
    """CP-style perturbation x control-expression x response-component factors."""
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        k = min(self.settings["interaction_rank"], self.pdim, self.rdim)
        self.perturbation = nn.Linear(self.pdim, k, bias=False)
        self.context = nn.Linear(self.cdim, k, bias=False)
        self.response = nn.Linear(k, self.rdim, bias=False)
        nn.init.zeros_(self.response.weight)

    def forward(self, b):
        _, cq = context_values(b, self.use_context)
        factor = 1 + self.context(cq)
        return dict(mean=self.response(self.perturbation(b["p"])*factor[:, None]))


class ConditionalNeuralProcess(ResidualModule):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        h = self.settings["hidden"]
        self.encoder = mlp(self.cdim+self.rdim, h, h, zero_output=False)
        self.decoder = mlp(h+self.cdim+self.pdim, h, 2*self.rdim)

    def forward(self, b):
        cs, cq = context_values(b, self.use_context)
        n = b["p"].shape[1]
        encoded = self.encoder(torch.cat([cs[:, None].expand(-1,n,-1,-1), b["source_y"]], -1))
        pooled = masked_mean(encoded, b["mask"])
        query = torch.cat([pooled, cq[:, None].expand(-1,n,-1), b["p"]], -1)
        mean, raw_variance = self.decoder(query).chunk(2, dim=-1)
        variance = 0.02 + F.softplus(raw_variance)
        return dict(mean=mean, variance=variance)


class LearnedMetricTransfer(ResidualModule):
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        k = self.settings["interaction_rank"]
        self.metric = nn.Linear(self.cdim, k, bias=False)
        self.perturbation = nn.Linear(self.pdim, k)

    def forward(self, b):
        cs, cq = context_values(b, self.use_context)
        direction = self.metric(cq[:, None]-cs)
        sensitivity = 0.1 + F.softplus(self.perturbation(b["p"]))
        distance = (sensitivity[:, :, None]*direction[:, None].square()).mean(-1)
        weights = (-distance).masked_fill(~b["mask"][:, None].bool(), -1e9).softmax(-1)
        weighted = (weights[:, :, :, None]*b["source_y"]).sum(2)
        # Uniform weights give exactly zero innovation: b16 remains the anchor.
        mean = weighted - masked_mean(b["source_y"], b["mask"])
        return dict(mean=mean, source_weights=weights)


class SharedKernelGP(ResidualModule):
    """Exact tiny GP over source contexts; three parameters shared across targets/factors."""
    def __init__(self, dims, use_context=True, **settings):
        super().__init__(dims, use_context, **settings)
        self.log_length = nn.Parameter(torch.tensor(0.0))
        self.log_noise = nn.Parameter(torch.tensor(-1.0))
        self.log_amplitude = nn.Parameter(torch.tensor(0.0))

    def forward(self, b):
        cs, cq = context_values(b, self.use_context)
        length = self.log_length.clamp(-3, 3).exp()
        noise = 0.01 + self.log_noise.clamp(-5, 3).exp()
        amplitude = 0.01 + self.log_amplitude.clamp(-5, 3).exp()
        dist = ((cs[:, :, None]-cs[:, None, :])**2).sum(-1)
        k = amplitude*torch.exp(-dist/(2*length**2))
        mask = b["mask"]
        k = k*mask[:, :, None]*mask[:, None, :]
        k = k + torch.diag_embed(mask*(noise+1e-5)+(1-mask))
        query_dist = ((cq[:, None]-cs)**2).sum(-1)
        kq = amplitude*torch.exp(-query_dist/(2*length**2))*mask
        chol = torch.linalg.cholesky(k)
        weights = torch.cholesky_solve(kq[:, :, None], chol).squeeze(-1)
        source_dev = b["source_y"] - masked_mean(b["source_y"], mask)[:, :, None]
        mean = (weights[:, None, :, None]*source_dev).sum(2)
        posterior = (amplitude-(weights*kq).sum(-1)).clamp_min(1e-8)
        variance = (posterior+noise)[:, None, None].expand_as(mean)
        return dict(mean=mean, variance=variance, posterior_variance=posterior,
                    observation_noise=noise, source_weights=weights,
                    kernel_parameters=torch.stack([length, noise, amplitude]))


MODULES = dict(pairwise=PairwiseDeltaNet, hypernetwork=LowRankHypernetwork,
               deepsets=ConditionalDeepSets, moe=MixtureOfExperts, tensor=LowRankTensor,
               cnp=ConditionalNeuralProcess, metric=LearnedMetricTransfer, gp=SharedKernelGP)


@dataclass
class Fitted:
    family: str
    module: ResidualModule
    history: list
    seconds: float
    seed: int
    shuffle_labels: bool
    shuffle_contexts: bool

    @property
    def params(self):
        return sum(p.numel() for p in self.module.parameters())

    def predict(self, batch):
        self.module.eval()
        with torch.no_grad():
            raw = self.module(batch)
            raw["mean"] = center_correction(raw["mean"])
            if "variance" in raw:
                # Marginal variance after removing the all-target mean under an
                # independent-target diagonal approximation. Baseline uncertainty
                # and gene directions outside the output basis are not modeled.
                n = raw["mean"].shape[1]
                v = raw["variance"]
                raw["variance"] = v*(1-2/n)+v.sum(1,keepdim=True)/(n*n)
            return {k: v.detach().numpy() for k,v in raw.items() if torch.is_tensor(v)}

    def save(self, path):
        torch.save(dict(family=self.family, state_dict=self.module.state_dict(),
                        dims=(self.module.pdim,self.module.rdim,self.module.cdim),
                        use_context=self.module.use_context, settings=self.module.settings,
                        history=self.history, seconds=self.seconds, seed=self.seed,
                        shuffle_labels=self.shuffle_labels, shuffle_contexts=self.shuffle_contexts,
                        params=self.params), path)

    @classmethod
    def load(cls, path):
        a = torch.load(path, map_location="cpu", weights_only=False)
        module = MODULES[a["family"]](a["dims"], a["use_context"], **a["settings"])
        module.load_state_dict(a["state_dict"])
        return cls(a["family"],module,a["history"],a["seconds"],a["seed"],
                   a["shuffle_labels"],a["shuffle_contexts"])


def fit(prepared, family, use_context=True, seed=0, shuffle_labels=False,
        shuffle_contexts=False, **settings):
    start = time.monotonic()
    cfg = {**DEFAULTS, **settings}
    torch.manual_seed(seed)
    module = MODULES[family](prepared.dims, use_context, **cfg)
    batch, labels = prepared.training_inputs(shuffle_labels, shuffle_contexts)
    optimizer = torch.optim.AdamW(module.parameters(), lr=cfg["lr"],
                                  weight_decay=cfg["weight_decay"])
    history = []
    for step in range(cfg["steps"]):
        module.train()
        optimizer.zero_grad(set_to_none=True)
        out = module(batch)
        mean = center_correction(out["mean"])
        error = mean[:, None] - labels
        mse = error.square().mean()
        magnitude = mean.square().mean()
        loss = mse + cfg["residual_penalty"]*magnitude
        if "pair_delta" in out:
            # Directly supervise every valid ordered source -> destination pair.
            # Shuffled residual labels also replace the destinations in this loss.
            true_destination = labels + batch["baseline"][:, None]
            true_difference = true_destination[:, :, :, None] - batch["source_y"][:, None]
            pair_error = out["pair_delta"][:, None]-true_difference
            pair_mse = (pair_error.square().mean((1,2,4))*batch["mask"]).sum()/batch["mask"].sum()
            loss = loss + 0.5*pair_mse
        if "variance" in out:
            # Shared squared-error objective for means; additionally fit a diagonal
            # uncertainty model to replicate residuals (not a claim of calibration).
            variance = out["variance"][:, None]
            uncertainty = 0.5*(error.detach().square()/variance + variance.log()).mean()
            loss = loss + 0.05*uncertainty
        loss = loss + out.get("regularizer", 0)
        if not torch.isfinite(loss):
            raise FloatingPointError(f"{family} nonfinite loss at step {step}")
        loss.backward()
        grad_norm = nn.utils.clip_grad_norm_(module.parameters(), 5.0)
        optimizer.step()
        if step == 0 or (step+1)%40 == 0 or step+1 == cfg["steps"]:
            history.append(dict(step=step+1, loss=float(loss.detach()), mse=float(mse.detach()),
                                residual_energy=float(magnitude.detach()),
                                grad_norm=float(grad_norm)))
    return Fitted(family,module,history,time.monotonic()-start,seed,
                  shuffle_labels,shuffle_contexts)


class FastPDS:
    """Exact cosine-rank algebra for baseline + low-rank correction; scoring only."""
    def __init__(self, prepared, truth):
        self.truth_norm = np.linalg.norm(truth,axis=1)
        self.base_dot = prepared.baseline @ truth.T
        self.correction_dot = prepared.scale_y * (prepared.basis_y.T @ truth.T)
        self.base_energy = np.sum(prepared.baseline**2,axis=1)
        self.base_projection = prepared.scale_y*(prepared.baseline @ prepared.basis_y)
        self.scale_y = prepared.scale_y

    def score(self, correction, strength):
        from jiang_benchmark import ranks_from_similarity
        c = np.asarray(correction,dtype=np.float64)*float(strength)
        dot = self.base_dot + c @ self.correction_dot
        energy = self.base_energy+2*np.sum(c*self.base_projection,axis=1)+self.scale_y**2*np.sum(c*c,axis=1)
        norms = np.sqrt(np.maximum(energy,0))[:,None]*self.truth_norm[None]
        sim = np.divide(dot,norms,out=np.zeros_like(dot),where=norms>0)
        return float(ranks_from_similarity(sim).mean())
