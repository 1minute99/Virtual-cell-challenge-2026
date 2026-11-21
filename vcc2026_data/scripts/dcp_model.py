"""The three small interaction families (Part 6) plus the shared trunk / adapter plumbing.

Shape of every model:

    z_context (64, shared space)  --> ContextAdapter  --> u_c (hidden)        [SHARED]
    z_pert    (rank, per dataset) --> PertAdapter[ds] --> u_p (hidden)        [per dataset]
                                        |
                                   interaction trunk                          [SHARED]
                                        |
                                    Head[ds]  --> latent (rank) --> @ basis_ds.T --> genes

The trunk is the only thing transferred between drugs and genes; everything touching a
dataset-specific coordinate system is an adapter. The final Head is zero-initialised so an
untrained model emits exactly zero residual, which is what makes `alpha=0 == refit` exact.
"""
from dataclasses import dataclass, field
import numpy as np
import torch
import torch.nn as nn


@dataclass
class Spec:
    kind: str = 'film'            # film | operator | mixture
    ctx_dim: int = 64
    hidden: int = 32
    rank: int = 16
    experts: int = 4              # mixture only
    operators: int = 4            # operator only: K
    op_rank: int = 4              # operator only: rank of each A_k
    weight_decay: float = 1e-3
    lr: float = 3e-3
    steps: int = 4000
    batch: int = 256
    patience: int = 8
    eval_every: int = 100
    seed: int = 20260920


class ContextAdapter(nn.Module):
    def __init__(self, ctx_dim, hidden):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(ctx_dim, hidden), nn.LayerNorm(hidden), nn.GELU())

    def forward(self, z):
        return self.net(z)


class PertAdapter(nn.Module):
    def __init__(self, rank, hidden):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(rank, hidden), nn.LayerNorm(hidden), nn.GELU())

    def forward(self, z):
        return self.net(z)


class FiLMTrunk(nn.Module):
    """h = (1 + gamma(u_c)) * u_p. Multiplicative only -- no additive beta by design."""
    def __init__(self, spec):
        super().__init__()
        self.gamma = nn.Linear(spec.hidden, spec.hidden)
        nn.init.zeros_(self.gamma.weight); nn.init.zeros_(self.gamma.bias)

    def forward(self, uc, up):
        return (1.0 + self.gamma(uc)) * up


class OperatorTrunk(nn.Module):
    """z_out = [I + sum_k phi_k(u_c) U_k V_k^T] u_p, each A_k low rank."""
    def __init__(self, spec):
        super().__init__()
        K, h, r = spec.operators, spec.hidden, spec.op_rank
        self.phi = nn.Linear(h, K)
        nn.init.zeros_(self.phi.weight); nn.init.zeros_(self.phi.bias)
        self.U = nn.Parameter(torch.randn(K, h, r) * (1.0 / np.sqrt(h)))
        self.V = nn.Parameter(torch.randn(K, h, r) * (1.0 / np.sqrt(h)))

    def forward(self, uc, up):
        w = self.phi(uc)                                   # (B,K)
        t = torch.einsum('bh,khr->bkr', up, self.V)        # (B,K,r)
        d = torch.einsum('bkr,khr->bkh', t, self.U)        # (B,K,h)
        return up + (w.unsqueeze(-1) * d).sum(1)


class MixtureTrunk(nn.Module):
    """h = sum_k softmax(Router(u_c))_k * Expert_k(u_p): interpolate among context archetypes."""
    def __init__(self, spec):
        super().__init__()
        self.router = nn.Linear(spec.hidden, spec.experts)
        self.experts = nn.ModuleList([nn.Linear(spec.hidden, spec.hidden)
                                      for _ in range(spec.experts)])

    def forward(self, uc, up):
        w = torch.softmax(self.router(uc), -1)
        e = torch.stack([m(up) for m in self.experts], 1)
        return (w.unsqueeze(-1) * e).sum(1)


TRUNKS = dict(film=FiLMTrunk, operator=OperatorTrunk, mixture=MixtureTrunk)


class InteractionModel(nn.Module):
    def __init__(self, spec, datasets):
        super().__init__()
        self.spec = spec
        self.ctx = ContextAdapter(spec.ctx_dim, spec.hidden)
        self.pert = nn.ModuleDict({d: PertAdapter(spec.rank, spec.hidden) for d in datasets})
        self.trunk = TRUNKS[spec.kind](spec)
        self.head = nn.ModuleDict({d: nn.Linear(spec.hidden, spec.rank) for d in datasets})
        for d in datasets:
            nn.init.zeros_(self.head[d].weight)
            nn.init.zeros_(self.head[d].bias)

    def forward(self, ds, zc, zp):
        return self.head[ds](self.trunk(self.ctx(zc), self.pert[ds](zp)))

    def trunk_parameters(self):
        return list(self.ctx.parameters()) + list(self.trunk.parameters())

    def adapter_parameters(self, ds=None):
        keys = list(self.pert.keys()) if ds is None else [ds]
        out = []
        for k in keys:
            out += list(self.pert[k].parameters()) + list(self.head[k].parameters())
        return out


def n_params(m):
    return int(sum(p.numel() for p in m.parameters() if p.requires_grad))


def add_dataset(model, name, device):
    """Attach a fresh adapter/head pair for a new dataset to an already-pretrained model."""
    s = model.spec
    model.pert[name] = PertAdapter(s.rank, s.hidden).to(device)
    h = nn.Linear(s.hidden, s.rank)
    nn.init.zeros_(h.weight); nn.init.zeros_(h.bias)
    model.head[name] = h.to(device)
    return model


# ------------------------------------------------------------------ training

class Batches:
    """One task = one dataset. Sampling is balanced across tasks by construction."""
    def __init__(self, tasks, rng):
        self.tasks = tasks          # name -> dict(zc=(N,ctx), zp=(N,rank), y=(N,rank), w=(N,), g=(N,))
        self.rng = rng

    def draw(self, name, size):
        t = self.tasks[name]
        i = self.rng.integers(0, len(t['y']), size=size)
        return i


def _j(v):
    """NaN/inf -> None, so every record here is writable with json allow_nan=False."""
    if v is None:
        return None
    v = float(v)
    return v if np.isfinite(v) else None


def _to(d, device):
    return {k: (torch.as_tensor(v, device=device) if isinstance(v, np.ndarray) else v)
            for k, v in d.items()}


def train(model, tasks, spec, device, weights=None, val=None, log=None,
          trunk_lr_scale=1.0, freeze_trunk=False, max_steps=None):
    """tasks: name -> dict(zc, zp, y, w). val: name -> same, held-out CONTEXTS only.

    weights: task sampling probabilities. Defaults to uniform over tasks (balanced), which is
    what keeps Tahoe from swamping GSE purely because it has 8x more conditions.
    """
    model.to(device)
    T = {k: _to(v, device) for k, v in tasks.items()}
    V = {k: _to(v, device) for k, v in (val or {}).items()}
    names = list(T)
    p = np.array([weights.get(n, 1.0) for n in names] if weights else [1.0] * len(names), float)
    p = p / p.sum()
    if freeze_trunk:
        for q in model.trunk_parameters():
            q.requires_grad_(False)
    groups = [dict(params=[q for q in model.adapter_parameters() if q.requires_grad], lr=spec.lr)]
    tp = [q for q in model.trunk_parameters() if q.requires_grad]
    if tp:
        groups.append(dict(params=tp, lr=spec.lr * trunk_lr_scale))
    opt = torch.optim.AdamW(groups, lr=spec.lr, weight_decay=spec.weight_decay)
    rng = np.random.default_rng(spec.seed)
    torch.manual_seed(spec.seed)
    steps = max_steps or spec.steps
    best, best_step, bad, curve = float('inf'), 0, 0, []
    best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}

    def evaluate(D):
        if not D:
            return float('nan')
        model.eval()
        tot, den = 0.0, 0.0
        with torch.no_grad():
            for n, t in D.items():
                pr = model(n, t['zc'], t['zp'])
                e = ((pr - t['y']) ** 2).mean(1) * t['w']
                tot += float(e.sum()); den += float(t['w'].sum())
        model.train()
        return tot / max(den, 1e-9)

    zero_val = None
    if V:
        zero_val = evaluate(V)
    model.train()
    for step in range(1, steps + 1):
        n = names[rng.choice(len(names), p=p)]
        t = T[n]
        i = torch.as_tensor(rng.integers(0, len(t['y']), size=min(spec.batch, len(t['y']))),
                            device=device)
        pr = model(n, t['zc'][i], t['zp'][i])
        loss = (((pr - t['y'][i]) ** 2).mean(1) * t['w'][i]).sum() / t['w'][i].sum().clamp_min(1e-9)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_([q for q in model.parameters() if q.requires_grad], 1.0)
        opt.step()
        if step % spec.eval_every == 0:
            v = evaluate(V)
            tr = evaluate(T)
            # None rather than NaN: these dicts are written with json allow_nan=False
            entry = dict(step=step, train=_j(tr), val=_j(v))
            curve.append(entry)
            if log:
                log(dict(entry))
            if V and v < best - 1e-9:
                best, best_step, bad = v, step, 0
                best_state = {k: q.detach().clone() for k, q in model.state_dict().items()}
            else:
                bad += 1
                if V and bad >= spec.patience:
                    break
    if V:
        model.load_state_dict(best_state)
    return dict(best_val=_j(best), best_step=best_step, zero_residual_val=_j(zero_val),
                beat_zero=bool(zero_val is not None and best < zero_val),
                params=n_params(model), steps_run=step, curve=curve)
