"""Conditional perturbation x context architectures on one shared interface, for the screen in
`arch_screen.py`. Every model consumes a batch dict with
    zc  (B, dc)   context vector: training-only PCA of the pooled control profile
    tok (B, T, d) context tokens: the same profile cut into T gene blocks (for attention/set models)
    zp  (B, 32)   perturbation programme: PCA of the refit prediction for the target
    pf  (B, 64)   frozen ESM-2 gene features of the target (for joint models that must not see B(g))
    x0  (B, K)    base state in the full-response basis (= main effect mu), for flow matching
and emits a latent that a frozen basis decodes to 5,592 genes:
    RESIDUAL family  -> rank-32 residual latent, prediction = base + alpha * residual
    JOINT family     -> rank-K full-response latent, prediction = decoded response (base unused)
Heads are zero-initialised where that makes "untrained == base" exact.
"""
from dataclasses import dataclass
import math
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F


def mlp(i, h, o, depth=2, act=nn.GELU, ln=True):
    layers, d = [], i
    for _ in range(depth - 1):
        layers += [nn.Linear(d, h)] + ([nn.LayerNorm(h)] if ln else []) + [act()]; d = h
    layers.append(nn.Linear(d, o))
    return nn.Sequential(*layers)


def zero_last(seq):
    last = [m for m in seq.modules() if isinstance(m, nn.Linear)][-1]
    nn.init.zeros_(last.weight); nn.init.zeros_(last.bias)
    return seq


class Base(nn.Module):
    family = 'residual'          # or 'joint'
    def loss(self, b, target):
        return F.mse_loss(self(b), target)
    def predict(self, b):
        return self(b)


# ---------------------------------------------------------------- residual family
class Additive(Base):
    """context -> residual, no perturbation input: the pure generic per-context shift."""
    def __init__(self, dc, h=32, r=32):
        super().__init__(); self.f = zero_last(mlp(dc, h, r))
    def forward(self, b): return self.f(b['zc'])


class FiLM(Base):
    """the current model: h = (1 + gamma(u_c)) * u_p, multiplicative only."""
    def __init__(self, dc, h=32, r=32):
        super().__init__()
        self.uc = mlp(dc, h, h, 1); self.up = nn.Sequential(nn.Linear(32, h), nn.LayerNorm(h), nn.GELU())
        self.gamma = nn.Linear(h, h); nn.init.zeros_(self.gamma.weight); nn.init.zeros_(self.gamma.bias)
        self.head = nn.Linear(h, r); nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
    def forward(self, b):
        uc = F.gelu(F.layer_norm(self.uc(b['zc']), (self.uc[-1].out_features,)))
        return self.head((1 + self.gamma(uc)) * self.up(b['zp']))


class ConcatMLP(Base):
    """F([h_g, c]) plain conditional MLP on the concatenation."""
    def __init__(self, dc, h=64, r=32):
        super().__init__(); self.f = zero_last(mlp(dc + 32, h, r, depth=3))
    def forward(self, b): return self.f(torch.cat([b['zp'], b['zc']], -1))


class FullFiLM(Base):
    """two FiLM blocks with gamma AND beta, on a 2-layer perturbation trunk."""
    def __init__(self, dc, h=64, r=32):
        super().__init__()
        self.c = mlp(dc, h, h, 2); self.l1 = nn.Linear(32, h); self.l2 = nn.Linear(h, h)
        self.g1 = nn.Linear(h, 2 * h); self.g2 = nn.Linear(h, 2 * h)
        for g in (self.g1, self.g2): nn.init.zeros_(g.weight); nn.init.zeros_(g.bias)
        self.head = nn.Linear(h, r); nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
    def forward(self, b):
        c = self.c(b['zc']); x = F.gelu(self.l1(b['zp']))
        g, be = self.g1(c).chunk(2, -1); x = F.gelu((1 + g) * x + be)
        x = F.gelu(self.l2(x)); g, be = self.g2(c).chunk(2, -1); x = (1 + g) * x + be
        return self.head(x)


class Gated(Base):
    """h = sigmoid(Wc c) * (Wp p): multiplicative gate on the perturbation representation."""
    def __init__(self, dc, h=64, r=32):
        super().__init__(); self.g = mlp(dc, h, h, 2); self.p = mlp(32, h, h, 2)
        self.head = nn.Linear(h, r); nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
    def forward(self, b): return self.head(torch.sigmoid(self.g(b['zc'])) * self.p(b['zp']))


class Bilinear(Base):
    """r_k = c^T A_k p over rank-r factorised A_k plus linear terms."""
    def __init__(self, dc, r=32, f=16):
        super().__init__(); self.U = nn.Linear(dc, r * f, bias=False); self.Vp = nn.Linear(32, r * f, bias=False)
        self.lin_p = nn.Linear(32, r); self.r, self.f = r, f
        nn.init.zeros_(self.lin_p.weight); nn.init.zeros_(self.lin_p.bias)
        self.scale = nn.Parameter(torch.zeros(1))
    def forward(self, b):
        u = self.U(b['zc']).view(-1, self.r, self.f); v = self.Vp(b['zp']).view(-1, self.r, self.f)
        return self.scale * (u * v).sum(-1) + self.lin_p(b['zp'])


class CrossAttn(Base):
    """perturbation queries attend over the context's gene-block tokens; multiplicative interaction."""
    def __init__(self, dc, dtok, h=64, heads=4, r=32):
        super().__init__(); self.h, self.heads = h, heads
        self.wq = nn.Linear(32, h, bias=False); self.wk = nn.Linear(dtok, h, bias=False); self.wv = nn.Linear(dtok, h, bias=False)
        self.wr = nn.Linear(h, h); self.wp = nn.Linear(32, h)
        self.head = nn.Linear(h, r); nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
    def forward(self, b):
        B, T, _ = b['tok'].shape; hd = self.h // self.heads
        q = self.wq(b['zp']).view(B, self.heads, 1, hd); k = self.wk(b['tok']).view(B, T, self.heads, hd).transpose(1, 2)
        v = self.wv(b['tok']).view(B, T, self.heads, hd).transpose(1, 2)
        w = torch.softmax(q @ k.transpose(-1, -2) / math.sqrt(hd), -1); c = (w @ v).reshape(B, self.h)
        return self.head(self.wr(c) * self.wp(b['zp']))


class LoRAHyper(Base):
    """context-conditioned low-rank adapter: W(c) = W0 + A(c) B on the perturbation map."""
    def __init__(self, dc, h=64, r=32, lr_rank=4):
        super().__init__(); self.W0 = nn.Linear(32, h); self.A = nn.Linear(dc, h * lr_rank); self.Bm = nn.Parameter(torch.randn(lr_rank, 32) / 32)
        nn.init.zeros_(self.A.weight); nn.init.zeros_(self.A.bias); self.k = lr_rank
        self.out = zero_last(mlp(h, h, r, 2))
    def forward(self, b):
        A = self.A(b['zc']).view(-1, self.W0.out_features, self.k)
        delta = torch.einsum('bhk,kp,bp->bh', A, self.Bm, b['zp'])
        return self.out(F.gelu(self.W0(b['zp']) + delta))


class FullHyper(Base):
    """context generates ALL weights of a 32->32->32 perturbation MLP."""
    def __init__(self, dc, h=32, r=32):
        super().__init__(); self.h = h; n = 32 * h + h + h * r + r
        self.gen = mlp(dc, 64, n, 2); nn.init.zeros_(self.gen[-1].weight); nn.init.zeros_(self.gen[-1].bias)
        self.bias0 = nn.Parameter(torch.zeros(n))
    def forward(self, b):
        w = self.gen(b['zc']) + self.bias0; h = self.h; B = w.shape[0]; i = 0
        W1 = w[:, i:i + 32 * h].view(B, 32, h); i += 32 * h; b1 = w[:, i:i + h]; i += h
        W2 = w[:, i:i + h * 32].view(B, h, 32); i += h * 32; b2 = w[:, i:]
        x = F.gelu(torch.einsum('bp,bph->bh', b['zp'], W1) + b1)
        return torch.einsum('bh,bhr->br', x, W2) + b2


class MoE(Base):
    """softmax router over context, K linear experts on the perturbation representation."""
    def __init__(self, dc, h=64, r=32, K=8):
        super().__init__(); self.router = mlp(dc, 32, K, 2); self.p = mlp(32, h, h, 2)
        self.experts = nn.ModuleList([nn.Linear(h, r) for _ in range(K)])
        for e in self.experts: nn.init.zeros_(e.weight); nn.init.zeros_(e.bias)
    def forward(self, b):
        w = torch.softmax(self.router(b['zc']), -1); u = self.p(b['zp'])
        return (w.unsqueeze(-1) * torch.stack([e(u) for e in self.experts], 1)).sum(1)


class TokenTransformer(Base):
    """set-to-vector Transformer over [perturbation token, context gene-block tokens]; proxy for
    a STATE-style set model (per-cell control sets exist for only 7/38 contexts and 5 of those
    files are the scrambled ones, so tokens are gene blocks of the pooled control profile)."""
    def __init__(self, dc, dtok, h=64, heads=4, layers=2, r=32):
        super().__init__(); self.pe = nn.Linear(32, h); self.te = nn.Linear(dtok, h)
        self.pos = nn.Parameter(torch.randn(1, 64, h) * 0.02)
        enc = nn.TransformerEncoderLayer(h, heads, 2 * h, dropout=0.0, batch_first=True, norm_first=True)
        self.tr = nn.TransformerEncoder(enc, layers); self.head = nn.Linear(h, r)
        nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
    def forward(self, b):
        x = torch.cat([self.pe(b['zp']).unsqueeze(1), self.te(b['tok'])], 1)
        x = x + self.pos[:, :x.shape[1]]; x = self.tr(x)
        return self.head(x[:, 0])


# ---------------------------------------------------------------- joint family
class JointMLP(Base):
    """full response from [ESM gene features, context]: never sees B(g)."""
    family = 'joint'
    def __init__(self, dc, K, h=256):
        super().__init__(); self.f = mlp(64 + dc, h, K, depth=3)
    def forward(self, b): return self.f(torch.cat([b['pf'], b['zc']], -1))


class JointMLPProg(Base):
    """full response from [refit programme, context]: sees B(g) as an input feature."""
    family = 'joint'
    def __init__(self, dc, K, h=128):
        super().__init__(); self.f = mlp(32 + dc, h, K, depth=3)
    def forward(self, b): return self.f(torch.cat([b['zp'], b['zc']], -1))


class CPA(Base):
    """CPA-style: z = enc_pert(p) + enc_ctx(c), decode; additive latent composition."""
    family = 'joint'
    def __init__(self, dc, K, h=64):
        super().__init__(); self.ep = mlp(32, h, h, 2); self.ec = mlp(dc, h, h, 2); self.dec = mlp(h, 128, K, 2)
    def forward(self, b): return self.dec(self.ep(b['zp']) + self.ec(b['zc']))


class Biolord(Base):
    """biolord-style: separate perturbation and context latents, decorrelation penalty, decoder
    on the concatenation."""
    family = 'joint'
    def __init__(self, dc, K, h=32):
        super().__init__(); self.ep = mlp(32, 64, h, 2); self.ec = mlp(dc, 64, h, 2); self.dec = mlp(2 * h, 128, K, 2)
    def latents(self, b): return self.ep(b['zp']), self.ec(b['zc'])
    def forward(self, b): zp, zc = self.latents(b); return self.dec(torch.cat([zp, zc], -1))
    def loss(self, b, target):
        zp, zc = self.latents(b); zpn = zp - zp.mean(0); zcn = zc - zc.mean(0)
        cov = (zpn.T @ zcn) / max(zp.shape[0] - 1, 1)
        return F.mse_loss(self.dec(torch.cat([zp, zc], -1)), target) + 0.1 * (cov ** 2).mean()


class CondFlow(Base):
    """conditional affine-coupling flow for y | (c, p) in the K-dim response basis; NLL training,
    prediction = inverse transform of z = 0 (the conditional mode/median)."""
    family = 'joint'
    def __init__(self, dc, K, h=64, layers=4):
        super().__init__(); self.K = K; self.cond = mlp(32 + dc, h, h, 2)
        self.nets = nn.ModuleList([mlp(K // 2 + h, h, 2 * (K - K // 2), 2) for _ in range(layers)])
        for n in self.nets: nn.init.zeros_(n[-1].weight); nn.init.zeros_(n[-1].bias)
    def _st(self, i, x1, c): s, t = self.nets[i](torch.cat([x1, c], -1)).chunk(2, -1); return torch.tanh(s), t
    def fwd(self, y, c):      # y -> z, logdet
        z, ld, k = y, 0.0, self.K // 2
        for i in range(len(self.nets)):
            a, bb = (z[:, :k], z[:, k:]) if i % 2 == 0 else (z[:, k:], z[:, :k])
            s, t = self._st(i, a, c); bb = bb * torch.exp(s) + t; ld = ld + s.sum(-1)
            z = torch.cat([a, bb], -1) if i % 2 == 0 else torch.cat([bb, a], -1)
        return z, ld
    def inv(self, z, c):
        k = self.K // 2
        for i in reversed(range(len(self.nets))):
            a, bb = (z[:, :k], z[:, k:]) if i % 2 == 0 else (z[:, k:], z[:, :k])
            s, t = self._st(i, a, c); bb = (bb - t) * torch.exp(-s)
            z = torch.cat([a, bb], -1) if i % 2 == 0 else torch.cat([bb, a], -1)
        return z
    def loss(self, b, target):
        c = self.cond(torch.cat([b['zp'], b['zc']], -1)); z, ld = self.fwd(target, c)
        return (0.5 * (z ** 2).sum(-1) - ld).mean() / self.K
    def predict(self, b):
        c = self.cond(torch.cat([b['zp'], b['zc']], -1)); return self.inv(torch.zeros(b['zp'].shape[0], self.K, device=c.device), c)
    def forward(self, b): return self.predict(b)


class FlowMatch(Base):
    """conditional flow matching from the base state x0 (main-effect latent) to the perturbed
    state x1 (observed latent); v(x_t, t | c, p) MLP; prediction integrates 16 Euler steps."""
    family = 'joint'
    def __init__(self, dc, K, h=128):
        super().__init__(); self.K = K; self.v = mlp(K + 1 + 32 + dc, h, K, 3)
    def field(self, x, t, b): return self.v(torch.cat([x, t, b['zp'], b['zc']], -1))
    def loss(self, b, target):
        x0 = b['x0']; t = torch.rand(x0.shape[0], 1, device=x0.device); xt = (1 - t) * x0 + t * target
        return F.mse_loss(self.field(xt, t, b), target - x0)
    def predict(self, b):
        x = b['x0']; n = 16
        for i in range(n):
            t = torch.full((x.shape[0], 1), i / n, device=x.device); x = x + self.field(x, t, b) / n
        return x
    def forward(self, b): return self.predict(b)


def build(name, dc, dtok, K):
    R = dict(additive=lambda: Additive(dc), film=lambda: FiLM(dc), concat_mlp=lambda: ConcatMLP(dc),
             full_film=lambda: FullFiLM(dc), gated=lambda: Gated(dc), bilinear=lambda: Bilinear(dc),
             cross_attn=lambda: CrossAttn(dc, dtok), lora_hyper=lambda: LoRAHyper(dc), full_hyper=lambda: FullHyper(dc),
             moe=lambda: MoE(dc), token_transformer=lambda: TokenTransformer(dc, dtok),
             joint_mlp=lambda: JointMLP(dc, K), joint_mlp_prog=lambda: JointMLPProg(dc, K), cpa=lambda: CPA(dc, K),
             biolord=lambda: Biolord(dc, K), cond_flow=lambda: CondFlow(dc, K), flow_match=lambda: FlowMatch(dc, K))
    return R[name]()


ARCHS = ['additive', 'film', 'concat_mlp', 'full_film', 'gated', 'bilinear', 'cross_attn', 'lora_hyper', 'full_hyper',
         'moe', 'token_transformer', 'joint_mlp', 'joint_mlp_prog', 'cpa', 'biolord', 'cond_flow', 'flow_match']


def train(model, T, V, device, lr=2e-3, wd=1e-3, batch=256, min_steps=2000, patience_evals=30, max_steps=20000,
          eval_every=100, seed=0, log=None, aug=lambda b: b):
    """generous early stopping: no stop before min_steps, patience 30 evals (3,000 steps), max 20k."""
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    model.to(device); opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    n = len(T['y']); best, best_step, bad, curve = float('inf'), 0, 0, []
    best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    def evaluate(D):
        model.eval(); out = 0.0
        with torch.no_grad():
            for s in range(0, len(D['y']), 4096):
                sl = slice(s, s + 4096); b = aug({k: v[sl] for k, v in D.items() if k != 'y'})
                out += float(((model.predict(b) - D['y'][sl]) ** 2).mean(1).sum())
        model.train(); return out / len(D['y'])
    zero_val = float((V['y'] ** 2).mean()) if model.family == 'residual' else float(((V['x0'] - V['y']) ** 2).mean())
    for step in range(1, max_steps + 1):
        i = torch.as_tensor(rng.integers(0, n, size=batch), device=device)
        b = aug({k: v[i] for k, v in T.items() if k != 'y'})
        loss = model.loss(b, T['y'][i]); opt.zero_grad(set_to_none=True); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
        if step % eval_every == 0:
            v = evaluate(V); curve.append(dict(step=step, train_loss=float(loss), val=v))
            if log: log(curve[-1])
            if v < best - 1e-9: best, best_step, bad = v, step, 0; best_state = {k: q.detach().clone() for k, q in model.state_dict().items()}
            else:
                bad += 1
                if step >= min_steps and bad >= patience_evals: break
    model.load_state_dict(best_state)
    return dict(best_val=best, best_step=best_step, zero_val=zero_val, beat_zero=best < zero_val, steps_run=step,
                params=int(sum(p.numel() for p in model.parameters())), curve=curve)
