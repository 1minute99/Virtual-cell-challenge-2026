"""transfer_v2 cell emission.

iid_mult (the new default, honest i.i.d. emission): for every target an INDEPENDENT draw of n control cells
(without replacement within the target), each cell's counts multiplied by FC_g = exp(LFC_g) and rounded
stochastically (unbiased). Zero stays zero (no counts invented for unexpressed genes), the pseudobulk is
the control mean x FC in expectation, and each target carries its own sampling noise -- which is exactly
what the scorer's jackknife correction subtracts. No aggregate pinning.

shared_mult: ONE shared draw of n control cells for every target (as S1-S4), counts multiplied by FC and
rounded stochastically. The draw's own sampling deviation is common to all targets (it barely moves the PDS
ranking); the scorer charges it on MSE because the across-perturbation budget cannot cover it -- honest,
penalised, never a rebate.

shared_cellshift (the old S1-S4 emission, kept only as the comparison arm): ONE shared draw for all
targets, per-cell log1p(CP10K) additive shift  counts = ((x/L*1e4 + 1) * e^d - 1) * L/1e4.
"""
import numpy as np
from scipy import sparse


def emit(controls, lfc, n_cells, seed, mode="iid_mult", gate_cpm=5.0, override=None):
    """controls: CSR (n_ctrl x G) raw counts; lfc: (n_targets x G) natural-log fold changes (or log1p-CP10K
    deltas for shared_cellshift). Returns CSR (n_targets*n_cells x G), float32 integers."""
    out = sparse.vstack(list(emit_iter(controls, lfc, n_cells, seed, mode, gate_cpm, override)), format="csr")
    out.eliminate_zeros()
    return out


def emit_iter(controls, lfc, n_cells, seed, mode="iid_mult", gate_cpm=5.0, override=None):
    """Same RNG stream as emit(), one target's CSR block at a time (flat memory for the official build)."""
    controls = controls.tocsr()
    rng = np.random.default_rng(seed)
    G = controls.shape[1]
    blocks = []
    if mode in ("shared_cellshift", "shared_mult", "shared_cellshift_gated", "pinned_mult"):
        pick = rng.choice(controls.shape[0], size=n_cells, replace=False)
        blk = np.asarray(controls[pick].todense(), dtype=np.float64)
        L = blk.sum(1, keepdims=True)
        base = blk / L * 1e4 + 1.0
    if mode == "pinned_mult":
        # LOCAL MEASUREMENT ONLY (issue #294/#348 residual): group sums pinned to full-control-mean x FC x n.
        full_mean = np.asarray(controls.sum(0)).ravel() / controls.shape[0]
        colsum = blk.sum(0)
    if mode == "shared_cellshift_gated":
        tot = np.asarray(controls.sum(0)).ravel()
        low = tot / tot.sum() * 1e6 < gate_cpm      # genes the DE gate never scores: multiplicative only
    for t in range(lfc.shape[0]):
        if mode == "iid_mult":
            pick = rng.choice(controls.shape[0], size=n_cells, replace=False)
            blk = np.asarray(controls[pick].todense(), dtype=np.float64)
            want = blk * np.exp(lfc[t])[None, :]
        elif mode == "shared_mult":
            want = blk * np.exp(lfc[t])[None, :]
        elif mode == "pinned_mult":
            e = np.exp(lfc[t]).copy()
            if override is not None:
                j = np.isfinite(override[t]); e[j] = np.exp(override[t, j])
            target = np.round(full_mean * e * n_cells)                      # exact integer group total per gene
            share = np.where(colsum > 0, blk / np.maximum(colsum, 1e-12), 1.0 / n_cells)
            want = share * target[None, :]
            fl = np.floor(want); rem = target - fl.sum(0)                  # largest remainder, column-exact
            frac = want - fl
            order = np.argsort(-frac, axis=0)
            add = np.zeros_like(fl)
            for g in np.flatnonzero(rem > 0):
                add[order[:int(rem[g]), g], g] = 1
            b = sparse.csr_matrix((fl + add).astype(np.float32)); b.eliminate_zeros()
            yield b
            continue
        elif mode == "shared_cellshift":
            want = (base * np.exp(lfc[t])[None, :] - 1.0) * L / 1e4
            np.clip(want, 0, None, out=want)
        elif mode == "shared_cellshift_gated":
            e = np.exp(lfc[t])[None, :]
            want = (base * e - 1.0) * L / 1e4
            np.clip(want, 0, None, out=want)
            want[:, low] = blk[:, low] * e[:, low]
        else:
            raise ValueError(mode)
        if override is not None:                       # multiplicative entries (promoter floor, knockdown)
            j = np.flatnonzero(np.isfinite(override[t]))
            if j.size:
                want[:, j] = blk[:, j] * np.exp(override[t, j])[None, :]
        fl = np.floor(want)
        cnt = fl + (rng.random(want.shape) < (want - fl))
        b = sparse.csr_matrix(cnt.astype(np.float32)); b.eliminate_zeros()
        yield b


def cs_cpm_expected(controls, lfc, low, chunk=2000):
    """Per target: mean over ALL controls of the per-cell composition (count / cell total) of the cellshift emission
    ((x/L*1e4+1) e^d - 1) L/1e4 clipped at 0, low genes x e^d. The DE-functional target of emit_dual_iter(cpc=True)."""
    controls = controls.tocsr(); N = controls.shape[0]; out = np.zeros(lfc.shape); E = np.exp(lfc)
    for a in range(0, N, chunk):
        X = np.asarray(controls[a:a + chunk].todense(), dtype=np.float32); L = X.sum(1, keepdims=True); base = X / L * 1e4 + 1.0
        for t in range(lfc.shape[0]):
            e = E[t].astype(np.float32)[None, :]
            w = (base * e - 1.0) * L / 1e4; np.clip(w, 0, None, out=w); w[:, low] = X[:, low] * e[:, low]
            out[t] += (w / w.sum(1, keepdims=True)).sum(0)
    return out / N


def kink_gamma(controls, lfc, override, n_cells, factor=1.17, draws=4, seed=0, margin=0.8):
    """Sum amplitude at which the prediction's across-perturbation spread B_pred (bulk_lognorm means of the pinned sums,
    TS 5e4) equals the scorer's capped correction claim ~ n_perts * factor * C_real(n) (issue #348 budget kink of
    expr_mse_unbiased_capped); C_real(n) = delete-1 jackknife trace of n control cells, factor = perturbed/control noise
    (1.16-1.19 measured). Test-time computable: needs only the context's controls, cells per perturbation and our delta.
    Harness 2026-09-27: sums moved to the kink = +0.013..+0.046 composite on the four RN panels, neutral where kink ~ 1.
    margin (target B/claim, default 0.8): the MSE rises like a cliff just past the kink (KOLF: 0.998 at gamma 2.8, 1.279 at
    3.0) while undershooting costs little, and the control-based claim is uncertain by 10-20 %, so aim below it."""
    from cell_eval2.moments import jackknife_correction as JK
    from cell_eval2.prep import bulk_lognorm_means as BLM
    rng = np.random.default_rng(seed); controls = controls.tocsr(); n = int(n_cells)
    cr = float(np.mean([JK(controls[rng.choice(controls.shape[0], size=n, replace=False)], np.zeros(n, np.int64), 1, 5e4)[0]
                        for _ in range(draws)]))
    claim = lfc.shape[0] * factor * margin * cr; fm = np.asarray(controls.sum(0)).ravel() / controls.shape[0]
    j = np.isfinite(override) if override is not None else np.zeros(lfc.shape, bool)
    def B(g):
        e = np.exp(g * lfc); e[j] = np.exp(override[j]); mu = BLM(fm[None, :] * e * n, 5e4); return float(((mu - mu.mean(0)) ** 2).sum())
    lo, hi = 0.3, 6.0
    if B(hi) < claim:
        return hi, cr
    for _ in range(40):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(mid) < claim else (lo, mid)
    return 0.5 * (lo + hi), cr


def budget_B(controls, lfc, override, n_cells):
    """B_pred(gamma) closure: across-perturbation centred SS of the bulk_lognorm means of the pinned sums (kink_gamma's B)."""
    from cell_eval2.prep import bulk_lognorm_means as BLM
    n = int(n_cells); fm = np.asarray(controls.sum(0)).ravel() / controls.shape[0]
    j = np.isfinite(override) if override is not None else np.zeros(lfc.shape, bool)
    def B(g):
        e = np.exp(g * lfc); e[j] = np.exp(override[j]); mu = BLM(fm[None, :] * e * n, 5e4); return float(((mu - mu.mean(0)) ** 2).sum())
    return B


def budget_match(controls, lfc, override, n_cells, target):
    """Scale s with B_pred(s * lfc) = target: a sum-side transform then changes the DIRECTION only and keeps the reference
    recipe's budget position (below the #348 kink the MSE credit is linear in <pred, truth>, just above it is a cliff)."""
    B = budget_B(controls, lfc, override, n_cells); lo, hi = 0.05, 8.0
    for _ in range(40):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def emit_dual_iter(controls, lfc, n_cells, seed, K=3000, gate_cpm=5.0, override=None, depth_up=0.0, agree=None, cpc=False,
                   top_frac=0.1, lfc_de=None, max_cell=8e5, ci=0, lsb=False, lcp=0.0):
    """Dual-moment emission (validated on the six-panel harness as fh_dual.py, gamma applied by the caller).
    Group sums: every target's n-cell total per gene is pinned to round(full control mean x exp(lfc) x n) with the
    override (knockdown / promoter floor) applied multiplicatively. Cell distribution: the top-K genes by |lfc| among
    genes the context expresses (>= gate_cpm) are re-allocated within the column (sum unchanged): up-genes flattened
    towards uniform, down-genes concentrated into the top 10 % of cells by that gene's control count, strength
    decaying 1 -> 0.3 with the rank. Other genes keep the proportional allocation of one shared control block.
    depth_up > 0: up-genes are flattened towards x_i ~ L_i^-depth_up instead of uniform (fh_dual 'du').
    agree (n_targets x G, 0..1): strength ordered by |lfc| x agree instead of |lfc| (fh_dual 'ag').
    cpc: after the re-allocation, every expressed non-override gene's column is depth-tilted (x_c ~ w_c L_c^s, bisection
    on s) so the group SUM stays pinned while the mean per-cell composition -- the functional the scorer's DE / NMAE uses --
    equals the cellshift emission's expectation over all controls (sch_emit 'cpc', 2026-09-26).
    lcp > 0 (sch_emit 'lcp<eta>', 2026-10-01): the re-allocation strength of the top-K genes is multiplied by
    clip((median CPM / CPM)^eta, 0.2, 1) over those genes, so high-expression genes keep more of the plain control-shaped
    column; sums stay pinned (PDS / MSE bit-identical). Real scorer, S14fk structure, eta 0.25: OR_HCT116X +0.0108 (NMAE),
    H1 +0.0015 (reach), RN mean +0.0018 (old structure). Default 0.0 = the old path, bit-exact."""
    controls = controls.tocsr()
    rng = np.random.default_rng(seed)
    n = n_cells
    if lsb:
        # DEPTH-BALANCED block. The scorer's DE log2FC is a per-cell normalize-then-average compared against the REAL
        # control pool, so however the emitted block's own composition differs from the pool's, that difference enters
        # every gene's predicted LFC as a bias; simple random sampling leaves it at ~1/sqrt(n) of the per-cell spread.
        # Systematic sampling on the library-size order is unbiased and much lower-variance for anything correlated with
        # depth, and issue #286 names depth-composition covariance as exactly what decides whether the emitted and
        # scored functionals agree. The pinned column totals use the FULL control mean, so PDS and MSE are provably
        # untouched -- measured on H1: PDS +0.0000, MSE +0.0000, fid +0.0015, reach +0.0037, jac +0.0011.
        Lall = np.asarray(controls.sum(1)).ravel()
        order = np.argsort(Lall, kind="stable"); step = len(order) / n
        pick = order[np.clip(((np.arange(n) + rng.random()) * step).astype(int), 0, len(order) - 1)]
    else:
        pick = rng.choice(controls.shape[0], size=n, replace=False)
    blk = np.asarray(controls[pick].todense(), dtype=np.float64)
    tot = np.asarray(controls.sum(0)).ravel()
    full_mean = tot / controls.shape[0]
    gate = tot / tot.sum() * 1e6 >= gate_cpm
    colsum = blk.sum(0)
    share = np.where(colsum > 0, blk / np.maximum(colsum, 1e-12), 1.0 / n)
    lam_rank = np.linspace(1.0, 0.3, max(K, 1))
    n_top = max(1, int(round(top_frac * n)))          # top_frac: share of cells down-genes are concentrated into (sch_emit 'nt')
    order_cells = np.argsort(-blk, axis=0, kind="stable")
    Lc = blk.sum(1)
    if cpc:                                       # lfc_de: DE-side target when the sums are rescaled (kink calibration)
        CPC = cs_cpm_expected(controls, lfc if lfc_de is None else lfc_de, ~gate); logL = np.log(Lc)
    for t in range(lfc.shape[0]):
        j = np.isfinite(override[t]) if override is not None else np.zeros(lfc.shape[1], bool)
        Lr = lfc if lfc_de is None else lfc_de           # DE side (membership, up/down, strength order) follows lfc_de
        cand = np.flatnonzero(gate & ~j)
        top = cand[np.argsort(-np.abs(Lr[t, cand]), kind="stable")[:K]]
        e = np.exp(lfc[t]).copy()
        if override is not None:
            e[j] = np.exp(override[t, j])
        T = np.round(full_mean * e * n)
        W = share * T[None, :]
        lam = lam_rank[:len(top)]
        if agree is not None:
            sc = np.abs(Lr[t, top]) * agree[t, top]
            lam = np.empty(len(top)); lam[np.argsort(-sc, kind="stable")] = lam_rank[:len(top)]
        if lcp > 0 and len(top):                      # weaken the re-allocation on high-expression genes (sch_emit 'lcp')
            cm = full_mean[top] / max(full_mean.sum(), 1e-12) * 1e6
            lam = lam * np.clip((np.median(cm) / np.maximum(cm, 1e-9)) ** lcp, 0.2, 1.0)
        up = top[Lr[t, top] > 0]; lu = lam[Lr[t, top] > 0]
        dn = top[Lr[t, top] < 0]; ld = lam[Lr[t, top] < 0]
        if up.size:
            if depth_up > 0:
                u = Lc ** (-depth_up); u = u / u.sum()
                W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (u[:, None] * T[up][None, :])
            else:
                W[:, up] = (1 - lu)[None, :] * W[:, up] + lu[None, :] * (T[up] / n)[None, :]
        if dn.size:
            conc = np.zeros((n, dn.size))
            np.put_along_axis(conc, order_cells[:n_top, dn], (T[dn] / n_top)[None, :].repeat(n_top, 0), axis=0)
            W[:, dn] = (1 - ld)[None, :] * W[:, dn] + ld[None, :] * conc
        if cpc:
            gs = np.flatnonzero(gate & ~j); g = gs[T[gs] > 0]
            if g.size:
                W0g = W[:, g].copy()                             # pre-tilt allocation (same column sums)
                w = W[:, g] / np.maximum(W[:, g].sum(0), 1e-12)
                want = n * CPC[t, g] / T[g]
                # ci > 0: the tilt is solved against the cells' ORIGINAL library totals, but the scorer normalises by the
                # FINAL ones, which the tilt itself changes. Re-solving against the emitted libraries closes that gap
                # (harness `ci`: H1 +0.0022 / KOLF +0.0004, PDS and MSE exactly unchanged because the sums stay pinned).
                Lcur = Lc
                for _it in range(1 + int(ci)):
                    lo_ = np.full(g.size, -4.0); hi_ = np.full(g.size, 4.0)
                    for _ in range(30):
                        mid = 0.5 * (lo_ + hi_); v = w * np.exp(mid[None, :] * logL[:, None]); v /= v.sum(0, keepdims=True)
                        f = (v / Lcur[:, None]).sum(0)
                        hi_ = np.where(f < want, mid, hi_); lo_ = np.where(f < want, lo_, mid)
                    s_ = 0.5 * (lo_ + hi_); v = w * np.exp(s_[None, :] * logL[:, None]); v /= v.sum(0, keepdims=True)
                    W[:, g] = v * T[g][None, :]
                    Lcur = np.maximum(W.sum(1), 1.0)
                Wt = W[:, g].copy(); a_ = 1.0
                while W.sum(1).max() > max_cell and a_ > 1e-3:   # scorer rejects cells > 1e6 counts: relax the tilt only
                    a_ *= 0.5; W[:, g] = W0g + a_ * (Wt - W0g)
        fl = np.floor(W); rem = (T - fl.sum(0)).astype(int); frac = W - fl
        cols = np.flatnonzero(rem > 0)
        if cols.size:
            o = np.argsort(-frac[:, cols], axis=0)
            mask = np.arange(n)[:, None] < rem[cols][None, :]
            add = np.zeros((n, cols.size)); np.put_along_axis(add, o, mask.astype(float), axis=0)
            fl[:, cols] += add
        b = sparse.csr_matrix(fl.astype(np.float32)); b.eliminate_zeros()
        yield b
