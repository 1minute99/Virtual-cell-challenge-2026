"""RQ2 Task 8 -- parameters of Proposition 2 estimated WITHOUT the gain curve (library + CLI).

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 8 (H2b), operational definitions 3-4.
Task 9 reuses the generic functions (participation_ratio, descriptor_pcs, halves_stats, aggregate,
bootstrap, c_star_*) on cell-level genetic data; only the loader differs.

QUANTITIES (definitions 3-4, unchanged):
  d        participation ratio (sum lam)^2 / sum lam^2 of the descriptor matrix (contexts x features),
           centred across contexts. lam = squared singular values of the centred matrix.
  Var(i)   noise-free interaction variance = <i^A, i^B> (cross-product of the two cell-level split-half
           residuals after removing mu, a_t, b_c, the additive fit done WITHIN each half).
  sig2_eps measurement-noise variance of the FULL-depth response = 1/2 * E||i^A - i^B||^2, rescaled from
           half to full depth (x 1/2: halving the cells doubles the sampling variance, A1 below).
  R2       (H2b only) adjusted R^2 of an OLS map from the top-k descriptor PCs (k = ceil(d) primary;
           floor(d), 2 ceil(d) sensitivity) to i across contexts, per drug, intercept included.
           Noise-free: fit on half A, score by the cross-product with half B,
               explained_t = <(H_t - J_t) i^A_t, i^B_t>,  total_t = <(I - J_t) i^A_t, i^B_t>
           (H_t = hat matrix of [1, PCs] on the contexts measuring drug t, J_t = mean operator; H_t is
           symmetric, so fit-on-A/score-on-B equals fit-on-B/score-on-A). Classical d.o.f. correction,
           pooled over drugs: R2_adj = 1 - [sum(total-explained)/sum((n_t-k-1)G)] / [sum total/sum((n_t-1)G)].
  C*_theory = d + 1 + d * sig2_res / (R2 * Var(i)),  sig2_res = (1 - R2) Var(i) + sig2_eps   (R2 = R2_adj, k = ceil d)
  C*_min    = d + 1 + d * sig2_eps / Var(i)
  UNITS: Var(i), sig2_eps are per-element averages, element = (context, drug, gene), in (log1p CP10K)^2,
  the response space of the model's target. Denominator = residual degrees of freedom of the additive
  fit (sum over elements of 1 - leverage), the same for both, so every ratio entering C* is
  denominator-free.

TAHOE INPUTS:
  descriptor  data/processed/drug_context_pretrain/contexts.npz['profile_adj'], Tahoe rows (48 x 5,955),
              identical to the Task 6 runner's context_descriptors() (checked in the tests; the loader is
              re-implemented here so that this module imports nothing from the Task 6 runner, which names the
              blinded curve file). Secondary d: the network's input z = (x - mean_train) @ B / s from
              dcp_data.context_space (top min(64, n-1) training PCs, RMS-scaled).
  responses   the model's target, rebuilt per half from RAW COUNTS: per Tahoe row (plate x sample x line,
              >= 25 cells, non-DMSO; 13,011 rows) log1p(1e4 * treated_sum / treated_native_library) minus
              log1p CP10K of the pooled DMSO_TF control of the same plate and line, on the 18,150-gene axis,
              rows averaged into (line, drug) with the dcp_prepare weight 1/(1/n_treated + 1/n_control)
              (computed from the half's own cell counts). The cells of EVERY group (each treated row and each
              DMSO well) are split into two equal random halves (floor/ceil), per split seed.
  cell source the 157 frozen Tahoe parquet shards listed in
              reports/context_signal_2026-09-15/tahoe/candidate_shard_plan.json under
              reports/context_signal_2026-09-15/tahoe/raw/ (override: --shard-root). If they are not ALL on
              disk with the frozen sizes, the split-half part STOPS (status 'blocked', exit code 3). No
              substitute (plate14 replicate rows, Poisson model) is used without approval.

ASSUMPTIONS (stated; checked where possible -- see knowledge/synthesis/rq2_params-assumptions.md):
  A1 sampling noise of a pseudobulk response scales as 1/(cells) (treated and control alike), so a
     half-depth difference has 4x the full-depth noise variance: sig2_full = (1/2) * [1/2 E||iA - iB||^2].
     Checked by simulation in tests (test_depth_rescaling_poisson).
  A2 noise is independent between the two halves (disjoint cells, disjoint control cells).
  A3 linear-Gaussian model of Proposition 2: i = f(x_c) + u, f linear in the top-k descriptor PCs.
  A4 the additive part is removed by an UNWEIGHTED two-way fit per gene; its leverage is the same for all
     genes (same observation mask), so the residual d.o.f. is exact for homoscedastic noise.
  A5 bootstrap over drugs and split seeds treats the additive fit and the descriptor as fixed (d is not
     resampled; its context-jackknife spread is reported separately).
  A6 per-element averages run over all 18,150 genes of the model's target axis.
BLINDING: this module reads no gain-vs-C result; tests/test_rq2_params.py greps it (and any module it imports
from output/) for the curve parquet, the dcp scaling outputs and other gain files.

CLI:  python rq2_params.py --dataset tahoe [--seeds 5] [--boot 2000] [--shard-root DIR]
  writes results/rq2_params_tahoe.json and results/rq2_h2b.json; exit 3 if the cell-level source is missing.
Environment: ~/miniconda3/envs/ct_transfer; run under nice -n 10, OMP/MKL threads 4.
"""
from __future__ import annotations

from pathlib import Path
import argparse, json, math, os, sys, time

import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC')
VCC = ROOT / 'vcc2026_data'
SCRIPTS = VCC / 'scripts'
PROC = VCC / 'data/processed/drug_context_pretrain'
TAHOE_META = VCC / 'reports/context_signal_2026-09-15/tahoe'
RESULTS = ROOT / 'output/data-analysis/results'

SPLIT_BASE = 20261108          # split-half RNG: default_rng([SPLIT_BASE, seed, group_id])
BOOT_SEED = 20261109
N_SPLIT_SEEDS = 5
N_BOOT = 2000
MIN_LINES = 3
DEPTH_SCALE = 0.5              # half depth -> full depth (A1)


# =========================================================================== d (definition 3)

def participation_ratio(X, center=True):
    """(sum lam)^2 / sum lam^2, lam = eigenvalues of the (centred) Gram/covariance of X (rows = contexts)."""
    X = np.asarray(X, np.float64)
    if center:
        X = X - X.mean(0)
    s = np.linalg.svd(X, compute_uv=False)
    lam = s ** 2
    return float(lam.sum() ** 2 / (lam ** 2).sum())


def descriptor_pcs(X):
    """PC scores (U*S) and eigenvalues of the context-centred descriptor (rows = contexts)."""
    Xc = np.asarray(X, np.float64) - np.asarray(X, np.float64).mean(0)
    U, s, _ = np.linalg.svd(Xc, full_matrices=False)
    r = int((s > s[0] * 1e-10).sum())
    return U[:, :r] * s[:r], s[:r] ** 2


def tahoe_descriptor():
    """(48 x 5,955) Tahoe descriptor, rows = sorted 'tahoe:<line>' keys (same as the Task 6 runner)."""
    z = np.load(PROC / 'contexts.npz', allow_pickle=True)
    names = np.array([str(x) for x in z['names']])
    m = np.array([str(x) for x in z['dataset']]) == 'tahoe'
    X = np.asarray(z['profile_adj'], np.float64)[m]
    nm = names[m]
    o = np.argsort(nm)
    return X[o], nm[o].tolist()


def _folds(names, k=6):        # = dcp_drug_lodo.folds (checked in tests)
    return {n: i % k for i, n in enumerate(sorted(names))}


def _subset(names, k):         # = dcp_drug_lodo.subset (seed-0 evenly spaced draw; checked in tests)
    names = sorted(names)
    if k >= len(names):
        return names
    step = len(names) / k
    return [names[min(len(names) - 1, int(i * step))] for i in range(k)]


def model_view_d(names, Cs=(6, 8, 12, 16, 24, 32, 40)):
    """Secondary d: participation ratio of the network's actual context input z (dcp_data.context_space)."""
    sys.path.insert(0, str(SCRIPTS))
    import dcp_data as DD
    out = {}
    cs = DD.context_space(list(names), dim=64)
    ti = [cs['index'][n] for n in names]
    Z = cs['Z'][ti].astype(np.float64)
    out['all48'] = dict(d=participation_ratio(Z), z_dim=int(cs['dim']), var_kept=cs['var_kept'],
                        note='PCA fitted on all 48 Tahoe contexts, k = min(64, 47); rank-preserving -> equals primary d')
    fm = _folds(names)
    per_fold = []
    for f in range(6):
        pool = sorted(n for n in names if fm[n] != f)
        csf = DD.context_space(pool, dim=64)
        Zall = csf['Z'][[csf['index'][n] for n in names]].astype(np.float64)
        Ztr = csf['Z'][[csf['index'][n] for n in pool]].astype(np.float64)
        per_fold.append(dict(fold=f, z_dim=int(csf['dim']), d_train40=participation_ratio(Ztr),
                             d_all48_projected=participation_ratio(Zall)))
    out['C40_folds'] = dict(per_fold=per_fold,
                            mean_d_train40=float(np.mean([p['d_train40'] for p in per_fold])),
                            mean_d_all48_projected=float(np.mean([p['d_all48_projected'] for p in per_fold])),
                            note='training pool of each of the 6 drug folds (= the C = 40 draw); held-out '
                                 'contexts projected on the training PCs (out-of-span part dropped)')
    vs_C = {}
    for C in Cs:
        ds = []
        for f in range(6):
            pool = sorted(n for n in names if fm[n] != f)
            tr = _subset(pool, C)
            csf = DD.context_space(tr, dim=64)
            ds.append(participation_ratio(csf['Z'][[csf['index'][n] for n in tr]].astype(np.float64)))
        vs_C[str(C)] = dict(mean=float(np.mean(ds)), min=float(np.min(ds)), max=float(np.max(ds)))
    out['train_subset_d_vs_C'] = dict(values=vs_C, note='DESCRIPTIVE: d of the C training descriptors '
                                      '(seed-0 evenly spaced draw, 6 folds); bounded by C - 1, so it is a '
                                      'finite-sample quantity, not the population d. Uses no response data.')
    return out


def d_jackknife(X):
    n = len(X)
    loo = np.array([participation_ratio(np.delete(X, i, 0)) for i in range(n)])
    se = float(np.sqrt((n - 1) / n * ((loo - loo.mean()) ** 2).sum()))
    return dict(loo_min=float(loo.min()), loo_max=float(loo.max()), jackknife_se=se)


# =========================================================================== split halves (generic)

def half_assignment(n_cells, seed, base=SPLIT_BASE):
    """Flat bool array over all cells (groups concatenated): True = half A. Exactly floor(n/2) per group in A.

    Group g's cells are identified by their ordinal 0..n_g-1 in a fixed reading order.
    """
    n_cells = np.asarray(n_cells, np.int64)
    out = np.zeros(int(n_cells.sum()), bool)
    off = 0
    for g, n in enumerate(n_cells):
        if n > 0:
            perm = np.random.default_rng([base, int(seed), g]).permutation(int(n))
            out[off:off + n] = perm < n // 2
        off += n
    return out


def log_cp10k(sums, libraries):
    sums = np.asarray(sums, np.float64)
    lib = np.asarray(libraries, np.float64)
    if np.any(lib <= 0):
        raise ValueError('every pseudobulk needs a positive library')
    return np.log1p(sums * 1e4 / lib[:, None])


def _additive_basis(ctx_idx, drug_idx):
    """Orthonormal basis of the column space of [1, drug dummies, context dummies] on the observed rows."""
    n = len(ctx_idx)
    uc, ci = np.unique(ctx_idx, return_inverse=True)
    ud, di = np.unique(drug_idx, return_inverse=True)
    X = np.zeros((n, 1 + len(ud) + len(uc)))
    X[:, 0] = 1
    X[np.arange(n), 1 + di] = 1
    X[np.arange(n), 1 + len(ud) + ci] = 1
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    r = int((s > s[0] * 1e-10).sum())
    return U[:, :r], r


def filter_min_lines(ctx_idx, drug_idx, min_lines=MIN_LINES):
    """Boolean keep-mask over observed rows: drugs measured in >= min_lines distinct contexts."""
    keep = np.ones(len(ctx_idx), bool)
    while True:
        cnt = {}
        for c, t in zip(ctx_idx[keep], drug_idx[keep]):
            cnt.setdefault(t, set()).add(c)
        bad = {t for t, s in cnt.items() if len(s) < min_lines}
        if not bad:
            return keep
        keep &= ~np.isin(drug_idx, list(bad))


def halves_stats(YA, YB, ctx_idx, drug_idx, pcs, ks, gene_chunk=2048):
    """Per-drug sufficient statistics from two split-half response matrices.

    YA, YB   (n_obs x G) responses of half A / half B, one row per observed (context, drug) element
    ctx_idx  context index of each row (indexes `pcs` rows); drug_idx: drug id of each row
    pcs      (n_contexts x r) descriptor PC scores (descriptor_pcs)
    ks       PC counts for the R^2 regressions
    Rows of drugs in < MIN_LINES contexts must be removed beforehand (filter_min_lines).
    Returns dict of per-drug arrays (aligned with `drugs`) plus meta.
    """
    ctx_idx, drug_idx = np.asarray(ctx_idx), np.asarray(drug_idx)
    Q, rank = _additive_basis(ctx_idx, drug_idx)
    lev = (Q ** 2).sum(1)
    drugs = np.unique(drug_idx)
    G = YA.shape[1]
    rows = [np.flatnonzero(drug_idx == t) for t in drugs]
    n_t = np.array([len(r) for r in rows])
    if n_t.min() < MIN_LINES:
        raise ValueError('drugs with < MIN_LINES contexts present; run filter_min_lines first')
    Hc = {}
    for k in ks:
        mats = []
        for r in rows:
            n = len(r)
            if n - k - 1 < 1:
                mats.append(None); continue
            Xd = np.column_stack([np.ones(n), pcs[ctx_idx[r], :k]])
            H = Xd @ np.linalg.pinv(Xd)
            mats.append(H - np.full((n, n), 1.0 / n))
        Hc[k] = mats
    T = len(drugs)
    cross = np.zeros(T); diff2 = np.zeros(T); tot = np.zeros(T)
    expl = {k: np.zeros(T) for k in ks}
    for lo in range(0, G, gene_chunk):
        hi = min(G, lo + gene_chunk)
        A = np.asarray(YA[:, lo:hi], np.float64); B = np.asarray(YB[:, lo:hi], np.float64)
        A = A - Q @ (Q.T @ A); B = B - Q @ (Q.T @ B)
        for j, r in enumerate(rows):
            a, b = A[r], B[r]
            cross[j] += float((a * b).sum())
            diff2[j] += float(((a - b) ** 2).sum())
            tot[j] += float(((a - a.mean(0)) * b).sum())
            for k in ks:
                if Hc[k][j] is not None:
                    expl[k][j] += float(((Hc[k][j] @ a) * b).sum())
        del A, B
    df = np.array([(1 - lev[r]).sum() for r in rows]) * G
    return dict(drugs=drugs, n_t=n_t, cross=cross, diff2=diff2, tot=tot, df=df, G=G,
                expl={int(k): v for k, v in expl.items()}, ks=[int(k) for k in ks],
                additive_rank=int(rank), n_obs=int(len(ctx_idx)), n_contexts=int(len(np.unique(ctx_idx))))


def c_star_theory(d, R2, var_i, sig2_eps):
    if not (R2 > 0 and var_i > 0):
        return math.inf
    sig2_res = (1 - R2) * var_i + sig2_eps
    return d + 1 + d * sig2_res / (R2 * var_i)


def c_star_min(d, var_i, sig2_eps):
    if not var_i > 0:
        return math.inf
    return d + 1 + d * sig2_eps / var_i


def _combine(S, wd, ws, d, k_primary, depth_scale):
    """Weighted pooled estimates; S = list (per seed) of halves_stats dicts with identical drug order.
    wd (B x T) drug multiplicities, ws (B x S) seed multiplicities. Returns dict of (B,) arrays."""
    stack = lambda key: np.stack([s[key] for s in S])                    # (S x T)
    W = np.einsum('bs,bt->bst', ws, wd)
    tot_sum = lambda M: np.einsum('bst,st->b', W, M)
    G = S[0]['G']
    n_t = S[0]['n_t']
    df = tot_sum(stack('df'))
    var_i = tot_sum(stack('cross')) / df
    sig2 = depth_scale * 0.5 * tot_sum(stack('diff2')) / df
    out = dict(var_i=var_i, sig2_eps=sig2, ratio=sig2 / var_i)
    tot = stack('tot')
    for k in S[0]['ks']:
        ok = (n_t - k - 1) >= 1
        ex = np.stack([s['expl'][k] for s in S])
        dfres = np.broadcast_to(np.where(ok, (n_t - k - 1) * G, 0.0), tot.shape)
        dftot = np.broadcast_to(np.where(ok, (n_t - 1) * G, 0.0), tot.shape)
        T_ = tot_sum(np.where(ok, tot, 0.0)); E_ = tot_sum(np.where(ok, ex, 0.0))
        out[f'R2_raw_k{k}'] = E_ / T_
        out[f'R2_adj_k{k}'] = 1 - ((T_ - E_) / tot_sum(dfres)) / (T_ / tot_sum(dftot))
        out[f'n_drugs_k{k}'] = np.full(len(df), int(ok.sum()))
    R2 = out[f'R2_adj_k{k_primary}']
    out['C_theory'] = np.array([c_star_theory(d, r, v, e) for r, v, e in zip(R2, var_i, sig2)])
    out['C_min'] = np.array([c_star_min(d, v, e) for v, e in zip(var_i, sig2)])
    for k in S[0]['ks']:
        out[f'C_theory_k{k}'] = np.array([c_star_theory(d, r, v, e) for r, v, e in
                                          zip(out[f'R2_adj_k{k}'], var_i, sig2)])
    return out


def aggregate(S, d, k_primary, depth_scale=DEPTH_SCALE, n_boot=N_BOOT, boot_seed=BOOT_SEED, alpha=0.05):
    """Point estimates (all drugs, all seeds) and percentile bootstrap CIs over drugs x split seeds."""
    for s in S[1:]:
        assert np.array_equal(s['drugs'], S[0]['drugs'])
    T, nS = len(S[0]['drugs']), len(S)
    point = _combine(S, np.ones((1, T)), np.ones((1, nS)), d, k_primary, depth_scale)
    rng = np.random.default_rng(boot_seed)
    wd = np.stack([np.bincount(rng.integers(0, T, T), minlength=T) for _ in range(n_boot)]).astype(float)
    ws = np.stack([np.bincount(rng.integers(0, nS, nS), minlength=nS) for _ in range(n_boot)]).astype(float)
    boot = _combine(S, wd, ws, d, k_primary, depth_scale)
    res = {}
    for key, v in point.items():
        b = boot[key]
        lo, hi = np.quantile(b, [alpha / 2, 1 - alpha / 2], method='inverted_cdf')   # inf-safe
        res[key] = dict(estimate=float(v[0]), ci95=[float(lo), float(hi)],
                        boot_frac_infinite=float(np.mean(~np.isfinite(b))))
    # per-seed point estimates (split-seed spread)
    res['per_seed'] = []
    for i in range(nS):
        ws1 = np.zeros((1, nS)); ws1[0, i] = 1
        p = _combine(S, np.ones((1, T)), ws1, d, k_primary, depth_scale)
        res['per_seed'].append({k: float(v[0]) for k, v in p.items()})
    res['bootstrap'] = dict(B=int(n_boot), seed=int(boot_seed), resampled=['drugs', 'split seeds'],
                            method='percentile (inverted_cdf)', fixed=['descriptor d', 'additive fit'])
    return res


# =========================================================================== SECONDARY sensitivities S1-S3
# Pre-specified 2026-11-05 (knowledge/_meta/log.md, "clarification", before any Task 7 result was seen by the main
# session); reported NEXT TO the primary H2b numbers, never replacing them.
#   S1  C* for the model's actual target b_c + i (b_c = line main effect, i = interaction):
#         signal  V_b + V_i;  noise  n_b = sig2_eps / nbar (b is an average over ~nbar drugs per line), n_i = sig2_eps
#         C*_S1      = d + 1 + d * [ (1-R2_b) V_b + n_b + (1-R2_i) V_i + n_i ] / [ R2_b V_b + R2_i V_i ]
#         C*_S1,min  = d + 1 + d * (n_b + n_i) / (V_b + V_i)
#       i.e. the excess term d*sig2_res/(R2 V) of Proposition 2 with residual and explained variances SUMMED over
#       the two independently learned components (the gain is the sum of the component gains), equivalently a
#       mixture of the component noise/signal ratios weighted by each component's (explained) variance.
#   S2  C*_b (line component alone): d + 1 + d [(1-R2_b) V_b + n_b] / (R2_b V_b);  C*_b,min = d + 1 + d n_b / V_b
#   S3  context-level bootstrap: resample the 48 lines with replacement (B = 1,000); in every draw recompute the
#       descriptor PCs, d and k = ceil(d) from the resampled rows, and all variance components / R^2 as weighted
#       sums of per-line sufficient statistics (additive fit and the half-specific b, i held fixed, as in A5).
# V_b noise-free = <(I-J) b^A, b^B> / ((L-1) G) from the two halves' line coefficients of the additive fit;
# R2_b = adjusted R^2 of b on [1, top-k PCs] across lines (fit on A, scored on B, as for i).

S3_BOOT = 1000
S3_SEED = 20261110


def halves_extra_stats(YA, YB, ctx_idx, drug_idx, gene_chunk=2048):
    """Per-drug cross matrices of the interaction residuals and per-line cross matrices of the line effects.

    M[t] = sum_g iA_g iB_g^T over the rows (contexts) of drug t; Dd[t] = diag sum_g (iA - iB)^2;
    Mb = sum_g bA_g bB_g^T (lines x lines), b = line coefficients of the additive fit, centred across lines;
    Db = diag sum_g (bA - bB)^2.  Everything the primary halves_stats returns is a function of these at unit
    weights (checked in tests), and they also support any line weighting (S3).
    """
    ctx_idx, drug_idx = np.asarray(ctx_idx), np.asarray(drug_idx)
    Q, rank = _additive_basis(ctx_idx, drug_idx)
    lev = (Q ** 2).sum(1)
    uc, ci = np.unique(ctx_idx, return_inverse=True)
    ud, di = np.unique(drug_idx, return_inverse=True)
    n = len(ctx_idx)
    X = np.zeros((n, 1 + len(ud) + len(uc)))
    X[:, 0] = 1
    X[np.arange(n), 1 + di] = 1
    X[np.arange(n), 1 + len(ud) + ci] = 1
    Pl = np.linalg.pinv(X)[1 + len(ud):]                  # (L x n): line coefficients (min-norm solution)
    Pl = Pl - Pl.mean(0)                                  # centred across lines -> identified contrasts
    rows = [np.flatnonzero(drug_idx == t) for t in ud]
    G = YA.shape[1]
    M = [np.zeros((len(r), len(r))) for r in rows]
    Dd = [np.zeros(len(r)) for r in rows]
    Mb = np.zeros((len(uc), len(uc))); Db = np.zeros(len(uc))
    for lo in range(0, G, gene_chunk):
        hi = min(G, lo + gene_chunk)
        A = np.asarray(YA[:, lo:hi], np.float64); B = np.asarray(YB[:, lo:hi], np.float64)
        bA, bB = Pl @ A, Pl @ B
        Mb += bA @ bB.T
        Db += ((bA - bB) ** 2).sum(1)
        A = A - Q @ (Q.T @ A); B = B - Q @ (Q.T @ B)
        for j, r in enumerate(rows):
            a, b = A[r], B[r]
            M[j] += a @ b.T
            Dd[j] += ((a - b) ** 2).sum(1)
        del A, B
    return dict(drugs=ud, rows_ctx=[ctx_idx[r] for r in rows], one_minus_lev=[1 - lev[r] for r in rows],
                M=M, Dd=Dd, lines=uc, Mb=Mb, Db=Db, G=int(G), n_seeds=1,
                drugs_per_line=np.bincount(ci, minlength=len(uc)))


def sum_extra_stats(E):
    """Sum per-seed extra stats (identical design) into one dict; n_seeds counts the summed seeds."""
    out = dict(E[0])
    for e in E[1:]:
        assert np.array_equal(e['drugs'], out['drugs']) and np.array_equal(e['lines'], out['lines'])
    out['M'] = [sum(e['M'][j] for e in E) for j in range(len(out['M']))]
    out['Dd'] = [sum(e['Dd'][j] for e in E) for j in range(len(out['Dd']))]
    out['Mb'] = sum(e['Mb'] for e in E); out['Db'] = sum(e['Db'] for e in E)
    out['n_seeds'] = len(E)
    return out


def weighted_descriptor(X, w):
    """PC scores of every context row and d for a weighted (= resampled-with-duplicates) descriptor."""
    X = np.asarray(X, np.float64); w = np.asarray(w, np.float64)
    xbar = w @ X / w.sum()
    Xc = X - xbar
    _, s, Vt = np.linalg.svd(np.sqrt(w)[:, None] * Xc, full_matrices=False)
    r = int((s > s[0] * 1e-10).sum())
    lam = s[:r] ** 2
    return Xc @ Vt[:r].T, float(lam.sum() ** 2 / (lam ** 2).sum())


def _wterms(Mm, w, Z):
    """trace(W (I-J_w) M) and trace(W (H_w - J_w) M) for a weighted OLS on design [1, Z]."""
    n = len(w)
    sw = w.sum()
    J = np.outer(np.ones(n), w / sw)
    Xd = np.column_stack([np.ones(n), Z])
    H = Xd @ np.linalg.pinv((Xd.T * w) @ Xd) @ (Xd.T * w)
    I = np.eye(n)
    tot = float((w[:, None] * (I - J) * Mm.T).sum())
    expl = float((w[:, None] * (H - J) * Mm.T).sum())
    return tot, expl


def _adj(tot, expl, dfres, dftot):
    return 1 - ((tot - expl) / dfres) / (tot / dftot)


def sensitivity_estimates(E, X, w=None, k=None, depth_scale=DEPTH_SCALE, min_lines=MIN_LINES):
    """Primary-equivalent and S1/S2 quantities from summed extra stats E, for context weights w (len = rows of X;
    None = all ones = the point estimate). k = None -> ceil(d) of the weighted descriptor."""
    nctx = len(X)
    w = np.ones(nctx) if w is None else np.asarray(w, float)
    pcs, d = weighted_descriptor(X, w)
    k = math.ceil(d) if k is None else int(k)
    G, S = E['G'], E['n_seeds']
    cross = diff2 = df = 0.0
    T_ = E_ = dfr = dft = 0.0
    n_drugs = 0
    for j in range(len(E['drugs'])):
        rc = E['rows_ctx'][j]
        wr = w[rc]
        m = wr > 0
        if len(np.unique(rc[m])) < min_lines:
            continue
        n_drugs += 1
        Mm = E['M'][j][np.ix_(m, m)]; wr = wr[m]
        cross += float((wr * np.diag(Mm)).sum())
        diff2 += float((wr * E['Dd'][j][m]).sum())
        df += float((wr * E['one_minus_lev'][j][m]).sum()) * G * S
        neff = wr.sum()
        if neff - k - 1 >= 1:
            tot, ex = _wterms(Mm, wr, pcs[rc[m], :k])
            T_ += tot; E_ += ex
            dfr += (neff - k - 1) * G * S; dft += (neff - 1) * G * S
    var_i = cross / df
    sig2 = depth_scale * 0.5 * diff2 / df
    R2_i = _adj(T_, E_, dfr, dft) if dft > 0 else float('nan')
    # line component
    lw = w[E['lines']]
    mL = lw > 0
    Mb = E['Mb'][np.ix_(mL, mL)]; lw_ = lw[mL]
    nL = lw_.sum()
    totb, exb = _wterms(Mb, lw_, pcs[E['lines'][mL], :k])
    var_b = totb / ((nL - 1) * G * S)
    R2_b = _adj(totb, exb, (nL - k - 1) * G * S, (nL - 1) * G * S) if nL - k - 1 >= 1 else float('nan')
    nbar = float((lw_ * E['drugs_per_line'][mL]).sum() / nL)
    n_b = sig2 / nbar
    n_b_emp = depth_scale * 0.5 * float((lw_ * E['Db'][mL]).sum()) / ((nL - 1) * G * S)
    out = dict(d=d, k=k, n_drugs=n_drugs, n_lines_eff=float(nL), var_i=var_i, sig2_eps=sig2, R2_adj_i=R2_i,
               var_b=var_b, R2_adj_b=R2_b, nbar_drugs_per_line=nbar, noise_b_model=n_b,
               noise_b_empirical=n_b_emp)
    out['C_theory'] = c_star_theory(d, R2_i, var_i, sig2)
    out['C_min'] = c_star_min(d, var_i, sig2)
    den = (R2_b * var_b if R2_b > 0 else 0.0) + (R2_i * var_i if R2_i > 0 else 0.0)
    res_ = (1 - R2_b) * var_b + n_b + (1 - R2_i) * var_i + sig2
    out['S1_C_theory'] = d + 1 + d * res_ / den if den > 0 else math.inf
    out['S1_C_min'] = d + 1 + d * (n_b + sig2) / (var_b + var_i)
    out['S2_C_b_theory'] = c_star_theory(d, R2_b, var_b, n_b)
    out['S2_C_b_min'] = c_star_min(d, var_b, n_b)
    # FLAGGED VARIANT (not pre-specified): b's noise from its own split halves. All drugs of a plate x line share
    # one DMSO control pseudobulk, so control noise does not average out over drugs and lands in b.
    den_e = den
    res_e = (1 - R2_b) * var_b + n_b_emp + (1 - R2_i) * var_i + sig2
    out['S1emp_C_theory'] = d + 1 + d * res_e / den_e if den_e > 0 else math.inf
    out['S1emp_C_min'] = d + 1 + d * (n_b_emp + sig2) / (var_b + var_i)
    out['S2emp_C_b_theory'] = c_star_theory(d, R2_b, var_b, n_b_emp)
    out['S2emp_C_b_min'] = c_star_min(d, var_b, n_b_emp)
    return out


def line_bootstrap(E, X, n_boot=S3_BOOT, seed=S3_SEED, alpha=0.05, log=None):
    """S3: resample contexts (lines) with replacement; d, PCs, k and every component recomputed per draw."""
    rng = np.random.default_rng(seed)
    n = len(X)
    draws = []
    for b in range(n_boot):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        draws.append(sensitivity_estimates(E, X, w))
        if log and (b + 1) % 100 == 0:
            log(json.dumps(dict(event='s3_boot', done=b + 1, total=n_boot)))
    keys = [k for k in draws[0] if isinstance(draws[0][k], float)]
    res = {}
    for key in keys + ['k']:
        v = np.array([dr[key] for dr in draws], float)
        fin = np.isfinite(v)
        lo, hi = np.quantile(np.where(fin, v, np.inf), [alpha / 2, 1 - alpha / 2], method='inverted_cdf')
        res[key] = dict(ci95=[float(lo), float(hi)], boot_frac_infinite=float(np.mean(~fin)),
                        boot_median=float(np.quantile(np.where(fin, v, np.inf), 0.5, method='inverted_cdf')))
    res['bootstrap'] = dict(B=int(n_boot), seed=int(seed), resampled=['contexts (lines)'],
                            recomputed=['descriptor centring + PCs', 'd', 'k = ceil(d)', 'Var(i)', 'sig2_eps',
                                        'R2 (i and b)', 'Var(b)', 'nbar'],
                            fixed=['additive fit and half-specific b, i (as A5)', 'split seeds pooled'],
                            method='percentile (inverted_cdf)')
    return res


def line_jackknife(E, X, z=1.959964):
    """S3b (FLAGGED alternative to the S3 line bootstrap): delete-1 jackknife over contexts (lines), with d, PCs,
    k = ceil(d) and every component recomputed without the deleted line. Resampling WITH replacement keeps only
    ~63 % distinct lines, which caps the rank of the descriptor and duplicates rows in the R^2 fits; on Tahoe that
    pulls d down and R^2_adj up so far that the S3 percentile CIs exclude the point estimates. Deleting one line
    keeps the statistic in its own regime. CI = point +- z * SE_jack (normal approximation)."""
    n = len(X)
    pt = sensitivity_estimates(E, X)
    loo = []
    for c in range(n):
        w = np.ones(n); w[c] = 0.0
        loo.append(sensitivity_estimates(E, X, w))
    res = {}
    for key, v0 in pt.items():
        if not isinstance(v0, (float, int)) or key == 'n_drugs':
            continue
        v = np.array([dr[key] for dr in loo], float)
        if not np.all(np.isfinite(v)):
            res[key] = dict(se=None, ci95=None, loo_frac_infinite=float(np.mean(~np.isfinite(v))))
            continue
        se = float(np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum()))
        res[key] = dict(se=se, ci95=[float(v0 - z * se), float(v0 + z * se)],
                        bias_jack=float((n - 1) * (v.mean() - v0)), loo_min=float(v.min()), loo_max=float(v.max()))
    res['method'] = dict(kind='delete-1 jackknife over the 48 lines', recomputed='d, PCs, k, all components',
                         fixed='additive fit, half-specific b and i (A5)', ci='normal, point +- 1.96 SE')
    return res


def save_extra_stats(E, path):
    """Cache of the summed sufficient statistics (small) so S1-S3 can be recomputed without the shards."""
    obj = dict(drugs=E['drugs'], lines=E['lines'], Mb=E['Mb'], Db=E['Db'], G=E['G'], n_seeds=E['n_seeds'],
               drugs_per_line=E['drugs_per_line'], n_drug=len(E['drugs']))
    for j in range(len(E['drugs'])):
        obj[f'M_{j}'] = E['M'][j]; obj[f'Dd_{j}'] = E['Dd'][j]
        obj[f'rc_{j}'] = E['rows_ctx'][j]; obj[f'oml_{j}'] = E['one_minus_lev'][j]
    np.savez_compressed(path, **obj)


def load_extra_stats(path):
    z = np.load(path)
    T = int(z['n_drug'])
    return dict(drugs=z['drugs'], lines=z['lines'], Mb=z['Mb'], Db=z['Db'], G=int(z['G']), n_seeds=int(z['n_seeds']),
                drugs_per_line=z['drugs_per_line'], M=[z[f'M_{j}'] for j in range(T)],
                Dd=[z[f'Dd_{j}'] for j in range(T)], rows_ctx=[z[f'rc_{j}'] for j in range(T)],
                one_minus_lev=[z[f'oml_{j}'] for j in range(T)])


# =========================================================================== Tahoe cell-level reader

def _tahoe_helpers():
    sys.path.insert(0, str(SCRIPTS))
    import simple_residual_prepare as SRP
    import context_signal_tahoe_effects as TE
    return SRP, TE


def tahoe_row_table():
    """Rebuild the 13,011-row Tahoe table structure from METADATA ONLY (no counts) + model indices."""
    SRP, _ = _tahoe_helpers()
    groups = SRP._tahoe_groups()
    sel = SRP.select_tahoe_rows(groups['sample'], groups['cell_line'], groups['plate'], groups['drug'],
                                groups['cells'].to_numpy(), SRP.MIN_CELLS)
    take = sel['indices']
    names = SRP._tahoe_line_names()
    line_names = [names.get(v, v) for v in groups.cell_line.iloc[take]]
    bg = np.array([SRP.biological_family(v) for v in line_names])
    drug = groups.drug.iloc[take].astype(str).to_numpy()
    ctx_names, pert_names = sorted(set(bg)), sorted(set(drug))
    ctrl = groups.drug.eq('DMSO_TF').to_numpy()
    ckeys = sorted(set(zip(groups.plate[ctrl], groups.cell_line[ctrl])))
    cidx = {k: i for i, k in enumerate(ckeys)}
    ctrl_groups = [np.flatnonzero(ctrl & groups.plate.eq(p).to_numpy() & groups.cell_line.eq(l).to_numpy())
                   for p, l in ckeys]
    row_ctrl = np.array([cidx[(groups.plate.iloc[i], groups.cell_line.iloc[i])] for i in take])
    return dict(groups=groups, take=take, is_repeat=sel['is_repeat'],
                ctx=np.array([ctx_names.index(b) for b in bg]), ctx_names=ctx_names,
                pert=np.array([pert_names.index(x) for x in drug]), pert_names=pert_names,
                row_ctrl=row_ctrl, ctrl_groups=ctrl_groups, ctrl_keys=ckeys,
                cells=groups['cells'].to_numpy(np.int64))


def check_row_table_vs_model(rt):
    """Compare the metadata-rebuilt table with the model's tahoe.npz (indices, repeat flag, weights)."""
    z = np.load(PROC / 'tahoe.npz', allow_pickle=True)
    nt = rt['cells'][rt['take']].astype(float)
    nc = np.array([rt['cells'][g].sum() for g in rt['ctrl_groups']], float)[rt['row_ctrl']]
    w = 1.0 / (1.0 / np.maximum(nt, 1) + 1.0 / np.maximum(nc, 1))
    return dict(rows=int(len(rt['take'])), rows_model=int(len(z['ctx'])),
                ctx_names_equal=[str(x) for x in z['ctx_names']] == rt['ctx_names'],
                pert_names_equal=[str(x) for x in z['pert_names']] == rt['pert_names'],
                ctx_equal=bool(np.array_equal(z['ctx'], rt['ctx'])),
                pert_equal=bool(np.array_equal(z['pert'], rt['pert'])),
                is_repeat_equal=bool(np.array_equal(z['is_repeat'], rt['is_repeat'])),
                weight_max_rel_diff=float(np.max(np.abs(z['weight'] - w) / w)),
                cells_total=int(rt['cells'].sum()))


def tahoe_source_status(shard_root=None):
    plan = json.loads((TAHOE_META / 'candidate_shard_plan.json').read_text())
    shards = plan['shards'] if isinstance(plan, dict) else plan
    root = Path(shard_root) if shard_root else TAHOE_META / 'raw'
    present = [s for s in shards if (root / s['path']).exists()
               and (root / s['path']).stat().st_size == int(s['bytes'])]
    return dict(shard_root=str(root), shards_expected=len(shards), shards_present=len(present),
                bytes_expected=int(sum(int(s['bytes']) for s in shards)),
                source_revision='2dc57900b7981cfcf5e211527169a0b006546a95',
                complete=len(present) == len(shards), shards=shards)


class HalfAccumulator:
    """Streams raw cell counts into per-group FULL sums and per-seed HALF-A sums (B = full - A)."""

    def __init__(self, n_cells_per_group, gene_lookup, n_genes, seeds, base=SPLIT_BASE):
        self.n = np.asarray(n_cells_per_group, np.int64)
        self.off = np.concatenate([[0], np.cumsum(self.n)])
        self.lookup = np.asarray(gene_lookup, np.int64)
        self.G = int(n_genes)
        self.seeds = list(seeds)
        self.inA = [half_assignment(self.n, s, base) for s in self.seeds]
        ng = len(self.n)
        self.full = np.zeros((ng, self.G), np.float32)
        self.A = [np.zeros((ng, self.G), np.float32) for _ in self.seeds]
        self.lib_full = np.zeros(ng); self.lib_A = [np.zeros(ng) for _ in self.seeds]
        self.seen = np.zeros(ng, np.int64)

    def add(self, codes, offsets, tokens, counts):
        codes = np.asarray(codes, np.int64)
        offsets = np.asarray(offsets, np.int64)
        lo, hi = int(offsets[0]), int(offsets[-1])
        tokens = np.asarray(tokens, np.int64)[lo:hi]
        counts = np.asarray(counts, np.float64)[lo:hi]
        offsets = offsets - lo
        if np.any(counts < 0) or np.any(np.abs(counts - np.rint(counts)) > 1e-5):
            raise ValueError('non-count expression values')
        kept = codes >= 0
        if not kept.any():
            return
        # ordinal of each kept cell within its group, in reading order
        kc = codes[kept]
        order = np.argsort(kc, kind='stable')
        sc = kc[order]
        start = np.r_[0, np.flatnonzero(np.diff(sc)) + 1]
        run = np.repeat(start, np.diff(np.r_[start, len(sc)]))
        ordk = np.empty(len(kc), np.int64)
        ordk[order] = self.seen[sc] + (np.arange(len(sc)) - run)
        self.seen += np.bincount(kc, minlength=len(self.n))
        if np.any(ordk >= self.n[kc]):
            raise ValueError('more cells than the frozen group cell count')
        flat = self.off[kc] + ordk
        lengths = np.diff(offsets)
        cell_of_tok = np.repeat(np.arange(len(codes)), lengths)
        lib = np.bincount(cell_of_tok, weights=counts, minlength=len(codes))
        uniq, loc = np.unique(kc, return_inverse=True)
        local = np.full(len(codes), -1, np.int64); local[kept] = loc
        tg = local[cell_of_tok]
        gi = self.lookup[tokens]
        keep = (tg >= 0) & (gi >= 0)
        idx = tg[keep] * self.G + gi[keep]
        ctok = cell_of_tok[keep]
        nb = len(uniq) * self.G
        self.full[uniq] += np.bincount(idx, weights=counts[keep], minlength=nb).reshape(len(uniq), self.G)
        self.lib_full += np.bincount(kc, weights=lib[kept], minlength=len(self.n))
        cellA_full = np.zeros(len(codes), bool)
        for s, m in enumerate(self.inA):
            cellA_full[:] = False
            cellA_full[np.flatnonzero(kept)] = m[flat]
            wA = cellA_full[ctok]
            self.A[s][uniq] += np.bincount(idx[wA], weights=counts[keep][wA],
                                           minlength=nb).reshape(len(uniq), self.G)
            self.lib_A[s] += np.bincount(kc[m[flat]], weights=lib[kept][m[flat]], minlength=len(self.n))

    def finish_check(self):
        if not np.array_equal(self.seen, self.n):
            raise ValueError('scanned cell counts differ from the frozen group table')
        if self.full.max() >= 2 ** 24:
            raise ValueError('float32 sums no longer exact')


def read_tahoe_shards(acc, key_to_group, shards, shard_root, TE, log=print):
    import pyarrow.parquet as pq
    t0 = time.monotonic()
    for i, sh in enumerate(shards, 1):
        pf = pq.ParquetFile(Path(shard_root) / sh['path'])
        seen = 0
        for batch in pf.iter_batches(batch_size=2048, columns=['genes', 'expressions', 'sample',
                                                               'cell_line_id', 'plate'], use_threads=False):
            if not all(p == sh['expected_plate'] for p in batch.column('plate').to_pylist()):
                raise ValueError('shard plate differs from the frozen plan')
            codes = np.array([key_to_group.get((s, c), -1) for s, c in
                              zip(batch.column('sample').to_pylist(), batch.column('cell_line_id').to_pylist())])
            ga, xa = batch.column('genes'), batch.column('expressions')
            offs = ga.offsets.to_numpy()
            if not np.array_equal(offs, xa.offsets.to_numpy()):
                raise ValueError('gene/count offsets disagree')
            offs, tok, val = TE.remove_cls(offs, ga.values.to_numpy(), xa.values.to_numpy())
            acc.add(codes, offs, tok, val)
            seen += len(batch)
        if seen != int(sh['expected_rows']):
            raise ValueError('shard row count differs from the frozen plan')
        if i == 1 or i % 10 == 0 or i == len(shards):
            log(json.dumps(dict(event='shard', done=i, total=len(shards),
                                secs=round(time.monotonic() - t0, 1))))
    acc.finish_check()


def tahoe_half_responses(acc, rt, s, with_full=False, row_chunk=1024):
    """Observed (context, drug) response matrices for split-seed index s.

    Returns ctx, drug, YA, YB[, YF]: one row per (context, drug) with >= 1 Tahoe row, float32 (n_obs x G);
    per Tahoe row: log1p CP10K(treated half) - log1p CP10K(pooled DMSO half of the same plate x line);
    rows averaged with weight 1/(1/n_treated + 1/n_control) of that half (dcp_prepare rule).
    """
    take, rc, cg = rt['take'], rt['row_ctrl'], rt['ctrl_groups']
    nA = np.array([acc.inA[s][acc.off[g]:acc.off[g + 1]].sum() for g in range(len(acc.n))])
    halves = {'A': (lambda idx: acc.A[s][idx].astype(np.float64), acc.lib_A[s], nA),
              'B': (lambda idx: acc.full[idx].astype(np.float64) - acc.A[s][idx], acc.lib_full - acc.lib_A[s],
                    acc.n - nA)}
    if with_full:
        halves['F'] = (lambda idx: acc.full[idx].astype(np.float64), acc.lib_full, acc.n)
    nP = len(rt['pert_names'])
    keys = rt['ctx'] * nP + rt['pert']
    uk, inv = np.unique(keys, return_inverse=True)
    out = []
    for h, (get, lib, ncell) in halves.items():
        # control pseudobulks only for the pools some row uses: an unused pool may be empty in a half
        # (plate14 x CVCL_1571: 2 wells x 1 cell -> floor(1/2) = 0 cells in half A)
        clog = np.full((len(cg), acc.G), np.nan)
        cn = np.zeros(len(cg))
        used = np.unique(rc)
        clog[used] = log_cp10k(np.stack([get(cg[u]).sum(0) for u in used]),
                               np.array([lib[cg[u]].sum() for u in used]))
        cn[used] = [ncell[cg[u]].sum() for u in used]
        w = 1.0 / (1.0 / np.maximum(ncell[take], 1) + 1.0 / np.maximum(cn[rc], 1))
        Y = np.zeros((len(uk), acc.G), np.float64)
        for lo in range(0, len(take), row_chunk):
            sl = slice(lo, lo + row_chunk)
            resp = log_cp10k(get(take[sl]), lib[take[sl]]) - clog[rc[sl]]
            np.add.at(Y, inv[sl], resp * w[sl, None])
        Y /= np.bincount(inv, weights=w)[:, None]
        out.append(Y.astype(np.float32))
        del Y
    return (uk // nP, uk % nP, *out)


# =========================================================================== CLI

def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    return x


def write_json(path, obj, keep_foreign_keys=False):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if keep_foreign_keys and path.exists():
        # keys added by others (e.g. the main session's H2b verdict) are carried over OPAQUELY: never parsed or
        # used here, only preserved so that a rerun does not delete them
        old = json.loads(path.read_text())
        obj = {**{k: v for k, v in old.items() if k not in obj}, **obj}
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(_jsonable(obj), indent=2) + '\n')
    os.replace(tmp, path)


def run_tahoe(a):
    t0 = time.monotonic()
    X, names = tahoe_descriptor()
    pcs, lam = descriptor_pcs(X)
    d = participation_ratio(X)
    ks = sorted({math.floor(d), math.ceil(d), 2 * math.ceil(d)})
    k_primary = math.ceil(d)
    res = dict(task='Task 8 (RQ2, H2b) -- Proposition-2 parameters on Tahoe, estimated without the gain curve',
               plan='input/ideas/zero-shot-limits-context-complexity-plan.md', created=time.strftime('%Y-%m-%d'),
               script='output/data-analysis/rq2_params.py')
    res['d'] = dict(primary=dict(value=d, matrix='contexts.npz profile_adj, 48 Tahoe contexts x 5,955 genes, '
                                 'centred across contexts', n_contexts=len(names), n_features=int(X.shape[1]),
                                 top_eig_frac=(lam[:10] / lam.sum()).tolist(), context_jackknife=d_jackknife(X)),
                    k_used=dict(primary=k_primary, sensitivity=[k for k in ks if k != k_primary]))
    print(json.dumps(dict(event='d', d=d, ks=ks)), flush=True)
    res['d']['secondary_model_input'] = model_view_d(names)
    print(json.dumps(dict(event='d_secondary', all48=res['d']['secondary_model_input']['all48']['d'])), flush=True)
    rt = tahoe_row_table()
    res['row_table_check'] = check_row_table_vs_model(rt)
    st = tahoe_source_status(a.shard_root)
    res['cell_level_source'] = {k: v for k, v in st.items() if k != 'shards'}
    res['units'] = ('Var(i), sig2_eps: per (line, drug, gene) element, (log1p CP10K)^2, denominator = residual '
                    'd.o.f. of the additive fit; R2 dimensionless; C* in contexts')
    h2b = dict(hypothesis='H2b: C*_theory lies within a factor of 2 of C*_emp (rejected if C*_theory/C*_emp '
                          'outside [0.5, 2]; if C*_emp censored at > 40, rejected only if C*_theory < 20)',
               verdict=None,
               verdict_note='The H2b verdict requires C*_emp from Task 7 (its H2a result file); it is NOT computed '
                            'here. The comparison is done by the main session.')
    if not st['complete']:
        msg = (f"BLOCKED: cell-level split halves impossible -- {st['shards_present']} of "
               f"{st['shards_expected']} frozen Tahoe shards on disk under {st['shard_root']} "
               f"({st['bytes_expected'] / 2**30:.1f} GiB expected); the processed table.npz is also absent. "
               "Per the task spec no substitute (plate14 replicate rows) is used without approval.")
        res['status'] = 'blocked'
        res['split_half'] = dict(status='blocked', reason=msg)
        res['C_star_min_lower_bound_from_d_only'] = dict(
            value=d + 1, note='C*_min >= d + 1 for any sig2_eps/Var(i) >= 0 (algebra, no noise estimate)')
        h2b.update(status='blocked', C_star_theory=None, reason=msg, d=d)
        write_json(a.out, res)
        write_json(a.h2b_out, h2b, keep_foreign_keys=True)
        print(json.dumps(dict(event='blocked', msg=msg)), flush=True)
        return 3
    # ---------------- full run (only reached when every shard is on disk)
    SRP, TE = _tahoe_helpers()
    genes, lookup = SRP._tahoe_gene_axis()
    g = rt['groups']
    key_to_group = {(r.sample, r.cell_line): int(r.group_id) for r in g.itertuples()}
    seeds = list(range(a.seeds))
    acc = HalfAccumulator(rt['cells'], lookup, len(genes), seeds)
    read_tahoe_shards(acc, key_to_group, st['shards'], st['shard_root'], TE)
    S = []
    EX = []
    for si in range(len(seeds)):
        if si == 0:   # full-depth reconstruction check against the model's own target rows
            ctx, drug, YA, YB, YF = tahoe_half_responses(acc, rt, si, with_full=True)
            z = np.load(PROC / 'tahoe.npz', allow_pickle=True)
            res['full_depth_check'] = dict(note='YF vs the model-target (ctx, drug) means rebuilt from tahoe.npz',
                                           **_full_check(z, ctx, drug, YF))
            del YF, z
        else:
            ctx, drug, YA, YB = tahoe_half_responses(acc, rt, si)
        keep = filter_min_lines(ctx, drug)
        S.append(halves_stats(YA[keep], YB[keep], ctx[keep], drug[keep], pcs, ks))
        if getattr(a, 'sensitivity', False):
            EX.append(halves_extra_stats(YA[keep], YB[keep], ctx[keep], drug[keep]))
        print(json.dumps(dict(event='seed', seed=si, n_obs=int(keep.sum()), drugs=len(S[-1]['drugs']),
                              secs=round(time.monotonic() - t0, 1))), flush=True)
        del YA, YB
    agg = aggregate(S, d, k_primary, n_boot=a.boot)
    s0 = S[0]
    res['status'] = 'complete'
    res['split_half'] = dict(seeds=seeds, split_base=SPLIT_BASE, depth_scale=DEPTH_SCALE,
                             n_drugs=len(s0['drugs']), n_lines=s0['n_contexts'], n_genes=s0['G'],
                             n_obs=s0['n_obs'], additive_rank=s0['additive_rank'],
                             estimates=agg)
    ct = agg['C_theory']
    h2b.update(status='complete', d=d, k_primary=k_primary,
               C_star_theory=dict(estimate=ct['estimate'], ci95=ct['ci95']),
               C_star_theory_sensitivity={f'k{k}': agg[f'C_theory_k{k}'] for k in ks},
               C_star_min=agg['C_min'])
    if getattr(a, 'sensitivity', False):
        E = sum_extra_stats(EX)
        pt = sensitivity_estimates(E, X)
        chk = dict(C_theory=pt['C_theory'] - ct['estimate'], C_min=pt['C_min'] - agg['C_min']['estimate'],
                   var_i=pt['var_i'] - agg['var_i']['estimate'], sig2_eps=pt['sig2_eps'] - agg['sig2_eps']['estimate'])
        print(json.dumps(dict(event='s3_start', point=pt, primary_route_diff=chk)), flush=True)
        save_extra_stats(E, RESULTS / 'logs' / 'rq2_params_tahoe_extra_stats.npz')
        boot = line_bootstrap(E, X, n_boot=a.line_boot, log=lambda m: print(m, flush=True))
        jk = line_jackknife(E, X)
        boot_valid = {key: (boot[key]['ci95'][0] <= pt[key] <= boot[key]['ci95'][1])
                      for key in ('d', 'R2_adj_i', 'C_theory', 'C_min', 'S1_C_theory', 'S2_C_b_theory')}
        res['sensitivity'] = dict(
            label='SECONDARY, pre-specified 2026-11-05 (log: clarification); does not change the H2b rule',
            definitions=dict(
                S1='C* for the model target b + i: C*_S1 = d+1+d[(1-R2_b)V_b + n_b + (1-R2_i)V_i + n_i]/[R2_b V_b + R2_i V_i]; '
                   'C*_S1,min = d+1+d(n_b+n_i)/(V_b+V_i); n_i = sig2_eps, n_b = sig2_eps/nbar',
                S2='C*_b = d+1+d[(1-R2_b)V_b + n_b]/(R2_b V_b); C*_b,min = d+1+d n_b/V_b',
                S3='line (context) bootstrap, d / PCs / k recomputed per draw; CIs for primary and S1/S2',
                emp_variant='FLAGGED, not pre-specified: S1emp/S2emp replace n_b = sig2_eps/nbar by the empirical '
                            'split-half noise of b (noise_b_empirical = 1/2 * 1/2 ||(I-J)(bA-bB)||^2/((L-1)G)); '
                            'shared DMSO control noise per plate x line does not average over drugs'),
            point=pt, point_route_minus_primary=chk, S3_line_bootstrap=boot,
            S3_validity=dict(point_inside_boot_ci=boot_valid,
                             note='FLAG: resampling lines with replacement leaves ~63 % distinct lines; d (rank-capped) '
                                  'is biased down and duplicated rows inflate R2_adj, so the S3 percentile CIs are not '
                                  'valid where the point estimate falls outside them. Use S3b.'),
            S3b_line_jackknife=jk)
        h2b['sensitivity_secondary'] = dict(
            S1_C_theory=dict(estimate=pt['S1_C_theory'], ci95_line_boot=boot['S1_C_theory']['ci95']),
            S1_C_min=dict(estimate=pt['S1_C_min'], ci95_line_boot=boot['S1_C_min']['ci95']),
            S2_C_b_theory=dict(estimate=pt['S2_C_b_theory'], ci95_line_boot=boot['S2_C_b_theory']['ci95']),
            S3_C_theory_ci95_line_boot=boot['C_theory']['ci95'], S3_C_min_ci95_line_boot=boot['C_min']['ci95'],
            S3b_C_theory_ci95_line_jackknife=jk['C_theory']['ci95'], S3b_C_min_ci95_line_jackknife=jk['C_min']['ci95'],
            S3b_S1_C_theory_ci95_line_jackknife=jk['S1_C_theory']['ci95'],
            S3_validity_point_inside_boot_ci=boot_valid,
            S1emp_C_theory=dict(estimate=pt['S1emp_C_theory'], ci95_line_boot=boot['S1emp_C_theory']['ci95']),
            S2emp_C_b_theory=dict(estimate=pt['S2emp_C_b_theory'], ci95_line_boot=boot['S2emp_C_b_theory']['ci95']),
            note='secondary; the H2b verdict uses the primary C*_theory. *emp = flagged variant with the '
                 'empirical split-half noise of b instead of the pre-specified sig2_eps/nbar')
    write_json(a.out, res)
    write_json(a.h2b_out, h2b, keep_foreign_keys=True)
    print(json.dumps(dict(event='done', secs=round(time.monotonic() - t0, 1))), flush=True)
    return 0


def _full_check(z, ctx, drug, YF):
    resp, zc, zp, w = z['resp'], z['ctx'], z['pert'], z['weight'].astype(np.float64)
    nP = len(z['pert_names'])
    pick = np.arange(0, len(ctx), max(1, len(ctx) // 200))
    errs = []
    for i in pick:
        m = (zc == ctx[i]) & (zp == drug[i])
        ref = (resp[m] * w[m][:, None]).sum(0) / w[m].sum()
        errs.append(float(np.abs(ref - YF[i]).max()))
    return dict(n_checked=len(pick), max_abs_diff=float(np.max(errs)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dataset', default='tahoe', choices=['tahoe'])
    ap.add_argument('--seeds', type=int, default=N_SPLIT_SEEDS)
    ap.add_argument('--boot', type=int, default=N_BOOT)
    ap.add_argument('--shard-root', default=None)
    ap.add_argument('--sensitivity', action='store_true',
                    help='also compute the SECONDARY sensitivities S1-S3 (results under the "sensitivity" key)')
    ap.add_argument('--line-boot', type=int, default=S3_BOOT)
    ap.add_argument('--out', default=str(RESULTS / 'rq2_params_tahoe.json'))
    ap.add_argument('--h2b-out', default=str(RESULTS / 'rq2_h2b.json'))
    a = ap.parse_args()
    sys.exit(run_tahoe(a))


if __name__ == '__main__':
    main()
