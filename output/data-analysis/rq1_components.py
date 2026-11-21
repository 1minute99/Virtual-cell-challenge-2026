"""RQ1 Task 12 -- noise-corrected variance components and ceilings on the 7 harness panels (library + CLI).

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 12, operational definitions 4 and 6, H1
(this script computes NO H1 quantity: no method prediction is read; Task 13 does that).

MODEL (Notation of the plan): delta_{c,t} = mu + a_t + b_c + i_{c,t} + eps_{c,t}, per gene.
For a destination panel c0 with target set T (panel targets measured, >= 20 cells, in >= 2 OTHER lines):
  q_t      = lambda + alpha_t, the OTHER lines' "average-line" additive profile of target t: an unweighted
             two-way fit delta_{k,t} = c + alpha_t + beta_k per gene over the other lines' (k, t) rows (t in T),
             completed over all lines and averaged over k (gauge-invariant).  Uses NO destination data.
  mu       = mean_t q_t                      (shared response = other lines' mean over T)
  a_t      = q_t - mu                        (conserved perturbation component, from the OTHER lines' truth)
  b_c0     = mean_t delta_{c0,t} - mu        (destination cell-line component; the oracle's knowledge of b)
  i_{c0,t} = delta_{c0,t} - mu - b_c0 - a_t  (line x perturbation interaction incl. the finite-K leakage of
                                              the other lines' interactions into a_t)
  Oracle-transfer predictor (Task 13): mu + a_t + b_c0.
Per-element (target x gene) noise-free quantities, from cell-level split halves A/B of EVERY line with SEPARATE
control halves (noise independent between halves => <X^A, X^B> is unbiased for ||X||^2 of the noise-free X):
  Var(a)  = <a^A, a^B>/N        Var(b) = <b^A, b^B>/G        Var(i) = <i^A, i^B>/N
  C_ai    = (<a^A,i^B> + <i^A,a^B>)/(2N)   (a and i are NOT orthogonal by construction: a comes from other lines)
  Var_rep = <(delta-mu)^A, (delta-mu)^B>/N  = Var(a) + Var(b) + Var(i) + 2 C_ai   (exact identity)
  sigma2_eps (destination, full panel depth) = (1/4) mean (delta^A - delta^B)^2     (A1: half = 2x full noise)
  sigma2_eps_i (definition-4 analogue)       = (1/4) mean (i^A - i^B)^2             (incl. other lines' noise)
  N = |T| * G.
CEILINGS (correlation scale; PRIMARY space = mu-centred response delta - mu over all target x gene elements):
  measurement ceiling  rel = Var_rep / (Var_rep + sigma2_eps),  ceiling_meas = sqrt(rel)   (max r of ANY
                       predictor against the full-depth noisy truth); also the split-half r_hh and its
                       Spearman-Brown value 2 r_hh/(1 + r_hh) as a check.
  r_max      = sqrt[(Var(a) + Var(b)) / Var_rep]                     (definition 6, as registered)
  r_oracle   = (Var(a) + Var(b) + C_ai) / sqrt[(Var(a) + Var(b)) Var_rep]   (noise-free r of mu + a + b itself)
  r_proj     = sqrt[((Var(a) + C_ai)^2 / Var(a) + Var(b)) / Var_rep]  (best rescaling of a and b; >= r_oracle)
  r_max_orth = sqrt[(Var(a) + Var(b)) / (Var(a) + Var(b) + Var(i))]   (Molina-Zhang convention: C_ai dropped)
  r_max * ceiling_meas = expected r of the oracle against the noisy full-depth truth.
  SECONDARY space 'pearson' = plain Pearson over elements (scalar grand mean removed; mu counts as transferable).
UNCERTAINTY: the Gram matrices of [delta_t; q_t; 1] across halves are averaged over 5 split seeds, then targets
are resampled with replacement (B = 2,000, percentile CIs); every quantity (centring, b, a, i) is recomputed per
draw by closed-form quadratic forms. The other lines' additive fit (q_t) is held fixed (as A5 of Task 8).

ASSUMPTIONS (stated; checked where possible -- knowledge/synthesis/rq1_components-assumptions.md):
  B1 Truth of the destination = the harness truth file (common.open_panel, which="truth") = what is scored; all
     its cells, controls = 'non-targeting'. Library = cell total over ALL panel genes (as Task 14/15).
  B2 Other lines = cell-level pseudobulks of the Task 9 cache (results/rq2_gene_cache, read-only; universe U)
     plus results/rq1_gene_cache (this task's extraction of the panel targets outside U, same readers, same
     gene axis, control split aligned with Task 9 and verified identical). Library = native cell total (Task 9 E2).
     A (line, target) is usable if the group has >= 20 cells and every half has a positive library.
  B3 Gene set G = official-axis genes in G_cap (results/g_cap.json; all 7 panels, >= 1 control UMI) AND on the
     cache axis (measured in every source). First panel column per official gene (as Task 15 A3).
  B4 Response = log1p CP10K(target pseudobulk) - log1p CP10K(control pseudobulk) of the SAME half (Task 9 / Task 14
     convention). FENG (pooled HipSci donors) = one line: per-donor responses averaged with weights
     1/(1/n_t + 1/n_c) (the half's own counts); for the KOLF destination kolf_2/kolf_3 donors are dropped (A3).
     CD4 halves are unit-level (Task 9 E4 deviation).
  B5 Noise independent between halves (disjoint cells, disjoint control halves); sampling noise ~ 1/cells (A1 of
     Task 8). The O(1/n) bias of log pseudobulks is common to both halves and not removed (small; flagged).
  B6 Other-line set PRIMARY 'pack' = the panel's ERA pack families (the pipeline's sources; RNJURKAT->JURKAT,
     H1->VCC_H1); SENSITIVITY 'harness6' = the other 6 harness lines; 'rn3' (RN panels only) = the other 3 of
     K562/RPE1/HepG2/Jurkat; 'rn4_insample' (RN panels only, CITE-CHECK, not a transfer ceiling) = the 4 RN lines
     with the destination INSIDE the additive fit, as in Molina & Zhang's 4-line decomposition (their Fig. 1c
     implies r_max = sqrt((27.8 + 29.4)/80.7) = 0.84).
  B8 The registered r_max formula equals the oracle's correlation only if a and i are orthogonal. With a from the
     OTHER lines they are not (C_ai != 0, e.g. a destination that responds more strongly than the other lines'
     average): then r_oracle (the oracle's own noise-free r) and r_proj (best rescaled oracle = the correlation
     'achievable by an oracle that knows a_t and b_c', definition 6 in words) differ from r_max. All three are
     reported; which one Task 13 uses is not decided here (flagged).
  B7 Unweighted additive fit; targets and lines enter equally (noise-free cross-products are unbiased anyway).
Seeds: destination split default_rng([20261112, seed, panel_id, group_id]); other lines: Task 9's splits;
bootstrap 20261113. No scorer runs, no method predictions are read.
CLI: python rq1_components.py run [--panels ...] [--boot 2000]   -> results/rq1_components.parquet, .json, npz
Environment: ~/miniconda3/envs/ct_transfer; nice -n 10; OMP/MKL threads 4.
"""
from __future__ import annotations

from pathlib import Path
import argparse, json, sys, time

import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC')
DA = ROOT / 'output/data-analysis'
RESULTS = DA / 'results'
T9CACHE = RESULTS / 'rq2_gene_cache'
T12CACHE = RESULTS / 'rq1_gene_cache'
OUTDIR = RESULTS / 'rq1_components'

DEST_BASE = 20261112
BOOT_SEED = 20261113
N_SEEDS = 5
N_BOOT = 2000
MIN_OTHER = 2
DEST_LINE = '__dest__'
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
DEST_SRC = {'H1': 'VCC_H1', 'OR_HCT116X': 'HCT116', 'KOLF': 'KOLF', 'RN_K562': 'K562', 'RN_HepG2': 'HEPG2',
            'RN_Jurkat': 'JURKAT', 'RN_RPE1': 'RPE1'}
FAMILY_SRC = {'RNJURKAT': 'JURKAT', 'H1': 'VCC_H1'}
FENG_DROP = {'KOLF': ('kolf_2', 'kolf_3')}
RN4 = ['K562', 'RPE1', 'HEPG2', 'JURKAT']
MZ_IMPLIED = dict(b=27.8, a=29.4, i=23.5, noise=19.3,
                  r_max=float(np.sqrt((27.8 + 29.4) / (27.8 + 29.4 + 23.5))),
                  source='input/notes/overlap-molina-zhang.md (Fig. 1c / Supp. Note 1)')


def log(**k):
    print(json.dumps(k, default=str), flush=True)


# =========================================================================== generic estimator (library)

def lcp(S, L):
    S = np.asarray(S, np.float64); L = np.asarray(L, np.float64)
    if np.any(L <= 0):
        raise ValueError('pseudobulk with non-positive library')
    return np.log1p(S * 1e4 / L[..., None])


def additive_design(line_idx, targ_idx, T, K):
    """pinv of [1, target dummies, line dummies] on the observed rows; checks connectivity (rank T + K - 1)."""
    n = len(line_idx)
    X = np.zeros((n, 1 + T + K))
    X[:, 0] = 1
    X[np.arange(n), 1 + targ_idx] = 1
    X[np.arange(n), 1 + T + line_idx] = 1
    rank = np.linalg.matrix_rank(X)
    if rank != T + K - 1:
        raise ValueError(f'additive design not connected: rank {rank} != T + K - 1 = {T + K - 1}')
    return np.linalg.pinv(X)


def target_profile(Y, Pinv, T, K):
    """q_t = c + alpha_t + mean_k beta_k (gauge-invariant 'average other line' profile), (T x G)."""
    coef = Pinv @ np.asarray(Y, np.float64)
    return coef[0] + coef[1:1 + T] + coef[1 + T:].mean(0)


def grams(dA, dB, qA, qB):
    """Cross- and same-half Gram matrices of Z^h = [delta_t (T rows); q_t (T rows); 1]  ((2T+1) square)."""
    G = dA.shape[1]
    ZA = np.vstack([dA, qA, np.ones((1, G))]); ZB = np.vstack([dB, qB, np.ones((1, G))])
    return dict(AB=ZA @ ZB.T, AA=ZA @ ZA.T, BB=ZB @ ZB.T, G=G, T=dA.shape[0])


def _blocks(M, T):
    d, q, o = slice(0, T), slice(T, 2 * T), 2 * T
    return dict(dd=M[d, d], dq=M[d, q], qd=M[q, d], qq=M[q, q], d1=M[d, o], one_d=M[o, d])


def _S(Mxy, W, T):
    """sum_t w_t <x_t^A - xbar^A, y_t^B - ybar^B>, xbar = sum w x / T; W (B x T) -> (B,)."""
    return W @ np.diag(Mxy) - np.einsum('bt,ts,bs->b', W, Mxy, W) / T


def _qf(M, W):
    return np.einsum('bt,ts,bs->b', W, M, W)


def components(gr, W=None):
    """All quantities for target weights W (B x T; default all ones). Returns dict of (B,) arrays."""
    T, G = gr['T'], gr['G']
    W = np.ones((1, T)) if W is None else np.atleast_2d(np.asarray(W, np.float64))
    N = T * G
    X = {h: _blocks(gr[h], T) for h in ('AB', 'AA', 'BB')}
    out = {}

    def parts(b):
        pp = b['dd'] - b['dq'] - b['qd'] + b['qq']                  # <p^A, p^B>, p = delta - q
        qp = b['qd'] - b['qq']; pq = b['dq'] - b['qq']
        Va = _S(b['qq'], W, T) / N
        Vi = _S(pp, W, T) / N
        Cai = 0.5 * (_S(qp, W, T) + _S(pq, W, T)) / N
        Vb = _qf(pp, W) / T ** 2 / G
        Vmu = _qf(b['qq'], W) / T ** 2 / G
        Vrep = (W @ np.diag(b['dd']) - _qf(b['dq'] + b['qd'], W) / T + _qf(b['qq'], W) / T) / N
        # Pearson space: scalar grand mean m^h of delta^h
        sA = b['d1']; sB = b['one_d']
        mA = (W @ sA) / N; mB = (W @ sB) / N
        VrepP = (W @ np.diag(b['dd']) - mB * (W @ sA) - mA * (W @ sB) + N * mA * mB) / N
        piP = (_qf(b['dd'], W) / T ** 2 - mB * (W @ sA) / T - mA * (W @ sB) / T + G * mA * mB) / G   # ||dbar - m||^2/G
        return dict(Va=Va, Vi=Vi, Cai=Cai, Vb=Vb, Vmu=Vmu, Vrep=Vrep, VrepP=VrepP, VtrP=piP + Va, pp=pp)

    ab = parts(X['AB']); aa = parts(X['AA']); bb = parts(X['BB'])
    Va, Vb, Vi, Cai, Vrep = ab['Va'], ab['Vb'], ab['Vi'], ab['Cai'], ab['Vrep']
    diff2 = np.diag(X['AA']['dd']) + np.diag(X['BB']['dd']) - 2 * np.diag(X['AB']['dd'])
    sig2 = (W @ diff2) / (4 * N)
    sig2_i = (_S(aa['pp'], W, T) + _S(bb['pp'], W, T) - 2 * _S(ab['pp'], W, T)) / (4 * N)
    rel = Vrep / (Vrep + sig2)
    r_hh = Vrep / np.sqrt(aa['Vrep'] * bb['Vrep'])
    relP = ab['VrepP'] / (ab['VrepP'] + sig2)
    Vt = Va + Vb
    with np.errstate(invalid='ignore', divide='ignore'):
        out.update(
            Var_mu=ab['Vmu'], Var_a=Va, Var_b=Vb, Var_i=Vi, C_ai=Cai, Var_rep=Vrep,
            identity_resid=Vrep - (Va + Vb + Vi + 2 * Cai),
            sigma2_eps=sig2, sigma2_eps_i=sig2_i, noise_to_rep=sig2 / Vrep,
            frac_a=Va / Vrep, frac_b=Vb / Vrep, frac_i=Vi / Vrep, frac_2Cai=2 * Cai / Vrep,
            frac_total_a=Va / (Vrep + sig2), frac_total_b=Vb / (Vrep + sig2), frac_total_i=Vi / (Vrep + sig2),
            frac_total_2Cai=2 * Cai / (Vrep + sig2), frac_total_noise=sig2 / (Vrep + sig2),
            reliability=rel, ceiling_meas=np.sqrt(rel), r_halfhalf=r_hh,
            reliability_SB=2 * r_hh / (1 + r_hh), ceiling_meas_SB=np.sqrt(2 * r_hh / (1 + r_hh)),
            r_max=np.sqrt(Vt / Vrep),
            r_max_orth=np.sqrt(Vt / (Vt + Vi)),
            r_oracle=(Vt + Cai) / np.sqrt(Vt * Vrep),
            r_proj=np.sqrt(((Va + Cai) ** 2 / Va + Vb) / Vrep),
            r_max_vs_noisy=np.sqrt(Vt / Vrep) * np.sqrt(rel),
            amp_slope_a=(Va + Cai) / Va,
            Var_rep_pearson=ab['VrepP'], Var_transferable_pearson=ab['VtrP'],
            r_max_pearson=np.sqrt(ab['VtrP'] / ab['VrepP']),
            r_oracle_pearson=(ab['VtrP'] + Cai) / np.sqrt(ab['VtrP'] * ab['VrepP']),
            reliability_pearson=relP, ceiling_meas_pearson=np.sqrt(relP),
        )
    return out


KEYS_CI = ['Var_mu', 'Var_a', 'Var_b', 'Var_i', 'C_ai', 'Var_rep', 'sigma2_eps', 'sigma2_eps_i', 'noise_to_rep',
           'frac_a', 'frac_b', 'frac_i', 'frac_2Cai', 'frac_total_a', 'frac_total_b', 'frac_total_i',
           'frac_total_2Cai', 'frac_total_noise', 'reliability', 'ceiling_meas', 'r_halfhalf', 'reliability_SB',
           'ceiling_meas_SB', 'r_max', 'r_max_orth', 'r_oracle', 'r_proj', 'r_max_vs_noisy', 'amp_slope_a', 'Var_rep_pearson',
           'Var_transferable_pearson', 'r_max_pearson', 'r_oracle_pearson', 'reliability_pearson',
           'ceiling_meas_pearson', 'identity_resid']


def mean_grams(grs):
    return dict(AB=np.mean([g['AB'] for g in grs], 0), AA=np.mean([g['AA'] for g in grs], 0),
                BB=np.mean([g['BB'] for g in grs], 0), G=grs[0]['G'], T=grs[0]['T'])


N_GENE_BLOCKS = 20
GENE_BLOCK_SEED = 20261114


def gene_blocks(G, n_blocks=N_GENE_BLOCKS, seed=GENE_BLOCK_SEED):
    """Fixed random partition of the gene columns into n_blocks blocks."""
    return np.array_split(np.random.default_rng(seed).permutation(G), n_blocks)


def block_grams(dA, dB, qA, qB, blocks):
    return [grams(dA[:, b], dB[:, b], qA[:, b], qB[:, b]) for b in blocks]


def sum_grams(grs):
    return dict(AB=sum(g['AB'] for g in grs), AA=sum(g['AA'] for g in grs), BB=sum(g['BB'] for g in grs),
                G=sum(g['G'] for g in grs), T=grs[0]['T'])


def gene_jackknife(blocks_gr, z=1.959964):
    """Delete-one-gene-block jackknife (captures gene sampling AND control noise, which is shared by all targets
    and invisible to the target bootstrap). Returns {key: (lo, hi, se)}."""
    tot = sum_grams(blocks_gr)
    J = len(blocks_gr)
    th = {k: [] for k in KEYS_CI}
    for g in blocks_gr:
        c = components(dict(AB=tot['AB'] - g['AB'], AA=tot['AA'] - g['AA'], BB=tot['BB'] - g['BB'],
                            G=tot['G'] - g['G'], T=tot['T']))
        for k in KEYS_CI:
            th[k].append(float(c[k][0]))
    full = {k: float(v[0]) for k, v in components(tot).items()}
    out = {}
    for k in KEYS_CI:
        v = np.array(th[k])
        se = float(np.sqrt((J - 1) / J * ((v - v.mean()) ** 2).sum()))
        out[k] = (full[k] - z * se, full[k] + z * se, se)
    return out


def bootstrap(gr, n_boot=N_BOOT, seed=BOOT_SEED, alpha=0.05, chunk=250):
    """Percentile CIs over targets (multinomial weights), all quantities recomputed per draw."""
    T = gr['T']
    rng = np.random.default_rng(seed)
    res = {k: [] for k in KEYS_CI}
    for lo in range(0, n_boot, chunk):
        b = min(chunk, n_boot - lo)
        W = np.stack([np.bincount(rng.integers(0, T, T), minlength=T) for _ in range(b)]).astype(np.float64)
        c = components(gr, W)
        for k in KEYS_CI:
            res[k].append(c[k])
    out = {}
    for k in KEYS_CI:
        v = np.concatenate(res[k])
        v = v[np.isfinite(v)]
        if len(v) == 0:
            out[k] = (float('nan'), float('nan'), 0); continue
        out[k] = (float(np.quantile(v, alpha / 2)), float(np.quantile(v, 1 - alpha / 2)), int(len(v)))
    return out


# =========================================================================== data: other lines (caches)

def gene_set():
    """Official-axis indices: G_cap AND the cache axis. Returns (official idx array, cache column positions)."""
    import common as CM
    off = CM.official_genes()
    first = {}
    for i, g in enumerate(off):
        first.setdefault(g, i)
    gc = json.loads((RESULTS / 'g_cap.json').read_text())['genes']
    gcap = {first[g] for g in gc}
    axis = json.loads((T9CACHE / 'plan.json').read_text())['gene_axis']
    ax12 = json.loads((T12CACHE / 'plan.json').read_text())['gene_axis']
    assert axis == ax12, 'Task 12 cache must share the Task 9 gene axis'
    pos = {g: j for j, g in enumerate(axis)}
    G = np.array(sorted(gcap & set(axis)))
    return G, np.array([pos[g] for g in G])


def _load_cache(path, cols):
    with np.load(path, allow_pickle=False) as z:
        groups = [str(x) for x in z['groups']]
        return dict(groups=groups, gi={g: i for i, g in enumerate(groups)}, n=z['n'], nA=z['nA'],
                    lib=z['lib'], libA=z['libA'], full=z['full'][:, cols], A=z['A'][:, :, cols].astype(np.float64))


def _resp(c, g, ctrl):
    """full (G,), A (S,G), B (S,G) responses of group g vs control group ctrl in cache c; None if a half is empty."""
    i, j = c['gi'][g], c['gi'][ctrl]
    libB_i = c['lib'][i] - c['libA'][:, i]; libB_j = c['lib'][j] - c['libA'][:, j]
    if np.any(c['libA'][:, i] <= 0) or np.any(libB_i <= 0) or np.any(c['libA'][:, j] <= 0) or np.any(libB_j <= 0):
        return None
    full = lcp(c['full'][i][None], c['lib'][i:i + 1])[0] - lcp(c['full'][j][None], c['lib'][j:j + 1])[0]
    A = lcp(c['A'][:, i], c['libA'][:, i]) - lcp(c['A'][:, j], c['libA'][:, j])
    Bi = c['full'][i][None] - c['A'][:, i]; Bj = c['full'][j][None] - c['A'][:, j]
    B = lcp(Bi, libB_i) - lcp(Bj, libB_j)
    w = dict(full=1 / (1 / c['n'][i] + 1 / c['n'][j]),
             A=1 / (1 / c['nA'][:, i] + 1 / c['nA'][:, j]),
             B=1 / (1 / (c['n'][i] - c['nA'][:, i]) + 1 / (c['n'][j] - c['nA'][:, j])))
    return dict(full=full, A=A, B=B, w=w)


def source_responses(src, targets, cols, drop_donors=()):
    """{canonical target: dict(full, A, B)} for the targets this source measures (>= 20 cells, both halves)."""
    tset = set(targets)
    out, dropped = {}, []
    paths = [T9CACHE / f'{src}.npz'] + ([T12CACHE / f'{src}.npz'] if (T12CACHE / f'{src}.npz').exists() else [])
    for path in paths:
        c = _load_cache(path, cols)
        if src == 'FENG':
            per = {}
            for g in c['groups']:
                line, t = g.split('|', 1)
                if t == 'CTRL' or t not in tset or line in drop_donors:
                    continue
                r = _resp(c, g, f'{line}|CTRL')
                if r is None:
                    dropped.append(g); continue
                per.setdefault(t, []).append(r)
            for t, rs in per.items():
                agg = {}
                for h in ('full', 'A', 'B'):
                    ws = np.array([r['w'][h] for r in rs], np.float64)          # (D,) or (D, S)
                    xs = np.stack([r[h] for r in rs])                            # (D, G) or (D, S, G)
                    agg[h] = (ws[..., None] * xs).sum(0) / ws.sum(0)[..., None]
                agg['n_donors'] = len(rs)
                out[t] = agg
        else:
            for g in c['groups']:
                if g == 'CTRL' or g not in tset or g in out:
                    continue
                r = _resp(c, g, 'CTRL')
                if r is None:
                    dropped.append(g); continue
                out[g] = r
        del c
    return out, dropped


# =========================================================================== data: destination panel

def dest_sums(panel, G_off, seeds=N_SEEDS, chunk=20000):
    import common as CM
    from scipy import sparse
    canon, _ = CM._canon()
    a = CM.open_panel(panel)
    try:
        lab = a.obs['target'].astype(str).to_numpy()
        om = CM.to_official(a.var_names.astype(str).to_numpy())
        pos = {}
        for j, o in enumerate(om):
            if o >= 0 and o not in pos:
                pos[int(o)] = j
        cols = np.array([pos[int(g)] for g in G_off])
        tg = sorted(set(lab) - {'non-targeting'})
        names = tg + ['CTRL']
        gid = {t: i for i, t in enumerate(tg)}; gid['non-targeting'] = len(tg)
        code = np.array([gid[x] for x in lab])
        ng, n = len(names), len(lab)
        pid = PANELS7.index(panel)
        inA = np.zeros((seeds, n), bool)
        for g in range(ng):
            ii = np.flatnonzero(code == g)
            for s in range(seeds):
                perm = np.random.default_rng([DEST_BASE, s, pid, g]).permutation(len(ii))
                inA[s, ii] = perm < len(ii) // 2
        full = np.zeros((ng, len(cols))); A = np.zeros((seeds, ng, len(cols)))
        lib = np.zeros(ng); libA = np.zeros((seeds, ng))
        for r0 in range(0, n, chunk):
            r1 = min(n, r0 + chunk)
            X = sparse.csr_matrix(a.X[r0:r1], dtype=np.float64)
            cl = np.asarray(X.sum(1)).ravel()
            Xg = X[:, cols].toarray()
            c = code[r0:r1]
            I = sparse.csr_matrix((np.ones(r1 - r0), (c, np.arange(r1 - r0))), shape=(ng, r1 - r0))
            full += I @ Xg; lib += np.bincount(c, weights=cl, minlength=ng)
            for s in range(seeds):
                m = inA[s, r0:r1]
                Ia = sparse.csr_matrix((np.ones(int(m.sum())), (c[m], np.flatnonzero(m))), shape=(ng, r1 - r0))
                A[s] += Ia @ Xg; libA[s] += np.bincount(c[m], weights=cl[m], minlength=ng)
        nn = np.bincount(code, minlength=ng)
        nA = np.stack([np.bincount(code[inA[s]], minlength=ng) for s in range(seeds)])
    finally:
        a.file.close()
    ctarget = [canon(t) for t in tg]
    if len(set(ctarget)) != len(ctarget):
        raise ValueError(f'{panel}: canonical target names collide')
    return dict(labels=tg, targets=ctarget, n=nn, nA=nA, full=full, A=A, lib=lib, libA=libA, cols=cols)


def dest_responses(ds):
    ng = len(ds['labels']) + 1; c = ng - 1
    full = lcp(ds['full'][:c], ds['lib'][:c]) - lcp(ds['full'][c:], ds['lib'][c:])
    A = lcp(ds['A'][:, :c], ds['libA'][:, :c]) - lcp(ds['A'][:, c:], ds['libA'][:, c:])
    B = (lcp(ds['full'][None, :c] - ds['A'][:, :c], ds['lib'][None, :c] - ds['libA'][:, :c])
         - lcp(ds['full'][None, c:] - ds['A'][:, c:], ds['lib'][None, c:] - ds['libA'][:, c:]))
    return full, A, B


# =========================================================================== panel driver

def other_sets(panel):
    import common as CM
    pj = json.loads((RESULTS / 'panels.json').read_text())['panels'][panel]
    excl = set(pj['excluded_sources']['all'])
    fam = [f for f in pj['pack_families'] if f not in excl]
    pack = [FAMILY_SRC.get(f, f) for f in fam]
    assert DEST_SRC[panel] not in pack
    sets = dict(pack=pack, harness6=[DEST_SRC[p] for p in PANELS7 if p != panel])
    if panel.startswith('RN_'):
        sets['rn3'] = [s for s in RN4 if s != DEST_SRC[panel]]
        sets['rn4_insample'] = sets['rn3'] + [DEST_LINE]     # cite-check only: destination INSIDE the fit (MZ-like)
    return sets


def run_set(panel, set_name, srcs, dest, resp_cache, G_off, n_boot, save_arrays):
    dfull, dA, dB = dest
    ds = resp_cache['__dest__']
    targets = ds['targets']
    avail = {s: resp_cache[s][0] for s in srcs if s != DEST_LINE}
    if DEST_LINE in srcs:                                   # in-sample variant: the destination is one of the lines
        avail[DEST_LINE] = {t: dict(full=dfull[j], A=dA[:, j], B=dB[:, j]) for j, t in enumerate(targets)}
    cnt = np.array([sum(1 for s in srcs if s != DEST_LINE and t in avail[s]) for t in targets])
    okd = (ds['n'][:-1] >= 2) & (ds['nA'][:, :-1] >= 1).all(0) & ((ds['n'][None, :-1] - ds['nA'][:, :-1]) >= 1).all(0)
    keep = np.flatnonzero((cnt >= MIN_OTHER) & okd)
    T = len(keep); K = len(srcs)
    if T < 10:
        log(event='skip', panel=panel, set=set_name, n_targets=int(T)); return None
    tk = [targets[i] for i in keep]
    rows = [(k, j) for j, t in enumerate(tk) for k, s in enumerate(srcs) if t in avail[s]]
    used_lines = sorted({srcs[k] for k, _ in rows})
    li = np.array([used_lines.index(srcs[k]) for k, _ in rows]); ti = np.array([j for _, j in rows])
    Kused = len(used_lines)
    Pinv = additive_design(li, ti, T, Kused)
    getr = lambda h, s=None: np.stack([avail[srcs[k]][tk[j]][h] if s is None else avail[srcs[k]][tk[j]][h][s]
                                       for k, j in rows])
    q_full = target_profile(getr('full'), Pinv, T, Kused)
    grs, per_seed, halves = [], [], []
    blocks = gene_blocks(dA.shape[2])
    bsum = None
    for s in range(N_SEEDS):
        qA = target_profile(getr('A', s), Pinv, T, Kused); qB = target_profile(getr('B', s), Pinv, T, Kused)
        bg = block_grams(dA[s][keep], dB[s][keep], qA, qB, blocks)
        g = sum_grams(bg)
        grs.append(g)
        bsum = bg if bsum is None else [dict(AB=x['AB'] + y['AB'], AA=x['AA'] + y['AA'], BB=x['BB'] + y['BB'],
                                             G=x['G'], T=x['T']) for x, y in zip(bsum, bg)]
        per_seed.append({k: float(v[0]) for k, v in components(g).items()})
        if save_arrays:
            halves.append((dA[s][keep].astype(np.float32), dB[s][keep].astype(np.float32),
                           qA.astype(np.float32), qB.astype(np.float32)))
    gm = mean_grams(grs)
    point = {k: float(v[0]) for k, v in components(gm).items()}
    ci = bootstrap(gm, n_boot=n_boot)
    jk = gene_jackknife([dict(AB=x['AB'] / N_SEEDS, AA=x['AA'] / N_SEEDS, BB=x['BB'] / N_SEEDS, G=x['G'], T=x['T'])
                         for x in bsum])
    lines_per_target = np.bincount(ti, minlength=T)
    meta = dict(panel=panel, other_set=set_name, n_targets=int(T), n_panel_targets=len(targets),
                n_genes=int(gm['G']), lines=used_lines, n_lines=Kused, n_rows=len(rows),
                lines_per_target_median=float(np.median(lines_per_target)),
                lines_per_target_min=int(lines_per_target.min()),
                dest_cells_per_target_median=float(np.median(ds['n'][:-1][keep])),
                dest_ctrl_cells=int(ds['n'][-1]),
                line_target_counts={s: int(sum(1 for t in tk if t in avail[s])) for s in used_lines})
    if save_arrays:
        mu = q_full.mean(0); b = dfull[keep].mean(0) - mu
        OUTDIR.mkdir(parents=True, exist_ok=True)
        np.savez(OUTDIR / f'{panel}_{set_name}.npz', targets=np.array(tk), labels=np.array([ds['labels'][i] for i in keep]),
                 genes_official_idx=G_off, panel_cols=ds['cols'], delta_full=dfull[keep].astype(np.float32),
                 q_full=q_full.astype(np.float32), mu=mu.astype(np.float32), b=b.astype(np.float32),
                 a_full=(q_full - mu).astype(np.float32),
                 dA=np.stack([h[0] for h in halves]), dB=np.stack([h[1] for h in halves]),
                 qA=np.stack([h[2] for h in halves]), qB=np.stack([h[3] for h in halves]),
                 note=np.array('dest truth halves (dA,dB) and other-line target profiles q (qA,qB) per split seed; '
                               'mu = mean_t q, a = q - mu, b = mean_t delta - mu (full depth); oracle = mu + a + b'))
    log(event='set_done', **{k: v for k, v in meta.items() if k != 'line_target_counts'},
        r_max=point['r_max'], ceiling_meas=point['ceiling_meas'], r_oracle=point['r_oracle'])
    return dict(meta=meta, point=point, ci=ci, jk=jk, per_seed=per_seed)


def run_panel(panel, G_off, cols, n_boot):
    t0 = time.time()
    ds = dest_sums(panel, G_off)
    dest = dest_responses(ds)
    log(event='dest', panel=panel, targets=len(ds['targets']), ctrl=int(ds['n'][-1]), secs=round(time.time() - t0, 1))
    sets = other_sets(panel)
    srcs_all = sorted(set().union(*sets.values()) - {DEST_LINE})
    rc = {'__dest__': ds}
    for s in srcs_all:
        rc[s] = source_responses(s, ds['targets'], cols, drop_donors=FENG_DROP.get(panel, ()) if s == 'FENG' else ())
        log(event='source', panel=panel, src=s, n_targets=len(rc[s][0]), dropped=len(rc[s][1]),
            secs=round(time.time() - t0, 1))
    out = {}
    for name, srcs in sets.items():
        r = run_set(panel, name, srcs, dest, rc, G_off, n_boot, save_arrays=(name == 'pack'))
        if r is not None:
            r['meta']['dropped_groups'] = {s: rc[s][1] for s in srcs if s != DEST_LINE and rc[s][1]}
            out[name] = r
    log(event='panel_done', panel=panel, secs=round(time.time() - t0, 1))
    return out


def to_rows(res):
    rows = []
    for panel, sets in res.items():
        for name, r in sets.items():
            m = r['meta']
            for k in KEYS_CI:
                ps = np.array([p[k] for p in r['per_seed']])
                lo, hi, nb = r['ci'][k]
                jlo, jhi, jse = r['jk'][k]
                rows.append(dict(panel=panel, other_set=name, primary=(name == 'pack'), quantity=k,
                                 estimate=r['point'][k], ci_lo=lo, ci_hi=hi, n_boot_finite=nb,
                                 genejk_lo=jlo, genejk_hi=jhi, genejk_se=jse,
                                 seed_min=float(ps.min()), seed_max=float(ps.max()),
                                 n_targets=m['n_targets'], n_genes=m['n_genes'], n_lines=m['n_lines'],
                                 lines=','.join(m['lines'])))
    return rows


def plot(tag=''):
    """Dot-and-interval chart per panel (primary 'pack' set): measurement ceiling, r_max, r_oracle, r_proj."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    df = pd.read_parquet(RESULTS / f'rq1_components{tag}.parquet')
    df = df[df.other_set == 'pack']
    series = [('ceiling_meas', 'measurement ceiling', '#2a78d6'), ('r_max', 'r_max (registered formula)', '#eb6834'),
              ('r_oracle', 'r_oracle (oracle mu+a+b)', '#1baf7a'), ('r_proj', 'r_proj (best-scaled oracle)', '#eda100')]
    panels = [p for p in PANELS7 if p in set(df.panel)]
    fig, ax = plt.subplots(figsize=(8.5, 4.2), dpi=150)
    off = np.linspace(-0.27, 0.27, len(series))
    for (q, lab, col), o in zip(series, off):
        d = df[df.quantity == q].set_index('panel').loc[panels]
        x = np.arange(len(panels)) + o
        ax.errorbar(x, d.estimate, yerr=[d.estimate - d.ci_lo, d.ci_hi - d.estimate], fmt='o', ms=5, color=col,
                    elinewidth=1.5, capsize=0, label=lab, mec='white', mew=1)
    ax.axhline(MZ_IMPLIED['r_max'], color='#52514e', lw=1, ls=(0, (4, 3)))
    ax.text(1.5, MZ_IMPLIED['r_max'] + 0.01, 'Molina & Zhang implied r_max 0.84 (in-sample, 4 lines)', ha='center',
            va='bottom', fontsize=8, color='#52514e')
    ax.set_xticks(np.arange(len(panels))); ax.set_xticklabels(panels, fontsize=9); ax.set_xlim(-0.5, len(panels) - 0.5)
    ax.set_ylabel('noise-corrected correlation (mu-centred)', fontsize=9)
    ax.set_ylim(0, 1.0); ax.grid(axis='y', color='#e6e5e1', lw=0.8); ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=8, frameon=False, ncol=2, loc='lower left')
    ax.set_title('Task 12: ceilings per harness panel (other lines = pack families; 95 % target-bootstrap CI)',
                 fontsize=10, loc='left')
    fig.tight_layout()
    fig.savefig(RESULTS / f'rq1_components{tag}.png'); fig.savefig(RESULTS / f'rq1_components{tag}.svg')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'plot'])
    ap.add_argument('--panels', nargs='*', default=PANELS7)
    ap.add_argument('--boot', type=int, default=N_BOOT)
    ap.add_argument('--tag', default='')
    a = ap.parse_args()
    sys.path.insert(0, str(DA))
    if a.cmd == 'plot':
        plot(a.tag); return
    t0 = time.time()
    G_off, cols = gene_set()
    log(event='genes', n=len(G_off))
    res = {}
    for p in a.panels:
        res[p] = run_panel(p, G_off, cols, a.boot)
    import pandas as pd
    df = pd.DataFrame(to_rows(res))
    tag = a.tag
    df.to_parquet(RESULTS / f'rq1_components{tag}.parquet', index=False)
    summary = dict(created=time.strftime('%Y-%m-%d %H:%M'), plan='input/ideas/zero-shot-limits-context-complexity-plan.md',
                   task=12, n_seeds=N_SEEDS, n_boot=a.boot, boot_seed=BOOT_SEED, dest_split_base=DEST_BASE,
                   n_genes=int(len(G_off)), molina_zhang_implied=MZ_IMPLIED,
                   panels={p: {n: dict(meta=r['meta'], point=r['point'],
                                       ci={k: r['ci'][k][:2] for k in KEYS_CI},
                                       gene_jackknife_ci={k: r['jk'][k][:2] for k in KEYS_CI}) for n, r in s.items()}
                           for p, s in res.items()},
                   runtime_s=round(time.time() - t0, 1))
    (RESULTS / f'rq1_components{tag}.json').write_text(json.dumps(summary, indent=1, default=float))
    plot(tag)
    log(event='done', secs=round(time.time() - t0, 1))


if __name__ == '__main__':
    main()
