"""Tests for output/data-analysis/rq1_components.py (Task 12, RQ1).

1. Known components (response level, balanced design, K other lines): Var(a), Var(b), Var(i), Var_rep,
   sigma2_eps, r_max and the measurement ceiling are recovered within 5 % of their population values
   (finite-K leakage of the other lines' interactions included analytically), with control noise shared by all
   targets of a line (separate control halves).
2. Unbalanced design: the split-half estimates match the noise-free plug-in (same decomposition on noise-free
   data) within 5 %.
3. Exact algebra: Var_rep = Var(a) + Var(b) + Var(i) + 2 C_ai; uniform bootstrap weights = point estimate;
   the quadratic-form bootstrap equals a brute-force recomputation on a resampled target set (fixed fit).
4. Gauge invariance of the target profile q_t (line order, design re-parametrisation); disconnected design raises.
5. Count level (Poisson cells -> log1p CP10K pseudobulks through the module's own response builders):
   Var_rep within 5 % of the noise-free value; sigma2_eps within 5 % of the variance between two independent
   full-depth replicates (depth rescaling A1).
6. Not-H1 guard: the module reads no prediction file and no cell-eval2 output.
"""
from pathlib import Path
import re, sys

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import rq1_components as R  # noqa: E402


def _simulate(seed=0, T=300, G=1500, K=6, sa=1.0, sb=0.6, si=0.5, sm=0.8, st=0.6, sc=0.3, miss=0.0):
    rng = np.random.default_rng(seed)
    L = K + 1                                   # line 0 = destination
    mu = sm * rng.standard_normal(G)
    a = sa * rng.standard_normal((T, G)); a -= a.mean(0)
    b = sb * rng.standard_normal((L, G))
    i = si * rng.standard_normal((L, T, G))
    clean = mu + a[None] + b[:, None] + i       # (L, T, G)
    # half-depth noise: target part var 2 st^2, control part (shared by all targets of a line) var 2 sc^2
    def half():
        return (np.sqrt(2) * st * rng.standard_normal((L, T, G))
                - np.sqrt(2) * sc * rng.standard_normal((L, 1, G)))
    YA, YB = clean + half(), clean + half()
    obs = np.ones((K, T), bool)
    if miss:
        obs = rng.random((K, T)) > miss
        for t in range(T):                       # keep >= 2 lines per target
            if obs[:, t].sum() < 2:
                obs[rng.choice(K, 2, replace=False), t] = True
    return dict(clean=clean, YA=YA, YB=YB, obs=obs, T=T, G=G, K=K, a=a, b=b, i=i, pars=dict(sa=sa, sb=sb, si=si, st=st, sc=sc))


def _estimate(sim, Y=None):
    T, K = sim['T'], sim['K']
    kk, tt = np.nonzero(sim['obs'])
    P = R.additive_design(kk, tt, T, K)
    YA = sim['YA'] if Y is None else Y
    YB = sim['YB'] if Y is None else Y
    qA = R.target_profile(YA[1:][kk, tt], P, T, K)
    qB = R.target_profile(YB[1:][kk, tt], P, T, K)
    gr = R.grams(YA[0], YB[0], qA, qB)
    return gr, {k: float(v[0]) for k, v in R.components(gr).items()}


def test_known_components_population_balanced():
    """Known (realised) components of the simulation, written out from the generating model, not from the module."""
    sim = _simulate()
    _, e = _estimate(sim)
    p, cl = sim['pars'], sim['clean']
    a, b, i = sim['a'], sim['b'], sim['i']
    ibar = i[1:].mean(0)                                   # other lines' mean interaction (T, G)
    a_hat = a + ibar - ibar.mean(0)                        # a from the OTHER lines (balanced fit), centred over t
    b_dest = b[0] - b[1:].mean(0) + i[0].mean(0) - ibar.mean(0)
    i_dest = (i[0] - i[0].mean(0)) - (ibar - ibar.mean(0))
    Va, Vb, Vi = np.mean(a_hat ** 2), np.mean(b_dest ** 2), np.mean(i_dest ** 2)
    Cai = np.mean(a_hat * i_dest)
    Vrep = np.mean((cl[0] - (cl[1:].mean(0)).mean(0)) ** 2)
    assert np.isclose(Vrep, Va + Vb + Vi + 2 * Cai)
    sig2 = p['st'] ** 2 + p['sc'] ** 2                     # full depth, target + control noise
    true = dict(Var_a=Va, Var_b=Vb, Var_i=Vi, Var_rep=Vrep, sigma2_eps=sig2,
                r_max=np.sqrt((Va + Vb) / Vrep), ceiling_meas=np.sqrt(Vrep / (Vrep + sig2)))
    for k, v in true.items():
        assert abs(e[k] / v - 1) < 0.05, (k, e[k], v)
    assert abs(e['C_ai'] - Cai) < 0.05 * Vrep
    # and the population expectations (finite-K leakage) within 10 % (finite G draw of b)
    K, T = sim['K'], sim['T']
    assert abs(Va / ((p['sa'] ** 2 + p['si'] ** 2 / K) * (1 - 1 / T)) - 1) < 0.10
    assert abs(Vi / (p['si'] ** 2 * (1 + 1 / K) * (1 - 1 / T)) - 1) < 0.10


def test_unbalanced_matches_noise_free_plugin():
    sim = _simulate(seed=3, miss=0.45)
    _, e = _estimate(sim)
    _, z = _estimate(sim, Y=sim['clean'])        # same decomposition on the noise-free tensor
    for k in ['Var_a', 'Var_b', 'Var_i', 'Var_rep', 'r_max', 'r_oracle', 'r_proj', 'r_max_pearson']:
        assert abs(e[k] / z[k] - 1) < 0.05, (k, e[k], z[k])
    assert abs(e['C_ai'] - z['C_ai']) < 0.05 * z['Var_rep']
    assert abs(z['sigma2_eps']) < 1e-12 and abs(e['sigma2_eps'] / 0.45 - 1) < 0.05   # 0.6^2 + 0.3^2


def test_identity_and_bootstrap_consistency():
    sim = _simulate(seed=5, T=60, G=300, miss=0.3)
    gr, e = _estimate(sim)
    assert abs(e['identity_resid']) < 1e-10 * e['Var_rep']
    c1 = R.components(gr, np.ones((1, gr['T'])))
    assert all(np.isclose(c1[k][0], e[k], equal_nan=True) for k in R.KEYS_CI)
    # quadratic-form bootstrap draw == brute force recomputation on the duplicated target set (fixed fit)
    rng = np.random.default_rng(1)
    idx = rng.integers(0, gr['T'], gr['T'])
    w = np.bincount(idx, minlength=gr['T']).astype(float)
    cw = R.components(gr, w[None])
    T, K = sim['T'], sim['K']
    kk, tt = np.nonzero(sim['obs'])
    P = R.additive_design(kk, tt, T, K)
    qA = R.target_profile(sim['YA'][1:][kk, tt], P, T, K)[idx]
    qB = R.target_profile(sim['YB'][1:][kk, tt], P, T, K)[idx]
    brute = R.components(R.grams(sim['YA'][0][idx], sim['YB'][0][idx], qA, qB))
    for k in ['Var_a', 'Var_b', 'Var_i', 'C_ai', 'Var_rep', 'sigma2_eps', 'sigma2_eps_i', 'r_max',
              'Var_rep_pearson', 'Var_transferable_pearson', 'r_halfhalf']:
        assert np.isclose(cw[k][0], brute[k][0], rtol=1e-8), (k, cw[k][0], brute[k][0])
    ci = R.bootstrap(gr, n_boot=200, seed=0)
    assert ci['r_max'][0] <= e['r_max'] <= ci['r_max'][1]


def test_gauge_invariance_and_connectivity():
    rng = np.random.default_rng(2)
    T, K, G = 20, 4, 7
    kk = np.repeat(np.arange(K), T); tt = np.tile(np.arange(T), K)
    keep = rng.random(len(kk)) > 0.3
    keep[:T] = True; keep[T:2 * T] = True          # lines 0, 1 measure everything (connected)
    kk, tt = kk[keep], tt[keep]
    Y = rng.standard_normal((len(kk), G))
    q = R.target_profile(Y, R.additive_design(kk, tt, T, K), T, K)
    perm = np.array([2, 0, 3, 1])                  # relabel lines
    q2 = R.target_profile(Y, R.additive_design(perm[kk], tt, T, K), T, K)
    assert np.allclose(q, q2)
    # balanced design: q_t = mean_k delta_kt
    kb = np.repeat(np.arange(K), T); tb = np.tile(np.arange(T), K); Yb = rng.standard_normal((K * T, G))
    qb = R.target_profile(Yb, R.additive_design(kb, tb, T, K), T, K)
    assert np.allclose(qb, Yb.reshape(K, T, G).mean(0))
    with pytest.raises(ValueError):               # two blocks of targets measured by disjoint lines
        R.additive_design(np.array([0, 1, 2, 3]), np.array([0, 0, 1, 1]), 2, 4)


def _poisson_line(rng, rates_t, rate_c, n_t, n_c):
    """Cells ~ Poisson(rate * size factor); returns per-group full / half-A sums and libraries."""
    T, G = rates_t.shape
    groups = list(range(T)) + ['C']
    full, A, lib, libA, n, nA = [], [], [], [], [], []
    for g in groups:
        r = rate_c if g == 'C' else rates_t[g]
        m = n_c if g == 'C' else n_t
        sf = rng.uniform(0.7, 1.3, m)
        X = rng.poisson(sf[:, None] * r[None])
        inA = rng.permutation(m) < m // 2
        full.append(X.sum(0)); A.append(X[inA].sum(0)); lib.append(X.sum()); libA.append(X[inA].sum())
        n.append(m); nA.append(int(inA.sum()))
    return dict(full=np.array(full, float), A=np.array(A, float)[None], lib=np.array(lib, float),
                libA=np.array(libA, float)[None], n=np.array(n), nA=np.array(nA)[None], labels=list(range(T)))


def test_count_level_end_to_end():
    rng = np.random.default_rng(11)
    T, G, K = 120, 400, 4
    base = np.exp(rng.normal(1.5, 1.0, G))                     # mean counts per cell
    a = 0.4 * rng.standard_normal((T, G)); b = 0.25 * rng.standard_normal((K + 1, G))
    i = 0.2 * rng.standard_normal((K + 1, T, G))
    logfc = a[None] + b[:, None] + i
    resp = {}
    for line in range(K + 1):
        ds = _poisson_line(rng, base[None] * np.exp(logfc[line]), base, 150, 3000)
        resp[line] = R.dest_responses(ds)                       # (full, A, B)
    # second independent full-depth destination replicate for the A1 check
    ds2 = _poisson_line(rng, base[None] * np.exp(logfc[0]), base, 150, 3000)
    full2 = R.dest_responses(ds2)[0]
    kk = np.repeat(np.arange(K), T); tt = np.tile(np.arange(T), K)
    P = R.additive_design(kk, tt, T, K)
    qA = R.target_profile(np.concatenate([resp[k + 1][1][0] for k in range(K)]), P, T, K)
    qB = R.target_profile(np.concatenate([resp[k + 1][2][0] for k in range(K)]), P, T, K)
    gr = R.grams(resp[0][1][0], resp[0][2][0], qA, qB)
    e = {k: float(v[0]) for k, v in R.components(gr).items()}
    # noise-free: expected log-profile responses (large-count limit) through the same decomposition
    clean = lambda line: (np.log1p(1e4 * (base * np.exp(logfc[line])) / (base * np.exp(logfc[line])).sum(1, keepdims=True))
                          - np.log1p(1e4 * base / base.sum())[None])
    q0 = R.target_profile(np.concatenate([clean(k + 1) for k in range(K)]), P, T, K)
    z = {k: float(v[0]) for k, v in R.components(R.grams(clean(0), clean(0), q0, q0)).items()}
    for k in ['Var_rep', 'Var_a', 'Var_i', 'r_max']:
        assert abs(e[k] / z[k] - 1) < 0.05, (k, e[k], z[k])
    rep_noise = 0.5 * np.mean((resp[0][0] - full2) ** 2)        # var between two full-depth replicates
    assert abs(e['sigma2_eps'] / rep_noise - 1) < 0.05, (e['sigma2_eps'], rep_noise)


def test_no_prediction_or_scorer_inputs():
    src = (HERE / 'rq1_components.py').read_text()
    for pat in [r'preds/', r'_prem3', r'_m3\.npz', r'cell[-_]eval', r'rq1_ceiling', r'sch_sumeval', r'arch_screen']:
        assert not re.search(pat, src), pat


def test_gene_jackknife_tracks_control_noise_in_b():
    """Control noise is shared by all targets of a line: the target bootstrap cannot see it in Var(b); the gene-block
    jackknife must. Fixed components, 24 noise replicates: jackknife SE within [0.6, 1.6] x the replicate SD."""
    base = _simulate(seed=21, T=100, G=800, sc=0.5)
    rng = np.random.default_rng(99)
    est, se_jk, se_bt = [], [], []
    for r in range(24):
        L, T, G = base['clean'].shape
        noise = lambda: (np.sqrt(2) * 0.6 * rng.standard_normal((L, T, G)) - np.sqrt(2) * 0.5 * rng.standard_normal((L, 1, G)))
        sim = dict(base, YA=base['clean'] + noise(), YB=base['clean'] + noise())
        T_, K = sim['T'], sim['K']
        kk, tt = np.nonzero(sim['obs'])
        P = R.additive_design(kk, tt, T_, K)
        qA = R.target_profile(sim['YA'][1:][kk, tt], P, T_, K); qB = R.target_profile(sim['YB'][1:][kk, tt], P, T_, K)
        bl = R.block_grams(sim['YA'][0], sim['YB'][0], qA, qB, R.gene_blocks(G))
        tot = R.sum_grams(bl)
        est.append(float(R.components(tot)['Var_b'][0]))
        se_jk.append(R.gene_jackknife(bl)['Var_b'][2])
        lo, hi, _ = R.bootstrap(tot, n_boot=200, seed=r)['Var_b']
        se_bt.append((hi - lo) / 3.92)
    sd = np.std(est, ddof=1)
    assert 0.6 < np.mean(se_jk) / sd < 1.6, (np.mean(se_jk), sd)
    assert np.mean(se_bt) < 0.5 * sd                 # documents the blind spot of the target bootstrap
