"""Tests for output/data-analysis/rq2_params.py (Task 8).

1. Blinding: rq2_params.py and every module it imports from output/ contain no reference to the Tahoe curve
   parquet, the dcp scaling outputs, or any other gain-vs-C result.
2. Estimators on synthetic data with known d, sig2_eps, Var(i), R^2: recovered within 10 %.
3. Depth rescaling (assumption A1) on a cell-level Poisson simulation: split-half sig2_eps (rescaled to full
   depth) matches the noise variance seen between two independent full-depth replicates within 10 %.
4. The streaming cell-level reader on a synthetic parquet shard in the Tahoe schema.
5. Real-data structure: descriptor identical to the Task 6 runner's; metadata-rebuilt row table identical to the
   model's tahoe.npz indices and weights.
"""
from pathlib import Path
import ast, math, re, sys

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import rq2_params as P  # noqa: E402

FORBIDDEN = [
    r"rq2_tahoe_curve[\w\-\*]*\.parquet",
    r"drug_context_pretrain_2026-09-20",
    r"scaling_table", r"scaling_rows", r"scaling/k",
    r"rq2_h2a", r"curves\.json", r"mse_gain", r"cos_gain", r"gain_vs_C",
]


def _output_imports(path):
    tree = ast.parse(Path(path).read_text())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name.split('.')[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split('.')[0])
    return [HERE / f'{n}.py' for n in sorted(names) if (HERE / f'{n}.py').exists()]


def test_blinding_grep():
    todo, seen, hits = [HERE / 'rq2_params.py'], set(), []
    while todo:
        f = todo.pop()
        if f in seen:
            continue
        seen.add(f)
        text = f.read_text()
        hits += [(f.name, pat) for pat in FORBIDDEN if re.search(pat, text)]
        todo += _output_imports(f)
    assert not hits, f'blinded inputs referenced: {hits}'


def test_blinding_vcc_helpers():
    """The vcc2026_data helper modules rq2_params imports must not touch gain files either."""
    for name in ('dcp_data', 'simple_residual_prepare', 'context_signal_tahoe_effects'):
        text = (P.SCRIPTS / f'{name}.py').read_text()
        assert not [p for p in FORBIDDEN if re.search(p, text)], name


# ----------------------------------------------------------------------------- synthetic estimator test

def _descriptor_with_spectrum(C, F, lam, rng):
    A = rng.standard_normal((C, len(lam)))
    A -= A.mean(0)
    U, _ = np.linalg.qr(A)
    V, _ = np.linalg.qr(rng.standard_normal((F, len(lam))))
    return 3.0 + U @ np.diag(np.sqrt(lam)) @ V.T          # offset: centring must remove it


def test_participation_ratio_known():
    rng = np.random.default_rng(0)
    lam = 0.8 ** np.arange(30)
    X = _descriptor_with_spectrum(48, 200, lam, rng)
    true = lam.sum() ** 2 / (lam ** 2).sum()
    assert abs(P.participation_ratio(X) - true) / true < 1e-8
    # PR is invariant to a rank-preserving projection + scalar (model-view equality)
    pcs, _ = P.descriptor_pcs(X)
    assert abs(P.participation_ratio(pcs / 7.3) - true) / true < 1e-8


def _synthetic_halves(seed, C=48, T=150, G=100, k_true=4, v_e=0.6, v_u=0.4, sig2_full=0.5, n_seeds=3,
                      miss=0.05):
    rng = np.random.default_rng(seed)
    lam = np.r_[np.array([6.0, 5.0, 4.0, 3.0]), 0.6 ** np.arange(40)]
    X = _descriptor_with_spectrum(C, 300, lam, rng)
    pcs, _ = P.descriptor_pcs(X)
    Z = pcs[:, :k_true] / pcs[:, :k_true].std(0)
    W = rng.standard_normal((T, k_true, G)) * np.sqrt(v_e / k_true)
    f = np.einsum('cj,tjg->ctg', Z, W)
    u = rng.standard_normal((C, T, G)) * np.sqrt(v_u)
    mu = rng.standard_normal(G) * 2
    a = rng.standard_normal((T, G)) * 1.5
    b = rng.standard_normal((C, G)) * 1.0
    truth = mu + a[None] + b[:, None] + f + u
    M = rng.random((C, T)) > miss
    ci, ti = np.nonzero(M)
    S = []
    for s in range(n_seeds):
        nz = np.random.default_rng([seed, s])
        half_sd = np.sqrt(2 * sig2_full)                     # half depth = 2 x full-depth variance
        YA = truth[ci, ti] + nz.standard_normal((len(ci), G)) * half_sd
        YB = truth[ci, ti] + nz.standard_normal((len(ci), G)) * half_sd
        keep = P.filter_min_lines(ci, ti)
        S.append(P.halves_stats(YA[keep], YB[keep], ci[keep], ti[keep], pcs, [k_true, 2 * k_true]))
    d_true = lam.sum() ** 2 / (lam ** 2).sum()
    return S, X, d_true, dict(var_i=v_e + v_u, sig2=sig2_full, R2=v_e / (v_e + v_u)), k_true


def test_estimators_recover_known_parameters():
    S, X, d_true, true, k = _synthetic_halves(1)
    d = P.participation_ratio(X)
    assert abs(d - d_true) / d_true < 1e-8
    agg = P.aggregate(S, d, k, n_boot=200)
    rel = lambda est, tr: abs(est - tr) / tr
    assert rel(agg['var_i']['estimate'], true['var_i']) < 0.10, agg['var_i']
    assert rel(agg['sig2_eps']['estimate'], true['sig2']) < 0.10, agg['sig2_eps']
    assert rel(agg[f'R2_adj_k{k}']['estimate'], true['R2']) < 0.10, agg[f'R2_adj_k{k}']
    # raw (unadjusted) R^2 is inflated by in-sample fitting of the unexplained part; adjustment must help
    assert agg[f'R2_raw_k{k}']['estimate'] > agg[f'R2_adj_k{k}']['estimate']
    ct_true = P.c_star_theory(d_true, true['R2'], true['var_i'], true['sig2'])
    assert rel(agg['C_theory']['estimate'], ct_true) < 0.10
    cm_true = P.c_star_min(d_true, true['var_i'], true['sig2'])
    assert rel(agg['C_min']['estimate'], cm_true) < 0.10
    lo, hi = agg['var_i']['ci95']
    assert lo <= agg['var_i']['estimate'] <= hi


def test_noise_free_cross_product_unaffected_by_noise_level():
    """Var(i) and R^2 must not move with the noise level (they are cross-products of independent halves)."""
    a = P.aggregate(_synthetic_halves(2, sig2_full=0.1)[0], 10.0, 4, n_boot=50)
    b = P.aggregate(_synthetic_halves(2, sig2_full=2.0)[0], 10.0, 4, n_boot=50)
    assert abs(a['var_i']['estimate'] - b['var_i']['estimate']) / a['var_i']['estimate'] < 0.10
    assert abs(a['R2_adj_k4']['estimate'] - b['R2_adj_k4']['estimate']) < 0.06


def test_c_star_algebra():
    d, R2, v, e = 7.5, 0.3, 2.0, 1.0
    res = (1 - R2) * v + e
    assert math.isclose(P.c_star_theory(d, R2, v, e), d + 1 + d * res / (R2 * v))
    assert math.isclose(P.c_star_theory(d, 1.0, v, e), P.c_star_min(d, v, e))   # R2 = 1 -> lower bound
    assert P.c_star_theory(d, 0.0, v, e) == math.inf
    assert P.c_star_theory(d, 0.5, v, e) > P.c_star_min(d, v, e)


# ----------------------------------------------------------------------------- A1: depth rescaling

def _simulate_counts(rng, mean_cells, n_cells):
    """Per cell Poisson counts; mean_cells (G,) expected counts per cell; plus 30 % extra 'native' library."""
    X = rng.poisson(mean_cells[None], size=(n_cells, len(mean_cells)))
    lib = X.sum(1) + rng.poisson(0.3 * mean_cells.sum(), n_cells)
    return X, lib


def test_depth_rescaling_poisson():
    rng = np.random.default_rng(5)
    C, T, G = 6, 25, 60
    base = np.exp(rng.uniform(np.log(0.05), np.log(20), G))
    eff = rng.standard_normal((C, T, G)) * 0.3
    ncell = rng.integers(40, 90, (C, T))
    nctrl = 300

    def dataset(rng):
        ctrl = [_simulate_counts(rng, base * np.exp(0.2 * c / C), nctrl) for c in range(C)]
        trt = {(c, t): _simulate_counts(rng, base * np.exp(0.2 * c / C + eff[c, t]), ncell[c, t])
               for c in range(C) for t in range(T)}
        return ctrl, trt

    def deltas(ctrl, trt, pick_c, pick_t):
        out = []
        for c in range(C):
            Xc, lc = ctrl[c]; mc = pick_c(c, len(lc))
            cl = P.log_cp10k(Xc[mc].sum(0)[None], [lc[mc].sum()])[0]
            for t in range(T):
                X, l = trt[(c, t)]; m = pick_t(c, t, len(l))
                out.append(P.log_cp10k(X[m].sum(0)[None], [l[m].sum()])[0] - cl)
        return np.array(out)

    all_ = lambda *a: slice(None)
    r1, r2 = dataset(rng), dataset(rng)
    emp = 0.5 * np.mean((deltas(*r1, lambda c, n: all_(), lambda c, t, n: all_())
                         - deltas(*r2, lambda c, n: all_(), lambda c, t, n: all_())) ** 2)
    ests = []
    for s in range(4):
        hc = lambda c, n: P.half_assignment([n], s, base=100 + c)
        ht = lambda c, t, n: P.half_assignment([n], s, base=1000 + c * T + t)
        dA = deltas(*r1, hc, ht)
        dB = deltas(*r1, lambda c, n: ~hc(c, n), lambda c, t, n: ~ht(c, t, n))
        ests.append(P.DEPTH_SCALE * 0.5 * np.mean((dA - dB) ** 2))
    est = float(np.mean(ests))
    assert abs(est - emp) / emp < 0.10, (est, emp)


# ----------------------------------------------------------------------------- reader on a synthetic shard

def _write_shard(path, cells, rng):
    import pyarrow as pa, pyarrow.parquet as pq
    genes, exprs, samples, lines = [], [], [], []
    for sample, line, toks, vals in cells:
        genes.append([1] + list(toks)); exprs.append([-2.0] + list(map(float, vals)))
        samples.append(sample); lines.append(line)
    tab = pa.table({'genes': pa.array(genes, pa.list_(pa.int64())),
                    'expressions': pa.array(exprs, pa.list_(pa.float32())),
                    'sample': samples, 'cell_line_id': lines, 'plate': ['plate3'] * len(cells)})
    pq.write_table(tab, path, row_group_size=7)


def test_half_accumulator_on_synthetic_shard(tmp_path):
    TE = P._tahoe_helpers()[1]
    rng = np.random.default_rng(3)
    groups = [('s1', 'L1'), ('s2', 'L1'), ('dmso', 'L1')]
    ncell = np.array([11, 26, 40])
    cells = []
    for g, (smp, line) in enumerate(groups):
        for _ in range(ncell[g]):
            toks = rng.choice(np.arange(2, 14), 6, replace=False)
            cells.append((smp, line, toks, rng.integers(1, 9, 6)))
    cells.append(('other', 'L9', [3], [5]))                        # not in any group -> ignored
    order = rng.permutation(len(cells)); cells = [cells[i] for i in order]
    _write_shard(tmp_path / 'x.parquet', cells, rng)
    lookup = np.full(14, -1); lookup[2:12] = np.arange(10)          # tokens 12, 13 outside the panel
    acc = P.HalfAccumulator(ncell, lookup, 10, seeds=[0, 1])
    k2g = {k: i for i, k in enumerate(groups)}
    P.read_tahoe_shards(acc, k2g, [dict(path='x.parquet', expected_plate='plate3', expected_rows=len(cells))],
                        tmp_path, TE, log=lambda *_: None)
    full = np.zeros((3, 10)); lib = np.zeros(3)
    ordinal = np.zeros(3, int); A = [np.zeros((3, 10)) for _ in range(2)]
    masks = [P.half_assignment(ncell, s) for s in (0, 1)]
    off = np.r_[0, np.cumsum(ncell)]
    for smp, line, toks, vals in cells:
        g = k2g.get((smp, line))
        if g is None:
            continue
        for tk, v in zip(toks, vals):
            lib[g] += v
            if lookup[tk] >= 0:
                full[g, lookup[tk]] += v
                for s in range(2):
                    if masks[s][off[g] + ordinal[g]]:
                        A[s][g, lookup[tk]] += v
        ordinal[g] += 1
    assert np.array_equal(acc.full, full) and np.allclose(acc.lib_full, lib)
    for s in range(2):
        assert np.array_equal(acc.A[s], A[s])
        assert [masks[s][off[g]:off[g + 1]].sum() for g in range(3)] == list(ncell // 2)
    # half responses: F equals the direct formula; A/B are built from disjoint cells
    rt = dict(take=np.array([0, 1]), row_ctrl=np.array([0, 0]), ctrl_groups=[np.array([2])],
              ctx=np.array([0, 0]), pert=np.array([0, 1]), pert_names=['d0', 'd1'])
    ctx, drug, YA, YB, YF = P.tahoe_half_responses(acc, rt, 0, with_full=True)
    ref = P.log_cp10k(full[:2], lib[:2]) - P.log_cp10k(full[2:3], lib[2:3])
    assert np.allclose(YF, ref, atol=1e-5)
    assert not np.allclose(YA, YB)


def test_half_responses_ignore_unused_empty_control_pool():
    """Regression (2026-11-05 full run): an UNUSED control pool whose half A is empty (1-cell wells) must not
    abort the run; used pools are unaffected."""
    rng = np.random.default_rng(5)
    ncell = np.array([30, 30, 40, 1, 1])                 # rows 0, 1; used ctrl = group 2; unused ctrl = 3 + 4
    acc = P.HalfAccumulator(ncell, np.arange(6), 6, seeds=[0])
    acc.full[:] = rng.integers(1, 50, (5, 6)); acc.lib_full[:] = acc.full.sum(1) + 10
    acc.A[0][:] = np.floor(acc.full / 2); acc.A[0][3:] = 0
    acc.lib_A[0][:] = acc.A[0].sum(1) + 5; acc.lib_A[0][3:] = 0
    rt = dict(take=np.array([0, 1]), row_ctrl=np.array([0, 0]), ctrl_groups=[np.array([2]), np.array([3, 4])],
              ctx=np.array([0, 1]), pert=np.array([0, 0]), pert_names=['d0'])
    ctx, drug, YA, YB, YF = P.tahoe_half_responses(acc, rt, 0, with_full=True)
    ref = P.log_cp10k(acc.full[:2], acc.lib_full[:2]) - P.log_cp10k(acc.full[2:3], acc.lib_full[2:3])
    assert np.allclose(YF, ref, atol=1e-5) and np.isfinite(YA).all() and np.isfinite(YB).all()


# ----------------------------------------------------------------------------- real-data structure

def test_descriptor_matches_task6_runner():
    import rq2_tahoe_curve as R6      # the TEST may import it; rq2_params may not (blinding)
    X6, n6 = R6.context_descriptors(return_names=True)
    X, n = P.tahoe_descriptor()
    assert n == n6 and np.array_equal(X, X6) and X.shape == (48, 5955)


def test_folds_and_subset_match_dcp():
    pytest.importorskip('torch')
    sys.path.insert(0, str(P.SCRIPTS))
    import dcp_drug_lodo as L
    _, names = P.tahoe_descriptor()
    assert L.folds(names) == P._folds(names)
    for C in (6, 8, 12, 16, 24, 32, 40, 47):
        assert L.subset(names[:40], C) == P._subset(names[:40], C)


def test_row_table_matches_model_target():
    chk = P.check_row_table_vs_model(P.tahoe_row_table())
    assert chk['rows'] == chk['rows_model'] == 13011
    assert chk['ctx_names_equal'] and chk['pert_names_equal'] and chk['ctx_equal'] and chk['pert_equal']
    assert chk['is_repeat_equal'] and chk['weight_max_rel_diff'] < 1e-5
    assert chk['cells_total'] == 4431233


# ----------------------------------------------------------------------------- SECONDARY sensitivities S1-S3

def _synthetic_with_b(seed=11, C=48, T=120, G=80, k_true=4, v_b=1.0, R2b=0.5, v_e=0.6, v_u=0.4, sig2_full=0.5,
                      n_seeds=2, miss=0.05):
    rng = np.random.default_rng(seed)
    lam = np.r_[np.array([6.0, 5.0, 4.0, 3.0]), 0.6 ** np.arange(40)]
    X = _descriptor_with_spectrum(C, 300, lam, rng)
    pcs, _ = P.descriptor_pcs(X)
    Z = pcs[:, :k_true] / pcs[:, :k_true].std(0)
    f = np.einsum('cj,tjg->ctg', Z, rng.standard_normal((T, k_true, G)) * np.sqrt(v_e / k_true))
    u = rng.standard_normal((C, T, G)) * np.sqrt(v_u)
    b = Z @ (rng.standard_normal((k_true, G)) * np.sqrt(R2b * v_b / k_true)) \
        + rng.standard_normal((C, G)) * np.sqrt((1 - R2b) * v_b)
    truth = rng.standard_normal(G) + rng.standard_normal((T, G))[None] + b[:, None] + f + u
    M = rng.random((C, T)) > miss
    ci, ti = np.nonzero(M)
    S, EX = [], []
    for s in range(n_seeds):
        nz = np.random.default_rng([seed, s])
        sd = np.sqrt(2 * sig2_full)
        YA = truth[ci, ti] + nz.standard_normal((len(ci), G)) * sd
        YB = truth[ci, ti] + nz.standard_normal((len(ci), G)) * sd
        keep = P.filter_min_lines(ci, ti)
        S.append(P.halves_stats(YA[keep], YB[keep], ci[keep], ti[keep], pcs, [k_true]))
        EX.append(P.halves_extra_stats(YA[keep], YB[keep], ci[keep], ti[keep]))
    nbar = len(ci) / C
    return S, P.sum_extra_stats(EX), X, k_true, dict(v_b=v_b, R2b=R2b, sig2=sig2_full, nbar=nbar)


def test_sensitivity_unit_weights_reproduce_primary():
    S, E, X, k, _ = _synthetic_with_b()
    d = P.participation_ratio(X)
    agg = P.aggregate(S, d, k, n_boot=20)
    pt = P.sensitivity_estimates(E, X, k=k)
    assert math.isclose(pt['d'], d, rel_tol=1e-10)
    for key, ref in [('var_i', 'var_i'), ('sig2_eps', 'sig2_eps'), ('R2_adj_i', f'R2_adj_k{k}'),
                     ('C_theory', 'C_theory'), ('C_min', 'C_min')]:
        assert math.isclose(pt[key], agg[ref]['estimate'], rel_tol=1e-8), key


def test_sensitivity_recovers_line_component():
    _, E, X, k, tr = _synthetic_with_b()
    pt = P.sensitivity_estimates(E, X, k=k)
    assert abs(pt['var_b'] - tr['v_b']) / tr['v_b'] < 0.25          # only 48 lines: sampling spread of b
    assert abs(pt['R2_adj_b'] - tr['R2b']) < 0.2
    # b noise: empirical split-half estimate agrees with the sig2_eps / nbar model
    assert abs(pt['noise_b_empirical'] - tr['sig2'] / tr['nbar']) / (tr['sig2'] / tr['nbar']) < 0.15
    assert abs(pt['noise_b_model'] - tr['sig2'] / tr['nbar']) / (tr['sig2'] / tr['nbar']) < 0.10
    # S1 / S2 algebra
    d, vb, vi, nb, e = pt['d'], pt['var_b'], pt['var_i'], pt['noise_b_model'], pt['sig2_eps']
    rb, ri = pt['R2_adj_b'], pt['R2_adj_i']
    assert math.isclose(pt['S1_C_min'], d + 1 + d * (nb + e) / (vb + vi))
    assert math.isclose(pt['S1_C_theory'],
                        d + 1 + d * ((1 - rb) * vb + nb + (1 - ri) * vi + e) / (rb * vb + ri * vi))
    assert math.isclose(pt['S2_C_b_theory'], P.c_star_theory(d, rb, vb, nb))
    ne = pt['noise_b_empirical']
    assert math.isclose(pt['S1emp_C_min'], d + 1 + d * (ne + e) / (vb + vi))
    assert math.isclose(pt['S2emp_C_b_theory'], P.c_star_theory(d, rb, vb, ne))
    # with b (cheap, low-noise signal) added, the minimum context count can only fall
    assert pt['S1_C_min'] < pt['C_min']


def test_weighted_descriptor_equals_duplicated_rows():
    rng = np.random.default_rng(4)
    X = _descriptor_with_spectrum(20, 60, 0.7 ** np.arange(12), rng)
    w = np.bincount(rng.integers(0, 20, 20), minlength=20).astype(float)
    _, d_w = P.weighted_descriptor(X, w)
    assert math.isclose(d_w, P.participation_ratio(np.repeat(X, w.astype(int), axis=0)), rel_tol=1e-9)
    _, d1 = P.weighted_descriptor(X, np.ones(20))
    assert math.isclose(d1, P.participation_ratio(X), rel_tol=1e-10)


def test_wterms_weighted_equals_duplicated():
    rng = np.random.default_rng(6)
    n, G = 9, 30
    a, b = rng.standard_normal((n, G)), rng.standard_normal((n, G))
    Z = rng.standard_normal((n, 2))
    w = np.array([1, 2, 0, 1, 3, 1, 1, 2, 1], float)
    m = w > 0
    tot, ex = P._wterms((a @ b.T)[np.ix_(m, m)], w[m], Z[m])
    idx = np.repeat(np.arange(n), w.astype(int))
    A, B, Zd = a[idx], b[idx], Z[idx]
    Xd = np.column_stack([np.ones(len(idx)), Zd])
    H = Xd @ np.linalg.pinv(Xd); J = np.full((len(idx),) * 2, 1 / len(idx))
    assert math.isclose(tot, float(((A - A.mean(0)) * B).sum()), rel_tol=1e-9)
    assert math.isclose(ex, float(((H - J) @ A * B).sum()), rel_tol=1e-9)


def test_line_bootstrap_runs_and_covers_point():
    _, E, X, k, _ = _synthetic_with_b(n_seeds=1)
    pt = P.sensitivity_estimates(E, X)
    bs = P.line_bootstrap(E, X, n_boot=60, seed=1)
    for key in ('d', 'var_i', 'var_b', 'C_min', 'S1_C_min'):
        lo, hi = bs[key]['ci95']
        assert lo <= pt[key] <= hi, key
    assert bs['d']['ci95'][1] > bs['d']['ci95'][0]                   # d actually varies across draws


def test_line_jackknife_and_cache_roundtrip(tmp_path):
    _, E, X, k, _ = _synthetic_with_b(n_seeds=1)
    jk = P.line_jackknife(E, X)
    pt = P.sensitivity_estimates(E, X)
    for key in ('d', 'var_i', 'var_b', 'C_min', 'S1_C_min'):
        lo, hi = jk[key]['ci95']
        assert lo < pt[key] < hi and jk[key]['se'] > 0, key
    P.save_extra_stats(E, tmp_path / 'e.npz')
    E2 = P.load_extra_stats(tmp_path / 'e.npz')
    pt2 = P.sensitivity_estimates(E2, X)
    assert all(math.isclose(pt[key], pt2[key], rel_tol=1e-12) for key in ('C_theory', 'S1_C_theory', 'var_b'))
