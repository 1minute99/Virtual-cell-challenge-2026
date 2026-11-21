"""Task 6 tests for rq2_tahoe_curve.py.

Reference: vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling_table.json, row
"47 requested / 40 trained": mean over 6 folds of the fold-pooled residual-MSE gain
+1.8985e-05 and total-response cosine gain +0.001936.

The full C = 40, seed 0, 6-fold run is too slow for an interactive test on the shared machine, so
the reproduction test reads the cache written BY THE WRAPPER ITSELF in one detached run:

    scripts/run_detached.sh rq2-t6-repro-c40s0 output/data-analysis/results/logs \
        ~/miniconda3/envs/ct_transfer/bin/python output/data-analysis/rq2_tahoe_curve.py \
        --C 40 --seed 0 --fold 0,1,2,3,4,5 \
        --out output/data-analysis/results/rq2_tahoe_curve_test_C40_seed0.parquet

DESCRIPTOR (Task 8 input; full definition in rq2_tahoe_curve's module docstring):
x_c = log1p(CP10K) pseudobulk of all untreated DMSO_TF control cells of Tahoe line c (pooled over
plates 3/6/14), on the 5,955-gene shared axis (library = sum over that axis), centred on the 48
Tahoe lines' mean -> context_descriptors() is 48 x 5,955. The network sees, per fold, the
projection of x_c onto the top min(64, n_train - 1) PCs of the training contexts (39 at C = 40),
divided by the training RMS.
"""
from pathlib import Path
import json, sys

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import rq2_tahoe_curve as T  # noqa: E402

REF = T.ROOT / 'vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling_table.json'
CACHE = HERE / 'results/rq2_tahoe_curve_test_C40_seed0.parquet'
INIT = HERE / 'results/rq2_tahoe_fold0_init_check.json'


def _ref_row():
    return next(r for r in json.load(open(REF)) if r['n_contexts_requested'] == 47)


def test_descriptor_shape_and_centring():
    X, names = T.context_descriptors(return_names=True)
    assert X.shape == (48, 5955)
    assert names == sorted(names) and all(n.startswith('tahoe:') for n in names)
    assert np.abs(X.mean(0)).max() < 1e-5          # per-dataset centring over the 48 lines


def test_folds_and_seed0_draw():
    keys = sorted(T.context_descriptors(return_names=True)[1])
    held_all = []
    for f in range(6):
        held, pool = T.fold_split(keys, f)
        assert len(held) == 8 and len(pool) == 40
        held_all += held
        assert T.draw_subset(pool, 40, 0, f) == pool          # "47 requested" == all 40
        assert T.draw_subset(pool, 47, 0, f) == pool
        a = T.draw_subset(pool, 12, 1, f)
        assert len(a) == 12 and set(a) <= set(pool)
        assert set(T.draw_subset(pool, 6, 1, f)) <= set(a)    # nested across C within a seed
        assert a == T.draw_subset(pool, 12, 1, f)             # deterministic
    assert sorted(held_all) == keys


@pytest.mark.skipif(not (CACHE.exists() and INIT.exists()), reason='run the C=40 seed-0 job + init check')
def test_fold0_reference_inside_init_distribution():
    """REPLACES the plan's 10 % test (logged deviation, 2026-11-05): the original fold-0 init is
    unrecoverable, so the acceptance test is folds 1-5 exact per fold (next test) AND the original
    fold-0 MSE gain inside the 5-95 % range of the fold-0 initialisation distribution
    (results/rq2_tahoe_fold0_init_check.json, 20 inits). The primary metric (definition 1) is the
    MSE gain; the cosine gain is printed as a secondary diagnostic only."""
    j = json.load(open(INIT))
    ref0 = j['fold0_reference']
    m = np.array([r['mse_gain'] for r in j['runs']]); c = np.array([r['cos_gain'] for r in j['runs']])
    lo, hi = np.percentile(m, [5, 95])
    clo, chi = np.percentile(c, [5, 95])
    print(f'fold-0 ref mse_gain {ref0["mse_gain"]:.3e} in [{lo:.3e}, {hi:.3e}]; '
          f'cos_gain {ref0["cos_gain"]:.5f} vs [{clo:.5f}, {chi:.5f}] (secondary, not gated)')
    assert len(m) >= 20
    assert lo <= ref0['mse_gain'] <= hi
    # the reference fold-0 value used here is the scaling/k47 fold-0 value
    k47 = {r['fold']: r for r in json.load(open(REF.parent / 'scaling/k47/results.json'))}
    assert abs(ref0['mse_gain'] - (k47[0]['no_context_resid_mse'] - k47[0]['correct_resid_mse'])) < 1e-12


@pytest.mark.skipif(not CACHE.exists(), reason='run the detached C=40 seed-0 job first')
def test_folds_1_to_5_reproduce_exactly():
    """Folds 1-5 have a recoverable initialisation (torch.manual_seed(20260920) state); each
    fold must match scaling/k47/results.json "exactly", i.e. to float precision (CPU vs the original
    device): fold-pooled residual MSE with context to 1e-5 relative (measured <= 2.7e-6) and,
    because the gain is a small difference of two MSEs, the gains to 1e-3 relative (measured
    <= 3.7e-4).
    Fold 0's original initialisation came from torch's per-process random default seed and is
    not recoverable (see rq2_tahoe_curve docstring and results/rq2_tahoe_fold0_init_check.json)."""
    ref = {r['fold']: r for r in json.load(open(REF.parent / 'scaling/k47/results.json'))}
    df = pd.read_parquet(CACHE)
    fp = T.fold_pooled(df[(df.C_requested == 40) & (df.seed == 0)]).set_index('fold')
    for f in range(1, 6):
        rm = ref[f]['no_context_resid_mse'] - ref[f]['correct_resid_mse']
        rc = ref[f]['correct_total_cos'] - ref[f]['no_context_total_cos']
        assert abs(fp.loc[f, 'mse_gain'] - rm) <= 1e-3 * abs(rm), (f, fp.loc[f, 'mse_gain'], rm)
        assert abs(fp.loc[f, 'cos_gain'] - rc) <= 1e-3 * abs(rc), (f, fp.loc[f, 'cos_gain'], rc)
        assert abs(fp.loc[f, 'resid_mse_nocontext'] - ref[f]['no_context_resid_mse']) < 1e-9
        assert abs(fp.loc[f, 'resid_mse_context'] / ref[f]['correct_resid_mse'] - 1) <= 1e-5


# ------------------------------------------------------------------ arms / selection (rework)
# One shared DrugData load; four short training runs at C = 12, seed 1, fold 0 (~2 s each).

@pytest.fixture(scope='module')
def tahoe():
    import torch
    torch.set_num_threads(4)
    import dcp_data as DD
    return DD.DrugData(T.DATASET)


@pytest.fixture(scope='module')
def both_rows(tahoe):
    return pd.DataFrame(T.run_one(tahoe, 12, 1, 0, arm='both'))


def test_arm_both_context_equals_arm_context(tahoe, both_rows):
    ctx = pd.DataFrame(T.run_one(tahoe, 12, 1, 0, arm='context'))
    assert list(ctx.context) == list(both_rows.context) and len(ctx) == 8
    np.testing.assert_allclose(both_rows.resid_mse_context, ctx.resid_mse_context, rtol=0, atol=1e-9)
    np.testing.assert_allclose(both_rows.total_cos_context, ctx.total_cos_context, rtol=0, atol=1e-9)
    assert (both_rows.best_step == ctx.best_step).all()
    print('max |diff| resid_mse_context:',
          float(np.abs(both_rows.resid_mse_context - ctx.resid_mse_context).max()))
    # derived columns
    r = both_rows
    np.testing.assert_allclose(r.ctx_pathway_gain, r.resid_mse_ablate - r.resid_mse_context)
    assert (r.resid_mse_zero == r.resid_mse_nocontext).all()
    for k in ('beat_zero', 'best_val', 'zero_residual_val', 'best_step', 'beat_zero_ablate',
              'best_val_ablate', 'zero_residual_val_ablate', 'best_step_ablate'):
        assert k in r and r[k].notna().all(), k
    assert (r.selection == 'frozen').all() and (r.arm == 'both').all()
    assert (r.ablate_z_maxabs < 1e-5).all()          # training-mean z is 0 after centring
    assert (r.params == r.params_ablate).all()       # identical network


def test_ablate_ignores_test_descriptors(tahoe, both_rows, monkeypatch):
    """Permute the held-out contexts' z rows (cyclic shift, every row moves): the context arm
    must change, the ablate arm must not move at all."""
    import dcp_drug_lodo as L
    held, _ = T.fold_split(tahoe.ctx_keys, 0)
    orig = L.build

    def build_perm(dss, train_names, rank, min_ctx=3):
        cs, fo = orig(dss, train_names, rank, min_ctx)
        ix = [cs['index'][k] for k in held]
        Z = cs['Z'].copy()
        Z[ix] = Z[np.roll(ix, 1)]
        return dict(cs, Z=Z), fo

    monkeypatch.setattr(L, 'build', build_perm)
    perm = pd.DataFrame(T.run_one(tahoe, 12, 1, 0, arm='both'))
    assert (perm.resid_mse_ablate.to_numpy() == both_rows.resid_mse_ablate.to_numpy()).all()
    assert (perm.best_step_ablate == both_rows.best_step_ablate).all()
    assert not np.allclose(perm.resid_mse_context, both_rows.resid_mse_context, rtol=0, atol=1e-12)


def test_select_include_zero_can_return_step0(tahoe, both_rows):
    """Frozen selection at C=12/seed 1/fold 0 never beats the zero residual on inner val
    (beat_zero False), so the include-zero rule must restore step 0 = the zero-residual model."""
    assert not both_rows.beat_zero.iloc[0]
    iz = pd.DataFrame(T.run_one(tahoe, 12, 1, 0, arm='context', include_zero=True))
    assert (iz.selection == 'include_zero').all()
    assert (iz.best_step == 0).all() and not iz.beat_zero.iloc[0]
    assert (iz.best_val == iz.zero_residual_val).all()
    np.testing.assert_allclose(iz.resid_mse_context, iz.resid_mse_nocontext, rtol=0, atol=1e-12)
    assert (iz.steps_run == 800).all()               # 8 non-improving evals, patience 8


def test_aggregate_is_unweighted(both_rows):
    agg = T.aggregate(both_rows)
    assert len(agg) == 1 and agg.n_contexts.iloc[0] == 8
    assert np.isclose(agg.mse_gain.iloc[0], both_rows.mse_gain.mean())
    assert np.isclose(agg.ctx_pathway_gain.iloc[0], both_rows.ctx_pathway_gain.mean())
    w = both_rows.n_pert.to_numpy(float)
    if np.ptp(w) > 0:   # weighted and unweighted differ -> confirm we did not weight
        assert not np.isclose(agg.mse_gain.iloc[0],
                              np.sum(both_rows.mse_gain * w) / w.sum(), rtol=1e-6, atol=0)
