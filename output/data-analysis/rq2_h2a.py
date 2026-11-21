"""RQ2 Task 7 -- H2a analysis of the Tahoe context-count curve (plan zero-shot-limits-context-complexity,
Task 7, amended by A1).

Inputs (written by rq2_tahoe_curve.py, never modified here):
  results/rq2_tahoe_curve.parquet                     confirmatory: C in {6,8,12,16,24,32,40} x seeds {1,2,3}
                                                      x folds 0-5, --arm both, frozen selection
  results/rq2_tahoe_curve_explor_includezero.parquet  EXPLORATORY: same grid, --select-include-zero
  vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling/k{6,12,24,47}/per_context.json
                                                      the pre-existing seed-0 curve (context arm only),
                                                      reported separately, never pooled into the CIs.
Output: results/rq2_h2a.json, results/rq2_h2a_curve.png.
This script does NOT read rq2_params_tahoe.json / rq2_h2b.json (H2b theory side is blinded).

QUANTITIES (per held-out context row; > 0 = the context helps):
  PRIMARY (A1)   ctx_pathway_gain     = resid_mse_ablate - resid_mse_context
  REGISTERED     mse_gain             = resid_mse_nocontext (analytic zero residual) - resid_mse_context
  secondary      ctx_pathway_cos_gain = total_cos_context - total_cos_ablate
                 cos_gain             = total_cos_context - total_cos_nocontext
                 mse_gain_ablate      = resid_mse_nocontext - resid_mse_ablate   (ablate arm vs zero)
AGGREGATION: unweighted mean over the 8 held-out contexts of a run (definition 1; = rq2_tahoe_curve.aggregate),
  then plain mean over the 18 (seed, fold) runs -> curve g(C).
C*_emp: zero crossing of g(C), linear interpolation in ln C. Operationalised as the FIRST upward crossing:
  the first C_j with g(C_j) > 0 whose predecessor has g(C_{j-1}) <= 0;
  ln C* = ln C_{j-1} + (0 - g_{j-1}) (ln C_j - ln C_{j-1}) / (g_j - g_{j-1}).
  g(6) > 0 -> left-censored "<= 6" (numeric 6); never > 0 up to C = 40 -> right-censored "> 40" (numeric +inf).
  For continuous gains "> 0" and ">= 0" coincide a.s.; strict ">" only matters for the include-zero
  exploratory arm, where gains are exactly 0 when both arms select step 0. Sensitivity: LAST upward
  crossing (the crossing after which g stays > 0), reported, flagged if it changes the verdict.
CI (A1): cluster bootstrap, B = 2,000, seed 20261107: draw 18 (seed, fold) runs with replacement, then for
  each drawn run 8 of its held-out contexts with replacement; the SAME run and context draws are used for
  every C (paired, nested subsets) and for every quantity. 95 % percentile intervals (gains: linear
  interpolation; C*: inverted-CDF quantile so censored +inf draws are handled). Sensitivity: fold-cluster
  bootstrap (draw 6 folds, then 8 contexts within the fold, applied to all 3 seeds of that fold, since the
  3 seeds of a fold share the same held-out contexts).
H2a RULE (plan + A1; thresholds unchanged):
  confirmed    if the 95 % CI of g(12) is entirely < 0 AND the lower 95 % bound of C*_emp >= 24
               (censored > 40 counts as >= 24);
  rejected     if the point estimate g(12) >= 0 OR the point estimate C*_emp < 24;
  inconclusive otherwise. (Rejection is checked first; a simultaneous confirm+reject is flagged.)

ASSUMPTIONS (checked below, results in the JSON 'checks' block):
  S1 grid complete: 7 C x 3 seeds x 6 folds x 8 contexts = 1,008 rows per file, no NaN in the gain columns;
  S2 pairing: within a (seed, fold) the held-out context set is identical across C;
  S3 nesting: within a (seed, fold) the training contexts at C are a subset of those at C' > C;
  S4 ablate arm constant: ablate_z_maxabs < 1e-5;
  S5 runs as exchangeable clusters: the 3 seeds of a fold share held-out contexts, so runs are NOT
     independent -> fold-cluster bootstrap reported as sensitivity;
  S6 the bootstrap captures seed (init + subset) variance through run resampling and context variance
     through within-run resampling; it does not capture variance from the fixed fold partition.
"""
from pathlib import Path
import json, sys, time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/wonmin-kim/AIVC')
RES = HERE / 'results'
CONF = RES / 'rq2_tahoe_curve.parquet'
EXPL = RES / 'rq2_tahoe_curve_explor_includezero.parquet'
SEED0 = ROOT / 'vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling'
OUT = RES / 'rq2_h2a.json'
FIG = RES / 'rq2_h2a_curve.png'
sys.path.insert(0, str(HERE))
import rq2_tahoe_curve as T  # noqa: E402  (aggregate() cross-check only)

CS = [6, 8, 12, 16, 24, 32, 40]
SEEDS = [1, 2, 3]
FOLDS = list(range(6))
B = 2000
BOOT_SEED = 20261107
QUANT = ['ctx_pathway_gain', 'mse_gain', 'ctx_pathway_cos_gain', 'cos_gain', 'mse_gain_ablate']
LNC = np.log(np.array(CS, float))


# ------------------------------------------------------------------ crossing

def crossing(g, last=False):
    """C*_emp of a curve g over CS (see docstring). Returns (value, label)."""
    g = np.asarray(g, float)
    if g[0] > 0:
        return float(CS[0]), f'<= {CS[0]} (left-censored)'
    ups = [j for j in range(1, len(g)) if g[j] > 0 and g[j - 1] <= 0]
    if last:
        ups = [j for j in ups if (g[j:] > 0).all()]
    if not ups:
        return float('inf'), f'> {CS[-1]} (censored)'
    j = ups[-1] if last else ups[0]
    x = LNC[j - 1] + (0 - g[j - 1]) * (LNC[j] - LNC[j - 1]) / (g[j] - g[j - 1])
    c = float(np.exp(x))
    return c, f'{c:.2f}'


def n_up_crossings(g):
    g = np.asarray(g, float)
    return int(sum(1 for j in range(1, len(g)) if g[j] > 0 and g[j - 1] <= 0) + (g[0] > 0))


def _f(x):
    return None if x is None or not np.isfinite(x) else float(x)


# ------------------------------------------------------------------ data

def cube(df, col):
    """(18 runs, 8 contexts, 7 C) array; runs ordered (seed, fold); contexts sorted by name."""
    A = np.full((len(SEEDS) * len(FOLDS), 8, len(CS)), np.nan)
    ctx_names = {}
    for r, (s, f) in enumerate([(s, f) for s in SEEDS for f in FOLDS]):
        for k, C in enumerate(CS):
            g = df[(df.seed == s) & (df.fold == f) & (df.C_requested == C)].sort_values('context')
            A[r, :len(g), k] = g[col].to_numpy(float)
            ctx_names.setdefault((s, f), []).append(tuple(g['context']))
    return A, ctx_names


def checks(df, name):
    out = {}
    n_exp = len(CS) * len(SEEDS) * len(FOLDS) * 8
    out['S1_rows'] = dict(expected=n_exp, got=int(len(df)),
                          nan_in_gains=int(df[QUANT].isna().to_numpy().sum()),
                          duplicates=int(df.duplicated(['C_requested', 'seed', 'fold', 'context']).sum()))
    out['S1_pass'] = bool(len(df) == n_exp and out['S1_rows']['nan_in_gains'] == 0
                          and out['S1_rows']['duplicates'] == 0)
    pair_ok, nest_ok = True, True
    for (s, f), g in df.groupby(['seed', 'fold']):
        sets = {C: frozenset(gg['context']) for C, gg in g.groupby('C_requested')}
        if len(set(sets.values())) != 1:
            pair_ok = False
        tr = {C: set(gg['train_contexts'].iloc[0].split(',')) for C, gg in g.groupby('C_requested')}
        n = {C: int(gg['n_train'].iloc[0]) for C, gg in g.groupby('C_requested')}
        for a, b in zip(CS[:-1], CS[1:]):
            if not tr[a] <= tr[b] or n[a] != a:
                nest_ok = False
    out['S2_pairing_pass'] = pair_ok
    out['S3_nesting_pass'] = nest_ok
    out['S4_ablate_z_maxabs'] = float(df['ablate_z_maxabs'].max())
    out['S4_pass'] = bool(out['S4_ablate_z_maxabs'] < 1e-5)
    out['selection_values'] = sorted(df['selection'].unique().tolist())
    out['arm_values'] = sorted(df['arm'].unique().tolist())
    # held-out contexts shared across seeds of the same fold (S5)
    shared = all(len({frozenset(df[(df.fold == f) & (df.seed == s)]['context']) for s in SEEDS}) == 1
                 for f in FOLDS)
    out['S5_seeds_of_a_fold_share_heldout'] = bool(shared)
    return out


# ------------------------------------------------------------------ bootstrap

def boot_indices(rng):
    R = rng.integers(0, 18, size=(B, 18))
    K = rng.integers(0, 8, size=(B, 18, 8))
    return R, K


def boot_curves(A, R, K):
    """A (18, 8, 7) -> (B, 7) bootstrap curves; same R, K for every quantity (paired)."""
    sub = A[R[:, :, None], K, :]                # (B, 18, 8, 7)
    return sub.mean(axis=(1, 2))


def boot_curves_fold(A, rng):
    """Fold-cluster sensitivity: draw 6 folds, then 8 contexts per drawn fold, shared by its 3 seeds."""
    A4 = A.reshape(len(SEEDS), len(FOLDS), 8, len(CS))       # (seed, fold, ctx, C)
    Fd = rng.integers(0, 6, size=(B, 6))
    Kd = rng.integers(0, 8, size=(B, 6, 8))
    sub = A4[:, Fd[:, :, None], Kd, :]                        # (seed, B, 6, 8, C)
    return sub.mean(axis=(0, 2, 3))


def q_lin(x, p):
    return float(np.quantile(x, p))


def q_cdf(x, p):
    return float(np.quantile(x, p, method='inverted_cdf'))


def summarise(A, bc, label):
    """Point curve, per-C CI, C*_emp point + CI, censoring shares."""
    point = A.mean(axis=(0, 1))                 # mean over contexts then runs (balanced)
    cstar, clab = crossing(point)
    cstar_last, clab_last = crossing(point, last=True)
    bs = np.array([crossing(b)[0] for b in bc])
    bs_last = np.array([crossing(b, last=True)[0] for b in bc])
    curve = []
    for k, C in enumerate(CS):
        curve.append(dict(C=C, gain=float(point[k]), ci95=[q_lin(bc[:, k], .025), q_lin(bc[:, k], .975)],
                          boot_mean=float(bc[:, k].mean()), boot_sd=float(bc[:, k].std(ddof=1)),
                          share_boot_gt0=float((bc[:, k] > 0).mean())))
    lo, hi = q_cdf(bs, .025), q_cdf(bs, .975)
    return dict(
        quantity=label, curve=curve, n_up_crossings_point=n_up_crossings(point),
        cstar_emp=_f(cstar), cstar_emp_label=clab,
        cstar_emp_ci95=[_f(lo) if np.isfinite(lo) else None, _f(hi) if np.isfinite(hi) else None],
        cstar_emp_ci95_label=[('> 40' if not np.isfinite(lo) else f'{lo:.2f}'),
                              ('> 40' if not np.isfinite(hi) else f'{hi:.2f}')],
        cstar_lower95_ge24=bool(lo >= 24),
        share_boot_censored_gt40=float(np.isinf(bs).mean()),
        share_boot_leftcensored_le6=float((bs == CS[0]).mean()),
        boot_cstar_median=_f(q_cdf(bs, .5)),
        sensitivity_last_crossing=dict(cstar_emp=_f(cstar_last), label=clab_last,
                                       ci95_lower=_f(q_cdf(bs_last, .025)),
                                       share_boot_censored_gt40=float(np.isinf(bs_last).mean())),
    )


def verdict(s):
    g12 = next(c for c in s['curve'] if c['C'] == 12)
    cpt = s['cstar_emp'] if s['cstar_emp'] is not None else float('inf')
    reject = bool(g12['gain'] >= 0 or cpt < 24)
    confirm = bool(g12['ci95'][1] < 0 and s['cstar_lower95_ge24'])
    v = 'rejected' if reject else ('confirmed' if confirm else 'inconclusive')
    return dict(verdict=v, gain_C12=g12['gain'], gain_C12_ci95=g12['ci95'],
                gain_C12_ci_entirely_below_0=bool(g12['ci95'][1] < 0),
                cstar_emp_point=s['cstar_emp_label'], cstar_emp_point_lt24=bool(cpt < 24),
                cstar_lower95=s['cstar_emp_ci95_label'][0], cstar_lower95_ge24=s['cstar_lower95_ge24'],
                confirm_conditions_met=confirm, reject_conditions_met=reject,
                conflict_flag=bool(confirm and reject))


def per_seed(A, label):
    A4 = A.reshape(len(SEEDS), len(FOLDS), 8, len(CS))
    seed_curves = A4.mean(axis=(1, 2))           # (3, 7)
    k12 = CS.index(12)
    out = dict(quantity=label, seeds={})
    for i, s in enumerate(SEEDS):
        c, lab = crossing(seed_curves[i])
        out['seeds'][str(s)] = dict(curve=[float(x) for x in seed_curves[i]], cstar_emp=_f(c),
                                    cstar_emp_label=lab)
    out['gain_C12_by_seed'] = [float(x) for x in seed_curves[:, k12]]
    out['between_seed_sd_gain_C12'] = float(seed_curves[:, k12].std(ddof=1))
    fold_curves = A4.mean(axis=(0, 2))
    out['between_fold_sd_gain_C12'] = float(fold_curves[:, k12].std(ddof=1))
    run_means = A[:, :, k12].mean(1)
    out['between_run_sd_gain_C12'] = float(run_means.std(ddof=1))
    out['within_run_context_sd_gain_C12'] = float(np.sqrt(np.mean(A[:, :, k12].var(axis=1, ddof=1))))
    return out


def run_freq(df):
    """best_step / beat_zero per run (one row per run), by C and arm."""
    r = df.drop_duplicates(['C_requested', 'seed', 'fold'])
    out = []
    for C, g in r.groupby('C_requested'):
        for arm, sfx in (('context', ''), ('ablate', '_ablate')):
            bs = g['best_step' + sfx].to_numpy()
            out.append(dict(C=int(C), arm=arm, n_runs=int(len(g)),
                            beat_zero_share=float(g['beat_zero' + sfx].mean()),
                            best_step_0_share=float((bs == 0).mean()),
                            best_step_100_share=float((bs == 100).mean()),
                            best_step_median=float(np.median(bs)), best_step_max=int(bs.max()),
                            steps_run_median=float(np.median(g['steps_run' + sfx]))))
    return out


def seed0_existing():
    """Pre-existing seed-0 curve (registered quantity only; context arm), unweighted context means."""
    rows = []
    for k, C in ((6, 6), (12, 12), (24, 24), (47, 40)):
        pc = json.load(open(SEED0 / f'k{k}' / 'per_context.json'))
        df = pd.DataFrame([dict(fold=p['fold'], context=p['context'],
                                mse_gain=p['no_context']['resid_mse'] - p['correct']['resid_mse'],
                                cos_gain=p['correct']['total_cos'] - p['no_context']['total_cos'])
                           for p in pc])
        fm = df.groupby('fold')[['mse_gain', 'cos_gain']].mean()
        rows.append(dict(C=C, k_requested=k, n_contexts=int(len(df)), mse_gain=float(fm['mse_gain'].mean()),
                         cos_gain=float(fm['cos_gain'].mean()),
                         folds_mse_gain_gt0=f"{int((fm['mse_gain'] > 0).sum())}/{len(fm)}"))
    g = np.array([r['mse_gain'] for r in rows])
    lnc = np.log([r['C'] for r in rows])
    ups = [j for j in range(1, 4) if g[j] > 0 and g[j - 1] <= 0]
    if ups:
        j = ups[0]
        c = float(np.exp(lnc[j - 1] - g[j - 1] * (lnc[j] - lnc[j - 1]) / (g[j] - g[j - 1])))
        lab = f'{c:.2f}'
    else:
        c, lab = None, '> 40 (censored)'
    return dict(note='original 2026-09-20 run, seed 0, evenly spaced subsets, C in {6,12,24,40} only, '
                     'context arm only (no ablation) -> registered quantity only; unweighted mean over '
                     'held-out contexts per fold, then over 6 folds. NOT pooled into the confirmatory CIs.',
                curve=rows, cstar_emp_registered=c, cstar_emp_registered_label=lab)


def plot(conf_s, expl_s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    col = {'ctx_pathway_gain': '#2a78d6', 'mse_gain': '#eb6834', 'mse_gain_ablate': '#1baf7a'}
    lab = {'ctx_pathway_gain': 'Context pathway (ablate - context) [primary, A1]',
           'mse_gain': 'Context net vs zero residual [registered]',
           'mse_gain_ablate': 'Ablated net vs zero residual'}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, S, title in ((axes[0], conf_s, 'Confirmatory (frozen selection), seeds 1-3 x 6 folds'),
                         (axes[1], expl_s, 'Exploratory: step 0 selectable')):
        ax.axhline(0, color='#52514e', lw=1)
        for q in ('mse_gain', 'mse_gain_ablate', 'ctx_pathway_gain'):
            c = S[q]['curve']
            x = [r['C'] for r in c]
            y = np.array([r['gain'] for r in c]) * 1e4
            lo = np.array([r['ci95'][0] for r in c]) * 1e4
            hi = np.array([r['ci95'][1] for r in c]) * 1e4
            ax.fill_between(x, lo, hi, color=col[q], alpha=.15, lw=0)
            ax.plot(x, y, color=col[q], lw=2, marker='o', ms=5, label=lab[q])
        ax.set_xscale('log')
        ax.set_xticks(CS); ax.set_xticklabels([str(c) for c in CS]); ax.minorticks_off()
        ax.set_xlabel('training contexts C (log scale)')
        ax.set_title(title, fontsize=10, color='#0b0b0b')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.grid(axis='y', color='#e6e5e0', lw=.6)
    axes[0].set_ylabel('residual-MSE gain, x 1e-4 (> 0: context helps)')
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], fontsize=8.5, frameon=False, loc='lower center', ncol=3)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(FIG, dpi=150)


def main():
    t0 = time.time()
    conf = pd.read_parquet(CONF)
    expl = pd.read_parquet(EXPL)
    res = dict(task='T7 H2a', plan='input/ideas/zero-shot-limits-context-complexity-plan.md (A1)',
               bootstrap=dict(B=B, seed=BOOT_SEED, scheme='cluster: 18 (seed,fold) runs w/ repl., then 8 '
                              'held-out contexts w/ repl. within each drawn run; same draws for all C and '
                              'all quantities', ci='95 % percentile'),
               crossing_rule='first upward crossing (g_{j-1} <= 0 < g_j), linear in ln C; '
                             'g(6) > 0 -> <= 6; never > 0 -> > 40 (censored)',
               inputs=dict(confirmatory=str(CONF.relative_to(ROOT)), exploratory=str(EXPL.relative_to(ROOT))))
    res['checks'] = dict(confirmatory=checks(conf, 'conf'), exploratory=checks(expl, 'expl'))
    assert set(conf['selection']) == {'frozen'} and set(expl['selection']) == {'include_zero'}

    rng = np.random.default_rng(BOOT_SEED)
    R, K = boot_indices(rng)
    rng_f = np.random.default_rng(BOOT_SEED + 1)
    out = {}
    for tag, df in (('confirmatory', conf), ('exploratory_include_zero', expl)):
        S, PS, SF = {}, {}, {}
        cubes = {q: cube(df, q)[0] for q in QUANT}
        for q in QUANT:
            assert np.isfinite(cubes[q]).all(), (tag, q)
            S[q] = summarise(cubes[q], boot_curves(cubes[q], R, K), q)
            PS[q] = per_seed(cubes[q], q)
        # fold-cluster sensitivity (primary + registered only)
        st = rng_f.bit_generator.state
        for q in ('ctx_pathway_gain', 'mse_gain'):
            rng_f.bit_generator.state = st                # same fold draws for both quantities
            SF[q] = summarise(cubes[q], boot_curves_fold(cubes[q], rng_f), q)
            SF[q] = dict(curve_ci95={str(c['C']): c['ci95'] for c in SF[q]['curve']},
                         cstar_emp_ci95=SF[q]['cstar_emp_ci95_label'],
                         share_boot_censored_gt40=SF[q]['share_boot_censored_gt40'], verdict=verdict(SF[q]))
        # cross-check against the runner's aggregate() (unweighted)
        agg = T.aggregate(df).groupby('C_requested')[['ctx_pathway_gain', 'mse_gain']].mean()
        for q in ('ctx_pathway_gain', 'mse_gain'):
            pt = np.array([c['gain'] for c in S[q]['curve']])
            assert np.allclose(pt, agg.loc[CS, q].to_numpy(), rtol=0, atol=1e-12), q
        out[tag] = dict(summary=S, per_seed=PS, fold_cluster_sensitivity=SF, run_frequencies=run_freq(df))
    res['checks']['aggregate_crosscheck'] = 'pass (point curves equal rq2_tahoe_curve.aggregate, atol 1e-12)'

    C_ = out['confirmatory']['summary']
    res['H2a_primary_amended'] = dict(
        quantity='ctx_pathway_gain = resid_mse_ablate - resid_mse_context (A1)',
        **verdict(C_['ctx_pathway_gain']),
        last_crossing_sensitivity=C_['ctx_pathway_gain']['sensitivity_last_crossing'],
        share_boot_censored_gt40=C_['ctx_pathway_gain']['share_boot_censored_gt40'],
        fold_cluster_sensitivity_verdict=out['confirmatory']['fold_cluster_sensitivity']['ctx_pathway_gain']['verdict']['verdict'])
    res['H2a_registered_quantity'] = dict(
        quantity='mse_gain = resid_mse_nocontext (analytic zero residual) - resid_mse_context, frozen protocol',
        confound='under the frozen selection the context network can never select the zero-residual step-0 '
                 'model, while the comparator IS the zero residual; part of any negative gain is the cost of '
                 'a forced non-zero residual, not of the context pathway (A1). Reported, not used for H2a.',
        **verdict(C_['mse_gain']),
        share_boot_censored_gt40=C_['mse_gain']['share_boot_censored_gt40'])
    res['H2b_note'] = ('C*_emp for H2b = primary cstar_emp above; if censored (> 40) H2b is "untestable '
                       '(censored)" per A1. H2b verdict is NOT computed here (blinded).')
    res['results'] = out
    res['seed0_existing'] = seed0_existing()
    res['runtime_secs'] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(res, indent=1, allow_nan=False))
    plot(out['confirmatory']['summary'], out['exploratory_include_zero']['summary'])

    # ---- console summary
    for tag in ('confirmatory', 'exploratory_include_zero'):
        print(f'== {tag}')
        for q in QUANT:
            s = out[tag]['summary'][q]
            print(f'  {q:22s} ' + '  '.join(f"{c['C']}:{c['gain']:+.2e}[{c['ci95'][0]:+.1e},{c['ci95'][1]:+.1e}]"
                                            for c in s['curve']))
            print(f'  {"":22s} C*={s["cstar_emp_label"]} CI={s["cstar_emp_ci95_label"]} '
                  f'cens>40={s["share_boot_censored_gt40"]:.3f} last={s["sensitivity_last_crossing"]["label"]}')
    print('H2a primary (A1):', json.dumps({k: v for k, v in res['H2a_primary_amended'].items()
                                            if k not in ('last_crossing_sensitivity',)}))
    print('registered     :', json.dumps({k: v for k, v in res['H2a_registered_quantity'].items()
                                            if k != 'confound'}))
    for q in ('ctx_pathway_gain', 'mse_gain'):
        p = out['confirmatory']['per_seed'][q]
        print(q, 'per-seed C*:', {s: v['cstar_emp_label'] for s, v in p['seeds'].items()},
              'SD_seed(g12)=%.2e' % p['between_seed_sd_gain_C12'])
    print('seed0:', json.dumps(res['seed0_existing']['curve']), res['seed0_existing']['cstar_emp_registered_label'])
    print('checks:', json.dumps(res['checks']))
    print('runtime', res['runtime_secs'])


if __name__ == '__main__':
    main()
