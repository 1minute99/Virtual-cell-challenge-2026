"""RQ2 Task 6 -- Tahoe context-count curve runner (wraps the frozen dcp_* protocol; edits nothing).

One call = one or more (C, seed, fold) runs. Each run trains the selected context model on C
Tahoe training contexts drawn from the fold's 40-context training pool and scores the 8 held-out
Tahoe contexts. One output row per (C, seed, fold, held-out context) is appended to `--out`.

PROTOCOL (confirmed from `vcc2026_data/scripts/dcp_drug_lodo.py`, `dcp_data.py`, `dcp_model.py`
and the seed-0 outputs in `reports/drug_context_pretrain_2026-09-20/scaling/k*/results.json`):

* Data: `DrugData('tahoe')` -- 48 contexts x 185 drugs x 18,150 genes, response
  log1p(CP10K treated) - log1p(CP10K plate-matched DMSO), condition weights from the npz.
  The scaling experiment was run with `--datasets tahoe` (only Tahoe; `n_train_contexts_all`
  equals the Tahoe count in every `scaling/k*/results.json` row).
* Folds: `dcp_drug_lodo.folds` = context i of the SORTED 48 dataset-qualified keys goes to fold
  i % 6 -> 6 folds x 8 held-out contexts, 40-context training pool per fold.
  NOTE: `context_splits.json` in the report directory does NOT hold these folds -- it holds the
  five gene-perturbation folds (K562/RPE1/HepG2/Jurkat/KOLF). The drug folds are the
  `folds()` rule above; this wrapper reproduces them by calling that same function.
* Subset draw:
    - seed 0 = the original draw: `dcp_drug_lodo.subset(pool, C)`, deterministic and EVENLY
      SPACED over the sorted 40-context pool (not random). The original arm "47 requested"
      hits `k >= len(names)` and returns the whole pool, hence "47 requested becomes 40 trained";
      C = 40 here is the identical draw.
    - seed s >= 1: `numpy.random.default_rng([20260920, s, fold]).permutation(pool)[:C]`, sorted.
      One permutation per (seed, fold), so the subsets are NESTED across C within a seed/fold
      (same spirit as the seed-0 evenly-spaced arms, which are "nested-ish").
* Per fold, fitted on the C training contexts only (dcp_data.DrugData.fold, min_ctx=3):
  perturbation main effect (weighted mean over training contexts), usable drugs (seen in >= 3
  training contexts), perturbation PCA basis rank 32, residual-response PCA basis rank 32,
  context PCA (see DESCRIPTOR). Inner validation for early stopping = the LAST 2 training
  contexts in dataset (sorted) order; the model is fitted on the other C - 2.
* Model: `dcp_model.InteractionModel`, Spec(kind='film', rank=32, hidden=32, weight_decay=1e-3,
  lr=3e-3, steps=4000 max, batch=256, eval_every=100, patience=8, AdamW, grad-clip 1.0),
  ctx_dim = context-PCA dim (= C - 1 for C <= 65). Zero-initialised head and FiLM gamma.
  Budget is identical across C (same step cap / patience / batch); the trained parameter count
  varies only through ctx_dim (4,576 params at C = 40). Checkpoint selection starts at step 100
  (the zero-residual step-0 model is never a candidate under the frozen rule). In 22 of the 24
  seed-0 scaling runs early stopping restored step 100; k6 folds 4 and 5 restored step 200
  (verified from scaling/k{6,12,24,47}/results.json). beat_zero (inner-val) was False in all 6
  runs at k6 and k12, 4/6 at k24 and 0/6 at k47.
* Seeds and initialisation:
    - seed 0 reproduces the original single-process run as far as it is recoverable. In it the
      model was CONSTRUCTED before `dcp_model.train` called torch.manual_seed(20260920). Folds
      1-5 were therefore initialised from the state left by the previous fold's
      torch.manual_seed(20260920) (nothing in training/evaluation/build consumes CPU torch
      RNG), which this wrapper reproduces exactly with torch.manual_seed(20260920) before
      construction. Fold 0 was initialised from torch's DEFAULT generator, which in the
      installed torch 2.14 (installed 2026-09-06, before the run) is seeded NON-deterministically
      per process (checked: torch.initial_seed() differs between fresh processes). The original
      fold-0 initialisation is therefore unrecoverable; seed 0 / fold 0 here uses
      torch.manual_seed(20260920) like the other folds (deviation, logged). Spec.seed = 20260920
      drives the minibatch sampler in all folds.
    - seed s >= 1: torch.manual_seed(20260920 + s) before construction and Spec.seed =
      20260920 + s (minibatch sampler), plus the subset draw above.
* Evaluation: `dcp_drug_lodo.evaluate` on the 8 held-out contexts, full 18,150-gene space, drugs
  that are `usable` and measured in that context. Arms: no_context (residual = 0, i.e. the
  training-context main effect alone), correct (own descriptor), wrong (another held-out
  context's descriptor, partner drawn with the original rng(20260920) sequence), zero_context
  (network fed z = 0; stored as resid_mse_zctx / total_cos_zctx).
  resid_mse_zero = the analytic zero residual (identical to resid_mse_nocontext, kept under both
  names). residual MSE = mean over drugs x genes of (pred_resid - (R - main))^2;
  total cosine = mean over drugs of cos(main + pred_resid, R).
  mse_gain = resid_mse_nocontext - resid_mse_context;  cos_gain = total_cos_context - total_cos_nocontext.
  The reference `scaling_table.json` numbers are the MEAN OVER FOLDS of the fold-POOLED values
  (pooled over held-out contexts, weighted by n_pert); `fold_pooled()` below reconstructs them
  from the per-context rows exactly.

DESCRIPTOR (the exact context input the model receives; Task 8 uses it):
  Base descriptor x_c (dimension D = 5,955 genes), built in `dcp_prepare.py`, stored as
  `data/processed/drug_context_pretrain/contexts.npz['profile_adj']`:
    1. raw UMI counts of ALL untreated DMSO_TF control cells of that Tahoe cell line, summed over
       every plate of the line (plate3/plate6 primary + plate14 repeat; both DMSO wells pooled),
    2. restricted to the 5,955 genes shared by Tahoe, GSE306429 and the gene-perturbation axis,
    3. log1p(1e4 * counts / sum of counts on that 5,955-gene axis)   (log1p CP10K pseudobulk),
    4. centred on the mean of the 48 Tahoe contexts (per-dataset centring).
  `context_descriptors()` returns this 48 x 5,955 matrix (rows = sorted Tahoe contexts).
  Per fold the network input is z_c = (x_c - mean_train x) @ B / s, where B = top-k right
  singular vectors (randomized SVD, random_state 0) of the training contexts' centred x,
  k = min(64, n_train - 1) (= C - 1; 39 at C = 40), and s = RMS of the training z entries.
  z_c is a linear, rank-preserving map of x_c on the training span, so the participation ratio
  of the training descriptors is unchanged by it; the 48 x 5,955 matrix is the object for d.

ARMS AND SELECTION (rework 2026-11-05, after the Task 6 quality review; DEFAULTS = frozen protocol):
  --arm context (default)  the frozen dcp protocol above, unchanged (bit-identical rows).
  --arm ablate             constant-z ablation: identical data, split, spec, network (same
                           ctx_dim/params), init seed, sampler seed, budget and early stopping, but
                           z_c is replaced by the training-context mean of z (= 0 up to float
                           error after the context-PCA centring; stored as ablate_z_maxabs) for
                           EVERY context, train, inner-val and held-out. -> resid_mse_ablate.
  --arm both               context then ablate on ONE shared build in one process; the context
                           columns are bit-identical to an --arm context run (tested).
  ctx_pathway_gain = resid_mse_ablate - resid_mse_context   (per held-out context; > 0 = the
  context pathway helps). mse_gain_ablate = resid_mse_nocontext - resid_mse_ablate.
  Per row, per arm: best_step, best_val, zero_residual_val, beat_zero (inner-val, from training),
  steps_run, train_secs (suffix _ablate for the ablate arm); run_secs = whole run incl. build.
  --select-include-zero    EXPLORATORY, default off: step 0 (zero-residual model) is a checkpoint
                           candidate (train_include_zero; otherwise identical to dcp_model.train).
                           Column selection = 'frozen' | 'include_zero'.
  aggregate(df): UNWEIGHTED mean over held-out contexts (definition 1); fold_pooled(df) is the
  n_pert-weighted dcp_report reconstruction, used only for the reproduction test.

CLI:
  python rq2_tahoe_curve.py --C 40 --seed 0 --fold 0,1,2,3,4,5 --out results/rq2_tahoe_curve.parquet
  --C/--seed/--fold accept comma lists (one data load amortised over many runs).
  --arm context|ablate|both, --select-include-zero (see above).
  --dry prints the draws and exits without training.
  --device cpu (default; the machine's GPU is shared).

REPRODUCTION (2026-11-05, CPU, 4 threads, nice 10, machine load ~35): C = 40, seed 0, folds 0-5
  -> folds 1-5 match scaling/k47 per fold to <= 1e-5 relative (per-context resid MSE, wrong
  partners 48/48); fold 0 differs only through its lost initialisation. 6-fold mean mse_gain
  1.643e-05 (ref 1.898e-05, -13.5 %), cos_gain 0.001683 (ref 0.001936, -13.1 %). Over 20 fold-0
  init seeds the init-only SD of the 6-fold mean is 2.1e-06 (11 % of the reference) and the
  original fold-0 value sits at the 90th (MSE) / 95th (cos) percentile
  (results/rq2_tahoe_fold0_init_check.json). Acceptance test (logged deviation replacing the 10 %
  test): folds 1-5 exact per fold (resid MSE <= 2.7e-6 rel., gains <= 3.7e-4 rel.) AND the original
  fold-0 MSE gain inside the 5-95 % range of that init distribution -> passes.
TIMING: data load 12 s once per process; one (C = 40, seed, fold) run 3-7 s wall-clock on CPU
  (build + train to early stop at step 900 + evaluate); 6 folds in one process 38 s. A C = 6,
  seed 1 smoke run took 29 s under load -- budget <= 30 s/run, i.e. 126 runs <= ~65 min.
  --arm both (load ~5-11): C = 6 4.6 s, C = 12 3.9 s, C = 40 6.5 s per (C, seed, fold) run
  (context train+eval 1.7-3.1 s, ablate 1.2-2.1 s, shared build the rest) -> ~1.7x one arm.

Environment: ~/miniconda3/envs/ct_transfer (python 3, torch 2.14, numpy 2.5, sklearn 1.9),
the env that ran the dcp scripts (dcp_chain.sh). Run under nice -n 10, OMP/MKL threads 4.

Assumptions (stated, checked in tests / the assumption note):
  A1 the frozen dcp_* modules and processed npz files are unchanged since 2026-09-20;
  A2 seed-0 folds 1-5 reproduce the original init exactly; fold 0's original init is lost (see above);
  A3 held-out contexts influence nothing but their own descriptor (inherited from dcp_data);
  A4 CPU vs the original device changes floats at <= ~3e-6 relative in resid MSE (tested);
  A5 the ablate arm sees no held-out descriptor (tested by permuting held-out z: resid_mse_ablate
     unchanged bit-for-bit).
"""
from pathlib import Path
import argparse, fcntl, json, os, sys, time

ROOT = Path('/home/wonmin-kim/AIVC')
SCRIPTS = ROOT / 'vcc2026_data/scripts'
sys.path.insert(0, str(SCRIPTS))

import numpy as np

BASE_SEED = 20260920
RANK, WD, KIND = 32, 1e-3, 'film'
DATASET = 'tahoe'
PROC = ROOT / 'vcc2026_data/data/processed/drug_context_pretrain'


# ------------------------------------------------------------------ descriptor

def context_descriptors(return_names=False):
    """(48 x 5,955) Tahoe context descriptor matrix x_c (see module docstring, DESCRIPTOR)."""
    z = np.load(PROC / 'contexts.npz', allow_pickle=True)
    names = np.array([str(x) for x in z['names']])
    m = np.array([str(x) for x in z['dataset']]) == DATASET
    X = np.asarray(z['profile_adj'], np.float64)[m]
    nm = names[m]
    o = np.argsort(nm)
    X, nm = X[o], nm[o]
    return (X, nm.tolist()) if return_names else X


# ------------------------------------------------------------------ draws

def fold_split(ctx_keys, fold):
    import dcp_drug_lodo as L
    fm = L.folds(ctx_keys)
    held = [k for k in ctx_keys if fm[k] == fold]
    pool = sorted(k for k in ctx_keys if fm[k] != fold)
    return held, pool


def draw_subset(pool, C, seed, fold):
    import dcp_drug_lodo as L
    pool = sorted(pool)
    if seed == 0:
        return sorted(L.subset(pool, C))
    perm = np.random.default_rng([BASE_SEED, seed, fold]).permutation(len(pool))
    return sorted(pool[i] for i in perm[:min(C, len(pool))])


def wrong_partners(ctx_keys, fold, nfold=6):
    """Replicates dcp_drug_lodo.main's rng(20260920) wrong-context draws up to `fold`."""
    import dcp_drug_lodo as L
    fm = L.folds(ctx_keys)
    rng = np.random.default_rng(BASE_SEED)
    for fi in range(fold + 1):
        held = [k for k in ctx_keys if fm[k] == fi]
        wrong_of = {}
        for c in held:
            alt = [x for x in held if x != c]
            wrong_of[c] = alt[rng.integers(len(alt))] if alt else c
    return wrong_of


def init_seed(seed, fold):
    return BASE_SEED if seed == 0 else BASE_SEED + seed


# ------------------------------------------------------------------ training (selection variants)

def train_include_zero(model, tasks, spec, device, val):
    """EXPLORATORY (--select-include-zero). dcp_model.train with ONE change: step 0 (the
    zero-initialised head, i.e. the zero-residual model) is a checkpoint candidate, so the
    running best starts at zero_val with best_step 0 instead of +inf. Optimiser, sampler, seeds,
    eval cadence, patience counter and tolerance are those of dcp_model.train (copied, since the
    frozen module must not be edited; called only with weights=None, no freezing, lr scale 1)."""
    import torch
    import dcp_model as M
    model.to(device)
    T = {k: M._to(v, device) for k, v in tasks.items()}
    V = {k: M._to(v, device) for k, v in (val or {}).items()}
    assert V, 'include-zero selection needs an inner validation set'
    names = list(T)
    p = np.full(len(names), 1.0 / len(names))
    groups = [dict(params=[q for q in model.adapter_parameters() if q.requires_grad], lr=spec.lr)]
    tp = [q for q in model.trunk_parameters() if q.requires_grad]
    if tp:
        groups.append(dict(params=tp, lr=spec.lr))
    opt = torch.optim.AdamW(groups, lr=spec.lr, weight_decay=spec.weight_decay)
    rng = np.random.default_rng(spec.seed)
    torch.manual_seed(spec.seed)

    def evaluate(D):
        model.eval()
        tot, den = 0.0, 0.0
        with torch.no_grad():
            for n, t in D.items():
                pr = model(n, t['zc'], t['zp'])
                e = ((pr - t['y']) ** 2).mean(1) * t['w']
                tot += float(e.sum()); den += float(t['w'].sum())
        model.train()
        return tot / max(den, 1e-9)

    zero_val = evaluate(V)
    best, best_step, bad, curve = zero_val, 0, 0, []          # <- the only change
    best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    model.train()
    for step in range(1, spec.steps + 1):
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
            curve.append(dict(step=step, train=M._j(tr), val=M._j(v)))
            if v < best - 1e-9:
                best, best_step, bad = v, step, 0
                best_state = {k: q.detach().clone() for k, q in model.state_dict().items()}
            else:
                bad += 1
                if bad >= spec.patience:
                    break
    model.load_state_dict(best_state)
    return dict(best_val=M._j(best), best_step=best_step, zero_residual_val=M._j(zero_val),
                beat_zero=bool(best < zero_val), params=M.n_params(model), steps_run=step,
                curve=curve)


# ------------------------------------------------------------------ one run

ARMS = ('context', 'ablate', 'both')


def ablate_space(cs, train_names):
    """Context space for the constant-z ablation: every context (train, inner-val, held-out) gets
    the SAME vector, the mean of the training contexts' z (= 0 up to float error, because the
    context PCA is centred on exactly these contexts). Nothing of a held-out descriptor enters."""
    zbar = cs['Z'][[cs['index'][k] for k in train_names]].astype(np.float64).mean(0)
    out = dict(cs)
    out['Z'] = np.tile(zbar.astype(np.float32), (cs['Z'].shape[0], 1))
    return out, float(np.abs(zbar).max())


def _fit_eval(d, fo, cs, fitc, inner, hi, wrong_of, spec, init, device, include_zero):
    import torch
    import dcp_model as M
    import dcp_drug_lodo as L
    t0 = time.monotonic()
    tasks = {DATASET: L.task_arrays(d, fo, cs, fitc)[0]}
    val = {DATASET: L.task_arrays(d, fo, cs, inner)[0]} if inner else {}
    torch.manual_seed(init)
    model = M.InteractionModel(spec, [DATASET])
    if include_zero:
        rec = train_include_zero(model, tasks, spec, torch.device(device), val)
    else:
        rec = M.train(model, tasks, spec, torch.device(device), val=val)     # frozen protocol
    met, pc = L.evaluate(model, d, fo, cs, hi, torch.device(device), wrong_of)
    return rec, met, pc, time.monotonic() - t0


def run_one(d, C, seed, fold, device='cpu', arm='context', include_zero=False):
    """arm='context' (default) is the frozen protocol. arm='ablate' trains the identical network
    (same data, split, spec, init seed, sampler seed, budget, early stopping) with z replaced by
    the constant training mean for every context. arm='both' runs both on one shared build."""
    from dataclasses import replace
    import dcp_model as M
    import dcp_drug_lodo as L
    assert arm in ARMS
    t0 = time.monotonic()
    held, pool = fold_split(d.ctx_keys, fold)
    train_names = draw_subset(pool, C, seed, fold)
    wrong_of = wrong_partners(d.ctx_keys, fold)
    cs, fo = L.build({DATASET: d}, train_names, RANK)
    fo = fo[DATASET]
    spec = M.Spec(kind=KIND, rank=RANK, weight_decay=WD,
                  seed=BASE_SEED if seed == 0 else BASE_SEED + seed)
    spec = replace(spec, ctx_dim=int(cs['Z'].shape[1]))
    tr = fo['train_ctx']
    inner = list(tr[-L.INNER_VAL:]) if len(tr) > L.INNER_VAL + 1 else []
    fitc = [c for c in tr if c not in inner]
    hi = [d.ctx_keys.index(c) for c in held]
    rest = (fitc, inner, hi, wrong_of, spec, init_seed(seed, fold), device, include_zero)
    res = {}
    if arm in ('context', 'both'):
        res['context'] = _fit_eval(d, fo, cs, *rest)
    if arm in ('ablate', 'both'):
        cs_ab, zmax = ablate_space(cs, train_names)
        res['ablate'] = _fit_eval(d, fo, cs_ab, *rest)
    secs = time.monotonic() - t0
    any_pc = next(iter(res.values()))[2]
    nan = float('nan')
    rows = []
    for k, p0 in enumerate(any_pc):
        r = dict(C_requested=int(C), n_train=len(tr), seed=int(seed), fold=int(fold),
                 context=p0['context'], n_pert=p0['n_pert'], arm=arm,
                 selection='include_zero' if include_zero else 'frozen',
                 resid_mse_nocontext=p0['no_context']['resid_mse'],
                 resid_mse_zero=p0['no_context']['resid_mse'],
                 total_cos_nocontext=p0['no_context']['total_cos'])
        if 'context' in res:
            rec, met, pc, ts = res['context']
            p = pc[k]; assert p['context'] == p0['context']
            r.update(resid_mse_context=p['correct']['resid_mse'],
                     total_cos_context=p['correct']['total_cos'],
                     resid_mse_wrong=p['wrong']['resid_mse'], total_cos_wrong=p['wrong']['total_cos'],
                     resid_mse_zctx=p['zero_context']['resid_mse'],
                     total_cos_zctx=p['zero_context']['total_cos'],
                     wrong_from=p['wrong_from'],
                     fold_resid_cos_context=met['correct']['resid_cos'],
                     fold_resid_cos_wrong=met['wrong']['resid_cos'],
                     best_step=rec['best_step'], steps_run=rec['steps_run'],
                     best_val=rec['best_val'], zero_residual_val=rec['zero_residual_val'],
                     beat_zero=rec['beat_zero'], params=rec['params'], train_secs=ts)
        else:
            r.update(resid_mse_context=nan, total_cos_context=nan)
        if 'ablate' in res:
            rec, met, pc, ts = res['ablate']
            p = pc[k]; assert p['context'] == p0['context']
            r.update(resid_mse_ablate=p['correct']['resid_mse'],
                     total_cos_ablate=p['correct']['total_cos'],
                     best_step_ablate=rec['best_step'], steps_run_ablate=rec['steps_run'],
                     best_val_ablate=rec['best_val'],
                     zero_residual_val_ablate=rec['zero_residual_val'],
                     beat_zero_ablate=rec['beat_zero'], params_ablate=rec['params'],
                     train_secs_ablate=ts, ablate_z_maxabs=zmax)
        else:
            r.update(resid_mse_ablate=nan, total_cos_ablate=nan)
        r.update(ctx_dim=spec.ctx_dim, torch_init_seed=init_seed(seed, fold),
                 sampler_seed=spec.seed, train_contexts=','.join(train_names),
                 device=device, run_secs=secs)
        r['mse_gain'] = r['resid_mse_nocontext'] - r['resid_mse_context']
        r['cos_gain'] = r['total_cos_context'] - r['total_cos_nocontext']
        r['mse_gain_ablate'] = r['resid_mse_nocontext'] - r['resid_mse_ablate']
        r['ctx_pathway_gain'] = r['resid_mse_ablate'] - r['resid_mse_context']
        r['ctx_pathway_cos_gain'] = r['total_cos_context'] - r['total_cos_ablate']
        rows.append(r)
    return rows


def aggregate(df, by=('C_requested', 'seed', 'selection')):
    """Frozen definition 1: UNWEIGHTED mean over held-out contexts (each held-out context counts
    once; never n_pert-weighted -- contrast fold_pooled(), which reproduces the n_pert-pooled
    dcp_report numbers). Returns one row per `by` group with the mean of every gain / MSE column
    present and n_contexts. Seed-and-fold means for C*_emp are a further plain mean over seeds."""
    df = df.copy()
    if 'selection' not in df:
        df['selection'] = 'frozen'
    cols = [c for c in ('mse_gain', 'cos_gain', 'ctx_pathway_gain', 'ctx_pathway_cos_gain',
                        'mse_gain_ablate', 'resid_mse_nocontext', 'resid_mse_context',
                        'resid_mse_ablate', 'total_cos_nocontext', 'total_cos_context',
                        'total_cos_ablate') if c in df]
    g = df.groupby(list(by))
    out = g[cols].mean()            # pandas mean skips NaN (arm not run)
    out['n_contexts'] = g.size()
    return out.reset_index()


def fold_pooled(df):
    """Per (C, seed, fold) pooled metrics exactly as dcp_drug_lodo/dcp_report computed them."""
    import pandas as pd
    def agg(g):
        w = g['n_pert'].to_numpy(float)
        f = lambda c: float(np.sum(g[c].to_numpy() * w) / w.sum())
        return pd.Series(dict(n_train=g['n_train'].iloc[0],
                              resid_mse_nocontext=f('resid_mse_nocontext'),
                              resid_mse_context=f('resid_mse_context'),
                              total_cos_nocontext=f('total_cos_nocontext'),
                              total_cos_context=f('total_cos_context')))
    out = df.groupby(['C_requested', 'seed', 'fold']).apply(agg, include_groups=False).reset_index()
    out['mse_gain'] = out['resid_mse_nocontext'] - out['resid_mse_context']
    out['cos_gain'] = out['total_cos_context'] - out['total_cos_nocontext']
    return out


def append_parquet(rows, path):
    import pandas as pd
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(path) + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        new = pd.DataFrame(rows)
        if path.exists():
            new = pd.concat([pd.read_parquet(path), new], ignore_index=True)
        tmp = path.with_suffix('.tmp.parquet')
        new.to_parquet(tmp, index=False)
        os.replace(tmp, path)


def _ints(s):
    return [int(x) for x in str(s).split(',') if x != '']


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--C', required=True, help='training-context count(s), comma list')
    ap.add_argument('--seed', required=True, help='seed(s), comma list; 0 = original draw')
    ap.add_argument('--fold', default='0,1,2,3,4,5', help='fold(s) 0-5, comma list')
    ap.add_argument('--out', default=str(ROOT / 'output/data-analysis/results/rq2_tahoe_curve.parquet'))
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--arm', default='context', choices=ARMS,
                    help="context (default, frozen protocol) | ablate (constant z) | both")
    ap.add_argument('--select-include-zero', action='store_true',
                    help='EXPLORATORY: step 0 (zero-residual model) is a checkpoint candidate')
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    import torch
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '4')))
    Cs, seeds, folds_ = _ints(a.C), _ints(a.seed), _ints(a.fold)
    assert all(0 <= f < 6 for f in folds_)
    import dcp_data as DD
    if a.dry:
        z = np.load(PROC / 'tahoe.npz', allow_pickle=True)
        keys = [f'{DATASET}:{c}' for c in z['ctx_names']]
        for C in Cs:
            for s in seeds:
                for f in folds_:
                    held, pool = fold_split(keys, f)
                    tn = draw_subset(pool, C, s, f)
                    print(json.dumps(dict(event='dry', C=C, seed=s, fold=f, n_train=len(tn),
                                          inner_val=tn[-2:], held=held, train=tn,
                                          torch_init_seed=init_seed(s, f))), flush=True)
        return
    t0 = time.monotonic()
    d = DD.DrugData(DATASET)
    print(json.dumps(dict(event='loaded', secs=round(time.monotonic() - t0, 1))), flush=True)
    for C in Cs:
        for s in seeds:
            for f in folds_:
                rows = run_one(d, C, s, f, a.device, arm=a.arm, include_zero=a.select_include_zero)
                append_parquet(rows, a.out)
                r0 = rows[0]
                ev = dict(event='run', C=C, seed=s, fold=f, arm=a.arm, selection=r0['selection'],
                          n_train=r0['n_train'], secs=round(r0['run_secs'], 1))
                # unweighted mean over held-out contexts (definition 1)
                for k in ('mse_gain', 'cos_gain', 'mse_gain_ablate', 'ctx_pathway_gain'):
                    v = np.array([r[k] for r in rows], float)
                    if np.isfinite(v).all():
                        ev[f'fold_{k}'] = float(v.mean())
                for k in ('best_step', 'beat_zero', 'train_secs', 'best_step_ablate',
                          'beat_zero_ablate', 'train_secs_ablate'):
                    if k in r0:
                        ev[k] = round(r0[k], 1) if 'secs' in k else r0[k]
                print(json.dumps(ev), flush=True)
    print(json.dumps(dict(event='done', secs=round(time.monotonic() - t0, 1))), flush=True)


if __name__ == '__main__':
    main()
