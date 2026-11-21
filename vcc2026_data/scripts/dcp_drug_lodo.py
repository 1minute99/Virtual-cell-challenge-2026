"""Parts 8 + 9: strict leave-drug-context-out validation and the small architecture sweep.

Folds hold out COMPLETE cell lines / cell types from both drug datasets at once. A held-out
context influences nothing: not training, not early stopping, not the context PCA, not the
perturbation main effect, not the response basis, not the architecture choice.

Four arms are scored on every held-out context:
    no_context   residual = 0, i.e. the perturbation main effect alone
    correct      the held-out context's own z_context
    wrong        a REAL different held-out context's z_context
    zero_context z_context = 0
"""
from pathlib import Path
from dataclasses import asdict, replace
import argparse, itertools, json, time
import numpy as np
import torch

import dcp_data as DD
import dcp_model as M

ROOT = Path('/home/wonmin-kim/AIVC/vcc2026_data')
OUT = ROOT / 'reports/drug_context_pretrain_2026-09-20'
NFOLD = 6
INNER_VAL = 2          # training contexts reserved for early stopping, per dataset


def folds(names, k=NFOLD):
    return {n: i % k for i, n in enumerate(sorted(names))}


def build(dss, train_names, rank, min_ctx=3):
    """Per-dataset fold objects + the shared context space, fitted on `train_names` only."""
    cs = DD.context_space(train_names, dim=64)
    out = {}
    for name, d in dss.items():
        tr = [i for i, k in enumerate(d.ctx_keys) if k in train_names]
        out[name] = d.fold(tr, rank, min_ctx)
        out[name]['train_ctx'] = tr
    return cs, out


def tensors(d, fo, cs, rows):
    zc = np.stack([cs['Z'][cs['index'][d.ctx_keys[c]]] for c, _ in rows]).astype(np.float32)
    zp = fo['zp'][rows[:, 1]].astype(np.float32)
    return zc, zp


def task_arrays(d, fo, cs, keep_ctx=None):
    rows = fo['rows']
    if keep_ctx is not None:
        m = np.isin(rows[:, 0], keep_ctx)
        rows, y, w = rows[m], fo['y'][m], fo['w'][m]
    else:
        y, w = fo['y'], fo['w']
    zc, zp = tensors(d, fo, cs, rows)
    w = (w / max(w.mean(), 1e-9)).astype(np.float32)
    return dict(zc=zc, zp=zp, y=y.astype(np.float32), w=w), rows


def evaluate(model, d, fo, cs, held, device, wrong_of):
    """Full-gene-space metrics for every arm on the held-out contexts of one dataset."""
    B = fo['resp_basis']
    res = {a: dict(se=0.0, n=0, cos=[], tot=[]) for a in
           ('no_context', 'correct', 'wrong', 'zero_context')}
    per_ctx = []
    for c in held:
        ps = np.flatnonzero(fo['usable'] & d.M[c])
        if len(ps) == 0:
            continue
        truth = (d.R[c, ps] - fo['main'][ps]).astype(np.float32)
        zp = torch.as_tensor(fo['zp'][ps], device=device)
        zs = dict(
            correct=cs['Z'][cs['index'][d.ctx_keys[c]]],
            wrong=cs['Z'][cs['index'][wrong_of[d.ctx_keys[c]]]],
            zero_context=np.zeros(cs['Z'].shape[1], np.float32))
        loc = {}
        for arm in res:
            if arm == 'no_context':
                pred = np.zeros_like(truth)
            else:
                zc = torch.as_tensor(np.tile(zs[arm], (len(ps), 1)), device=device)
                with torch.no_grad():
                    lat = model(d.name, zc, zp).cpu().numpy()
                pred = (lat @ B.T).astype(np.float32)
            se = float(((pred - truth) ** 2).sum())
            res[arm]['se'] += se; res[arm]['n'] += truth.size
            if arm != 'no_context':
                num = (pred * truth).sum(1)
                den = np.linalg.norm(pred, axis=1) * np.linalg.norm(truth, axis=1)
                res[arm]['cos'] += list(num / np.maximum(den, 1e-9))
            full_t = d.R[c, ps]
            full_p = fo['main'][ps] + pred
            num = (full_p * full_t).sum(1)
            den = np.linalg.norm(full_p, axis=1) * np.linalg.norm(full_t, axis=1)
            res[arm]['tot'] += list(num / np.maximum(den, 1e-9))
            loc[arm] = dict(resid_mse=se / truth.size,
                            total_cos=float(np.mean(num / np.maximum(den, 1e-9))))
        per_ctx.append(dict(dataset=d.name, context=d.ctx_names[c], n_pert=int(len(ps)),
                            wrong_from=wrong_of[d.ctx_keys[c]], **{k: v for k, v in loc.items()}))
    out = {}
    for a, v in res.items():
        out[a] = dict(resid_mse=v['se'] / max(v['n'], 1),
                      resid_cos=float(np.mean(v['cos'])) if v['cos'] else float('nan'),
                      total_cos=float(np.mean(v['tot'])) if v['tot'] else float('nan'))
    return out, per_ctx


def specs():
    out = []
    for rank, wd in itertools.product((16, 32), (1e-4, 1e-3)):
        out.append(M.Spec(kind='film', rank=rank, weight_decay=wd))
        out.append(M.Spec(kind='operator', rank=rank, weight_decay=wd, operators=4, op_rank=4))
        out.append(M.Spec(kind='operator', rank=rank, weight_decay=wd, operators=4, op_rank=8))
        out.append(M.Spec(kind='mixture', rank=rank, weight_decay=wd, experts=4))
        out.append(M.Spec(kind='mixture', rank=rank, weight_decay=wd, experts=8))
    return out


def tag(s):
    x = f'{s.kind}_r{s.rank}_wd{s.weight_decay:g}'
    if s.kind == 'mixture':
        x += f'_e{s.experts}'
    if s.kind == 'operator':
        x += f'_K{s.operators}x{s.op_rank}'
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--out', default=str(OUT / 'drug_lodo'))
    ap.add_argument('--specs', default='')
    ap.add_argument('--datasets', default='tahoe,gse')
    ap.add_argument('--context-limit', type=int, default=0,
                    help='cap on Tahoe TRAINING contexts, for the scaling experiment')
    ap.add_argument('--label', default='')
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device)
    dss = {n: DD.DrugData(n) for n in args.datasets.split(',')}
    fm = {n: folds(d.ctx_keys) for n, d in dss.items()}
    allspecs = [s for s in specs() if not args.specs or tag(s) in args.specs.split(',')]
    rng = np.random.default_rng(20260920)

    rows, perctx, curves = [], [], []
    for fi in range(NFOLD):
        held = {n: [k for k in d.ctx_keys if fm[n][k] == fi] for n, d in dss.items()}
        train_names = sorted({k for n, d in dss.items() for k in d.ctx_keys
                              if fm[n][k] != fi})
        if args.context_limit:
            th = [c for c in train_names if c.startswith('tahoe:')]
            keep = set(subset(th, args.context_limit))
            train_names = sorted([c for c in train_names
                                  if not c.startswith('tahoe:') or c in keep])
        # wrong-context partner: a REAL different held-out context
        wrong_of = {}
        for n, hs in held.items():
            pool = [c for nn, h2 in held.items() for c in h2]
            for c in hs:
                alt = [x for x in pool if x != c]
                wrong_of[c] = alt[rng.integers(len(alt))] if alt else c
        cache = {}
        for spec in sorted(allspecs, key=lambda x: x.rank):
            t0 = time.monotonic()
            if spec.rank not in cache:
                cache.clear()          # one rank in memory at a time; the cubes are large
                cache[spec.rank] = build(dss, train_names, spec.rank)
            cs, fo = cache[spec.rank]
            spec_f = replace(spec, ctx_dim=int(cs['Z'].shape[1]))
            # nested inner validation on TRAINING contexts only
            tasks, val = {}, {}
            for n, d in dss.items():
                tr = fo[n]['train_ctx']
                inner = list(tr[-INNER_VAL:]) if len(tr) > INNER_VAL + 1 else []
                fitc = [c for c in tr if c not in inner]
                tasks[n] = task_arrays(d, fo[n], cs, fitc)[0]
                if inner:
                    val[n] = task_arrays(d, fo[n], cs, inner)[0]
            model = M.InteractionModel(spec_f, list(dss))
            rec = M.train(model, tasks, spec_f, device, val=val,
                          log=lambda e, f=fi, s=spec: curves.append(
                              dict(fold=f, spec=tag(s), **e)))
            for n, d in dss.items():
                hi = [d.ctx_keys.index(c) for c in held[n]]
                met, pc = evaluate(model, d, fo[n], cs, hi, device, wrong_of)
                rows.append(dict(fold=fi, spec=tag(spec), kind=spec.kind, rank=spec.rank,
                                 weight_decay=spec.weight_decay, experts=spec.experts,
                                 op_rank=spec.op_rank, dataset=n,
                                 n_train_contexts=len(fo[n]['train_ctx']),
                                 n_contexts_requested=args.context_limit or 0,
                                 label=args.label,
                                 n_train_contexts_all=len(train_names),
                                 params=rec['params'], best_step=rec['best_step'],
                                 best_val=rec['best_val'], beat_zero=rec['beat_zero'],
                                 basis_var_kept=fo[n]['basis_var_kept'],
                                 ctx_var_kept=cs['var_kept'], **{f'{a}_{k}': v
                                 for a, m2 in met.items() for k, v in m2.items()}))
                perctx += [dict(fold=fi, spec=tag(spec), **p) for p in pc]
            json.dump(rows, open(out / 'results.json', 'w'), default=float)
            json.dump(perctx, open(out / 'per_context.json', 'w'), default=float)
            json.dump(curves, open(out / 'curves.json', 'w'), default=float)
            r = rows[-len(dss):]
            print(json.dumps(dict(event='cfg', fold=fi, spec=tag(spec), params=r[0]['params'],
                                  secs=round(time.monotonic() - t0, 1),
                                  **{f"{x['dataset']}_d_cos": round(
                                      x['correct_resid_cos'], 4) for x in r},
                                  **{f"{x['dataset']}_cmw": round(
                                      x['correct_resid_cos'] - x['wrong_resid_cos'], 4)
                                      for x in r})), flush=True)
    print(json.dumps(dict(event='done', n=len(rows))), flush=True)


def subset(names, k):
    """Deterministic representative subset: evenly spaced over the sorted context list.

    Evenly spaced rather than random so that the 6/12/24/47 arms are nested-ish and a single
    unlucky draw cannot produce the scaling curve.
    """
    names = sorted(names)
    if k >= len(names):
        return names
    step = len(names) / k
    return [names[min(len(names) - 1, int(i * step))] for i in range(k)]


if __name__ == '__main__':
    main()
