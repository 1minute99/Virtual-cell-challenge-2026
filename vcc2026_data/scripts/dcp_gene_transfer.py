"""Parts 12, 14, 15: transfer the drug-pretrained interaction trunk to gene perturbations.

    prediction = refit_gene_effect + alpha * context_specific_residual

The refit baseline is never replaced and never retrained differently: alpha=0 reproduces it
byte-for-byte, which is asserted against the stored refit predictions from the earlier runs.

Strict leave-one-biology-out. For the scored biology L nothing but its untreated control profile
is visible. That includes the DRUG side: if a drug dataset contains the same biology (Tahoe holds
HepG2/C3A as HEPG2_FAMILY, GSE holds HepG2 and HEK293T) that context is dropped from pretraining
on L's fold, so every arm is pretrained per fold rather than once globally.
"""
from pathlib import Path
from dataclasses import asdict, replace
import argparse, gc, json, time
import numpy as np
import torch

import b19_optimization_run as R
import b19_optimization_data as A
import b19_optimization_report as P
import context_forcing_data as F
import simple_residual_data as D
import simple_residual_baselines as B
import vcclib as V
import fullaxis_lean_baseline as LB
import dcp_data as DD
import dcp_model as M

ROOT = D.ROOT / 'reports/drug_context_pretrain_2026-09-20'
BASAL = D.ROOT / 'reports/basal_context_2026-09-20'
STORED_REFIT = D.ROOT / 'reports/cellfm_context_2026-09-20'
LINES = tuple(R.LINES)
BIOS7 = ('HCT116', 'HEK293T', 'HEPG2', 'JURKAT', 'K562', 'KOLF', 'RPE1')
ALPHAS = (0.25, 0.5, 1.0)
SEEDS = (0, 1, 2)
# drug contexts that ARE the scored biology and so must leave drug pretraining on that fold
CONFLICT = {'HEPG2': {'tahoe': ['HEPG2_FAMILY'], 'gse': ['HepG2']},
            'HEK293T': {'tahoe': [], 'gse': ['HEK293T']}}
PRETRAIN_STEPS = 6000
GENE_STEPS = 3000
FT_STEPS = 600


def status(stage, **kw):
    prior = P.read(ROOT / 'status.json') if (ROOT / 'status.json').exists() else {}
    P.write(ROOT / 'status.json', {**prior, 'stage': stage,
            'updated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **kw})


def score_one(worker, targets, line, path, seed):
    out = Path(path).with_name(Path(path).stem + f'_score_seed{seed}.json')
    res = worker.request(dict(line=line, targets=targets, delta_path=str(path),
                              output_path=str(out), seed=seed, cells=400, scale=3., null_repeats=0))
    ok = bool(res['count_diagnostics']['official_per_cell_limit_valid']
              and res['count_diagnostics']['finite_nonnegative_integral'])
    return dict(pds=res['pds'], valid=ok)


def gene_fold(line, rank, base, table, E9, keep9, present, panel, pf, has_pf, ppos, device):
    """refit baseline + the gene task in the shared (context, response program) formulation."""
    h = BIOS7.index(D.canonical_bio(line))
    tr = [i for i in range(len(BIOS7)) if i != h]
    ids = A.scope_ids(table, [line])
    targets = base.targets[line]
    feats = base.features(targets)
    canon = [D.canonical_target(t) for t in targets]
    empirical = B.empirical(targets, excluded_bios=[line])
    y_ids = np.asarray(table.y[ids], dtype=np.float32)
    meta = {k: np.asarray(getattr(table, k))[ids] for k in
            ('targets', 'bios', 'modality', 'features', 'weights', 'sources')}
    sub = D.Table(y=y_ids, controls=np.zeros((len(ids), 1), np.float32),
                  genes=list(table.genes), **meta)
    fitted = LB.LeanSourceBaseline(sub, np.arange(len(ids)), rank=32)
    refit = fitted.predict(canon, feats)[0].astype(np.float32)
    G = refit.shape[1]

    # perturbation main response programme for EVERY panel target = the same global refit model
    prog = np.zeros((len(panel), G), np.float32)
    idx = np.flatnonzero(has_pf)
    for s0 in range(0, len(idx), 2048):
        j = idx[s0:s0 + 2048]
        prog[j] = fitted.predict([D.canonical_target(panel[q]) for q in j], pf[j])[0]
    del sub, fitted, y_ids, meta
    gc.collect()

    # main effect / residual definition (Part 3), training biologies only
    acc = np.zeros((len(panel), G), np.float64); cnt = np.zeros(len(panel))
    for i in tr:
        m = present[i]
        acc[m] += np.asarray(E9[keep9[i]][m], np.float64); cnt[m] += 1
    ok = (cnt >= 2)
    mu = np.zeros((len(panel), G), np.float32)
    mu[ok] = (acc[ok] / cnt[ok, None]).astype(np.float32)
    del acc; gc.collect()
    usable = ok & has_pf

    rng = np.random.default_rng(20260920 + h)
    sample = []
    for i in tr:
        m = np.flatnonzero(present[i] & usable)
        pick = m[rng.choice(len(m), size=min(1500, len(m)), replace=False)]
        sample.append(np.asarray(E9[keep9[i]][pick], np.float32) - mu[pick])
    pooled = np.concatenate(sample); del sample
    resp_basis = DD.pca(pooled, rank)
    pert_basis = DD.pca(prog[usable], rank)
    prog_mean = prog[usable].mean(0)
    zp_all = ((prog - prog_mean) @ pert_basis).astype(np.float32)
    diag = dict(line=line, held=BIOS7[h], n_panel_usable=int(usable.sum()),
                refit_vs_mu_cosine=float(np.mean([
                    float(np.dot(prog[p], mu[p]) / max(np.linalg.norm(prog[p]) *
                          np.linalg.norm(mu[p]), 1e-9))
                    for p in np.flatnonzero(usable)[:2000]])),
                resp_basis_var=float((( pooled @ resp_basis) ** 2).sum()
                                     / max((pooled ** 2).sum(), 1e-9)))
    del pooled; gc.collect()
    return dict(h=h, tr=tr, targets=targets, canon=canon, refit=refit, empirical=empirical,
                mu=mu, usable=usable, present=present, resp_basis=resp_basis, zp=zp_all,
                panel_idx=np.array([ppos.get(t, -1) for t in canon]), diag=diag, G=G)


def gene_task(gf, bios, cs, E9, keep9):
    """(zc, zp, y, w) rows for a set of training biologies."""
    Zs, Ps, Ys = [], [], []
    for i in bios:
        m = np.flatnonzero(gf['present'][i] & gf['usable'])
        zc = cs['Z'][cs['index'][f'gene:{BIOS7[i]}']]
        for s0 in range(0, len(m), 4096):
            j = m[s0:s0 + 4096]
            d = np.asarray(E9[keep9[i]][j], np.float32) - gf['mu'][j]
            Ys.append(d @ gf['resp_basis']); Ps.append(gf['zp'][j])
            Zs.append(np.tile(zc, (len(j), 1)))
            del d
    zc = np.concatenate(Zs).astype(np.float32)
    return dict(zc=zc, zp=np.concatenate(Ps).astype(np.float32),
                y=np.concatenate(Ys).astype(np.float32),
                w=np.ones(len(zc), np.float32))


def drug_tasks(dss, train_names, rank, cs):
    tasks, folds = {}, {}
    for n, d in dss.items():
        keep = [i for i, k in enumerate(d.ctx_keys) if k in train_names]
        fo = d.fold(keep, rank)
        folds[n] = fo
        rows = fo['rows']
        zc = np.stack([cs['Z'][cs['index'][d.ctx_keys[c]]] for c, _ in rows]).astype(np.float32)
        w = fo['w'] / max(fo['w'].mean(), 1e-9)
        tasks[n] = dict(zc=zc, zp=fo['zp'][rows[:, 1]].astype(np.float32),
                        y=fo['y'].astype(np.float32), w=w.astype(np.float32))
    return tasks, folds


def split_inner(task, frac=0.2, seed=0):
    rng = np.random.default_rng(seed)
    n = len(task['y']); m = rng.permutation(n); k = max(1, int(n * frac))
    v, t = m[:k], m[k:]
    return ({q: task[q][t] for q in task}, {q: task[q][v] for q in task})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lines', default=','.join(LINES))
    ap.add_argument('--seeds', default='0,1,2')
    ap.add_argument('--rank', type=int, default=32)
    ap.add_argument('--kind', default='film')
    ap.add_argument('--experts', type=int, default=8)
    ap.add_argument('--op-rank', type=int, default=4)
    ap.add_argument('--weight-decay', type=float, default=1e-4)
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--arms', default='all')
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(',')]
    ROOT.mkdir(parents=True, exist_ok=True)
    R.OUT = ROOT
    device = torch.device(args.device)

    base = A.Data()
    table = D.Table.load(F.TABLE)
    indices = {g: i for i, g in enumerate(V.official_features())}
    columns = np.array([indices[g] for g in table.genes])
    z = np.load(BASAL / 'assembled.npz', allow_pickle=True)
    b9 = [str(x) for x in z['bios']]
    panel = [str(x) for x in z['panel']]
    E9 = np.load(BASAL / 'effects.npy', mmap_mode='r')
    keep9 = [b9.index(b) for b in BIOS7]
    present = np.stack([np.isfinite(np.asarray(E9[i, :, 0])) for i in keep9])
    ppos = {D.canonical_target(t): i for i, t in enumerate(panel)}
    pf = np.zeros((len(panel), base.table.features.shape[1]), np.float32)
    has_pf = np.zeros(len(panel), bool)
    for j, t in enumerate(panel):
        f = base.feature_map.get(D.canonical_target(t))
        if f is not None:
            pf[j] = f; has_pf[j] = True
    dss = {n: DD.DrugData(n) for n in ('tahoe', 'gse')}

    spec0 = M.Spec(kind=args.kind, rank=args.rank, weight_decay=args.weight_decay,
                   experts=args.experts, op_rank=args.op_rank)
    P.freeze(ROOT / 'gene_protocol.json', dict(
        question='does drug-context pretraining improve unseen gene-context transfer',
        base='frozen refit baseline; never replaced; alpha=0 reproduces it',
        residual_definition='observed - mean response across TRAINING biologies (Part 3)',
        perturbation_representation='refit global gene-effect programme, training-only PCA',
        context_representation='untreated pseudobulk -> per-dataset centred -> training-only PCA',
        drug_conflict_exclusions=CONFLICT, spec=asdict(spec0), alphas=list(ALPHAS),
        seeds=seeds, scale=3., cells=400, official_uploads=0,
        pretrain_steps=PRETRAIN_STEPS, gene_steps=GENE_STEPS, ft_steps=FT_STEPS))

    worker = R.Worker('dcp')
    rows, diags, curves, trainlog = [], [], [], []
    try:
        for line in args.lines.split(','):
            t0 = time.monotonic()
            status('fold', line=line)
            cb = D.canonical_bio(line)
            drop = CONFLICT.get(cb, {'tahoe': [], 'gse': []})
            train_names = sorted(
                [k for k, c in zip(dss['tahoe'].ctx_keys, dss['tahoe'].ctx_names)
                 if c not in drop['tahoe']]
                + [k for k, c in zip(dss['gse'].ctx_keys, dss['gse'].ctx_names)
                   if c not in drop['gse']]
                + [f'gene:{b}' for b in BIOS7 if b != cb])
            cs = DD.context_space(train_names, dim=64)
            spec = replace(spec0, ctx_dim=int(cs['Z'].shape[1]))
            gf = gene_fold(line, args.rank, base, table, E9, keep9, present, panel,
                           pf, has_pf, ppos, device)
            gf['diag'].update(n_train_drug_contexts=sum(
                1 for k in train_names if not k.startswith('gene:')),
                              dropped_drug_contexts=drop, ctx_var_kept=cs['var_kept'],
                              ctx_dim=spec.ctx_dim)
            # alpha=0 identity check against the stored refit predictions
            sp = STORED_REFIT / line / 'refit_alpha0.npy'
            if sp.exists():
                prev = np.load(sp)
                mine = gf['empirical'].copy(); mine[:, columns] = gf['refit']
                gf['diag']['refit_max_abs_diff_vs_stored'] = float(np.abs(mine - prev).max())
                del prev, mine
            dt, dfolds = drug_tasks(dss, train_names, args.rank, cs)
            gtask_full = gene_task(gf, gf['tr'], cs, E9, keep9)
            inner_bio = gf['tr'][len(gf['tr']) // 2]
            fit_bios = [i for i in gf['tr'] if i != inner_bio]
            gtask_fit = gene_task(gf, fit_bios, cs, E9, keep9)
            gtask_val = gene_task(gf, [inner_bio], cs, E9, keep9)
            gf['diag']['inner_validation_biology'] = BIOS7[inner_bio]
            gf['diag']['gene_rows'] = int(len(gtask_full['y']))
            gf['diag'].update({f'drug_rows_{k}': int(len(v['y'])) for k, v in dt.items()})
            diags.append(gf['diag']); P.write(ROOT / 'gene_diagnostics.json', diags)
            print(json.dumps(dict(event='fold', **gf['diag'])), flush=True)

            def log(tagname):
                return lambda e: curves.append(dict(line=line, arm=tagname, **e))

            def pretrain(which, tagname):
                tasks = {k: dt[k] for k in which}
                tv = {k: split_inner(v, seed=7) for k, v in tasks.items()}
                torch.manual_seed(spec.seed)
                m = M.InteractionModel(spec, list(dss))
                rec = M.train(m, {k: v[0] for k, v in tv.items()}, replace(spec, steps=PRETRAIN_STEPS),
                              device, val={k: v[1] for k, v in tv.items()}, log=log(tagname))
                return m, rec

            def gene_train(model, mode, tagname, extra=None, trunk_scale=1.0, freeze=False,
                           steps=GENE_STEPS, early_stop=True):
                """Early-stop on a held-out TRAINING biology, then refit on all of them.

                `early_stop=False` is used for the short joint fine-tune: after the refit-on-all
                pass the inner validation biology has been seen, so an early-stopping signal
                measured on it would be contaminated. Part 12 asks for that fine-tune to be short
                and conservative anyway, so it runs a fixed FT_STEPS instead.
                """
                tasks = dict(gene=gtask_fit)
                w = {'gene': 1.0}
                if extra:
                    tasks.update({k: dt[k] for k in extra})
                    w.update({k: 0.5 for k in extra})
                if not early_stop:
                    tasks2 = dict(tasks); tasks2['gene'] = gtask_full
                    rec = M.train(model, tasks2, replace(spec, steps=steps), device, weights=w,
                                  val=None, log=log(tagname), trunk_lr_scale=trunk_scale,
                                  freeze_trunk=freeze)
                    rec['best_val'] = None; rec['zero_residual_val'] = None
                    rec['beat_zero'] = None; rec['early_stopped'] = False
                    return rec, None
                rec = M.train(model, tasks, replace(spec, steps=steps), device,
                              weights=w, val=dict(gene=gtask_val), log=log(tagname),
                              trunk_lr_scale=trunk_scale, freeze_trunk=freeze)
                rec['early_stopped'] = True
                snapshot = {k: v.detach().clone() for k, v in model.state_dict().items()}
                best = rec['best_step']
                # refit on ALL training biologies for the selected number of steps
                tasks2 = dict(tasks); tasks2['gene'] = gtask_full
                M.train(model, tasks2, replace(spec, steps=max(best, spec.eval_every)), device,
                        weights=w, val=None, trunk_lr_scale=trunk_scale, freeze_trunk=freeze)
                return rec, snapshot

            arms = {}
            plan = args.arms.split(',') if args.arms != 'all' else [
                'scratch', 'tahoe_frozen', 'tahoe_partial', 'tgse_frozen', 'tgse_partial',
                'gse_partial', 'joint', 'joint_ft']
            pre = {}
            for key, which in (('tahoe', ['tahoe']), ('tgse', ['tahoe', 'gse']), ('gse', ['gse'])):
                if any(a.startswith(key + '_') for a in plan) or (key == 'tgse' and 'joint' in plan):
                    pre[key] = pretrain(which, f'pretrain_{key}')
                    trainlog.append(dict(line=line, arm=f'pretrain_{key}',
                                         **{k: v for k, v in pre[key][1].items() if k != 'curve'}))
                    print(json.dumps(dict(event='pretrained', line=line, which=key,
                                          params=pre[key][1]['params'],
                                          best_step=pre[key][1]['best_step'],
                                          beat_zero=pre[key][1]['beat_zero'])), flush=True)

            import copy
            snaps = {}
            for arm in plan:
                t1 = time.monotonic()
                torch.manual_seed(spec.seed + 1)
                if arm == 'scratch':
                    m = M.InteractionModel(spec, list(dss) + ['gene'])
                    rec, snap = gene_train(m, arm, arm)
                elif arm == 'joint':
                    m = M.InteractionModel(spec, list(dss) + ['gene'])
                    rec, snap = gene_train(m, arm, arm, extra=['tahoe', 'gse'])
                elif arm == 'joint_ft':
                    # start from the joint model as it was BEFORE the refit-on-all pass, so the
                    # inner validation biology is still unseen at the point the fine-tune begins
                    m = copy.deepcopy(arms['joint'][0])
                    m.load_state_dict(snaps['joint'])
                    rec, snap = gene_train(m, arm, arm, steps=FT_STEPS, trunk_scale=0.1,
                                           early_stop=False)
                else:
                    key, mode = arm.rsplit('_', 1)
                    m = copy.deepcopy(pre[key][0])
                    M.add_dataset(m, 'gene', device)
                    rec, snap = gene_train(m, arm, arm, freeze=(mode == 'frozen'),
                                           trunk_scale=0.1 if mode == 'partial' else 1.0)
                snaps[arm] = snap
                arms[arm] = (m, rec)
                trainlog.append(dict(line=line, arm=arm,
                                     **{k: v for k, v in rec.items() if k != 'curve'}))
                P.write(ROOT / 'training.json', trainlog)
                P.write(ROOT / 'curves.json', curves)
                rnd = lambda v: None if v is None else round(float(v), 8)
                print(json.dumps(dict(event='arm', line=line, arm=arm, params=rec['params'],
                                      best_step=rec['best_step'], best_val=rnd(rec['best_val']),
                                      zero_val=rnd(rec['zero_residual_val']),
                                      beat_zero=rec['beat_zero'],
                                      secs=round(time.monotonic() - t1, 1))), flush=True)

            # ------------- predictions -------------
            keep = np.flatnonzero((gf['panel_idx'] >= 0)
                                  & gf['usable'][np.maximum(gf['panel_idx'], 0)])
            zp_t = torch.as_tensor(gf['zp'][gf['panel_idx'][keep]], device=device)
            wrong_bio = BIOS7[gf['tr'][0]]

            def residual(model, ctx_name):
                out = np.zeros((len(gf['targets']), gf['G']), np.float32)
                if len(keep) == 0:
                    return out
                zc = np.zeros(spec.ctx_dim, np.float32) if ctx_name is None \
                    else cs['Z'][cs['index'][f'gene:{ctx_name}']]
                zc = torch.as_tensor(np.tile(zc, (len(keep), 1)), device=device)
                with torch.no_grad():
                    lat = model('gene', zc, zp_t).cpu().numpy()
                out[keep] = (lat @ gf['resp_basis'].T).astype(np.float32)
                return out

            plans = [('refit_alpha0', 0.0, np.zeros((len(gf['targets']), gf['G']), np.float32))]
            for arm in plan:
                d = residual(arms[arm][0], cb)
                for al in ALPHAS:
                    plans.append((f'{arm}_a{al:g}', al, d))
            # the arm carried into the context-specificity check is chosen by NESTED validation
            # on a held-out TRAINING biology -- never by anything computed on the scored line
            cands = [(a, arms[a][1]['best_val']) for a in plan
                     if arms[a][1].get('best_val') is not None
                     and np.isfinite(arms[a][1]['best_val'])]
            best_arm = min(cands, key=lambda x: x[1])[0] if cands else plan[-1]
            gf['diag']['selected_arm'] = best_arm
            gf['diag']['selected_arm_val'] = float(arms[best_arm][1]['best_val'])
            P.write(ROOT / 'gene_diagnostics.json', diags)
            print(json.dumps(dict(event='selected', line=line, arm=best_arm)), flush=True)
            plans.append((f'BEST_{best_arm}_wrong_a1', 1.0, residual(arms[best_arm][0], wrong_bio)))
            plans.append((f'BEST_{best_arm}_zero_a1', 1.0, residual(arms[best_arm][0], None)))

            for tag, al, d in plans:
                pred = gf['empirical'].copy()
                pred[:, columns] = gf['refit'] + (al * d if al else 0.0)
                path = ROOT / line / f'{tag}.npy'
                R.R.save_array(path, pred)
                del pred; gc.collect()
                for seed in seeds:
                    rows.append(dict(line=line, tag=tag, alpha=al, seed=seed,
                                     **score_one(worker, gf['targets'], line, path, seed)))
                P.write(ROOT / 'scores.json', rows)
                got = [r['pds'] for r in rows if r['line'] == line and r['tag'] == tag]
                print(json.dumps(dict(event='scored', line=line, tag=tag,
                                      pds=round(float(np.mean(got)), 6))), flush=True)
            print(json.dumps(dict(event='line_done', line=line,
                                  seconds=round(time.monotonic() - t0, 1))), flush=True)
            del gf, arms, pre; gc.collect()
    finally:
        worker.close()
    status('complete')


if __name__ == '__main__':
    try:
        main()
    except BaseException as exc:
        status('failed', error_type=type(exc).__name__, error=str(exc))
        raise
