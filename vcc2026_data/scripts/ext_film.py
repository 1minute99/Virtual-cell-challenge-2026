"""Part 4 of the Feng+CD4 experiment: the small FiLM context residual trained on gene data only,
with the extended set of gene contexts.

    prediction = refit_v1 + alpha * residual(context, perturbation programme)

Everything is the `scratch` arm of `dcp_gene_transfer.py` (same Spec: FiLM, rank 32, hidden 32,
wd 1e-3, lr 3e-3, batch 256, 3,000 steps, eval every 100, patience 8, early stop on ONE inner
training biology then refit on all training contexts for the selected number of steps), except:
  * gene contexts come from `ext_effects.npy` / `ext_gene_contexts.npz` (7 old biologies + 12 CD4
    donor x state + 19 HipSci lines) instead of the 7-biology basal panel;
  * the base is the refit fitted on the same extended table (`refit_v1`), so FiLM_ext vs refit_v1
    is a paired comparison of the residual alone.
Held-out line L: its context contributes only its control profile. On the KOLF fold the two
kolf-background HipSci lines are held out too; on the HepG2 fold the drug contexts of HepG2 leave
the context PCA (as in dcp). The inner early-stopping biology is the SAME one dcp used per line.

--curve trains the identical model on the first n training contexts of a fixed order (old
biologies, then CD4 donor x state, then HipSci lines) for n in 6/12/24/all, and reports residual
cosine (correct / wrong / zero context), total-response cosine gain and harness PDS (seed 0).
"""
from pathlib import Path
from dataclasses import replace
import argparse, gc, json, sys, time
import numpy as np, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dcp_gene_transfer as T
import dcp_data as DD
import dcp_model as M
import b19_optimization_run as R
import b19_optimization_data as A
import b19_optimization_report as P
import context_forcing_data as F
import simple_residual_data as D
import vcclib as V
import fullaxis_lean_baseline as LB
import ext_baselines as XB

ROOT = D.ROOT / 'reports/ext_contexts_2026-09-21'
BIOS7 = T.BIOS7
CONFLICT_GENE = {'KOLF': ['FENGLKOLF2', 'FENGLKOLF3']}
ALPHAS = (0.25, 0.5, 1.0)


def load_assets():
    z = np.load(ROOT / 'ext_gene_contexts.npz', allow_pickle=True)
    bios = [str(x) for x in z['bios']]
    return dict(bios=bios, panel=[str(x) for x in z['panel']], adj=np.asarray(z['profile_adj'], np.float64),
                E=np.load(ROOT / 'ext_effects.npy', mmap_mode='r'))


def context_space(train_gene_bios, adj_gene, gene_bios, drop_drug, dim=64):
    """Training-only PCA of the gene contexts' table-axis control profiles. Gene contexts only:
    the drug-context profiles live on another gene axis and the old gene-side profiles they were
    paired with were scrambled (CONTROLS_COMMON_BUG.md). dim = min(64, n_train - 1)."""
    keys = [f'gene:{b}' for b in gene_bios]
    idx = {k: i for i, k in enumerate(keys)}
    Xp = np.asarray(adj_gene, np.float64)
    tr = np.array([idx[f'gene:{b}'] for b in train_gene_bios])
    cm = Xp[tr].mean(0)
    Bm = DD.pca(Xp[tr] - cm, min(dim, len(tr) - 1))
    Z = ((Xp - cm) @ Bm).astype(np.float32)
    s = float(np.sqrt(np.mean(Z[tr] ** 2))); Z = Z / max(s, 1e-9)
    return dict(Z=Z, index=idx, dim=int(Bm.shape[1]), n_train=len(tr), n_drug=0)


def fold(line, S, X, base, table, rank, device):
    cb = D.canonical_bio(line)
    held = {cb, *CONFLICT_GENE.get(cb, [])}
    bios = X['bios']
    tr_all = [i for i, b in enumerate(bios) if b not in held]
    h = bios.index(cb)
    ids = F.permitted_ids(table, [line])
    tb = np.asarray(table.bios)[ids].astype(str)
    ids = ids[~np.isin(tb, list(held))]
    targets = base.targets[line]; feats = base.features(targets)
    canon = [D.canonical_target(t) for t in targets]
    y = np.asarray(table.y[ids], np.float32)
    meta = {k: np.asarray(getattr(table, k))[ids] for k in ('targets', 'bios', 'modality', 'features', 'weights', 'sources')}
    sub = D.Table(y=y, controls=np.zeros((len(ids), 1), np.float32), genes=list(table.genes), **meta)
    fitted = LB.LeanSourceBaseline(sub, np.arange(len(ids)), rank=32)
    refit = fitted.predict(canon, feats)[0].astype(np.float32)
    G = refit.shape[1]
    panel = X['panel']
    pf = np.zeros((len(panel), base.table.features.shape[1]), np.float32); has_pf = np.zeros(len(panel), bool)
    for j, t in enumerate(panel):
        f = base.feature_map.get(D.canonical_target(t))
        if f is not None:
            pf[j] = f; has_pf[j] = True
    prog = np.zeros((len(panel), G), np.float32)
    idx = np.flatnonzero(has_pf)
    for s0 in range(0, len(idx), 2048):
        j = idx[s0:s0 + 2048]
        prog[j] = fitted.predict([D.canonical_target(panel[q]) for q in j], pf[j])[0]
    del sub, fitted, y, meta; gc.collect()
    E = X['E']; present = np.isfinite(np.asarray(E[:, :, 0]))
    ppos = {D.canonical_target(t): i for i, t in enumerate(panel)}
    return dict(h=h, tr_all=tr_all, held=held, targets=targets, canon=canon, refit=refit, G=G,
                prog=prog, has_pf=has_pf, present=present, panel_idx=np.array([ppos.get(t, -1) for t in canon]),
                rows=int(len(ids)))


def prepare(gf, tr, X, rank, seed):
    """main effect, residual basis and perturbation programme from TRAINING contexts `tr` only."""
    E = X['E']; panel_n = E.shape[1]; G = gf['G']
    acc = np.zeros((panel_n, G), np.float64); cnt = np.zeros(panel_n)
    for i in tr:
        m = gf['present'][i]; acc[m] += np.asarray(E[i][m], np.float64); cnt[m] += 1
    ok = cnt >= 2
    mu = np.zeros((panel_n, G), np.float32); mu[ok] = (acc[ok] / cnt[ok, None]).astype(np.float32)
    del acc
    usable = ok & gf['has_pf']
    rng = np.random.default_rng(seed)
    sample = []
    for i in tr:
        m = np.flatnonzero(gf['present'][i] & usable)
        if len(m) == 0:
            continue
        pick = m[rng.choice(len(m), size=min(1500, len(m)), replace=False)]
        sample.append(np.asarray(E[i][pick], np.float32) - mu[pick])
    pooled = np.concatenate(sample); del sample
    resp_basis = DD.pca(pooled, rank)
    var = float(((pooled @ resp_basis) ** 2).sum() / max((pooled ** 2).sum(), 1e-9)); del pooled
    pert_basis = DD.pca(gf['prog'][usable], rank)
    zp = ((gf['prog'] - gf['prog'][usable].mean(0)) @ pert_basis).astype(np.float32)
    return dict(mu=mu, usable=usable, resp_basis=resp_basis, zp=zp, var=var, n_usable=int(usable.sum()))


def gene_task(gf, pr, ctx, cs, X):
    Zs, Ps, Ys = [], [], []
    for i in ctx:
        m = np.flatnonzero(gf['present'][i] & pr['usable'])
        zc = cs['Z'][cs['index'][f"gene:{X['bios'][i]}"]]
        for s0 in range(0, len(m), 4096):
            j = m[s0:s0 + 4096]
            d = np.asarray(X['E'][i][j], np.float32) - pr['mu'][j]
            Ys.append(d @ pr['resp_basis']); Ps.append(pr['zp'][j]); Zs.append(np.tile(zc, (len(j), 1)))
    zc = np.concatenate(Zs).astype(np.float32)
    return dict(zc=zc, zp=np.concatenate(Ps).astype(np.float32), y=np.concatenate(Ys).astype(np.float32),
                w=np.ones(len(zc), np.float32))


def train(gf, pr, tr, inner, cs, X, spec, device, log):
    fit = [i for i in tr if i != inner]
    t_full, t_fit, t_val = (gene_task(gf, pr, c, cs, X) for c in (tr, fit, [inner]))
    torch.manual_seed(spec.seed + 1)
    m = M.InteractionModel(spec, ['tahoe', 'gse', 'gene'])
    rec = M.train(m, dict(gene=t_fit), replace(spec, steps=T.GENE_STEPS), device, weights={'gene': 1.0},
                  val=dict(gene=t_val), log=log)
    M.train(m, dict(gene=t_full), replace(spec, steps=max(rec['best_step'], spec.eval_every)), device,
            weights={'gene': 1.0}, val=None)
    info = dict(best_step=rec['best_step'], best_val=rec['best_val'], zero_residual_val=rec['zero_residual_val'],
                beat_zero=rec['beat_zero'], params=rec['params'], rows=int(len(t_full['y'])),
                inner=X['bios'][inner], n_train_contexts=len(tr))
    return m, info


def residual(model, gf, pr, spec, zc_vec, device):
    keep = np.flatnonzero((gf['panel_idx'] >= 0) & pr['usable'][np.maximum(gf['panel_idx'], 0)])
    out = np.zeros((len(gf['targets']), gf['G']), np.float32)
    zp_t = torch.as_tensor(pr['zp'][gf['panel_idx'][keep]], device=device)
    zc = torch.as_tensor(np.tile(zc_vec.astype(np.float32), (len(keep), 1)), device=device)
    with torch.no_grad():
        lat = model('gene', zc, zp_t).cpu().numpy()
    out[keep] = (lat @ pr['resp_basis'].T).astype(np.float32)
    return out, keep


def truth(gf, pr, X):
    """held-out context's observed residual (vs training mu) on the line's panel targets, for
    the cheap residual-cosine diagnostics; never used in any fit."""
    keep = np.flatnonzero((gf['panel_idx'] >= 0) & pr['usable'][np.maximum(gf['panel_idx'], 0)])
    obs = np.asarray(X['E'][gf['h']][gf['panel_idx'][keep]], np.float32)
    ok = np.isfinite(obs[:, 0])
    return keep[ok], obs[ok] - pr['mu'][gf['panel_idx'][keep[ok]]], pr['mu'][gf['panel_idx'][keep[ok]]]


def cos(a, b):
    return float(np.dot(a.ravel(), b.ravel()) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-9))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lines', default=','.join(T.LINES))
    ap.add_argument('--seeds', default='0,1,2')
    ap.add_argument('--curve', action='store_true', help='context-count scaling instead of the main arms')
    ap.add_argument('--device', default='cuda')
    a = ap.parse_args()
    device = torch.device(a.device); seeds = [int(s) for s in a.seeds.split(',')]
    R.OUT = ROOT
    X = load_assets(); base = A.Data(); table = D.Table.load(ROOT / 'table_v1')
    indices = {g: i for i, g in enumerate(V.official_features())}
    columns = np.array([indices[g] for g in table.genes])
    spec0 = M.Spec(kind='film', rank=32, weight_decay=1e-3)
    rows_path = ROOT / 'scores.json'; rows = json.loads(rows_path.read_text()) if rows_path.exists() else []
    tag_prefix = 'curve' if a.curve else 'film'
    diag_path = ROOT / ('curve_diagnostics.json' if a.curve else 'film_diagnostics.json')
    diags = json.loads(diag_path.read_text()) if diag_path.exists() else []
    curves = []
    worker = R.Worker('extfilm')
    try:
        for line in a.lines.split(','):
            t0 = time.monotonic(); cb = D.canonical_bio(line)
            gf = fold(line, None, X, base, table, 32, device)
            old_tr = [BIOS7.index(b) for b in BIOS7 if b != cb]
            inner_bio = BIOS7[old_tr[len(old_tr) // 2]]
            inner = X['bios'].index(inner_bio)
            wrong_bio = BIOS7[old_tr[0]]
            drop_drug = T.CONFLICT.get(cb, {'tahoe': [], 'gse': []})
            drop = {'tahoe': drop_drug['tahoe'], 'gse': drop_drug['gse']}
            # fixed context order for the curve: old biologies, CD4 donor x state, HipSci lines
            order = ([i for i in gf['tr_all'] if X['bios'][i] in BIOS7]
                     + sorted(i for i in gf['tr_all'] if X['bios'][i].startswith('CD4D'))
                     + sorted(i for i in gf['tr_all'] if X['bios'][i].startswith('FENGL')))
            sizes = [6, 12, 24, len(order)] if a.curve else [len(order)]
            refit_v1 = np.load(ROOT / line / 'refit_v1.npy') if (ROOT / line / 'refit_v1.npy').exists() else None
            for n in sizes:
                tr = order[:n]
                if inner not in tr:
                    tr = [inner] + tr[:-1]
                train_bios = [X['bios'][i] for i in tr]
                cs = context_space(train_bios, X['adj'], X['bios'], drop, dim=64)
                spec = replace(spec0, ctx_dim=cs['dim'])
                pr = prepare(gf, tr, X, 32, 20260921 + gf['h'] + 1000 * n)
                log = lambda e, n=n: curves.append(dict(line=line, n=n, **e))
                m, info = train(gf, pr, tr, inner, cs, X, spec, device, log)
                zc = lambda b: cs['Z'][cs['index'][f'gene:{b}']]
                res = {k: residual(m, gf, pr, spec, v, device)[0] for k, v in
                       [('correct', zc(cb)), ('wrong', zc(wrong_bio)), ('zero', np.zeros(spec.ctx_dim, np.float32))]}
                keep_t, tru, mu_t = truth(gf, pr, X)
                cheap = {k: cos(res[k][keep_t], tru) for k in res}
                base_pred = mu_t; tot = {k: cos(base_pred + res[k][keep_t], base_pred + tru) - cos(base_pred, base_pred + tru)
                                         for k in res}
                d = dict(line=line, train_contexts=train_bios,
                         wrong=wrong_bio, ctx_dim=cs['dim'], n_drug_ctx=cs['n_drug'], panel_usable=pr['n_usable'],
                         resp_basis_var=pr['var'], refit_rows=gf['rows'], **info,
                         resid_cos=cheap, total_cos_gain=tot,
                         correct_minus_wrong=cheap['correct'] - cheap['wrong'],
                         residual_norm_ratio=float(np.linalg.norm(res['correct']) / np.linalg.norm(gf['refit'])))
                # predictions on the refit_v1 base (identical fit; assert) and harness PDS
                pred_base = np.load(ROOT / line / 'mean_v1.npy') if (ROOT / line / 'mean_v1.npy').exists() else None
                if pred_base is None:
                    raise FileNotFoundError('run ext_baselines.py first (mean_v1 / refit_v1)')
                pred_base = pred_base.copy(); pred_base[:, columns] = gf['refit']
                if refit_v1 is not None:
                    d['max_abs_diff_refit_vs_refit_v1'] = float(np.abs(pred_base - refit_v1).max())
                plans = []
                if a.curve:
                    plans += [(f'{tag_prefix}{n}_a0.25', 0.25, res['correct'], [0])]
                else:
                    plans += [(f'film_ext_a{al:g}', al, res['correct'], seeds) for al in ALPHAS]
                    plans += [(f'film_ext_wrong_a1', 1.0, res['wrong'], seeds), (f'film_ext_zero_a1', 1.0, res['zero'], seeds),
                              (f'film_ext_wrong_a0.25', 0.25, res['wrong'], seeds), (f'film_ext_zero_a0.25', 0.25, res['zero'], seeds)]
                pds = {}
                for tag, al, r, sd in plans:
                    pred = pred_base.copy(); pred[:, columns] = gf['refit'] + al * r
                    path = ROOT / line / f'{tag}.npy'; np.save(path, pred.astype(np.float32))
                    for seed in sd:
                        if any(x['line'] == line and x['tag'] == tag and x['seed'] == seed for x in rows):
                            continue
                        rows.append(dict(line=line, tag=tag, alpha=al, seed=seed, **XB.score_one(worker, gf['targets'], line, path, seed)))
                        P.write(rows_path, rows)
                    pds[tag] = float(np.mean([x['pds'] for x in rows if x['line'] == line and x['tag'] == tag]))
                d['pds'] = pds
                diags.append(d); P.write(diag_path, diags); P.write(ROOT / f'{tag_prefix}_curves.json', curves)
                print(json.dumps(dict(event='arm', line=line, n=len(tr), best_step=info['best_step'], beat_zero=info['beat_zero'],
                                      resid_cos=cheap, cmw=d['correct_minus_wrong'], pds=pds,
                                      secs=round(time.monotonic() - t0, 1)), default=float), flush=True)
                del m, pr, cs; gc.collect()
            del gf; gc.collect()
    finally:
        worker.close()
    print(json.dumps(dict(event='done')), flush=True)


if __name__ == '__main__':
    main()
