"""S4: refit on the extended table (v2, contexts nested in CD4/FENG) + 0.25 x FiLM context
residual trained on all 38 gene contexts, per official context A/B/C.

Authorised by the user 2026-09-22 ("해봐") after `reports/ext_contexts_2026-09-21/FINDINGS.md`
ranked film_v2 first on the local official-preset composite (+0.0247 vs refit, 15/15 seeds, 5/5
lines). Mirrors exactly what was scored: the residual's perturbation programme and training data
come from the v1 table (as in `ext_film.py`), the base is the v2-table refit (as in film_v2), alpha
0.25 (nested choice on 4/5 lines). Context profiles for A/B/C are pooled official control counts
-> log1p CP10K on the 5,592 table genes, projected with the PCA fitted on the 38 training contexts.
`--verify K562` reproduces the stored leave-one-biology-out residual with the same code path.
"""
from pathlib import Path
import pathlib
from dataclasses import replace
import argparse, gc, json, sys, time
import numpy as np, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ext_film as XF
import ext_baselines as XB
import dcp_data as DD
import dcp_model as M
import b19_optimization_data as A
import b19_optimization_report as P
import context_forcing_data as F
import simple_residual_data as D
import vcclib as V
import fullaxis_lean_baseline as LB

import os
ROOT = D.ROOT / 'reports/ext_contexts_2026-09-21'
# transfer_v2 copy of ext_submit_build.py: output dir, context names and unseen-target handling are configurable so
# the S4 refit+FiLM part can be rebuilt for a new panel/contexts (final phase) WITHOUT touching the S4 files.
OUT = pathlib.Path(os.environ['TV2_S4_OUT'])
CTX = tuple(os.environ.get('TV2_CTX', 'A,B,C').split(',')); ALPHA = 0.25; BIOS7 = XF.BIOS7
ALLOW_UNSEEN = os.environ.get('TV2_ALLOW_UNSEEN', '1') == '1'
OUT.mkdir(parents=True, exist_ok=True)


def log(**kw):
    print(json.dumps(kw, default=float), flush=True)


def refit_on(table, ids, targets, feats, canon, rank=32):
    y = np.asarray(table.y[ids], np.float32)
    meta = {k: np.asarray(getattr(table, k))[ids] for k in ('targets', 'bios', 'modality', 'features', 'weights', 'sources')}
    sub = D.Table(y=y, controls=np.zeros((len(ids), 1), np.float32), genes=list(table.genes), **meta)
    fitted = LB.LeanSourceBaseline(sub, np.arange(len(ids)), rank=rank)
    common, seen = fitted.predict(canon, feats)
    if not seen.all():
        if not ALLOW_UNSEEN:
            raise ValueError(f'{int((~seen).sum())} targets unseen by the baseline')
        log(event='unseen_targets_zeroed', n=int((~seen).sum()))
        common = np.where(seen[:, None], common, 0.0)
    del sub, y, meta; gc.collect()
    return common.astype(np.float32), fitted


def full_fold(X, base, table1, device):
    """ext_film.fold with nothing held out."""
    ids = F.permitted_ids(table1, [])
    targets = V.official_targets(); feats = base.features(targets)
    canon = [D.canonical_target(t) for t in targets]
    refit, fitted = refit_on(table1, ids, targets, feats, canon)
    G = refit.shape[1]; panel = X['panel']
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
    del fitted; gc.collect()
    E = X['E']; present = np.isfinite(np.asarray(E[:, :, 0]))
    ppos = {D.canonical_target(t): i for i, t in enumerate(panel)}
    return dict(h=None, tr_all=list(range(len(X['bios']))), held=set(), targets=targets, canon=canon, refit=refit,
                G=G, prog=prog, has_pf=has_pf, present=present,
                panel_idx=np.array([ppos.get(t, -1) for t in canon]), rows=int(len(ids)), feats=feats)


def context_space_q(train_bios, adj, bios, extra, dim=64):
    """XF.context_space plus projection of new gene-assay contexts (raw table-axis profiles)."""
    keys = [f'gene:{b}' for b in bios]; idx = {k: i for i, k in enumerate(keys)}
    Xp = np.asarray(adj, np.float64); tr = np.array([idx[f'gene:{b}'] for b in train_bios])
    cm = Xp[tr].mean(0); Bm = DD.pca(Xp[tr] - cm, min(dim, len(tr) - 1))
    Z = ((Xp - cm) @ Bm).astype(np.float32); s = float(np.sqrt(np.mean(Z[tr] ** 2))); Z = Z / max(s, 1e-9)
    ref = XF.context_space(train_bios, adj, bios, {}, dim); np.testing.assert_array_equal(Z, ref['Z'])
    Zq = {k: (((v - cm) @ Bm).astype(np.float32) / max(s, 1e-9)) for k, v in extra.items()}
    return dict(Z=Z, index=idx, dim=int(Bm.shape[1]), Zq=Zq)


def official_profile(c, genes):
    """pooled official control counts -> log1p CP10K over all official features -> table genes."""
    import anndata as ad, scipy.sparse as sp
    a = ad.read_h5ad(V.OFFICIAL_DIR / f'context_{c}.h5ad')
    if [str(x) for x in a.var_names] != list(V.official_features()):
        raise ValueError(f'context {c}: var order differs from gene_names.csv')
    csr = sp.csr_matrix(a.X)
    if not np.all(csr.data == np.floor(csr.data)) or csr.data.min() < 0:
        raise ValueError(f'context {c}: not integral counts')
    tot = np.asarray(csr.sum(0)).ravel().astype(np.float64)
    raw = np.log1p(1e4 * tot / tot.sum())
    gi = {g: i for i, g in enumerate(V.official_features())}
    prof = raw[[gi[g] for g in genes]].astype(np.float32)
    n = int(csr.shape[0]); del a, csr; gc.collect()
    return prof, n


def verify(line, X, base, table1, device, spec0):
    gf = XF.fold(line, None, X, base, table1, 32, device)
    cb = D.canonical_bio(line)
    old_tr = [BIOS7.index(b) for b in BIOS7 if b != cb]
    inner = X['bios'].index(BIOS7[old_tr[len(old_tr) // 2]])
    order = ([i for i in gf['tr_all'] if X['bios'][i] in BIOS7] + sorted(i for i in gf['tr_all'] if X['bios'][i].startswith('CD4D'))
             + sorted(i for i in gf['tr_all'] if X['bios'][i].startswith('FENGL')))
    tr = order; n = len(order)
    cs = XF.context_space([X['bios'][i] for i in tr], X['adj'], X['bios'], {}, 64)
    spec = replace(spec0, ctx_dim=cs['dim'])
    pr = XF.prepare(gf, tr, X, 32, 20260921 + gf['h'] + 1000 * n)
    m, info = XF.train(gf, pr, tr, inner, cs, X, spec, device, lambda e: None)
    res, keep = XF.residual(m, gf, pr, spec, cs['Z'][cs['index'][f'gene:{cb}']], device)
    indices = {g: i for i, g in enumerate(V.official_features())}; columns = np.array([indices[g] for g in table1.genes])
    stored = (np.load(ROOT / line / 'film_ext_a1.npy') - np.load(ROOT / line / 'refit_v1.npy'))[:, columns]
    out = dict(line=line, max_abs_diff_vs_stored_residual=float(np.abs(res - stored).max()),
               corr_vs_stored=float(np.corrcoef(res.ravel(), stored.ravel())[0, 1]),
               residual_norm=float(np.linalg.norm(res)), **info)
    log(event='verify', **out); return out


def full(X, base, device, spec0):
    OUT.mkdir(parents=True, exist_ok=True)
    table1 = D.Table.load(ROOT / 'table_v1'); table2 = D.Table.load(ROOT / 'table_v2')
    indices = {g: i for i, g in enumerate(V.official_features())}; columns = np.array([indices[g] for g in table1.genes])
    gf = full_fold(X, base, table1, device)
    log(event='refit_v1_full', rows=gf['rows'], targets_in_panel=int(((gf['panel_idx'] >= 0)).sum()))
    inner = X['bios'].index('JURKAT')
    tr = gf['tr_all']; n = len(tr)
    extra = {}; ctx_diag = {}
    for c in CTX:
        extra[c], ncell = official_profile(c, table1.genes); ctx_diag[c] = dict(cells=ncell)
    cs = context_space_q([X['bios'][i] for i in tr], X['adj'], X['bios'], extra, 64)
    Zg = cs['Z']; C = np.asarray(np.load(ROOT / 'ext_gene_contexts.npz', allow_pickle=True)['profile_raw'], np.float64)
    for c in CTX:
        d = np.linalg.norm(Zg - cs['Zq'][c][None], axis=1); cc = [float(np.corrcoef(extra[c], C[i])[0, 1]) for i in range(len(C))]
        ctx_diag[c].update(nearest_context=X['bios'][int(d.argmin())], nn_distance=float(d.min()),
                           gene_nn_median=float(np.median(np.sort(np.linalg.norm(Zg[:, None] - Zg[None], axis=-1), 1)[:, 1])),
                           norm_z=float(np.linalg.norm(cs['Zq'][c])), gene_norm_z_median=float(np.median(np.linalg.norm(Zg, axis=1))),
                           profile_corr_median=float(np.median(cc)), profile_corr_max=float(np.max(cc)),
                           profile_top_genes=[table1.genes[i] for i in np.argsort(-extra[c])[:5]])
    log(event='context_space', dim=cs['dim'], contexts=ctx_diag)
    spec = replace(spec0, ctx_dim=cs['dim'])
    pr = XF.prepare(gf, tr, X, 32, 20260921 + 7 + 1000 * n)
    curve = []
    m, info = XF.train(gf, pr, tr, inner, cs, X, spec, device, lambda e: curve.append(e))
    log(event='trained', **info)
    torch.save(dict(state_dict=m.state_dict(), spec=spec.__dict__, resp_basis=pr['resp_basis'], info=info), OUT / 'S4_film_model.pt')
    json.dump(curve, open(OUT / 'S4_training_curve.json', 'w'), indent=1)
    res = {c: XF.residual(m, gf, pr, spec, cs['Zq'][c], device)[0] for c in CTX}
    res['zero'] = XF.residual(m, gf, pr, spec, np.zeros(spec.ctx_dim, np.float32), device)[0]
    for b in BIOS7:
        res[f'gene:{b}'] = XF.residual(m, gf, pr, spec, cs['Z'][cs['index'][f'gene:{b}']], device)[0]
    keep = np.flatnonzero((gf['panel_idx'] >= 0) & pr['usable'][np.maximum(gf['panel_idx'], 0)])
    # base: refit on the v2 table, all rows; outside the table genes the v2 mean transfer (as scored)
    ids2 = F.permitted_ids(table2, [])
    refit2, _ = refit_on(table2, ids2, gf['targets'], gf['feats'], gf['canon'])
    emp2, nv = XB.empirical(gf['targets'], [], 'v2')
    base2 = emp2.copy(); base2[:, columns] = refit2
    np.save(OUT / 'S4_refit_v2.delta.npy', base2.astype(np.float32))
    s1 = np.load(D.ROOT / 'reports/submission_2026-09-20/S1_refit.delta.npy')
    same = s1.shape == base2.shape
    log(event='refit_v2_full', rows=int(len(ids2)), voters_outside=nv, corr_vs_S1_refit=float(np.corrcoef(base2.ravel(), s1.ravel())[0, 1]) if same else None,
        norm=float(np.linalg.norm(base2)), norm_S1=float(np.linalg.norm(s1)) if same else None)
    cos = lambda a, b: float(np.dot(a.ravel(), b.ravel()) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-9))
    diagnostics = []
    outside = np.ones(base2.shape[1], bool); outside[columns] = False
    for c in CTX:
        delta = base2.copy(); delta[:, columns] += ALPHA * res[c]
        np.testing.assert_array_equal(delta[:, outside], base2[:, outside])
        if delta.shape != (300, 18533) or not np.isfinite(delta).all():
            raise ValueError(f'invalid S4 delta {c}')
        np.save(OUT / f'S4_ctx{c}.delta.npy', delta.astype(np.float32))
        diagnostics.append(dict(context=c, targets_with_residual=int(len(keep)),
                                residual_norm=float(np.linalg.norm(ALPHA * res[c])), refit_norm=float(np.linalg.norm(refit2)),
                                ratio=float(np.linalg.norm(ALPHA * res[c]) / np.linalg.norm(refit2)),
                                cos_to_other_contexts={o: cos(res[c], res[o]) for o in CTX if o != c},
                                cos_to_zero_context=cos(res[c], res['zero']),
                                cos_to_gene_contexts={b: cos(res[c], res[f'gene:{b}']) for b in BIOS7}, **ctx_diag[c]))
        log(event='s4_built', **diagnostics[-1])
    P.freeze(OUT / 'build_receipt.json', dict(
        authorised_by_user='2026-09-22: "해봐" -- submit film_v2 (FiLM residual on the extended refit)',
        source_experiment=str(ROOT), local_official_composite=dict(refit_cur=0.2755, refit_v2=0.2815, film_v2=0.3001,
                                                                     delta_vs_refit=0.0247, seeds='15/15', lines='5/5'),
        alpha=ALPHA, spec=spec.__dict__, training=info,
        refit_v2_rows=int(len(ids2)), refit_v1_rows_for_programme=gf["rows"], contexts=diagnostics,
        scale_for_submission=3.0, official_upload=False))
    log(event='build_complete')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', default=''); ap.add_argument('--full', action='store_true')
    ap.add_argument('--device', default='cuda'); a = ap.parse_args()
    device = torch.device(a.device); spec0 = M.Spec(kind='film', rank=32, weight_decay=1e-3)
    X = XF.load_assets(); base = A.Data()
    if a.verify:
        OUT.mkdir(parents=True, exist_ok=True)
        r = verify(a.verify, X, base, D.Table.load(ROOT / 'table_v1'), device, spec0)
        json.dump(r, open(OUT / f'verify_{a.verify}.json', 'w'), indent=2)
        if not (r['max_abs_diff_vs_stored_residual'] < 1e-3 and r['corr_vs_stored'] > 0.999):
            raise SystemExit(f'VERIFY FAILED: {r}')
    if a.full:
        full(X, base, device, spec0)


if __name__ == '__main__':
    main()
