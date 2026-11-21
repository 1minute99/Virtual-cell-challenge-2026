"""Architecture screen for perturbation x context prediction on the extended gene contexts.

Fixed: held-out K562, seed 0, 35 training contexts (34 for fitting + JURKAT for early stopping),
identical inputs, bases and evaluation for every architecture. Residual family: prediction =
refit_v2 + alpha * residual (alpha 1 for the official runs). Joint family: prediction = decoded
full response, refit unused. Generous early stopping: >= 2,000 steps, patience 30 evals (3,000
steps), max 20,000 steps. Outputs harness PDS (seed 0) and deltas for correct / wrong / zero
context; `arch_preset.sh` runs the official vcc2026 preset on them.
"""
from pathlib import Path
import argparse, gc, json, sys, time
import numpy as np, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ext_film as XF
import ext_baselines as XB
import arch_zoo as Z
import dcp_data as DD
import b19_optimization_run as R
import b19_optimization_data as A
import b19_optimization_report as P
import simple_residual_data as D
import vcclib as V

EXT = D.ROOT / 'reports/ext_contexts_2026-09-21'
OUT = D.ROOT / 'reports/arch_screen_2026-09-22'
SEED = 0; K_FULL = 64; T_TOK = 32
LINE = 'K562'   # overridden by --line


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--archs', default=','.join(Z.ARCHS)); ap.add_argument('--device', default='cuda')
    ap.add_argument('--line', default='K562'); ap.add_argument('--out', default='')
    a = ap.parse_args(); device = torch.device(a.device)
    global LINE, OUT
    LINE = a.line; OUT = Path(a.out) if a.out else OUT
    OUT.mkdir(parents=True, exist_ok=True); (OUT / LINE).mkdir(exist_ok=True); R.OUT = OUT
    X = XF.load_assets(); base = A.Data(); table1 = D.Table.load(EXT / 'table_v1')
    indices = {g: i for i, g in enumerate(V.official_features())}; columns = np.array([indices[g] for g in table1.genes])
    gf = XF.fold(LINE, None, X, base, table1, 32, device)
    cb = D.canonical_bio(LINE); bios = X['bios']
    tr = gf['tr_all']; inner = bios.index('JURKAT'); fit = [i for i in tr if i != inner]
    old_tr = [b for b in XF.BIOS7 if b != cb]; wrong_bio = old_tr[0]     # as in ext_film: first old training biology
    cs = XF.context_space([bios[i] for i in tr], X['adj'], bios, {}, 64)
    pr = XF.prepare(gf, tr, X, 32, 20260922)
    E = X['E']; G = gf['G']
    # full-response basis (rank 64) from training effects, same sampling as prepare
    rng = np.random.default_rng(20260922 + 1); sample = []
    for i in tr:
        m = np.flatnonzero(gf['present'][i] & pr['usable']); pick = m[rng.choice(len(m), size=min(1500, len(m)), replace=False)]
        sample.append(np.asarray(E[i][pick], np.float32))
    full_basis = DD.pca(np.concatenate(sample), K_FULL); del sample
    # context tokens: centred pooled profile cut into 32 gene blocks
    prof = np.asarray(np.load(EXT / 'ext_gene_contexts.npz', allow_pickle=True)['profile_raw'], np.float32)
    prof = prof - prof[tr].mean(0); dtok = int(np.ceil(prof.shape[1] / T_TOK)); pad = T_TOK * dtok - prof.shape[1]
    TOK = torch.as_tensor(np.pad(prof, ((0, 0), (0, pad))).reshape(len(bios), T_TOK, dtok) / (prof.std() + 1e-6), device=device)
    TOK = torch.cat([TOK, torch.zeros(1, T_TOK, dtok, device=device)])          # index len(bios) = zero context
    Zc = torch.as_tensor(np.vstack([cs['Z'], np.zeros((1, cs['dim']), np.float32)]), device=device)
    pf_panel = np.zeros((len(X['panel']), 64), np.float32)
    for j, t in enumerate(X['panel']):
        f = base.feature_map.get(D.canonical_target(t))
        if f is not None: pf_panel[j] = f
    def rows(ctxs):
        zc, ctx, zp, pf, yr, yf, x0 = [], [], [], [], [], [], []
        for i in ctxs:
            m = np.flatnonzero(gf['present'][i] & pr['usable'])
            for s0 in range(0, len(m), 4096):
                j = m[s0:s0 + 4096]; e = np.asarray(E[i][j], np.float32)
                yr.append((e - pr['mu'][j]) @ pr['resp_basis']); yf.append(e @ full_basis); x0.append(pr['mu'][j] @ full_basis)
                zp.append(pr['zp'][j]); pf.append(pf_panel[j]); ctx.append(np.full(len(j), i))
        t = lambda arr: torch.as_tensor(np.concatenate(arr).astype(np.float32), device=device)
        return dict(ctx=torch.as_tensor(np.concatenate(ctx), device=device), zp=t(zp), pf=t(pf), yr=t(yr), yf=t(yf), x0=t(x0))
    Tfit, Tval = rows(fit), rows([inner])
    def aug(b):
        b = dict(b); b['zc'] = Zc[b['ctx']]; b['tok'] = TOK[b['ctx']]; return b
    # held-out K562 rows on the harness targets
    keep = np.flatnonzero((gf['panel_idx'] >= 0) & pr['usable'][np.maximum(gf['panel_idx'], 0)]); pj = gf['panel_idx'][keep]
    H = dict(zp=torch.as_tensor(pr['zp'][pj], device=device), pf=torch.as_tensor(pf_panel[pj], device=device),
             x0=torch.as_tensor((pr['mu'][pj] @ full_basis).astype(np.float32), device=device))
    keep_t, tru, mu_t = XF.truth(gf, pr, X)
    base_delta = np.load(EXT / LINE / 'refit_v2.npy')
    rows_path = OUT / 'scores.json'; scores = json.loads(rows_path.read_text()) if rows_path.exists() else []
    results = json.loads((OUT / 'results.json').read_text()) if (OUT / 'results.json').exists() else {}
    worker = R.Worker('arch')
    cos = lambda u, v: float(np.dot(u.ravel(), v.ravel()) / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-9))
    def score(tag, seed=SEED):
        if any(r['line'] == LINE and r['tag'] == tag and r['seed'] == seed for r in scores): 
            return [r['pds'] for r in scores if r['line'] == LINE and r['tag'] == tag and r['seed'] == seed][0]
        r = dict(line=LINE, tag=tag, alpha=1.0, seed=seed, **XB.score_one(worker, gf['targets'], LINE, OUT / LINE / f'{tag}.npy', seed))
        scores.append(r); P.write(rows_path, scores); return r['pds']
    try:
        np.save(OUT / LINE / 'base_refit_v2.npy', base_delta); results.setdefault('base_refit_v2', {})['pds'] = score('base_refit_v2')
        P.write(OUT / 'results.json', results)
        print(json.dumps(dict(event='base', pds=results['base_refit_v2']['pds'], fit_rows=int(len(Tfit['zp'])), val_rows=int(len(Tval['zp'])),
                              ctx_dim=cs['dim'], tokens=[T_TOK, dtok], K_full=K_FULL, held_targets=int(len(keep)))), flush=True)
        for arch in a.archs.split(','):
            if arch in results and 'pds_correct' in results[arch]:
                continue
            t0 = time.monotonic(); torch.manual_seed(SEED)
            model = Z.build(arch, cs['dim'], dtok, K_FULL); res_family = model.family == 'residual'
            T = {k: v for k, v in Tfit.items() if k not in ('yr', 'yf')}; T['y'] = Tfit['yr'] if res_family else Tfit['yf']
            Vv = {k: v for k, v in Tval.items() if k not in ('yr', 'yf')}; Vv['y'] = Tval['yr'] if res_family else Tval['yf']
            curve = []
            rec = Z.train(model, T, Vv, device, seed=SEED, aug=aug, log=lambda e: curve.append(e))
            model.eval(); out = dict(arch=arch, family=model.family, params=rec['params'], best_step=rec['best_step'], steps_run=rec['steps_run'],
                                     best_val=rec['best_val'], zero_val=rec['zero_val'], beat_zero=rec['beat_zero'], secs=round(time.monotonic() - t0, 1))
            json.dump(curve, open(OUT / f'curve_{arch}.json', 'w'))
            preds = {}
            for name, ci in (('correct', bios.index(cb)), ('wrong', bios.index(wrong_bio)), ('zero', len(bios))):
                b = dict(H); b['ctx'] = torch.full((len(keep),), ci, device=device, dtype=torch.long); b = aug(b)
                with torch.no_grad(): lat = model.predict(b).cpu().numpy()
                full = np.zeros((len(gf['targets']), G), np.float32)
                if res_family:
                    full[keep] = lat @ pr['resp_basis'].T; resid = full
                else:
                    full[keep] = lat @ full_basis.T; resid = full.copy(); resid[keep] -= pr['mu'][pj]
                preds[name] = (full, resid)
                out[f'resid_cos_{name}'] = cos(resid[keep_t], tru)
                out[f'total_cos_gain_{name}'] = cos(mu_t + resid[keep_t], mu_t + tru) - cos(mu_t, mu_t + tru)
                for al in ((1.0, 0.25) if res_family else (1.0,)):
                    delta = base_delta.copy()
                    if res_family: delta[:, columns] = gf_base_cols(base_delta, columns) + al * full
                    else: delta[:, columns] = full
                    if not np.isfinite(delta).all(): raise ValueError(f'{arch}/{name}: non-finite delta')
                    tag = f'{arch}_{name}' + ('' if al == 1.0 else f'_a{al:g}')
                    np.save(OUT / LINE / f'{tag}.npy', delta.astype(np.float32))
                    out[f'pds_{name}' + ('' if al == 1.0 else f'_a{al:g}')] = score(tag)
            out['residual_norm_ratio'] = float(np.linalg.norm(preds['correct'][1]) / np.linalg.norm(gf_base_cols(base_delta, columns)))
            out['correct_minus_wrong_cos'] = out['resid_cos_correct'] - out['resid_cos_wrong']
            results[arch] = out; P.write(OUT / 'results.json', results)
            print(json.dumps(dict(event='arch', **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()})), flush=True)
            del model; gc.collect(); torch.cuda.empty_cache()
    finally:
        worker.close()
    print(json.dumps(dict(event='done')), flush=True)


def gf_base_cols(base_delta, columns):
    return base_delta[:, columns]


if __name__ == '__main__':
    main()
