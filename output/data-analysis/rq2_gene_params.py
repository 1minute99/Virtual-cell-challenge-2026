"""RQ2 Task 9 -- Proposition-2 parameters on the public GENETIC screens and the H2c verdict.

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md -- H2c, operational definitions 3-5, Task 9, Amendments
A2 (C_public counted two ways) and A3. Pre-specifications in knowledge/_meta/log.md (2026-11-05): d_gene over the
Count-1 contexts (capped at C - 1, reported as conservative); H2c against BOTH counts; secondary C*_min with the Count-2
donor-line contexts; secondary sensitivities S1 (target b + i), S2 (C*_b), S3 (line bootstrap recomputing d) that do
NOT change the verdict rule.

Inputs: results/rq2_gene_cache/<SRC>.npz written by rq2_gene_extract.py (cell-level split halves, 5 split seeds,
separate control halves), results/families.json (Count 1), results/rq2_params_gene.json['count2'] (Count 2, written
by rq2_gene_count2.py BEFORE this script ran). Estimators are the Task 8 library rq2_params.py (unchanged):
participation_ratio, descriptor_pcs, filter_min_lines, halves_stats, aggregate (B = 2,000 over targets x split seeds),
c_star_min; S1-S3 via halves_extra_stats / sensitivity_estimates / line_bootstrap.

ASSUMPTIONS (stated here, checked where possible; note: knowledge/synthesis/rq2_gene_params-assumptions.md)
  G1 One context per Count-1 line (families.json line_to_families): K562 (GWPS), RPE1, HepG2, Jurkat (Replogle/Nadig),
     HCT116, HEK293T (Orion), CD4 (3 states pooled), KOLF2.1J, HipSci pool (Feng targeted screen, all 19 lines pooled),
     DLD-1, H1 = H1GW (the counted H1 family; VCC2025-H1 enters only the sensitivities).
  G1b The HipSci pool's response for a target = the 1/(1/n_t + 1/n_c)-weighted mean of the WITHIN-LINE responses
     (each line against its own control halves), so that the line mix of target cells vs control cells cannot create
     a spurious response; its descriptor = control-weighted mean of the lines' control profiles. The registry-style
     summed pseudobulk is a sensitivity. Per (line, target) groups below 20 cells were not extracted.
  G2 Response of (line, target) in half h = log1p CP10K(target cells of h) - log1p CP10K(control cells of h of the same
     line); native-library denominators; no cross-half baseline sharing. A target enters a line with >= 20 cells; the
     analysis keeps targets present in >= 3 lines (rq2_params.filter_min_lines).
  G3 Gene axis = official-axis genes measured in every source (rq2_gene_extract E1); per-element units
     (line, target, gene) in (log1p CP10K)^2, as on Tahoe.
  G4 Descriptor = full-depth control pseudobulk (log1p CP10K) of each line on the gene axis (definition 3);
     d = participation ratio, capped at C - 1.
  G5 sig2_eps rescaled half -> full depth by x 1/2 (rq2_params A1), exact for cell-level halves of floor/ceil size;
     CD4 halves are pseudobulk-unit halves balanced by cell count (deviation; sensitivity drops CD4).
  G6 R^2 / C*_theory are reported only if >= 1 target has n_t - k - 1 >= 1 with k = ceil(d); otherwise 'not estimable'.
Seeds: split SPLIT_BASE 20261108 (extract), bootstrap BOOT_SEED 20261109, S3 seed 20261110. Deterministic.
Run: OMP_NUM_THREADS=4 nice -n 10 python rq2_gene_params.py
"""
from __future__ import annotations
import json, math, sys, time, datetime
from pathlib import Path
import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC'); HERE = ROOT / 'output/data-analysis'
sys.path.insert(0, str(HERE))
import rq2_params as RP

RES = HERE / 'results'; CACHE = RES / 'rq2_gene_cache'
N_SEEDS = 5
MAX_IMBALANCE = 0.25        # |n - 2 nA| / n in every seed (only binds for CD4 unit-level halves; G5)
PRIMARY = ['K562', 'RPE1', 'HEPG2', 'JURKAT', 'HCT116', 'HEK293T', 'CD4', 'KOLF', 'FENG', 'DLD1', 'H1GW']
LINE_NAME = {'K562': 'K562', 'RPE1': 'hTERT-RPE1', 'HEPG2': 'HepG2', 'JURKAT': 'Jurkat', 'HCT116': 'HCT116',
             'HEK293T': 'HEK293T', 'CD4': 'primary CD4+ T', 'KOLF': 'KOLF2.1J iPSC', 'FENG': 'HipSci iPSC pool',
             'DLD1': 'DLD-1', 'H1GW': 'H1 hESC', 'VCC_H1': 'H1 hESC (VCC2025)'}
T0 = time.time()


def log(**k):
    print(json.dumps(dict(k, secs=round(time.time() - T0, 1)), default=str), flush=True)


# ------------------------------------------------------------------ per-context group sums

_cache = {}


def load(src):
    if src not in _cache:
        z = np.load(CACHE / f'{src}.npz')
        _cache[src] = {k: z[k] for k in z.files}
    return _cache[src]


def context_sums(spec):
    """spec: ('src', SRC) | ('feng_line', LINE) | ('sum', [spec, ...]).
    Returns dict target -> index, arrays full (T+1 x G), A (S x T+1 x G), lib, libA, n, nA; last row = control."""
    kind = spec[0]
    if kind == 'sum':
        parts = [context_sums(s) for s in spec[1]]
        tg = sorted(set().union(*[set(p['targets']) for p in parts]))
        G = parts[0]['full'].shape[1]
        out = dict(targets=tg, full=np.zeros((len(tg) + 1, G)), A=np.zeros((N_SEEDS, len(tg) + 1, G), np.float32),
                   lib=np.zeros(len(tg) + 1), libA=np.zeros((N_SEEDS, len(tg) + 1)), n=np.zeros(len(tg) + 1),
                   nA=np.zeros((N_SEEDS, len(tg) + 1)))
        ix = {t: i for i, t in enumerate(tg)}
        for p in parts:
            rows = [ix[t] for t in p['targets']] + [len(tg)]
            for key in ('full', 'lib', 'n'):
                np.add.at(out[key], rows, p[key])
            for key in ('A', 'libA', 'nA'):
                out[key][:, rows] += p[key]
        return out
    if kind == 'src':
        z = load(spec[1]); g = [str(x) for x in z['groups']]
        tsel = [i for i, x in enumerate(g) if x != 'CTRL']; csel = [g.index('CTRL')]
    elif kind == 'feng_line':
        z = load('FENG'); g = [str(x) for x in z['groups']]; pre = spec[1] + '|'
        tsel = [i for i, x in enumerate(g) if x.startswith(pre) and not x.endswith('|CTRL')]
        csel = [g.index(pre + 'CTRL')]
    elif kind == 'feng_pool':                                   # sum over (all or given) Feng lines per target
        lines = spec[1]
        return context_sums(('sum', [('feng_line', l) for l in lines]))
    else:
        raise ValueError(kind)
    sel = tsel + csel
    names = [g[i].split('|', 1)[-1] for i in tsel]
    return dict(targets=names, full=z['full'][sel], A=z['A'][:, sel], lib=z['lib'][sel],
                libA=z['libA'][:, sel], n=z['n'][sel], nA=z['nA'][:, sel])


def feng_lines():
    z = load('FENG')
    return sorted({str(x).split('|')[0] for x in z['groups']})


def _table(spec, min_cells):
    c = context_sums(spec)
    n = c['n'][:-1]; okh = (c['nA'][:, :-1] >= 1).all(0) & ((c['n'][:-1] - c['nA'][:, :-1]) >= 1).all(0)
    imb = (np.abs(c['n'][:-1] - 2 * c['nA'][:, :-1]) / np.maximum(c['n'][:-1], 1)).max(0)
    keep = np.flatnonzero((n >= min_cells) & okh & (imb <= MAX_IMBALANCE))
    desc = RP.log_cp10k(c['full'][-1:], c['lib'][-1:])[0]
    return dict(c=c, keep=keep, n_dropped_imbalance=int(((n >= min_cells) & okh & (imb > MAX_IMBALANCE)).sum()),
                targets=[c['targets'][i] for i in keep], desc=desc, n_ctrl=float(c['n'][-1]), cells=n[keep])


def context_tables(specs, min_cells=20):
    """Per context: targets kept (>= min_cells, both halves non-empty and balanced in every seed) and the control
    descriptor. spec ('merge', [spec, ...]) = one line measured in several datasets: responses are formed WITHIN each
    dataset (own control halves) and averaged per target with weight 1/(1/n_target + 1/n_control); descriptor =
    control-cell-weighted mean of the datasets' control log1p CP10K profiles."""
    T = {}
    for name, spec in specs.items():
        if spec[0] == 'merge':
            parts = [_table(p, min_cells) for p in spec[1]]
            tg = sorted(set().union(*[set(p['targets']) for p in parts]))
            wc = np.array([p['n_ctrl'] for p in parts])
            cells = np.zeros(len(tg)); ix = {t: i for i, t in enumerate(tg)}
            for p in parts:
                for t, n in zip(p['targets'], p['cells']):
                    cells[ix[t]] += n
            T[name] = dict(parts=parts, targets=tg, desc=(wc[:, None] * np.stack([p['desc'] for p in parts])).sum(0) / wc.sum(),
                           n_ctrl=float(wc.sum()), cells=cells,
                           n_dropped_imbalance=int(sum(p['n_dropped_imbalance'] for p in parts)))
        else:
            T[name] = _table(spec, min_cells)
    _cache.clear()
    return T


def _resp(t, seed):
    c, keep = t['c'], t['keep']
    sA = c['A'][seed].astype(np.float64); sF = c['full']
    lA = c['libA'][seed]; lF = c['lib']
    a = RP.log_cp10k(sA[keep], lA[keep]) - RP.log_cp10k(sA[-1:], lA[-1:])
    b = RP.log_cp10k(sF[keep] - sA[keep], lF[keep] - lA[keep]) - RP.log_cp10k(sF[-1:] - sA[-1:], lF[-1:] - lA[-1:])
    w = 1.0 / (1.0 / c['n'][keep] + 1.0 / c['n'][-1])
    return a, b, w


def build_rows(T, seed):
    """Rows (context, target) with half-A / half-B responses for one split seed (separate control halves)."""
    ctx, drug, YA, YB = [], [], [], []
    names = list(T)
    for ci, name in enumerate(names):
        t = T[name]
        if 'parts' in t:
            ix = {x: i for i, x in enumerate(t['targets'])}; G = len(t['desc'])
            a = np.zeros((len(ix), G)); b = np.zeros((len(ix), G)); W = np.zeros(len(ix))
            for p in t['parts']:
                pa, pb, pw = _resp(p, seed); r = [ix[x] for x in p['targets']]
                a[r] += pw[:, None] * pa; b[r] += pw[:, None] * pb; W[r] += pw
            a /= W[:, None]; b /= W[:, None]
        else:
            a, b, _ = _resp(t, seed)
        YA.append(a); YB.append(b); ctx += [ci] * len(t['targets']); drug += t['targets']
    return np.array(ctx), np.array(drug), np.vstack(YA), np.vstack(YB)


# ------------------------------------------------------------------ one configuration

def run_config(label, specs, n_boot=RP.N_BOOT, extras=True, desc_override=None, s3=False):
    t0 = time.time()
    T = context_tables(specs)
    names = list(T)
    X = np.stack([T[n]['desc'] for n in names]) if desc_override is None else desc_override(T)
    d = RP.participation_ratio(X)
    pcs, lam = RP.descriptor_pcs(X)
    k = math.ceil(d)
    S, E = [], []
    keep_rows = None
    for s in range(N_SEEDS):
        ctx, drug, YA, YB = build_rows(T, s)
        if keep_rows is None:
            keep_rows = RP.filter_min_lines(ctx, drug)
        c2, d2, A2, B2 = ctx[keep_rows], drug[keep_rows], YA[keep_rows], YB[keep_rows]
        S.append(RP.halves_stats(A2, B2, c2, d2, pcs, [k]))
        if extras:
            E.append(RP.halves_extra_stats(A2, B2, c2, d2))
        del YA, YB, A2, B2
    agg = RP.aggregate(S, d, k, n_boot=n_boot)
    n_r2 = int(agg[f'n_drugs_k{k}']['estimate'])
    # coverage after the >= 3-line filter
    ctx, drug = c2, d2
    cov = {}
    for ci, n in enumerate(names):
        m = ctx == ci
        tset = set(drug[m]); cells = [T[n]['cells'][T[n]['targets'].index(t)] for t in tset]
        cov[n] = dict(line=LINE_NAME.get(n, n), n_targets_used=int(m.sum()),
                      cells_per_target_median=float(np.median(cells)) if cells else None,
                      cells_per_target_min=float(np.min(cells)) if cells else None,
                      control_cells=T[n]['n_ctrl'], dropped_half_imbalance=T[n]['n_dropped_imbalance'])
    out = dict(label=label, contexts=names, C=len(names), d=d, d_cap=len(names) - 1, eigen_share=(lam / lam.sum()).tolist(),
               d_jackknife=RP.d_jackknife(X), n_targets=int(len(np.unique(drug))), n_rows=int(len(drug)),
               lines_per_target=dict(zip(*[x.tolist() for x in np.unique(np.bincount(np.unique(drug, return_inverse=True)[1]), return_counts=True)])),
               n_genes=int(S[0]['G']), additive_rank=S[0]['additive_rank'], coverage=cov,
               var_i=agg['var_i'], sig2_eps=agg['sig2_eps'], ratio_noise_interaction=agg['ratio'], C_min=agg['C_min'],
               per_seed=[dict(var_i=p['var_i'], sig2_eps=p['sig2_eps'], C_min=p['C_min']) for p in agg['per_seed']],
               k=k, n_targets_R2_estimable=n_r2, bootstrap=agg['bootstrap'])
    if n_r2 >= 1:
        out['R2_adj'] = agg[f'R2_adj_k{k}']; out['C_theory'] = agg['C_theory']
        out['R2_note'] = (f'estimable on {n_r2} targets with n_t - k - 1 >= 1 (k = {k}); residual d.o.f. per target '
                          f'<= C - k - 1 = {len(names) - k - 1}: unstable, descriptive only')
    else:
        out['R2_adj'] = 'not estimable'; out['C_theory'] = 'not estimable'
        out['R2_note'] = f'no target has n_t - k - 1 >= 1 with k = ceil(d) = {k} and C = {len(names)}'
    if extras:
        ES = RP.sum_extra_stats(E)
        se = RP.sensitivity_estimates(ES, X)
        out['S1_S2_point'] = se
        # consistency: primary point from the extra stats equals halves_stats
        out['check_extra_vs_primary'] = dict(var_i=[se['var_i'], agg['var_i']['estimate']],
                                             sig2=[se['sig2_eps'], agg['sig2_eps']['estimate']])
        if s3:
            out['S3_line_bootstrap'] = RP.line_bootstrap(ES, X, log=lambda m: print(m, flush=True))
    out['secs'] = round(time.time() - t0, 1)
    log(event='config', label=label, C=len(names), d=round(d, 3), n_targets=out['n_targets'],
        var_i=agg['var_i']['estimate'], sig2=agg['sig2_eps']['estimate'], C_min=agg['C_min']['estimate'],
        C_min_ci=agg['C_min']['ci95'])
    out['_S'], out['_k'] = S, k
    return out, T


def with_descriptor(label, base_out, T, desc_fn, n_boot=RP.N_BOOT):
    """Same rows / halves statistics as base_out; only the descriptor (hence d) changes -> C*_min re-aggregated."""
    X = desc_fn(T)
    d = RP.participation_ratio(X)
    agg = RP.aggregate(base_out['_S'], d, base_out['_k'], n_boot=n_boot)
    out = dict(label=label, C=base_out['C'], contexts=base_out['contexts'], d=d, d_cap=base_out['C'] - 1,
               n_features=int(X.shape[1]), var_i=agg['var_i'], sig2_eps=agg['sig2_eps'],
               ratio_noise_interaction=agg['ratio'], C_min=agg['C_min'],
               note='rows, targets, halves and Var(i), sig2_eps identical to the primary; only d recomputed')
    log(event='config', label=label, d=round(d, 3), C_min=agg['C_min']['estimate'], C_min_ci=agg['C_min']['ci95'])
    return out


def ctrl_pc_descriptor(kpc):
    """Descriptor restricted to the top-k PCs of the pooled control cells (<= 2,000 per line; log1p CP10K)."""
    def f(T):
        samp = []
        for n in T:
            z = load(n)
            s = z['ctrl_sample'].astype(np.float64); L = z['ctrl_sample_lib']
            ok = L > 0
            samp.append(np.log1p(s[ok] * 1e4 / L[ok, None]))
        Z = np.vstack(samp); mu = Z.mean(0)
        _, _, Vt = np.linalg.svd(Z - mu, full_matrices=False)
        D = np.stack([T[n]['desc'] for n in T])
        return (D - mu) @ Vt[:kpc].T
    return f


def verdict(cmin, cpub):
    est, lo = cmin['estimate'], cmin['ci95'][0]
    if est <= cpub:
        return 'rejected'
    if lo > cpub:
        return 'confirmed'
    return 'inconclusive'


def main():
    gene_json = RES / 'rq2_params_gene.json'
    G = json.loads(gene_json.read_text())
    assert 'count2' in G, 'Count 2 must be written before any C*_min (A2)'
    c1 = int(json.loads((RES / 'families.json').read_text())['C_public']); c2 = int(G['count2']['C_public_count2'])
    plan = json.loads((CACHE / 'plan.json').read_text())
    fl = feng_lines()
    base = {s: ('src', s) for s in PRIMARY if s != 'FENG'}
    base['FENG'] = ('merge', [('feng_line', l) for l in fl])     # pool = within-line responses averaged (G1b)
    base = {s: base[s] for s in PRIMARY}
    res = {}
    res['primary'], T = run_config('primary (Count-1 contexts, H1 = H1GW)', base, s3=True)
    sens = {}
    for kpc in (8, 24, 50):
        sens[f'ctrlPC_top{kpc}'] = with_descriptor(f'descriptor = projection on the top {kpc} pooled-control-cell PCs',
                                                   res['primary'], T, ctrl_pc_descriptor(kpc))
    sens['no_KOLF'], _ = run_config('KOLF excluded', {k: v for k, v in base.items() if k != 'KOLF'}, extras=False)
    hm = dict(base); hm['H1GW'] = ('merge', [('src', 'H1GW'), ('src', 'VCC_H1')])
    sens['H1GW_VCCH1_one_line'], _ = run_config('H1GW + VCC2025-H1 as one H1 line (responses within dataset, averaged per target)', hm, extras=False)
    hs = dict(base); hs['VCC_H1'] = ('src', 'VCC_H1')
    sens['VCCH1_separate_context'], _ = run_config('VCC2025-H1 as a 12th context (datasets-as-contexts reading)', hs, extras=False)
    fs = dict(base); fs['FENG'] = ('feng_pool', fl)
    sens['FENG_pooled_pseudobulk'], _ = run_config('HipSci pool as one summed pseudobulk (registry-style; line mix of targets != of controls)', fs, extras=False)
    sens['no_CD4'], _ = run_config('CD4 excluded (unit-level halves)', {k: v for k, v in base.items() if k != 'CD4'}, extras=False)
    res['sensitivity'] = sens
    for o in [res['primary'], *sens.values()]:
        o.pop('_S', None); o.pop('_k', None)
    # secondary: Count-2 donor-line contexts (pooled HipSci context replaced by the qualifying non-kolf lines)
    q = G['count2']['qualifying_nonkolf_lines']
    sec = {k: v for k, v in base.items() if k != 'FENG'}
    for l in q:
        sec['FENG:' + l] = ('feng_line', l)
    res['secondary_count2_contexts'], _ = run_config(f'Count-2 contexts ({len(sec)}): 10 lines + {len(q)} HipSci donor lines', sec, s3=True)
    res['secondary_count2_contexts'].pop('_S', None); res['secondary_count2_contexts'].pop('_k', None)
    # H2c
    P = res['primary']['C_min']
    v1, v2 = verdict(P, c1), verdict(P, c2)
    if v1 == 'confirmed' and v2 == 'confirmed':
        h2c = 'confirmed'
    elif v1 == 'confirmed':
        h2c = 'holds for distinct cell types, not for donor lines'
    else:
        h2c = v1
    h = dict(hypothesis='H2c: C*_min = d + 1 + d sig2_eps / Var(i) on public genetic screens exceeds C_public',
             rule='rejected if point C*_min <= C_public; confirmed if the 95 % lower bound > C_public; else inconclusive; '
                  'A2: confirmed only if it holds under both counts; Count 1 only -> "holds for distinct cell types, not for donor lines"',
             C_public_count1=c1, C_public_count2=c2,
             C_min=P['estimate'], C_min_ci95=P['ci95'], d=res['primary']['d'], d_cap=res['primary']['d_cap'],
             var_i=res['primary']['var_i']['estimate'], sig2_eps=res['primary']['sig2_eps']['estimate'],
             ratio=res['primary']['ratio_noise_interaction']['estimate'],
             verdict_count1=v1, verdict_count2=v2, verdict=h2c,
             secondary_count2_contexts=dict(C=res['secondary_count2_contexts']['C'], d=res['secondary_count2_contexts']['d'],
                                            C_min=res['secondary_count2_contexts']['C_min'],
                                            verdict_vs_count2_descriptive=verdict(res['secondary_count2_contexts']['C_min'], c2)),
             sensitivities={k: dict(C=v['C'], d=v['d'], C_min=v['C_min']['estimate'], ci95=v['C_min']['ci95'],
                                    verdict_count1=verdict(v['C_min'], c1), verdict_count2=verdict(v['C_min'], c2))
                            for k, v in sens.items()},
             S1_S2_S3=dict(point=res['primary']['S1_S2_point'], line_bootstrap=res['primary'].get('S3_line_bootstrap'),
                           note='secondary (pre-specified 2026-11-05); does not change the verdict rule'),
             computed=datetime.date.today().isoformat())
    (RES / 'rq2_h2c.json').write_text(json.dumps(RP._jsonable(h), indent=1))
    G.update(RP._jsonable(dict(plan=dict(n_universe=plan['n_universe'], n_axis=plan['n_axis'],
                                         universe_rule='targets with >= 20 cells in >= 1 of the RPE1/HepG2/Jurkat essential screens'),
                               **res)))
    G['runtime_s'] = round(time.time() - T0, 1)
    gene_json.write_text(json.dumps(G, indent=1))
    log(event='done', verdict=h2c, v1=v1, v2=v2)


if __name__ == '__main__':
    main()
