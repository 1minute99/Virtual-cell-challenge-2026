"""Part 3 of the Feng+CD4 experiment: the exact current baselines with and without the new
contexts, strict leave-one-scored-biology-out, scored on the five-line harness (PDS) here and by
the official vcc2026 preset downstream (`ext_preset_score.py`).

  mean_cur   B.empirical over the six frozen effects stores (one vote per store)
  mean_v2    same, with pooled CD4 replaced by the donor x state store and the Feng line store
             added -- still one vote per STORE (7 voters): "finer data, same weighting"
  mean_v1    each new context is its own voter (5 old stores + 12 CD4 + 19 Feng = 36 voters):
             "more contexts"
  refit_cur  LeanSourceBaseline rank 32 on the frozen table (== stored refit_alpha0, asserted)
  refit_v1/v2  the same fit on the extended tables (`ext_assemble.py`)

Held-out line L: every row whose canonical biology is L is excluded from every fit and vote; on
the KOLF fold that includes the kolf_2 / kolf_3 HipSci lines (labelled KOLF).
"""
from pathlib import Path
import argparse, gc, json, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import b19_optimization_run as R
import b19_optimization_data as A
import b19_optimization_report as P
import context_forcing_data as F
import simple_residual_data as D
import simple_residual_baselines as B
import vcclib as V
import fullaxis_lean_baseline as LB

ROOT = D.ROOT / 'reports/ext_contexts_2026-09-21'
DCP = D.ROOT / 'reports/drug_context_pretrain_2026-09-20'
LINES = tuple(R.LINES)
NEW = dict(cd4_Rest='cd4_target_condition_donor_all_Rest', cd4_Stim8hr='cd4_target_condition_donor_all_Stim8hr',
           cd4_Stim48hr='cd4_target_condition_donor_all_Stim48hr', feng='feng_targetedscreen_target_line',
           feng_fit='feng_genomewidescreen_fitnessgenes_target')
KOLF_LINES = ('kolf_2', 'kolf_3')
CONFLICT_GENE = {'KOLF': ['FENGLKOLF2', 'FENGLKOLF3']}   # v1 labels of the kolf-background lines


def new_meta(store):
    p = pd.read_csv(D.ROOT / 'effects' / NEW[store] / 'pairs.csv', low_memory=False)
    gm = V.GeneMap()
    p['canonical'] = [gm.canon(x) for x in p.target_symbol.astype(str)]
    if store.startswith('cd4'):
        p['bio'] = 'CD4'
        p['ctx'] = [f'cd4donor/{d}/{s}' for d, s in zip(p.donor, p.stimulation)]
    elif store == 'feng_fit':
        p['bio'] = 'FENG'; p['ctx'] = 'fengpooled/fitness'
    else:
        p['bio'] = ['KOLF' if str(d) in KOLF_LINES else 'FENG' for d in p.donor]
        p['ctx'] = [f'fengline/{d}' for d in p.donor]
    return p


def voter(p, delta, ids, canonical):
    """weighted per-target mean over the rows `ids`, exactly as B.empirical does per store."""
    dest = np.full((len(canonical), delta.shape[1]), np.nan, np.float32)
    if not len(ids):
        return dest
    sub = p.iloc[ids]; val = delta[ids]
    np_ = sub.n_cells_pert.to_numpy(float); nc = sub.n_cells_control.to_numpy(float)
    w = 1 / (1 / np.maximum(np_, 1e-9) + 1 / np.maximum(nc, 1e-9))
    own_all = sub.canonical.to_numpy()
    for j, t in enumerate(canonical):
        own = own_all == t
        if not own.any():
            continue
        y = val[own]; ok = np.isfinite(y); weights = ok * w[own, None]
        den = weights.sum(0); num = (np.nan_to_num(y) * weights).sum(0)
        dest[j] = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
    return dest


def empirical(targets, excluded_bios, variant):
    banned = {D.canonical_bio(x) for x in excluded_bios}
    gm = V.GeneMap()
    canonical = [gm.canon(t) for t in targets]
    stack = []
    old = [s for s in B.STORES if not (variant != 'cur' and s == 'cd4_target_condition_all')]
    for store in old:
        p = B.metadata(store)
        ids = np.flatnonzero(~p.bio.isin(banned) & p.canonical.isin(canonical))
        z = np.load(D.ROOT / 'effects' / store / 'delta.npy', mmap_mode='r')
        stack.append(voter(p, z, ids, canonical))
    if variant != 'cur':
        cd4_parts = []
        for store in NEW:
            p = new_meta(store)
            z = np.load(D.ROOT / 'effects' / NEW[store] / 'delta.npy', mmap_mode='r')
            keep = ~p.bio.isin(banned) & p.canonical.isin(canonical)
            if variant == 'v2':
                ids = np.flatnonzero(keep)
                if not len(ids):
                    continue
                v = voter(p, z, ids, canonical)
                # v2 keeps the CURRENT weighting: CD4 is ONE voter (its three condition stores are
                # averaged into one profile, as the pooled cd4_target_condition_all store was)
                (cd4_parts if store.startswith('cd4') else stack).append(v)
            else:
                for c in sorted(set(p.ctx)):
                    ids = np.flatnonzero(keep & (p.ctx == c))
                    if len(ids):
                        stack.append(voter(p, z, ids, canonical))
        if cd4_parts:
            a3 = np.stack(cd4_parts); d3 = np.isfinite(a3).sum(0)
            stack.append(np.divide(np.nansum(a3, 0), d3, out=np.full(a3.shape[1:], np.nan, np.float32),
                                   where=d3 > 0).astype(np.float32))
    a = np.stack(stack)
    den = np.isfinite(a).sum(0)
    out = np.divide(np.nansum(a, axis=0), den, out=np.zeros(a.shape[1:], np.float32),
                    where=den > 0).astype(np.float32)
    return out, len(stack)


def refit(table, line, targets, feats, canon, rank=32):
    ids = F.permitted_ids(table, [line])
    conf = CONFLICT_GENE.get(D.canonical_bio(line), [])
    if conf:
        ids = ids[~np.isin(np.asarray(table.bios)[ids].astype(str), conf)]
    held = {D.canonical_bio(x) for x in np.asarray(table.bios)[ids]}
    assert D.canonical_bio(line) not in held and not (held & set(conf))
    y = np.asarray(table.y[ids], np.float32)
    meta = {k: np.asarray(getattr(table, k))[ids] for k in
            ('targets', 'bios', 'modality', 'features', 'weights', 'sources')}
    sub = D.Table(y=y, controls=np.zeros((len(ids), 1), np.float32), genes=list(table.genes), **meta)
    fitted = LB.LeanSourceBaseline(sub, np.arange(len(ids)), rank=rank)
    common, seen = fitted.predict(canon, feats)
    del sub, fitted, y, meta; gc.collect()
    return common.astype(np.float32), int(len(ids)), len(held), bool(seen.all())


def score_one(worker, targets, line, path, seed):
    out = Path(path).with_name(Path(path).stem + f'_score_seed{seed}.json')
    res = worker.request(dict(line=line, targets=targets, delta_path=str(path),
                              output_path=str(out), seed=seed, cells=400, scale=3., null_repeats=0))
    return dict(pds=res['pds'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lines', default=','.join(LINES))
    ap.add_argument('--arms', default='mean_cur,mean_v2,mean_v1,refit_cur,refit_v2,refit_v1')
    ap.add_argument('--seeds', default='0,1,2')
    ap.add_argument('--no-score', action='store_true')
    a = ap.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True); R.OUT = ROOT
    seeds = [int(s) for s in a.seeds.split(',')]
    base = A.Data()
    indices = {g: i for i, g in enumerate(V.official_features())}
    columns = np.array([indices[g] for g in base.table.genes])
    tables = {'cur': base.table}
    for v in ('v1', 'v2'):
        if (ROOT / f'table_{v}' / 'y.npy').exists():
            tables[v] = D.Table.load(ROOT / f'table_{v}')
    rows_path = ROOT / 'scores.json'
    rows = json.loads(rows_path.read_text()) if rows_path.exists() else []
    diag_path = ROOT / 'baseline_diagnostics.json'
    diags = json.loads(diag_path.read_text()) if diag_path.exists() else []
    worker = None if a.no_score else R.Worker('ext')
    try:
        for line in a.lines.split(','):
            targets = base.targets[line]
            canon = [D.canonical_target(t) for t in targets]
            feats = base.features(targets)
            (ROOT / line).mkdir(exist_ok=True)
            emp = {}
            for arm in a.arms.split(','):
                kind, variant = arm.split('_')
                t0 = time.monotonic()
                if variant not in emp:
                    emp[variant], nv = empirical(targets, [line], variant)
                    if variant == 'cur':
                        ref = B.empirical(targets, excluded_bios=[line])
                        np.testing.assert_allclose(emp['cur'], ref, atol=1e-6)
                    d = dict(line=line, arm=f'mean_{variant}', voters=nv,
                             norm=float(np.linalg.norm(emp[variant])))
                else:
                    d = dict(line=line, arm=f'mean_{variant}')
                if kind == 'mean':
                    pred = emp[variant]
                else:
                    if variant not in tables:
                        print(json.dumps(dict(event='skip', line=line, arm=arm, reason='no table')), flush=True)
                        continue
                    common, nrows, nbios, seen = refit(tables[variant], line, targets, feats, canon)
                    pred = emp[variant].copy(); pred[:, columns] = common
                    d = dict(line=line, arm=arm, rows=nrows, bios=nbios, all_seen=seen,
                             norm=float(np.linalg.norm(pred)))
                    if variant == 'cur':
                        prev = np.load(DCP / line / 'refit_alpha0.npy')
                        d['max_abs_diff_vs_dcp_refit'] = float(np.abs(pred - prev).max())
                if not np.isfinite(pred).all() or pred.shape != (len(targets), 18533):
                    raise ValueError(f'invalid delta {line}/{arm}')
                path = ROOT / line / (f'{arm}_check.npy' if arm == 'refit_cur' else f'{arm}.npy')
                np.save(path, pred.astype(np.float32))
                d['seconds_build'] = round(time.monotonic() - t0, 1)
                if worker is not None and arm != 'refit_cur':
                    for seed in seeds:
                        if any(r['line'] == line and r['tag'] == arm and r['seed'] == seed for r in rows):
                            continue
                        rows.append(dict(line=line, tag=arm, alpha=0.0, seed=seed,
                                         **score_one(worker, targets, line, path, seed)))
                        P.write(rows_path, rows)
                    d['pds'] = float(np.mean([r['pds'] for r in rows if r['line'] == line and r['tag'] == arm]))
                elif arm == 'refit_cur':
                    d['note'] = 'check only; official refit_cur = bit-identical dcp refit_alpha0 copied by ext_chain.sh'
                diags.append(d); P.write(diag_path, diags)
                print(json.dumps(dict(event='arm', **d)), flush=True)
            del emp; gc.collect()
    finally:
        if worker is not None:
            worker.close()
    print(json.dumps(dict(event='done')), flush=True)


if __name__ == '__main__':
    main()
