"""Official vcc2026 preset scorer for the Feng+CD4 context experiment (copy of dcp_preset_score with SRC redirected; EXT_SRC env overrides). Part 17: score the drug-context-pretraining candidates with the OFFICIAL cell-eval2 vcc2026
preset, not with PDS.

This reuses the real bundles already built in reports/preset_rescore_2026-09-20/bundles (the
`cell-eval2 baseline` 0-end and the `prep-real-bundle` 1-end for each of the five lines), so the
normalisation anchors are identical to the previous official re-scoring and the two experiments
are directly comparable.

Arms are named by prediction stem inside reports/drug_context_pretrain_2026-09-20/<line>/.
Every run asserts that the preset's own pds_cosine equals the harness PDS stored in scores.json,
so a silent mismatch in target ordering or reconstruction fails loudly.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse, gc, json, os, subprocess, sys, time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

_SCRIPTS = str(Path(__file__).resolve().parent)
sys.path.insert(0, _SCRIPTS)
import context_forcing_score as C
# importing context_forcing_score removes the scripts dir from sys.path again
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)
import simple_residual_full_score_batch as S

E = C.E
SRC = C.ROOT / (os.environ.get('EXT_SRC') or 'reports/ext_contexts_2026-09-21')
ROOT = SRC / 'preset'
BUNDLES = C.ROOT / 'reports/preset_rescore_2026-09-20/bundles'
CE = Path(sys.executable).with_name('cell-eval2')
LINES = ('K562', 'RPE1', 'HepG2', 'Jurkat', 'KOLF')


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    S._write_json(path, value)


def annotated(matrix, labels, genes, prefix):
    return ad.AnnData(X=matrix,
                      obs=pd.DataFrame({'target': pd.Categorical(labels)},
                                       index=[f'{prefix}_{i}' for i in range(len(labels))]),
                      var=pd.DataFrame(index=list(genes)))


def build_prediction(data, scaled, targets, seed, cells=400):
    labels, bulk, diagnostics = E._reconstructed_bulk(data.query_controls, scaled, targets, cells, seed)
    obj = SimpleNamespace(controls=data.query_controls, cells=cells, targets=targets, genes=data.genes)
    matrix = E.MockEvaluator._prediction_matrix(obj, scaled, seed)
    result = annotated(sparse.vstack([matrix, data.query_controls], format='csr'),
                       [*np.repeat(targets, cells), *[E.CONTROL] * data.query_controls.shape[0]],
                       data.genes, 'pred')
    actual = E.pseudobulk_bulk_lognorm(result, 'target', bulk_target_sum=E.BULK_TARGET_SUM)
    np.testing.assert_array_equal(actual[0], labels)
    err = float(np.max(np.abs(actual[1] - bulk)))
    if err > 1e-12:
        raise AssertionError(f'reconstruction bulk mismatch {err}')
    diagnostics['bulk_max_abs_error'] = err
    return result, diagnostics


def stored_pds(line, tag, seed):
    # the gated/ensemble arms are scored by dcp_gate.py into a separate file
    rows = json.loads((SRC / 'scores.json').read_text())
    gp = SRC / 'gate_scores.json'
    if gp.exists():
        rows = rows + json.loads(gp.read_text())
    hit = [r for r in rows if r['line'] == line and r['tag'] == tag and r['seed'] == seed]
    if len(hit) != 1:
        raise ValueError(f'expected one stored score for {line}/{tag}/{seed}, got {len(hit)}')
    return float(hit[0]['pds'])


def real_path(line):
    p = C.ROOT / 'reports/b19_sweep2_2026-09-18/full_metrics' / line / 'real_disjoint.h5ad'
    if not p.exists():
        raise FileNotFoundError(f'missing frozen real panel for {line}: {p}')
    return p


def phase_run(line, arms, seeds):
    bundle_dir = BUNDLES / line
    if not (bundle_dir / 'bundle').exists():
        raise FileNotFoundError(f'no prebuilt real bundle for {line}: {bundle_dir}. '
                                f'Run preset_rescore.py --phase bundle first.')
    worker = C.ContextForcingScoreWorker(ROOT)
    data = worker._load_line(line)
    targets0 = json.loads(
        (C.ROOT / 'reports/context_forcing_2026-09-16/protocol.json').read_text())['final_targets'][line]
    try:
        for tag in arms:
            delta_file = SRC / line / f'{tag}.npy'
            if not delta_file.exists():
                print(json.dumps(dict(event='missing', line=line, arm=tag)), flush=True)
                continue
            scaled, audit = C.prepare_delta(np.load(delta_file), targets0, data, 3.)
            targets = audit['sorted_targets']
            for seed in seeds:
                out = ROOT / 'runs' / line / tag / f'seed{seed}'
                if (out / 'receipt.json').exists():
                    continue
                out.mkdir(parents=True, exist_ok=True)
                t0 = time.monotonic()
                pred_path = out / 'prediction.h5ad'
                if not pred_path.exists():
                    pred, diag = build_prediction(data, scaled, targets, seed)
                    pred.write_h5ad(pred_path, compression='gzip')
                    save(out / 'reconstruction.json', diag)
                    del pred; gc.collect()
                cmd = [str(CE), 'run', '-ap', str(pred_path), '-ar', str(real_path(line)),
                       '--preset', 'vcc2026', '--pert-col', 'target', '--control', E.CONTROL,
                       '--input-type', 'counts', '-o', str(out / 'cell_eval2'),
                       '--cache-real', str(bundle_dir / 'cache_real'),
                       '--set', 'de.backend=pdex', '--set', 'num_threads=2',
                       '--set', 'gather_threads=2']
                if not (out / 'cell_eval2/agg_results.csv').exists():
                    with (out / 'scorer.log').open('a') as g:
                        subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT,
                                       env=S.thread_environment(), check=True)
                agg = pd.read_csv(out / 'cell_eval2/agg_results.csv', index_col=0)
                raw = {k: float(agg.loc['mean', k]) for k in S.RAW_SIX}
                want = stored_pds(line, tag, seed)
                err = abs(raw['pds_cosine'] - want)
                if not np.isfinite(err) or err > 1e-9:
                    raise RuntimeError(f'PDS parity failed {line}/{tag}/seed{seed}: '
                                       f'{raw["pds_cosine"]} vs stored {want} (err {err})')
                save(out / 'receipt.json', dict(status='complete', line=line, arm=tag, seed=seed,
                                                cells=400, scale=3, raw_six=raw, stored_pds=want,
                                                pds_parity_error=err,
                                                delta_sha256=S.sha256(delta_file),
                                                seconds=round(time.monotonic() - t0, 1)))
                try:
                    pred_path.unlink()
                except FileNotFoundError:
                    pass
                print(json.dumps(dict(event='run_done', line=line, arm=tag, seed=seed,
                                      pds=round(raw['pds_cosine'], 6),
                                      seconds=round(time.monotonic() - t0, 1))), flush=True)
            del scaled; gc.collect()
    finally:
        worker.close()


def phase_score(line, arms, seeds):
    bundle = BUNDLES / line / 'bundle'
    for tag in arms:
        for seed in seeds:
            out = ROOT / 'runs' / line / tag / f'seed{seed}'
            agg = out / 'cell_eval2/agg_results.csv'
            if not agg.exists() or (out / 'scored.csv').exists():
                continue
            cmd = [str(CE), 'score', '--user-agg', str(agg), '--real-bundle', str(bundle),
                   '-o', str(out / 'scored.csv')]
            with (out / 'score.log').open('a') as g:
                subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT,
                               env=S.thread_environment(), check=True)
            print(json.dumps(dict(event='scored', line=line, arm=tag, seed=seed)), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=['run', 'score', 'all'], default='all')
    ap.add_argument('--lines', default=','.join(LINES))
    ap.add_argument('--arms', default='')
    ap.add_argument('--arm-map', default='', help='JSON from dcp_select.py: per-line tag lists')
    ap.add_argument('--seeds', default='0,1,2')
    a = ap.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    seeds = [int(s) for s in a.seeds.split(',')]
    amap = json.loads(Path(a.arm_map).read_text())['arm_map'] if a.arm_map else None
    for line in a.lines.split(','):
        arms = amap[line] if amap else a.arms.split(',')
        if a.phase in ('run', 'all'):
            phase_run(line, arms, seeds)
        if a.phase in ('score', 'all'):
            phase_score(line, arms, seeds)
    print(json.dumps(dict(event='phase_complete', phase=a.phase)), flush=True)


if __name__ == '__main__':
    main()
