"""Arm C: installed six-metric evaluation of sweep-2 outer seed-0 deltas with the frozen disjoint controls. Diagnostic only."""
import argparse,gc,json,subprocess,sys
from pathlib import Path
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import context_forcing_score as C
sys.path.insert(0,str(Path(__file__).resolve().parent))
import simple_residual_full_score_batch as S
E=C.E
ROOT=C.ROOT/'reports/b19_sweep2_2026-09-18'
TAGS=('b19_reference','selected','refit','ensemble_k3','ensemble_k5')

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    S._write_json(path,value)

def annotated(matrix,labels,genes,prefix):
    return ad.AnnData(X=matrix,obs=pd.DataFrame({'target':pd.Categorical(labels)},index=[f'{prefix}_{i}' for i in range(len(labels))]),var=pd.DataFrame(index=list(genes)))

def prediction(data,scaled,targets,cells=400,seed=0):
    labels,bulk,diagnostics=E._reconstructed_bulk(data.query_controls,scaled,targets,cells,seed)
    obj=SimpleNamespace(controls=data.query_controls,cells=cells,targets=targets,genes=data.genes)
    matrix=E.MockEvaluator._prediction_matrix(obj,scaled,seed)
    result=annotated(sparse.vstack([matrix,data.query_controls],format='csr'),[*np.repeat(targets,cells),*[E.CONTROL]*data.query_controls.shape[0]],data.genes,'pred')
    actual=E.pseudobulk_bulk_lognorm(result,'target',bulk_target_sum=E.BULK_TARGET_SUM)
    np.testing.assert_array_equal(actual[0],labels)
    error=float(np.max(np.abs(actual[1]-bulk)))
    if error>1e-12:raise AssertionError(f'Reconstruction bulk mismatch {error}')
    diagnostics['bulk_max_abs_error']=error
    return result,diagnostics

def replace_real_controls(real,data,targets):
    if tuple(real.var_names)!=tuple(data.genes):raise ValueError('Real gene axis differs')
    labels=real.obs.target.astype(str).to_numpy();keep=np.isin(labels,targets)
    return annotated(sparse.vstack([sparse.csr_matrix(real.X[keep]),data.reference_controls],format='csr'),[*labels[keep],*[E.CONTROL]*data.reference_controls.shape[0]],data.genes,'real')

def delta_targets(line,delta):
    manifest=delta.with_suffix('.manifest.json')
    if manifest.exists():return json.loads(manifest.read_text())['targets']
    return json.loads((C.ROOT/'reports/context_forcing_2026-09-16/protocol.json').read_text())['final_targets'][line]

def run(line,tag):
    rows=json.loads((ROOT/'outer'/line/'scores.json').read_text());row=next(r for r in rows if r['seed']==0 and r['tag']==tag)
    out=ROOT/'full_metrics'/line/tag;out.mkdir(parents=True,exist_ok=True)
    if not row['valid']:
        save(out/'receipt.json',dict(status='ineligible_invalid_counts',line=line,tag=tag));return
    delta=Path(row['prediction_path']);targets=delta_targets(line,delta)
    request=dict(line=line,tag=tag,delta=str(delta),delta_sha256=S.sha256(delta),expected_pds=row['pds'],seed=0,cells=400,scale=3,adapter_sha256=S.sha256(Path(__file__)))
    if (out/'receipt.json').exists():
        receipt=json.loads((out/'receipt.json').read_text())
        if receipt['status']!='complete' or any(receipt[k]!=v for k,v in request.items()):raise ValueError('Completed full-metric request differs')
        for path,key in [(out/'prediction.h5ad','prediction_sha256'),(ROOT/'full_metrics'/line/'real_disjoint.h5ad','real_sha256'),(out/'cell_eval2/agg_results.csv','aggregate_sha256')]:
            if S.sha256(path)!=receipt[key]:raise ValueError('Completed full-metric artifact changed')
        S.require_pds_parity(receipt['raw_six']['pds_cosine'],row['pds'])
        return
    if (out/'request.json').exists() and json.loads((out/'request.json').read_text())!=request:raise ValueError('Full-metric request changed')
    save(out/'request.json',request)
    worker=C.ContextForcingScoreWorker(ROOT);data=worker._load_line(line)
    scaled,audit=C.prepare_delta(np.load(delta),targets,data,3.);targets=audit['sorted_targets']
    real_path=ROOT/'full_metrics'/line/'real_disjoint.h5ad'
    if not real_path.exists():
        original=ad.read_h5ad(C.KOLF_REAL if line=='KOLF' else C.RN_MOCK/line/'panel300/real.h5ad')
        real=replace_real_controls(original,data,targets);del original
        truth=E.pseudobulk_bulk_lognorm(real,'target',bulk_target_sum=E.BULK_TARGET_SUM)
        expected=C._subset_truth(data.real_bulk,targets,data.reference_controls)
        np.testing.assert_array_equal(truth[0],expected[0]);np.testing.assert_allclose(truth[1],expected[1],rtol=0,atol=1e-12)
        real.write_h5ad(real_path,compression='gzip');del real;gc.collect()
    pred_path=out/'prediction.h5ad'
    if not pred_path.exists():
        pred,diagnostics=prediction(data,scaled,targets)
        pred.write_h5ad(pred_path,compression='gzip');save(out/'reconstruction.json',diagnostics);del pred;gc.collect()
    worker.close();del data;gc.collect()
    command=S.score_command(Path(sys.executable).parent/'cell-eval2',pred_path,real_path,out/'cell_eval2')
    if not (out/'cell_eval2/agg_results.csv').exists():
        with (out/'scorer.log').open('a') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=S.thread_environment(),check=True)
    agg=pd.read_csv(out/'cell_eval2/agg_results.csv',index_col=0);metrics={k:float(agg.loc['mean',k]) for k in S.RAW_SIX}
    error=S.require_pds_parity(metrics['pds_cosine'],row['pds'])
    missing=[k for k,v in metrics.items() if not np.isfinite(v)]
    metrics={k:(v if np.isfinite(v) else None) for k,v in metrics.items()}
    save(out/'receipt.json',dict(status='complete',**request,raw_six=metrics,missing_metrics=missing,pds_parity_error=error,command=command,prediction_sha256=S.sha256(pred_path),real_sha256=S.sha256(real_path),aggregate_sha256=S.sha256(out/'cell_eval2/agg_results.csv'),note='Local raw metrics; not calibrated official overall. Seed 0 diagnostic only; never a selector.'))
    print(json.dumps(dict(event='full_metrics_complete',line=line,tag=tag,metrics=metrics)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--line',required=True);p.add_argument('--tag',choices=TAGS,required=True);a=p.parse_args();run(a.line,a.tag)
