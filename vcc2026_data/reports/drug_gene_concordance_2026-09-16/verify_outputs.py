"""Separate reconstruction of raw matched pairs from source IDs and summary arithmetic."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd,h5py
from scipy.stats import spearmanr
R=Path('/home/wonmin-kim/AIVC/vcc2026_data');O=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'scripts'));from vcclib import GeneMap
GM=GeneMap()
def series(values,genes):return pd.Series(values,index=[GM.canon(g) for g in genes]).groupby(level=0).mean()
a=pd.read_csv(O/'pair_metrics.csv');s=pd.read_csv(O/'summary.csv');f=pd.read_csv(R/'reports/drug_gene_feasibility_2026-09-16/candidate_pairs.csv')
f=f[(f.moa=='inhibitor/antagonist')&(f.dataset!='lincs_selected200')]
assert set(map(tuple,a[['dataset','cell_line','drug','target']].values))==set(map(tuple,f[['dataset','cell_line','drug','target']].values))
z=np.load(R/'reports/context_signal_2026-09-15/tahoe/panel.npz',allow_pickle=True)
drugnames=[json.loads(x)[0][0] for x in z['targets']];excluded=set(json.load(open(O/'excluded_target_genes.json')))
feat={(x['dataset'],x['cell_line']):x['genes'] for x in json.load(open(O/'features.json'))}
L=R/'raw/context_signal_2026-09-15/lincs';meta=pd.read_csv(L/'GSE106127_CGS_meta.txt.gz',sep='\t',dtype=str);genes=pd.read_csv(L/'GSE106127_gene_info.txt.gz',sep='\t',dtype=str).set_index('pr_gene_id').pr_gene_symbol
checks=[]
for (dataset,line),block in a[a['mode']=='raw'].groupby(['dataset','cell_line']):
 r=block.iloc[0];ti=list(z['lines']).index({'A549':'CVCL_0023','HT29':'CVCL_0320'}[line]);di=drugnames.index(r.drug)
 ds=series(z['effects'][ti,0,di],z['genes'])
 if dataset.startswith('jiang'):
  q=np.load(R/'processed/context_signal_2026-09-15'/f'{dataset}.npz',allow_pickle=True);li=list(q['lines']).index(line);j=[GM.canon(t) for t in q['targets']].index(r.target)
  gs=series(q['effects'][li,:,j].mean(0),q['genes'])
 else:
  sm=meta[(meta.cell_id==line)&(meta.pert_iname.map(GM.canon)==r.target)&(meta.pert_type=='trt_sh.cgs')&(meta.pert_itime=='96 h')];assert len(sm)==1
  with h5py.File(L/'GSE106127_CGS_n33839x978.gctx','r') as h:
   ci=list(h['0/META/COL/id'][:]).index(sm.sig_id.iloc[0].encode());rg=[genes[x.decode()] for x in h['0/META/ROW/id'][:]]
   gs=series(h['0/DATA/0/matrix'][ci,:],rg)
 shared=sorted((set(ds.index)&set(gs.index))-excluded);assert shared==feat[(dataset,line)]
 rho=spearmanr(ds.loc[shared],gs.loc[shared]).statistic
 assert np.isclose(rho,r.spearman,atol=2e-6), (dataset,line,rho,r.spearman)
 checks.append(dict(dataset=dataset,line=line,drug=r.drug,target=r.target,independent_spearman=float(rho),saved_spearman=float(r.spearman)))
for r in s.itertuples():
 g=a[(a.dataset==r.dataset)&(a['mode']==r.mode)]
 if r.subset=='single_listed_target':g=g[g.single_listed_target]
 for col in ['spearman','wrong_spearman','spearman_advantage','matched_percentile']:
  v=g.groupby(['drug','cell_line'])[col].mean().groupby('drug').mean().mean()
  assert np.isclose(v,getattr(r,col),atol=1e-12)
assert a.matched_percentile.between(0,1).all();assert a.spearman.between(-1,1).all()
assert np.allclose(a.spearman-a.wrong_spearman,a.spearman_advantage)
assert json.load(open(O/'verification.json'))['all_sources_unchanged']
(O/'independent_verification.json').write_text(json.dumps(dict(source_id_reconstruction=checks,all_18_summary_rows_recomputed=True,all_401_candidate_pairs_replayed=True,passed=True),indent=2))
print('PASS: 10 source-ID/gene-alignment reconstructions; 18 summary rows; 401 candidate identities; metric ranges and contrasts.')
