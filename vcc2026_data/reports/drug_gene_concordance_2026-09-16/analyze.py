import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import sys,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy.stats import rankdata,spearmanr,pearsonr
ROOT=Path('/home/wonmin-kim/AIVC/vcc2026_data');OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from vcclib import GeneMap,HGNC_TSV
GM=GeneMap();canon=lambda x:GM.canon(str(x).strip())
F=ROOT/'reports/drug_gene_feasibility_2026-09-16';T=ROOT/'reports/context_signal_2026-09-15/tahoe';L=ROOT/'raw/context_signal_2026-09-15/lincs';P=ROOT/'processed/context_signal_2026-09-15'
inputs={}
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def track(p):
 p=Path(p);inputs[str(p)]=digest(p);return p
def collapse(x,genes):
 names=[canon(g) for g in genes];out=sorted(set(names));indices={}
 for i,g in enumerate(names):indices.setdefault(g,[]).append(i)
 result=x[..., [indices[g][0] for g in out]].copy()
 for j,g in enumerate(out):
  if len(indices[g])>1:result[...,j]=x[...,indices[g]].mean(-1)
 return result,np.array(out)
def unit(x,center=False,rank=False):
 x=np.asarray(x,dtype=float)
 if rank:x=rankdata(x,axis=-1,method='average');center=True
 if center:x=x-x.mean(-1,keepdims=True)
 n=np.linalg.norm(x,axis=-1,keepdims=True)
 assert (n>0).all(), 'undefined correlation: zero vector'
 return x/n
track(HGNC_TSV);track(OUT/'ANALYSIS_PLAN.md')
cand=pd.read_csv(track(F/'candidate_pairs.csv'));cand=cand[(cand.moa=='inhibitor/antagonist')&(cand.dataset!='lincs_selected200')].copy()
ann=pd.read_csv(track(F/'drug_annotations.csv')).fillna('')
ann['target_set']=ann.targets.map(lambda s:set(str(s).split(';'))-{''})
lookup=ann.set_index('drug').to_dict('index')
excluded=set().union(*(lookup[d]['target_set'] for d in cand.drug.unique()))
z=np.load(track(T/'panel.npz'),allow_pickle=True)
TD,tgenes=collapse(z['effects'][:,0],z['genes']);TR,_=collapse(z['biological_repeat_effects'],z['genes'])
repmap={int(j):i for i,j in enumerate(z['biological_repeat_indices'])}
lineids={'A549':'CVCL_0023','HT29':'CVCL_0320'};tid={k:list(z['lines']).index(v) for k,v in lineids.items()}
records=[];repeats=[];features=[];checks=[]
def process(name,line,G,Grep,gtargets,ggenes):
 G,ggenes=collapse(G,ggenes)
 if Grep is not None:Grep,_=collapse(Grep,original_genes[name])
 assert len(gtargets)==len(set(gtargets)),(name,line,'duplicated canonical targets')
 common=sorted((set(ggenes)&set(tgenes))-excluded);assert len(common)>500
 gindex={g:i for i,g in enumerate(ggenes)};dindex={g:i for i,g in enumerate(tgenes)}
 gi=[gindex[g] for g in common];di=[dindex[g] for g in common]
 G=G[:,gi];D=TD[tid[line]][:,di];R=TR[tid[line]][:,di]
 if Grep is not None:Grep=Grep[:,:,gi]
 sub=cand[(cand.dataset==name)&(cand.cell_line==line)]
 universe=sorted(sub.target.unique());ix=[gtargets.index(g) for g in universe]
 features.append(dict(dataset=name,cell_line=line,n_genes=len(common),genes=common))
 for mode in ['raw','generic_removed']:
  GG=G if mode=='raw' else G-G.mean(0)
  DD=D if mode=='raw' else D-D.mean(0)
  YY=GG[ix]
  rho=unit(DD,rank=True)@unit(YY,rank=True).T
  pear=unit(DD,center=True)@unit(YY,center=True).T
  cos=unit(DD)@unit(YY).T
  for r in sub.itertuples():
   a=lookup[r.drug];i=int(a['drug_index']);j=universe.index(r.target)
   wrong=np.array([v not in a['target_set'] for v in universe]);assert wrong.any();assert not wrong[j]
   assert not (set(np.array(universe)[wrong])&a['target_set'])
   null=rho[i,wrong];value=float(rho[i,j]);x,y=DD[i],YY[j];mask=(x!=0)&(y!=0)
   records.append(dict(dataset=name,cell_line=line,drug=r.drug,target=r.target,mode=mode,n_genes=len(common),single_listed_target=r.n_annotated_drug_targets==1,spearman=value,pearson=float(pear[i,j]),cosine=float(cos[i,j]),sign_agreement=float(np.mean(np.sign(x[mask])==np.sign(y[mask]))),wrong_spearman=float(null.mean()),spearman_advantage=value-float(null.mean()),matched_percentile=float((np.sum(value>null)+.5*np.sum(value==null))/len(null)),n_wrong=len(null)))
  r=next(sub.itertuples());i=int(lookup[r.drug]['drug_index']);j=universe.index(r.target)
  assert np.isclose(rho[i,j],spearmanr(DD[i],YY[j]).statistic,atol=1e-12)
  assert np.isclose(pear[i,j],pearsonr(DD[i],YY[j]).statistic,atol=1e-12)
  assert np.isclose(cos[i,j],np.dot(DD[i],YY[j])/np.linalg.norm(DD[i])/np.linalg.norm(YY[j]),atol=1e-12)
  checks.append(dict(dataset=name,line=line,mode=mode,scipy_spearman_pearson_and_direct_cosine=True))
 # Raw repeat comparisons; identical feature and pair support for each row.
 repeat_corr={'cross_primary':unit(D,rank=True)@unit(G[ix],rank=True).T, 'cross_drug_repeat':unit(R,rank=True)@unit(G[ix],rank=True).T}
 if Grep is not None:
  for k in range(2):repeat_corr['cross_gene_rep'+str(k+1)]=unit(D,rank=True)@unit(Grep[k,ix],rank=True).T
 for r in sub.itertuples():
  i=int(lookup[r.drug]['drug_index']);j=gtargets.index(r.target)
  base=dict(dataset=name,cell_line=line,drug=r.drug,target=r.target,single_listed_target=r.n_annotated_drug_targets==1,n_genes=len(common),cross_primary=float(spearmanr(D[i],G[j]).statistic))
  if Grep is not None:
   base.update(genetic_repeat_rho=float(spearmanr(Grep[0,j],Grep[1,j]).statistic),cross_gene_rep1=float(spearmanr(D[i],Grep[0,j]).statistic),cross_gene_rep2=float(spearmanr(D[i],Grep[1,j]).statistic))
  if i in repmap:
   dr=R[repmap[i]]
   base.update(drug_repeat_rho=float(spearmanr(D[i],dr).statistic),cross_drug_repeat=float(spearmanr(dr,G[j]).statistic))
   if Grep is not None:
    base.update(cross_both_rep1=float(spearmanr(dr,Grep[0,j]).statistic),cross_both_rep2=float(spearmanr(dr,Grep[1,j]).statistic))
  uj=universe.index(r.target);wrong=np.array([v not in lookup[r.drug]['target_set'] for v in universe])
  for key,matrix in repeat_corr.items():
   if key=='cross_drug_repeat' and i not in repmap:continue
   row=repmap[i] if key=='cross_drug_repeat' else i
   actual=float(matrix[row,uj]);null=float(matrix[row,wrong].mean())
   assert np.isclose(actual,base[key],atol=1e-10)
   base[key+'_wrong']=null;base[key+'_advantage']=actual-null
  repeats.append(base)
 print(name,line,'pairs',len(sub),'genes',len(common),flush=True)
original_genes={}
for path in sorted(P.glob('jiang_*.npz')):
 name=path.stem;q=np.load(track(path),allow_pickle=True);original_genes[name]=q['genes']
 for line in lineids:
  li=list(q['lines']).index(line);gr=q['effects'][li]
  process(name,line,gr.mean(0),gr,[canon(t) for t in q['targets']],q['genes'])
m=pd.read_csv(track(L/'GSE106127_CGS_meta.txt.gz'),sep='\t',dtype=str)
gi=pd.read_csv(track(L/'GSE106127_gene_info.txt.gz'),sep='\t',dtype=str)
gmap=dict(zip(gi.pr_gene_id,gi.pr_gene_symbol))
with h5py.File(track(L/'GSE106127_CGS_n33839x978.gctx'),'r') as h:
 ids=[x.decode() for x in h['0/META/COL/id'][:]];idx={v:i for i,v in enumerate(ids)}
 genes=[gmap[x.decode()] for x in h['0/META/ROW/id'][:]]
 for line in lineids:
  sub=m[(m.cell_id==line)&(m.pert_type=='trt_sh.cgs')&(m.pert_itime=='96 h')].copy()
  sub['row']=sub.sig_id.map(idx);sub=sub.sort_values('row')
  G=h['0/DATA/0/matrix'][sub.row.to_numpy(),:];assert np.isfinite(G).all()
  process('lincs_full_local_CGS',line,G,None,[canon(t) for t in sub.pert_iname],genes)
a=pd.DataFrame(records);rr=pd.DataFrame(repeats)
assert len(a)==2*len(cand)==802
assert len(rr)==len(cand)==401
assert a.drop(columns=['dataset','cell_line','drug','target','mode']).notna().all().all()
a.to_csv(OUT/'pair_metrics.csv',index=False);rr.to_csv(OUT/'repeat_metrics.csv',index=False)
(OUT/'features.json').write_text(json.dumps(features));(OUT/'excluded_target_genes.json').write_text(json.dumps(sorted(excluded)))
metrics=['spearman','pearson','cosine','sign_agreement','wrong_spearman','spearman_advantage','matched_percentile']
sums=[];drugs=[];line_rows=[]
for (name,mode),g in a.groupby(['dataset','mode']):
 for subset in ['all_annotated_inhibitors','single_listed_target']:
  s=g if subset=='all_annotated_inhibitors' else g[g.single_listed_target]
  if s.empty:continue
  byline=s.groupby(['drug','cell_line'])[metrics].mean().reset_index()
  bydrug=byline.groupby('drug')[metrics].mean()
  for drug,v in bydrug.iterrows():drugs.append(dict(dataset=name,mode=mode,subset=subset,drug=drug,**v.to_dict()))
  for line,sl in byline.groupby('cell_line'):line_rows.append(dict(dataset=name,mode=mode,subset=subset,cell_line=line,n_drugs=len(sl),**sl[metrics].mean().to_dict()))
  rng=np.random.default_rng(20260916);v=bydrug[metrics].to_numpy();boot=v[rng.integers(len(v),size=(2000,len(v)))].mean(1)
  row=dict(dataset=name,mode=mode,subset=subset,n_pairs=len(s),n_drugs=len(bydrug),n_targets=s.target.nunique(),n_genes=int(s.n_genes.iloc[0]),**bydrug.mean().to_dict())
  for k,metric in enumerate(metrics):row[metric+'_lo'],row[metric+'_hi']=np.quantile(boot[:,k],[.025,.975])
  sums.append(row)
pd.DataFrame(sums).to_csv(OUT/'summary.csv',index=False);pd.DataFrame(drugs).to_csv(OUT/'by_drug.csv',index=False);pd.DataFrame(line_rows).to_csv(OUT/'by_line.csv',index=False)
# Each repeat comparison uses exactly the rows for which that repeat exists.
rs=[]
for name,g in rr.groupby('dataset'):
 for repeat_kind,needed in [('genetic',['genetic_repeat_rho','cross_gene_rep1','cross_gene_rep2']),('drug',['drug_repeat_rho','cross_primary','cross_drug_repeat'])]:
  if any(k not in g for k in needed):continue
  s=g.dropna(subset=needed)
  if s.empty:continue
  needed=needed+[k+suffix for k in needed if k.startswith('cross_') for suffix in ['_wrong','_advantage']]
  v=s.groupby(['drug','cell_line'])[needed].mean().groupby('drug').mean()
  rs.append(dict(dataset=name,repeat_kind=repeat_kind,n_pairs=len(s),n_drugs=len(v),**v.mean().to_dict()))
pd.DataFrame(rs).to_csv(OUT/'repeat_summary.csv',index=False)
unchanged={p:digest(p)==v for p,v in inputs.items()};assert all(unchanged.values())
(OUT/'verification.json').write_text(json.dumps(dict(real_metric_crosschecks=checks,source_sha256=inputs,all_sources_unchanged=all(unchanged.values()),expected_pairs=401,observed_pairs=len(rr),nulls_exclude_all_listed_drug_targets=True,finite_metrics=True,script_sha256=digest(__file__)),indent=2))
print(pd.DataFrame(sums).query("subset=='all_annotated_inhibitors'")[['dataset','mode','n_drugs','n_genes','spearman','wrong_spearman','spearman_advantage','matched_percentile']].to_string(index=False))
