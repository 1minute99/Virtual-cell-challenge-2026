"""Full-population, excluded-biology preparation for the b19 follow-up."""
import gc
import json
from pathlib import Path
import numpy as np
import context_forcing_data as F
import context_forcing_model as M
import context_forcing_runner as R
import simple_residual_data as D
import simple_residual_baselines as B
import vcclib as V

OUT=D.ROOT/'reports/b19_optimization_2026-09-16'
LINES=F.LINES


def read(path):return json.loads(Path(path).read_text())


def freeze(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if read(path)!=value:raise ValueError(f'Frozen artifact differs: {path}')
    else:D.write_json(path,value)


def scope_ids(table,excluded):
    if not set(excluded)<=set(LINES) or len(set(excluded))!=len(excluded):raise ValueError('Invalid excluded scope')
    return F.permitted_ids(table,excluded)


def decode(prepared,global_latent,context_latent,columns,baseline,beta=1.,gamma=1.):
    if not np.isfinite([beta,gamma]).all() or min(beta,gamma)<0:raise ValueError('Invalid branch strength')
    if beta==0 and gamma==0:return baseline.copy()
    return F.full_delta(prepared,beta*global_latent+gamma*context_latent,columns,
                        baseline.shape[1],empirical=baseline)


class Data:
    def __init__(self):
        self.table=D.Table.load(F.TABLE)
        self.protocol=read(F.OUT/'protocol.json')
        self.targets=self.protocol['final_targets'];self.inner_targets=self.protocol['inner_targets']
        indices={g:i for i,g in enumerate(V.official_features())}
        self.columns=np.array([indices[g] for g in self.table.genes])
        with np.load(F.TABLE/'gene_feature_table.npz') as z:
            self.feature_map={D.canonical_target(str(t)):f.copy() for t,f in zip(z['targetlist'],z['features'])}
        self.queries={};self.available={}
        for line in LINES:
            path=F.OUT/'query_features'/f'{line}.npz'
            with np.load(path) as q:
                self.queries[line]=q['control_log1p_cp10k'][self.columns].copy()
                self.available[line]=q['available_mask'][self.columns].copy()
        self.baseline_cache={}

    def features(self,targets):return np.stack([self.feature_map[D.canonical_target(t)] for t in targets])

    def folder(self,excluded):return OUT/'folds'/('__'.join(sorted(excluded)) if excluded else 'all')

    def prepared(self,excluded,rank=32):
        ids=scope_ids(self.table,excluded);folder=self.folder(excluded);folder.mkdir(parents=True,exist_ok=True)
        scope=dict(excluded=sorted(excluded),ids_sha256=M._array_sha(ids.astype(np.int64)),rows=len(ids),
                   targets=len(set(self.table.targets[ids])),modality='genetic',target_subset=False)
        freeze(folder/'scope.json',scope)
        path=folder/f'prepared_rank{rank}.npz'
        if path.exists():
            receipt=read(path.with_suffix('.json'))
            if D.sha(path)!=receipt['sha256']:raise ValueError('Cached preparation changed')
            prepared=F.Prepared.load(path)
            np.testing.assert_array_equal(prepared.ids,ids)
            if prepared.metadata['rank']!=rank:raise ValueError('Cached rank differs')
            return prepared
        old=F.OUT/'full/folds'/excluded[0] if len(excluded)==1 else None
        if rank==32 and old is not None and (old/'residual.prepared.npz').exists():
            source=old/'residual.prepared.npz';prepared=F.Prepared.load(source)
            np.testing.assert_array_equal(prepared.ids,ids)
            lineage=dict(reused_preparation=str(source),source_sha256=D.sha(source))
        else:
            anchor_path=old/'anchors.npy' if old is not None and (old/'anchors.npy').exists() else folder/'anchors.npy'
            if not anchor_path.exists():
                print(json.dumps(dict(event='anchors',excluded=list(excluded),rows=len(ids))),flush=True)
                anchors=D.training_anchors(self.table,ids,rank=32);R.save_array(anchor_path,anchors);del anchors;gc.collect()
            freeze(folder/'anchors.json',dict(path=str(anchor_path),sha256=D.sha(anchor_path),
                ids_sha256=scope['ids_sha256'],baseline_rank=32,excludes_destination_biology=True))
            print(json.dumps(dict(event='prepare',excluded=list(excluded),rank=rank,rows=len(ids))),flush=True)
            prepared=F.prepare(self.table,ids,route='residual',rank=rank,anchors=np.load(anchor_path,mmap_mode='r'))
            lineage=dict(anchor_receipt=str(folder/'anchors.json'))
        prepared.save(path)
        freeze(path.with_suffix('.json'),dict(sha256=D.sha(path),scope=scope,metadata=prepared.metadata,**lineage))
        return prepared

    def baseline(self,excluded,line,inner=True):
        targets=self.inner_targets[line] if inner else self.targets[line]
        path=self.folder(excluded)/'query'/f'{line}_{"inner" if inner else "full"}.npz'
        identity=dict(excluded=sorted(excluded),line=line,targets=targets,baseline_rank=32)
        if path.exists():
            meta=read(path.with_suffix('.json'))
            if meta['identity']!=identity or meta['sha256']!=D.sha(path):raise ValueError('Baseline identity differs')
            with np.load(path) as z:return {k:z[k] for k in z.files}
        key=tuple(sorted(excluded))
        if key not in self.baseline_cache:
            self.baseline_cache.clear();self.baseline_cache[key]=D.SourceBaseline(self.table,scope_ids(self.table,excluded),rank=32)
        empirical=B.empirical(targets,excluded_bios=excluded);refit=empirical.copy()
        refit[:,self.columns]=self.baseline_cache[key].predict([D.canonical_target(t) for t in targets],self.features(targets))[0]
        result=dict(refit=refit,mean=empirical,b16=B.legacy_b16(empirical,excluded_bios=excluded,held_targets=()),no_change=np.zeros_like(refit))
        path.parent.mkdir(parents=True,exist_ok=True);np.savez(path,**result)
        freeze(path.with_suffix('.json'),dict(identity=identity,sha256=D.sha(path)))
        return result

    def query(self,prepared,line,targets,baseline,context_line=None):
        line=context_line or line
        profile=R.impute_query_controls(self.queries[line],self.available[line],prepared.control_mean)
        h,a=prepared.query(np.broadcast_to(profile,(len(targets),len(self.columns))),baseline[:,self.columns])
        return self.features(targets),h,a

    def truth(self,line,targets):return R.Experiment.truth(self,line,targets)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--scope',required=True);p.add_argument('--rank',type=int,default=32);args=p.parse_args()
    data=Data();excluded=args.scope.split(',') if args.scope!='all' else []
    prep=data.prepared(excluded,args.rank)
    for line in excluded:data.baseline(excluded,line,inner=len(excluded)==2)
    print(json.dumps(dict(event='preparation_complete',scope=excluded,rank=args.rank,rows=len(prep.ids))),flush=True)
