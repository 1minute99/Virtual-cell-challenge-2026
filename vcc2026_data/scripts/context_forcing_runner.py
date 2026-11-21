"""Frozen, resumable whole-line experiment for forced and anchor-free prediction."""
from __future__ import annotations
import argparse
from dataclasses import asdict, replace
import gc
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
import context_forcing_data as F
import context_forcing_model as M
import simple_residual_data as D
import simple_residual_baselines as B
import vcclib as V

OUT = F.OUT
VCCPY = '/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python'
LINES = F.LINES


def candidates():
    base = M.Config(steps=800, coverage_pass=True)
    return {
        'res_standard': replace(base, energy=.5, weight_decay=.03),
        'res_relaxed': base,
        'res_no_anchor_input': replace(base, use_baseline=False),
        'res_anchor_dropout': replace(base, baseline_dropout=.5),
        'res_context_only': replace(base, family='interaction', use_baseline=False),
        'res_contrast1': replace(base, family='interaction', use_baseline=False, pair_weight=1.),
        'res_contrast4': replace(base, family='interaction', use_baseline=False, pair_weight=4.),
        'direct_gene': replace(base, route='direct', use_baseline=False, mode='none'),
        'direct_context': replace(base, route='direct', use_baseline=False),
        'direct_contrast1': replace(base, route='direct', use_baseline=False, pair_weight=1.),
        'direct_contrast4': replace(base, route='direct', use_baseline=False, pair_weight=4.),
    }


def common_only(delta, columns):
    out = np.zeros_like(delta)
    out[:, np.asarray(columns)] = delta[:, np.asarray(columns)]
    return out


def impute_query_controls(profile, available, source_mean):
    """Missing assay coordinates receive permitted-source means, never zeros/truth."""
    profile=np.asarray(profile,dtype=np.float32).copy()
    available=np.asarray(available,dtype=bool)
    source_mean=np.asarray(source_mean,dtype=np.float32)
    if profile.shape!=available.shape or profile.shape!=source_mean.shape:
        raise ValueError('Query coverage shape mismatch')
    profile[~available]=source_mean[~available]
    if not np.isfinite(profile).all(): raise ValueError('Nonfinite query control profile')
    return profile


def select_complete(frame, outer, names):
    expected = set(LINES)-{outer}
    scores = []
    for name in names:
        part = frame[(frame.outer == outer) & (frame.name == name)]
        if len(part) != len(expected) or set(part.inner) != expected:
            raise ValueError(f'Incomplete inner evaluation: {outer}/{name}')
        if not part.valid.all() or not np.isfinite(part.pds).all():
            raise ValueError(f'Invalid inner evaluation: {outer}/{name}')
        scores.append(dict(name=name, pds=float(part.pds.mean()), minimum=float(part.pds.min())))
    if not scores:
        raise ValueError('No eligible configurations')
    return sorted(scores, key=lambda x:(-x['pds'], -x['minimum'], x['name']))[0]


def native_pds(pred, truth):
    p,t = np.asarray(pred,dtype=float),np.asarray(truth,dtype=float)
    den = np.linalg.norm(p,axis=1)[:,None]*np.linalg.norm(t,axis=1)[None,:]
    sim = np.divide(p@t.T,den,out=np.zeros((len(p),len(t))),where=den>0)
    own = np.diag(sim)[:,None]
    return float(np.mean(((sim<own).sum(1)+.5*((sim==own).sum(1)-1))/(len(p)-1)))


def array_sha(x):
    x=np.ascontiguousarray(x)
    return hashlib.sha256(str((x.shape,x.dtype.str)).encode()+x.tobytes()).hexdigest()


def save_array(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if not np.array_equal(np.load(path), value):
            raise ValueError(f'Refusing to replace different array: {path}')
        return
    with path.with_suffix('.tmp').open('wb') as f:
        np.save(f,value)
    path.with_suffix('.tmp').replace(path)


def write_once(path, value):
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError(f'Frozen record differs: {path}')
    else:
        D.write_json(path,value)


class Scorer:
    def __init__(self):
        self.log=(OUT/'scorer.log').open('a')
        self.process=subprocess.Popen([VCCPY,str(D.ROOT/'scripts/context_forcing_score.py')],
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,bufsize=1)
    def request(self, request):
        self.process.stdin.write(json.dumps(request)+'\n'); self.process.stdin.flush()
        while True:
            line=self.process.stdout.readline()
            if not line:
                raise RuntimeError('Scorer exited; inspect scorer.log')
            try: result=json.loads(line)
            except json.JSONDecodeError:
                self.log.write(line); self.log.flush(); continue
            if 'ok' not in result:
                self.log.write(line); continue
            if not result['ok']:
                raise RuntimeError(result)
            return result['result']
    def close(self):
        if self.process.poll() is None:
            self.request({'command':'quit'}); self.process.wait()
        self.log.close()


class Experiment:
    def __init__(self):
        OUT.mkdir(parents=True,exist_ok=True)
        self.table=D.Table.load(F.TABLE)
        official=V.official_features(); index={g:i for i,g in enumerate(official)}
        self.columns=np.array([index[g] for g in self.table.genes])
        z=np.load(F.TABLE/'gene_feature_table.npz')
        self.feature_map={D.canonical_target(str(t)):f for t,f in zip(z['targetlist'],z['features'])}
        self.pilot=F.pilot_targets(self.table)
        self.configs=candidates()
        self.scorer=Scorer()
        self.queries={}
        self.available={}
        for line in LINES:
            path=OUT/'query_features'/f'{line}.npz'
            self.scorer.request(dict(command='export_query',line=line,output_path=str(path)))
            with np.load(path) as q:
                self.available[line]=q['available_mask'][self.columns]
                self.queries[line]=q['control_log1p_cp10k'][self.columns].astype(np.float32)
        mock=json.loads((D.OUT/'mock/selection.json').read_text())
        kolf=json.loads((D.OUT/'kolf_full_axis_2026-09-16/axis.json').read_text())
        self.targets={line:sorted(kolf['targets'] if line=='KOLF' else mock['targets']) for line in LINES}
        self.small_targets={line:sorted(np.random.default_rng(7301).choice(ts,min(64,len(ts)),replace=False))
                            for line,ts in self.targets.items()}
        self.freeze()

    def freeze(self):
        protocol=dict(version=1, configs={k:asdict(v) for k,v in self.configs.items()},
            lines=list(LINES),pilot_targets=self.pilot,inner_targets=self.small_targets,
            final_targets=self.targets,table_manifest_sha256=D.sha(F.TABLE/'manifest.json'),
            common_axis_sha256=D.sha(F.TABLE/'genes.json'),
            missing_query_coordinates={line:np.asarray(self.table.genes)[~mask].tolist() for line,mask in self.available.items()},
            missing_query_policy='Impute unavailable control coordinates with permitted-fold source control mean. They remain absent from the native scoring axis. KOLF support is 5583 of 5592 common genes.',
            selection='Within each outer exclusion, native PDS shortlist two per route across four whole-line inner exclusions; installed count PDS across all four inner lines selects one per route. Outer responses never select.',
            pilot_steps=800,full_steps=2500,coverage_pass=True,seeds=[0,1,2],alpha=1.,
            inner_cells=64,full_cells=400,reconstruction_scale=3.,null_repeats=32,
            scopes='Genetic rows only; entire outer/inner biology excluded across studies and stimuli. Outer query controls only. Seen-target transfer, not gene-holdout generalization.',
            support='Full-axis residuals preserve source baseline outside 5592 coordinates; direct outputs zero outside; support-matched baseline/residual scores also reported.',
            controls='Full separately refitted mode none and context_shuffle, seed0, same selected config. If gene-only wins these are explicit identical aliases. Same-checkpoint zero and next-line swapped query additionally scored.',
            oracle='Evaluation-only projected observed effect/residual; diagnostic raw MSE energy, never an input or selection candidate.',
            prospective='Protocol frozen before new experiment scores; biological panels have been studied previously and are not an untouched external validation set.',
            scripts={p.name:D.sha(p) for p in [Path(__file__),D.ROOT/'scripts/context_forcing_data.py',
                  D.ROOT/'scripts/context_forcing_model.py',D.ROOT/'scripts/context_forcing_score.py']})
        write_once(OUT/'protocol.json',protocol)

    def features(self, targets):
        return np.stack([self.feature_map[D.canonical_target(t)] for t in targets])

    def truth(self,line,targets):
        # Access only after a prediction is frozen; this function is never used by fit/prepare.
        sources=self.table.sources.astype(str)
        source=(np.char.startswith(sources,'xl/kolf:') if line=='KOLF'
                else sources == 'rn/'+line)
        result=[]
        for target in targets:
            ids=np.flatnonzero(source & (self.table.targets==D.canonical_target(target)))
            if not len(ids): raise ValueError(f'No evaluation truth: {line}/{target}')
            result.append(np.average(self.table.y[ids],axis=0,weights=self.table.weights[ids]))
        return np.stack(result)

    def folder(self,excluded,stage):
        return OUT/stage/'folds'/'__'.join(sorted(excluded))

    def preparation(self, excluded, stage):
        folder=self.folder(excluded,stage); folder.mkdir(parents=True,exist_ok=True)
        ids=F.permitted_ids(self.table,excluded,self.pilot if stage=='pilot' else None)
        signature=dict(ids_sha256=array_sha(ids),excluded=sorted(excluded),rows=len(ids),
                       targets=len(set(self.table.targets[ids])),modality='genetic',stage=stage)
        write_once(folder/'scope.json',signature)
        ready={}
        anchors=None
        for route in ('direct','residual'):
            path=folder/f'{route}.prepared.npz'
            if path.exists():
                ready[route]=F.Prepared.load(path)
                if not np.array_equal(ready[route].ids,ids): raise ValueError('Cached preparation ids differ')
                continue
            if route=='residual':
                cache=folder/'anchors.npy'
                if cache.exists(): anchors=np.load(cache,mmap_mode='r')
                else:
                    old=(D.OUT/'finish02/full_confirmation/folds'/D.canonical_bio(excluded[0])/'training_anchors.npy'
                         if stage=='full' and len(excluded)==1 and excluded[0]!='KOLF' else Path('/nonexistent'))
                    if old.exists():
                        old_ids=D.scope(self.table,excluded_bios=excluded)
                        keep=self.table.modality[old_ids]==0
                        if not np.array_equal(old_ids[keep],ids): raise ValueError('Old anchor scope mismatch')
                        old_array=np.load(old,mmap_mode='r')
                        if old_array.shape != (len(old_ids),len(self.columns)): raise ValueError('Old anchor shape mismatch')
                        anchors=np.asarray(old_array[keep])
                        D.write_json(folder/'anchor_reuse.json',dict(path=str(old),sha256=D.sha(old),
                            assertion='Genetic anchors fitted independently of drug modality, exact full excluded-biology ids'))
                    else: anchors=D.training_anchors(self.table,ids,rank=32)
                    save_array(cache,anchors)
            print(json.dumps(dict(event='prepare',stage=stage,excluded=excluded,route=route,rows=len(ids))),flush=True)
            ready[route]=F.prepare(self.table,ids,route=route,rank=32,anchors=anchors if route=='residual' else None)
            ready[route].save(path)
            D.write_json(folder/f'{route}.preparation.json',ready[route].metadata)
        del anchors; gc.collect()
        return ready,ids

    def baselines(self, excluded, ids, targets, folder):
        path=folder/'baselines.npz'
        if path.exists():
            with np.load(path) as z: return {k:z[k] for k in z.files}
        features=self.features(targets); canonical=[D.canonical_target(t) for t in targets]
        mean=B.empirical(targets,excluded_bios=excluded)
        b16=B.legacy_b16(mean,excluded_bios=excluded,held_targets=())
        refit=mean.copy(); fit=D.SourceBaseline(self.table,ids,rank=32)
        refit[:,self.columns]=fit.predict(canonical,features)[0]
        result=dict(mean=mean,b16=b16,refit=refit,no_change=np.zeros_like(mean))
        folder.mkdir(parents=True,exist_ok=True); np.savez(path,**result)
        return result

    def fitted(self,folder,prepared,config,seed):
        if prepared.metadata['route']!=config.route or prepared.y.shape[1]!=config.rank:
            raise ValueError('Prepared route/rank differs from requested model')
        path=folder/'model.pt'; folder.mkdir(parents=True,exist_ok=True)
        if path.exists():
            model=M.Fitted.load(path,device='cuda')
            if model.config!=config or model.seed!=seed or model.receipt['training_ids_sha256']!=array_sha(prepared.ids):
                raise ValueError(f'Cached model scope/config differs: {path}')
            return model
        model=M.fit(prepared,config,seed=seed)
        model.save(path)
        D.write_json(folder/'receipt.json',dict(config=asdict(config),**model.receipt))
        print(json.dumps(dict(event='fit',path=str(folder),seconds=model.receipt['seconds'],
                             rows=model.receipt['training_rows'],unique=model.receipt['unique_rows_seen'])),flush=True)
        return model

    def predict(self,model,prep,line,targets,refit,zero=False,context_line=None):
        query=context_line or line
        profile=impute_query_controls(self.queries[query],self.available[query],prep.control_mean)
        c=np.broadcast_to(profile,(len(targets),len(self.columns)))
        h,a=prep.query(c,refit[:,self.columns] if prep.metadata['route']=='residual' else None)
        latent=model.predict(self.features(targets),h,a,zero_context=zero)
        return F.full_delta(prep,latent,self.columns,
                            empirical=refit if prep.metadata['route']=='residual' else None)

    def score(self,line,targets,delta_path,output,seed=0,inner=False):
        result=self.scorer.request(dict(line=line,targets=targets,delta_path=str(delta_path),
            output_path=str(output),seed=seed,cells=64 if inner else 400,scale=3.,null_repeats=0 if inner else 32))
        valid=result['count_diagnostics']['official_per_cell_limit_valid']
        if not valid: raise ValueError(f'Invalid reconstructed output: {delta_path}')
        return dict(pds=result['pds'],valid=valid,null_mean=result['random_null']['mean'],
                    score_path=str(output))

    def pilot_screen(self):
        rows=[]
        for excluded in itertools.combinations(LINES,2):
            ready,ids=self.preparation(excluded,'pilot'); folder=self.folder(excluded,'pilot')
            base={line:self.baselines(excluded,ids,self.small_targets[line],folder/'query'/line) for line in excluded}
            for name,config in self.configs.items():
                model=self.fitted(folder/name,ready[config.route],config,0)
                for line in excluded:
                    path=folder/name/f'{line}.delta.npy'
                    if not path.exists():
                        save_array(path,self.predict(model,ready[config.route],line,self.small_targets[line],base[line]['refit']))
                    prediction=np.load(path)[:,self.columns]
                    truth=self.truth(line,self.small_targets[line])
                    rows.append(dict(outer=next(x for x in excluded if x!=line),inner=line,name=name,
                        route=config.route,pds=native_pds(prediction,truth),valid=True,delta_path=str(path)))
                del model
            pd.DataFrame(rows).to_csv(OUT/'pilot/native_scores.csv',index=False)
            del ready; gc.collect()
        frame=pd.DataFrame(rows)
        shortlist={}
        for outer in LINES:
            shortlist[outer]={}
            for route in ('residual','direct'):
                names=[k for k,v in self.configs.items() if v.route==route]
                ranked=[]
                for name in names:
                    ranked.append(select_complete(frame,outer,[name]))
                ranked.sort(key=lambda x:(-x['pds'],-x['minimum'],x['name']))
                shortlist[outer][route]=[x['name'] for x in ranked[:2]]
        write_once(OUT/'pilot/shortlist.json',shortlist)
        count_rows=[]
        for outer, routes in shortlist.items():
            for route,names in routes.items():
                for name in names:
                    for inner in LINES:
                        if inner==outer: continue
                        row=frame[(frame.outer==outer)&(frame.inner==inner)&(frame.name==name)].iloc[0]
                        path=Path(row.delta_path)
                        scored=self.score(inner,self.small_targets[inner],path,path.with_suffix('.score.json'),inner=True)
                        count_rows.append(dict(outer=outer,inner=inner,name=name,route=route,**scored))
        count=pd.DataFrame(count_rows); count.to_csv(OUT/'pilot/count_scores.csv',index=False)
        selection={outer:{route:select_complete(count,outer,names)
                          for route,names in routes.items()} for outer,routes in shortlist.items()}
        write_once(OUT/'selection.json',selection)
        print(json.dumps(dict(event='selection',selection=selection)),flush=True)

    def prepare_full(self):
        for line in LINES:
            ready,ids=self.preparation((line,),'full')
            self.baselines((line,),ids,self.targets[line],self.folder((line,),'full')/'query'/line)
            del ready; gc.collect()

    def full(self, lines=LINES):
        selection=json.loads((OUT/'selection.json').read_text()); rows=[]
        diagnostics=[]
        for line in lines:
            excluded=(line,); ready,ids=self.preparation(excluded,'full')
            folder=self.folder(excluded,'full'); targets=self.targets[line]
            base=self.baselines(excluded,ids,targets,folder/'query'/line)
            predictions={}
            for name,delta in base.items():
                predictions[name]=delta
                if name!='no_change': predictions[name+'_common']=common_only(delta,self.columns)
            for route in ('residual','direct'):
                selected=selection[line][route]['name']; baseconfig=replace(self.configs[selected],steps=2500)
                prep=ready[route]
                for mode,seed in [('correct',s) for s in (0,1,2)]+[('none',0),('context_shuffle',0)]:
                    # The selected gene-only control remains context-free in the correct run.
                    actual_mode=baseconfig.mode if mode=='correct' or baseconfig.mode=='none' else mode
                    config=replace(baseconfig,mode=actual_mode)
                    alias=(baseconfig.mode=='none' and mode!='correct')
                    fit_mode='correct' if alias else mode
                    model=self.fitted(folder/f'{route}_{fit_mode}_seed{seed}',prep,config,seed)
                    delta=self.predict(model,prep,line,targets,base['refit'])
                    name=f'{route}_{mode}_seed{seed}'
                    predictions[name]=delta
                    if route=='residual': predictions[name+'_common']=common_only(delta,self.columns)
                    if mode=='correct' and seed==0:
                        zero=self.predict(model,prep,line,targets,base['refit'],zero=True)
                        predictions[f'{route}_zero_query']=zero
                        if route=='residual': predictions[f'{route}_zero_query_common']=common_only(zero,self.columns)
                        swapped_line=LINES[(LINES.index(line)+1)%len(LINES)]
                        swapped=self.predict(model,prep,line,targets,base['refit'],context_line=swapped_line)
                        predictions[f'{route}_swapped_query']=swapped
                        if route=='residual': predictions[f'{route}_swapped_query_common']=common_only(swapped,self.columns)
                    # Held-out MSE and projection ceilings only after fitting and fixed selection.
                    truth=self.truth(line,targets)
                    effect=(delta-base['refit'])[:,self.columns] if route=='residual' else delta[:,self.columns]
                    truth_objective=truth-base['refit'][:,self.columns] if route=='residual' else truth
                    truth_objective-=truth_objective.mean(0)
                    projection=(truth_objective@prep.output_basis)@prep.output_basis.T
                    zero_query=self.predict(model,prep,line,targets,base['refit'],zero=True)
                    diagnostics.append(dict(line=line,route=route,mode=mode,seed=seed,selected=selected,
                        actual_context_mode=actual_mode,identical_gene_only_alias=alias,
                        output_norm=float(np.linalg.norm(effect)),baseline_norm=float(np.linalg.norm(base['refit'][:,self.columns])),
                        zero_query_difference=float(np.linalg.norm(delta[:,self.columns]-zero_query[:,self.columns])),
                        mse=float(np.mean((effect-truth_objective)**2)),zero_mse=float(np.mean(truth_objective**2)),
                        oracle_projection_mse=float(np.mean((projection-truth_objective)**2)),
                        oracle_energy_fraction=float(np.sum(projection**2)/np.sum(truth_objective**2)),
                        initial_train_mse=model.receipt['trace'][0]['train_mse'],
                        final_train_mse=model.receipt['trace'][-1]['train_mse'],
                        training_rows=len(ids),training_targets=len(set(self.table.targets[ids])),
                        unique_rows_seen=model.receipt['unique_rows_seen']))
                    del model
            for name,delta in predictions.items():
                path=folder/'predictions'/f'{name}.npy'; save_array(path,delta)
                seeds=[int(name.split('seed')[1].split('_')[0])] if '_seed' in name else ([0] if '_query' in name else [0,1,2])
                for seed in seeds:
                    output=folder/'scores'/f'{name}.reconstruction_seed{seed}.json'
                    scored=self.score(line,targets,path,output,seed=seed)
                    rows.append(dict(line=line,name=name,seed=seed,**scored))
                    print(json.dumps(dict(event='score',line=line,name=name,seed=seed,pds=scored['pds'])),flush=True)
            pd.DataFrame([r for r in rows if r['line']==line]).to_csv(folder/'count_scores.csv',index=False)
            pd.DataFrame([r for r in diagnostics if r['line']==line]).to_csv(folder/'diagnostics.csv',index=False)
            del ready,predictions,base; gc.collect()
        if set(lines)==set(LINES):
            self.collect()

    def collect(self):
        score_frames=[];diagnostic_frames=[]
        for line in LINES:
            folder=self.folder((line,),'full')
            scores=pd.read_csv(folder/'count_scores.csv')
            diag=pd.read_csv(folder/'diagnostics.csv')
            if len(scores)!=42 or len(diag)!=10 or set(scores.line)!={line} or set(diag.line)!={line}:
                raise ValueError(f'Incomplete full fold: {line}')
            score_frames.append(scores);diagnostic_frames.append(diag)
        scores=pd.concat(score_frames,ignore_index=True);diag=pd.concat(diagnostic_frames,ignore_index=True)
        scores.to_csv(OUT/'full/count_scores.csv',index=False)
        diag.to_csv(OUT/'full/diagnostics.csv',index=False)
        D.write_json(OUT/'research_status.json',dict(status='experiments_complete_review_pending',
            full_score_rows=len(scores),full_model_records=len(diag),official_submission='none in this follow-up',
            selection=str(OUT/'selection.json'),count_scores=str(OUT/'full/count_scores.csv')))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--stage',choices=('pilot','prepare-full','full','collect','all'),default='all')
    parser.add_argument('--lines',nargs='+',choices=LINES,default=list(LINES))
    args=parser.parse_args()
    exp=Experiment()
    try:
        if args.stage in ('pilot','all'): exp.pilot_screen()
        if args.stage=='prepare-full': exp.prepare_full()
        if args.stage in ('full','all'): exp.full(tuple(args.lines))
        if args.stage=='collect': exp.collect()
    finally: exp.scorer.close()


if __name__=='__main__': main()
