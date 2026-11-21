"""Nested convergence, branch and capacity experiments; historical assets are read-only."""
from dataclasses import asdict,replace
import argparse
import gc
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import numpy as np
import torch
import b19_optimization_data as A
import b19_optimization_model as O
import context_forcing_model as M
import context_forcing_runner as R
import context_forcing_data as F
import simple_residual_data as D

OUT=A.OUT;LINES=A.LINES;STEPS=(2500,10000,25000)
STRENGTHS=tuple(itertools.product((0.,.5,1.),(0.,.25,.5,1.)))
VCCPY='/home/wonmin-kim/miniconda3/envs/vcc2026/bin/python'


def select(rows,outer):
    part=[r for r in rows if r['outer']==outer];expected=set(LINES)-{outer};ranked=[]
    for name in sorted(set(r['name'] for r in part)):
        rr=[r for r in part if r['name']==name]
        if len(rr)!=4 or {r['inner'] for r in rr}!=expected:raise ValueError(f'Incomplete or contaminated inner set: {outer}/{name}')
        if not all(r['valid'] and r['pds'] is not None and math.isfinite(r['pds']) for r in rr):continue
        spec=rr[0].get('spec',{})
        ranked.append(dict(name=name,pds=statistics.mean(r['pds'] for r in rr),minimum=min(r['pds'] for r in rr),
            spec=spec,beta=rr[0].get('beta',1.),gamma=rr[0].get('gamma',1.),
            by_line={r['inner']:r['pds'] for r in rr}))
    if not ranked:raise ValueError('No valid complete source-inner candidate')
    return sorted(ranked,key=lambda x:(-x['pds'],-x['minimum'],x['spec'].get('steps',0),
        x['spec'].get('rank',0),x['spec'].get('hidden',0),x['name']))[0]


def promising(candidate,reference):
    return all(x is not None and math.isfinite(x) for x in candidate+reference) and statistics.mean(candidate)-statistics.mean(reference)>=.001 and sum(a>b for a,b in zip(candidate,reference))>=3


def key(spec):return hashlib.sha256(json.dumps(asdict(spec),sort_keys=True).encode()).hexdigest()[:16]
def name(spec,beta=1.,gamma=1.):return f'{spec.method}_{spec.schedule}_s{spec.steps}_r{spec.rank}_h{spec.hidden}_c{spec.context_penalty:g}_b{beta:g}_g{gamma:g}'


def freeze_protocol():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=[Path(__file__),Path(A.__file__),Path(O.__file__),Path(M.__file__),Path(F.__file__),Path(D.__file__),Path(R.__file__),
           D.ROOT/'scripts/simple_residual_baselines.py',
           D.ROOT/'scripts/context_forcing_score.py',D.ROOT/'scripts/context_forcing_fast_worker.py',D.ROOT/'scripts/context_forcing_fast_bulk.py']
    value=dict(lines=list(LINES),steps=list(STEPS),schedules=['constant','cosine'],horizon=25000,
        base=asdict(O.TrainSpec()),branch_methods=['joint_total','joint_split','staged'],context_penalties=[.5,2.],
        strengths=[list(x) for x in STRENGTHS],inner_targets=A.read(F.OUT/'protocol.json')['inner_targets'],
        inner_population='all permitted genetic rows, no pilot target subset',inner_cells=400,outer_cells=400,
        reconstruction_scale=3.,outer_seeds=[0,1,2],null_repeats=32,baseline_rank=32,
        stage_order=['convergence','branches','rank64','conditional_rank128_width128','outer'],
        native_shortlist=2,capacity_gate='mean gain >= .001 and >=3/4 inner line wins',
        deployment_gate='mean outer PDS gain >= .001, >=3/5 line wins, >=9/15 paired wins, all valid',
        scientific_scope='Seen-target genetic transfer, source-only nested selection; historically studied panels.',
        official_uploads=0,baseline_preserved='b19',
        source_hashes={str(p):D.sha(p) for p in paths},
        table_manifest_sha256=D.sha(F.TABLE/'manifest.json'),old_protocol_sha256=D.sha(F.OUT/'protocol.json'),
        input_hashes={str(p):D.sha(p) for p in [F.TABLE/'genes.json',F.TABLE/'gene_feature_table.npz',
            *[F.OUT/'query_features'/f'{line}.npz' for line in LINES]]})
    A.freeze(OUT/'protocol.json',value)


class Worker:
    def __init__(self,tag):
        self.log=(OUT/f'{tag}.worker.log').open('a')
        env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
        self.process=subprocess.Popen([VCCPY,str(D.ROOT/'scripts/context_forcing_fast_worker.py'),'--report-root',str(OUT)],
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,bufsize=1,env=env)
    def request(self,request):
        self.process.stdin.write(json.dumps(request)+'\n');self.process.stdin.flush()
        while True:
            line=self.process.stdout.readline()
            if not line:raise RuntimeError('Score worker exited')
            try:r=json.loads(line)
            except json.JSONDecodeError:self.log.write(line);continue
            if 'ok' not in r:self.log.write(line);continue
            if not r['ok']:raise RuntimeError(r)
            return r['result']
    def close(self):
        if self.process.poll() is None:self.request({'command':'quit'});self.process.wait()
        self.log.close()


class Experiment:
    def __init__(self,tag):
        freeze_protocol();self.data=A.Data();self.worker=Worker(tag)
    def close(self):self.worker.close()

    def trained(self,scope,spec,seed=0,points=None):
        points=list(points or [spec.steps]);folder=self.data.folder(scope)/'training'/f'{key(spec)}_seed{seed}'
        prep=self.data.prepared(scope,spec.rank)
        identity=dict(spec=asdict(spec),seed=seed,
            ids_sha256=M._array_sha(prep.ids.astype(np.int64)),prepared_sha256=D.sha(self.data.folder(scope)/f'prepared_rank{spec.rank}.npz'))
        if (folder/'request.json').exists():
            prior=A.read(folder/'request.json')
            if any(prior[k]!=v for k,v in identity.items()):raise ValueError('Training request differs')
            if not set(points)<=set(prior['points']):raise ValueError('Requested checkpoints not in original trajectory')
            if (folder/'COMPLETE.json').exists():
                complete=A.read(folder/'COMPLETE.json')
                if any(D.sha(folder/f'step{step}'/'model.pt')!=digest for step,digest in complete['hashes'].items()):raise ValueError('Completed checkpoint hash differs')
                models={step:M.Fitted.load(folder/f'step{step}'/'model.pt') for step in points}
                for step,model in models.items():model.checkpoint_path=folder/f'step{step}'/'model.pt'
                return prep,models
        elif len(points)==1 and spec.method=='joint_total' and spec.steps<max(STEPS):
            longer=replace(spec,steps=max(STEPS));prefix=self.data.folder(scope)/'training'/f'{key(longer)}_seed{seed}'/f'step{spec.steps}'/'model.pt'
            if prefix.exists():
                model=M.Fitted.load(prefix)
                if model.receipt['train_spec']!=asdict(longer) or model.receipt['checkpoint_step']!=spec.steps or model.receipt['training_ids_sha256']!=identity['ids_sha256']:
                    raise ValueError('Trajectory prefix identity differs')
                complete=A.read(prefix.parent.parent/'COMPLETE.json')
                if D.sha(prefix)!=complete['hashes'][str(spec.steps)]:raise ValueError('Prefix checkpoint hash differs')
                model.checkpoint_path=prefix
                return prep,{spec.steps:model}
        A.freeze(folder/'request.json',dict(**identity,points=points))
        paths={p:folder/f'step{p}'/'model.pt' for p in points}
        if all(p.exists() for p in paths.values()):
            models={step:M.Fitted.load(path) for step,path in paths.items()}
        else:
            print(json.dumps(dict(event='fit',scope=scope,spec=asdict(spec),seed=seed,rows=len(prep.ids))),flush=True)
            models=O.fit(prep,spec,seed=seed,checkpoints=points)
            for step,model in models.items():
                path=paths[step];path.parent.mkdir(parents=True,exist_ok=True)
                if path.exists():
                    old=M.Fitted.load(path)
                    if any(not torch.equal(t,old.module.state_dict()[k]) for k,t in model.module.state_dict().items()):raise ValueError('Checkpoint replay differs')
                    models[step]=old
                else:model.save(path);D.write_json(path.with_name('receipt.json'),model.receipt)
        for step,model in models.items():
            if model.receipt['train_spec']!=asdict(spec) or model.receipt['checkpoint_step']!=step or model.receipt['training_ids_sha256']!=M._array_sha(prep.ids.astype(np.int64)):
                raise ValueError('Cached training identity mismatch')
            model.checkpoint_path=paths[step]
        A.freeze(folder/'COMPLETE.json',dict(hashes={str(step):D.sha(path) for step,path in paths.items()}))
        return prep,models

    def prediction(self,scope,spec,model,prep,line,beta=1.,gamma=1.,seed=0,inner=True,zero=False,context_line=None):
        targets=self.data.inner_targets[line] if inner else self.data.targets[line]
        baseline=self.data.baseline(scope,line,inner)['refit']
        query=self.data.query(prep,line,targets,baseline,context_line)
        g,c=O.predict_components(model,*query,zero=zero)
        pred=A.decode(prep,g,c,self.data.columns,baseline,beta,gamma)
        tag=('zero' if zero else ('query_'+context_line if context_line else 'correct'))
        folder=OUT/'predictions'/('inner' if inner else 'outer')/'__'.join(sorted(scope))/line/f'{key(spec)}_s{seed}_{tag}'
        path=folder/f'b{beta:g}_g{gamma:g}.npy';R.save_array(path,pred)
        A.freeze(path.with_suffix('.manifest.json'),dict(checkpoint=str(model.checkpoint_path),checkpoint_sha256=D.sha(model.checkpoint_path),
            beta=beta,gamma=gamma,scope=sorted(scope),line=line,targets=targets,seed=seed,tag=tag,sha256=D.sha(path)))
        return path

    def score(self,line,path,seed=0,inner=True):
        targets=self.data.inner_targets[line] if inner else self.data.targets[line]
        output=path.with_name(path.stem+f'_score_seed{seed}.json')
        request=dict(line=line,targets=targets,delta_path=str(path),output_path=str(output),seed=seed,cells=400,scale=3.,null_repeats=0 if inner else 32)
        try:
            result=self.worker.request(request)
            valid=bool(result['count_diagnostics']['official_per_cell_limit_valid'] and result['count_diagnostics']['finite_nonnegative_integral'])
            return dict(pds=result['pds'],valid=valid,null_mean=result['random_null']['mean'],score_path=str(output),
                        score_sha256=D.sha(output),prediction_path=str(path),prediction_sha256=D.sha(path))
        except RuntimeError as exc:
            payload=exc.args[0]
            allowed={'type':'ValueError','message':'reconstruction has integer counts outside exact float32 range'}
            if not isinstance(payload,dict) or payload.get('error')!=allowed:raise
            failure=output.with_name(output.stem+'_unscorable.json')
            A.freeze(failure,dict(status='unscorable',request=request,error=allowed,prediction_sha256=D.sha(path)))
            return dict(pds=None,valid=False,null_mean=None,score_path=str(failure),score_sha256=D.sha(failure),prediction_path=str(path),prediction_sha256=D.sha(path))

    def convergence(self,scope):
        rows=[]
        for schedule in ('constant','cosine'):
            spec=O.TrainSpec(schedule=schedule)
            prep,models=self.trained(scope,spec,points=STEPS)
            for step,model in models.items():
                candidate=replace(spec,steps=step)
                for line in scope:
                    outer=next(x for x in scope if x!=line)
                    path=self.prediction(scope,candidate,model,prep,line)
                    row=dict(outer=outer,inner=line,name=name(candidate),spec=asdict(candidate),beta=1.,gamma=1.,**self.score(line,path))
                    rows.append(row);D.write_json(self.data.folder(scope)/'convergence.json',rows)
                    print(json.dumps(dict(event='convergence_score',outer=outer,inner=line,name=row['name'],pds=row['pds'])),flush=True)

    def collect_convergence(self):
        rows=[]
        for scope in itertools.combinations(LINES,2):rows.extend(A.read(self.data.folder(scope)/'convergence.json'))
        if len(rows)!=120:raise ValueError('Incomplete convergence grid')
        picks={line:select(rows,line) for line in LINES}
        A.freeze(OUT/'convergence/selection.json',picks);D.write_json(OUT/'convergence/scores.json',rows)
        print(json.dumps(dict(event='convergence_selected',selections=picks)),flush=True)

    def evaluate_candidate(self,outer,spec,beta=1.,gamma=1.,native=False):
        rows=[]
        for inner in LINES:
            if inner==outer:continue
            scope=sorted((outer,inner));prep,models=self.trained(scope,spec);model=models[spec.steps]
            path=self.prediction(scope,spec,model,prep,inner,beta,gamma)
            if native:
                truth=self.data.truth(inner,self.data.inner_targets[inner]);pds=R.native_pds(np.load(path)[:,self.data.columns],truth)
                score=dict(pds=pds,valid=True,prediction_path=str(path))
            else:score=self.score(inner,path)
            rows.append(dict(outer=outer,inner=inner,name=name(spec,beta,gamma),spec=asdict(spec),beta=beta,gamma=gamma,**score))
        return rows

    def branches(self,outer):
        pick=A.read(OUT/'convergence/selection.json')[outer];base=O.TrainSpec(**pick['spec'])
        specs=[base]+[replace(base,method=m,context_penalty=p) for m in ('joint_split','staged') for p in (.5,2.)]
        native=[];count=[];shortlist={}
        # Train one model per scope/method, decode all branch weights without refitting.
        for spec in specs:
            part=[]
            for inner in LINES:
                if inner==outer:continue
                scope=sorted((outer,inner));prep,models=self.trained(scope,spec);model=models[spec.steps]
                truth=self.data.truth(inner,self.data.inner_targets[inner])
                for beta,gamma in STRENGTHS:
                    path=self.prediction(scope,spec,model,prep,inner,beta,gamma)
                    part.append(dict(outer=outer,inner=inner,name=name(spec,beta,gamma),spec=asdict(spec),beta=beta,gamma=gamma,
                        pds=R.native_pds(np.load(path)[:,self.data.columns],truth),valid=True,prediction_path=str(path)))
            native.extend(part)
            ranks=[select([r for r in part if r['name']==n],outer) for n in sorted(set(r['name'] for r in part))]
            ranks.sort(key=lambda r:(-r['pds'],-r['minimum'],r['name']))
            chosen={r['name'] for r in ranks[:2]}
            if spec==base:chosen.update(name(spec,b,g) for b,g in [(0.,0.),(1.,1.),(1.,0.),(0.,1.)])
            shortlist[key(spec)]=sorted(chosen)
            for row in part:
                if row['name'] in chosen:
                    count.append({**row,**self.score(row['inner'],Path(row['prediction_path']))})
            D.write_json(OUT/'branches'/outer/'native.json',native);D.write_json(OUT/'branches'/outer/'scores.json',count)
        A.freeze(OUT/'branches'/outer/'shortlist.json',shortlist)
        A.freeze(OUT/'branches'/outer/'selection.json',select(count,outer))
        print(json.dumps(dict(event='branch_selection',outer=outer,selection=select(count,outer))),flush=True)

    def capacity(self,outer):
        chosen=A.read(OUT/'branches'/outer/'selection.json');spec=O.TrainSpec(**chosen['spec']);beta=chosen['beta'];gamma=chosen['gamma']
        baseline_rows=[r for r in A.read(OUT/'branches'/outer/'scores.json') if r['name']==chosen['name']]
        all_rows=list(baseline_rows);rank64=[]
        if beta==0 and gamma==0:
            gates=dict(rank64=False,rank128=False,width128=False,reason='Source-inner selection chose baseline-only; no output capacity can change a zero residual')
        else:
            rank64=self.evaluate_candidate(outer,replace(spec,rank=64),beta,gamma);all_rows.extend(rank64)
            order=sorted(set(LINES)-{outer});ordered=lambda rr:[next(r['pds'] if r['valid'] else None for r in rr if r['inner']==l) for l in order]
            rank_gate=promising(ordered(rank64),ordered(baseline_rows))
            reference=O.TrainSpec(steps=2500)
            reference_rows=[r for r in A.read(OUT/'convergence/scores.json') if r['outer']==outer and r['name']==name(reference)]
            width_gate=rank_gate or promising(ordered(baseline_rows),ordered(reference_rows))
            gates=dict(rank64=True,rank128=rank_gate,width128=width_gate,rank64_by_line=ordered(rank64),rank32_by_line=ordered(baseline_rows))
            A.freeze(OUT/'capacity'/outer/'gates.json',gates)
            if rank_gate:all_rows.extend(self.evaluate_candidate(outer,replace(spec,rank=128),beta,gamma))
            if width_gate:
                best=select(all_rows,outer);wide=replace(O.TrainSpec(**best['spec']),hidden=128)
                all_rows.extend(self.evaluate_candidate(outer,wide,best['beta'],best['gamma']))
        A.freeze(OUT/'capacity'/outer/'gates.json',gates)
        D.write_json(OUT/'capacity'/outer/'scores.json',all_rows)
        A.freeze(OUT/'capacity'/outer/'selection.json',select(all_rows,outer))
        print(json.dumps(dict(event='capacity_selection',outer=outer,gates=gates,selection=select(all_rows,outer))),flush=True)

    def outer(self,line):
        convergence=A.read(OUT/'convergence/selection.json')[line]
        branch=A.read(OUT/'branches'/line/'selection.json');final=A.read(OUT/'capacity'/line/'selection.json')
        reference=dict(spec=asdict(O.TrainSpec(steps=2500)),beta=1.,gamma=1.)
        choices=dict(b19_reference=reference,duration=convergence,branches=branch,selected=final)
        A.freeze(OUT/'outer'/line/'frozen_choices.json',choices)
        rows=[];scope=[line];base=self.data.baseline(scope,line,False)
        for seed in (0,1,2):
            for tag,choice in choices.items():
                spec=O.TrainSpec(**choice['spec']);prep,models=self.trained(scope,spec,seed);model=models[spec.steps]
                path=self.prediction(scope,spec,model,prep,line,choice['beta'],choice['gamma'],seed,False)
                rows.append(dict(line=line,seed=seed,tag=tag,name=name(spec,choice['beta'],choice['gamma']),spec=asdict(spec),
                    beta=choice['beta'],gamma=choice['gamma'],**self.score(line,path,seed,False)))
                D.write_json(OUT/'outer'/line/'scores.json',rows)
            for tag,prediction in base.items():
                path=OUT/'outer'/line/'baselines'/f'{tag}.npy';R.save_array(path,prediction)
                rows.append(dict(line=line,seed=seed,tag=tag,**self.score(line,path,seed,False)))
                D.write_json(OUT/'outer'/line/'scores.json',rows)
        controls=[];spec=O.TrainSpec(**final['spec'])
        for mode in ('none','context_shuffle','zero_query','wrong_query'):
            candidate=replace(spec,mode=mode) if mode in ('none','context_shuffle') else spec
            prep,models=self.trained(scope,candidate,0);model=models[candidate.steps]
            query=LINES[(LINES.index(line)+1)%len(LINES)] if mode=='wrong_query' else None
            path=self.prediction(scope,candidate,model,prep,line,final['beta'],final['gamma'],0,False,mode=='zero_query',query)
            controls.append(dict(line=line,seed=0,mode=mode,**self.score(line,path,0,False)))
            D.write_json(OUT/'outer'/line/'controls.json',controls)
        print(json.dumps(dict(event='outer_complete',line=line,records=len(rows),controls=len(controls))),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['freeze','convergence','select_convergence','branches','capacity','outer'],required=True)
    p.add_argument('--scope');p.add_argument('--line',choices=LINES);args=p.parse_args()
    if args.stage=='freeze':freeze_protocol();return
    x=Experiment(args.stage+'_'+str(args.line or args.scope))
    try:
        if args.stage=='convergence':x.convergence(args.scope.split(','))
        elif args.stage=='select_convergence':x.collect_convergence()
        else:getattr(x,args.stage)(args.line)
    finally:x.close()


if __name__=='__main__':main()
