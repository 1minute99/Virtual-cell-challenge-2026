"""Read-only experiment aggregation and source-inner deployment consensus."""
import argparse,csv,json,math,statistics
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'reports/b19_optimization_2026-09-16'
LINES=('K562','RPE1','HepG2','Jurkat','KOLF')

def read(p):return json.loads(p.read_text())
def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def freeze(p,v):
    if p.exists():
        if read(p)!=v:raise ValueError(f'Frozen value changed: {p}')
    else:write(p,v)
def csvout(p,rows):
    if not rows:return
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def deployment_choice(picks):
    if set(picks)!=set(LINES):raise ValueError('All source-inner choices required')
    groups=defaultdict(list)
    for line,p in picks.items():groups[p['name']].append((line,p))
    ranked=[]
    for name,g in groups.items():
        p=g[0][1]
        ranked.append(dict(name=name,spec=p['spec'],beta=p['beta'],gamma=p['gamma'],occurrences=len(g),source_inner_mean=statistics.mean(x['pds'] for _,x in g),supporting_lines=sorted(l for l,_ in g)))
    return sorted(ranked,key=lambda x:(-x['occurrences'],-x['source_inner_mean'],x['name']))[0]
def promotion(rows):
    want={(l,s,t) for l in LINES for s in range(3) for t in ('b19_reference','selected')}
    subset=[r for r in rows if r['tag'] in ('b19_reference','selected')]
    actual={(r['line'],r['seed'],r['tag']) for r in subset}
    complete=actual==want and len(subset)==len(want)
    valid=complete and all(r['valid'] and r['pds'] is not None and math.isfinite(r['pds']) for r in subset)
    if not valid:return dict(passed=False,complete=complete,all_valid=False,reason='Missing or invalid paired evidence')
    values={(r['line'],r['seed'],r['tag']):r['pds'] for r in subset}
    diffs={(l,s):values[l,s,'selected']-values[l,s,'b19_reference'] for l in LINES for s in range(3)}
    byline={l:statistics.mean(diffs[l,s] for s in range(3)) for l in LINES}
    gain=statistics.mean(byline.values());wins=sum(x>0 for x in byline.values());paired=sum(x>0 for x in diffs.values())
    return dict(passed=gain>=.001 and wins>=3 and paired>=9,complete=True,all_valid=True,mean_pds_gain=gain,line_wins=wins,paired_wins=paired,line_gains=byline,criterion='gain >= .001; >=3/5 line wins; >=9/15 paired wins; all valid')
def main(stage):
    if stage=='freeze_deployment':
        picks={l:read(ROOT/'capacity'/l/'selection.json') for l in LINES}
        freeze(ROOT/'deployment_choice.json',deployment_choice(picks));return
    rows=[];controls=[]
    for l in LINES:
        rows+=read(ROOT/'outer'/l/'scores.json');controls+=read(ROOT/'outer'/l/'controls.json')
    if len(rows)!=120 or len(controls)!=20:raise ValueError('Outer evaluation incomplete')
    gate=promotion(rows);freeze(ROOT/'promotion.json',gate)
    flat=[{k:r.get(k) for k in ('line','seed','tag','name','beta','gamma','pds','valid','null_mean','prediction_path','score_path')} for r in rows]
    csvout(ROOT/'outer_scores.csv',flat);csvout(ROOT/'context_controls.csv',[{k:r.get(k) for k in ('line','seed','mode','pds','valid','null_mean','prediction_path')} for r in controls])
    summary=[]
    for tag in sorted(set(r['tag'] for r in rows)):
        part=[r for r in rows if r['tag']==tag];valid=[r for r in part if r['valid'] and r['pds'] is not None]
        byline={l:statistics.mean(r['pds'] for r in valid if r['line']==l) for l in LINES if any(r['line']==l for r in valid)}
        summary.append(dict(tag=tag,mean_pds=statistics.mean(byline.values()) if len(valid)==15 else None,valid_runs=len(valid),total_runs=len(part),**{l:byline.get(l) for l in LINES}))
    csvout(ROOT/'model_comparison.csv',summary);write(ROOT/'outer_summary.json',dict(comparison=summary,promotion=gate))
    conv=read(ROOT/'convergence/scores.json');cs=[]
    for schedule in ('constant','cosine'):
        for steps in (2500,10000,25000):
            part=[r for r in conv if r['spec']['schedule']==schedule and r['spec']['steps']==steps]
            cs.append(dict(schedule=schedule,steps=steps,mean_pds=statistics.mean(r['pds'] for r in part),n=len(part)))
    csvout(ROOT/'convergence_comparison.csv',cs)
    traces=[]
    for p in ROOT.glob('folds/*/training/*/step*/receipt.json'):
        v=read(p)
        traces.append(dict(scope=p.parents[3].name,checkpoint=str(p),method=v['train_spec']['method'],schedule=v['train_spec']['schedule'],steps=v['checkpoint_step'],rank=v['train_spec']['rank'],source_mse=v['source_main_mse'],source_penalized_loss=v['source_penalized_loss'],optimizer_steps=v['optimizer_steps']))
    csvout(ROOT/'training_diagnostics.csv',traces)
    print(json.dumps(dict(promotion=gate,comparison=summary)),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['freeze_deployment','collect'],required=True);a=p.parse_args();main(a.stage)
