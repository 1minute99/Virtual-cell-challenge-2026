"""Auditable b19 training trajectories and separately trained residual branches."""
from dataclasses import dataclass, asdict, replace
import copy
import hashlib
import math
import os
from pathlib import Path
import time
import numpy as np
import torch
from torch import nn
import context_forcing_model as M
import simple_residual_data as D


@dataclass(frozen=True)
class TrainSpec:
    rank:int=32
    hidden:int=64
    context_dim:int=16
    steps:int=25000
    batch_size:int=256
    schedule:str='constant'
    horizon:int=25000
    lr:float=.003
    min_lr:float=.00003
    method:str='joint_total'
    global_penalty:float=.5
    context_penalty:float=.5
    weight_decay:float=.03
    mode:str='correct'

    def __post_init__(self):
        if self.method not in ('joint_total','joint_split','staged') or self.schedule not in ('constant','cosine'):
            raise ValueError('Unknown training method/schedule')
        if min(self.rank,self.hidden,self.context_dim,self.steps,self.batch_size,self.horizon)<=0:
            raise ValueError('Training dimensions and budgets must be positive')
        if not 0<self.min_lr<=self.lr or min(self.global_penalty,self.context_penalty,self.weight_decay)<0:
            raise ValueError('Invalid optimization setting')
        self.config()

    def config(self):
        return M.Config(rank=self.rank,hidden=self.hidden,context_dim=self.context_dim,
            steps=self.steps,batch_size=self.batch_size,lr=self.lr,energy=self.global_penalty,
            weight_decay=self.weight_decay,coverage_pass=True,mode=self.mode)


def learning_rate(spec,step):
    if spec.schedule=='constant':return spec.lr
    return spec.min_lr+(spec.lr-spec.min_lr)*.5*(1+math.cos(math.pi*min(max(step,0),spec.horizon)/spec.horizon))


def components(net,features,controls,baseline,zero_context=False):
    gene=net._gene_hidden(features,baseline)
    context=net.encode_context(controls,zero_context=zero_context or net.config.mode=='none')
    interaction=net.interaction_head(net.interaction_body(torch.cat([gene*context,context],dim=-1)))
    return net.global_head(gene),interaction


def branch_hash(net,global_branch):
    h=hashlib.sha256()
    for name,p in net.named_parameters():
        if name.startswith(('gene.','global_head.'))==global_branch:
            h.update(name.encode());h.update(p.detach().cpu().numpy().tobytes())
    return h.hexdigest()


def predict_components(model,features,controls,baseline,zero=False):
    device=next(model.module.parameters()).device
    with torch.no_grad():
        arrays=[torch.as_tensor(np.array(x,copy=True),dtype=torch.float32,device=device) for x in (features,controls,baseline)]
        g,c=components(model.module,*arrays,zero_context=zero)
    return g.cpu().numpy(),c.cpu().numpy()


def fit(prepared,spec,seed=0,checkpoints=None,out=None,device=None):
    start=time.monotonic();config=spec.config()
    features,controls,baseline,target,weights,bios,pairs,ids=M._validate_prepared(prepared,config)
    points=sorted(set(checkpoints or [spec.steps]))
    if points[-1]!=spec.steps or points[0]<=0:raise ValueError('Checkpoints must end at the requested budget')
    if spec.method=='staged' and points!=[spec.steps]:raise ValueError('Staged fits require a single final checkpoint')
    n=len(ids);coverage_batches=math.ceil(n/spec.batch_size)
    split=spec.steps//2
    if spec.method=='staged' and spec.steps-split<coverage_batches:
        raise ValueError('Staged context phase cannot cover all rows within budget')
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
    torch.use_deterministic_algorithms(True);torch.backends.cuda.matmul.allow_tf32=False
    torch.manual_seed(seed)
    if torch.cuda.is_available():torch.cuda.manual_seed_all(seed)
    torch.set_num_threads(2);device=device or ('cuda' if torch.cuda.is_available() else 'cpu')
    mapping={};controls=np.asarray(controls,dtype=np.float32)
    if spec.mode=='none':controls=np.zeros_like(controls)
    elif spec.mode=='context_shuffle':controls,mapping=M._derange_groups(controls,bios,np.random.default_rng(38111+seed))
    net=M.Net(features.shape[1],target.shape[1],controls.shape[1],config).to(device)
    f,c,a,y=[torch.as_tensor(np.array(x,copy=True),dtype=torch.float32,device=device) for x in (features,controls,baseline,target)]
    rng=np.random.default_rng(19317+seed);coverage_rng=np.random.default_rng(829317+seed)
    probe=np.sort(np.random.default_rng(619317+seed).choice(n,size=min(n,512),replace=False));pi=torch.tensor(probe,device=device)
    probabilities=weights/weights.sum();covered=np.zeros(n,bool);gseen=np.zeros(n,bool);cseen=np.zeros(n,bool)
    counts=dict(optimizer_steps=0,coverage_batches=0,coverage_examples=0,sampled_batches=0,sampled_examples=0)
    phase='global' if spec.method=='staged' else 'joint';trace=[];result={};freeze_hash=None
    initial_context_hash=branch_hash(net,False)

    def optimizer_for(which):
        for name,p in net.named_parameters():
            p.grad=None
            is_global=name.startswith(('gene.','global_head.'))
            p.requires_grad_(which=='joint' or (is_global if which=='global' else not is_global))
        return torch.optim.AdamW([p for p in net.parameters() if p.requires_grad],lr=spec.lr,weight_decay=spec.weight_decay)
    optimizer=optimizer_for(phase)

    def update(rows,coverage,clock):
        ix=torch.as_tensor(rows,dtype=torch.long,device=device)
        for group in optimizer.param_groups:group['lr']=learning_rate(spec,clock)
        optimizer.zero_grad(set_to_none=True)
        if spec.method=='joint_total':
            prediction=net(f[ix],c[ix],a[ix]);loss=torch.mean((prediction-y[ix])**2)+spec.global_penalty*torch.mean(prediction**2)
        else:
            g,h=components(net,f[ix],c[ix],a[ix]);prediction=g if phase=='global' else g+h
            loss=torch.mean((prediction-y[ix])**2)+spec.global_penalty*torch.mean(g**2)
            if phase!='global':loss=loss+spec.context_penalty*torch.mean(h**2)
        if not torch.isfinite(loss):raise FloatingPointError('Nonfinite training objective')
        loss.backward();nn.utils.clip_grad_norm_(net.parameters(),5.);optimizer.step()
        covered[rows]=True
        if phase in ('global','joint'):gseen[rows]=True
        if phase in ('context','joint'):cseen[rows]=True
        counts['optimizer_steps']+=1
        key='coverage' if coverage else 'sampled';counts[key+'_batches']+=1;counts[key+'_examples']+=len(rows)

    def measure(step):
        with torch.no_grad():
            g,h=components(net,f[pi],c[pi],a[pi]);prediction=g if phase=='global' else g+h
            mse=torch.mean((prediction-y[pi])**2)
            penalty=(spec.global_penalty*torch.mean(prediction**2) if spec.method=='joint_total' else
                     spec.global_penalty*torch.mean(g**2)+(spec.context_penalty*torch.mean(h**2) if phase!='global' else 0))
            trace.append(dict(step=step,phase=phase,optimizer_steps=counts['optimizer_steps'],lr=optimizer.param_groups[0]['lr'],
                probe_mse=float(mse),probe_penalized_loss=float(mse+penalty),global_rms=float(torch.sqrt(torch.mean(g**2))),
                context_rms=float(torch.sqrt(torch.mean(h**2)))))

    def snapshot(step):
        totals=np.zeros(4,dtype=np.float64)
        with torch.no_grad():
            for lo in range(0,n,4096):
                hi=min(lo+4096,n);g,h=components(net,f[lo:hi],c[lo:hi],a[lo:hi]);prediction=g+h
                main=torch.mean((prediction-y[lo:hi])**2,dim=1)
                penalty=(spec.global_penalty*torch.mean(prediction**2,dim=1) if spec.method=='joint_total' else
                         spec.global_penalty*torch.mean(g**2,dim=1)+spec.context_penalty*torch.mean(h**2,dim=1))
                values=torch.stack([main,main+penalty,torch.mean(g**2,dim=1),torch.mean(h**2,dim=1)],dim=1).cpu().numpy()
                totals+=np.sum(values.astype(np.float64)*probabilities[lo:hi,None],axis=0)
        receipt=dict(**counts,train_spec=asdict(spec),checkpoint_step=step,training_rows=n,
            training_ids_sha256=M._array_sha(ids.astype(np.int64)),unique_rows_seen=int(covered.sum()),
            global_unique_rows_seen=int(gseen.sum()),context_unique_rows_seen=int(cseen.sum()),
            examples_seen=counts['coverage_examples']+counts['sampled_examples'],seed=seed,mode=spec.mode,
            context_shuffle=mapping,coverage_pass=True,trace=list(trace),parameters=sum(p.numel() for p in net.parameters()),
            global_at_freeze_sha256=freeze_hash,global_final_sha256=branch_hash(net,True),
            context_before_sha256=initial_context_hash,context_final_sha256=branch_hash(net,False),
            seconds=time.monotonic()-start,schedule_horizon=spec.horizon,route='residual',family='mlp',
            source_main_mse=float(totals[0]),source_penalized_loss=float(totals[1]),
            source_global_rms=float(np.sqrt(totals[2])),source_context_rms=float(np.sqrt(totals[3])))
        model=M.Fitted(copy.deepcopy(net).cpu().eval(),config,seed,receipt);result[step]=model
        if out is not None:
            folder=Path(out)/f'step{step}';folder.mkdir(parents=True,exist_ok=True)
            if (folder/'model.pt').exists():raise FileExistsError('Refusing to replace a checkpoint')
            model.save(folder/'model.pt');D.write_json(folder/'receipt.json',receipt)
        if spec.method=='staged' and receipt['global_final_sha256']!=freeze_hash:raise AssertionError('Frozen global branch changed')

    measure(0)
    order=coverage_rng.permutation(n)
    for lo in range(0,n,spec.batch_size):update(order[lo:lo+spec.batch_size],True,0)
    measure(0)
    context_order=None
    for step in range(spec.steps):
        if spec.method=='staged' and step==split:
            freeze_hash=branch_hash(net,True);phase='context';optimizer=optimizer_for(phase)
            context_order=coverage_rng.permutation(n)
        in_context_coverage=spec.method=='staged' and split<=step<split+coverage_batches
        if in_context_coverage:
            lo=(step-split)*spec.batch_size;rows=context_order[lo:lo+spec.batch_size]
        else:rows=rng.choice(n,size=spec.batch_size,replace=True,p=probabilities)
        update(rows,in_context_coverage,step)
        if step+1 in points or (step+1)%max(1,spec.steps//10)==0:measure(step+1)
        if step+1 in points:snapshot(step+1)
    return result
