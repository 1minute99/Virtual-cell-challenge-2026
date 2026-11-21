#!/usr/bin/env python3
"""Small, nested Jiang LOCO sweep for eight context-dependent residual models.

The official b16 checkpoint and prior preprocessing/scoring implementations are
read-only references. Each permitted-source denoiser is fixed across models,
hyperparameters, seeds and controls. Inner checkpoints may be reused when the
allowed source set is identical; no outer response enters fitting or selection.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import pathlib
import shutil
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

import context_architectures as A
import jiang_benchmark as J

ROOT = pathlib.Path(__file__).resolve().parents[1]
RANKS = (8, 16)
STRENGTHS = (0.0, 0.25, 0.5, 1.0)
SEEDS = (0, 1, 2)
SOURCE_FILES = ("scripts/context_architectures.py", "scripts/run_context_architectures.py",
                "scripts/context_interaction.py", "scripts/jiang_benchmark.py",
                "tests/test_context_architectures.py")
REFERENCE_FILES = ("scripts/context_interaction.py", "scripts/jiang_benchmark.py",
                   "processed/cross_source_residual/crsres128/model.npz")


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = pathlib.Path(path)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(value,indent=2,allow_nan=True)+"\n")
    temp.replace(path)


def extra_metrics(pred, reps, baseline):
    truth = reps.mean(0)
    def agreement(a,b):
        aa = a-a.mean(1,keepdims=True)
        bb = b-b.mean(1,keepdims=True)
        den = np.linalg.norm(aa,axis=1)*np.linalg.norm(bb,axis=1)
        return float(np.divide(np.sum(aa*bb,axis=1),den,
                     out=np.zeros(len(aa)),where=den>0).mean())
    pn,tn = np.linalg.norm(pred,axis=1),np.linalg.norm(truth,axis=1)
    return dict(row_pearson=agreement(pred,truth),
                pearson_rep1=agreement(pred,reps[0]), pearson_rep2=agreement(pred,reps[1]),
                rmse=float(np.sqrt(np.mean((pred-truth)**2))),
                mae=float(np.mean(np.abs(pred-truth))),
                relative_l2_error=float(np.linalg.norm(pred-truth)/max(np.linalg.norm(truth),1e-12)),
                magnitude_ratio=float(np.mean(pn/np.maximum(tn,1e-12))),
                magnitude_log_error=float(np.mean(np.abs(np.log((pn+1e-12)/(tn+1e-12))))),
                residual_l2_ratio=float(np.linalg.norm(pred-baseline)/max(np.linalg.norm(baseline),1e-12)))


class Experiment:
    def __init__(self,args):
        self.args = args
        self.out = args.out
        self.out.mkdir(parents=True,exist_ok=True)
        self.panels, self.prepared, self.scorers = {}, {}, {}
        self.new_fits = 0
        self.training_seconds = 0.0
        self.start = time.monotonic()
        self.source_hashes = {p:sha(ROOT/p) for p in SOURCE_FILES}
        self.reference_hashes = {p:sha(ROOT/p) for p in REFERENCE_FILES}
        protocol = dict(families=list(A.FAMILIES),ranks=list(RANKS),strengths=list(STRENGTHS),
            seeds=list(SEEDS),training={**A.DEFAULTS,"steps":args.steps},
            stimuli=["IFNB","INS","TGFB","TNFA"],min_cells=25,official_only=True,
            baseline="Unchanged Engine rank128/alpha1/strength1; frozen for each allowed source set",
            official_b16="Frozen checkpoint unchanged; Jiang uses existing permitted-source denoiser",
            target_scope="Same shared Jiang targets and official measured genes; no unseen-target claim",
            context="Training-only control PCA rank2; no cell-line IDs; destination controls only",
            training_views="For each allowed destination: all other sources and one deterministic source removed",
            replicates="Two separate labels, equal-weight mean squared error; source input averages both",
            objective="Residual MSE + 0.5 residual energy + weight decay; pairwise difference auxiliary 0.5; uncertainty NLL auxiliary 0.05",
            centering="Every new correction centered over the target panel, identically across architectures",
            tuning="Inner whole-cell-line PDS only; ties prefer smaller residual strength then rank",
            ablations="Independently selected no-context; same rank/strength no-context; zero query; every wrong source query context",
            shuffles="One fixed training-residual-label shuffle and training-context shuffle; 32 prediction-target derangements",
            stability="Selected context and matched no-context retrained at seeds1,2 without outer retuning",
            source_hashes=self.source_hashes,reference_hashes=self.reference_hashes,
            research=["https://arxiv.org/abs/1609.09106", "https://arxiv.org/abs/1703.06114",
                      "https://proceedings.mlr.press/v80/garnelo18a.html",
                      "https://gaussianprocess.org/gpml/chapters/RW2.pdf"],
            limitations="Exploratory benchmark already studied; six independent cell lines; 23 line-stimulus folds are not 23 independent biologies")
        spec = self.out/"spec.json"
        if spec.exists():
            if json.loads(spec.read_text()) != protocol:
                raise RuntimeError("Existing experiment spec differs; use a fresh output directory")
        else:
            write_json(spec,protocol)
            snapshot = self.out/"sources"
            snapshot.mkdir()
            for p in SOURCE_FILES:
                shutil.copy2(ROOT/p,snapshot/pathlib.Path(p).name)
            write_json(snapshot/"manifest.json",self.source_hashes)
        (self.out/"cache").mkdir(exist_ok=True)
        (self.out/"folds").mkdir(exist_ok=True)
        (self.out/"features").mkdir(exist_ok=True)
        self.state("running")

    def state(self,stage,**more):
        completed = list((self.out/"folds").glob("*/*.metrics.json"))
        write_json(self.out/"run_state.json",dict(stage=stage,elapsed_seconds=time.monotonic()-self.start,
            completed_model_folds=len(completed),expected_model_folds=8*23,
            new_fit_count=self.new_fits,training_seconds_this_invocation=self.training_seconds,**more))

    def panel(self,stim):
        if stim not in self.panels:
            panel = J.load_panel(stim,25,official_only=True)
            self.panels[stim] = panel
            write_json(self.out/f"{stim}.audit.json",panel.audit)
            for line in panel.lines:
                previous = ROOT/"reports/context_interaction_2026-09-15"/f"{stim}.{line}.prediction.npz"
                saved = np.load(previous)
                np.testing.assert_array_equal(saved["targets"],panel.targets)
                np.testing.assert_array_equal(saved["genes"],panel.genes)
        return self.panels[stim]

    def features(self,stim,train,rank):
        key = (stim,tuple(sorted(train)),rank)
        if key not in self.prepared:
            panel = self.panel(stim)
            p = A.prepare(panel.effects,panel.controls,key[1],rank)
            other_key = (stim,key[1],16 if rank==8 else 8)
            if other_key in self.prepared:
                np.testing.assert_array_equal(p.baseline,self.prepared[other_key].baseline)
            self.prepared[key] = p
        return self.prepared[key]

    def fitted(self,stim,p,family,use_context=True,seed=0,shuffle_labels=False,shuffle_contexts=False):
        sources = "".join(map(str,p.train))
        name = f"{stim}.s{sources}.r{p.rank}.{family}.c{int(use_context)}.seed{seed}.label{int(shuffle_labels)}.ctx{int(shuffle_contexts)}.pt"
        path = self.out/"cache"/name
        if path.exists():
            return A.Fitted.load(path)
        model = A.fit(p,family,use_context=use_context,seed=seed,
                      shuffle_labels=shuffle_labels,shuffle_contexts=shuffle_contexts,steps=self.args.steps)
        model.save(path)
        self.new_fits += 1
        self.training_seconds += model.seconds
        return model

    def inner_score(self,stim,p,query_index,correction,strength):
        key = (stim,p.train,p.rank,query_index)
        if key not in self.scorers:
            panel = self.panel(stim)
            self.scorers[key] = A.FastPDS(p,panel.effects[query_index].mean(0))
        return self.scorers[key].score(correction,strength)

    def select(self,stim,outer,family):
        panel = self.panel(stim)
        allowed = tuple(i for i in range(len(panel.lines)) if i != outer)
        rows=[]
        for inner in allowed:
            train = tuple(i for i in allowed if i != inner)
            assert outer not in train and inner not in train
            for rank in RANKS:
                p = self.features(stim,train,rank)
                query = p.query(panel.controls[inner])
                for context in (False,True):
                    fitted = self.fitted(stim,p,family,use_context=context)
                    correction = fitted.predict(query)["mean"][0]
                    for strength in STRENGTHS:
                        rows.append(dict(inner_line=panel.lines[inner],context=context,rank=rank,
                            strength=strength,pds=self.inner_score(stim,p,inner,correction,strength)))
        scores={}
        for r in rows:
            scores.setdefault((r["context"],r["rank"],r["strength"]),[]).append(r["pds"])
        scores={k:float(np.mean(v)) for k,v in scores.items()}
        def choose(context):
            key=max((k for k in scores if k[0]==context),key=lambda k:(scores[k],-k[2],-k[1]))
            return dict(rank=key[1],strength=key[2],inner_pds=scores[key])
        return choose(True),choose(False),rows

    def run_fold(self,stim,outer,family):
        panel=self.panel(stim)
        line=panel.lines[outer]
        folder=self.out/"folds"/family
        folder.mkdir(exist_ok=True)
        prefix=folder/f"{stim}.{line}"
        metrics_path=prefix.with_suffix(prefix.suffix+".metrics.json")
        if metrics_path.exists():
            return
        start=time.monotonic()
        correct_choice,no_choice,inner_rows=self.select(stim,outer,family)
        train=tuple(i for i in range(len(panel.lines)) if i!=outer)
        p=self.features(stim,train,correct_choice["rank"])
        no_p=self.features(stim,train,no_choice["rank"])
        np.testing.assert_array_equal(p.baseline,no_p.baseline)
        previous=np.load(ROOT/"reports/context_interaction_2026-09-15"/f"{stim}.{line}.prediction.npz")
        np.testing.assert_allclose(p.baseline,previous["baseline"],atol=2e-7,rtol=2e-6)
        for pp in (p,no_p):
            feature_path=self.out/"features"/f"{stim}.{line}.r{pp.rank}.npz"
            if not feature_path.exists():
                pp.save_features(feature_path,panel.targets,panel.genes)

        # Outer responses enter only this evaluation section, after selection.
        truth_reps=panel.effects[outer]
        perms=J.derangements(len(panel.targets),seed=0,count=32)
        baseline_pds=float(J.pds(p.baseline,truth_reps.mean(0)).mean())
        mean_pds=float(J.pds(p.original_mean,truth_reps.mean(0)).mean())
        rows,latents={},{}
        target_scores={}
        diagnostics={}

        def record(name,pred,model=None,features=p,strength=0,raw=None,seed=0):
            assert np.isfinite(pred).all()
            m,target=J.metrics(pred,truth_reps,perms)
            m.update(extra_metrics(pred,truth_reps,p.baseline))
            for label,base in (("denoiser",baseline_pds),("mean",mean_pds)):
                den=base-m["shuffle"]
                m[f"gain_vs_{label}_above_own_shuffle"]=(m["pds"]-m["shuffle"])/den if abs(den)>1e-12 else np.nan
            m.update(params=model.params if model else 0,rank=features.rank,strength=float(strength),seed=seed,
                     context_dim=features.dims[2],fit_seconds=model.seconds if model else 0)
            if raw is not None:
                latents[name]=raw["mean"][0]
                if "variance" in raw:
                    # Save the uncertainty of the unscaled latent residual model,
                    # and selected-correction variance separately. This does not
                    # include uncertainty in b16 or omitted gene directions.
                    var=raw["variance"][0].astype(np.float64)*features.scale_y**2
                    latents[name+"__residual_variance_unscaled"]=var
                    latents[name+"__correction_variance_selected"]=var*strength**2
                    actual=(truth_reps-features.baseline[None])@features.basis_y
                    error=actual-raw["mean"][0][None]*features.scale_y
                    m["unshrunk_latent_nll"]=float(0.5*np.mean(np.log(2*np.pi*np.maximum(var,1e-12))[None]+error**2/np.maximum(var,1e-12)[None]))
                    m["unshrunk_latent_90_coverage"]=float(np.mean(np.abs(error)<=1.644854*np.sqrt(var)[None]))
                    m["latent_predictive_variance_mean"]=float(var.mean())
                diagnostics[name]={k:np.asarray(v).tolist() for k,v in raw.items()
                                   if k in ("kernel_parameters","expert_weights","source_weights","posterior_variance","observation_noise")}
            rows[name]=m
            target_scores[name]=target.tolist()

        record("mean",p.original_mean)
        record("denoiser",p.baseline)
        query=p.query(panel.controls[outer])
        correct=self.fitted(stim,p,family)
        raw=correct.predict(query)
        strength=correct_choice["strength"]
        pred=p.decode(raw["mean"][0],strength)
        record("correct_context",pred,correct,p,strength,raw)
        record("correct_unshrunk",p.decode(raw["mean"][0],1),correct,p,1,raw)
        # First real fold verifies the low-rank fast evaluator against full gene scoring.
        exact=A.FastPDS(p,truth_reps.mean(0)).score(raw["mean"][0],strength)
        assert abs(exact-rows["correct_context"]["pds"])<1e-12

        matched=self.fitted(stim,p,family,use_context=False)
        matched_raw=matched.predict(query)
        record("matched_no_context",p.decode(matched_raw["mean"][0],strength),matched,p,strength,matched_raw)
        record("matched_no_context_unshrunk",p.decode(matched_raw["mean"][0],1),matched,p,1,matched_raw)
        selected_no=self.fitted(stim,no_p,family,use_context=False)
        no_raw=selected_no.predict(no_p.query(panel.controls[outer]))
        record("selected_no_context",no_p.decode(no_raw["mean"][0],no_choice["strength"]),
               selected_no,no_p,no_choice["strength"],no_raw)
        zero_raw=correct.predict(p.query(panel.controls[outer],zero_context=True))
        record("zero_query_context",p.decode(zero_raw["mean"][0],strength),correct,p,strength,zero_raw)
        for wrong in train:
            wrong_raw=correct.predict(p.query(panel.controls[wrong]))
            record("wrong_context:"+panel.lines[wrong],p.decode(wrong_raw["mean"][0],strength),
                   correct,p,strength,wrong_raw)
        for name,kwargs in (("training_labels_shuffled",dict(shuffle_labels=True)),
                            ("training_contexts_shuffled",dict(shuffle_contexts=True))):
            null_model=self.fitted(stim,p,family,**kwargs)
            null_raw=null_model.predict(query)
            record(name,p.decode(null_raw["mean"][0],strength),null_model,p,strength,null_raw)
        for seed in SEEDS[1:]:
            for use_context,name in ((True,"correct_context"),(False,"matched_no_context")):
                other=self.fitted(stim,p,family,use_context=use_context,seed=seed)
                other_raw=other.predict(query)
                record(name+f":seed{seed}",p.decode(other_raw["mean"][0],strength),other,p,strength,other_raw,seed)

        np.savez_compressed(str(prefix)+".predictions.npz",prediction=pred.astype(np.float32),
                            targets=np.asarray(panel.targets),genes=np.asarray(panel.genes),**latents)
        write_json(str(prefix)+".diagnostics.json",diagnostics)
        write_json(str(prefix)+".inner.json",inner_rows)
        elapsed=time.monotonic()-start
        result=dict(family=family,stimulus=stim,line=line,train_lines=[panel.lines[i] for i in train],
                    correct_choice=correct_choice,no_context_choice=no_choice,
                    feature_file=f"features/{stim}.{line}.r{p.rank}.npz",metrics=rows,
                    targets=panel.targets,target_pds=target_scores,wall_seconds=elapsed)
        write_json(metrics_path,result)
        print(json.dumps(dict(family=family,stimulus=stim,line=line,
              denoiser=baseline_pds,correct=rows["correct_context"]["pds"],
              no_context=rows["selected_no_context"]["pds"],rank=p.rank,strength=strength,
              params=correct.params,fold_seconds=round(elapsed,2),new_fits=self.new_fits)),flush=True)
        self.state("running",active_family=family,active_fold=f"{stim}/{line}")

    def collect(self):
        rows,selection,targets=[],[],[]
        for path in sorted((self.out/"folds").glob("*/*.metrics.json")):
            a=json.loads(path.read_text())
            common={k:a[k] for k in ("family","stimulus","line")}
            for name,m in a["metrics"].items():
                rows.append(dict(**common,variant=name,**m))
                targets.extend(dict(**common,variant=name,target=t,pds=v)
                               for t,v in zip(a["targets"],a["target_pds"][name]))
            selection.append(dict(**common,**{"correct_"+k:v for k,v in a["correct_choice"].items()},
                **{"no_context_"+k:v for k,v in a["no_context_choice"].items()},wall_seconds=a["wall_seconds"]))
        pd.DataFrame(rows).to_csv(self.out/"fold_metrics.csv",index=False)
        pd.DataFrame(selection).to_csv(self.out/"selection.csv",index=False)
        pd.DataFrame(targets).to_csv(self.out/"target_metrics.csv",index=False)

    def run(self):
        for family in self.args.families:
            for stim in ("IFNB","INS","TGFB","TNFA"):
                panel=self.panel(stim)
                for outer in range(len(panel.lines)):
                    self.run_fold(stim,outer,family)
            self.collect()
            print(json.dumps(dict(family_complete=family,elapsed_seconds=time.monotonic()-self.start)),flush=True)
        self.collect()
        after={p:sha(ROOT/p) for p in REFERENCE_FILES}
        assert after==self.reference_hashes,"Frozen b16 or evaluation reference changed"
        complete=len(list((self.out/"folds").glob("*/*.metrics.json")))==8*23
        self.state("complete" if complete else "partial",frozen_reference_verified=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=pathlib.Path,default=ROOT/"reports/context_architectures_2026-09-15")
    ap.add_argument("--families",nargs="+",choices=A.FAMILIES,default=list(A.FAMILIES))
    ap.add_argument("--steps",type=int,default=160)
    ap.add_argument("--threads",type=int,default=4)
    args=ap.parse_args()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    with threadpool_limits(args.threads):
        experiment=Experiment(args)
        try:
            experiment.run()
        except Exception as exc:
            experiment.state("failed",error=repr(exc))
            raise


if __name__=="__main__":
    main()
