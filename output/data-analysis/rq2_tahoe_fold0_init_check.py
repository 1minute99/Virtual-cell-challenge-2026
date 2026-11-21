"""Task 6 diagnostic: how much does the (unrecoverable) fold-0 initialisation move the seed-0 C=40 gain?

The original seed-0 run initialised fold 0 from torch's per-process random default seed (see
rq2_tahoe_curve.py docstring). Folds 1-5 reproduce exactly. This script reruns fold 0 at C = 40
with the original subset draw and minibatch sampler (Spec.seed 20260920) and 20 different torch
INITIALISATION seeds (20260920 + 1000 + k, k = 0..19), and reports where the original fold-0 gain
lies in that distribution and what 6-fold mean it implies with the reproduced folds 1-5.
Seeds fixed: init seeds listed above; sampler 20260920. Output JSON in results/.
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rq2_tahoe_curve as T
import dcp_data as DD
import pandas as pd

REFDIR = T.ROOT / 'vcc2026_data/reports/drug_context_pretrain_2026-09-20'
d = DD.DrugData('tahoe')
inits = [T.BASE_SEED + 1000 + k for k in range(20)]
out = []
orig = T.init_seed
for s in inits:
    T.init_seed = lambda seed, fold, s=s: s
    rows = T.run_one(d, 40, 0, 0, 'cpu')
    fp = T.fold_pooled(pd.DataFrame(rows)).iloc[0]
    out.append(dict(init_seed=s, mse_gain=float(fp.mse_gain), cos_gain=float(fp.cos_gain),
                    secs=rows[0]['run_secs']))
    print(json.dumps(out[-1]), flush=True)
T.init_seed = orig
ref = {r['fold']: r for r in json.load(open(REFDIR / 'scaling/k47/results.json'))}
r0m = ref[0]['no_context_resid_mse'] - ref[0]['correct_resid_mse']
r0c = ref[0]['correct_total_cos'] - ref[0]['no_context_total_cos']
cache = pd.read_parquet(HERE / 'results/rq2_tahoe_curve_test_C40_seed0.parquet')
fp = T.fold_pooled(cache[(cache.C_requested == 40) & (cache.seed == 0)]).set_index('fold')
rest_m = float(fp.loc[1:5, 'mse_gain'].sum()); rest_c = float(fp.loc[1:5, 'cos_gain'].sum())
m = np.array([o['mse_gain'] for o in out]); c = np.array([o['cos_gain'] for o in out])
summary = dict(
    fold0_reference=dict(mse_gain=r0m, cos_gain=r0c),
    fold0_init_distribution=dict(mse_mean=float(m.mean()), mse_sd=float(m.std(ddof=1)),
                                 mse_min=float(m.min()), mse_max=float(m.max()),
                                 cos_mean=float(c.mean()), cos_sd=float(c.std(ddof=1)),
                                 ref_mse_percentile=float((m < r0m).mean()),
                                 ref_cos_percentile=float((c < r0c).mean())),
    implied_6fold_mean=dict(mse_mean=float((rest_m + m.mean()) / 6),
                            mse_sd_from_fold0_init=float(m.std(ddof=1) / 6),
                            cos_mean=float((rest_c + c.mean()) / 6),
                            cos_sd_from_fold0_init=float(c.std(ddof=1) / 6)),
    runs=out)
json.dump(summary, open(HERE / 'results/rq2_tahoe_fold0_init_check.json', 'w'), indent=2)
print(json.dumps({k: v for k, v in summary.items() if k != 'runs'}, indent=2))
