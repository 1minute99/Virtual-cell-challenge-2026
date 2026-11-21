---
title: "RQ2 Task 6 — Tahoe context-curve runner: assumption checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ2 Task 6 — Tahoe context-curve runner: assumption checks

Script: `output/data-analysis/rq2_tahoe_curve.py` (wraps `dcp_drug_lodo` / `dcp_data` / `dcp_model`, unchanged).
Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` Task 6, operational definitions 1–3.

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Fold rule | 6 folds × 8 held-out, pool 40 | `dcp_drug_lodo.folds` (sorted key i % 6): 6 × 8, pool 40 | pass |
| Fold source | `context_splits.json` | that file holds the 5 GENE folds, not the drug folds | deviation (logged) |
| Seed-0 subset draw | reproduce "47 requested → 40 trained" | `dcp_drug_lodo.subset`, evenly spaced, k ≥ 40 returns the whole pool | pass |
| Folds 1–5 reproduce | per-fold gain = scaling/k47 | ≤ 1e-5 relative; no-context MSE < 1e-9 abs; wrong partners 48/48 | pass |
| Fold 0 init recoverable | yes | no: torch 2.14 default generator is seeded randomly per process; original fold-0 init lost | fail (logged) |
| 6-fold mean MSE gain within 10 % | 1.8985e-05 ± 10 % | 1.6429e-05 (−13.5 %) | fail |
| 6-fold mean cos gain within 10 % | 0.001936 ± 10 % | 0.001683 (−13.1 %) | fail |
| Init noise vs tolerance | init SD ≪ 10 % | fold-0 init alone gives SD 2.1e-06 on the 6-fold mean (11 % of ref); ref fold 0 at 90th/95th pct | tolerance tighter than init noise |
| Device | result device-independent | CPU reproduces the original (probably CUDA) folds 1–5 to 1e-5 | pass |
| Budget identical across C | same step cap/patience | steps 4000, eval 100, patience 8; best_step 100 in 22/24 seed-0 runs, 200 in k6 folds 4–5 (earlier text "all 24 → 100" corrected) | pass |
| Replacement acceptance test (logged deviation) | folds 1–5 exact AND ref fold-0 MSE gain in init 5–95 % | resid MSE ≤ 2.7e-6 rel., gains ≤ 3.7e-4 rel.; ref −5.30e-6 ∈ [−4.07e-5, 1.65e-6] | pass |
| Same, cosine gain (secondary, not gated) | inside 5–95 % | ref −0.00185 just above the 95th pct (−0.00195) of 20 inits | note |
| `--arm both` context columns = `--arm context` | equal to 1e-9 | bit-identical (max diff 0.0) at C=12/s1/f0 | pass |
| Default path unchanged by the rework | = cached C=40 s0 rows | fold 1 resid MSE max diff 0.0, best_step equal | pass |
| Ablate arm ignores held-out z | permuting held-out z leaves resid_mse_ablate unchanged | bit-identical; context arm changes | pass |
| Ablate constant = training mean z ≈ 0 | max abs < 1e-5 | yes (float32 centring error) | pass |
| `--select-include-zero` can pick step 0 | best_step 0 when nothing beats zero | C=12/s1/f0: best_step 0, pred = zero residual, 800 steps | pass |
| `aggregate()` unweighted | mean over held-out contexts | equals plain mean; differs from n_pert-weighted | pass |

## Descriptor (for Task 8)

`context_descriptors()` → 48 × 5,955: log1p CP10K pseudobulk of all untreated DMSO_TF cells of each Tahoe line
(plates 3/6/14 pooled), shared 5,955-gene axis, library = axis sum, centred on the Tahoe 48-line mean
(`contexts.npz['profile_adj']`). The network receives its projection onto the top min(64, n_train−1) training PCs, RMS-scaled.

## Evidence base

`vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling/k47/{results,per_context}.json`, `scaling_table.json`;
`output/data-analysis/results/rq2_tahoe_curve_test_C40_seed0.parquet`, `results/rq2_tahoe_fold0_init_check.json`.

## Rework (2026-11-05): ablation arm and selection option

Constant-z ablation (`--arm ablate|both`) and exploratory `--select-include-zero` added; defaults unchanged.
Column change: the network-with-z=0 arm is now `resid_mse_zctx`/`total_cos_zctx`; `resid_mse_zero` now holds the analytic
zero residual (= `resid_mse_nocontext`), as the review requested. Only the test cache carried the old meaning; no Task 7 data exists.
Smoke (C=12, seed 1, fold 0, unweighted mean over 8 held-out contexts): gain vs zero −6.70e-4 (context) / −7.2e-5 (ablate);
ctx_pathway_gain −5.98e-4 (negative in 8/8 contexts); best_step 100 (context) vs 500 (ablate); beat_zero False in both arms.
Single seed/fold — illustrative only, not evidence for H2a.

## What stays open

- Whether the parent accepts the 10 % reproduction gate as failed-for-a-known-reason (lost fold-0 init) or re-specifies it
  as per-fold reproduction of folds 1–5. The hypothesis is not touched.
- Initialisation noise at C = 40 is of the same order as the gain; Task 7's 3 seeds × 6 folds average over it, but
  the bootstrap over held-out contexts in H2a does not capture init variance.
