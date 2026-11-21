---
title: "RQ3 Task 15 — split-half residual-energy estimator (definition 7): assumptions and checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ3 Task 15: split-half residual-energy estimator

Script: `output/data-analysis/rq3_energy.py` (library + CLI for Task 16). Tests:
`output/data-analysis/tests/test_rq3_energy.py` (6 pass; log `output/data-analysis/results/logs/rq3_energy_tests.log`).
Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` — Task 15, operational definition 7, H3a and its
falsification criterion. No amendment (A1–A3) concerns RQ3; A3's truth paths are used via `common.py`.
**No H3a quantity was computed on real data** (Task 16). Upstream: Task 14 (`rq3_build_preds.py`) writes
`results/preds/<panel>_{prem3,m3}.npz` (keys `delta`, `targets`, `genes` = truth var_names) and
`<panel>_ctrl_split_seed<k>.npz` (keys `C_pc`, `C_a`, `C_b`, positions among the control cells in file order, rule
`common.split_halves(n_ctrl, seed, 3)`); this module uses the same rule and the same keys.

## Assumptions stated
- A1 Truth = harness truth file (`which="truth"`); controls = `non-targeting`; split indices = file-order control positions.
- A2 pb(S) = log1p(1e4 · Σ counts / Σ library), library over all panel genes (as `common.pseudobulk` and Task 14's delta).
- A3 G∩ = `results/g_cap.json` (6,014); first panel column per official gene (1 duplicate on H1).
- A4 Δpred aligned to G∩ by official symbol; `--ctrl-ab` uses Task 14's C_a ∪ C_b rebuild (independent of C_pc).
- A5 Three-way split rule `common.split_halves(n_ctrl, seed, 3)`; existing files reused and validated.
- A6 Log-pseudobulk bias O(1/n), common to both halves, negligible (test i).
- A7 Targets are the bootstrap unit; split seeds averaged, not pooled.
- A8 Target halves seeded by (split seed, crc32(target)); targets with < 2 cells dropped.
- A9 Randomized SVD equals exact SVD subspace (test v: min cosine of principal angles > 0.999).

## Checks (synthetic Poisson data, 300 genes, 4 latent control factors)
| Check | Expected | Got | Pass |
|---|---|---|---|
| (i) zero residual: E[numerator], E[denominator] | ≈ 0 | \|z\| < 4 on U_dest, U_j, Haar, denominator | pass |
| (ii) residual on dest control span | F(U_dest) ≫ F(U_perm); D_other CI > 0 | F_dest > 0.6, > 5× F_perm, > 5× F_other; D CIs > 0; F_haar ≈ k/G | pass |
| (iii) isotropic residual, 16 replicates, 100-cell thirds | shared baseline biased > 0; three-way ≈ 0 | shared D_perm +0.052 (z +5.8), three-way −0.054 (z −2.1) | pass |
| (iv) negative control, Δpred = 0 | numerator-unit D ≈ 0 | \|N_other\|, \|N_perm\| < 4 SE; scaled fraction < 50 % | pass |
| (v) P_U projector; rotation invariance; split/loader | exact | exact | pass |
| smoke: H1 controls → U_dest, U_perm (k = 8) on G∩ | orthonormal; U_dest > U_perm variance on held-out C_a | 0.041 vs 0.002 | pass |
| end-to-end `run()` negctrl on RN_RPE1 (2 seeds, temp dir, discarded) | runs; finite | 45 rows, finite | pass |

## Findings that matter for Task 16 (flagged, not acted on)
1. **The target bootstrap misses baseline noise.** The cross-term ⟨P ε_Ca, P ε_Cb⟩ is shared by all targets, so
   resampling targets cannot see it. Synthetic, true D = 0: with 100-cell thirds the spread across replicates is
   2–4× the bootstrap SE. With 1,000-cell thirds it is about 1.4×, and with 4,000-cell thirds the two match
   (0.025 vs 0.030). H1 (12,000 controls, 4,000-cell thirds) is fine, so the registered H3a test is unaffected. RN
   panels (800 controls, 266-cell thirds) will have CIs that are too narrow. Suggested fix: a two-level bootstrap
   (targets × control cells within C_a/C_b) as a labelled sensitivity analysis.
2. **The ratio has a small negative bias with small thirds.** The ratio-form D shows a slight negative bias
   (−0.04, z −2 over 48 replicates, 100-cell thirds; none at 1,000-cell thirds). The numerator itself is unbiased.
   The bias makes D conservative.
3. **The negative-control ratio is ill-posed.** With Δpred = 0 on control cells, both the numerator and the
   denominator of F have expectation 0, so the ratio-form negative-control D is unstable. The module reports the
   ratio form and also the numerator-unit N per target, scaled by the observed per-target denominator
   (`negctrl_vs_observed`). Which form the H3a 50 % rule uses must be decided before Task 16 runs. Log the
   decision as an amendment or deviation.
4. **Δpred depends on the full control mean.** Task 14's Δpred uses the mean of all destination controls (fm),
   including C_pc. The `--ctrl-ab` option uses Task 14's per-seed C_a ∪ C_b rebuild.

## Runtime estimate (Task 16)
One pass per panel covers all k (8/24/50 are nested) and all 5 seeds. Measured on the loaded machine (4 threads):
- Reads: 12,000 controls 4 s; 300 targets about 80 s (from 20 targets in 5 s).
- PCs: U_dest at k = 50 on 4,000 cells 2 s; each permutation 1.4 s (20 × 5 seeds ≈ 2.5 min); other-panel PCs
  about 10 s per panel, cached in `results/rq3_cache/`.
- Total: about 4–5 min per large panel for real mode and about the same for negctrl; under 1 min per RN panel.
  Estimate ≈ 30–60 min for 7 panels × both modes (+ m3 / ctrlAB variants if run).

## Evidence base
Synthetic tests only; H1 controls (smoke); RN_RPE1 controls (negctrl dry run, values discarded).

## What stays open
Items 1 and 3 above need a parent/user decision before Task 16.
