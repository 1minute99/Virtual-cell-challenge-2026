---
title: "RQ3 Task 14 — frozen pre-M3 / M3 predictions on the 7 harness panels: assumptions and checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ3 Task 14: frozen pre-M3 / M3 predictions

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 14, definition 7 and Amendment A3.
Script: `output/data-analysis/rq3_build_preds.py`, run as the subcommands `feng` → `s4kolf` → `build` → `kolf_nofilm`.
Results are in `output/data-analysis/results/preds/`:
- `<panel>_{prem3,m3}.npz`
- `<panel>_ctrl_split_seed{0..4}.npz`
- `<panel>_prem3_ctrlAB_seed{0..4}.npz`
- `KOLF_{prem3,m3}_noFiLM.npz` and `KOLF_noFiLM.json`
- `rq3_build_preds_summary.json`
- `kolf_fix/` (the kolf-free FENG sources, the KOLF S4, and their JSON reports)

Logs are in `output/data-analysis/results/logs/` (`paper-t14-build.log`, `rq3_build_preds_{feng,s4kolf,kolf_nofilm}.log`).
Seeds: rng 0 throughout, as in the frozen builder; control splits use seeds 0–4. No `cell-eval2` was run.

## Assumptions (stated in the script header)

1. **The "frozen recipe" on a harness panel = the harness M1 plus F2's sum-side settings.**
   - *Harness M1:* the ERA pack `base_total`, built with the destination's own family left out as listed in
     `results/panels.json`.
   - *F2 settings, in builder order:*
     1. ew (1 + cpm/50)^-0.25. F2 has `ew_adaptive false`, so ew is always on.
     2. Centring 0.5.
     3. Guarded removal of the top 3 library PCs (λ 0.75, guard 0.15/0.25).
     4. M2.
     5. M3 (in the M3 arm only).
     6. Budget match to the base budget (kink 1.17, margin 0.8, floor 1).
   - *F2 M1 items that cannot be reproduced on the harness:* the `cens` pseudocount members, the DLD1 family, promoter-off
     fusion, knn_noK, the learned residual 0.25, same-line replacement, and the somatic S4noFK branch.
2. **M2 = F2's `--sum-tamp 0.5`.**
   - The rank of the cross-family agreement (`tv2_tamp_prep.agreement` verbatim) over the panel's fusion families, raised to 0.5.
   - The historical harness M2 (S16/S18 learned rho_hat^0.75) is used only as a reference arm.
   - Spearman between the two rankings: 0.73–0.86.
3. **M3 = builder `sum.ctrl_defl`** (k 8, η 0.75, weaker half by agreement^0.75).
   - V is computed from 6,000 random (rng 0) destination controls, in log1p CP10K on the official axis.
   - Panel genes that map to an already-used official gene are dropped: 1 gene on H1, 1 on OR_HCT116X.
4. **Deliverable delta** = log1p CP10K(fm·e^L) − log1p CP10K(fm).
   - Override entries are applied multiplicatively.
   - No rounding (the noise-free expectation).
   - The gene axis is the truth file's var_names.
5. **The offline scorer (pack moments) is exact here**, because pre-M3 and M3 are matched to the same budget.
   - On RN_HepG2, `sch_sumeval.evaluate` agrees with the pack evaluator to 4 decimals (PDS 0.9298, MSE 0.3191).

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| FENG rebuild with all lines vs the harness stores (fp, fc, ncell, center) | identical | max rel 0 / 6.6e-16 | yes |
| KOLF S4: unmodified refit vs the saved refit_v2 | identical | max abs 0 | yes |
| KOLF M1 rebuild with the original sources and S4 vs the pack | base_total < 1e-4 | 0.0; V3 abs-dot 1.0; Bt rel 2.5e-12; Od equal | yes |
| Recomputed budget Bt vs the pack (6 non-KOLF panels) | rel < 1e-5 | 0 | yes |
| Control mean vs the pack `fm` | equal | equal | yes |
| M3 − pre-M3 PDS > 0 | 7/7 | 7/7 | yes |

## M3 sanity (offline, exact for this budget-matched pair)

| panel | pre-M3 PDS | M3 PDS | ΔPDS | ΔMSE | ref S18-M2 ΔPDS |
|---|---|---|---|---|---|
| H1 | 0.9092 | 0.9104 | +0.0012 | +0.0019 | +0.0008 |
| OR_HCT116X | 0.7328 | 0.7434 | +0.0105 | +0.0008 | +0.0111 |
| KOLF (fixed) | 0.6366 | 0.6459 | +0.0092 | −0.0012 | +0.0018 |
| RN_K562 | 0.9055 | 0.9274 | +0.0219 | +0.0038 | +0.0189 |
| RN_HepG2 | 0.9077 | 0.9298 | +0.0221 | +0.0084 | +0.0198 |
| RN_Jurkat | 0.9205 | 0.9284 | +0.0079 | +0.0065 | +0.0089 |
| RN_RPE1 | 0.9147 | 0.9335 | +0.0187 | +0.0028 | +0.0175 |

- The result agrees with the reference: 7/7 positive.
- MSE rises slightly on 6/7 panels, as expected for a PDS lever.
- Nothing was tuned.

## KOLF leak fix (A3): rebuild

**FENG source stores.** kolf_2 and kolf_3 were removed at cell level, using the per-(target × line) raw-UMI sums.

| arm | perturbed cells, before → after | control cells, before → after |
|---|---|---|
| FENG_T | 1,153,623 → 1,020,967 (−11.5 %) | 8,241 → 7,266 |
| FENG_GWF | 104,130 → 93,856 (−9.9 %) | 218,616 → 193,407 |
| FENG_GWN | 115,124 → 90,014 (−21.8 %) | 281,334 → 236,721 |

- These stores feed every FENG arm: the fusion family column, the first part, the neighbour library, the library PCs and
  the M2 agreement.
- KOLF-panel targets with ≥ 20 cells: GWN 97 → 59, GWF 24 → 24, T 5 → 5.

**S4.** The refit drops the 3,719 pooled Feng rows (xl/feng, fengpooled/fitness) and the pooled empirical stores; the
per-line non-kolf rows are kept.
- The fixed S4 film has cosine 0.818 with the original.
- The base_total has cosine 0.938 with the original.
- The agreement ranks have Spearman 0.986 with the original.

**Remaining exposure.** The saved FiLM residual is still used, with weight 0.25.
- *Sensitivity arm with that weight set to 0* (`KOLF_*_noFiLM.npz`):
  - lfc cosine 0.9992 and rel L2 0.040 against the main build, for both pre-M3 and M3.
  - Pseudobulk delta cosine 0.9985 and rel L2 0.054.
  - PDS: pre-M3 0.6339 (main 0.6366), M3 0.6436 (main 0.6459), ΔPDS +0.0097.

## Control dependence (definition 7)

**The pre-M3 prediction depends on the controls beyond the baseline profile.**
- *Profile uses:* fm (the pseudobulk baseline and the budget match) and the ew cpm.
- *Beyond the profile:* the budget Bt comes from `kink_gamma`'s C_real(n), the delete-1 jackknife trace of n random
  control cells. This is a cell-level variance quantity.
- *How it enters:* it acts as one global scale on the prediction.
- *Source-side parts* (M2, fusion, S4, knn, library PCs) use no destination control cells.

**Size of the effect (seed 0, rebuilt with C_a ∪ C_b only).** Relative change ‖ΔL‖/‖L‖:

| panel | total | through the budget only |
|---|---|---|
| H1 | 0.08 % | — |
| OR_HCT116X | 0.22 % | — |
| KOLF | 0.78 % | 0.76 % |
| RN_K562 | 0.47 % | — |
| RN_HepG2 | 0.46 % | — |
| RN_Jurkat | 1.98 % | 1.87 % |
| RN_RPE1 | 1.21 % | 1.08 % |

**Therefore** `<panel>_prem3_ctrlAB_seed{0..4}.npz` are provided. In them, every control-derived quantity (fm, cpm,
C_real/Bt, budget) comes from C_a ∪ C_b only. They are the primary input for Task 16.

## Evidence base

- `results/preds/rq3_build_preds_summary.json`
- `kolf_fix/{feng_rebuild,s4kolf}.json`
- `KOLF_noFiLM.json`
- The ERA packs (`vcc2026_data/reports/era_2026-10-07/packs`)
- `results/panels.json` (exclusions, leak record)

## What stays open

- The M1 is the harness M1, not the full F2 M1 (Assumption 1).
- The KOLF FiLM residual exposure: bounded by the noFiLM arm, with PDS moving by 0.0027.
- The KOLF S4 loses the non-kolf cells of the pooled genome-wide Feng arms: only per-line rows exist for the targeted arm.
- `TV2_BASE_TABLE` stub used for the S4 step: the base table's y.npy was deleted on 10/11. Parity with the saved
  refit_v2 is exact.
