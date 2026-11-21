---
title: "RQ2 Task 11 — transfer-pipeline source-count curve and LOFO value: assumptions and checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ2 Task 11: source-count curve PDS(k) and the LOFO family value (descriptive)

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 11 (descriptive, no hypothesis). It extends
`vcc2026_data/scripts/ml_ctx_curve.py`, whose older S10 recipe gave PDS ≈ 0.796 + 0.048 ln k on H1 (AGENTS.md).
It builds on [[rq3_build_preds-assumptions]] (the Task 14 frozen recipe) and feeds [[rq4_h4-assumptions]] (Task 18's per-context value).

- **Script:** `output/data-analysis/rq2_source_curve.py`. Subcommands: `libfam`, `curve P…`, `collect`, `spot_prep`, `spot_build` and `spot_real`.
- **Runners:**
  - `rq2_source_curve_run.sh` runs one subcommand in a scope capped at MemoryHigh 7G / MemoryMax 8G, nice 10, with 4 threads.
  - `rq2_source_curve_chain.sh` runs the panels.
  - `rq2_source_curve_spot.sh` is the real-scorer queue, unit `aivc-paper-t11-spot`.
- **Outputs:**
  - `results/rq2_source_curve.parquet`: one row per panel × subset × arm, 2,523 rows.
  - `results/rq2_source_curve.json`: curves, fits, LOFO values, budget flags, parity checks and the spot check.
  - `results/rq2_source_curve.png`
  - Per-subset arrays (pds_t, num, dist) in `results/rq2_source_curve/sub/<P>/`.
  - Per-family library arrays in `results/rq2_source_curve/libfam/`.
- **Logs:** `results/logs/paper-t11-*.log`.

## Design

The eligible families for a panel are the panel's M1 fusion families after the destination exclusions (`results/panels.json`
`pack_families`). That gives:

| panel | K |
|---|---|
| H1 | 9 |
| OR_HCT116X | 7 |
| KOLF | 9 |
| RN_K562, RN_HepG2, RN_Jurkat, RN_RPE1 | 8 each |

For each k = 1…K:
- If there are ≤ 20 subsets, all of them are used. Otherwise R = 20 distinct subsets are drawn uniformly (rng `[0, panel index, k]`).
- k = K−1 is therefore exactly the leave-one-family-out (LOFO) set.
- 841 subset builds in total.

## Assumptions

1. **M1 + M2 on a subset S is the Task 14 pre-M3 recipe rebuilt on S:**
   - fused(S) = `ml_h1fuse.fixed_fusion` over the 8 min_cells members (from `h1fuse/<P>_data.npz`, the arrays the ERA pack
     was built from), with the weights outside S set to 0;
   - the neighbour library and the library PCs V3 come from both recipe parts fused over S only (exact per-family
     decomposition; the fusion is linear before the clip);
   - M2 = rank(cross-family agreement over S)^0.5;
   - F2 sum side, then the budget match.
2. **S4 is held fixed** at its frozen leave-destination-out fit (KOLF: the A3 kolf-free S4).
   - So the primary curve prices a family through fusion + neighbour library + library PCs + M2, *conditional on S4*.
   - Sensitivity arm `noS4`: M1 without S4, i.e. a pure source-count curve.
3. **M2 is off** (constant amplitude) when no target has ≥ 2 measuring families in S (always at k = 1). This is the only change to the verbatim code.
4. **Budget.** Primary arm `own`: the recipe's own budget Bt(S). Arm `fixBt`: matched to the full-set Bt(K), which is budget-invariant.
5. **Offline scorer** = the ERA pack evaluator, using the subset's own override Od.
6. **Fits** (descriptive): the subsets overlap, so the points are not independent.
   - a + b ln k (OLS) vs a − c e^(−k/τ).
   - Gaussian AIC; |ΔAIC| < 2 counts as indistinguishable.
   - Fitted on the subset-level points (primary) and on the k-means (sensitivity).
7. **LOFO value** = mean over f and t of pds_t(K) − pds_t(K∖f).
   - 95 % CI from B = 2,000 bootstraps that resample families and targets independently (rng 0).
   - Converted to a sum-side composite with the vcc2026 preset constants (PDS and MSE legs only; the DE side is not rebuilt).

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Full-set rebuild vs ERA pack (6 non-KOLF panels): fused, knn, base_total | max abs < 1e-4 | ≤ 1.5e-7 | yes |
| V3 abs-dot / Od / Bt rel | > 0.999 / equal / < 1e-5 | 1.0 / equal / ≤ 1.3e-11 | yes |
| k = K pre-M3 lfc vs Task 14 `<P>_prem3.npz` (7 panels, incl. kolf-free KOLF) | rel < 1e-4, PDS identical | ≤ 4.8e-8, identical PDS and MSE | yes |
| M2 rank vs Task 14; agreement vs verbatim function | equal | equal (7/7) | yes |
| Budget invariance (`own` vs `fixBt`) | — | Bt changes only on H1 (median \|ΔBt\| 7 %, max 69 %) and RN_K562 (median 5 %, max 33 %); kink ≥ 1 elsewhere → Bt = claim, invariant; max \|PDS own − fixBt\| 0.0019 (H1), 0.0009 (RN_K562) | flagged |
| Real scorer at k = K (proxy: Task 18 `noM3` = identical sums, lfc rel 4.6e-8) | \|ΔPDS\| ≤ 0.002 | −0.0003 (MSE exact) | yes |
| Real scorer at k = 1 and k = 5 (H1, median-PDS subset) | \|ΔPDS\| ≤ 0.002 | **queued** (`aivc-paper-t11-spot`, waits for the T18/T10/T17-rest queues) | pending |

## Results (offline PDS, primary arm `own`; fits on subsets)

| panel | K | PDS k=1 → K | b (ln k) | AIC verdict (subsets / k-means) | LOFO ΔPDS [95 % CI] | `noS4` b |
|---|---|---|---|---|---|---|
| H1 | 9 | 0.873 → 0.909 | +0.0172 | indist. (+1.1) / indist. (−0.5) | +0.0021 [−0.0003, +0.0059] | +0.0705 |
| OR_HCT116X | 7 | 0.712 → 0.733 | +0.0129 | indist. / saturating (−3.5; τ at bound 200 on subsets) | +0.0013 [−0.0047, +0.0096] | +0.0711 |
| KOLF | 9 | 0.601 → 0.637 | +0.0170 | indist. / indist. | +0.0022 [−0.0022, +0.0066] | +0.0280 |
| RN_K562 | 8 | 0.894 → 0.906 | +0.0046 | indist. / indist. | +0.0013 [−0.0011, +0.0034] | +0.0560 |
| RN_HepG2 | 8 | 0.898 → 0.908 | +0.0043 | log-linear (+2.0) / indist. | +0.0001 [−0.0049, +0.0053] | +0.0636 |
| RN_Jurkat | 8 | 0.897 → 0.920 | +0.0133 | indist. / indist. | +0.0017 [−0.0027, +0.0079] | +0.0687 |
| RN_RPE1 | 8 | 0.898 → 0.915 | +0.0094 | indist. / saturating (−5.2) | +0.0005 [−0.0015, +0.0036] | +0.0633 |

**With S4 fixed.**
- The curve is shallow: b = 0.004–0.017 per ln k.
- No LOFO CI excludes 0.
- The LOFO value is carried by K562: its per-family ΔPDS is +0.014 on H1, +0.019 on OR_HCT116X, +0.014 on RN_HepG2 and +0.021 on RN_Jurkat. Most other families are ≈ 0 or negative (RPE1 −0.012 on RN_HepG2).

**The S4-free curve is 4–8× steeper.** It reproduces the size of the old ml_ctx_curve slope (b ≈ 0.07 vs 0.048 there). On the
k-means it is mostly saturating (ΔAIC −4.7 to −21).

**H1 per-context value (Task 18 input):**
- ΔPDS +0.0021 [−0.0003, +0.0059];
- sum-side composite +0.00086 [−0.00045, +0.00255];
- PDS-only composite +0.00079.

Task 18's `rq4_sameline.py` hard-codes `per_context_value=None` and rewrites `rq4_h4.json` on every verdict call. The value
is therefore kept in `rq2_source_curve.json` (`per_context_value_H1`) and NOT written into `rq4_h4.json`.

## Deviations (logged)

- **Bisection range widened.** `bmatch` falls back to a log-space bisection on [1e-4, 100] at an edge. This is flagged per row; it never triggered.
- **Zero-coverage subsets in `noS4`.** On OR_HCT116X (HEPG2, RNJURKAT, RPE1 and their combinations; 7 rows) and KOLF (5
  rows), the subset covers no panel target. There `noS4` is the zero delta (overrides only) and is flagged `zero_delta`. The
  first chain run crashed on these (NaN in the PC fraction) and was restarted with this rule; completed builds were kept.
- **Spot-check DE side.** The spot-check packages reuse Task 18's full-set DE delta. Only PDS/MSE (the pinned sums) are
  compared; the DE metrics of those runs are not interpretable.

## Evidence base

- `results/rq2_source_curve.json`, `results/rq2_source_curve.parquet` and `results/rq2_source_curve.png`
- the per-subset npz files and `sub/<P>/parity.json`
- `results/rq4/pkg/noM3/result.json` (the k = K real-scorer proxy)
- the Task 14 preds

## What stays open

- **The k = 1 and k = 5 real-scorer spot checks** (H1, the budget-changing panel) are queued. Until they land:
  - Absolute offline values are trustworthy for the 5 budget-invariant panels and for k = K on H1.
  - H1 and RN_K562 curves at k < K are provisionally *relative*.
- **S4 is not refitted per subset**, so the primary LOFO value is a lower bound on a family's full value to the pipeline.
- **Subset points overlap**, so the AIC comparisons are descriptive only.
