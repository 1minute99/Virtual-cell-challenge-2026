---
title: "RQ2 Task 7 — Tahoe context curve and H2a: assumption checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ2 Task 7: Tahoe context curve and H2a, assumption checks

Scripts: `output/data-analysis/rq2_tahoe_curve.py` (runner, unchanged in this task) and `output/data-analysis/rq2_h2a.py`
(analysis). Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 7, as amended by A1.
Results: `output/data-analysis/results/rq2_h2a.json`, `rq2_h2a_curve.png`. The sweep data are in `rq2_tahoe_curve.parquet`
(confirmatory) and `rq2_tahoe_curve_explor_includezero.parquet` (exploratory). The H2b theory files were not read.

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| S1 grid complete | 7 C × 3 seeds × 6 folds × 8 contexts = 1,008 rows, no NaN, no duplicates | 1,008 / 0 / 0 in both files | pass |
| S2 pairing | same held-out set across C within (seed, fold) | yes | pass |
| S3 nesting | train set at C ⊂ train set at C′ > C; n_train = C | yes | pass |
| S4 ablate constant | ablate_z_maxabs < 1e-5 | 4.5e-8 | pass |
| S5 run independence | runs independent | no: the 3 seeds of a fold share the same 8 held-out contexts | fold-cluster sensitivity reported (verdict unchanged) |
| aggregate() cross-check | point curve = runner `aggregate()` (unweighted) | equal (atol 1e-12) | pass |
| Runner tests | 8/8 | 8/8 pass | pass |
| Analysis reproducible | identical JSON on rerun (bootstrap seed 20261107) | identical | pass |

## Choices made in the analysis (not fixed by the plan text)

- **Zero crossing.** The first upward crossing is used: the first C with g > 0 whose predecessor has g ≤ 0, interpolated
  linearly in ln C. The last crossing (after which g stays > 0) is reported as a sensitivity check. The two are equal for
  the primary and registered curves (37.42 and 37.68), so this choice does not affect the verdict.
- **Percentile CIs.** For C\*, the quantiles use the inverted CDF, with censored draws set to +∞. A lower bound of 33.26
  means that 97.5 % of the draws cross at or above 33.26 or never cross.
- **Rule order.** Rejection is checked first. A case where both the confirm and the reject conditions hold is flagged.
  No such conflict occurred.

## Caveats

- **Stopping step.** The context network almost always restores step 100, the first eligible checkpoint: 94–100 % of
  runs at every C. The ablated network stops later (median 100–300). So at small C the context network moves at most
  100 optimiser steps from zero before inner validation stops it. The primary gain therefore compares two
  early-stopped networks that stop at different points. It is not a comparison of two networks trained to convergence.
- **Weak confound in the ablate arm.** The ablated network is itself worse than the zero residual at small C
  (−1.85e-4 at C = 6, −2.6e-5 at C = 12; CI below 0 up to C = 16). This is the forced-non-zero-selection cost that
  A1 describes. It is about 6 % of the registered gain at C = 12, so the primary and registered curves nearly coincide.
- **No-harm effect of step-0 selection (exploratory).** When step 0 can be selected, both arms choose it in 18 of 18
  runs at C ≤ 12, so the gain there is exactly 0. The context arm chooses step 0 in 72 % of runs at C = 24 and in
  50 % at C = 32, and the gain there is still slightly negative. The crossing is at 36.2, CI [31.5, > 40]. The negative
  gain at C ≤ 12 under the frozen protocol is therefore the cost of a selection rule that cannot decline to use the
  context pathway. Inner validation would decline it if allowed. The break-even C of about 37 is the same under both
  selection rules.
- **Non-monotone curve at small C.** The primary gain is −3.9e-4 at C = 6 and −4.9e-4 at C = 8. Their CIs overlap.
- **Variance.** At C = 12, the between-seed SD of the gain is 9.3e-5 and the between-fold SD is 9.7e-5. The SD of the
  18 run means is 1.8e-4, and the mean within-run SD across contexts is 3.0e-4. Seed-level C\*: 39.3 / 37.0 / 36.7.
- **Fold partition.** The bootstrap does not capture variance from the fixed fold partition (6 folds of 8). With
  fold-level resampling the CIs widen: C\* lower bound 31.1, and 22 % of draws are censored. The verdict is unchanged.
- **Seed 0.** The seed-0 curve (C = 6/12/24/40, context arm only) is reported separately and not pooled. Unweighted
  registered gains are −4.57e-4, −3.88e-4, −1.00e-4 and +1.84e-5, with C\* = 36.96. This agrees with the new seeds.

## Evidence base

- `output/data-analysis/results/rq2_h2a.json` (all numbers)
- `vcc2026_data/reports/drug_context_pretrain_2026-09-20/scaling/k{6,12,24,47}/per_context.json` (seed 0)
- Runner assumptions: `knowledge/synthesis/rq2_tahoe_curve-assumptions.md`

## What stays open

- Whether the paper presents the early-stopping caveat as part of the mechanism (the context pathway overfits within
  100 steps at small C) or as a protocol limitation.
- H2b uses the primary C\*_emp = 37.42, CI [33.26, > 40], with 8.75 % of draws censored. The point estimate is not
  censored, so the "untestable (censored)" clause in A1 does not apply. The H2b verdict is computed independently.
