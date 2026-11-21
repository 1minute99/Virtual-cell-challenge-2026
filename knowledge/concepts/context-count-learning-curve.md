---
title: "Context-count learning curve (performance vs number of training contexts)"
type: concept
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [context-count, learning-curve, break-even, mean-baseline, h2a, h2b]
sources: ["[[pan-2026-global-cell-state]]", "[[ji-2024-universal-prediction]]", "[[dibaeinia-2026-context-scale]]", "[[molina-2026-response-decomposition]]", "[[rosenfeld-2021-risks-irm]]"]
relations:
  - target: ji-2024-universal-prediction
    type: builds-on
    confidence: extracted
    because: "The only published version of this curve with a context-free comparator on the same axes: «evaluating Prophet with 80%, 50%, 30%, 20%, 10%, and 5% of the perturbations or cell lines included in the training dataset» (PDF p. 32; figs. S6–S7, PDF pp. 50–52)."
  - target: simple-perturbation-baselines
    type: builds-on
    confidence: inferred
    because: "The curve is only interpretable against a context-free comparator, and the mean estimator / Perturb Mean is that comparator in every instance recorded here."
  - target: rosenfeld-2021-risks-irm
    type: builds-on
    confidence: extracted
    because: "Theoretical floor for the count axis, by analogy only: «the optimal invariant predictor is the global minimum of the IRM objective if and only if E > de» (Corollary 5.2, p. 6), population-level, linear, binary classification."
---

# Context-count learning curve

## Definition

A plot of held-out predictive performance against **the number (or share) of distinct cellular contexts seen in
training**, with the perturbation set and the evaluation protocol held fixed. It is informative only when a
**context-free comparator** is drawn on the same axes: without one, a rising curve shows that more data help, not
that modelling the context helps. The point where the context-using curve crosses the context-free one is the
**break-even context count C\***, the quantity H2a and H2b of
`input/ideas/zero-shot-limits-context-complexity-plan.md` are about.

## Where it comes from

Not a named concept in the literature; this page collects the instances. The project's own formulation (gain of a
context network over a constant-descriptor ablation, as a function of C) is set by Amendment A1 of the plan.

## How it is applied

- **[[ji-2024-universal-prediction]] (Prophet), the prior empirical evidence.** Cell states (and, separately,
  perturbations) are subsampled to 5 / 10 / 20 / 30 / 50 / 80 % of a dataset and the model is plotted against the
  mean estimator (PDF p. 32; figs. S6, S7, PDF pp. 50–52). Crossings read off the figures, in training cell lines:
  GDSC ≈ 170 of 682, LINCS ≈ 6 of 57, and **no crossing** on CTRP (≤ 658), GDSCcomb (≤ 88) or SCORE (≤ 199). With
  pretraining the GDSC and LINCS crossings fall below 5 % and GDSC² crosses at ≈ 11 lines, while SCORE still never
  crosses. The spread of crossings within one paper and one protocol is the central fact: **the break-even is a
  property of the dataset and the pipeline, not a universal constant.**
- **[[dibaeinia-2026-context-scale]]**, as recorded on that page: DEG-F1 rises with the number of training contexts
  per perturbation, but the curve carries no context-free baseline, so it has no crossing.
- **[[molina-2026-response-decomposition]]**, as recorded on that page: one versus three source cell lines for a
  model that is itself context-free (Fig. 3a) — a two-point count manipulation without a context-aware arm.
- **[[rosenfeld-2021-risks-irm]], the theoretical floor (by analogy).** In linear-Gaussian invariant prediction, a
  pooled objective (ERM) cannot be beaten by the invariance objective until the number of environments exceeds the
  environmental dimension, E > d_e (Theorem 5.1, Corollary 5.2, p. 6); below it, d_e + 1 − E environmental directions
  stay unidentified (Lemma C.2, pp. 14–15). Their synthetic curve (d_c = 3, d_e = 6; Fig. C.1, p. 18) has a pooled
  comparator and switches after E = d_e. Population-level (infinite samples, p. 4), classification only, no noise
  term, and no count threshold at all for non-linear models (Theorem 6.1, p. 8).
- **The project**: a seven-point curve C ∈ {6, 8, 12, 16, 24, 32, 40} on Tahoe-100M with a matched constant-descriptor
  ablation, giving C*_emp = 37.4 [33.3, > 40] and a negative gain at C = 12.

## Criticism

- **The comparator decides the answer.** Prophet's mean estimator is a different model class set against a trained
  transformer; a crossing then confounds the value of the context pathway with the capacity of the network. A
  matched ablation (identical data, budget, early stopping, selection, seed) is the only version that isolates the
  context pathway — the reason for the plan's Amendment A1.
- **Shares hide counts.** A curve over a *fraction* of contexts cannot be compared across datasets or carried into a
  design recommendation until it is converted with the dataset's total.
- **The crossing moves with the pipeline.** Pretraining alone moves Prophet's GDSC crossing from ≈ 25 % to below 5 %
  ([[ji-2024-universal-prediction]], figs. S6A vs S7A). Any reported C* must name its protocol.
- **A curve can fail to cross for reasons unrelated to count.** On SCORE the mean estimator is never beaten even at
  ≈ 199 training lines, and the authors attribute this to replicate unreliability (PDF p. 9; Table S4, PDF p. 59).
  Non-crossing is therefore not evidence of a high C* unless reliability is measured too.
- **Readout matters.** Every published instance but the project's own uses a scalar endpoint (viability, IC50,
  CERES, DEG-F1), not a transcriptome-wide residual.

## Not to be confused with

- **Scaling with the number of perturbations** (the other axis of the same figures), which is a different and
  generally easier problem — Prophet beats the mean estimator on the perturbation axis where it fails on the
  cell-line axis ([[ji-2024-universal-prediction]], figs. S6B, S3D).
- **Variance-reduction curves.** Prophet's fig. S5A / S7C show the *standard deviation* of R² falling with the number
  of training conditions (PDF pp. 49, 51). That is precision, not a break-even.
- **The theoretical break-even C\*** of the plan's Proposition 2 (rejected as a quantitative predictor, H2b:
  C*_theory = 230 against C*_emp = 37.4), and the identifiability bound of Proposition 1 (a context model can span
  only rank ≤ C − 1 of the interaction), neither of which has any counterpart in the published curves. Proposition 1's
  nearest published precedent is the E > d_e threshold of [[rosenfeld-2021-risks-irm]] — an analogy (same
  "E − 1 degrees of freedom" counting), not a source of the rank bound, which that paper does not state.

- **Calibration curves over target-context perturbations** (shots, not contexts). [[pan-2026-global-cell-state]]
  sweep K = 5…200 measured perturbations **inside** the held-out context and fit a ridge map from their
  context-free functional embeddings: accuracy rises from the train-mean baseline level at K = 5 to a plateau at
  K ≈ 100 in K562, RPE1 and H1 hESC (pp. 11, 15; Figs. 4D–E, 6G). The number of **contexts** is never varied there —
  it is fixed at the seven Perturb-seq datasets available, and the paper says why: «publicly available Perturb-seq
  collections currently cover a limited number of cell-type contexts» (p. 10). The same within-context saturation
  appears in [[molina-2026-response-decomposition]] at 100–150 target perturbations.

## Relevance to this project

This is the object of RQ2. The literature supplies the shape (performance rises with context count; a context-free
baseline can hold out for a long time; whether it is ever overtaken is dataset dependent) and nothing about the
position. The project's contribution is the position — a measured break-even count with an interval, on a
transcriptomic destination, against a matched ablation.

## Sources

- [[ji-2024-universal-prediction]] — PDF p. 32; figs. S6–S7, pp. 50–52; fig. S3D, pp. 45–46; p. 9; Table S4, p. 59.
- [[dibaeinia-2026-context-scale]] — context curve without a context-free baseline.
- [[molina-2026-response-decomposition]] — 1 vs 3 source lines, context-free model.
- [[rosenfeld-2021-risks-irm]] — Theorem 5.1 / Corollary 5.2, p. 6; Lemma C.2, pp. 14–15; Fig. C.1, p. 18; Theorem 6.1, p. 8.
