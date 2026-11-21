---
title: "Simple perturbation baselines (no change, additive, mean, linear GWPᵀ + b)"
type: concept
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
tags: [baseline, linear-model, mean-baseline, additive-model, benchmark, h4]
sources: ["[[pan-2026-global-cell-state]]", "[[ahlmanneltze-2025-linear-baselines]]", "[[piao-2026-adaptive-contexts]]", "[[ji-2024-universal-prediction]]", "[[vinastorne-2026-systema]]"]
relations:
  - target: ahlmanneltze-2025-linear-baselines
    type: builds-on
    confidence: extracted
    because: "The linear model W = argmin ||Y_train − (GWPᵀ + b)||² and the mean baseline b are defined there (eq. 1, p. 1660; eq. 3, PDF p. 6)."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "The mean baseline b is the within-line average response, i.e. an estimate of the template μ + α_c; reader's mapping, not stated in either source."
---

# Simple perturbation baselines

## Definition

The family of deliberately simple predictors against which deep perturbation models are benchmarked. As defined by
[[ahlmanneltze-2025-linear-baselines]]:

- **No change**: predict the control expression y∅ (PDF p. 6).
- **Additive** (double perturbations): ŷ_add = y_A + y_B − y∅, the sum of single-perturbation log fold changes (eq. 2,
  PDF p. 6).
- **Mean**: b = the mean expression across the training perturbations of the same line (p. 1660; PDF p. 6), taken from
  Kernfeld et al. and Csendes et al. (p. 1660).
- **Linear model**: Y_train ≈ G W Pᵀ + b, with G = top K = 10 PCs of Y_train (genes × PCs), P = rows of G for the
  perturbed genes, W from the ridge normal equations with λ = 0.1 (eqs. 1, 3; p. 1660; PDF pp. 6–7). G or P can be swapped for
  external embeddings. The additive model is a special case (G = Y_single, binary P, W = I, b = −y∅; PDF p. 7).

## Where it comes from

It is a reaction to foundation-model claims (scGPT, scFoundation, GEARS) that pretrained representations predict
unseen perturbations. The authors attribute the idea to a reviewer of an earlier paper who suggested comparing
foundation models with simple linear models (PDF p. 7, Acknowledgements). The mean baseline was added during revision
after two preprints proposed it (p. 1660).

## How it is applied

On Norman (K562 CRISPRa doubles), Replogle K562 / RPE1 and Adamson (CRISPRi singles), no deep model consistently beats
these baselines (pp. 1657, 1660; Figs. 1a, 2a). The only consistent improvement is the linear model with a perturbation
embedding P pretrained on the other Replogle line (Fig. 2c, p. 1659; p. 1660). The mean baseline alone reaches Pearson
delta ≈ 0.40 / 0.68 / 0.73 on K562 / RPE1 / Adamson (ED Fig. 8a, PDF p. 15; read from plot).

**Across lines** ([[piao-2026-adaptive-contexts]]). A "Perturb Mean", the per-perturbation average response over the
source lines, is the cross-line analogue of the mean baseline. Training on RPE1 + HepG2 + Jurkat and testing on K562,
it scores Pearson-∆ 0.650 and DE Spearman LFC 0.817. Every learned model is lower: Linear-GenePT 0.515, AdaPert 0.460,
TxPert 0.405, STATE 0.309, MORPH −0.113 (Table 13, p. 22; Fig. 13, p. 21). Within lines, the source's "Context Mean"
(same numbers as its "Matching-Mean") gives 0.387 (K562) and 0.629 (RPE1), against AdaPert 0.620 / 0.676 (Tables 9, 12,
pp. 21–22).

**Across unseen cell lines, with a context-count curve** ([[ji-2024-universal-prediction]]). Prophet's "mean
estimator" is the same estimator on the zero-shot axis, and the source defines it exactly: «The mean baseline
estimates that the effect of an perturbation in an unseen cell state is the average value of that perturbation in all
previously seen states. Conversely, it estimates the effect of an unseen perturbation on a cell state is the average
of all seen perturbations in that same cell state» (fig. S1C legend, PDF p. 43; "an perturbation" sic). The figure
panel glosses the first case as "Mean baseline = assumes drug behaves the same across all cell lines". The main text
calls it «a difficult baseline in drug response prediction» (PDF p. 6), and the paper shows how difficult: on the
unseen-cell-line axis at full training it beats Prophet on CTRP (R² 0.47 vs 0.45), PRISM (0.62 vs 0.59) and SCORE
(0.84 vs 0.82), losing only on GDSC and GDSCcomb (fig. S3D, PDF pp. 45–46). Plotted against the share of training
cell lines it is not overtaken at all on CTRP, GDSCcomb or SCORE up to 80 % of lines (figs. S6–S7, PDF pp. 50–52; see
[[context-count-learning-curve]]). Where replicates are unreliable the authors make it a ceiling: «the mean estimator
represents an upper bound on predictive performance» (PDF p. 9).

**Ten datasets, control vs perturbed reference** ([[vinastorne-2026-systema]]). The "perturbed mean" (μ_pert, the
mean of all perturbed training cells, cell-weighted) and the "matching mean" for pairs (½(μ(X) + μ(Y)), with unseen
singles replaced by μ_pert; PDF p. 11) match or beat CPA, GEARS and scGPT on control-referenced PearsonΔ on all ten
datasets (e.g. 0.70 Adamson, 0.55 Replogle RPE1, 0.32 Replogle K562; Fig. 1c, p. 1052). The matching mean beats GEARS
by 11 % on unseen Norman pairs (p. 1051). The source explains the result by [[systematic-variation]]. Under its
perturbed-centroid reference the perturbed mean falls to near zero (PearsonΔ −0.11 to 0.15) and every method scores low (Fig. 4c, p. 1056).

**A context-conditioned model losing to the train-mean, in a genetic cross-context setting**
([[pan-2026-global-cell-state]]). Their mixture-of-experts maps a context-free functional embedding of the perturbed
gene through a cell-state-gated linear map (p. 10). Trained on seven Perturb-seq knockdown datasets and evaluated on
a held-out context (K562, RPE1 or H1 hESC), its **zero-shot predictions score below the train-mean baseline** — "the
mean value of shift vectors with the same perturbation across the training dataset", the baseline the authors take
from @ahlmanneltze-2025-linear-baselines (p. 11). Medians read off the violins: RPE1 ≈ 0.00 (zero-shot) vs ≈ 0.15
(baseline); K562 ≈ 0.04 vs ≈ 0.10; hESC ≈ 0.03 vs ≈ 0.05 (Figs. 4D–E p. 37, Fig. 6G p. 40). Fitting a ridge model on
as few as **K = 10 measured target-context perturbations** flips the ordering (RPE1 ≈ 0.44, hESC ≈ 0.09), with a
plateau at K ≈ 100 (pp. 11, 15). The metric is per-perturbation Pearson on pseudobulk shift vectors with no
centroid subtraction, so the numbers are not comparable with the perturbed-reference scores above.

- [[li-2024-systematic-comparison]]: AdditiveMean and Linear beat the complex models on seen combinations (seen1/seen2;
  p. 9), and no method generalises consistently to unseen cell types (p. 2).
- [[adduri-2025-state]] claim STATE beats the baselines, but on the genetic Replogle–Nadig data it only *matches* the
  perturbation mean on discrimination (−3 %; Fig. 2D, p. 8; «State matched the performance of the perturbation mean
  baseline», p. 10), with the context mean the strongest baseline on Pearson Δ (+5 %, p. 10; Fig. 2E).

## Criticism

- The Perturb Mean needs the perturbation to be measured in the source lines, so it is zero-shot in the cell line but
  not in the perturbation ([[piao-2026-adaptive-contexts]], p. 21). The same holds for Prophet's mean estimator on
  the unseen-cell-state axis ([[ji-2024-universal-prediction]], PDF p. 43).
- A mean baseline is a *different estimator class*, not an ablation. When it is set against a trained network, a
  crossing of the two curves confounds the value of the context pathway with the capacity of the model — the reason
  for Amendment A1 of `input/ideas/zero-shot-limits-context-complexity-plan.md`
  ([[ji-2024-universal-prediction]], PDF pp. 6, 31).
- *Reader's caveats (not from the source):* all four baselines except "no change" need measured perturbations of the
  destination line, so they define a within-line, unseen-perturbation benchmark. In an unseen-cell-line (zero-shot)
  setting the mean and linear baselines do not exist as defined; their analogue is a transferred template.
- Control-referenced Pearson delta rewards the mean baseline for recovering the shared response; how much of the
  baselines' score is "template" is not quantified in the source.
- The source's own limits: four datasets, all cancer lines, original QC kept (p. 1660).

## Not to be confused with

- **The project's transfer baseline** (plan H4, `input/ideas/zero-shot-limits-context-complexity-plan.md`): the
  prediction of an architecture *before* its modelling addition, in the unseen-line setting. It is not the
  within-line mean.
- **Systema-style common-response removal** ([[vinastorne-2026-systema]], [[systematic-variation]]): subtracts a centroid inside the metric or as
  a correction, rather than predicting it.

## Relevance to this project

The literature precedent for H4's "data beats modelling" framing, and the reason any ceiling-based evaluation must
state how much a perturbation-agnostic template already scores (see [[zero-shot-transfer-ceiling]]).

## Sources

- [[ahlmanneltze-2025-linear-baselines]] — eqs. 1–3, pp. 1657, 1660; PDF pp. 6–7; Figs. 1–2; ED Fig. 8.
- [[piao-2026-adaptive-contexts]] — cross-line Perturb Mean (Fig. 13, p. 21; Table 13, p. 22); within-line Context
  Mean (Tables 9–12, pp. 21–22).
- [[ji-2024-universal-prediction]] — the mean estimator on the unseen-cell-line axis (fig. S1C legend, PDF p. 43;
  text PDF pp. 6, 9; fig. S3D, PDF pp. 45–46; figs. S6–S7, PDF pp. 50–52).
- [[vinastorne-2026-systema]] — perturbed mean / matching mean on ten datasets (PDF p. 11; p. 1051; Figs. 1c, 4c).
- [[li-2024-systematic-comparison]]: pp. 2, 9.
- [[adduri-2025-state]]: pp. 8, 10; Fig. 2D–J.
