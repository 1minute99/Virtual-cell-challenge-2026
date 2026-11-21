---
title: "Systematic variation (and the perturbed-centroid reference)"
type: concept
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
tags: [evaluation-metrics, systematic-variation, perturbed-reference, centroid-accuracy, benchmark, h3b]
sources: ["[[vinastorne-2026-systema]]", "[[molina-2026-response-decomposition]]"]
relations:
  - target: vinastorne-2026-systema
    type: builds-on
    confidence: extracted
    because: "Concept, measure and perturbed-centroid reference are defined there (pp. 1050, 1054–1055; PDF p. 11)."
---

# Systematic variation

## Definition

Coined by [[vinastorne-2026-systema]]: «the consistent transcriptional differences between perturbed and control
cells arising from selection biases or confounders» (p. 1050). It is a shift that most perturbations of a dataset
share relative to the controls. Its named causes are a functionally narrow perturbation panel, unmeasured covariates
(cell-cycle phase, chromatin, target efficiency) and responses that are «biological in origin but systematic in
effect», such as stress or cell-cycle arrest (p. 1051).

**Measure.** For each perturbation X, take the shift from the control centroid to the centroid of X-cells, and take
the average perturbation effect as the shift from the control centroid to the centroid of all perturbed cells. The
dataset score is the cosine between the two, averaged over perturbations (p. 1054; Fig. 3a).

**Perturbed-centroid reference (the "Systema" fix).** Replace the control centroid in Δ_X and Δ̂_X by
O_pert = (1/|𝒫_train|) Σ_{X∈𝒫_train} O(X), the unweighted mean of the training perturbations' ground-truth centroids
(PDF p. 11). The same vector is subtracted from truth and prediction before the Pearson correlation is computed.

**Centroid accuracy.** The share of other perturbations Y whose true centroid lies farther (Euclidean) from the
predicted centroid of X than X's own true centroid does (PDF p. 12).

## Where it comes from

It replaced the GEARS/scGPT practice of control-referenced PearsonΔ / PearsonΔ20, after a perturbed-mean baseline was
found to match or beat GEARS, scGPT and CPA on ten datasets (p. 1051).

## How it is applied

- **Ten datasets.** Systematic variation is high in Adamson and Norman, higher in Replogle RPE1 than in K562, and lowest
  in Frangieh (Fig. 3b, p. 1054). It correlates 0.91–0.95 with methods' PearsonΔ across datasets (p. 1055). Under
  O_pert the mean cosine falls from 0.41 to −0.06, all scores drop, and the correlation with systematic variation
  falls to 0.21 (p. 1055).
- **Zero-shot transfer.** [[molina-2026-response-decomposition]] use perturbed-reference Pearson as their headline
  metric, but for zero-shot they take the *observed target* centroid as the reference (their p. 15), because no
  destination training perturbations exist.
- **This project (H3b).** The "Systema-style" arm subtracts η × a target mean from predictions on the weak half
  (log 2026-11-05, Task 17). That is a prediction-side analogue, not the metric-side reference.

- [[li-2024-systematic-comparison]] name the same object as a decomposition term: G_c, "the systematic shift shared
  across perturbations within the same cell state" (p. 30), and recommend centred metrics (p. 12); their Pearson Δ
  only subtracts the control (p. 25), with no Systema-style correction.

## Criticism

- From the source itself (PDF p. 12): reference-sensitive metrics are noisy for weak perturbations, invariant to
  strength, and depend on, or are even undefined for, the chosen reference.
- Two different means are in play: the cosine measure uses the cell-weighted mean of all perturbed cells (p. 1054),
  while O_pert is perturbation-weighted (PDF p. 11).
- *Reader's caveat:* O_pert removes only a rank-one common translation. A shared response that varies in amplitude
  across perturbations, such as a template scaled per perturbation, survives it.

## Not to be confused with

- **The perturbed-mean baseline** μ_pert ([[simple-perturbation-baselines]]): a *predictor*, cell-weighted, not a
  reference.
- **The template / μ + α_c** of [[perturbation-response-decomposition]]: a per-line component estimated across lines,
  with a cell-line × perturbation interaction beyond it.
- **[[control-subspace-projection]] / M3 control-PC deflation**: these remove a k-dimensional subspace estimated from
  control cells, not one centroid.

## Relevance to this project

Source of H3b's fourth alternative. It also explains why control-referenced correlations (H1) should be read against
what a template-only predictor scores, and why H2a's residual MSE (reference-insensitive) avoids this inflation.

## Sources

- [[vinastorne-2026-systema]]: pp. 1050–1051, 1054–1055; PDF pp. 11–12; Figs. 3–5.
- [[molina-2026-response-decomposition]]: zero-shot adaptation of the perturbed reference (p. 15 there).
- [[li-2024-systematic-comparison]]: pp. 12, 25, 30.
