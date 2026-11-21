---
title: "Directional cell quota and control-pool floor"
type: concept
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [sample-size, angular-error, control-pool, sequencing-depth, noise-floor]
sources: ["[[huuloc-2026-cell-quota]]"]
---

# Directional cell quota and control-pool floor

## Definition

Introduced by [[huuloc-2026-cell-quota]]: the number of treated cells n_t* needed so that the direction u = v/‖v‖ of a
perturbation's centroid displacement v = μ_t − μ_c is estimated within an angular tolerance θ* (p. 3):

n_t* = 1 / (m²θ*² / tr(PΣP) − 1/n_c), with P = I − uuᵀ, m = ‖v‖, Σ the per-cell within-condition covariance and n_c the
control-pool size.

Only noise perpendicular to the effect, tr(PΣP), enters (pp. 2–3). The **control-pool floor**
m_min = √(tr(PΣP)/(n_c θ*²)) is the magnitude below which n_t* = ∞: the control's own sampling noise already rotates the
direction by more than θ* (p. 3).

## Where it comes from

It replaces isotropic power calculations and detection-oriented tools (scPower, bulk replicate planning; Hotelling
non-centrality vᵀΣ⁻¹v) that ignore direction and the control arm ([[huuloc-2026-cell-quota]], pp. 2, 14). Its basis
is a first-order delta expansion of the normalisation map, valid for ρ² = m²/(uᵀSu) well above 1 (pp. 2–3, 14).

## How it is applied

Each condition of an atlas is classified as over-sampled, treated-depth-limited, ghost (n* > 50,000), pool-limited
(m < m_min) or not detectable (pp. 4, 7). Results at θ* = 0.1 rad (Table 1, p. 8): Tahoe-100M 55.4 % pool-limited and
11.4 % over-sampled; genome-wide Orion CRISPRi 14.7–35.3 % detectable and almost none over-sampled; TRADE
essential-gene CRISPRi mostly detectable but depth-limited at 48–85 cells per gene (pp. 10–11).

## Criticism

- The noise model is within-condition only. Between-well, plate and cell-cycle variance is excluded, so n* is a lower
  bound where those dominate (p. 14).
- The control-pool term has not been tested on replicate wells (p. 14).
- θ* = 0.1 rad is a convention and the fractions move strongly with it (Table 2, p. 9).
- It certifies a pseudobulk centroid direction in a 50-PC embedding, not gene-level or distributional predictions
  (pp. 3, 14).

## Not to be confused with

- **Control-PC deflation / [[control-subspace-projection]].** "Control-aware" here means the 1/n_c sampling term, not
  projecting responses onto control principal components. P removes the effect direction, not control axes.
- **Reliability / [[zero-shot-transfer-ceiling]].** The quota concerns whether a measured response is resolved, not how
  well a predictor can match it.
- **Context sample complexity (C*).** The quota counts cells per condition; the project's C* counts training cell
  lines.

## Relevance to this project

It is the cells-axis complement to the project's contexts axis. It explains why Tahoe responses referenced to a
plate × line DMSO pool carry shared control noise (see [[rq2_params-assumptions]]) and why genome-wide CRISPRi
destinations have small reproducible signal (H1 ceilings, H4 framing). Overlap per hypothesis is in the source page's
focus block.

## Sources

- [[huuloc-2026-cell-quota]] — pp. 1–4, 7–14; Tables 1–4.
