---
title: "Zero-shot transfer ceiling (r_max)"
type: concept
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [ceiling, noise-correction, zero-shot, h1]
sources: ["[[molina-2026-response-decomposition]]", "[[huuloc-2026-cell-quota]]", "[[piao-2026-adaptive-contexts]]", "[[ji-2024-universal-prediction]]"]
relations:
  - target: perturbation-response-decomposition
    type: builds-on
    confidence: inferred
    because: "r_max is defined from Var(a), Var(b) and the reproducible variance of the four-way split (plan, operational definition 6)."
  - target: molina-2026-response-decomposition
    type: builds-on
    confidence: inferred
    because: "Uses Molina & Zhang's components; the source itself never states a ceiling."
---

# Zero-shot transfer ceiling (r_max)

## Definition

The project defines this term; it is not taken from the literature. Operational definition 6 of the plan
(`input/ideas/zero-shot-limits-context-complexity-plan.md`) defines r_max as the noise-corrected correlation achievable
by an oracle that knows a_t and b_c exactly but not i_{c,t}: r_max = √[(Var(a) + Var(b)) / Var_reproducible]. The
attained fraction of a method is r_ours / r_max (hypothesis H1).

## Where it comes from

It combines two existing ideas:

- the component split of the [[perturbation-response-decomposition]] ([[molina-2026-response-decomposition]]), which
  says what is transferable (a, b) and what is not (i);
- reliability ceilings of the form |r| ≤ √ρ (@wang-2026-reliable-perturbations, not yet ingested).

## How it is applied

The project has not computed it yet (Tasks 12–13). A **value implied by Molina & Zhang's numbers** (a reader's
derivation; the source states none): (27.8 + 29.4) / 80.7 ≈ 0.71 of reproducible energy, so r ≈ 0.84 on their
4-line, 638-perturbation tensor (p. 2, Fig. 1c; Supp. Table 1, p. 34). For a centroid-subtracted metric the analogue
on the template-removed residual is ≈ 0.75 (Table 1, p. 34). Details are in the source page's focus item 4.

## Criticism

- Energy shares include μ. They are tensor-level, not the per-perturbation mean correlations most papers report.
- The ceiling inherits the split-half independence assumption.
- Template share varies from 23 % to 65 % across lines (Table 2, p. 35), so r_max is destination-specific. A single
  number misleads.

## Not to be confused with

- **Measurement / reliability ceiling** (√ρ), which limits *any* prediction because of noise alone.
- **Per-condition resolvability** ([[directional-cell-quota]]). [[huuloc-2026-cell-quota]] show that many measured
  responses are not resolved at all (Tahoe 55.4 % control-pool-limited; genome-wide CRISPRi 65–85 % below the
  detection floor; Table 1, p. 8). That lowers Var_reproducible and makes r_max noisy on such destinations; it is not
  itself a prediction ceiling.
- **A reliability ceiling stated in words** ([[ji-2024-universal-prediction]], PDF p. 9). Prophet's authors say that
  where «replicates of the same experimental condition are less similar to each other than to replicates of other
  random conditions … the mean estimator represents an upper bound on predictive performance», and give the
  replicate correlations for SCORE (random 0.481 / 0.674 after aggregation vs within-cell-line 0.615 / 0.820;
  Table S4, PDF p. 59). This is the same intuition as r_max but is never quantified: no reproducible-variance
  denominator, no split-half estimator, no fraction attained, and it is computed on CERES viability scores rather
  than expression. It bounds *the context-free baseline's superiority*, not the achievable correlation.

- **Within-line "upper reference"** ([[molina-2026-response-decomposition]], Fig. 8e, p. 32). This is an empirical
  benchmark, and it is exceeded in one case (target-only 30 % MLP in HepG2, 0.402 > 0.372; Table 22, p. 41).

- **"Upper bound" as a baseline label** ([[piao-2026-adaptive-contexts]], Table 13, p. 22). There the context-free
  Perturb Mean is called «a strong upper bound when perturbation effects are conserved across cell lines». It is an
  empirical predictor (Pearson-∆ 0.650, 3 lines → K562), not a noise-corrected ceiling. It is closer to the
  transferable part a_t that r_max assumes known.

## Relevance to this project

This is the core quantity of RQ1 / H1. The paper's contribution is to make it operational and per-destination, and to
place real methods on it. Molina & Zhang supply the ingredients and an implied value, which must be cited as such.

## Sources

- [[molina-2026-response-decomposition]]
- [[huuloc-2026-cell-quota]] — noise-floor context only (pp. 7, 10–11).
- [[piao-2026-adaptive-contexts]] — heuristic "upper bound" label only (p. 22).
- [[ji-2024-universal-prediction]] — verbal noise ceiling and SCORE replicate correlations (PDF p. 9; Tables S4–S5,
  PDF pp. 59–60).
