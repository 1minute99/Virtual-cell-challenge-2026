---
title: "Perturbation response decomposition (μ + α_c + β_p + γ)"
type: concept
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
tags: [variance-decomposition, split-half, crispri, notation]
sources: ["[[pan-2026-global-cell-state]]", "[[molina-2026-response-decomposition]]", "[[ahlmanneltze-2025-linear-baselines]]"]
---

# Perturbation response decomposition

## Definition

[[molina-2026-response-decomposition]] write the pseudobulk response of perturbation p in cell line c at gene g as
δ(c,p,g) = μ(g) + α_c(g) + β_p(g) + γ(c,p,g). The four terms are a global mean, a cell-line effect, a perturbation
effect conserved across lines, and a line × perturbation interaction. They are computed as row, column and
interaction means of a **balanced** tensor (4 lines × 638 perturbations × 2,000 genes), with zero-sum constraints, so
the four component tensors are mutually orthogonal and their sums of squares add up (pp. 20). μ + α_c is called the
**cell-line template** (p. 21).

**Noise correction (split-half).** Cells in each condition are split into two halves and the decomposition is
recomputed on each. The reproducible energy of a component is the cross-half Frobenius inner product ⟨Z⁽¹⁾, Z⁽²⁾⟩, and
noise is the mean half-energy minus total reproducible signal. Results are averaged over 50 splits (pp. 21–22, 28).

**Project notation.** μ, α_c, β_p and γ correspond to the plan's μ, b_c, a_t and i_{c,t}; ε is the plan's
measurement-noise term (plan `input/ideas/zero-shot-limits-context-complexity-plan.md`, "Notation").

## Where it comes from

The source presents it as a diagnostic that shows *which* response components a model recovers, as an alternative to
aggregate scores that reward the template (pp. 1–2, 9). Algebraically it is the classical additive two-way
(row / column / interaction) decomposition, applied per gene. The source cites no statistical precedent for it.

## How it is applied

- On the observed Replogle–Nadig data ([[replogle-nadig-crispri-collection]]) the full-response energy splits into
  template 27.8 %, conserved 29.4 %, interaction 23.5 % and noise 19.3 % (p. 2, Fig. 1c). After projective template
  removal the split is conserved 35.0 %, interaction 27.7 % and noise 37.4 % (Supp. Table 1, p. 34).
- To attribute predictions, the source correlates δ̂ with the true β_p and γ (r_β, r_γ; p. 27). Zero-shot models
  align with β but never with γ (Table 25, p. 43).

- A parallel, *diagnostic* decomposition in [[li-2024-systematic-comparison]] (preprint of
  @li-2026-systematic-comparison): observed pseudobulk profiles Y_{c,p} = B_c + G_c + S_p + I_{p1,p2} + C_{c,p} + ε,
  with B_c the control baseline, G_c the shift shared across perturbations within a cell state, S_p the
  perturbation-specific response, I_{p1,p2} the non-additive term for combinations, C_{c,p} the cell-state-by-
  perturbation interaction, and ε an "empirical noise floor" from split halves of control cells (pp. 30–31). Unlike
  Molina & Zhang's it is not an orthogonal partition (fractions "not constrained to sum to 100%", p. 31) and the
  components are not noise-corrected; it serves to explain why metric families disagree (p. 11; Fig. 5).

## Criticism

- The authors limit its scope themselves: it is «a diagnostic description of response variation rather than a complete
  mechanistic model» on four essential-enriched pseudobulk screens (p. 9).
- *Reader's caveats (not from the source):* the shares are fractions of uncentred energy, and the template term
  includes μ. The split-half estimator assumes the two halves have independent noise (p. 21); random cell halves can
  share clones (see @huang-2026-shared-clones, not yet ingested). The balanced design keeps only perturbations present
  in all lines (638), which biases the sample towards shared essential genes.

## Not to be confused with

- **Within-line template projection** (δ = τ·T_c + ε⊥; pp. 19–20, 24). This projects each response on a single
  per-line mean direction. It does not use the cross-line α / β / γ split.
- **Perturbed-reference (centroid-subtracted) metrics** in the Systema style ([[vinastorne-2026-systema]]; [[systematic-variation]]). These subtract
  a centroid inside the metric but do not decompose anything.

- **The "mean" baseline** of [[ahlmanneltze-2025-linear-baselines]] ([[simple-perturbation-baselines]]): the
  average response over a line's training perturbations (p. 1660), i.e. an estimate of the template μ + α_c computed
  within one line. It predicts no perturbation-specific part, yet no deep model consistently beats it (p. 1660). The
  mapping to the template is a reader's inference; that source decomposes nothing.

- **A multiplicative factorisation of the same object** ([[pan-2026-global-cell-state]], p. 10; Methods pp. 57–58):
  Δŷ = Σ_k π_k(z) W_k v_g — a **context-independent** functional embedding v_g of the perturbed gene times a
  **cell-state-conditioned** mixture of K = 8 linear maps, with v_g ∈ ℝ¹⁶ and 18,108 readout genes. It separates
  the same two things (what the perturbation is vs how a context reads it), but as a model, not as a variance
  accounting: no component shares, no noise term and no split-half estimator appear in that paper, so none of its
  numbers can be set next to the 27.8 / 29.4 / 23.5 / 19.3 % shares above.

## Relevance to this project

This is the plan's notation, and its split-half cross-product is operational definition 4 (Var(i) = ⟨i_A, i_B⟩). Its
components are the ingredients of the [[zero-shot-transfer-ceiling]] (H1). The project applies it per destination
panel and to transfer residuals, which the source does not do.

## Sources

- [[molina-2026-response-decomposition]] — pp. 2, 19–22, 34; Fig. 1a, 1c.
- [[ahlmanneltze-2025-linear-baselines]] — mean baseline as within-line template analogue (inferred); p. 1660.
- [[li-2024-systematic-comparison]]: pp. 11, 30–31; Fig. 5.
