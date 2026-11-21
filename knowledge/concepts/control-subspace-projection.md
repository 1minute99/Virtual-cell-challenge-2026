---
title: "Control-subspace projection (response energy on control PCs)"
type: concept
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [control-pcs, subspace, h3a]
sources: ["[[molina-2026-response-decomposition]]", "[[huuloc-2026-cell-quota]]"]
---

# Control-subspace projection

## Definition

This is the fraction of a response matrix's energy captured by the top-k principal axes of unperturbed (control)
variation: f_ctrl = ‖X V_k V_kᵀ‖²_F / ‖X‖²_F, with V_k the leading right singular vectors of the control matrix
([[molina-2026-response-decomposition]], p. 25). In the four-line analysis the control matrix is centred
non-targeting control cells of the same line (Supp. Note 3, p. 35). In the genome-wide K562 analysis it is the set of
non-significant perturbations (energy-test p ≥ 0.5), used as a control proxy (p. 25).

## Where it comes from

In [[molina-2026-response-decomposition]] it tests whether perturbation responses use the same axes as homeostatic
cell-state variation (p. 35). Their reciprocal check projects controls onto the perturbation subspace (p. 25).

## How it is applied

The 100 control PCs capture 92–98 % of the template but only 34.5–39.3 % of the template-removed residual. At k = 5
the control PCs capture 7.1–13.4 % of the residual, against 31.2–34.9 % for the residual's own PCs, and random bases
capture about 0.3 % (Tables 6, 7, 10, pp. 36–37; Fig. 1i, p. 5).

## Criticism

- The only null is a random orthonormal basis. There is no other-line control basis and no gene-wise permuted
  controls.
- No split-half noise control is applied to the fractions, so projected energy includes measurement noise.
- Table 7 and Table 10 of the source disagree at k = 5 (10.9 / 7.4 / 10.5 / 11.1 % vs 13.4 / 7.1 / 11.2 / 11.9 %).
- The Methods (§4.9.6) describe the non-significant-perturbation proxy, not the non-targeting-control PCA that
  Supp. Note 3 uses.

## Not to be confused with

- The project's **RQ3 estimator F(U)** (plan definition 7). It applies a three-way control split and a cross-half
  inner product to the **transfer residual** (truth − pre-M3 prediction) and compares destination, other-line and
  permuted bases.
- **The "control-aware" cell quota** ([[directional-cell-quota]], [[huuloc-2026-cell-quota]]). There the control
  enters as sampling variance 1/n_c, not as a subspace. Its one related observation: on chemical atlases drug
  displacements «do not concentrate on the high-variance axes» of the pooled within-condition covariance
  (pp. 13–14).
- **Control expression as model input** (Supp. Note 16, Tables 24, 26). That asks whether descriptors help
  prediction, not where a residual lies.

## Relevance to this project

This is the response-level precedent for H3a, and it must be cited together with the three differences listed above.
The source's reading that control PCs do not capture variance beyond the template (p. 3) must be reconciled with
H3a's claim about transfer residuals.

## Sources

- [[molina-2026-response-decomposition]] — pp. 2, 25, 35–37; Fig. 1i; Fig. 5e (p. 30).
- [[huuloc-2026-cell-quota]] — contrast only; pp. 3, 13–14.
