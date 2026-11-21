---
title: "Prophet (universal cellular-phenotype predictor, Theis lab)"
type: entity
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [model, transformer, perturbation-prediction, viability, ic50, transfer-learning, context-count]
sources: ["[[ji-2024-universal-prediction]]"]
---

# Prophet

## Location and authority data

Software, so no site or person authority ID applies. Code: github.com/theislab/prophet (PDF p. 30 of
[[ji-2024-universal-prediction]]). Primary publication: @ji-2024-universal-prediction (bioRxiv
10.1101/2024.08.12.607533; this wiki reads **v3**, posted 23 September 2025). Helmholtz Munich / TUM / LMU, Theis lab.

## What one sees

A supervised transformer that treats any perturbation experiment as a triplet — cell state, perturbation(s),
phenotype — and regresses **one scalar** outcome (PDF pp. 2–3, 33–34). 8 encoder blocks, 8 heads, 512-dim embedding,
20.1 M parameters (Table S3, PDF p. 59). Cell states enter as CCLE bulk-RNA profiles reduced to 300 PCs over 1,406
lines; genetic perturbations as co-expression + Enformer sequence embeddings of the target gene; small molecules as
RDKit2D + Morgan fingerprints; phenotypes as learnable embeddings (PDF pp. 34–35). Pretrained on 4.7 M experiments
(26,000 perturbations, 1,550 cell states), then fine-tuned per dataset (PDF p. 6). Vector phenotypes are handled by
splitting them into per-feature scalar tasks (PDF p. 42).

## History of research

Posted as a preprint in August 2024 (v1); the version in the library is v3 (September 2025), which renumbered the
supplementary figures and deleted the v2 sentence about "the amount of data required to outperform the mean
estimation" — see the v2 → v3 table on [[ji-2024-universal-prediction]].

## The evidence

- Beats linear regression, MLP, Random Forest and the mean estimator on most benchmarks; hit rates +63 % over the
  best simpler architecture (PDF p. 6), and wins 69 % of 64 stratified tasks (PDF p. 9).
- **Against the context-free mean estimator on the unseen-cell-line axis it is much weaker**: at full training it
  wins on GDSC (R² 0.75 vs 0.71) and GDSCcomb (0.35 vs 0.25) and loses on CTRP (0.45 vs 0.47), PRISM (0.59 vs 0.62)
  and SCORE (0.82 vs 0.84) (fig. S3D, PDF pp. 45–46). Spearman with the mean subtracted is ≈ 0 on most
  unseen-cell-line rows (fig. S5D, PDF p. 49).
- Its context-count curves (figs. S6–S7, PDF pp. 50–52) are the evidence base of
  [[context-count-learning-curve]].
- Validated prospectively: 4 of 5 ChEMBL compounds selected in silico reduced viability selectively in melanoma
  lines (PDF pp. 14–15).

## Contested

- **Comparison fairness.** The mean estimator is a closed-form average, not an ablation of Prophet under a matched
  budget, so a crossing of the two curves is a model-class crossing, not a context-pathway break-even
  ([[ji-2024-universal-prediction]], Kernthese 6).
- **Scope of the readout.** Every endpoint is a scalar (IC50, CERES, one L1000 gene, one morphology feature), so
  claims of "transcriptomics" rest on 42–50 L1000 genes (Table S1, PDF p. 57).
- **Noise.** «Prophet is not designed to correct for intra-study experimental noise» (PDF p. 19); where replicates
  are unreliable the authors concede the mean estimator is an upper bound (PDF p. 9).

## Relevance to this project

Named in Amendment A3 of `input/ideas/zero-shot-limits-context-complexity-plan.md` as **the** prior empirical
evidence for the H2a context-count curve. It establishes that such a curve exists and that its crossing is dataset
dependent; it leaves the break-even count, the matched ablation, the transcriptomic readout and all theory open.

## Open / to verify

- Whether a later version (v4+) restores the deleted "amount of data required" sentence or adds the crossings as
  numbers. The crossings recorded here are visual readings of the v3 figure pages.
- Whether the published (peer-reviewed) version keeps the v3 figure numbering.

## Sources

- [[ji-2024-universal-prediction]] — PDF pp. 2–3, 6, 9, 19, 32–37, 42–43, 45–46, 49–52, 57, 59.
