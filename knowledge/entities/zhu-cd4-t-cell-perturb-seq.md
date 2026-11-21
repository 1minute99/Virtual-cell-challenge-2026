---
title: "Zhu et al. 2025 — genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells"
type: entity
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [dataset, perturb-seq, crispri, cd4-t-cells, primary-cells]
---

# Zhu et al. 2025 CD4+ T-cell CRISPRi Perturb-seq

## What it is

A genome-scale CRISPRi Perturb-seq dataset of primary human CD4+ T cells (Zhu et al., 2025), as described secondhand
by [[dibaeinia-2026-context-scale]]: 22 million cells, approximately 12,000 single-gene knockdowns, four healthy donors
and three activation timepoints (resting, 8 h, 48 h post-stimulation) (p. 4; App. C.1, p. 13). After the source's
filtering (knockdown ratio < 0.5 in ≥ 2 donors per timepoint; ≥ 100 cells), ≈ 18 M cells and 10,570 perturbations
remain (App. C.3, p. 13). The original Zhu et al. paper is not yet in the library; all facts here are as reported by
Dibaeinia et al.

## Relevance to this project

- Listed in the plan's dataset citations to resolve ("CD4 T-cell genome-scale CRISPRi";
  `input/ideas/zero-shot-limits-context-complexity-plan.md`, Data Sources).
- For C_public (H2c, definition 5): one cell type (primary CD4+ T cells) — counts once under Count 1; its 4 donors are
  donors, not cell lines (A2).
- Its 12 donor × timepoint contexts are the testbed of the Dibaeinia context-coverage analysis (Fig. 4, p. 6).

## Open / to verify

- Bibliographic record of Zhu et al. 2025 and its accession; whether its per-donor responses meet the A2 Count-2 rule.

## Sources

- [[dibaeinia-2026-context-scale]]
