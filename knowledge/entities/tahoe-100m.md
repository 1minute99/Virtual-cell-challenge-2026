---
title: "Tahoe-100M (chemical single-cell perturbation atlas)"
type: entity
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
tags: [dataset, chemical-perturbation, tahoe]
sources: ["[[molina-2026-response-decomposition]]", "[[huuloc-2026-cell-quota]]", "[[huang-2025-x-atlas-orion]]"]
---

# Tahoe-100M

## Location and authority data

This is a dataset, so no authority identifier applies. Primary publication: @zhang-2026-tahoe (not yet ingested).

## What one sees

A giga-scale chemical perturbation atlas. Per the literaturguide entry it has about 50 lines and 1,100 drug-dose
conditions. The plan uses 48 Tahoe contexts for H2a/H2b. [[molina-2026-response-decomposition]] used only **5 lines**
(C32, HOP62, HepG2/C3A, Hs 766T, PANC-1), 268 compounds and 3 doses. Of 4,020 line × compound × dose conditions,
3,662 have ≥ 20 cells. Pooled DMSO serves as the control, without plate or batch matching (p. 17).

[[huuloc-2026-cell-quota]] give the full design (p. 3): 100,648,790 cells (95,624,334 after QC), 379 small
molecules × 3 doses, **50** cancer lines, 14 96-well plates (each well pools all lines under one treatment); 56,827
drug × dose × line conditions with a median of 1,296 cells; 2–3 DMSO_TF vehicle wells per plate, pooled to one vehicle
centroid per plate × line (median about 3,113 cells, serving a median of 94 conditions; p. 4). 86.3 % of drug–dose
combinations occupy one plate, so replicates are largely absent (p. 10).

## History of research

In [[molina-2026-response-decomposition]], one-shot transfer of the same compound's response from other lines (low-rank
Ridge, mean r = 0.86) far beats static chemical descriptors (Morgan r = 0.11). The result survives template centring
and collapses under compound-label permutation (pp. 3–4, 17–18; Fig. 3d–e).

## The evidence

- Against its plate × line vehicle at θ* = 0.1 rad, 55.4 % of conditions are control-pool-limited, 5.6 % not
  detectable and only 11.4 % over-sampled; the median vehicle floor m_min = 1.22 exceeds the median bias-corrected
  magnitude 1.096 ([[huuloc-2026-cell-quota]], p. 7; Table 1, p. 8). Within-condition per-cell variance in a
  50-PC embedding: σ² = 0.957 (p. 8).
- Magnitude is set jointly by compound and line (homoharringtonine 5 µM: m = 14.9 in NCI-H460, 3.0 in NCI-H661;
  pp. 8–9).
- Depth: median 2,298 UMIs per cell per [[huang-2025-x-atlas-orion]] (p. 5; Fig. 2E), against 19,416 for
  [[x-atlas-orion]] HCT116 (8.45× deeper).

- Primary source, ingested from the preprint [[zhang-2025-tahoe]] (bioRxiv v3, 2025-05-10; preprint of
  @zhang-2026-tahoe): 50 cancer lines pooled in every well ("Mosaic" design; p. 3, Fig. 1A), 379 drugs × 3 doses
  (0.05 / 0.5 / 5 µM; p. 5, Fig. 4B), 24 h, 14 plates, 100,648,790 cells before and 95,624,334 after QC (pp. 3, 21);
  lines assigned by SNP demultiplexing (demuxlet, > 98 % accuracy, > 90 % singlets; p. 19). Only **47** lines are
  analysed: NCI-H661, NCI-H596 and NCI-H2122 are excluded for low counts (p. 19); 52,886 conditions, median 1,287
  cells, 1,135 drug-dose combinations after QC (p. 7). Mean 2,288 and median 1,890 transcripts per cell (p. 3).
  Plate 14 replicates plate 6: matched-condition median r 0.975, unmatched 0.915 (p. 7; Fig. 3D). Two DMSO wells per
  plate (H11/H12; p. 18). No modelling benchmark and no noise or variance decomposition (scVI only; pp. 20–21).
- [[adduri-2025-state]] hold out 5 Tahoe lines (C32, HOP62, HepG2/C3A, Hs 766T, PANC-1; Fig. S3, p. 62).

## Contested

- Reference choice changes the estimand: a per-line-mean reference removes the floor (m_min ≈ 0.05) but measures
  "difference from the average drug" ([[huuloc-2026-cell-quota]], p. 9).
- Molina & Zhang's Tahoe result is **not** zero-shot in the context sense: the target compound's response in source
  lines is an input (p. 17).

## Relevance to this project

This is the H2a/H2b testbed (learned context correction vs context-free, C ∈ {6 … 40}). The source does not address
it. The source's 5-line use offers no context-count curve. The project's H2 response is matched to the plate × line
DMSO pool, the reference that [[huuloc-2026-cell-quota]] show to be noise-limiting (p. 7); this is the mechanism
behind the T8 shared-control-noise finding in [[rq2_params-assumptions]].

## Open / to verify

- The exact number of lines in the project's Tahoe build (48 in the plan). The atlas has 50 lines, of which the
  preprint analyses 47 ([[zhang-2025-tahoe]], p. 19). The project's 48 contain NCI-H661 but neither NCI-H596 nor
  NCI-H2122 (checked 2026-11-06 in `vcc2026_data/data/processed/drug_context_pretrain/contexts.npz`: 48 `tahoe:`
  names), i.e. the preprint's 47 + NCI-H661, assuming HEPG2_FAMILY = HepG2/C3A (not stated in any source).
- Counts that the preprint does not confirm: 56,827 conditions / median 1,296 cells (preprint: 52,886 / 1,287 on 47
  lines, p. 7), "2–3" DMSO wells per plate (preprint: two, p. 18), median 2,298 UMIs per cell (preprint: mean 2,288,
  median 1,890, p. 3), 86.3 % of combinations on one plate (not stated). The 2.87 noise/interaction ratio is the
  project's own estimate (the preprint has no decomposition).

## Sources

- [[molina-2026-response-decomposition]] — pp. 3–4, 10, 17–18.
- [[huuloc-2026-cell-quota]] — pp. 3–4, 7–10; Tables 1–2.
- [[zhang-2025-tahoe]] — pp. 1, 3, 5, 7, 18–21; Fig. 3D.
- [[adduri-2025-state]] — Fig. S3, p. 62.
- [[huang-2025-x-atlas-orion]] — depth comparator only; p. 5.
