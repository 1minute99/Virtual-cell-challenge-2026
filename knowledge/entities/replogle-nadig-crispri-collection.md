---
title: "Replogle–Nadig four-line CRISPRi Perturb-seq collection (K562, RPE1, HepG2, Jurkat)"
type: entity
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
tags: [dataset, crispri, perturb-seq]
sources: ["[[pan-2026-global-cell-state]]", "[[molina-2026-response-decomposition]]", "[[huuloc-2026-cell-quota]]", "[[ahlmanneltze-2025-linear-baselines]]", "[[piao-2026-adaptive-contexts]]", "[[vinastorne-2026-systema]]", "[[huang-2025-x-atlas-orion]]"]
---

# Replogle–Nadig four-line CRISPRi collection

## Location and authority data

This is a dataset, so no ORCID, GND or Wikidata identifier applies. Primary publications: @replogle-2022-genome-scale
(K562, RPE1) and @nadig-2025-transcriptome-wide (HepG2, Jurkat). Both are ingested from their open preprints
(2026-11-06): [[nadig-2024-transcriptome-wide]] (bioRxiv v1) and [[replogle-2021-genome-scale]] (bioRxiv v3).

## What one sees

These are four essential-scale CRISPRi Perturb-seq screens. As processed by [[molina-2026-response-decomposition]]
(p. 10), after excluding perturbations with < 30 cells:

| Line | Cells | NT controls | Perturbed genes |
|---|---|---|---|
| K562 | 188,590 | 10,691 | 1,383 |
| RPE1 | 173,737 | 11,485 | 1,499 |
| HepG2 | 96,616 | 4,976 | 1,340 |
| Jurkat | 184,470 | 12,013 | 1,537 |

The four lines share 6,642 genes and 638 perturbations (pp. 11, 20). RPE1 is the only non-transformed line (p. 2).

## History of research

The K562 and RPE1 screens come from @replogle-2022-genome-scale; HepG2 and Jurkat come as processed matrices from
@nadig-2025-transcriptome-wide (p. 10). [[molina-2026-response-decomposition]] use the four lines as the testbed for
the [[perturbation-response-decomposition]] and for 1- vs 3-source cross-line transfer.

## The evidence

- Noise-corrected shares: template 27.8 %, conserved 29.4 %, interaction 23.5 %, noise 19.3 % (p. 2).
- Template share of response energy: 23 % (K562), 65 % (RPE1), 49 % (HepG2), 28 % (Jurkat) (Table 2, p. 35).
- [[huuloc-2026-cell-quota]] reprocess the HepG2 and Jurkat screens ("TRADE essential-gene screens", GSE264667;
  2,393 genes, median 48–85 cells per gene, NTC 4,380 / 11,514 cells; pp. 3–4, Table 1 p. 8). Knockdowns are strong
  (median m 2.25 / 1.84) and 60.9 % / 70.6 % detectable, but only 4.4 % / 1.2 % of detectable knockdowns are deep
  enough to resolve their direction at 0.1 rad: the screens are depth-limited (p. 11; Table 4).

- [[ahlmanneltze-2025-linear-baselines]] use the Replogle K562 and RPE1 essential screens **as provided by GEARS**
  (5,000 genes; 162,264 / 161,423 cells; 1,087 / 1,534 single perturbations; ED Fig. 1, PDF p. 8) for an unseen-
  perturbation benchmark. A linear model whose perturbation embedding is pretrained on the *other* line is the only
  consistent winner over the mean baseline (p. 1660; Fig. 2c), and the two lines' baseline expression correlates 0.81
  (ED Fig. 9a, PDF p. 17).

- [[piao-2026-adaptive-contexts]] (AdaPert, ICML 2026) use all four lines with perturbation-level splits (K562 5,000
  HVGs, 1,088 entries; the others 3,352 HVGs; Table 4, p. 13). In the one cross-line test (RPE1 + HepG2 + Jurkat →
  K562), the context-free per-perturbation mean of the source lines beats every model (Pearson-∆ 0.650 vs ≤ 0.515;
  Table 13, p. 22). The paper calls all four screens «genome-wide» (p. 12) and does not cite Nadig et al.

- [[vinastorne-2026-systema]] use Replogle **genome-wide K562** and **RPE1** in a ten-dataset unseen-perturbation
  benchmark (scGPT processing: the 1,973 strong-phenotype perturbations, 100 cells per perturbation, 2,500 controls;
  1,813 / 1,410 perturbations after processing; PDF p. 12; Fig. 1b, p. 1052). RPE1 carries more
  [[systematic-variation]] than K562 (Fig. 3b, p. 1054): 46 % vs 25 % of perturbed vs control cells are in G1, which the
  authors attribute to CIN-induced arrest in p53-positive RPE1; K562 is p53-negative and shifts much less (p. 1053).

- [[huang-2025-x-atlas-orion]] use Replogle **genome-wide K562** as the benchmark for the [[x-atlas-orion]] atlas: median
  11,499 UMIs per cell (Orion HCT116 is 1.68× deeper), NTC batch-to-batch Spearman median 0.967 (vs 0.993 / 0.988 for
  Orion HCT116 / HEK293T; MAD 9.09e-3 vs 2.87e-3 / 3.58e-3), and 54.30 % of expressed-gene perturbations discernible
  by an L1 classifier (vs 48.36 % / 40.37 %) (pp. 5–8; Fig. 2F, PDF p. 7; Fig. 3A).

- [[pan-2026-global-cell-state]] use Replogle **K562 and RPE1** as the only genome-scale contexts of their
  mixture-of-experts training set (3 of 7 Perturb-seq datasets; p. 10) and hold each out in turn for few-shot
  evaluation (p. 11). Two reproducibility facts from there: only **608 perturbed genes** are concordant at Pearson
  r > 0.4 between the **two genome-scale K562 screens (day 6 vs day 8)** (pp. 10, 58), and K562 perturbation effect
  sizes are significantly larger than in H1 hESCs at matched knockdown (Fig. 6D, p. 40). Pathway-level rewiring
  between the two lines is explicit: mitochondrial transcription/translation knockdown activates ISR in K562 but not
  RPE1, and p53 responds broadly in RPE1 but is absent in p53-deficient K562 (p. 13).

- [[replogle-2021-genome-scale]] (preprint of @replogle-2022-genome-scale) is the primary source of the K562 and RPE1
  screens: K562 genome-wide (all expressed genes, day 8; 9,866 genes / 10,673 perturbations), K562 essential (day 6;
  2,057 / 2,176) and RPE1 essential (day 7; 2,393 / 2,549) (pp. 5–6; Fig. 1A, p. 32 gives slightly different n);
  CRISPRi with two sgRNAs per construct, dCas9-KRAB in K562 and ZIM3-KRAB in RPE1 (pp. 5–6); > 2.5 M cells, median
  > 100 cells per perturbation (p. 6), median ~10,000 UMIs per cell (Fig. 1A, p. 32); median knockdown 85.5 % (K562)
  and 91.6 % (RPE1) (p. 6). 31.1 % of K562 genome-wide perturbations (2,987 / 9,608) pass the energy-distance test
  against 1.9 % of non-targeting controls, explicitly a lower bound (p. 7). K562 vs RPE1 is compared qualitatively only
  (e.g. CIN-driven G1/G0 arrest in p53-positive RPE1, little shift in K562; p. 15). Day 8 = genome-wide and day 6 =
  essential-only, so "two genome-scale K562 screens (day 6 vs day 8)" is imprecise.
- [[nadig-2024-transcriptome-wide]] (preprint of @nadig-2025-transcriptome-wide) is the primary source of the
  **HepG2- and Jurkat-Essential** screens: new dual-sgRNA CRISPRi screens of 2,393 essential genes, 10x 3′ v3, 56 GEM
  groups; 262,956 Jurkat and 145,473 HepG2 cells retained after QC (p. 25); Zim3-dCas9 in Jurkat, KOX1 dCas9-KRAB in
  HepG2 (pp. 23–24). Median cells per perturbation (Table 1, p. 9): HepG2 45, RPE1 72, Jurkat 83, K562-Ess 121,
  K562-GW 178. Only K562-GW is genome-wide; the other four screens are essential-scale. Noise-corrected (TRADE)
  cross-line correlations of perturbation effects: medians 0.40 (K562/RPE1) to 0.75 (HepG2/RPE1) over 1,660
  perturbations; 56 % of perturbations correlate highly across all lines, 44 % only within groups of similar lines
  (pp. 13–14; Fig. 4C–D). Raw vs noise-corrected replicate correlation between two K562 screens: 0.16 vs 0.90 (p. 13).
- [[adduri-2025-state]] (preprint of @adduri-2026-state) uses the four lines as "Replogle–Nadig" (4 lines, 1,677
  perturbations in Table 1, p. 31; the text says 2,024, p. 7). Its leave-one-line-out benchmark is **few-shot**: 30 %
  of the held-out line's perturbations are added to training (p. 9).

## Contested

- With only four lines, any context-count conclusion from this collection is limited to 1 vs 3 sources (p. 15).
- Resolved (2026-11-06): STATE's leave-one-line-out on these lines is few-shot, 30 % of the target line's
  perturbations in training ([[adduri-2025-state]], p. 9).
- Cell counts differ by processing: Molina & Zhang's post-filter HepG2 96,616 / Jurkat 184,470 vs the preprint's
  QC-retained 145,473 / 262,956 ([[nadig-2024-transcriptome-wide]], p. 25); median cells per gene 48–85 in
  [[huuloc-2026-cell-quota]] vs 45 (HepG2) / 83 (Jurkat) in the preprint's Table 1 (p. 9).

## Relevance to this project

These lines feed the harness's RN panels (e.g. RN_K562, whose truth was built from K562 GWPS per plan Amendment A3)
and count towards C_public for H2c.

## Multiple roles

The collection is a data source for the M1 registry, a destination for harness panels, and the substrate of the
source paper's variance shares.

## Open / to verify

- Whether the project's gene space and filtering reproduce Molina & Zhang's shares.

## Sources

- [[molina-2026-response-decomposition]] — pp. 2, 10–11, 20, 35.
- [[replogle-2021-genome-scale]] — pp. 5–7, 15; Fig. 1A, p. 32.
- [[nadig-2024-transcriptome-wide]] — pp. 7, 9, 13–14, 23–25.
- [[adduri-2025-state]] — pp. 7, 9, 31.
- [[huuloc-2026-cell-quota]] — HepG2 / Jurkat only; pp. 3–4, 8, 11.
- [[ahlmanneltze-2025-linear-baselines]] — K562 / RPE1 only (GEARS processing); p. 1660; PDF pp. 6, 8, 17.
- [[piao-2026-adaptive-contexts]] — all four lines; pp. 12–13, 21–22.
- [[vinastorne-2026-systema]] — K562 GW / RPE1 only; pp. 1052–1055; PDF p. 12.
- [[huang-2025-x-atlas-orion]] — K562 GW as comparator only; PDF pp. 5–8.
- [[pan-2026-global-cell-state]] — K562 and RPE1 as MoE training/held-out contexts; pp. 10–13, 58; Fig. 6D p. 40.
