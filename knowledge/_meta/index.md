---
title: "Wiki Index"
type: synthesis
created: 2026-11-05
updated: 2026-11-06
status: review
author: llm
---

# Wiki index

As of: 2026-11-06

## Entities

- [[replogle-nadig-crispri-collection]] — K562/RPE1/HepG2/Jurkat CRISPRi Perturb-seq (primary: Replogle and Nadig preprints)
- [[tahoe-100m]] — chemical perturbation atlas; H2a/H2b testbed (primary: Tahoe preprint)
- [[x-atlas-orion]] — genome-wide CRISPRi Perturb-seq, HCT116 + HEK293T (Fix-Cryopreserve); OR_HCT116X source
- [[h1-hesc-crispri-perturb-seq]] — H1 hESC CRISPRi Perturb-seq (GSE295214; registry family H1GW)
- [[zhu-cd4-t-cell-perturb-seq]] — primary CD4+ T-cell genome-scale CRISPRi, 4 donors × 3 timepoints
- [[prophet-universal-predictor]] — Prophet (Ji et al.): transformer over cell line × perturbation × readout

## Concepts

- [[perturbation-response-decomposition]] — μ + α_c + β_p + γ with split-half noise correction (plan notation)
- [[zero-shot-transfer-ceiling]] — r_max for H1 (project-defined; implied ≈0.84 from Molina & Zhang)
- [[control-subspace-projection]] — response energy on control PCs; H3a response-level precedent
- [[context-count-learning-curve]] — accuracy vs number of training contexts; crossing of a context-free baseline (H2a)
- [[simple-perturbation-baselines]] — no-change, additive, mean, linear baselines (H4 context)
- [[systematic-variation]] — Systema: shared perturbed-vs-control shift; centred metrics (H3b alternative)
- [[directional-cell-quota]] — cells needed to resolve a response direction; control-pool term (H1/H4)

## Sources

- [[molina-2026-response-decomposition]] — Molina & Zhang 2026, Perturbation response decomposition (deep)
- [[huuloc-2026-cell-quota]] — Huu Loc & Van 2026, closed-form control-aware cell quota (standard)
- [[dibaeinia-2026-context-scale]] — Dibaeinia et al. 2026, Virtual cells need context, not just scale (standard)
- [[ahlmanneltze-2025-linear-baselines]] — Ahlmann-Eltze et al. 2025, linear baselines (Nat. Methods; standard)
- [[piao-2026-adaptive-contexts]] — Piao et al. 2026, AdaPert (ICML; standard)
- [[ji-2024-universal-prediction]] — Ji et al. 2024/25, Prophet (bioRxiv v3; deep) — prior context-count evidence
- [[rosenfeld-2021-risks-irm]] — Rosenfeld et al. 2021, The risks of IRM (ICLR; E > d_e theorem; standard)
- [[vinastorne-2026-systema]] — Viñas Torné et al. 2026, Systema (standard; data/benchmark focus)
- [[huang-2025-x-atlas-orion]] — Huang et al. 2025, X-Atlas/Orion (standard; data focus)
- [[pan-2026-global-cell-state]] — Pan et al. 2026, global cell-state manifold + H1 hESC screen (standard; data focus)
- [[adduri-2025-state]] — Adduri et al. 2025, State (bioRxiv v2; preprint of adduri-2026-state): few-shot (30 %) context benchmark, mixed genetic-data gain; zero-shot via Tahoe pre-training (standard)
- [[li-2024-systematic-comparison]] — Li et al. 2024, systematic comparison of perturbation models (bioRxiv v3; preprint of li-2026-systematic-comparison): baselines win on seen combinations; diagnostic non-orthogonal decomposition (standard)
- [[replogle-2021-genome-scale]] — Replogle et al. 2021, genome-scale Perturb-seq, K562 GW/essential + RPE1 essential (bioRxiv v3; preprint of replogle-2022-genome-scale; standard)
- [[nadig-2024-transcriptome-wide]] — Nadig et al. 2024, TRADE transcriptome-wide impact; HepG2/Jurkat essential screens (bioRxiv v1; preprint of nadig-2025-transcriptome-wide; standard)
- [[zhang-2025-tahoe]] — Zhang et al. 2025, Tahoe-100M atlas: 50 lines (47 analysed), 379 drugs × 3 doses, Mosaic/SNP demux (bioRxiv v3; preprint of zhang-2026-tahoe; standard)

## Synthesis

- [[results-summary]] — All pre-registered verdicts with deciding numbers, exploratory results, deviations
- [[related-work-positioning]] — What prior work shows per hypothesis and what this paper adds
- [[theory-context-complexity]] — Theory: identifiability (Prop. 1, cf.) and the rejected OLS break-even (Prop. 2); consequences for H2a–c
- [[rq2_tahoe_curve-assumptions]] — Task 6 Tahoe runner: protocol, descriptor, review notes
- [[rq2_h2a-assumptions]] — Task 7 Tahoe context curve + H2a verdict: checks, caveats (early stopping, ablate confound)
- [[rq2_params-assumptions]] — Task 8 theory-parameter estimators (d, σ²_ε, Var(i), R²)
- [[task04-freeze-inputs-assumptions]] — Task 4 input freeze and source registry
- [[e4_official-assumptions]] — Task 19 official VCC record: trajectory, 129/1,635, single-change gains (M2 +0.0106), M3 bundle not attributable, shrinkage
- [[appx-drug-gene-bound]] — Task 20 appendix: Cauchy–Schwarz bound on target-identity drug→CRISPRi transfer (ρ* ≈ 0.06; Jiang = CRISPRi)
- [[rq3_matchedn-exploratory]] — EXPLORATORY matched-n sensitivity for H3a/H3b: H1 D_other halves (+0.009), RN_K562/RN_Jurkat flip positive; H3b M3 − M3_other unchanged (+0.0058)
- [[rq1_components-assumptions]] — Task 12 RQ1 noise-corrected variance components and ceilings (7 harness panels)
- [[rq1_h1-assumptions]] — Task 13 RQ1 methods vs transfer ceilings and the H1 verdict
