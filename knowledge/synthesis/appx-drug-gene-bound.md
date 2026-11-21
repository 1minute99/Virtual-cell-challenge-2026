---
title: "Appendix — upper bound on target-identity-mediated drug → genetic-perturbation transfer"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# Appendix: drug → gene bound (Task 20)

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 20.

**Data source.** The frozen exploratory analysis in `vcc2026_data/reports/drug_gene_concordance_2026-09-16/`:
- inputs: `RESULTS.md` (Korean), `ANALYSIS_PLAN.md`, `pair_metrics.csv`, `repeat_metrics.csv`, `summary.csv`
- scope: A549 and HT-29; Tahoe inhibitors against two genetic sources:
  - Jiang CRISPRi, four stimuli: 110 conditions
  - LINCS shRNA: 291 pairs

**Code / outputs**
- `output/data-analysis/appx_drug_gene.py` (seed 20261105, B = 2,000)
- `output/data-analysis/results/appx_drug_gene.json` and `appx_drug_gene.csv`
- `output/data-analysis/results/appx_drug_gene.png`
- log: `results/logs/appx_drug_gene.log`

**Explicit exclusion.** The earlier gene-side drug-pretraining experiment is **invalid and is neither cited nor used**:
- experiment: `vcc2026_data/reports/drug_context_pretrain_2026-09-20/`
- why it is invalid: its gene context vectors were scrambled. Sorted gene names were placed over source-order counts
  for 5 of 7 biologies. See AGENTS.md, "BUG FOUND 2026-09-21", and that experiment's FINDINGS.md caveat ("CAVEAT
  added 2026-09-21").
- the script reads nothing from that directory, and no conclusion here depends on it.

## Modality of the "Jiang" screen: CRISPRi
Jiang et al. 2025 [@jiang-2025-pathway-signatures], Nat Cell Biol 27:505–517, doi:10.1038/s41556-025-01622-z. The
data come from GEO GSE281048 and Zenodo 14518762.
- The GEO series description states: "we modified each of these lines to express a CRISPR interference (CRISPRi)
  dCas9-KRAB-MeCP2 cassette". The guides are "from the Dolcetto genome-wide CRISPRi library".
- The paper's Fig. 2 is titled "Mixscale for weighting and analysing CRISPRi perturbation effects".

So this is knock-down, not knockout. LINCS is shRNA knock-down (96 h); the Tahoe drug treatments are 24 h.

## Derivation (Cauchy–Schwarz on the noise-corrected concordance)
Take one triple: drug d, annotated target t, line l. The genes are the samples, in rank space.

**Model.** Let the drug response be D = s_D + e_D and the genetic response be G = s_G + e_G. The noise terms e are
independent of each other and of the signals s.

**Decomposition.** Write s_G = β s_D + u with ⟨u, s_D⟩ = 0.

**Transfer predictor.** A target-identity-mediated transfer predictor ("use the response of a drug that hits t as
the prediction for knocking down t") has the form ĝ = a s_D + w, where w carries no target-specific information:
⟨w, s_G⟩ = ⟨w, s_D⟩ = 0 once the generic response shared with wrong targets is removed. Then

  corr(ĝ, s_G) = a⟨s_D, s_G⟩ / (‖ĝ‖ ‖s_G‖) ≤ |⟨s_D, s_G⟩| / (‖s_D‖ ‖s_G‖) = |ρ\*|.

The inequality holds because ‖ĝ‖ ≥ |a| ‖s_D‖ (w ⊥ s_D), which is Cauchy–Schwarz. Equality holds iff w = 0.

**Estimating ρ\*.** The target-identity part of the concordance is Δ = ρ(D_t, G_t) − mean over wrong targets t′ of
ρ(D_t, G_t′). The wrong targets exclude every annotated target of d. Spearman's disattenuation then gives

  **ρ\* = Δ / √(R_D R_G)**,

where:
- R_G = 2r_G/(1+r_G), the Spearman–Brown correction for the Rep1/Rep2 mean of Jiang. r_G is the Rep1 vs Rep2
  correlation of the same knock-down.
- R_D = r_D, the Tahoe primary-vs-repeat correlation. Only the primary treatment is used.

**Derived quantities.**
- **ρ\*²** is the largest share of the *reproducible* genetic-response variance that this transfer can explain.
- **ρ\*·√R_G** is the largest correlation it can reach against the *observed* genetic truth.

**Check of the attenuation model.** Averaging two genetic replicates should multiply Δ by √(2/(1+r̄_G)) = 1.353.
The observed ratio is 1.354, pooled over the four Jiang stimuli (Δ_avg / Δ_single-rep). The same ratio for the
matched ρ is 1.27.

## Results (raw effects, primary; equal weight per drug; 95 % percentile bootstrap)
| comparison | drugs | Δ_obs | R_D | R_G | **ρ\*** [95 % CI] | ρ\*² | bound vs observed truth |
|---|---|---|---|---|---|---|---|
| Jiang IFNB | 7 | +0.0146 | 0.204 | 0.174 | 0.078 [−0.214, 0.218] | 0.006 | 0.032 |
| Jiang INS | 10 | −0.0048 | 0.204 | 0.172 | −0.026 [−0.055, 0.002] | — | — |
| Jiang TGFB | 16 | +0.0046 | 0.204 | 0.129 | 0.028 [−0.021, 0.088] | 0.001 | 0.010 |
| Jiang TNFA | 13 | +0.0313 | 0.204 | 0.206 | **0.153 [0.036, 0.253]** | 0.023 | 0.069 |
| **Jiang, 4 stimuli (equal)** | — | — | — | — | **0.058 [−0.019, 0.103]** | 0.003 | — |
| LINCS shRNA | 66 | +0.0108 | 0.292 | n/a | drug-side only: 0.020 [0.006, 0.036] | — | — |

**Headline.** For a target-identity-mediated drug → CRISPRi transfer predictor, the noise-corrected correlation
ceiling is **ρ\* ≈ 0.06** (pooled over the four stimuli), so at most **0.3 %** of the reproducible knock-down
response is transferable. The most permissive single condition is TNFA, with ρ\* = 0.15 and an upper 95 % limit of
0.25. Even that limit explains at most **6.4 %** of the reproducible variance.

**LINCS.** There is no shRNA replicate, so R_G is unidentified. The drug-side-only value 0.020 is a floor of the
bound. If R_G were 0.3 the bound would be 0.037; if R_G were 0.1 it would be 0.063.

**Sensitivity (generic response removed; raw reliabilities reused).**

| comparison | ρ\* [95 % CI] |
|---|---|
| Jiang, pooled | 0.099 [0.026, 0.150] |
| TNFA | 0.215 [0.081, 0.337] |
| LINCS, drug-side only | 0.019 [0.004, 0.036] |

## Assumption checks
| check | expected | got | pass |
|---|---|---|---|
| Reproduce `summary.csv` Δ and ρ with the D5 aggregation | 0 | 1.5e-16 | ✓ |
| Genetic replicate ρ available for every Jiang pair | yes | yes | ✓ |
| Tahoe repeat drugs | — | 7 drugs on the Jiang axis, 26 on the LINCS axis; r_D spread across stimuli for the same drug×line ≤ 0.0015 | ✓ (small n) |
| LINCS genetic replicate | — | absent | bound unidentified (D3) |
| Attenuation-model ratio | 1.353 | 1.354 | ✓ (noisy) |
| Bootstrap draws with R_D·R_G ≤ 0 | 0 | 0 | ✓ |

## Evidence base
- Source analysis: `vcc2026_data/reports/drug_gene_concordance_2026-09-16/`. It had a pre-written analysis plan
  and was independently verified (`verification.json`, `independent_verification.json`).
- Its own plan declined a reliability correction because shared controls and study differences complicate the
  assumptions. **This appendix applies the correction anyway**, as the Task 20 spec requires, and states the
  direction of the remaining bias below.
- Modality: [@jiang-2025-pathway-signatures], GEO GSE281048.
- Related project pages: [[perturbation-response-decomposition]], [[zero-shot-transfer-ceiling]],
  [[tahoe-100m]].

## What stays open
- **Bias direction.** If replicates share non-targeting controls, r_G and r_D are inflated. That makes ρ\* too small,
  so the bound is anti-conservative. The effect is not identifiable from these files.
- **Small samples.** R_D rests on 7 Tahoe drugs on the Jiang axis. IFNB has only 4 candidate targets and 7 drugs,
  hence its very wide CI.
- **Approximation.** Rank-space disattenuation approximates the Pearson result. The Spearman–Brown step assumes the
  two Jiang replicates are exchangeable.
- **Scope of the bound.** It covers predictors whose target-specific content is the matched drug signature. A
  learned gene × gene map trained across many targets is not covered.
- **Comparability gaps.** Same-named lines are not the same cultures. Time points differ: 24 h drug vs 24 h
  stimulus + CRISPRi vs 96 h shRNA. Drug target annotations are the publisher's and have not been validated
  pharmacologically.
