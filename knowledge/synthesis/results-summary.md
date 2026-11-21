---
title: "Results summary: pre-registered verdicts, effect sizes, deviations"
type: synthesis
created: 2026-11-06
updated: 2026-11-06
status: review
author: mixed
---

# Results summary

Plan (pre-registered 2026-11-05, Amendments A1–A3): `input/ideas/zero-shot-limits-context-complexity-plan.md`.
Every number traces to a file in `output/data-analysis/results/`; the full chronology (verdicts, deviations,
clarifications, reviews) is in `knowledge/_meta/log.md`. Theory: [[theory-context-complexity]]. Positioning:
[[related-work-positioning]].

## Verdicts

| Hypothesis | Verdict | Deciding numbers | File |
|---|---|---|---|
| H2a (Tahoe break-even, A1 context-pathway gain) | **confirmed** | gain at C = 12: −4.20e-4 [−5.15, −3.32]e-4; C*_emp = 37.4 [33.3, > 40] | `rq2_h2a.json` |
| H2b (Prop. 2 predicts C*_emp within ×2) | **rejected** | C*_theory 230 [192, 287] / 37.4 = 6.2; C*_min 80.6 > C*_emp | `rq2_h2b.json` |
| H2c (genetic C*_min > C_public, both counts) | **holds for distinct cell types, not for donor lines** | C*_min 16.45 [16.09, 16.82] vs Count 1 = 11, Count 2 = 25 | `rq2_h2c.json` |
| H1 (pipeline ≥ 0.7 of r_max on the H1 panel) | **untestable (as registered)** | oracle r 0.490 [0.465, 0.521] > registered r_max 0.423 → ceiling estimator inconsistent | `rq1_h1.json` |
| H3a (transfer residual on destination control PCs, H1, k = 8) | **confirmed** | D_other 0.0182 [0.0092, 0.0267]; D_perm 0.0782 [0.0653, 0.0922]; neg. control ≤ 2 % | `rq3_h3a.json` |
| H3b (M3 beats five alternatives) | **inconclusive** | 7-panel ΔPDS > 0 vs all five, CIs exclude 0; H1 real scorer: Systema-style +0.0076 over M3 | `rq3_h3b.json` |
| H4 (H1GW beats the largest modelling gain) | **confirmed** | +0.0690 composite (w = 1) vs M2 +0.0065 (≈ ×10.6); robust to the seed-mean architecture gain +0.0075 (≈ ×9.2) | `rq4_h4.json`, `rq2_arch/rq2_arch_seeds_summary.json` |

## Exploratory and secondary results (labelled as such in the manuscript)

- **Ceiling (H1 panel):** against the bounded ceiling r_proj = 0.492 the pipeline attains 0.773 [0.722, 0.810];
  transferable share of the remaining error < 0.5 (bias bound 0.39–0.46). Other panels with the registered r_max:
  KOLF, RN_Jurkat, RN_RPE1 "model-limited"; RN_K562, RN_HepG2 inconclusive; OR not estimable. Measurement ceilings
  0.34 (OR) – 0.88 (H1).
- **Matched-n RQ3:** H1 D_other halves to 0.0093 [0.0003, 0.0176] (CI includes 0 at common n = 800 / 267); D_other > 0
  on 7/7 panels once PC-estimation noise is matched (registered RN_K562 / RN_Jurkat negatives were an artefact);
  M3 vs other-line PCs unchanged (+0.0058, 7/7).
- **Theory sensitivities:** target b + i gives C*_theory 192, C*_min 65 — still far above 37.4.
- **Source-count curves (transfer pipeline, S4 held fixed):** slopes +0.004…+0.017 PDS per ln k; leave-one-family-out
  values small, CIs include 0 everywhere; H1 per-context composite value +0.00086 [−0.00045, +0.00255].
- **Official record:** 38 submissions; end-of-validation rank 129 / 1,635 (top 7.9 %), best 34 / 1,244 (top 2.7 %);
  only M2 is a single-change official gain (+0.0106); local → official shrinkage ratio 0.053 (S4, the only valid case).
- **Drug → gene bound (appendix):** reliability-corrected ρ* = 0.058 [−0.019, 0.103] (≈ 0.3 % of reproducible CRISPRi
  variance); TNFA 0.153.
- **Learned-architecture seed replication (Task 10):** seeds 0–2, six architectures, K562 + HepG2 (seed-0 gate passed,
  |Δ| 1.7e-16). No learned architecture beats the refit/fusion baseline on the two-line mean with a CI excluding 0;
  additive beats it on K562 (+0.0240 [0.0048, 0.0433]) and loses on HepG2 (−0.0090 [−0.0132, −0.0047]). Updated best
  two-line gain +0.0075 [−0.0010, +0.0161] (additive); H4 stays confirmed (H1GW ≈ 9.2× it; ≈ 2.9× the single-line
  +0.0240). LoRA hypernetwork stays seed-0 only.

## Deviations that the manuscript must state

A1 (H2a/H2b quantity: context-pathway gain against a matched constant-descriptor ablation; cluster bootstrap; censored
H2b = untestable) — approved before Task 7 data, one disclosed smoke run; A2 (C_public two counts); A3 (truth paths,
KOLF leak fixed by rebuild, 38 submissions, framing). Task 6 reproduction criterion replaced (initialisation noise).
Tahoe cell-level shards and other inputs deleted by an earlier cleanup were restored from public sources (sha256).
The analysing agent's "confirmed" reading of H1 was not adopted.

## Evidence base

Result files listed above; per-task assumption notes [[rq2_h2a-assumptions]], [[rq2_params-assumptions]],
[[rq2_gene_params-assumptions]], [[rq1_components-assumptions]], [[rq1_h1-assumptions]], [[rq3_h3a-assumptions]],
[[rq3_h3b-assumptions]], [[rq3_matchedn-exploratory]], [[rq4_h4-assumptions]], [[rq2_source_curve-assumptions]],
[[e4_official-assumptions]], [[appx-drug-gene-bound]]. Independent quality reviews ran for Tasks 6, 8, 15–16.

## What stays open

- The empirical break-even rests on one learner family and one drug dataset; no public genetic dataset has enough
  contexts to measure the genetic break-even directly.
- H1 and H3b are not decided by the registered rules; the manuscript reports them as untestable / inconclusive.
