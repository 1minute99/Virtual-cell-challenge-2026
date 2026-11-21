---
title: "Related-work positioning: what is already shown, and what this paper adds"
type: synthesis
created: 2026-11-05
updated: 2026-11-06
status: stable
author: mixed
---

# Related-work positioning

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`. Theory: [[theory-context-complexity]]. Per-source
overlap tables are on each source page (rows H1–H4). Sources marked † are cited from their bibliographic record and
the project's logged facts only; their full texts are paywalled and not yet ingested.

## 1. Learned models do not beat simple baselines (background, not our claim)

- [[ahlmanneltze-2025-linear-baselines]]: deep models do not beat mean, linear or additive baselines; a
  perturbation-agnostic mean already reaches high Pearson-Δ. Two lines, no unseen-line task. Our addition: the
  context-count axis and a matched ablation that isolates the context pathway.
- [[piao-2026-adaptive-contexts]] (AdaPert): context-conditioned models lose to a context-free cross-line mean when
  trained on C = 3 lines (Table 13). Same direction as our negative context gain at small C, one point only.
- [[vinastorne-2026-systema]]: control-referenced metrics are inflated by the shared perturbed-vs-control shift;
  average-effect baselines rival deep models. Source of our "Systema-style" H3b arm, which is an *analogue* (we subtract
  the mean predicted response from predictions; Systema changes the reference inside the metric).
- The VCC 2026 description [@adduri-2026-virtual-cell-challenge]† defines the zero-shot unseen-context benchmark.
  STATE ([[adduri-2025-state]], preprint of @adduri-2026-state) is benchmarked mainly on an *underrepresented-context*
  (few-shot) task: on the genetic Replogle–Nadig data each line is held out with 30 % of its perturbations kept in
  training (p. 9). There its gain is mixed, not uniformly small: parity with the perturbation mean on discrimination
  (−3 %; «State matched the performance of the perturbation mean baseline», p. 10), +5 % Pearson-Δ with the context
  mean the strongest baseline, losses on DE overlap (−21 %) and fold change (−12 %), gains on DE AUPRC (+184 %) and
  effect size (+22 %) (Fig. 2D–J, p. 8), against «notably larger» gains on Tahoe and Parse (p. 11). Its zero-shot
  unseen-line task reports only effect-size ranking on Replogle–Nadig (≈ 0.57 vs ≈ 0.47 for the perturbation mean,
  +23 %; Fig. 3G, p. 12; read from plot).

## 2. Context diversity, not scale (closest framing)

- [[dibaeinia-2026-context-scale]] (ICML 2026 position): the failure mode of virtual cells is insufficient coverage of
  biological contexts, not model expressivity; no baseline crossing and no context count (Fig. 4). Our addition: a
  number — the matched-ablation break-even ≈ 37 contexts on 48 drug-treated lines — and the 11 distinct public genetic
  cell types to compare it with.
- [[ji-2024-universal-prediction]] (Prophet): the prior empirical evidence for H2a. Learning curves over the share of
  seen cell lines against a context-free mean, crossing at dataset-dependent points (≈ 6 lines on LINCS, ≈ 170 on GDSC,
  never on CTRP / GDSCcomb / SCORE; v3 Fig. S6A). Differences: scalar viability/IC50 readouts, not transcriptomes; the
  baseline is the mean vs a separately trained model, not a matched ablation; no theory. The spread across datasets
  itself argues against a single closed-form C*, consistent with our rejected H2b.
- [[rosenfeld-2021-risks-irm]]: environment-count identifiability (IRM optimal only with more environments than the
  environmental-feature dimension; linear, infinite-sample, binary, noise-free). Our Proposition 1 is cited as its
  context-count analogue ("cf."); our noise-dependent break-even (Proposition 2) was tested and rejected quantitatively.

## 3. Response decomposition and ceilings

- [[molina-2026-response-decomposition]]: noise-corrected four-component decomposition (shared, conserved target,
  context, interaction); the interaction is never recovered zero-shot; the observed template and template-removed
  residual are projected onto the same line's control PCs with a random-basis null. Their implied in-sample ceiling ≈ 0.84
  is an in-sample fit (our in-sample analogue 0.73–0.80, lower out of sample). Our additions: (i) the object analysed is
  the *transfer-prediction residual*, with other lines' control PCs and permuted controls as nulls and a three-way
  control split for noise control (H3a); (ii) five corrections compared (H3b); (iii) methods placed on a ceiling scale
  (H1, which turned out untestable as registered because the registered r_max is internally inconsistent).
- [[huuloc-2026-cell-quota]]: a closed-form cell quota per perturbation direction and a sampling-noise floor; most
  genome-wide CRISPRi knockdowns sit below the detection floor. Unit is cells, not contexts — complementary to RQ2;
  supports our large noise/interaction ratios (Tahoe 2.87, genetic 2.04).
- Systematic benchmark ([[li-2024-systematic-comparison]], preprint of @li-2026-systematic-comparison): 13 models
  plus mean/additive/linear/MLP baselines on 25 datasets and three tasks; additive and linear baselines win on seen
  combinations (p. 9) and no method generalises consistently to unseen cell types (p. 2). Its *diagnostic*,
  non-orthogonal energy decomposition of observed profiles, Y_{c,p} = B_c + G_c + S_p + I_{p1,p2} + C_{c,p} + ε (control
  baseline, shared within-state shift, perturbation-specific, combination non-additivity, cell-state × perturbation
  interaction, and ε an «empirical noise floor» from split halves of control cells; fractions «not constrained to sum
  to 100%»; pp. 30–31), explains why metric families disagree (p. 11; Fig. 5). It gives no noise-corrected
  components or ceiling, no context-count analysis and no correction comparison.

## 4. Data resources and same-line data

- [[pan-2026-global-cell-state]]: source of H1GW (genome-scale CRISPRi Perturb-seq in H1 hESCs, GSE295214); qualitative
  analogue of H4 (target-context data outperform modelling).
- [[huang-2025-x-atlas-orion]]: HCT116 and HEK293T screens (two of the 11 Count-1 cell types) and the noise inputs for
  the OR_HCT116X panel.
- Replogle ([[replogle-2021-genome-scale]], preprint of @replogle-2022-genome-scale): K562 genome-wide and essential,
  RPE1 essential screens (pp. 5–6). Nadig ([[nadig-2024-transcriptome-wide]], preprint of
  @nadig-2025-transcriptome-wide): HepG2 and Jurkat essential screens (p. 7; p. 25). Tahoe-100M ([[zhang-2025-tahoe]],
  preprint of @zhang-2026-tahoe): 50 lines, of which the preprint analyses 47 (p. 19), 379 drugs × 3 doses (p. 5); the
  project's 48-line drug data of the H2a/H2b analysis (local 4.43 M-cell subset) is the preprint's 47 plus NCI-H661
  ([[tahoe-100m]]).
  Jiang [@jiang-2025-pathway-signatures]†: CRISPRi under cytokine stimulation, used for the drug→gene bound (appendix).

## 5. One-paragraph positioning for the introduction

Simple baselines rival learned models for perturbation prediction ([[ahlmanneltze-2025-linear-baselines]],
[[vinastorne-2026-systema]], [[piao-2026-adaptive-contexts]]), and context coverage has been argued to be the binding
constraint ([[dibaeinia-2026-context-scale]]), with dataset-dependent learning curves for scalar readouts
([[ji-2024-universal-prediction]]). What is missing is (a) a measured, matched-ablation break-even in the number of
contexts for transcriptome-wide responses, (b) the corresponding count of public genetic contexts, (c) a test of whether
a closed-form identifiability/noise argument predicts that break-even (it does not: ×6 too pessimistic), and (d) a
noise-controlled account of where zero-shot transfer error lives relative to the destination's own control state, with a
principled correction compared against strong alternatives.

## Evidence base

Fifteen ingested source pages (each with an H1–H4 overlap table): [[ahlmanneltze-2025-linear-baselines]],
[[piao-2026-adaptive-contexts]], [[vinastorne-2026-systema]], [[dibaeinia-2026-context-scale]],
[[ji-2024-universal-prediction]], [[rosenfeld-2021-risks-irm]], [[molina-2026-response-decomposition]],
[[huuloc-2026-cell-quota]], [[pan-2026-global-cell-state]], [[huang-2025-x-atlas-orion]], and — from open preprints
of the paywalled published versions, claims quoted from the preprint only — [[adduri-2025-state]],
[[li-2024-systematic-comparison]], [[replogle-2021-genome-scale]], [[nadig-2024-transcriptome-wide]],
[[zhang-2025-tahoe]]; two sources (†) cited from bibliographic records and logged facts only (VCC 2026, no preprint;
Jiang). Literature guide:
`input/bibliography/literaturguide.md` (102 verified entries); overlap note `input/notes/overlap-molina-zhang.md`.

## What stays open

- STATE's genetic-data gain and the Science Advances decomposition are now verified from the preprints (pp. above);
  the remaining † sources (VCC 2026 — the Cell original is free to read but bot-blocked; Jiang) must be ingested or
  their claims dropped before drafting cites them for specific numbers.
- Preprint-based pages quote the preprint; before submission, check the cited numbers against the published versions
  (STATE: Cell 2026; Li: Science Advances 2026; Tahoe: Cell 2026, retitled; Nadig: Nature Genetics 2025, retitled).
- Prophet's v2 → v3 changes: the paper must quote v3 (the version on disk) only.
