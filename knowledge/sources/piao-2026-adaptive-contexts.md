---
title: "Piao et al. 2026 — AdaPert: adaptive perturbation-conditioned contexts"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: piao-2026-adaptive-contexts
depth: standard
tags: [perturbation-prediction, knowledge-graph, mean-collapse, unseen-perturbation, cross-cell-line, crispri, baselines]
relations:
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "In the only unseen-cell-line test (3 lines → K562), the context-free Perturb Mean beats all ten models: «Perturb Mean (dashed line) achieves the highest Pearson-∆, reflecting that the majority of perturbation effects are conserved across cell lines» (Fig. 13, p. 21; Table 13, p. 22: 0.650 vs AdaPert 0.460)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "«We evaluate on four publicly available genome-wide CRISPRi Perturb-seq datasets spanning the human cell lines K562, RPE1, JURKAT, and HEPG2» (p. 12; Table 4, p. 13)."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "Indirect only: a per-perturbation average over other lines (≈ μ + a_t) wins cross-line, with 59 % of K562 perturbations at Pearson-∆ ≥ 0.7 (Fig. 13, p. 21). Reader's mapping onto the conserved component; the source never decomposes responses."
  - target: zero-shot-transfer-ceiling
    type: mentions
    confidence: inferred
    because: "The source calls Perturb Mean «a strong upper bound when perturbation effects are conserved across cell lines» (Table 13, p. 22). That is a heuristic label for an empirical baseline, not a noise-corrected ceiling."
---

# Learning Adaptive Perturbation-Conditioned Contexts for Robust Transcriptional Response Prediction

## Bibliographic Details

Piao, Yinhua, Hyomin Kim, Seonghwan Kim, Yunhak Oh, Junhyeok Jeon, Sang-Yeon Hwang, Jaechang Lim, Woo Youn Kim,
Chanyoung Park, and Sungsoo Ahn. 2026. "Learning Adaptive Perturbation-Conditioned Contexts for Robust
Transcriptional Response Prediction." In *Proceedings of the 43rd International Conference on Machine Learning
(ICML 2026)*, Seoul. PMLR 306 (p. 1, footer). arXiv:2602.18885v2 [cs.CE], 5 July 2026 (p. 1). KAIST / HITS, Seoul.
Read from `<library>/pdf/piao-2026-adaptive-contexts.pdf` (34 pp.; printed page numbers = PDF page numbers). No code
or data availability statement, no hyperparameters and no seed count are given anywhere in pp. 1–34.
Citekey: @piao-2026-adaptive-contexts.

**Terminology warning.** "Context" in the title means a **perturbation-conditioned subgraph of a knowledge graph**
(z_context, eq. 13, p. 5), not a cellular context. The paper is an unseen-perturbation method paper. Cross-cell-line
transfer appears only in one appendix experiment (App. C, pp. 18–25).

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; 1 Introduction (mean collapse, Fig. 1) | 1–2 |
| — | 2 Related Works (data-driven, knowledge-driven, LLM priors) | 2–3 |
| ✓ F1 | 3.1 Problem definition: conditional autoencoder F: (X_c, p) → X̂_p (eqs. 1–4) | 3–4 |
| — | 3.2 Overview; 3.3 Perturbation-conditioned subgraph extraction (eqs. 5–13, Fig. 3) | 4–5 |
| ✓ F1 | 3.4 Adaptive learning: reconstruction, non-DEG Huber, DEG alignment losses (eqs. 14–21) | 5–6 |
| ✓ F1 | 4.1 Experiment setup (unseen-perturbation, four lines); 4.2 Baselines; 4.3 Metrics (cell-eval) | 6 |
| ✓ F1 | 5.1–5.3 Main results (Tables 1–3) | 6–8 |
| ✓ F1 | 5.4–5.5 Ablations of z_context and L_non; λ_non sweep (Figs. 4–5) | 8 |
| ✓ F1 | 6 Conclusion; Impact statement | 8–9 |
| — | References | 9–11 |
| ✓ F1 | App. A.1 Data statistics (Tables 4–5, Fig. 6) | 12–14 |
| — | App. A.2 STRING knowledge graph (Table 6, Fig. 7); A.3 GenePT descriptions (Table 7) | 14–16 |
| ✓ F1 | App. B.1 Baselines; B.2 Metric definitions (Pearson-∆, PDS, DES, DE-Spearman, direction) | 16–17 |
| — | App. B.3 Inference detail (Fig. 8) | 18 |
| ✓ F1 | App. C Per-line results (Figs. 9–12, Tables 9–12) and **cross-cell-line 3 → K562** (Fig. 13, Table 13) | 18–22 |
| ✓ F1 | App. C Systema evaluation (Figs. 14–18, Tables 14–18) | 22–25 |
| — | App. D DEG-threshold sensitivity (Fig. 19) | 25–26 |
| — | App. E KG variants: source, density, fixed vs adaptive, threshold (Figs. 20–24) | 26–29 |
| ✓ F1 | App. F Mean-collapse analysis (Fig. 25) | 29–30 |
| — | App. G Pathway case study (Figs. 26–28) | 30–33 |
| — | App. H Additional related works | 33–34 |

## Kernthesen

### 1. ⭐ Mean collapse is the main failure mode, and the authors blame the loss, not the data (pp. 1–2, 4; Fig. 1, p. 1; eq. 4)

Models «achieving high correlation by predicting the global average expression rather than perturbation-specific
responses» (p. 1). Because only 1–10 % of genes respond strongly (p. 2), the MSE splits into a DEG term and a much
larger non-DEG term. «Since the second term dominates, minimizing mean squared error is driven mainly by
non-responsive genes, encouraging predictions to shrink toward zero» (p. 4, eq. 4). They say: «We trace mean collapse
to how perturbation effects are modeled and supervised, rather than to data limitations» (p. 2). In Fig. 1 (UQCRB,
114 DEGs), DEG correlation rises from 0.617 (TxPert) to 0.689 (+11 %) (p. 1).

### 2. ⭐ The model conditions on control expression and a perturbation embedding; "context" means a knowledge-graph subgraph (pp. 3–5; Figs. 2–3, 8; eqs. 1–13)

AdaPert is a conditional autoencoder. It learns F: (X_c, p) → X̂_p, with z_c = ENC_θ(X_c) for the control profile
(eq. 1, p. 3), a perturbation embedding, and X̂_p = DEC_φ(z_c, z_p) (eq. 3). The perturbation side has three steps.
(i) STRING PPI node embeddings come from message passing (eqs. 5–6). (ii) A text embedding of the gene's NCBI / GPT-4o
description is used (GenePT; eqs. 7–8, p. 4; App. A.3, p. 15). (iii) Perturbation-conditioned node scores go through
Gumbel-softmax sampling and a threshold T (eqs. 9–12, p. 5). The selected nodes are sum-pooled into z_context
(eq. 13). At inference, «The decoder takes the concatenation of basal state, context vector, and semantic embedding»
(Fig. 8, p. 18). **No cell-line descriptor exists** apart from whatever the control encoder picks up from X_c.

### 3. A three-term loss separates signal from noise using DEG labels from the training data (pp. 5–6; eqs. 14–21)

DEGs are genes with q < 0.05 in a per-gene test on the training data (eq. 14, p. 5). The loss is
L_total = L_recon + λ_non·L_non + λ_align·L_align (eq. 21, p. 6). L_recon is the full-profile MSE (eq. 15). L_non is a
Huber penalty on predicted changes ΔX̂ = X̂_p − X_c for non-DEGs, with δ «set proportional to the empirical standard
deviation of non-DEG effects» (eqs. 16–17, p. 5). L_align is a cosine loss between z_context and a projection of the
DEG-masked true effect y⁽ᵖ⁾ (eqs. 18–20, p. 6).

### 4. ⭐ The main evaluation is unseen perturbations within one line, on four CRISPRi lines (p. 6; App. A.1, pp. 12–13; Table 4)

«All experiments are conducted under the unseen perturbation setting, where perturbations in the test set are not
observed during training» (p. 6). Each line is split at the perturbation level: K562 734 / 82 / 272 train/val/test
perturbations on 5,000 HVGs; RPE1 771 / 308 / 464; Jurkat 745 / 298 / 448; HepG2 826 / 330 / 497, all on 3,352 HVGs
(Table 4, p. 13). The prediction target is «differential expression relative to matched non-targeting controls»
(p. 12). The metrics come from cell-eval: Pearson-∆, PDS (L1, as in the Virtual Cell Challenge), DES, DE-Spearman and
direction match (p. 6; App. B.2, pp. 16–17).

### 5. AdaPert beats TxPert by small margins within a line (Table 1, p. 7; Tables 9–12, pp. 21–22)

Pearson-∆ for AdaPert vs TxPert is 0.619 vs 0.580 (K562) and 0.674 vs 0.655 (RPE1). PDS is 0.711 vs 0.665 (K562)
and 0.663 vs 0.618 (RPE1) (Table 1, p. 7). On RPE1, MORPH has the higher PDS (0.688). Models without a knowledge graph
have PDS 0.50–0.52, «close to chance level» (p. 6). In App. C the margin over the next-best method in Pearson-∆ is
+5.5 % (K562), +2.9 % (Jurkat), +4.3 % (HepG2) and +3.2 % (RPE1) (Figs. 9–12, pp. 19–20). Results are «mean ± std
over all test perturbations» (Table 1). No seeds or significance tests are reported.

### 6. ⭐ The within-line margin over simple baselines is large in K562 and small in RPE1 (Tables 9–12, pp. 21–22)

Pearson-∆ for Context Mean / Linear (GenePT) / AdaPert is 0.387 / 0.463 / 0.620 (K562), 0.398 / 0.420 / 0.500
(Jurkat), 0.385 / 0.400 / 0.469 (HepG2) and 0.629 / 0.633 / 0.676 (RPE1). On RPE1, Context Mean has DE direction match
0.843 vs AdaPert 0.867 and DE Spearman LFC 0.706 vs 0.723 (Table 12, p. 22). "Context Mean" is never defined in the
text. Its numbers match a within-line training-mean baseline (compare ED Fig. 8a of
[[ahlmanneltze-2025-linear-baselines]], about 0.40 / 0.68 on K562 / RPE1). STATE is weak here (0.247 K562, 0.189
Jurkat).

### 7. Gains concentrate in DEG-aware metrics and small-effect perturbations (Tables 2–3, pp. 7–8)

On K562, DES@50 is 0.263 vs TxPert 0.220, and direction match is 0.867 vs 0.843 (Table 2, p. 7). For the smallest
effect tertile, Pearson-∆ is about the same (0.462 vs 0.457, +1.1 %) but DES is +28 % (0.115 vs 0.090), and PDS is
slightly lower (0.737 vs 0.740) (Table 3, p. 8). Medium: +15 % / +28 % / +12 %. Large: +4.0 % / +19 % / +11 % (p. 8).
Tertiles are formed on the number of DEGs (BH-Wilcoxon p < 0.05; p. 12).

### 8. The ablations support both components, and the non-DEG penalty matters most for weak effects (Figs. 4–5, p. 8)

«Removing z_context consistently degrades performance across all comparisons» (p. 8). Removing L_non hurts
small- and medium-effect perturbations most. Raising λ_non from 0 to 0.01 helps those groups, larger effects prefer
a smaller λ_non, and λ_non = 0.1 degrades every group (Fig. 5, p. 8). These ablations are K562 only and within-line.

### 9. ⭐⭐ In the only unseen-cell-line test, a context-free per-perturbation mean beats every model (Fig. 13, p. 21; Table 13, p. 22)

Models are trained on RPE1 + HepG2 + Jurkat and tested on held-out K562 (p. 18). Pearson-∆: **Perturb Mean 0.650**,
Linear-GenePT 0.515, Linear-OneHot 0.442, AdaPert 0.460, TxPert 0.405, STATE 0.309, scVI 0.129, scGPT 0.078,
CPA 0.027, GEARS −0.002, MORPH −0.113 (Table 13, p. 22). Perturb Mean also wins Discrim. L1 (0.909 vs AdaPert 0.692),
Overlap@50 (0.348 vs 0.167), direction match (0.837 vs 0.703) and DE Spearman LFC (0.817 vs 0.416). It loses only on
MAE-∆ (0.065 vs 0.035). Perturb Mean is a «perturbation-aware training-set average» (Fig. 13, p. 21), i.e. the mean
response of the same perturbation in the three source lines. It uses no K562 information.

### 10. ⭐ The authors attribute this to conserved effects and say the model is not built for cross-line transfer (pp. 21–22)

«Perturb Mean (dashed line) achieves the highest Pearson-∆, reflecting that the majority of perturbation effects are
conserved across cell lines (59% of test perturbations have Perturb Mean ∆ ≥ 0.7)» (Fig. 13, p. 21). Linear-GenePT is
the best learned model, «likely benefiting from LLM-derived embeddings that encode cell-type-independent functional
priors» (p. 21). AdaPert is «the best performance among deep-learning models» (p. 21). «We note that A DA P ERT is
designed for the unseen-perturbation setting; the cross-cell-line setting introduces additional challenges related to
cell-state variation and imbalanced perturbation conservation across training cell lines» (p. 21). Table 13 calls
Perturb Mean «a strong upper bound when perturbation effects are conserved across cell lines» (p. 22).

### 11. ⭐ The Systema cross-line table reports AdaPert as best but leaves out Perturb Mean and TxPert (Fig. 18, p. 24; Table 18, p. 25)

Under the Systema protocol, AdaPert has P∆(Pert) 0.783, P∆20(Pert) 0.733 and centroid accuracy 0.689, against STATE
0.718 / 0.658 / 0.622 (Table 18, p. 25). The text concludes that this «demonstrates A DA P ERT's ability to generalize
perturbation predictions to entirely unseen cellular contexts with different gene-regulation programs» (Fig. 18,
p. 24). Table 18 lists ten methods. **Perturb Mean and TxPert are not in it**, so the baseline that won Table 13 is
not scored under Systema. The only mean-type entry is "Matching-Mean" (P∆(Pert) 0.132, P∆(Ctrl) 0.053), which is not
perturbation-specific (Kernthese 12). STATE is also not second on P∆20(Pert): Linear has 0.666 (Table 18). **Table 18
does not agree with Table 13** for the same setting. Control-referenced Pearson-∆ is 0.509 vs 0.460 (AdaPert),
0.358 vs 0.027 (CPA), 0.254 vs −0.002 (GEARS) and 0.377 vs 0.309 (STATE). Unlike every within-line table, the
perturbation-referenced scores exceed the control-referenced ones (AdaPert 0.783 vs 0.509). The paper does not
explain either point.

### 12. Within-line Systema scores agree with the main tables (Figs. 14–17, pp. 23–24; Tables 14–17, pp. 24–25)

AdaPert has the best perturbation-referenced Pearson-∆ on every line: K562 0.534, Jurkat 0.345, HepG2 0.394,
RPE1 0.379 (pp. 23–25). STATE and Matching-Mean are negative on K562 (−0.054, −0.119; Fig. 14, p. 23). The
perturbation-referenced values are below the control-referenced ones (K562 0.534 vs 0.618; Table 14, p. 24), as
expected when the shared response is removed. "Matching-Mean" has the same control-referenced values as "Context
Mean" in Tables 9–10 (K562 0.3865 vs 0.3867; Jurkat 0.3975 vs 0.3975), so it is the within-line training mean.

### 13. The best knowledge graph has wide gene coverage and moderate density (App. E, Figs. 20–23, pp. 27–28)

STRING PPI is best (Pearson-∆ 0.6188), then GO (0.6047). Co-expression graphs built from control cells (about 5K
genes) lose to curated graphs (18.2K nodes) at every density (Fig. 20, p. 27). Top-10 edges per gene is best, and a
perturbed-genes-only graph beats top-20 and top-50 (Fig. 22, p. 28). Learned subgraphs beat fixed 1-hop and 2-hop sums
(Fig. 23, p. 28).

### 14. The DEG definition does not drive the results (App. D, Fig. 19, p. 26)

Stricter FDR or |LFC| cutoffs shrink training DEGs (K562 641 → 229) but change Pearson-∆ by < 1.3 % (K562), 2.1 %
(Jurkat) and 0.7 % (RPE1) (p. 26).

### 15. The mean-collapse reduction is a few percent of MSE (App. F, Fig. 25, p. 30)

Against TxPert, rank-weighted MSE improves by 3.0 % on DEGs and 1.8 % on non-DEGs. Unweighted MSE improves by 4.6 %
and 4.7 % (Fig. 25, p. 30). The paper reads the order Linear → TxPert → AdaPert as progressive mitigation (p. 30).

### 16. Pathway-level gains are reported for K562-relevant programs (App. G, Figs. 26–28, pp. 31–33)

AdaPert beats TxPert on 9 of 12 K562-relevant Hallmark pathways (Fig. 27, p. 32). For the HIRA knockdown, predicted
and true NES correlate at r = 0.53 over 44 pathways (Fig. 28, p. 33). The authors claim AdaPert «preferentially
captures cell-type-specific biological programs» (p. 33). This is a within-line claim about K562 and is not tested
across lines.

## Figures and tables

- Fig. 1 (p. 1): mean collapse on UQCRB, TxPert vs AdaPert DEG correlation.
- Figs. 2–3 (pp. 3–4): architecture and subgraph extraction. Fig. 8 (p. 18): inference pipeline.
- Tables 1–3 (pp. 7–8): within-line main results, DEG-aware metrics, effect-size strata. Figs. 4–5 (p. 8): ablations.
- Tables 4–5 (p. 13): splits and effect-size strata for the four lines.
- Tables 9–12 (pp. 21–22): all metrics per line, **including Context Mean and linear baselines**.
- **Fig. 13 (p. 21) and Table 13 (p. 22): cross-cell-line 3 → K562, Perturb Mean best.**
- Figs. 14–18, Tables 14–18 (pp. 23–25): Systema evaluation. **Table 18 (p. 25) is the cross-line Systema table
  without Perturb Mean or TxPert.**
- Figs. 19–28 (pp. 26–33): DEG thresholds, KG variants, mean-collapse error split, pathway case study.

## What this source does not address

- **No context count beyond one split.** Cross-line transfer is a single experiment: 3 training lines → K562
  (pp. 18, 21). There is no other held-out line, no leave-one-line-out rotation, no sweep over the number or diversity
  of training lines, and no theory of how many contexts are needed.
- **No Tahoe and no chemical perturbations.** Only the four CRISPRi lines are used (p. 12).
- **No ceiling.** There is no noise-corrected or split-half reliability bound and no variance decomposition. "Upper
  bound" (p. 22) is a label for the Perturb Mean baseline.
- **No few-shot arm.** No target-line perturbation measurements are used in the cross-line test, so the paper has no
  zero-shot vs calibrated contrast.
- **No control-state descriptor beyond the encoder.** Control expression enters only through ENC_θ(X_c) (eq. 1, p. 3).
  There are no control PCs, no deflation, no common-response removal and no residual analysis. Systema appears only
  as an evaluation protocol (p. 22).
- **No uncertainty.** No seeds, confidence intervals or significance tests (Tables 1–18). The mean ± std is over test
  perturbations (Table 1, p. 7).
- **No limitations section** and no code or data statement.

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 2, 4, 6, 9, 10, 11** (and 12 for the Systema metric, 5 for the size of modelling gains). The
hypotheses are those of the plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top;
Amendments A1–A3 at the end).

### What this focus needs beyond the Kernthesen

1. **How it adapts to an unseen line.** Not by design. The only route for cell-line information is the control
   encoder z_c = ENC_θ(X_c) (eq. 1, p. 3; Fig. 8, p. 18). The paper does not say explicitly that K562 control cells
   are the encoder input in the cross-line test, but the architecture requires a control profile (eq. 3), so this is
   the most likely reading (inferred). No target perturbation measurements are used. The cross-line test is
   therefore **zero-shot in the cell line but not in the perturbation**: the same perturbations are measured in the
   source lines, which is what makes the Perturb Mean baseline possible (p. 21). The K562 test set has n = 1,087
   perturbations (Fig. 18, p. 24), i.e. the whole K562 set (734 + 82 + 272 = 1,088 entries including the control
   entry; Table 4, p. 13).
2. **Context count.** C = 3 training lines, 1 destination, one split (p. 18). There is no ablation over the number or
   diversity of lines and no theory. The genes differ by line (5,000 HVGs in K562, 3,352 in the others; Table 4,
   p. 13), and the paper does not say which gene space is used for the cross-line test.
3. **Baselines and margins.** The within-line learned-model gain over Context Mean is +0.23 Pearson-∆ in K562 but
   only +0.05 in RPE1 (Tables 9, 12). Across lines, the context-free Perturb Mean is +0.19 Pearson-∆ above AdaPert and
   +0.40 above it in DE Spearman LFC (Table 13, p. 22). So every context-aware model in the paper does worse than
   simply averaging the perturbation's measured response over the three source lines. This includes STATE, MORPH and
   TxPert, which have explicit context machinery.
4. **Internal inconsistencies to keep in mind when citing.** (a) The cross-line Systema Table 18 omits Perturb Mean
   and TxPert, and its control-referenced numbers do not match Table 13 (Kernthese 11). Cite Table 13 for the
   cross-line result, not Fig. 18's "entirely unseen cellular contexts" claim. (b) Table 1 and Table 9 give TxPert
   K562 Pearson-∆ as 0.580 and 0.588. (c) Test counts differ slightly between Table 4 and Figs. 10–12 (e.g. RPE1 464
   vs 461). (d) p. 6 says all four lines are «from Replogle et al.», while p. 12 says only K562 and RPE1 are. Jurkat
   and HepG2 (@nadig-2025-transcriptome-wide) are not cited. (e) Table 4 calls the screens «genome-wide», but the
   K562 set with 1,087 perturbations and 5,000 HVGs matches the GEARS-processed essential screen used by
   [[ahlmanneltze-2025-linear-baselines]] (ED Fig. 1, PDF p. 8 there).
5. **Literaturguide pre-read note** ("H4 (H1-adjacent): AdaPert, ICML 2026; mean collapse and noise;
   unseen-perturbation focus"): correct as far as it goes. It misses the one result that matters most for this
   project: the cross-line appendix experiment in which the context-free Perturb Mean beats all models (Kernthese 9).
   That result is H1 / H2a-adjacent rather than H4.
6. **Reader's interpretation (not in the source).** The cross-line result fits the project's account. With C = 3
   training lines, Proposition 1 caps a learned context model at rank ≤ C − 1 = 2 of the interaction. That is far
   below the project's empirical break-even C*_emp ≈ 37 contexts on Tahoe (H2a). A context-free transfer of the
   conserved component is the strongest zero-shot predictor in exactly that regime. The authors explain the shortfall
   differently: «cell-state variation and imbalanced perturbation conservation across training cell lines» and a
   model «designed for the unseen-perturbation setting» (p. 21). They never mention the number of training lines.
   The source has nothing on the rejected linear-Gaussian C* formula (H2b).

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partial (qualitative only)** | Cross-line, a context-free per-perturbation mean from other lines (≈ μ + a_t) scores Pearson-∆ 0.650, with 59 % of perturbations at ≥ 0.7. All learned models score lower (0.460 AdaPert, 0.515 Linear-GenePT). So published models do not even recover the transferable part, which is the "model-limited" pattern H1 is designed to detect. There is no r_max, no noise correction, no split-half and no error split. "Upper bound" is a heuristic label. | Fig. 13 p. 21; Table 13 p. 22 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **partial (adjacent, genetic, C = 3)** | Context-conditioned models lose to a context-free mean at C = 3 training lines (Table 13), which points the same way as a negative context gain below break-even. But this is not a learned correction on top of a context-free model, not Tahoe, not a C sweep, and has no zero crossing. | p. 18; Table 13 p. 22 |
| **H2b** (theory predicts the break-even within ×2) | **none** | No sample-complexity theory, no d, σ²_ε or Var(i). | — |
| **H2c** (genetic C*_min > C_public) | **none (tangential)** | It uses the same four public CRISPRi lines (3 for training) but counts no C_public and estimates no C*_min. | p. 12; Table 4 p. 13 |
| **H3a** (transfer residual on destination top-8 control PCs > other-line / permuted PCs) | **none** | Control expression is only an encoder input. There are no control PCs, no transfer residual and no projection. | eq. 1 p. 3; Fig. 8 p. 18 |
| **H3b** (M3 beats five PC-deflation / common-response corrections on PDS) | **none (metric overlap only)** | Systema perturbation-referenced Pearson-∆ and centroid accuracy are used for evaluation only. No correction is compared. PDS (L1, VCC definition) is reported. | p. 22; Tables 14–18 pp. 24–25; App. B.2 p. 17 |
| **H4** (adding same-line genome-wide data beats the largest modelling gain) | **none (tangential)** | No data-addition experiment. Within-line modelling gains over the next-best model are +2.9 to +5.5 % relative Pearson-∆ (Figs. 9–12), a useful yardstick for "largest modelling gain". | Table 1 p. 7; Figs. 9–12 pp. 19–20 |

### Direct quotes (for the draft)

> "Perturb Mean (dashed line) achieves the highest Pearson-∆, reflecting that the majority of perturbation effects are
> conserved across cell lines (59% of test perturbations have Perturb Mean ∆ ≥ 0.7)." (Fig. 13, p. 21)

> "We note that A DA P ERT is designed for the unseen-perturbation setting; the cross-cell-line setting introduces
> additional challenges related to cell-state variation and imbalanced perturbation conservation across training cell
> lines." (Fig. 13, p. 21)

> "Perturb Mean is a perturbation-aware training-set average that serves as a strong upper bound when perturbation
> effects are conserved across cell lines." (Table 13, p. 22)

> "We trace mean collapse to how perturbation effects are modeled and supervised, rather than to data limitations."
> (p. 2)

> "All experiments are conducted under the unseen perturbation setting, where perturbations in the test set are not
> observed during training." (p. 6)

### Examples & illustrations (for later drafting)

- 3 lines → K562: Perturb Mean 0.650 vs best deep model 0.460 vs STATE 0.309 vs MORPH −0.113 Pearson-∆ (Table 13,
  p. 22). This is a clean published case of a context-free transfer beating every context-aware model at small C.
- RPE1 within-line: Context Mean 0.629 vs AdaPert 0.676 (Table 12, p. 22). Even within a line, the mean baseline
  carries most of the score on template-heavy lines (compare RPE1's 65 % template share in
  [[molina-2026-response-decomposition]]).

### Boundary of this focus

This focus does not cover the knowledge-graph design, the loss engineering or the pathway case study (Kernthesen 3,
7–8, 13–16), apart from their size as modelling gains. A reader looking for a context-count curve, a ceiling, control
PCs or a data-scaling experiment will not find one in this source.

## Mentioned entities

- Concepts: [[simple-perturbation-baselines]], [[perturbation-response-decomposition]], [[zero-shot-transfer-ceiling]]
- Datasets: [[replogle-nadig-crispri-collection]]

## Connections

- Supports: [[simple-perturbation-baselines]] (the cross-line Perturb Mean result, Table 13, p. 22).
- Consistent with (inferred): [[perturbation-response-decomposition]]: conserved effects dominate cross-line
  (Fig. 13, p. 21), as in [[molina-2026-response-decomposition]] (zero-shot predictions align with β, never with γ).
- Mentions: [[zero-shot-transfer-ceiling]] ("upper bound" used heuristically, p. 22).
- Same within-line benchmark lineage as [[ahlmanneltze-2025-linear-baselines]] (GEARS-style K562 / RPE1 splits). The
  paper does not cite it.
- Cites (not yet ingested): @adduri-2026-state (as Adduri et al. 2025), @wenkel-2026-txpert (as Wenkel et al. 2025),
  @replogle-2022-genome-scale. It uses the Systema protocol (@vinastorne-2026-systema) without listing it in the
  references (pp. 9–11, 22).
