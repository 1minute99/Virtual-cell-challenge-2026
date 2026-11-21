---
title: "Molina & Zhang 2026 — Perturbation response decomposition"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: molina-2026-response-decomposition
depth: deep
tags: [perturbation-prediction, variance-decomposition, crispri, cross-cell-line, zero-shot, split-half, control-subspace]
relations:
  - target: perturbation-response-decomposition
    type: supports
    confidence: extracted
    because: "«We introduce a decomposition framework that explicitly separates transcriptional responses into global, perturbation-specific, cell-line-specific, and perturbation-by-cell-line interaction components» (p. 1)."
  - target: control-subspace-projection
    type: supports
    confidence: extracted
    because: "«we fitted PCA on centered non-targeting control cells from each cell line and projected the template and residual responses onto the resulting control subspace» (p. 35, Supp. Note 3)."
  - target: zero-shot-transfer-ceiling
    type: builds-on
    confidence: inferred
    because: "Supplies the noise-corrected variance components from which the project's r_max is built; the source itself never states a ceiling (no ceiling anywhere in pp. 1–44)."
---

# Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts

## Bibliographic Details

Molina, Alexis, and Xinyi Zhang. 2026. "Perturbation response decomposition enables biologically aligned generalization
to unseen perturbations and cellular contexts." bioRxiv preprint, posted 27 July 2026 (AITHYRA, Vienna).
doi:10.64898/2026.07.24.740459. CC BY 4.0. Code: github.com/xinyizhanglab/perturbation-decomposition (p. 28).
Read from `<library>/pdf/molina-2026-response-decomposition.pdf` (44 pp.; printed page numbers = PDF page numbers).
Citekey: @molina-2026-response-decomposition.

**Figure numbering.** The PDF numbers its supplementary figures 5–9; the text calls them Supplementary Fig. 1–5
(Supp. Fig. 1 = Fig. 5, p. 30; Supp. Fig. 2 = Fig. 6, p. 31; Supp. Fig. 3 = Fig. 7, p. 31; Supp. Fig. 4 = Fig. 8,
p. 32; Supp. Fig. 5 = Fig. 9, p. 33). Below, "Fig. N" is the PDF number.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; 1 Introduction | 1–2 |
| ✓ F1 | 2.1 Decomposition of perturbation response reveals predictability of response components (Fig. 1, p. 5) | 2 |
| — | 2.2 Response alignment makes gene coessentiality predictive of perturbation-specific effects (Fig. 2, p. 6) | 2–3 |
| ✓ F1 | 2.3 Mapping perturbation-aligned axes across cell lines enables generalization to unseen cell lines (Fig. 3, p. 7) | 3–4 |
| ✓ F1 | 2.4 Our decomposition framework reflects the predictability of the response components (Fig. 4, p. 8) | 4 |
| ✓ F1 | 3 Discussion (incl. limitations) | 9 |
| ✓ F1 | 4.1 Data resources and preprocessing (four CRISPRi screens, genome-wide K562, Norman, Tahoe, DepMap, gene selection) | 10–11 |
| ✓ F1 | 4.2 Response construction and evaluation (pseudobulk δ, splits, perturbed-reference Pearson) | 11–12 |
| — | 4.3 DepMap response mapping and within-cell-line prediction | 12–13 |
| — | 4.4 Model-specific predictors (scGPT, MORPH, linear probing, STATE) | 13–15 |
| ✓ F1 | 4.5 Cross-cell-line prediction tasks (one-/three-source, double-unseen, calibration, basis coverage) | 15–16 |
| — | 4.6 Combinatorial perturbation benchmark (Norman) | 16–17 |
| ✓ F1 | 4.7 Chemical perturbation transfer in Tahoe-100M | 17–19 |
| ✓ F1 | 4.8 Response decomposition (template projection, balanced additive decomposition, split-half partition) | 19–22 |
| ✓ F1 | 4.9 Analysis of decomposed spaces (enrichment, essential vs non-essential, control vs perturbation subspace, overlap) | 22–26 |
| ✓ F1 | 4.10 Prediction attribution | 26–27 |
| ✓ F1 | 4.11 Statistical analysis and reporting (50 split-half resamples) | 27–28 |
| — | References | 29 |
| ✓ F1 | Supplementary Figs. 5–9 (= Supp. Fig. 1–5) | 30–33 |
| ✓ F1 | Supp. Note 1: sensitivity of decomposition to template removal (Table 1) | 34 |
| ✓ F1 | Supp. Note 2: template and residual geometry (Tables 2–4) | 34–35 |
| ✓ F1 | Supp. Note 3: projection onto control expression subspaces (Tables 5–10) | 35–37 |
| — | Supp. Notes 4–8: DepMap similarity, Ridge/MLP maps, rank ablation, learned projections, MORPH stages (Tables 11–16) | 36–39 |
| ✓ F1 | Supp. Note 9: direct source→target with STATE and MORPH (Table 17) | 39–40 |
| ✓ F1 | Supp. Notes 10–11: three-source and one-source transfer, double-unseen (Tables 18–19) | 40 |
| ✓ F1 | Supp. Notes 12–13: target properties; response-basis gaps (Tables 20–21) | 40–41 |
| ✓ F1 | Supp. Note 14: target-only and within-cell-line reference baselines (Table 22) | 41 |
| ✓ F1 | Supp. Note 15: target calibration saturation | 42 |
| ✓ F1 | Supp. Note 16: prediction decomposition (Tables 23–29) | 42–44 |

## Kernthesen

### 1. ⭐⭐ Perturbation responses split into four orthogonal components with measurable, noise-corrected shares (p. 2; Fig. 1c, p. 5; Supp. Table 1, p. 34)

The authors write the pseudobulk response as δ(c,p,g) = μ(g) + α_c(g) + β_p(g) + γ(c,p,g): a global mean, a
cell-line effect, a perturbation effect conserved across lines, and a line × perturbation interaction (Fig. 1a, p. 5;
eqs. pp. 20). On the four-line Replogle–Nadig tensor, after split-half noise correction: «reproducible variance was
distributed across the cell-line-specific component (the cell-line template; 27.8%), the conserved
perturbation-specific component (29.4%), and the cell-line-by-perturbation interaction (23.5%), with the remaining
19.3% attributable to measurement noise» (p. 2). Supp. Table 1 gives 27.9 / 29.4 / 23.5 / 19.2 % (p. 34). The
"template" share is μ + α_c together (Fig. 1c legend "Cell-line template (μ + α)", p. 5; S_template = S_M + S_A,
p. 21). These are shares of total observed response **energy** Q (uncentred sums of squares), not variances
(f_Z = S_Z / Q, p. 22).

### 2. ⭐⭐ A cross-half inner product separates reproducible signal from sampling noise (pp. 21–22, 28)

Within each line × perturbation condition the cells are split into two disjoint halves, the full decomposition is
recomputed on each half, and the reproducible signal of each component is the cross-half Frobenius inner product
S_Z = ⟨Z⁽¹⁾, Z⁽²⁾⟩_F (p. 21). Noise is the remainder S_noise = Q − S_signal, with Q the mean of the two halves' energies
(p. 21). Estimates are averaged over 50 random splits (p. 28). The assumption is stated: «Under the assumption that
sampling noise was independent between halves, the cross-half product retained response structure reproduced in both
estimates while suppressing noise specific to either half» (p. 21). Its interpretation: «The reproducible portion of
γc,p,g represents target-specific interaction structure that is consistently detected across independent cell subsets
and is therefore, in principle, available for model recovery» (p. 22).

### 3. ⭐ The interaction survives every way of removing the template (p. 2; Supp. Note 1, p. 34)

The residual decomposition is computed three ways. ANOVA accounting and centroid removal give almost the same result:
conserved 40.7–40.8 %, interaction 32.6 %, noise 26.6–26.7 % of residual energy. Projective template removal, the
stricter method used in the main analysis, gives conserved 35.0 %, interaction 27.7 %, noise 37.4 % (Table 1, p. 34).
«Projective removal reduced the total residual signal and increased the relative contribution of measurement noise …
However, it did not eliminate either the conserved perturbation component or the interaction component» (p. 34).

### 4. Template prominence varies widely across lines; residual dimensionality does not (p. 2; Supp. Note 2, pp. 34–35)

The template captures about 23 % (K562), 65 % (RPE1), 49 % (HepG2) and 28 % (Jurkat) of full response energy (Table 2,
p. 35). The effective rank of the full response is 10.8–25.6, but after template removal it is 31.1–37.1 in all four
lines, with 5–7 PCs needed for 50 % of residual variance (Table 3, p. 35). The authors conclude that apparent
differences in response complexity «largely reflect template magnitude rather than differences in core perturbation
biology» (p. 2), and that standard losses and full-response metrics therefore reward template recovery (pp. 2, 35).

### 5. The template is a growth-arrest / stress program, not a lineage program (p. 2; Fig. 1f–g, p. 5; pp. 35)

Template enrichment shows depletion of proliferation programs (E2F, MYC, G2M) and enrichment of p53 / stress pathways,
with line-specific differences (Fig. 1f). Proliferation-shutdown cosines are 0.270–0.477 and p53-arrest cosines
0.194–0.382 (Table 4, p. 35). Baseline identity axes from non-targeting controls explain only 0.9 % (HepG2) to 7.4 %
(K562) of template variance (Fig. 1g, p. 5). Cell-type signatures are «largely orthogonal to the average perturbation
response» (p. 35).

### 6. Essential perturbations load about twice as strongly on the template (p. 2; Fig. 1h; pp. 22–23)

In the genome-wide K562 screen, essential perturbations have about two-fold higher template contribution than
significant non-essential ones (Mann–Whitney P = 2.8e-71, Fig. 1h, p. 5). The difference persists after matching
response norms (Fig. 5g, p. 30), and the essential and non-essential templates have concordant Hallmark enrichment
(Fig. 5h). So essentials differ in template strength, not template identity (p. 30).

### 7. ⭐⭐ Control PCs contain the template but not the template-removed residual (p. 2; Fig. 1i; Supp. Note 3, pp. 35–37)

PCA on centred non-targeting control cells of each line: the first 100 control PCs capture 92.0 / 98.3 / 96.8 / 95.2 %
of the K562 / RPE1 / HepG2 / Jurkat template (Table 6, p. 36), but only 36.9 / 34.5 / 39.3 / 38.6 % of residual
variance (Table 7, p. 36). At k = 5 the control PCs capture 13.4 / 7.1 / 11.2 / 11.9 % of residual variance, against
34.4 / 32.1 / 34.9 / 31.2 % for the residual's own PCs (Table 10, p. 37). The only null is a random basis: «Random
bases captured less than 0.5% in all cases» on full deltas and «approximately 0.3%» on residuals (p. 36). The
conclusion is that «control variation captures the template well, while the residual response requires
perturbation-induced directions that are poorly represented by the leading control manifold» (p. 36).
*Inconsistency in the source:* Table 7 lists the first 5 control PCs at 10.9 / 7.4 / 10.5 / 11.1 %, which differs from
Table 10's k = 5 row.

### 8. The in-control part of the residual is compact, the outside part is high-rank (p. 36; Table 8, p. 37)

The residual inside the 100-PC control subspace needs about 3 PCs for 50 % of its variance. The part outside (58–65 %)
needs 46–51 PCs. The authors read this as «a few shared state-like modes» versus «many perturbation-induced programs»
(pp. 36–37). In the genome-wide K562 variant, the "control" subspace is built from non-significant perturbations
(energy-test p ≥ 0.5) as a proxy. A Marchenko–Pastur threshold σ(√n + √d) with σ = 1/√n_cells then separates
structured from noise variation in the orthogonal part (p. 25).

### 9. Raw DepMap coessentiality is not aligned with the residual response geometry (pp. 2–3; Supp. Note 4, pp. 36–38)

Spearman correlations between DepMap cosine similarity and response similarity are 0.210–0.366 on full deltas but
−0.021 to −0.038 after template removal (Table 11, p. 37). Even the top DepMap-similarity decile has mean residual
correlations of only 0.0099–0.0171 (Table 12), and k-means clusters in DepMap space are only weakly coherent
(Table 13, p. 38).

### 10. A supervised linear or MLP map from DepMap matches MORPH within a line (p. 3; Supp. Note 5, p. 38)

Ridge (α = 10) and a two-layer MLP on DepMap PCA50 reach perturbed-reference Pearson 0.311–0.494. MORPH reaches
0.305–0.459, and is significantly lower in RPE1 (Table 14, p. 38; Fig. 2a, p. 6). «the recoverable signal arises
primarily from combining an external functional prior with response alignment, rather than from conditional
generative architecture» (p. 3).

### 11. The useful prior signal is low-rank but not in the top DepMap axis (p. 3; Supp. Notes 6–7, pp. 38–39)

In K562, PC1 alone gives 0.082, 50 PCs give 0.366, 200 PCs give 0.377, and the full rank gives 0.298 (Table 15,
p. 38). Performance saturates between 20 and 50 PCs in all lines (Fig. 2d). The learned Ridge map rotates DepMap into
response geometry: 93.6 % of the leading predictive mode's energy is in the residual, and the post-map similarity
reaches Spearman ≈ 0.25 with residual similarity (p. 39).

### 12. MORPH's encoder discards predictive information that the raw prior had (p. 3; Supp. Note 8, p. 39)

Linear probes on MORPH representation stages: raw DepMap gives 0.308–0.471, the perturbation encoder 0.249–0.419, and
the post-attention layers do not recover the loss (Table 16, p. 39). Adding the cell embedding before attention adds
nothing (p. 39).

### 13. For combinations, MORPH relies on measured singles and loses to DepMap pair features when both genes are unseen (p. 3; Fig. 2c; pp. 16–17)

On Norman (72 DepMap-covered test doubles: 13 / 44 / 15 in the 0/2, 1/2 and 2/2 unseen strata, p. 16), MORPH matches
additive baselines when both singles were seen. It underperforms Ridge/MLP on summed DepMap pair features when both are
unseen (Fig. 2c, p. 6).

### 14. ⭐⭐ Existing models do not transfer the perturbation-specific signal zero-shot (p. 3; Supp. Note 9, pp. 39–40; Tables 28–29)

«these models cannot predict perturbation-specific components in held-out cell lines beyond the cell-line template»
(p. 3). Direct K562↔RPE1 transfer gives STATE 0.046 / 0.022 and MORPH −0.013 / −0.041 perturbed-reference Pearson.
With 30 % of target perturbations added, STATE rises to 0.532 / 0.458 and MORPH to 0.226 / 0.116 (Table 17, p. 40).
In three-source zero-shot, STATE's corr(δ̂, β) is 0.043–0.120 and corr(δ̂, γ) is −0.021 to 0.014 (Table 28, p. 44).

### 15. ⭐⭐ A context-free DepMap map transfers zero-shot, and more source lines help (pp. 3, 15; Fig. 3a, p. 7; Supp. Notes 10–11, p. 40)

Three-source → one-target zero-shot, with no cell-line descriptor: MLP 0.226 / 0.271 / 0.340 / 0.241 and Ridge
0.118 / 0.208 / 0.287 / 0.167 for K562 / RPE1 / HepG2 / Jurkat (Table 18, p. 40). One-source → one-target is lower,
e.g. MLP 0.128 for RPE1→K562 (Table 19; all 12 directed pairs in Fig. 8a, p. 32). «Incorporation of more source cell
lines improved generalization performance for both models and in all held-out cell lines» (p. 3). For MLP the 1→3
improvement is significant in every target (Welch's t over five folds, p < 10⁻¹⁰; Fig. 3a legend, p. 7). This is the
only context-count manipulation in the paper. The model is explicitly context-free: «Ridge and MLP received no
explicit cell-line descriptor and consequently estimated a source-averaged response mapping» (p. 15).

### 16. Double-unseen transfer stays positive but is lower (p. 3; Fig. 3b; p. 40)

With 30 % of shared perturbation genes removed from all sources, MLP falls to 0.164 / 0.231 / 0.288 / 0.194 and Ridge
to 0.102 / 0.211 / 0.276 / 0.161 (Table 18, p. 40). Seeing the same gene in source lines therefore still carries
response information (p. 40).

### 17. Source response bases span little of the target residual, and this does not predict transfer (p. 3; Fig. 3c; Supp. Note 13, p. 41)

At rank 20, source residual bases capture 19.8–29.8 % of target residual variance (gap 0.702–0.802; Table 21, p. 41).
Coverage does not track MLP transfer: Spearman ρ = 0.24, p = 0.46 (Fig. 8d, p. 32). «response-program transfer and
DepMap-to-response prediction are distinct mechanisms» (p. 41).

### 18. Target-line properties, not source identity, set the difficulty (Supp. Note 12, pp. 40–41)

The zero-shot ordering HepG2 > RPE1 > Jurkat > K562 follows target response geometry. Template-dominated, lower-rank
targets are easier because «a larger fraction of its response can be explained by shared or conserved components»
(pp. 40–41; Table 20, p. 41). Table 20 also lists an "interaction variance" per target (7.19 / 10.29 / 7.40 / 5.69)
without units.

### 19. ⭐⭐ Zero-shot predictions align with the conserved component but never with the interaction (p. 4; Fig. 4b, p. 8; Table 25, p. 43)

«Models aligned positively with the conserved component but not with the target interaction» (Fig. 4b legend, p. 8).
For three-source zero-shot DepMap models, corr(δ̂, β) is 0.049–0.110 (Ridge) and 0.184–0.221 (MLP). The "R²(β)"
column reads 0.10–0.13 (Ridge) and 0.25–0.39 (MLP), and «correlations with γ(c, p, g) were negative across all targets»
(p. 43). R²(γ) ranges from −0.02 to −0.29 (Table 25). *Inconsistency:* the main text reports the R²(β) range as
«r=0.25–0.39 versus 0.10–0.13» (p. 4), and R²(β) / R²(γ) are not defined in the Methods.

### 20. ⭐ Control expression, as input or for fine-tuning, does not supply the missing context (p. 4; Supp. Note 16, pp. 42–44)

Within a line, control-expression descriptors alone give perturbed-reference r between −0.031 and 0.081, and adding
them to DepMap changes nothing (Table 24, p. 42). Across lines, concatenating target or source control descriptors with
DepMap leaves performance and component alignment «essentially unchanged» (Table 26, p. 43). MORPH fine-tuned on target
controls stays near zero (Table 29, p. 44). «The control cell state therefore is not predictive of the different
perturbation response components in the held-out cellular context» (p. 4).

### 21. ⭐⭐ Target perturbation measurements recover the interaction; pooling with sources can hide it (p. 4; Fig. 4c–d; Table 27, p. 44)

Training only on 30 % of target perturbations flips corr(δ̂, γ) from negative to positive in K562, HepG2 and Jurkat
(e.g. K562 −0.031 → 0.039), while corr(δ̂, β) «changed little or decreased». «Pooling source and target data produced
little improvement relative to zero shot … source responses can dominate simple pooled models and obscure the weaker
target calibration signal» (p. 43). For the DepMap MLP, pooled source + 30 % target lifts HepG2 0.340 → 0.386 but
leaves RPE1 flat (Table 18, p. 40).

### 22. ⭐ The interaction is "learnable from currently available data", but only with target measurements (p. 4; Fig. 4d)

«comparable performance for the interaction component can be achieved when disjoint sets of perturbations in the
source and target cell lines are used in training, indicating that there are learnable and generalizable patterns of
cell line-by-perturbation interactions» (p. 4). The paper also says the interaction «cannot be predicted without
information from the target perturbation context» (p. 1) and «is hard to predict in unseen cell lines but is learnable
given the current data» (p. 4). The learnability claim refers to target-calibrated training, not zero-shot.

### 23. ⭐ About 100–150 target perturbations calibrate the target geometry (Supp. Note 15, p. 42; Fig. 8c, p. 32)

Ridge coordinate mapping rises from near zero with one calibration perturbation to perturbed-reference r ≈ 0.38–0.40
after 100–150 target perturbations, then gains «less than 0.01 improvement per 50 added perturbations» (p. 42). The
authors attribute this to a low-dimensional target response space (p. 42).

### 24. Within-line training is an "upper reference", not a strict bound (Supp. Note 14, p. 41; Fig. 8e, p. 32)

Full within-cell-line training gives MLP 0.356 / 0.458 / 0.372 / 0.318 and Ridge 0.343 / 0.469 / 0.393 / 0.309
(Table 22, p. 41). «target perturbation measurements remain the most direct source of target response information.
DepMap supports cross-context and double-unseen extrapolation, but it does not fully replace target response
observations» (p. 41). The bound is not strict: target-only 30 % MLP in HepG2 (0.402) exceeds within-line MLP (0.372),
and the paper says «strongest reference for most cell lines» (p. 41).

### 25. For drugs, measured responses in other lines beat static chemical descriptors (pp. 3–4; Fig. 3d–e; pp. 17–18)

On Tahoe-100M (5 lines, 268 compounds × 3 doses; 3,662 of 4,020 conditions with ≥ 20 cells, p. 17), low-rank
source-response Ridge reaches mean r = 0.86 against 0.11 for Morgan fingerprints (p. 3). This persists under
leave-one-MoA-out, survives template centring, and collapses under compound-label permutation (100 permutations,
p. 18; Fig. 3e). The source-response representation averages the same compound-dose in the other lines (p. 17), so
this is few-shot ("one-shot") transfer, not zero-shot: «the compound was excluded from model training across all cell
lines» (p. 3), but its response in source lines is an input.

### 26. Pathway-level recovery separates models that aggregate scores conflate (p. 4; Fig. 3f, p. 7)

Against DepMap Ridge, the MLP gains most on signaling, inflammatory and stress programs. MORPH and STATE are worse
across all pathways (Fig. 3f).

### 27. Within-line models mostly reproduce the template (p. 4; Supp. Note 16, Table 23, p. 42)

Template-coefficient recovery is r = 0.62–0.75, residual recovery 0.20–0.26. Ridge puts 55 / 94 / 82 / 72 % of its
predicted variance on the template direction in K562 / RPE1 / HepG2 / Jurkat (p. 4).

### 28. Aggregate metrics are not enough; evaluate components and separate zero-shot from target-calibrated settings (pp. 1, 9)

«Performance before and after removal of shared response structure should be considered separately, and recovery of
perturbation-specific and interaction components should not be inferred from full-response metrics alone. Cross-context
studies should also distinguish zero-shot transfer from settings that use measured responses in the target context,
even when only a limited calibration set is available» (p. 9).

### 29. The decomposition is diagnostic, on four essential-enriched pseudobulk screens (p. 9)

«Our primary analyses use pseudobulk transcriptional responses from four CRISPRi screens enriched for essential-gene
perturbations … The additive decomposition is intended as a diagnostic description of response variation rather than a
complete mechanistic model» (p. 9). Regulatory state, lineage, knockdown efficiency, temporal dynamics or single-cell
distributions «may contain context information not represented in baseline expression» (p. 9).

## Figures and tables

- Fig. 1 (p. 5): decomposition schematic (a), task axes (b), **noise-corrected shares (c)**, template fraction (d),
  effective rank (e), template GSEA (f), identity axes vs template (g), essential vs non-essential (h),
  **template/residual captured by 5 and 100 control PCs (i)**.
- Fig. 3 (p. 7): **1 vs 3 source lines (a)**, double-unseen (b), one-to-one transfer and basis coverage (c), Tahoe (d–e),
  pathway recovery (f).
- Fig. 4 (p. 8): within-line attribution (a), **zero-shot β/γ attribution (b)**, target calibration (c–d).
- Fig. 8 (p. 32): 12 directed one-to-one transfers (a), STATE/MORPH with target info (b), **calibration saturation
  curves (c)**, basis coverage vs transfer (d), **reference baselines (e)**.
- Supp. Tables 1 (p. 34), 6–10 (pp. 36–37), 17–29 (pp. 40–44) hold the numbers quoted above.

## What this source does not address

- **No ceiling.** It never turns the variance components into an achievable-correlation bound or places a method as
  a fraction of one. Its only "upper reference" is empirical (within-line training, p. 41). The bound is the project's
  [[zero-shot-transfer-ceiling]]. Reliability-style ceilings (√ρ) are in @wang-2026-reliable-perturbations, which is
  not yet ingested.
- **No context-count curve.** The only manipulation is 1 vs 3 source lines, on 4 lines, for a context-free model
  (p. 15). There is no learned context correction, no break-even, and no sample-complexity argument in the number of
  contexts. Tahoe is used with 5 lines only (p. 17).
- **No transfer-residual analysis.** The control-PC projection (see [[control-subspace-projection]]) is applied to the
  observed response and its template-removed residual, not to truth − prediction. Prediction attribution correlates
  δ̂ with the true β and γ (p. 27) but never decomposes the prediction error.
- **No other-line or permuted-control null and no split-half control on the projection fractions.** The only null is a
  random basis (p. 36). Split-half noise correction is applied to the four-way partition only (pp. 21–22).
- **No correction method** (PC deflation, common-response removal) is proposed or compared. Template removal appears
  only as an analysis choice (Supp. Note 1).
- **No single-cell distributional evaluation.** Everything is pseudobulk (p. 9). There is no CRISPR-KO or non-cancer
  stem-cell destination: RPE1 is the only non-transformed line (p. 2), and the data are the
  [[replogle-nadig-crispri-collection]].

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 1, 2, 3, 7, 14, 15, 19, 20, 21, 22, 23, 24** (and 25 for Tahoe, 4 and 18 for per-destination
ceilings). The hypotheses are those of the plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4
at the top; Amendments A1–A3 at the end).

### What this focus needs beyond the Kernthesen

1. **Context count.** The genetic analyses use **4 cell lines** (K562, RPE1, HepG2, Jurkat; p. 10). The variance
   partition uses a balanced tensor of **4 × 638 perturbations × 2,000 HVGs** (p. 20), with HVGs taken once from the
   6,642 shared genes on pooled cells without a batch key (p. 11). Tahoe uses **5 lines** (p. 17). The only
   context-count contrast is 1 vs 3 sources. The one-source arm in Table 19 has 4 directed pairs, while Fig. 8a shows
   all 12, so the 1-vs-3 comparison mixes source identity with source count (p. 15).
2. **Baselines and nulls.** The headline metric is perturbed-reference Pearson (training-centroid subtraction, after
   Systema @vinastorne-2026-systema; p. 12). **In zero-shot the reference is the observed target centroid:**
   «Zero-shot perturbed-reference evaluation used the observed target perturbation centroid as an evaluation-only
   reference» (p. 15). Their cross-line numbers therefore never penalise a wrong template (b_c in project notation).
   The comparators are a DepMap Ridge baseline (Fig. 3f), within-line and target-only-30 % references (p. 41),
   STATE / MORPH / scGPT, and, for Tahoe, a mean-of-source-lines direct transfer plus a compound-label permutation
   null (pp. 18). The projection null is a random basis (p. 36). There is no context-free-vs-context-aware comparison,
   because every DepMap model is context-free.
3. **Raw responses vs residuals.** Every decomposition and projection is computed on the **observed** responses
   («All decompositions were computed from observed pseudobulk responses and used exclusively for post hoc
   analysis», p. 19). Models are attributed by corr(δ̂, β) and corr(δ̂, γ) (p. 27). No residual after transfer is
   analysed.
4. **Implied ceiling (a reader's derivation, not in the source).** (27.8 + 29.4) / 80.7 = 0.709 of reproducible
   energy is template + conserved, so √0.709 ≈ **0.84** (Supp. Table 1 numbers give the same 0.84). Caveats: (a) the
   template term includes the grand mean μ (p. 21), so this is an energy share, not Var(a) + Var(b); (b) it is
   tensor-level, not the per-perturbation mean correlation the paper reports; (c) the denominator excludes noise.
   Against the observed (noisy) truth it is √0.572 ≈ 0.76. (d) For the centroid-subtracted metric the paper actually
   uses, the analogue is conserved / (conserved + interaction) on the residual: 40.7 / 73.3 → **√0.555 ≈ 0.75**
   (ANOVA accounting) or 35.0 / 62.7 → √0.558 ≈ 0.75 (projective; Table 1, p. 34). The best observed zero-shot MLP
   (0.226–0.340, Table 18) sits far below either value, but those numbers are not noise-corrected.
5. **Overlap-note verification.** `input/notes/overlap-molina-zhang.md` was checked claim by claim against the PDF and
   is correct on every number and quote it gives. Five refinements: (i) its k = 5 residual numbers come from Table 10;
   Table 7 of the same note disagrees (Kernthese 7). (ii) "STATE/MORPH recover neither β nor γ" is slightly too
   strong: STATE has small positive corr(δ̂, β) of 0.043–0.120 (Table 28). (iii) The "upper reference" is not strict
   (Kernthese 24). (iv) The main-text "r = 0.25–0.39" for β alignment is the Table 25 R²(β) column (Kernthese 19).
   (v) The note read the HTML; equations, Tables 1–29 and Figs. 5–9 are now verified from the PDF.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partly shown (ingredients only)** | The noise-corrected components behind r_max, and the cross-half estimator that is the project's definition 4. Qualitatively, zero-shot models recover part of β and none of γ. No r_max, no attained fraction, no split of the remaining error into transferable vs non-transferable parts, and no frozen-pipeline or VCC-scorer evaluation. An r_max ≈ 0.84 is implied (focus item 4) and must be cited as "implied by" this paper. The per-destination template share (23–65 %) means r_max will differ by destination. | p. 2, Fig. 1c, pp. 21–22, Table 1 p. 34, Table 25 p. 43, Table 2 p. 35 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **not addressed** | Tahoe is used with 5 lines and few-shot source-response transfer only. There is no learned context correction, no gain-vs-C curve and no zero crossing. | pp. 17–18, Fig. 3d–e |
| **H2b** (theory predicts the break-even within ×2) | **not addressed** | No theory of context sample complexity and no d, σ²_ε/Var(i) or R² break-even. | — |
| **H2c** (genetic C*_min > C_public) | **not addressed (adjacent evidence)** | 1 → 3 sources improve a context-free map in all 4 targets (consistent with more sources sharpening β). γ is never recovered zero-shot, and the paper says it «cannot be predicted without information from the target perturbation context». The claim that γ is «learnable from currently available data» refers to target-calibrated training (30 % of target perturbations), so it does not contradict H2c. There is no C_public count and no C*_min. | p. 1, p. 3, p. 4, p. 15, Fig. 3a, Fig. 4d, Tables 18–19, 25 |
| **H3a** (transfer residual: destination top-8 control PCs > other lines' PCs and > permuted PCs, split-half) | **partly shown (response-level precedent only)** | It projects the **observed** template and template-removed residual onto the **same line's** control PCs: 34.5–39.3 % of the residual at 100 PCs and 7.1–13.4 % at 5 PCs (random basis about 0.3 %). There is no transfer residual, no other-line PCs, no permuted controls, no split-half on the fractions and no k = 8. Its §2.3 reading («the top PCs of the control cells do not capture variance beyond cell-line template», p. 3) must be reconciled with H3a, because a transfer residual still contains the mis-predicted template. | p. 2, Fig. 1i, pp. 35–37 Tables 6–10, p. 25 |
| **H3b** (M3 beats five PC-deflation / common-response corrections on PDS) | **not addressed** | No correction is compared. Template removal appears only as a decomposition choice (Supp. Note 1), and Systema-style centring only inside the metric. | p. 34, p. 12 |
| **H4** (adding same-line genome-wide data beats the largest modelling gain) | **partly shown (qualitative analogue)** | Target perturbation data outperform modelling: STATE 0.046 → 0.532 with 30 % target perturbations, and target-only 30 % beats zero-shot and recovers γ. Calibration saturates at about 100–150 perturbations, and pooled source data can mask target signal, which matters for H4's weight w. This is not on the VCC composite, not with H1GW, not a w-sweep, and not compared with a frozen pipeline's modelling gains. | Table 17 p. 40, Table 27 p. 44, p. 42, p. 43, Table 22 p. 41 |

### Direct quotes (for the draft)

> "after noise correction, we found that reproducible variance was distributed across the cell-line-specific component
> (the cell-line template; 27.8%), the conserved perturbation-specific component (29.4%), and the
> cell-line-by-perturbation interaction (23.5%), with the remaining 19.3% attributable to measurement noise" (p. 2)

> "The reproducible portion of γc,p,g represents target-specific interaction structure that is consistently detected
> across independent cell subsets and is therefore, in principle, available for model recovery." (p. 22)

> "Ridge and MLP received no explicit cell-line descriptor and consequently estimated a source-averaged response
> mapping for perturbations represented in multiple source contexts." (p. 15)

> "Zero-shot perturbed-reference evaluation used the observed target perturbation centroid as an evaluation-only
> reference. This centroid was used only after inference and was not available during model fitting or checkpoint
> selection." (p. 15)

> "Thus, control variation and perturbation-specific responses are not unrelated random directions, but they organize
> expression variation into different subspaces." (p. 36)

### Examples & illustrations (for later drafting)

- K562 ↔ RPE1 direct transfer: STATE 0.046 → 0.532 and MORPH −0.013 → 0.226 when 30 % of target perturbations are
  added (Table 17, p. 40). This is the clearest "data beats model" case for H4 framing.
- RPE1 vs K562: template 65 % vs 23 % of response energy (Table 2, p. 35). This illustrates why a ceiling must be per
  destination.
- K562 → RPE1 vs three sources → RPE1: MLP 0.211 vs 0.271 (Tables 18–19, p. 40), the 1-vs-3 context contrast.

### Boundary of this focus

This focus does not cover the DepMap-alignment and combinatorial results (Kernthesen 9–13) or the template biology
(5–6), apart from their relevance to ceilings. A reader looking for a context-count curve, a break-even or a
correction comparison will not find one here.

## Mentioned entities

- Concepts: [[perturbation-response-decomposition]], [[control-subspace-projection]], [[zero-shot-transfer-ceiling]]
- Datasets: [[replogle-nadig-crispri-collection]], [[tahoe-100m]]

## Connections

- Supports / introduces: [[perturbation-response-decomposition]] (p. 1), [[control-subspace-projection]] (pp. 35–36).
- Supplies the ingredients for: [[zero-shot-transfer-ceiling]] (inferred; no ceiling stated in the source).
- Cites (not yet ingested): @replogle-2022-genome-scale, @nadig-2025-transcriptome-wide, @adduri-2026-state,
  @zhang-2026-tahoe, @vinastorne-2026-systema, @ahlmanneltze-2025-linear-baselines.
- Earlier reading: `input/notes/overlap-molina-zhang.md` (verified; see focus item 5).
