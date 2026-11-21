---
title: "Ji et al. 2024/2025 — Prophet, a universal cellular-phenotype predictor"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: ji-2024-universal-prediction
depth: deep
tags: [perturbation-prediction, context-count, mean-baseline, zero-shot, unseen-cell-line, viability, ic50, transfer-learning, noise-ceiling, h2a]
relations:
  - target: context-count-learning-curve
    type: supports
    confidence: extracted
    because: "The only published curve of model performance against the share of training cell states, plotted against a context-free mean baseline: «evaluating Prophet with 80%, 50%, 30%, 20%, 10%, and 5% of the perturbations or cell lines included in the training dataset» (PDF p. 32; figs. S6, S7, PDF pp. 50–51)."
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "Defines the cross-context mean estimator exactly: «The mean baseline estimates that the effect of an perturbation in an unseen cell state is the average value of that perturbation in all previously seen states» (Fig. S1C legend, PDF p. 43)."
  - target: zero-shot-transfer-ceiling
    type: supports
    confidence: extracted
    because: "States a noise-induced ceiling in words: «Under such circumstances, the mean estimator represents an upper bound on predictive performance» (PDF p. 9), with the replicate correlations behind it in Tables S4–S5 (PDF pp. 59–60)."
  - target: dibaeinia-2026-context-scale
    type: supports
    confidence: inferred
    because: "Reader's comparison: both report performance rising with the number of training contexts, Prophet against a context-free baseline and Dibaeinia et al. without one; neither paper cites the other (Prophet's references, PDF pp. 21–30)."
  - target: prophet-universal-predictor
    type: builds-on
    confidence: extracted
    because: "Introduces the model and its three-axis experiment representation: «Prophet represents an experiment as the combination of three fundamental elements — the cellular state, the perturbations applied, and the phenotype measured» (PDF p. 2)."
---

# Scalable and universal prediction of cellular phenotypes enables in silico experiments (Prophet)

## Bibliographic Details

Ji, Yuge, Alejandro Tejada-Lapuerta, Niklas A. Schmacke, Zihe Zheng, Xinyue Zhang, Simrah Khan, Ina Rothenaigner,
Juliane Tschuck, Kamyar Hadian, Veit Hornung, and Fabian J. Theis. "Scalable and universal prediction of cellular
phenotypes enables in silico experiments." bioRxiv preprint, doi:10.1101/2024.08.12.607533. CC-BY 4.0.
Helmholtz Munich / TUM / LMU (Theis lab). Code: github.com/theislab/prophet (PDF p. 30).
Citekey: @ji-2024-universal-prediction.

**Version read: v3, the version posted 23 September 2025** (the running head of every page of
`<library>/pdf/ji-2024-universal-prediction.pdf` reads "this version posted September 23, 2025"). The BibTeX entry
carries `year = {2024}` for the v1 posting (12 August 2024). The pre-read note in
`input/bibliography/literaturguide.md` (lines ~13–36) was written from the **v2 full-text HTML**; every place where
v3 differs is flagged below under "v2 → v3 differences" and in the focus block. This page is written from v3 only.

**Pagination.** 71 PDF pages. The printed numbering restarts: main text PDF pp. 1–30 = printed 1–30; from PDF p. 31
(Methods) the printed numbers restart at 1, so PDF p. 43 prints "13", PDF p. 50 prints "20", PDF p. 57 prints "27".
**All anchors on this page are PDF pages**, written "PDF p. N" where confusion is possible. The BibTeX entry has no
`pages` field, so no printed range is being contradicted.

**Figure numbering shifted between v2 and v3.** The learning curves the literature guide calls "fig. S4 / S5" are
**fig. S6 (Prophet-individual) and fig. S7 (fine-tuned Prophet)** in v3. v3's fig. S4 is a bar chart of models per
benchmark setting, and v3's fig. S5 is "Mean baseline and evaluation schema". Below, "fig. SN" is always the **v3**
number.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | PDF pp. |
|---|---|---|
| ✓ F1 | Abstract; Introduction | 1–2 |
| ✓ F1 | Results: A unified framework for experimental cell biology (datasets, architecture); Fig. 1 legend | 3–6 |
| ✓ F1 | Results: Transfer learning improves in silico experimental prediction (baselines, mean estimator, scaling, SCORE/LINCS failures); Fig. 2 legend | 6–9 |
| — | Results: Prophet learns a universal landscape of biological experiments (activity fingerprints, MoA clusters, BIRC5, EGFR); Fig. 3 legend | 10–13 |
| — | Results: Virtual screening and in-vitro validation of melanoma-selective compounds; Fig. 4 legend | 13–15 |
| — | Results: Zero-shot in vitro validation under experimental constraints; Fig. 5 legend | 15–17 |
| — | Results: Prophet predicts in vivo phenotypes (zebrafish); Fig. 6 legend | 17–19 |
| ✓ F1 | Discussion, incl. Limitations | 19–20 |
| — | References | 21–30 |
| — | Acknowledgments, Funding, Author contributions, Competing interests, Data availability | 30 |
| ✓ F1 | Methods: External datasets; Baseline models; Data splits (incl. the 5–80 % grid); Performance evaluation (hit rate) | 31–33 |
| ✓ F1 | Methods: Overview of Prophet; Tokenization; Cell state representations; Treatment representations | 33–35 |
| — | Methods: Phenotype representations; In vivo phenotype prediction | 35–36 |
| ✓ F1 | Methods: Model training (incl. "capacity constrained by data availability"); Uncertainty | 36–37 |
| — | Methods: In silico screening compound selection | 37 |
| ✓ F1 | Supplementary Methods: Prophet-individual training; Comparison of Prophet-individual training; embeddings | 37–38 |
| — | Supplementary Methods: AUPRC, activity fingerprints, MoA annotation, plate scoring, wet-lab protocols, zebrafish KO space | 38–41 |
| ✓ F1 | Supplementary Notes: prompt syntax; usage of "phenotype"; input/output clarification (univariate scalar regression) | 41–42 |
| ✓ F1 | Fig. S1 (datasets; **mean-baseline definition**, C) | 43 |
| — | Fig. S2 (architecture) | 44–45 |
| ✓ F1 | Fig. S3 (hard partitions; **D = R² table per dataset and split**) | 45–46 |
| ✓ F1 | Fig. S4 (model comparison bars, Wilcoxon) | 47–48 |
| ✓ F1 | Fig. S5 (**A = R² variance vs sampled space**; B–E precision-recall, radar, Spearman table, HCC827) | 49 |
| ✓ F1 | **Fig. S6 (Prophet-individual vs mean as a function of the fraction of cell states / treatments in training)** | 50 |
| ✓ F1 | **Fig. S7 (fine-tuned Prophet vs mean, same two axes; C = R² SD vs # conditions)** | 51–52 |
| — | Figs. S8–S11 (clusters, spaces, melanoma, plate selection) | 52–55 |
| ✓ F1 | Fig. S12 (pretraining-configuration robustness, scaling) | 56 |
| ✓ F1 | Table S1 (dataset overview: readout type, # cell states, # perturbations) | 57 |
| — | Tables S2–S3 (hyperparameters) | 58–59 |
| ✓ F1 | **Tables S4–S5 (SCORE replicate correlations, cell-line and perturbation axis)** | 59–60 |
| — | Tables S6–S11 (MoA enrichment, compound annotations, cell lines, selected compounds, vendors) | 60–71 |

## Kernthesen

### 1. ⭐ An experiment is a triplet (cell state, perturbation, phenotype), and the model is a scalar regressor (PDF pp. 2–3, 33–34, 42)

«Prophet represents an experiment as the combination of three fundamental elements — the cellular state, the
perturbations applied, and the phenotype measured — and learns predictive mappings across this space» (p. 2). Formally
X_i = {X_env, X_treat, X_phe} is tokenised into a common space, a <CLS> token is prepended, and a feed-forward head
maps the transformer output to **one scalar** Y_i (pp. 33–34). Vector phenotypes are split: «Prophet is a univariate
regression model, so it operates on scalar targets. For high-dimensional data such as transcriptomics or Cell Painting
features … we decompose the task into multiple scalar regression problems—one per feature» (p. 42).

### 2. Architecture and training scale (PDF pp. 3, 6, 36; Table S3, p. 59)

8 transformer encoder blocks, 8 attention heads, FFN width 1,024, 512-dim experiment embedding, **20.1 M parameters**
(Table S3). Pretrained on **4.7 million experimental conditions, 26,000 perturbations, 1,550 cell states** across
seven screens, then fine-tuned per dataset for 20,000 steps at batch 256 (pp. 6, 36). Dataset-specific
("Prophet-individual") models are much smaller (2–4 blocks, 1–4 heads) to avoid overfitting (p. 36).

### 3. ⭐⭐ The readouts are viability, essentiality, L1000 fluorescence, morphology and cell-type proportions — there is no Perturb-seq and no transcriptome-wide readout (PDF p. 3; Table S1, p. 57)

Table S1 lists ten resources with their readout type: Horlbeck (CRISPR combination, viability, 1 cell state),
GDSC (drug singleton, ln(IC50), 682 cell states), GDSC2 (drug combination, ln(IC50), 110), JUMP/cpg-0000 (drug+CRISPR,
200 Cell Painting features, 2), **LINCS (drug singleton, "RNA (fluorescence)", 42 readouts, 57 cell states, 7,306
perturbations)**, PRISM (drug singleton, log2FC, 556), **SCORE (CRISPR singleton, CERES essentiality, 249 cell states,
15,573 perturbations)**, Saunders/zebrafish (cell-type proportions, 620), CTRP (drug singleton, IC50, 822).
*Inconsistency in the source:* the Table S1 legend says LINCS has «a total of 50 selected genes» while the table's
readout column says 42. Either way the transcriptomic arm is a few dozen L1000 genes, not a transcriptome.

### 4. ⭐ The two held-out axes are "unseen cell states" and "unseen perturbations"; the cell-state axis is the unseen-cell-line zero-shot task (PDF pp. 6, 31–32; fig. S1C, p. 43)

«For unseen cell states, the test dataset consists of experiments in entirely new cell states» (p. 31). Test = 20 % of
the perturbations or cell lines, validation = 10 % (p. 31). Evaluation is 5-fold cross-validation per seed, repeated
for 3 seeds (fig. S1C legend, p. 43). The axis is therefore the project's axis — an unseen cell line — but at the
granularity of a CCLE cell line, with a **scalar** phenotype per experiment.

### 5. ⭐⭐ EXACT definition of the mean estimator / mean baseline (fig. S1C legend, PDF p. 43; text p. 6)

> «The mean baseline estimates that the effect of an perturbation in an unseen cell state is the average value of that
> perturbation in all previously seen states. Conversely, it estimates the effect of an unseen perturbation on a cell
> state is the average of all seen perturbations in that same cell state. Naturally, it changes depending on the test
> set used.» (fig. S1C legend, p. 43; "an perturbation" sic)

The figure panel itself labels the two cases: **"Mean baseline = assumes drug behaves the same across all cell lines"**
(cell-state axis) and "Mean baseline = assumes cell state response is not drug-specific" (perturbation axis)
(fig. S1C, p. 43, read from the figure). In the main text it is «the mean estimator, which summarizes the global effect
of perturbations across prior experiments (fig. S1C, Methods) and is used as a difficult baseline in drug response
prediction» (p. 6).

**On the unseen-cell-state axis this averages over the seen cell lines for the same perturbation.** In the project's
notation that is exactly the context-free model μ + a_t: one number per perturbation, no cell-state input. It is
therefore the same *object* as our context-free arm — but see Kernthese 6 for why it is not the same *comparison*.

### 6. ⭐⭐ The mean estimator is a different estimator class, not an ablation of Prophet (PDF pp. 6, 31, 43)

Nothing in the Methods matches the mean baseline to Prophet in capacity, budget, early stopping or seed: it is a
closed-form average computed on the training split, compared against a trained transformer (pp. 6, 31). The baseline
family (linear regression, MLP, Random Forest, mean) is a **model comparison**, not a context-pathway ablation. The
project's A1 amendment exists precisely because that difference confounds the comparison
(`input/ideas/zero-shot-limits-context-complexity-plan.md`, Amendment A1).

### 7. ⭐⭐ The context-count manipulation: 5 / 10 / 20 / 30 / 50 / 80 % of cell lines or perturbations in training (PDF p. 32)

> «Additionally, we also assess the model's performance as it is exposed to progressively more perturbations or cell
> lines (fig. S4,5), evaluating Prophet with 80%, 50%, 30%, 20%, 10%, and 5% of the perturbations or cell lines
> included in the training dataset.» (p. 32)

*The in-text pointer "(fig. S4,5)" was not updated when the figures were renumbered in v3: the curves are figs. S6 and
S7 (pp. 50–51), not S4/S5.* Five datasets have enough cell states for the sweep: CTRP, GDSC, GDSC2, LINCS, SCORE
(p. 9). The curve is over a **share**, and the panel titles give the totals, so absolute counts follow (Kernthese 8).

### 8. ⭐⭐ Fig. S6A (PDF p. 50) — Prophet-individual vs the mean baseline along the **cell-line** axis: crossings read off the figure

Legend: «Fig. S6. Prophet-individual models' performance as a function of the number of conditions seen during
training. The more the experimental space is mapped during training, the better Prophet's performance. However, the
concrete pattern is dataset dependent as different datasets capture different amounts of variation. A) Performance of
Prophet-individual models along the cell line axis. B) … along the perturbation axis.» (p. 50)

Axes: x = "Fraction of cell states seeing in training" (sic), ticks at 0.25/0.50/0.75, points at 0.05/0.10/0.20/0.30/
0.50/0.80; y = **R²**. Two lines (dark = Prophet, tan = mean) with shaded bands and open circles for the 3 seeds × 5
folds. **All readings below are visual readings of the rendered figure page, not tabulated numbers** (the paper gives
no table for these curves):

| Panel (n cell states) | mean baseline | Prophet-individual | Crossing (share → approx. training lines) |
|---|---|---|---|
| CTRP (n = 822) | flat ≈ 0.45–0.47 | −0.3 at 0.05 → 0.36 (0.10) → 0.40 (0.20) → 0.42 (0.30) → 0.43 (0.50) → ≈ 0.45 (0.80) | **no crossing**; only converges at 0.80 (≈ 658 lines) |
| GDSC (n = 682) | flat ≈ 0.705–0.71 | 0.655 (0.05) → 0.68 (0.10) → 0.695 (0.20) → 0.725 (0.30) → 0.735 (0.50) → 0.745 (0.80) | **≈ 0.25 → ≈ 170 lines** |
| GDSCcomb (n = 110) | 0.05 (0.05) → 0.25 (0.30) → 0.27 (0.50) → 0.24 (0.80) | −0.4 (0.05) → 0.05 (0.20) → 0.12 (0.30) → 0.17 (0.50) → 0.20 (0.80) | **no crossing** up to 0.80 (≈ 88 lines); still below |
| LINCS (n = 57) | 0.40 (0.05) → 0.50 (0.10) → 0.57 (0.20) → 0.62 (0.30) → 0.645 (0.50) → 0.67 (0.80) | 0.30 (0.05) → 0.50 (0.10) → 0.59 (0.20) → 0.645 (0.30) → 0.665 (0.50) → 0.71 (0.80) | **≈ 0.10 → ≈ 6 lines**; a clear margin only from ≈ 0.3 (≈ 17 lines) |
| SCORE (n = 249) | 0.805 (0.05) → 0.835 (0.20) → 0.842 (0.50) → 0.842 (0.80) | 0.80 (0.05) → 0.811 (0.20) → 0.812 (0.50) → 0.820 (0.80) | **no crossing**; the mean wins at every fraction incl. 0.80 (≈ 199 lines) |

So of five datasets, **two cross (GDSC, LINCS) and three never do** within the measured range — the "varies per
dataset" statement of v2, now only visible in the figure.

### 9. ⭐ Fig. S6B (PDF p. 50) — the same curves along the **perturbation** axis

x = "Fraction of treatments seeing in training". CTRP (n = 545 perturbations): mean flat ≈ 0.10, Prophet −0.12 (0.05)
→ crosses at ≈ 0.5 → 0.19 (0.80), with a large dip at 0.30. GDSC (n = 224): mean ≈ 0.05–0.10, Prophet crosses at
≈ 0.10–0.15 and plateaus ≈ 0.25. GDSCcomb (n = 61): crosses at ≈ 0.5–0.6, reaching ≈ 0.33 vs 0.22. LINCS (n = 7,306):
mean ≈ 0.895 flat, Prophet 0.79 (0.05) → 0.88 → converges at ≈ 0.895 at 0.80, **never above**. SCORE (n = 15,573):
mean ≡ 0.00 at every fraction, Prophet 0.19 (0.05) → 0.30 (0.80), **above everywhere**. (Visual readings.)

### 10. ⭐⭐ Fig. S7A (PDF pp. 51–52) — fine-tuned (pretrained) Prophet along the cell-line axis: pretraining moves the crossing left but does not rescue SCORE

Legend: «Fig. S7. Prophet models' performance as a function of the number of conditions seen during fine tuning. It
can be seen that performance is more robust when fine tuning. As a remark, all models here are fine tuned for just
20,000 steps, independently of the size of the dataset. A) Performance of Prophet along the cell line axis. B) …
along the perturbation axis. C)» (pp. 51–52; the C) sentence is unfinished in the source).

Visual readings: **CTRP** mean ≈ 0.46 flat, Prophet 0.29 (0.05) → 0.43 (0.20) → ≈ 0.45 (0.80) — still **no crossing**.
**GDSC** Prophet ≈ 0.735 already at 0.05 against a mean of ≈ 0.705 — **above at every fraction, i.e. the crossing has
moved below 0.05 (< ≈ 34 lines)**. **GDSC² (n = 110)** Prophet ≈ 0.05–0.17 at 0.05, 0.35 at 0.20 against a mean of
≈ 0.22 — crossing at **≈ 0.10 (≈ 11 lines)**, where the untuned model never crossed. **LINCS** Prophet ≈ 0.60 vs mean
0.43 at 0.05 — above throughout (crossing < 0.05, < 3 lines). **SCORE** Prophet 0.811 vs mean 0.805 at 0.05, then the
mean overtakes at ≈ 0.10–0.15 and stays above (0.843 vs 0.818 at 0.80), with the Prophet curve *falling* after 0.50.

### 11. ⭐ Fig. S7B — fine-tuned, perturbation axis

CTRP crosses at ≈ 0.4; GDSC is above throughout (≈ 0.33 vs 0.08); GDSC² jumps above at ≈ 0.20 (≈ 0.45 vs 0.10);
**LINCS stays below the mean at every fraction** (0.885 vs 0.895); SCORE is above everywhere (mean ≡ 0). (Visual.)

### 12. ⭐⭐ What v3 actually says in the Results about scaling — the v2 "amount of data required" sentence is gone (PDF p. 9)

> «Finally, we assessed how Prophet's predictive uncertainty and performance scale with the number of conditions
> observed during training (Methods). We fine-tuned models on increasingly more training data, using datasets with
> enough cell states (CTRP, GDSC, GDSC2, LINCS, SCORE). This revealed that confidence in predictions increased as more
> perturbations and cell states were incorporated (fig. S5A). Performance similarly improved with additional data
> (fig. S6, S7).» (p. 9)

There is **no** sentence in v3 reading "The amount of data required to outperform the mean estimation varies per
dataset". The claim is now carried only by the figures. See "v2 → v3 differences" below.

### 13. ⭐⭐ The failure cases and the noise ceiling (PDF p. 9)

> «Despite its strong performance, Prophet encounters limitations in certain scenarios. For example, predicting the
> effects of observed perturbations in untested cell lines within SCORE, or forecasting the impact of untested small
> molecules in known cell lines within LINCS, can be challenging. Further analysis revealed that the correlation
> between experimental replicates is a key determinant of predictability (Table S4, Methods). In some cases,
> replicates of the same experimental condition are less similar to each other than to replicates of other random
> conditions, reflecting high experimental noise and batch effects. Under such circumstances, the mean estimator
> represents an upper bound on predictive performance.» (p. 9)

This is the paper's whole theory of when the mean baseline cannot be beaten: **replicate reliability, not the number
of contexts**. Note that v3 says "encounters limitations", not v2's "encounters challenges in outperforming the mean
estimator", and adds the LINCS unseen-drug case.

### 14. ⭐ Table S4 (PDF p. 59) — why the SCORE **cell-line** axis fails: the cell-line signal does not survive aggregation

Pearson correlations on SCORE along the cell-line axis: between random replicates 0.481 (nothing aggregated) and
**0.674** (perturbational replicates aggregated); within a cell line 0.615 (nothing aggregated) and 0.820
(perturbational replicates aggregated); across cell lines 0.625 / 0.755. The legend draws the conclusion:
«the Pearson correlation between random replicates once the perturbational replicates are aggregated (0.674) is
bigger than the Pearson correlation between replicates of a same cell line (0.615). This means that the cell
line-specific information gets lost once the replicates are aggregated and, hence, one can not distinguish a random
replicate from a real one.» (Table S4 legend, p. 59)

### 15. ⭐ Table S5 (PDF p. 60) — and why the **perturbation** axis works on the same data

Within perturbation 0.125 (nothing aggregated) / 0.208 (cell-line replicates aggregated) against random replicates
0.029 / 0.046. «This means that the perturbation-specific information does not get lost once the replicates are
aggregated and, hence, one can distinguish a random replicate from a real one. Therefore, the mean estimator can be
outperformed since the experiment retains biological signals.» (Table S5 legend, p. 60). Tables S4/S5 together are a
reliability argument of exactly the shape the project's H1 needs — but computed on CERES viability scores, not on
expression.

### 16. ⭐ Fig. S3D (PDF pp. 45–46) — the per-dataset R² table, Prophet vs the mean column

Read from the figure (columns: Prophet | Prophet-individual | mean | RF | MLP | LR): LINCS-unseen-cl 0.68 | 0.71 |
**0.67**; PRISM-unseen-cl 0.59 | 0.60 | **0.62**; CTRP-unseen-cl 0.45 | 0.46 | **0.47**; GDSC-unseen-cl 0.75 | 0.74 |
**0.71**; GDSCcomb-unseen-cl 0.35 | 0.20 | **0.25**; SCORE-unseen-cl 0.82 | 0.82 | **0.84**. On the unseen-cell-line
axis Prophet beats the mean on 3 of 6 datasets at the full training fraction and loses on CTRP, PRISM and SCORE.
On the unseen-perturbation rows Prophet is at or above the mean nearly everywhere (e.g. SCORE 0.27 vs −0.00,
GDSC 0.33 vs 0.09, JUMP 0.90 vs 0.89).

### 17. ⭐ Fig. S5D (PDF p. 49) — Spearman "with mean subtracted" is ≈ 0 on most unseen-cell-line rows

The heat map reports Spearman correlation **with the mean baseline's value subtracted**: LINCS-cell-line 0.03,
PRISM-cell-line 0.00, CTRP-cell-line 0.00, GDSC-cell-line 0.04 and 0.12, SCORE-cell-line 0.00 — against clearly
positive perturbation-axis rows (CTRP 0.21, GDSC 0.23 and 0.21, SCORE 0.35, ShifrutMarson 0.16). Read from the
figure. The cell-line axis is where Prophet is closest to the mean estimator.

### 18. The headline metric is a hit rate, not a correlation (PDF pp. 32–33)

«Of these, hit rate is the most relevant and interpretable»: the overlap of the true top-5 % and bottom-5 % extrema
with the predicted ones, averaged, **minus the chance level 0.05** (p. 32). R² and Spearman are secondary (p. 32).
Hit rate rewards extreme-response ranking, which is a drug-screening objective, not a transcriptome-reconstruction
objective.

### 19. Transfer-learning gains over the single-dataset model (PDF pp. 6, 8–9)

Prophet-individual already beats the simple baselines: «improving hit rates by an average of 63% over the best simpler
architecture and capturing up to 33% of all hits» (p. 6). Joint GDSC2+PRISM training improves unseen perturbations and
unseen cell states by 3 % and 28 % (p. 6). Pretraining + fine-tuning improves GDSC2 R² «by 70% for unseen
perturbations and 30% for unseen cell states compared to Prophet-individual» (p. 8). Fine-tuned Prophet «consistently
matched or outperformed all baseline models … across every benchmark» (p. 6).

### 20. ⭐ Unseen cell states are the harder axis (PDF p. 6)

«Hit rates for unseen perturbations were generally higher than for unseen cell states (Fig. 2E, fig. S5C), reflecting
the relative difficulty of extrapolating to new cellular contexts» (p. 6). The same ordering is visible in Fig. S3A–C
(pp. 45–46) and in the figures read above.

### 21. Robustness over 64 stratified tasks (PDF p. 9)

«Across 64 such stratified tasks, Prophet outperformed Prophet-individual and baseline models in 69% of cases, with
differences in the remaining tasks being minimal (≤ 0.015 hit ratio, fig. S3)» (p. 9). The stratifications include
unseen tissues, unseen scaffolds, unseen drug/gene/cell-line clusters, and both-unseen (pp. 32).

### 22. ⭐ The cell-state representation is bulk CCLE expression, PCA-reduced to 300 dimensions (PDF pp. 35, 42)

«we used bulk RNA-seq from the Cancer Cell Line Encyclopedia (CCLE) … which covers 1406 cell lines. To reduce training
time, we reduced the dimensionality of both sources using PCA to 300 components» (p. 35). Consequence stated in the
notes: «For unseen cell lines, predictions are limited to those for which transcriptomic profiles are available in the
CCLE dataset» (p. 42). The context descriptor is therefore a 300-dim basal-expression vector, the same family of
descriptor the project's d is estimated from — but its dimension is a preprocessing choice, never analysed.

### 23. ⭐ The only quantitative data-scaling law in the paper is about **variance**, not bias (fig. S5A p. 49; fig. S7C p. 51)

«A) The variance in the R2 score obtained by Prophet is reduced as more parts of the experimental space - more cell
states and more perturbations - are sampled during training» (fig. S5 legend, p. 49). Fig. S7C plots "R2 standard
deviation" against "# conditions in training" on log–log axes for both axes (cell states, interventions), a declining
cloud from ≈ 10⁰ at ~3 conditions to ≈ 10⁻³ at ~10⁴ (read from the figure). No functional form is fitted and no
break-even is derived.

### 24. ⭐ Capacity is data-limited, not architecture-limited (PDF pp. 36–37; fig. S12, p. 56)

«We also observed that scaling beyond ~20M parameters in the pre-trained setting offered no performance benefit
(fig. S12). This, taking into account the previously mentioned overfitting of the large model in the dataset-specific
case, suggests that model capacity is currently constrained by data availability rather than architectural
limitations» (pp. 36–37).

### 25. ⭐ Stated limitations — and the absence of any noise model (PDF p. 19)

«Limitations of the current approach include the simplified assumption of the three-axis experimental space (Fig. 1).
Other elements, such as temporal dynamics, dosage, or technical variations, are also part of each experiment …
**Prophet is not designed to correct for intra-study experimental noise, assuming that data preprocessing and quality
control have already been applied.** Additionally, it does not currently implement formal Bayesian uncertainty
estimation; ensemble-based approximations of predictive uncertainty provide only a rough estimate, with correlated
members» (p. 19). Uncertainty = the variance of six models differing only in the validation seed (p. 37).

### 26. ⭐⭐ There is no theory anywhere in the paper (PDF pp. 1–71)

No proposition, no lemma, no identifiability argument, no descriptor-dimension term, no noise term, no formula for how
many contexts are needed. The words "identifiab*", "break-even", "sample complexity" and "how many cell lines" do not
occur. The only predictions about data amount are the empirical curves (Kernthesen 8–11) and the variance law
(Kernthese 23).

### 27. Applications not in this focus (PDF pp. 10–19)

Activity-based fingerprints over 979 cell lines × 8,620 compounds (p. 11); MoA clusters with Fisher-exact enrichment
(Table S6, p. 60); a 1.9 M-compound ChEMBL in silico screen in 1,406 cell lines (2.7 billion experiments, 0.12 % seen
in training) yielding 123 candidates, 5 tested, 4 validated in A375/SK-MEL-2 vs HeLa (pp. 14–15); a pre-built-plate
zero-shot screen with median Spearman 0.53 against 0.07 for PRISM measurements of the same compounds (p. 16); and
zebrafish cell-type proportions from 23 knockouts trained from scratch (pp. 17–19).

## Figures and tables

- **Fig. S1C (p. 43)** — the 5-fold CV schematic and **the mean-baseline definition** (both axes), plus the in-figure
  statement "Mean baseline = assumes drug behaves the same across all cell lines".
- **Fig. S3D (pp. 45–46)** — R² per dataset × split for Prophet / Prophet-individual / mean / RF / MLP / LR.
- **Fig. S5A (p. 49)** — R² variance shrinks as more of the experimental space is sampled. **Fig. S5D** — Spearman
  with the mean subtracted.
- **Fig. S6 (p. 50)** — ⭐ Prophet-individual vs mean as a function of the fraction of cell states (A) and treatments
  (B) seen in training. The core figure for this focus.
- **Fig. S7 (pp. 51–52)** — ⭐ the same for the pretrained, fine-tuned Prophet; **C** = R² SD vs # conditions.
- **Fig. S12 (p. 56)** — no benefit beyond ≈ 20 M parameters.
- **Table S1 (p. 57)** — readout type, # cell states, # perturbations, sparsity per dataset.
- **Tables S4–S5 (pp. 59–60)** — ⭐ SCORE replicate correlations along the cell-line and perturbation axes.

## What this source does not address

- **No break-even count and no theory.** The curves are plotted against a *share* of cell states, never converted to a
  number of contexts, never fitted, and never accompanied by a formula, an identifiability condition or a noise term
  (Kernthese 26). Our C* formula has no counterpart here, neither to confirm nor to contradict.
- **No matched context-free ablation.** The comparator is a closed-form mean, not the same network with its context
  pathway removed under an identical budget (Kernthese 6).
- **No transcriptome.** Every readout is a scalar: IC50, viability, CERES, 42–50 L1000 genes, Cell Painting features,
  cell-type proportions (Kernthese 3). There is no CRISPRi Perturb-seq, no Tahoe, no gene-level log-fold-change vector,
  and no DE or perturbation-discrimination metric.
- **No noise-corrected ceiling.** Replicate correlations exist for SCORE only (Tables S4–S5) and are used verbally
  («the mean estimator represents an upper bound»), never as a quantitative bound on attainable correlation, and never
  as a denominator for a fraction-attained statement.
- **No decomposition and no residual analysis.** Responses are never split into transferable and interaction parts,
  predictions are never attributed to components, and the prediction error is never projected onto control PCs.
- **No correction method.** Nothing comparable to PC deflation or common-response removal is proposed or compared.
- **No scorer-level composite.** The metrics are hit rate, R² and Spearman per dataset; there is no VCC-style
  composite and no data-vs-modelling budget comparison.

## v2 → v3 differences (verified against the PDF; the literature guide quotes v2)

`input/bibliography/literaturguide.md` (lines ~13–36) was written from the v2 full-text HTML. Checked line by line
against v3:

| Guide's v2 quote / claim | Status in v3 (PDF read 2026-11-05) |
|---|---|
| «we trained Prophet-individual models to predict the effect of either unseen treatments or unseen cell states, including 5-80% of the data along the respective axis. The amount of data required to outperform the mean estimation varies per dataset … (fig. S4, S5)» | **Deleted.** No such sentence anywhere in v3. The 5–80 % grid survives only in Methods as «evaluating Prophet with 80%, 50%, 30%, 20%, 10%, and 5% of the perturbations or cell lines included in the training dataset» (p. 32), and the Results now read «Performance similarly improved with additional data (fig. S6, S7)» (p. 9). **The headline sentence of the STOP-FLAG can no longer be quoted from the version of record.** |
| «Fig. S4A legend: Performance of Prophet along the cell line axis.» | **Renumbered.** That exact wording is now **fig. S7A** (fine-tuned Prophet, p. 51); the Prophet-individual curve is **fig. S6A**, legended «Performance of Prophet-individual models along the cell line axis» (p. 50). v3's fig. S4 is a bar chart of models per benchmark setting (pp. 47–48). |
| «Fig. S5A defines the baseline: The mean baseline estimates that the effect of an intervention in an unseen cell state is the average value of that intervention in all previously seen states.» | **Moved and reworded.** The definition is in the **fig. S1C legend (p. 43)**, with "intervention" replaced throughout by "perturbation" («the effect of an perturbation in an unseen cell state», sic). v3's fig. S5A is a different panel (R² variance vs sampled space, p. 49). The substance is unchanged, and the guide's reading of it as our context-free model μ + a_t is correct. |
| «Prophet encounters challenges in outperforming the mean estimator in certain scenarios, such as predicting the effects of observed treatments in unobserved cell lines within SCORE» | **Reworded and broadened.** v3: «Prophet encounters limitations in certain scenarios. For example, predicting the effects of observed perturbations in untested cell lines within SCORE, or forecasting the impact of untested small molecules in known cell lines within LINCS, can be challenging» (p. 9). The explicit "in outperforming the mean estimator" clause is gone; a second failure case (LINCS, unseen drugs) is added. |
| «They attribute this to low replicate correlation (Table S4), not to the number of contexts.» | **Confirmed verbatim in v3** (p. 9; Tables S4–S5, pp. 59–60), including «the mean estimator represents an upper bound on predictive performance». |
| «It is a share-of-lines curve (5–80 %), without a reported break-even count C*.» | **Confirmed.** No count anywhere. The absolute numbers in Kernthese 8 are this page's conversion using the panel-title totals. |
| «There is no theory: nothing on identifiability rank ≤ C−1, descriptor dimension d, or a noise term.» | **Confirmed** (Kernthese 26). |
| «The data are mostly viability/IC50 (GDSC, CTRP, PRISM), CRISPR viability (SCORE) and L1000.» | **Confirmed** (Table S1, p. 57), with the detail that LINCS is 42 (legend: 50) selected genes measured as fluorescence, not a transcriptome. |
| «The figure image could not be retrieved … the per-dataset crossing points are **unread**. The user must look at fig. S4A manually.» | **Now read** — see Kernthesen 8–11. This closes the open action item of the STOP-FLAG; no manual inspection by the user is required for the crossings, only a check of the visual readings. |

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 22, 23, 24, 25, 26.** The hypotheses are
those of the plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 near the top; Amendments A1–A3
at the end; A3 names this paper as the prior empirical evidence for the H2a curve).

### What this focus needs beyond the Kernthesen

1. **Is Prophet's axis our axis?** Yes for the manipulated variable, no for the readout. The "unseen cell state" split
   holds out whole CCLE cell lines and predicts a perturbation's effect in them (pp. 31–32), which is our zero-shot
   destination setting. But the predicted quantity is one scalar per (line, perturbation, phenotype) — IC50, CERES,
   one L1000 gene, one Cell Painting feature — never a gene-expression vector (Kernthese 3, p. 42). Our H2a quantity
   is a transcriptome-wide held-out residual MSE on Tahoe. The two curves measure different things on the same axis.
2. **Is the mean baseline context-free in our sense?** Yes as an estimator (average over the seen cell lines for the
   same perturbation, fig. S1C legend, p. 43 — in our notation μ + a_t, with the in-figure gloss "assumes drug behaves
   the same across all cell lines"), **no as a comparison**: it is a closed-form average set against a trained
   transformer, not the identical network with a constant descriptor under the same budget, early stopping, selection
   and seed. That is precisely the confound Amendment A1 was written to remove, so Prophet's crossings are
   *model-class* crossings, not *context-pathway* break-evens. Our C*_emp and theirs are not the same quantity.
3. **Absolute counts at the crossings.** The paper plots shares; the panel titles give the totals (Table S1, p. 57
   agrees), so, as visual readings of fig. S6A / S7A (pp. 50–51): Prophet-individual crosses at ≈ 0.25 × 682 ≈ **170**
   GDSC lines and ≈ 0.10 × 57 ≈ **6** LINCS lines, and never crosses on CTRP (≤ 658 lines), GDSCcomb (≤ 88) or SCORE
   (≤ 199). Fine-tuned Prophet crosses below 0.05 on GDSC (< ≈ 34) and LINCS (< ≈ 3), at ≈ 0.10 on GDSC² (≈ 11), never
   on CTRP, and **never on SCORE, where the gap widens with more lines**. The spread 6 → > 658 across five datasets of
   one paper is itself the strongest published argument that a break-even is dataset-specific and must be measured.
4. **The genetic datapoint.** SCORE is the one CRISPR screen on the cell-line axis: 249 lines, 15,573 single
   knockouts, CERES essentiality (Table S1, p. 57). Prophet is below the mean at **every** training fraction including
   80 % (≈ 199 lines), and the authors trace it to replicate reliability — random replicates correlate 0.674 after
   aggregation against 0.615 within a cell line (Table S4, p. 59). This is adjacent to H2c (genetic screens sitting
   below break-even) but it is confounded: with 199 training lines the count is far above any plausible C*, so what
   fails there is the *signal*, not the *context budget*. Cited as H2c evidence, it must be cited with that caveat.
5. **Does it already contain our noise lower bound?** No. «Under such circumstances, the mean estimator represents an
   upper bound on predictive performance» (p. 9) is a verbal noise ceiling and the closest published statement to our
   H1 framing, but there is no r_max, no reproducible-variance denominator, no split-half estimator and no
   fraction-attained. Tables S4–S5 are the ingredients of a reliability ceiling on one dataset, never assembled.
6. **What Prophet does *not* do that the guide's near-miss list does.** Unlike [[dibaeinia-2026-context-scale]], whose
   context curve has no context-free baseline at all, Prophet's curve has one on both axes — that is what makes it
   the prior evidence for H2a rather than a near-miss. Unlike [[molina-2026-response-decomposition]] (1 vs 3 source
   lines), it varies the count over six levels. Neither of those papers, and not this one, de-confounds the context
   pathway from the model class.

### Overlap with H1–H4

| Hyp. | Verdict | What Prophet already shows exactly | Anchors (PDF pages) |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partial (verbal ceiling only)** | States that measurement noise can make the context-free mean unbeatable — «the mean estimator represents an upper bound on predictive performance» — and gives the replicate correlations behind that for one dataset (SCORE: random 0.481/0.674 vs within-line 0.615/0.820). No r_max, no reproducible-variance denominator, no split-half correction, no attained fraction, no error decomposition, and the metric is scalar R²/hit rate rather than noise-corrected correlation on a transcriptome. | p. 9; Tables S4–S5, pp. 59–60; fig. S3D, pp. 45–46 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **overlap — this is the prior empirical evidence** | Six training fractions (5–80 %) of cell states on five datasets, each plotted against the context-free mean estimator, with dataset-dependent crossings: GDSC ≈ 0.25 (≈ 170 lines) and LINCS ≈ 0.10 (≈ 6 lines) for Prophet-individual; no crossing on CTRP, GDSCcomb or SCORE up to 80 %; pretraining moves GDSC and LINCS below 0.05 and GDSC² to ≈ 0.10. **Still new for us:** a matched context-free ablation (identical data, budget, early stopping, selection, seed — Amendment A1) instead of a different estimator class; a transcriptome-wide residual-MSE gain on Tahoe instead of scalar R²; an actual break-even *count* with an interval, C*_emp = 37.4 [33.3, > 40]; and a negative gain demonstrated at a specific C = 12. | p. 32; figs. S6 (p. 50), S7 (pp. 51–52); Table S1, p. 57 |
| **H2b** (theory predicts the break-even within ×2) | **none** | No theory of any kind: no identifiability condition, no descriptor dimension d, no σ²_res, no R², no formula, no fitted curve shape (Kernthese 26). The 300-dim CCLE PCA descriptor is a preprocessing choice that is never analysed (p. 35). Prophet is in fact *negative* evidence for any single formula, since its five crossings span ≈ 6 to > 658 training lines under one protocol. Our rejected linear-Gaussian C* has no counterpart to compare against here. | pp. 35, 42; figs. S6–S7 |
| **H2c** (genetic C*_min > C_public) | **partial, adjacent and confounded** | SCORE, the one genetic cell-line-axis screen (249 lines, 15,573 KOs, CERES), is the clearest case where the mean estimator is never beaten — at 80 % of lines the mean still leads 0.842 vs 0.820, and fig. S3D gives 0.84 vs 0.82 at full training. But the cause the authors give is replicate noise, not context count, and with ≈ 199 training lines the count cannot be the binding constraint. Also: viability, not transcriptome; no C_public and no C*_min anywhere. | p. 9; fig. S6A, p. 50; fig. S3D, pp. 45–46; Tables S1 p. 57, S4 p. 59 |
| **H3a** (transfer residual on destination top-8 control PCs) | **none** | No residual is ever formed or analysed; no PCA of controls; no permutation null; no split-half. The CCLE profile enters only as a 300-dim input descriptor (p. 35). | p. 35 |
| **H3b** (M3 beats five PC-deflation / common-response corrections) | **none** | No correction method is proposed or compared. The paper explicitly declines the noise question: «Prophet is not designed to correct for intra-study experimental noise» (p. 19). | p. 19 |
| **H4** (same-line data beats the largest modelling gain) | **partial (qualitative, same direction)** | Says in so many words that data, not architecture, is the binding constraint — «model capacity is currently constrained by data availability rather than architectural limitations» — and shows that more training conditions along either axis raise performance while scaling past ≈ 20 M parameters does not. But there is no composite scorer, no w-sweep, no same-line-vs-modelling budget comparison, and the "data" added is other contexts rather than destination measurements. | pp. 36–37; fig. S12, p. 56; figs. S6–S7 |

### Direct quotes (for the draft)

> "The mean baseline estimates that the effect of an perturbation in an unseen cell state is the average value of that
> perturbation in all previously seen states. Conversely, it estimates the effect of an unseen perturbation on a cell
> state is the average of all seen perturbations in that same cell state. Naturally, it changes depending on the test
> set used." (fig. S1C legend, PDF p. 43)

> "Additionally, we also assess the model's performance as it is exposed to progressively more perturbations or cell
> lines (fig. S4,5), evaluating Prophet with 80%, 50%, 30%, 20%, 10%, and 5% of the perturbations or cell lines
> included in the training dataset." (PDF p. 32)

> "This revealed that confidence in predictions increased as more perturbations and cell states were incorporated
> (fig. S5A). Performance similarly improved with additional data (fig. S6, S7)." (PDF p. 9)

> "The more the experimental space is mapped during training, the better Prophet's performance. However, the concrete
> pattern is dataset dependent as different datasets capture different amounts of variation." (fig. S6 legend, PDF p. 50)

> "Despite its strong performance, Prophet encounters limitations in certain scenarios. For example, predicting the
> effects of observed perturbations in untested cell lines within SCORE, or forecasting the impact of untested small
> molecules in known cell lines within LINCS, can be challenging. Further analysis revealed that the correlation
> between experimental replicates is a key determinant of predictability (Table S4, Methods). … Under such
> circumstances, the mean estimator represents an upper bound on predictive performance." (PDF p. 9)

> "This means that the cell line-specific information gets lost once the replicates are aggregated and, hence, one can
> not distinguish a random replicate from a real one." (Table S4 legend, PDF p. 59)

> "Hit rates for unseen perturbations were generally higher than for unseen cell states … reflecting the relative
> difficulty of extrapolating to new cellular contexts." (PDF p. 6)

> "Prophet is a univariate regression model, so it operates on scalar targets. For high-dimensional data such as
> transcriptomics or Cell Painting features … we decompose the task into multiple scalar regression problems—one per
> feature." (PDF p. 42)

> "This … suggests that model capacity is currently constrained by data availability rather than architectural
> limitations." (PDF pp. 36–37)

### Examples & illustrations (for later drafting)

- **The headline figure for H2a's framing**: fig. S6A, GDSC panel (p. 50) — a flat mean baseline at R² ≈ 0.71 and a
  Prophet curve rising through it at ≈ 25 % of 682 cell lines. This is the published picture our Tahoe curve
  reproduces with a matched ablation and an explicit count.
- **The counter-example in the same figure**: the SCORE panel (p. 50) — the mean baseline is above Prophet at every
  fraction and the two *diverge*, with Table S4 (p. 59) supplying the reliability explanation. The best single
  illustration that a context curve can never cross for reasons that have nothing to do with the number of contexts.
- **The pretraining shift**: GDSC crosses at ≈ 25 % without pretraining (fig. S6A) and below 5 % with it (fig. S7A).
  A crossing point is a property of the training pipeline, not only of the data — which is why our C*_emp must be
  reported with its protocol.
- **Scalar-vs-vector**: Prophet predicts TNFRSF21 mRNA as one scalar task and beats the mean by 25 % on it (p. 9,
  Fig. 2H) — the closest the paper comes to a transcriptomic result, and it is one gene.

### Boundary of this focus

This focus ignores the drug-discovery half of the paper — activity-based fingerprints, MoA clustering, the 1.9 M
ChEMBL screen, the melanoma and plate validations, and the zebrafish work (Kernthese 27, pp. 10–19) — except where it
fixes what the readouts are. A reader looking here for a break-even count, a formula for how many contexts are needed,
a noise-corrected ceiling, a transcriptome-wide evaluation or any correction method will not find one: Prophet
supplies the empirical shape of the context curve and nothing about its position.

## Mentioned entities

- Models / software: [[prophet-universal-predictor]]
- Concepts: [[context-count-learning-curve]], [[simple-perturbation-baselines]], [[zero-shot-transfer-ceiling]]
- Datasets named in the source but without wiki pages: SCORE (CRISPR, CERES, 249 lines), LINCS/L1000, GDSC, GDSC2,
  CTRP, PRISM, JUMP Cell Painting, Horlbeck, Shifrut, Saunders zebrafish, CCLE/DepMap (Table S1, PDF p. 57).

## Connections

- Introduces / is the subject of: [[prophet-universal-predictor]] (PDF pp. 2–3).
- Prior empirical evidence for: [[context-count-learning-curve]] (figs. S6–S7, PDF pp. 50–52) — the only published
  context-count curve plotted against a context-free baseline.
- Defines the cross-context member of: [[simple-perturbation-baselines]] (fig. S1C legend, PDF p. 43).
- Supplies a verbal noise ceiling for: [[zero-shot-transfer-ceiling]] (PDF p. 9; Tables S4–S5, pp. 59–60).
- Compared in this project with: [[dibaeinia-2026-context-scale]] (context curve without a context-free baseline),
  [[molina-2026-response-decomposition]] (1 vs 3 source lines, context-free model), [[piao-2026-adaptive-contexts]]
  (Perturb Mean as a cross-line context-free baseline) — none of them cites Prophet and Prophet cites none of them.
- Cites (not ingested): @replogle-2022-genome-scale (ref. 74); Norman et al. 2019 (ref. 73), Behan et al. SCORE
  (ref. 3), Subramanian et al. LINCS/L1000 (ref. 20), Corsello et al. PRISM (ref. 1), Barretina et al. CCLE
  (ref. 65) — none of these has a key in the project bib (PDF pp. 21–30).
- Project use: the plan's Amendment A3 (`input/ideas/zero-shot-limits-context-complexity-plan.md`) fixes the framing —
  Prophet is prior empirical evidence; the contribution claimed is the explicit break-even with its noise lower bound,
  its measurement, and the de-confounded learning curve.
