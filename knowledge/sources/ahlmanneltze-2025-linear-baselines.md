---
title: "Ahlmann-Eltze, Huber & Anders 2025 — Deep-learning perturbation prediction does not yet beat linear baselines"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: ahlmanneltze-2025-linear-baselines
depth: standard
tags: [perturbation-prediction, benchmark, linear-baseline, mean-baseline, additive-model, foundation-models, genetic-interactions, crispri, crispra]
relations:
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "Introduces the linear model argmin_W ||Y_train − (GWPᵀ + b)||² (eq. 1, p. 1660) and reports «None of the deep learning models was able to consistently outperform the mean prediction or the linear model» (p. 1660)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "Uses «two CRISPR interference datasets by Replogle et al. obtained with K562 and RPE1 cells» (p. 1660), as provided by GEARS (PDF p. 6)."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "Indirect only: the hard-to-beat 'mean' baseline is the within-line average response (b, p. 1660), i.e. the template μ + α_c of the decomposition; the source never decomposes responses."
---

# Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines

## Bibliographic Details

Ahlmann-Eltze, Constantin, Wolfgang Huber, and Simon Anders. 2025. "Deep-learning-based gene perturbation effect
prediction does not yet outperform simple linear baselines." *Nature Methods* 22 (8): 1657–1661.
doi:10.1038/s41592-025-02772-6. Brief Communication; received 11 Oct 2024, accepted 24 June 2025, published online
4 Aug 2025 (p. 1657). CC BY 4.0 (p. 1661). Code: github.com/const-ae/linear_perturbation_prediction-Paper, archived on Zenodo
doi:10.5281/zenodo.14832393 (PDF p. 7). Read from `<library>/pdf/ahlmanneltze-2025-linear-baselines.pdf` (21 pp.).
Citekey: @ahlmanneltze-2025-linear-baselines.

**Page numbering.** The PDF is the version of record. "p. 1657"–"p. 1661" are the printed journal pages (= PDF
pp. 1–5: main text and references). The online Methods (PDF pp. 6–7), Extended Data Figs. 1–10 (PDF pp. 8–18; ED
Fig. 8 spans PDF pp. 15–16) and the reporting summary (PDF pp. 19–21, image-only, not read in detail) carry no printed
page numbers, so they are cited as "PDF p. N".

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. (printed / PDF) |
|---|---|---|
| ✓ F1 | Abstract; motivation and model list (scGPT, scFoundation, GEARS, CPA; repurposed scBERT, Geneformer, UCE with linear decoder) | 1657 |
| ✓ F1 | Double-perturbation benchmark (Norman K562 CRISPRa): no-change and additive baselines; Fig. 1a–b | 1657–1658 |
| ✓ F1 | Genetic-interaction detection and classification (buffering / synergistic / opposite); Fig. 1c–f | 1657–1660 |
| ✓ F1 | Single unseen-perturbation benchmark (Replogle K562 / RPE1, Adamson): linear model (eq. 1), mean baseline; Fig. 2a | 1660, 1659 |
| ✓ F1 | Pretrained embeddings in the linear model; P pretrained on the other Replogle line; Fig. 2c | 1660, 1659 |
| ✓ F1 | Discussion: earlier "easy" comparisons, related benchmarks, limitations | 1660 |
| — | References 1–29 | 1660–1661 |
| ✓ F1 | Methods: data, software versions, per-model protocol, linear decoder | PDF 6 |
| ✓ F1 | Methods: double-perturbation setup, additive model (eq. 2), interaction calling (locfdr) and classes | PDF 6 |
| ✓ F1 | Methods: single-perturbation setup, linear-model fitting (eq. 3), embeddings, additive model as special case | PDF 6–7 |
| ✓ F1 | Methods: metrics (L2, Pearson delta, FDP/TPR curves, eqs. 4–6); data and code availability | PDF 7 |
| ✓ F1 | ED Fig. 1 dataset overview (table, UMAP, target-gene knockdown) | PDF 8 |
| ✓ F1 | ED Figs. 2, 8: alternative metrics (Pearson delta; gene-set size n) | PDF 9, 15–16 |
| — | ED Figs. 3–6: locfdr null, PR/ROC curves, recurring genes, HBG2/HBZ | PDF 10–13 |
| ✓ F1 | ED Fig. 7: per-gene variation of predictions vs truth | PDF 14 |
| ✓ F1 | ED Fig. 9: K562 ↔ RPE1 transfer accuracy vs cross-line differential expression | PDF 17 |
| ✓ F1 | ED Fig. 10: compute cost | PDF 18 |
| — | Reporting summary | PDF 19–21 |

## Kernthesen

### 1. ⭐⭐ None of seven deep models beat deliberately simple baselines (p. 1657; p. 1660)

Five foundation models (scGPT, scFoundation, and the repurposed scBERT, Geneformer, UCE) and two other deep models
(GEARS, CPA) are compared against «deliberately simple baselines for predicting transcriptome changes after single or
double perturbations. None outperformed the baselines» (abstract, p. 1657). The conclusion is about the goal, not the
architectures: «As our deliberately simple baselines are incapable of representing realistic biological complexity,
yet were not outperformed by the foundation models, we conclude that the latter's goal of providing a generalizable
representation of cellular states and predicting the outcome of not-yet-performed experiments is still elusive»
(p. 1660). The authors state they ran each model «as much as possible, with their default parameters», asked the original
authors to review their code, and report CPA's original code because the authors' fix performed worse (PDF p. 6).

### 2. ⭐⭐ Double perturbations: every model has a higher error than the additive baseline (p. 1657; Fig. 1a–b, p. 1658; eq. 2, PDF p. 6)

On Norman et al. (K562, CRISPR activation; 100 single and 124 double perturbations, 19,264 genes; p. 1657), models are
fine-tuned on all 100 singles and 62 doubles and tested on the other 62 doubles, over five random splits (p. 1657; PDF p. 6).
The two baselines use no double-perturbation data (p. 1657): **no change** predicts the control y∅, and **additive**
predicts ŷ_add = y_A + y_B − y∅, the sum of the single-perturbation log fold changes (eq. 2, PDF p. 6). «All models had a
prediction error substantially higher than the additive baseline» (p. 1657). Error is L2 on the 1,000 most highly
expressed genes; Pearson delta and other gene subsets give «the same overall result» (p. 1657; ED Fig. 2, PDF p. 9). In the
CEBPE + CEBPB example the additive model has L2 4.7 / Pearson delta 0.93, scGPT 7.1 / 0.53, CPA 14.8 / 0.14, and
Geneformer* 2.9 / 0.93 (Fig. 1b, p. 1658), so the ranking is a mean over perturbations, not a per-example rule.

### 3. ⭐ No model detects genetic interactions better than "no change" (pp. 1657, 1659; Fig. 1c–f, p. 1658; PDF p. 6)

Interactions are double-perturbation values that deviate from the additive expectation beyond an Efron empirical null
(locfdr): 5,035 of a potential 124,000 (124 doubles × 1,000 genes) at FDR 5 % (p. 1657; PDF p. 6; ED Fig. 3, PDF p. 10). Of all
read-out values, 2.3 % are buffering, 0.6 % synergistic and none opposite (PDF p. 6). «None of the models was better than
the 'no change' baseline» on TPR vs false-discovery proportion (p. 1657; Fig. 1c), with the same ranking under PR and ROC
(ED Fig. 4, PDF p. 11). Models «mostly predicted buffering interactions» and synergy predictions were rare and rarely
correct (p. 1659; Fig. 1f).

### 4. ⭐ Deep-model predictions barely vary across perturbations (pp. 1659–1660; ED Fig. 7, PDF p. 14)

«for most genes, the predictions of scGPT, UCE and scBERT did not vary across perturbations, and those of GEARS and
scFoundation varied considerably less than the ground truth» (p. 1660; ED Fig. 7 on Norman and Replogle K562). All models
except Geneformer and scFoundation predicted LFC ≈ 0 for the HBG2/HBZ double despite strong single effects (pp. 1659–1660;
ED Figs. 5–6). This is the "mean collapse" behaviour that later papers name.

### 5. ⭐⭐ The linear baseline: Y_train ≈ G W Pᵀ + b, fitted in closed form (eq. 1, p. 1660; eq. 3, PDF p. 6; PDF p. 7)

Each read-out gene gets a K-dimensional vector (rows of G) and each perturbation an L-dimensional vector (rows of P);
Y_train has one row per read-out gene and one column per perturbation pseudobulk, and W = argmin ‖Y_train − (GWPᵀ + b)‖²
with **b the vector of row means of Y_train** (eq. 1, p. 1660). By default G = the top K = 10 principal components of
Y_train, and P = the rows of G for the genes perturbed in training; W = (GᵀG + λI)⁻¹Gᵀ(Y_train − b)P(PᵀP + λI)⁻¹ with
ridge λ = 0.1 «for numerical stability» (eq. 3, PDF p. 6; K = 10, PDF p. 7). Prediction for test perturbations is Ŷ = G W P̃ᵀ + b,
with P̃ the rows of G for the test targets (PDF p. 6). External G or P «simply replaced the estimate from the training
data» (PDF p. 7). Consequence: the default model can only predict perturbations whose target gene is itself a read-out
gene, which restricts evaluation to 73 (Adamson), 398 (K562) and 629 (RPE1) perturbations (PDF p. 6).

### 6. ⭐⭐ The "mean" baseline is not consistently beaten either (p. 1660; Fig. 2a, p. 1659; ED Fig. 8, PDF pp. 15–16; PDF p. 6)

«We also included an even simpler baseline, b, the mean across the perturbations in the training set», following
Kernfeld et al. and Csendes et al. (p. 1660); «The 'mean' model calculated the mean of the expression values in the
training data» (PDF p. 6). «None of the deep learning models was able to consistently outperform the mean prediction or
the linear model» (p. 1660; Fig. 2a; ED Fig. 8). The task is unseen single perturbations within one line, using GEARS'
"simulation" split repeated twice (PDF p. 6). scFoundation is excluded (gene vocabulary mismatch) and CPA is excluded
(not designed for unseen perturbations) (p. 1660).

### 7. ⭐ The additive model is a special case of the linear model (PDF p. 7)

«The additive model is a special case of the linear model (equation (1)) where the gene embedding is simply the
single perturbation data … (G = Y_single); the perturbation embedding P is a binary coding … and W is an identity
matrix and b = −y∅» (PDF p. 7). So both baselines live in one bilinear family.

### 8. Atlas-pretrained embeddings help a little; random embeddings almost as much (p. 1660; Fig. 2c, p. 1659)

Plugging scFoundation's or scGPT's gene embedding G, or a GO-derived P from GEARS, into the linear model performs «as
well or better than scGPT and GEARS with their in-built decoders» (p. 1660). The scFoundation/scGPT-G models beat the
mean baseline but «did not consistently outperform the linear model using G and P from the training data» (p. 1660).
«pretraining on the single-cell atlas data provided only a small benefit over random embeddings» (p. 1660).

### 9. ⭐⭐ The only consistent winner uses perturbation data from another cell line (p. 1660; Fig. 2c, p. 1659; PDF p. 7)

«The approach that did consistently outperform all other models was a linear model with P pretrained on the Replogle
data (using the K562 cell line data as pretraining for the Adamson and RPE1 data and the RPE1 cell line for the K562
data)» (p. 1660). P is a 10-dimensional PCA «on the columns of the matrix with the perturbation means of the reference
data» (PDF p. 7). In Fig. 2c its error relative to the mean baseline is roughly 0.85–0.89 in all three datasets (read
from the forest plot, p. 1659; the paper gives no table). Conclusion: «pretraining on perturbation data increased
predictive performance» (p. 1660). **This is not zero-shot cell-line transfer:** G, W and b are still fitted on the
target line's training perturbations (eqs. 1, 3; PDF p. 7); only the perturbation embedding comes from the other line.

### 10. ⭐ Cross-line transfer degrades with K562–RPE1 expression differences (p. 1660; ED Fig. 9, PDF p. 17)

«The predictions were more accurate for genes that were more similar between K562 and RPE1» (p. 1660). ED Fig. 9: baseline
mean expression of shared genes correlates 0.81 between RPE1 and K562 (a); per-read-out-gene absolute error rises with
the genes' cross-line expression difference (slopes 0.08 for K562 → RPE1 and 0.04 for RPE1 → K562) (b); per-perturbation
Pearson delta falls with the target gene's cross-line difference (slopes −0.19 and −0.15) (c). This is the paper's only
quantitative statement about cross-cell-line similarity.

### 11. Metrics: L2 on 1,000 control-expressed genes and Pearson delta (PDF p. 7; ED Figs. 2, 8)

L2 is √Σ_g (ŷ_g − y_g)² over the 1,000 genes most highly expressed in the control (the paper also calls it RMSE).
Pearson delta is cor(ŷ − y∅, y − y∅), after Cui et al.; «Unlike the L2 distance, the Pearson delta metric does not
penalize predictions that are consistently too small or too large in amplitude» (PDF p. 7). Ranking genes by differential
expression «is only possible if access to the ground truth is available and can thus not be applied in real-world use
cases» (ED Figs. 2b, 8b; PDF pp. 9, 16).

### 12. The mean baseline already scores high Pearson delta, and it differs by line (ED Fig. 8a, PDF p. 15)

Read from the beeswarm means (no table): the mean baseline reaches Pearson delta ≈ 0.40 (Replogle K562), ≈ 0.68
(Replogle RPE1) and ≈ 0.73 (Adamson); scGPT is lower (≈ 0.31 / 0.59 / 0.59); the linear model is ≈ 0.44 / 0.68 / 0.72.
A predictor that ignores perturbation identity therefore gets most of the attainable control-referenced correlation,
more so in RPE1 than in K562. The source does not comment on this.

### 13. Earlier positive comparisons were "easy" (p. 1660)

CPA was never designed for unseen perturbations, and «The linear model used in scGPT's benchmark appears to have been
set up such that it reverts to predicting no change over the control condition for any unseen perturbation» (p. 1660).
The result agrees with Kernfeld, Csendes, PertEval-scFM, Bendidi, PerturBench, Li L., Li C. and Wong et al. (p. 1660;
refs. 15–16, 22–27, p. 1661).

### 14. Deep models cost far more compute for no gain (p. 1660; PDF p. 6; ED Fig. 10, PDF p. 18)

The baselines win «despite significant computational expenses for fine-tuning the deep learning models» (p. 1660).
Fine-tuning was capped at 3 days, so scFoundation trained for five epochs (PDF p. 6); ED Fig. 10 gives elapsed time and
peak memory on Norman.

### 15. Self-declared limits: four datasets, cancer lines, no QC repair (p. 1660; PDF p. 6)

«One limitation of our benchmark is that we used only four datasets» (chosen because GEARS, scGPT and scFoundation
used them); «all datasets are from cancer cell lines»; and the authors «did not attempt to improve the original
quality control, for example, by excluding perturbations that did not affect the expression of their own target gene»
(p. 1660). UCE was given ground-truth expression rows for the perturbed genes, accepting possible test leakage in its
favour (PDF p. 6).

## Figures and tables

- Fig. 1 (p. 1658): **double-perturbation L2 per model vs no-change and additive (a)**, CEBPE + CEBPB example (b),
  interaction TPR vs FDP (c), interaction classes (d–e), predicted vs observed deviation from additive (f).
- Fig. 2 (p. 1659): **single unseen-perturbation L2 on Replogle K562 / RPE1 and Adamson (a)**, linear-model schematic
  Y ≈ GWPᵀ (b), **error relative to the mean baseline with pretrained G / P, incl. "LM with P from Replogle" (c)**.
- ED Fig. 1 (PDF p. 8): dataset table — Norman (K562) 19,264 genes × 84,143 cells, 124 doubles + 100 singles; Replogle
  K562 5,000 × 162,264, 1,087 singles; Replogle RPE1 5,000 × 161,423, 1,534 singles; Adamson (K562) 5,060 × 65,899,
  81 singles.
- ED Fig. 8 (PDF pp. 15–16): **Pearson delta per model for single perturbations (a)**; error vs number of read-out genes (b).
- ED Fig. 9 (PDF p. 17): **K562 ↔ RPE1 transfer accuracy vs cross-line differential expression**.
- *Inconsistencies in the source:* (i) ED Fig. 10 gives Norman as «19 624 genes and 81 143 cells grouped into 225
  conditions» (PDF p. 18), against 19,264 genes (p. 1657) and 84,143 cells (ED Fig. 1, PDF p. 8). (ii) Fig. 2's caption counts
  134 / 210 / 24 unseen perturbations (K562 / RPE1 / Adamson; p. 1659), the Methods 398 / 629 / 73 (PDF p. 6); the ratio is
  ≈ 1/3 but the relation is not explained. (iii) ED Fig. 9's caption speaks of «122 double perturbations» for a K562 /
  RPE1 single-perturbation analysis (PDF p. 17). (iv) Fig. 2c's caption dims a point range «if the confidence interval
  contains 0», but the plotted quantity is a ratio centred on 1 (p. 1659).

## What this source does not address

- **No ceiling of any kind.** There is no noise / replicate ceiling, no split-half or replicate reliability, and no
  statement of how close any model, or the baselines, are to what the data allow. Pearson delta and L2 are raw.
- **No unseen-cell-line task.** Every benchmark is within one line (unseen doubles or unseen single perturbations of
  the same line). The only cross-line element is the perturbation embedding P taken from the other Replogle line,
  with G, W and b fitted on the target line (Kernthese 9). There is no context count, no context descriptor and no
  learned context correction.
- **Two cell lines in total** (K562 in three datasets, RPE1 in one; ED Fig. 1, PDF p. 8), all cancer-derived (p. 1660).
  Tahoe, HepG2, Jurkat and non-cancer lines are not used.
- **No decomposition** of responses into shared vs perturbation-specific parts, and no control-PC or residual
  analysis. That the mean baseline is strong is shown, not explained.
- **No data-scaling experiment** (more perturbations, more cells, more lines) — only the qualitative contrast that
  perturbation-data pretraining helps and atlas pretraining barely does (p. 1660).
- **No VCC-style composite** (no DE-based or distributional metric); single-cell distributions are reduced to
  per-condition means (p. 1660).

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 1, 2, 5, 6, 7, 9, 10, 12, 15** (and 4 for mean collapse, 11 for metrics). The hypotheses are those
of the plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3 at the
end).

### What this focus needs beyond the Kernthesen

1. **Exact baselines, in project notation.** (a) *No change*: ŷ = y∅, the destination control (p. 1657; PDF p. 6). (b)
   *Additive* (doubles only): ŷ = y_A + y_B − y∅ (eq. 2, PDF p. 6). (c) *Mean*: b = the row mean of Y_train, i.e. the
   average response over the **same line's** training perturbations (p. 1660; PDF p. 6). In the plan's notation
   (δ = μ + a_t + b_c + i + ε) this is an estimate of μ + b_c (the within-line template of
   [[perturbation-response-decomposition]]) plus the mean of a_t + i over training targets. (d) *Linear model*:
   b + G W Pᵀ with K = 10 PCs of Y_train and ridge λ = 0.1 (eqs. 1, 3; p. 1660; PDF pp. 6–7), i.e. the template plus a rank-10
   perturbation-specific term learned within the line. All four need measured perturbations of the **destination**
   line, except "no change". In a zero-shot unseen-line setting the mean and linear baselines are not available as
   defined; their analogue is a transferred template. The project's "transfer baseline" (plan H4) is therefore not the
   same object as this paper's mean baseline.
2. **Context count.** Two cell lines (K562 in Norman, Replogle K562 and Adamson; RPE1 in Replogle RPE1; ED Fig. 1,
   PDF p. 8). No result is reported as a function of the number of lines, perturbations or cells.
3. **Cross-cell-line transfer.** Only Kernthese 9: P pretrained on the other Replogle line (K562 → RPE1 / Adamson,
   RPE1 → K562), with G, W, b fitted on target-line training perturbations. It is the single consistent winner,
   relative error ≈ 0.85–0.89 of the mean baseline (Fig. 2c, p. 1659, read from plot). Accuracy falls with cross-line
   expression differences (ED Fig. 9, PDF p. 17; overall expression r = 0.81). This is few-shot, target-calibrated
   transfer of the perturbation embedding, not transfer to an unseen line.
4. **Noise / replicate ceiling.** None. No split-half, no replicate correlation, no reliability bound (pp. 1657–1661; PDF pp. 6–7
   searched). The only noise model is the Efron empirical null for calling interactions in doubles (PDF p. 6).
5. **Reader's observation (not in the source).** The mean baseline's Pearson delta is ≈ 0.40 in Replogle K562 and
   ≈ 0.68 in Replogle RPE1 (ED Fig. 8a, PDF p. 15). The ordering matches the template share of response energy in
   [[molina-2026-response-decomposition]] (K562 23 %, RPE1 65 %; Table 2, p. 35): where the template is larger, a
   perturbation-agnostic predictor gets more control-referenced correlation for free. This is consistent with the
   project's per-destination [[zero-shot-transfer-ceiling]], but the two papers use different gene sets and
   processing, so it is an illustration, not evidence.
6. **Literaturguide pre-read verified.** `input/bibliography/literaturguide.md` says «H4: deep models do not beat
   linear/additive baselines». Correct (pp. 1657, 1660), with two refinements: the single-perturbation comparison is
   within-line unseen perturbations, not unseen lines; and the only consistent gain came from *other-line
   perturbation data*, which is the part closest to H4.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **none (background only)** | No ceiling, no noise correction, no error split. It shows that a perturbation-agnostic mean already reaches high Pearson delta (≈ 0.40–0.73), which motivates scoring against a ceiling rather than raw correlation, but it does not do so. | ED Fig. 8a PDF p. 15; PDF p. 7 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **none** | No unseen-line task, no context correction, no Tahoe, two lines. The spirit — added model complexity does not beat a simple predictor — parallels the project's negative FiLM gain at C = 12, but on a different axis (perturbations, not contexts). | pp. 1657, 1660 |
| **H2b** (theory predicts the break-even within ×2) | **none** | No sample-complexity theory. Nothing here bears on the project's C*_theory = 230 vs C*_emp = 37.4. | — |
| **H2c** (genetic C*_min > C_public) | **none** | No context count. It contributes two of the public lines (Replogle K562, RPE1) as data only. | ED Fig. 1 PDF p. 8 |
| **H3a** (transfer residual on destination control PCs) | **none** | No residual or control-PC analysis. Its "mean collapse" finding (predictions barely vary across perturbations) is about prediction variance, not residual geometry. | p. 1660; ED Fig. 7 PDF p. 14 |
| **H3b** (M3 beats five PC-deflation / common-response corrections) | **none** | No correction method. The bias term b is part of the model, not a correction. | eq. 1 p. 1660 |
| **H4** (adding same-line genome-wide data beats the largest modelling gain) | **partial (qualitative precedent, different setting)** | Architecture gains are nil: no deep model beats mean / linear / additive baselines. The one consistent gain comes from adding perturbation data (P from the other Replogle line), and atlas pretraining adds little over random embeddings — "data, not modelling". Differences from H4: within-line unseen perturbations, not an unseen line; the added data are from a *different* line (closer to the plan's related-line cases KOLF ← H1GW than to H1 ← H1GW); L2 / Pearson delta, not the VCC composite; no dose sweep (w). | pp. 1657, 1660; Fig. 2c p. 1659; ED Fig. 9 PDF p. 17 |

### Direct quotes (for the draft)

> "Here we compared five foundation models and two other deep learning models against deliberately simple baselines
> for predicting transcriptome changes after single or double perturbations. None outperformed the baselines" (p. 1657)

> "None of the deep learning models was able to consistently outperform the mean prediction or the linear model"
> (p. 1660)

> "The approach that did consistently outperform all other models was a linear model with P pretrained on the
> Replogle data (using the K562 cell line data as pretraining for the Adamson and RPE1 data and the RPE1 cell line for
> the K562 data)." (p. 1660)

> "Together, these results suggest that pretraining on the single-cell atlas data provided only a small benefit over
> random embeddings, but pretraining on perturbation data increased predictive performance." (p. 1660)

> "Unlike the L2 distance, the Pearson delta metric does not penalize predictions that are consistently too small or
> too large in amplitude and, thus, prioritizes correct prediction of the direction of the expression change." (PDF p. 7)

### Examples & illustrations (for later drafting)

- CEBPE + CEBPB double: additive L2 4.7 vs scGPT 7.1 vs CPA 14.8 (Fig. 1b, p. 1658) — the clearest single picture of
  "simple beats deep".
- "LM with P from Replogle": the only model below the mean baseline in all three datasets (relative error ≈ 0.85–0.89,
  Fig. 2c, p. 1659) — data from another line beat every architecture. Useful as the literature precedent for H4.
- Mean baseline Pearson delta ≈ 0.40 (K562) vs ≈ 0.68 (RPE1) (ED Fig. 8a, PDF p. 15) — how much a perturbation-agnostic
  template already scores, and how much that varies by line.

### Boundary of this focus

This focus does not cover the genetic-interaction analysis (Kernthese 3) beyond noting that it uses an empirical
null, nor the per-model engineering details (PDF p. 6). A reader looking for a ceiling, a context-count curve, a zero-shot
cell-line benchmark or a correction comparison will not find one in this source.

## Mentioned entities

- Concepts: [[simple-perturbation-baselines]], [[perturbation-response-decomposition]], [[zero-shot-transfer-ceiling]]
- Datasets: [[replogle-nadig-crispri-collection]] (Replogle K562 and RPE1 only)

## Connections

- Introduces / supports: [[simple-perturbation-baselines]] (eqs. 1–3, p. 1660; PDF pp. 6–7).
- Uses data from: [[replogle-nadig-crispri-collection]] (K562, RPE1 essential screens as provided by GEARS; PDF p. 6).
- Indirectly consistent with: [[perturbation-response-decomposition]] (the mean baseline is the within-line template;
  inferred) and the per-destination [[zero-shot-transfer-ceiling]] (focus item 5; reader's observation).
- Cited by: [[molina-2026-response-decomposition]], [[dibaeinia-2026-context-scale]].
- Cites (not yet ingested): @replogle-2022-genome-scale, @kernfeld-2025-expression-forecasting,
  @csendes-2025-foundation-models, @wenteler-2024-perteval, @bendidi-2024-transcriptomics-foundation,
  @wu-2025-perturbench, and Li L. et al. (2024 preprint; probably the preprint of @li-2026-systematic-comparison —
  title match, not verified); also Norman et al. 2019, Adamson et al. 2016, Li C. et al., Wong et al. 2025 (refs. 11,
  13–16, 22–27, pp. 1660–1661).
