---
title: "Adduri et al. 2025 (preprint) — State: predicting perturbation responses across contexts"
type: source
created: 2026-11-06
updated: 2026-11-06
status: review
author: llm
bibkey: adduri-2025-state
based_on: preprint
depth: standard
tags: [perturbation-prediction, virtual-cell, state, set-transformer, cell-embedding, context-generalization, few-shot, zero-shot, mean-baseline, cell-eval, tahoe, replogle-nadig, preprint]
relations:
  - target: simple-perturbation-baselines
    type: contradicts
    confidence: extracted
    because: "Claims to overturn the baseline result: «it is the first model to consistently outperform simple linear baselines in generalizing perturbation effects across cellular contexts» (pp. 4–5). On Replogle–Nadig, however, it only «matched the performance of the perturbation mean baseline» on discrimination (p. 10)."
  - target: ahlmanneltze-2025-linear-baselines
    type: builds-on
    confidence: extracted
    because: "Uses its model as a comparator: «a simple linear model (Ahlmann-Eltze et al., 2024)» (p. 7), formalised as Y ≈ GKPᵀ + b1ᵀ (eq. 49, p. 40)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "«Replogle-Nadig consists of 2,024 genetic perturbations applied to four distinct cell lines after filtering perturbations with low on-target efficacy» (p. 7). Table 1 gives 1,677 perturbations and 4 contexts (p. 31)."
  - target: tahoe-100m
    type: cites
    confidence: extracted
    because: "«Tahoe-100M includes perturbation responses from 50 diverse cancer cell lines treated under 1,138 conditions involving 380 distinct drug perturbations» (p. 7). It is also the ST pre-training corpus for the zero-shot task (p. 14)."
  - target: context-count-learning-curve
    type: supports
    confidence: inferred
    because: "Asserts that gains over mean baselines are larger on datasets with more contexts (Tahoe, Parse) than on genetic data «conducted in just a few cell lines» (p. 11). This is a cross-dataset contrast confounded with modality, size and effect strength; there is no sweep over the number of contexts."
---

# Predicting cellular responses to perturbation across diverse contexts with State

> Preprint version (adduri-2025-state, bioRxiv v2, posted 2025-07-10) of adduri-2026-state; claims quoted from the preprint.

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

## Bibliographic Details

Adduri, Abhinav K., Dhruv Gautam, Beatrice Bevilacqua, Alishba Imran, Rohan Shah, Mohsen Naghipourfar, Noam Teyssier,
Rajesh Ilango, et al. (27 authors; corresponding author Yusuf H. Roohani, Arc Institute). 2025. "Predicting cellular
responses to perturbation across diverse contexts with State." bioRxiv preprint, doi:10.1101/2025.06.26.661135,
version 2 posted 10 July 2025 (v1 27 June 2025). CC BY 4.0. Code: github.com/ArcInstitute/State, Cell-load,
Cell-eval, State-reproduce. Weights: huggingface.co/arcinstitute/SE-600M, ST-Tahoe, ST-Parse (p. 48).
Read from `<library>/pdf/adduri-2025-state.pdf` (69 pp.). Citekey: @adduri-2025-state. Published version:
@adduri-2026-state (Cell 189(19):5914–5931.e20, 2026). Do not cite it for the numbers below.

**Page numbers.** The printed page footers match the PDF page numbers, so "p. N" below is both. The main text runs
pp. 1–19, Methods pp. 20–41, OT theory pp. 42–45, interpretability pp. 45–47, statements pp. 48–49, references
pp. 50–59 and supplementary figures pp. 60–69.

**Figure values.** The paper gives no results tables. All baseline comparisons are bar charts with a printed
"% vs best baseline" label (Fig. 2D–J, p. 8; Fig. 3G–H, p. 12; Figs. S6–S8, pp. 65–67). The printed percentages are
quoted exactly. Absolute values marked "≈" were **read off the bars** of 300-dpi renders and are accurate to about
±0.02. Do not cite them as reported numbers.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; 1 Introduction (two noise sources; Eq. 1 heterogeneity model; SE + ST overview; "first to outperform linear baselines") | 1–5 |
| ✓ F1 | 2.1 Building the State Transition model (cell sets, MMD, set-size ablation; Fig. 1, p. 6) | 5–7 |
| ✓ F1 | 2.2 State outperforms baselines across cell contexts (underrepresented-context task; Tahoe / Parse / Replogle–Nadig; Fig. 2, pp. 8–9) | 7–11 |
| ✓ F1 | 2.3 State embeddings enhance zero-shot prediction (SE, embedding benchmark, Tahoe-pretrained zero-shot on 5 query datasets; Fig. 3, pp. 12–13) | 11–15 |
| — | 2.4 Cell-type-specific responses (Tahoe held-out lines, Trametinib in C32, drug-similarity ARI, survival; Fig. 4, pp. 16–17) | 15–17 |
| ✓ F1 | 3 Discussion (incl. limitations) | 17–19 |
| ✓ F1 | 4.1 Data generation process (Eqs. 2–3) | 20 |
| ✓ F1 | 4.2 Task description: 4.2.1 underrepresented-context task (α = 0.30); 4.2.2 zero-shot context task | 20–22 |
| ✓ F1 | 4.3 ST model (cell sets, control pairing, encoders, residual transformer, MMD energy kernel, embedding-space training) | 22–26 |
| — | 4.4 SE model (ESM-2 gene embeddings, top-2,048 expression set, soft-bin expression encoding, CLS + dataset token, dual-axis loss) | 26–30 |
| ✓ F1 | 4.5 Datasets (Table 1 ST datasets; knockdown filtering; Table 2 SE corpus) | 31–32 |
| — | 4.6 Training (ST hyperparameters, Table 3; components, Table 4; fine-tuning; SE hyperparameters and training) | 32–34 |
| ✓ F1 | 4.7 Evaluation: Cell-Eval metrics; DE metrics; embedding metrics; other evaluations; baseline models (4.7.5) | 34–41 |
| — | 5 State Transition as Optimal Transport (Theorems 2–3) | 42–45 |
| — | 6 Latent space of State (attention heads; K562/Jurkat) | 45–47 |
| — | Code availability, acknowledgments, contributions, competing interests | 48–49 |
| — | 7 References | 50–59 |
| ✓ F1 | Supplementary Figs. S1–S10 (S3 Tahoe holdout lines; S4–S5 extra DE metrics; S6 held-out tissue; S7–S8 zero-shot; S9 interpretability; S10 distance choice) | 60–69 |

## Kernthesen

### 1. Two noise sources mask perturbation effects and explain why models fail to generalise (pp. 2–4)

The authors attribute the failure of deep models to beat linear models across contexts to «two major sources of noise
that mask true perturbation effects in single-cell perturbation datasets: biological heterogeneity within the studied
population that is not explained by experimental covariates, and technical or experimental variation across
different perturbation datasets» (p. 2). They formalise this as X̂_p ∼ T̂_p(D_basal) + H(D_basal) + ε, with the
perturbation effect, the basal-population heterogeneity and independent technical noise (Eq. 1, p. 4; Eqs. 2–3,
p. 20). ST is meant to model T̂_p + H and SE to absorb ε (p. 20).

### 2. ST models perturbation as a transformation of a set of control cells, trained with MMD (pp. 5–7, 22–25)

Cells are grouped by context, perturbation and batch into sets of size S. Each perturbed set is paired with a set of
S control cells from the same line, optionally the same batch (eqs. 10–11, p. 22). A bidirectional transformer
without positional encoding maps the control set to the perturbed set as a residual, O = H + f_ST(H) (eq. 17, p. 24).
The loss is the squared MMD with the energy-distance kernel k(u, v) = −‖u − v‖₂ (eq. 19, p. 25). On Tahoe the
validation loss improves with set size up to 256 cells. The full model beats a mean-pooling "pseudobulk" variant, a
set-size-1 variant and a variant without self-attention (Fig. 1D, p. 6; p. 38).

### 3. ⭐ Context enters ST only through the control cells; there is no cell-line token (pp. 22–24)

«The final input to ST is constructed by summing the control cell embeddings with the perturbation and batch
embeddings» (p. 24): H = H_cell + H_pert + H_batch (eq. 16). Perturbations and batches are one-hot encoded («For
State, we use one-hot encodings of perturbation and batch», p. 23). The context label ℓ is used only to form sets and
pick matched controls (pp. 22–23). The batch encoder is switched off for Replogle–Nadig (Table 3, p. 33). Therefore
the model's only information about a destination line is its basal control population. In both tasks «Unperturbed
cells from all contexts including the underrepresented context are available to the model during training»
(pp. 21, 22).

### 4. SE is a 600M-parameter cell embedding trained on 167 million observational cells (pp. 13, 26–34)

SE is a 16-layer, h = 2,048 transformer encoder (p. 33). Its input is the cell's top 2,048 expressed genes, each
represented by an ESM-2 protein embedding plus a soft-binned expression encoding (pp. 26–28). The output is a CLS
embedding concatenated with a 10-d projection of a dataset token (eq. 31, p. 28). Training uses a dual-axis loss
(within cell across genes, and within gene across cells) plus classification of the dataset of origin over about
14,000 datasets (pp. 29–30). The corpus is 14,420 AnnData files from scBaseCount, CZ CELLxGENE and Tahoe-100M (Table 2,
p. 32; p. 34). To avoid leakage, only 20 Tahoe lines are used, separate from the five held-out Tahoe lines (p. 32).

### 5. ⭐⭐ The headline benchmark is few-shot: 30 % of the held-out context's perturbations are in training (pp. 7, 9, 20–21)

The main comparison (Fig. 2) is an «underrepresented context generalization task (Section 4.2.1), in which each model
had access to 30% of perturbations in the test context during training» (p. 7). It is formalised with α = 0.30 and a
support set P_sup ⊂ P_c* (eqs. 4–6, p. 21). The splits differ by dataset. For Tahoe, five lines were chosen from a
PCA of pseudobulks (C32, HOP62, HepG2/C3A, Hs 766T, PANC-1; Fig. S3, p. 62). Parse holds out 4 of 12 donors.
Replogle–Nadig holds out one line at a time: «models were trained on the remaining three cell lines plus an
additional 30% of perturbations randomly sampled from the test cell line» (p. 9). A further Tahoe split holds out
each of 11 organs (p. 9; Fig. S6, p. 65), again with 30 % of test perturbations in training (p. 65).

### 6. The mean baselines are context-free averages with a target anchor (pp. 7, 39–40)

"Perturbation mean" is the target line's control mean plus the perturbation's offset averaged over training lines:
x̂ = μ_t^ctrl + δ_p. «Thus, controls are reproduced exactly, while every non-control perturbation receives the same
global shift» (p. 39). "Context mean" returns the average perturbed expression of the target line over all training
perturbations, i.e. the same prediction for every perturbation (p. 39). In the few-shot task it is estimated from
the 30 % support perturbations. The linear baseline is the low-rank ridge model of Ahlmann-Eltze et al. added to each
cell's control profile (eqs. 48–52, p. 40). The other comparators are scVI, CPA and scGPT (p. 40). All models predict
the top 2,000 HVGs (p. 7).

### 7. ⭐ On Tahoe-100M and Parse-PBMC, State beats every baseline on all six headline metrics (pp. 10–11; Fig. 2D–J, p. 8)

The printed gains over the best baseline for Tahoe / Parse are: perturbation discrimination +54 % / +29 %; Pearson-Δ
+63 % / +47 %; DE AUPRC +36 % / +29 %; log-fold-change Spearman +10 % / +22 %; DE overlap +96 % / +43 %;
effect-size Spearman +70 % / +53 % (Fig. 2, p. 8). On Tahoe State reaches ≈ 0.97 discrimination against
perturbation-mean ≈ 0.62, and ≈ 0.80 effect-size Spearman, «approaching an absolute correlation of 0.8» (p. 11). The
text calls the discrimination gains «an absolute improvement of 54% and 29%» (p. 10), but the bars show they are
relative gains (≈ 0.97 / 0.62).

### 8. ⭐⭐ On Replogle–Nadig the picture is mixed: parity on discrimination, losses on DE overlap and fold change (pp. 10–11; Fig. 2D–J, p. 8)

The printed gains over the best baseline on Replogle–Nadig are: discrimination **−3 %** (State ≈ 0.38 vs
perturbation-mean ≈ 0.39); Pearson-Δ **+5 %** (≈ 0.41 vs context-mean ≈ 0.39); AUPRC **+184 %** (≈ 0.47 vs ≈ 0.17);
log-fold-change Spearman **−12 %** (≈ 0.44 vs context-mean ≈ 0.50); DE overlap **−21 %** (≈ 0.20 vs linear
≈ 0.25); effect-size Spearman **+22 %** (≈ 0.70 vs scVI ≈ 0.58) (Fig. 2, p. 8). In the authors' words: «On genetic
perturbation datasets, State matched the performance of the perturbation mean baseline and significantly outperforms
all other baseline models» (p. 10), and «Linear model outperformed on DE overlap by 20% and the context mean baseline
outperformed on fold change prediction by 12%» (p. 11). They attribute the small Pearson-Δ gain to «the subtler
effects of perturbations» (p. 10).

### 9. ⭐ The authors themselves say that gains are larger where there are more contexts (p. 11)

«the performance gains over mean baselines were notably larger on Tahoe-100M, which includes 100 million cells
spanning thousands of perturbations across dozens of baseline cellular contexts, and Parse-PBMC, which includes
10 million cells across 12 donors and 18 cell types, as compared to the genome-scale genetic perturbation datasets
conducted in just a few cell lines» (p. 11). They read this as «State's ability to leverage data scale and context
diversity more effectively than existing models» (p. 11). The comparison is across datasets that differ in modality,
cell count and effect size at once. No experiment varies the number of contexts within one dataset.

### 10. Effect-size ranking is the metric State is best at, including on genetic data (pp. 11, 18)

The number of DE genes per perturbation is ranked well even when the identities of the DE genes are not recovered:
«even when the specific genes predicted to be differentially expressed after genetic perturbation do not always match
those observed experimentally (Fig. 2I), State can still accurately estimate the overall size of the DE gene set»
(p. 11). The authors propose DE-set size as «a more robust and complementary indicator of model performance in
genetic perturbation experiments» (p. 11).

### 11. SE embeddings separate perturbations better than other foundation-model embeddings (pp. 13–14; Fig. 3D, p. 12)

On two perturbation datasets not seen in SE training (Srivatsan, 3 lines; Replogle–Nadig, 4 lines), a classifier on
SE embeddings separates perturbation labels better than one on scFoundation, scGPT, UCE or TranscriptFormer
embeddings, and better than one on raw expression (intrinsic). «State embeddings also led to over a 6% absolute
improvement in downstream perturbation classification accuracy compared to all baselines» (extrinsic, p. 14). The
authors infer that «SE is, in some cases, capable of denoising Perturb-seq data» (p. 14).

### 12. ⭐⭐ True zero-shot context transfer exists only as a secondary task with Tahoe pre-training (pp. 14, 21–22; Fig. 3E–G)

In the zero-shot task, ST is pre-trained on Tahoe-100M and fully fine-tuned on a query dataset with one context held
out. It is then evaluated on that context with no perturbation data from it (eqs. 7–8, pp. 21–22). The authors
describe this as «focusing on zero-shot context-level generalization within a dataset rather than zero-shot dataset
transfer» (p. 14). The task uses five query datasets with 2,102 perturbations and 40 contexts (p. 14). The main-text
metric is effect-size Spearman (Fig. 3G). For Replogle–Nadig the printed gain is **+23 %** (ST+SE ≈ 0.57 vs ST+HVG
≈ 0.48 and perturbation-mean ≈ 0.47). For Parse it is +17 % (≈ 0.78). For Parse-PBMC and Replogle–Nadig the text
says «using the State embedding achieved more than 17% improvement with an absolute Spearman correlation greater than
0.5» (p. 14). The other two genetic datasets (Jiang, 30 contexts; McFaline-Figueroa, 3 contexts) have baselines «near
zero» (p. 14).

### 13. ⭐ In the zero-shot supplement, Replogle–Nadig's printed gains are over ST+HVG, not over the mean baseline (Fig. S7, p. 66)

In Fig. S7 the percentages are «change between ST + SE and ST + HVG» (p. 66). For Replogle–Nadig they are +45.9 %
discrimination, +58.4 % Pearson, −7.6 % DE overlap, +49.3 % fold change, +24.3 % effect size and +11.9 % AUPRC. The
**perturbation-mean bars are drawn alongside them, and on three metrics they are at or above ST+SE**: discrimination
≈ 0.37 vs ≈ 0.35, Pearson ≈ 0.25 vs ≈ 0.18, fold-change Spearman ≈ 0.33 vs ≈ 0.22. DE overlap is about 0.10 for all
three. ST+SE leads only on effect size (≈ 0.58 vs ≈ 0.47) and AUPRC (≈ 0.18 vs ≈ 0.08). The text does not mention
this.

### 14. The claim that pre-training helps on genetic data is not borne out for Replogle–Nadig effect size (pp. 14–15; Fig. S8A, p. 67)

The text says that «ST+SE model performance also consistently benefited from pre-training even for datasets without
drug-based perturbations (e.g. genetic or signaling perturbations)» (pp. 14–15). Fig. S8A prints **−14.8 %** for
Replogle (ST+SE pre-trained vs not pre-trained). Read off the bars, ST+SE without pre-training is ≈ 0.67, ST+SE
pre-trained ≈ 0.57, ST+HVG without pre-training ≈ 0.65, ST+HVG pre-trained ≈ 0.46 and perturbation-mean ≈ 0.47. The
Fig. 3G / S7 values (ST+SE ≈ 0.57, ST+HVG ≈ 0.46–0.48) therefore appear to be the pre-trained arms. A non-pre-trained
HVG model would score about as well as the headline ST+SE. This is a reader inference from matching bar heights; the
captions do not say which arms Fig. 3G uses.

### 15. Cell-type-specific DE genes are recovered on held-out Tahoe lines (pp. 15–17; Fig. 4)

On the five held-out Tahoe lines (few-shot split), State beats both mean baselines on DE overlap and fold-change
Spearman across perturbations (Fig. 4B–C). For 0.05 µM Trametinib it recovered «over 30% of the top DE genes, while
both baseline methods identified fewer than 7%» in every line, about 20 more line-specific DE genes than the
baselines. On average only 16 of the top-100 DE genes are shared between any two lines (p. 15). The caption text
gives 0.05 µM (Fig. 4, p. 16) and the main text 0.5 µM (p. 15), an inconsistency in the source.

### 16. Predicted profiles preserve drug-similarity structure and predict survival (p. 17; Fig. 4G–H)

The drug-similarity clustering from State predictions has ARI «of over 0.7, compared to less than 0.1 for the
perturbation mean baseline». A survival regressor on predicted profiles reaches Pearson 0.52, against 0.2 for
context-mean and 0.31 for perturbation-mean inputs (p. 17; SVR with 35 held-out drugs, 10 repeats, pp. 38–39).

### 17. Cell-Eval: discrimination uses Manhattan distance, and Pearson-Δ uses absolute deltas (pp. 34–37; Fig. S10)

Perturbation discrimination is a normalised rank of the true perturbation among all true profiles, reported as
1 − 2·PDisc so that random = 0 (pp. 34–35). It uses Manhattan rather than cosine distance because cosine «reflects
excessive similarity» (Fig. S10, p. 69; p. 34). Pearson-Δ is computed on «∆t = |p̄t − c̄|, the element-wise absolute
difference between the mean expression of perturbed cells (p̄t) and control cells (c̄)» (p. 35). DE genes come from
Wilcoxon tests with BH correction at 0.05 (p. 36). Effect size is the Spearman correlation of DE-gene counts (eq. 47,
p. 37).

### 18. ST contains the optimal-transport map, which Jacobian constraints make unique (pp. 5, 18, 42–45)

Under regularity conditions and asymptotically, «the unique continuous OT map between unperturbed and perturbed cell
populations lies within the solution family of State» (p. 5; Theorem 2). With a symmetric positive-semidefinite
Jacobian constraint the family contains only that map (Theorem 3, p. 45). The authors note that the unconstrained
family «contains infinite additional elements», so OT «is not guaranteed to be learned» (p. 45).

### 19. Attention heads split into heterogeneity-sensitive and pseudobulk-like ones (pp. 45–47; Fig. S9)

On Replogle–Nadig, most of the change in hidden state happens in ST's final layer. Its heads are either
heterogeneity-sensitive or heterogeneity-insensitive. The insensitive ones «mimic operations akin to pseudobulking …
or PerturbMean-style computations» (p. 46). A mixed 50 % K562 / 50 % Jurkat input set polarises the heads (Fig. S9D,
p. 68). SE places K562 and Jurkat close together without labels (Fig. S9A; p. 46).

### 20. Stated limits: no held-out-dataset test, and subtle metrics depend on data size and quality (pp. 18–19)

«Although State has been trained and evaluated on over 70 cell contexts, we have not tested its usability on entirely
held-out datasets for which no cell contexts have been seen during training» (p. 19). «While State shows strong
zero-shot performance when predicting perturbation effect sizes, metrics that reflect subtler changes, such as the
accuracy of individual DE genes, are more sensitive to dataset size and quality» (p. 18). The authors expect larger
perturbation atlases to improve these metrics (p. 19).

### 21. Genetic data are filtered to strong knockdowns before training (p. 31)

Perturbations are kept if mean residual target expression is ≤ 0.30. Cells are kept if their own residual expression
is ≤ 0.50, and perturbations with fewer than 30 valid cells are dropped (p. 31). The Replogle–Nadig count after
filtering is given as 2,024 in the text (p. 7) and 1,677 in Table 1 (p. 31) and the figure labels (Fig. 2, p. 8;
Fig. 3, p. 12). The preprint never names the four lines. Only K-562 and Jurkat are named (p. 46).

## Figures and tables

- **Fig. 1** (p. 6): A noise sources; B multi-scale architecture; C ST schematic; D validation loss vs cell-set size
  and ablations (pseudobulk, set size 1, no self-attention).
- **Fig. 2** (pp. 8–9): A underrepresented-context task and the two mean baselines; B datasets; C Cell-Eval; **D–J
  bar charts per dataset (Tahoe / Parse / Replogle–Nadig rows) with % vs best baseline**; K effect-size confusion
  heatmaps.
- **Fig. 3** (pp. 12–13): A SE + ST; D embedding benchmark on Srivatsan + Replogle–Nadig; F zero-shot fine-tuning
  scheme; **G zero-shot effect-size Spearman on 5 query datasets**; H Parse zero-shot, all metrics.
- **Fig. 4** (pp. 16–17): Tahoe cell-type-specific DE, Trametinib/C32, similarity ARI, survival.
- **Table 1** (p. 31): ST datasets. Replogle–Nadig 624,158 cells / 1,677 perturbations / 4 contexts; Jiang 30
  contexts; Srivatsan 3; McFaline-Figueroa 3; Tahoe 50; Parse 12 donors / 18 cell types.
- **Table 2** (p. 32): SE corpus (71.7 M scBaseCount, 59.2 M CELLxGENE, 36.2 M Tahoe training cells).
- **Table 3** (p. 33): ST for Replogle–Nadig has set size 32, hidden 128, GPT-2 backbone and 10 M parameters, against
  244 M for Tahoe and Parse.
- **Figs. S4–S5** (pp. 63–64): AUPRC, overlap@k, precision@k. **Fig. S6** (p. 65): held-out tissue. **Fig. S7**
  (p. 66): zero-shot, 6 metrics × 5 datasets, % = ST+SE vs ST+HVG. **Fig. S8** (p. 67): effect of pre-training.

## What this source does not address

- **No context-count manipulation.** The number of training contexts is never varied within a dataset. The
  "more contexts → larger gain" reading (p. 11) is a cross-dataset contrast. The sweeps in the paper are over cell-set
  size (Fig. 1D) and pre-training (Fig. S8), not over lines or perturbations.
- **No noise ceiling or reliability correction.** Performance is reported relative to the best baseline, never
  relative to a split-half or replicate ceiling. There is no variance decomposition into conserved and interaction
  parts (cf. [[molina-2026-response-decomposition]]).
- **No residual analysis** and no analysis of where errors sit relative to the destination's control state. Controls
  are only used as the model input (Kernthese 3).
- **No systematic-variation correction** of the metrics. Pearson-Δ is taken against the control (on absolute deltas),
  which is the setting critiqued in [[vinastorne-2026-systema]].
- **No hESC, no CRISPR-KO, no VCC.** The preprint does not mention the Virtual Cell Challenge (see
  @adduri-2026-virtual-cell-challenge).
- **No zero-shot test on a held-out dataset** (p. 19), and the zero-shot task always fine-tunes on the query dataset's
  other contexts.

<!-- ═══ FOKUS ═══ -->

## Focus: facts and claims the manuscript cites; exactly what overlaps with H1–H4 — 2026-11-06

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

**Carries Kernthesen 3, 5, 6, 8, 9, 12, 13, 14, 20, 21** (and 7 for the Tahoe/Parse contrast). The hypotheses are
those of `input/ideas/zero-shot-limits-context-complexity-plan.md`.

### What this focus needs beyond the Kernthesen

1. **Verification of "STATE's own gain on the genetic Replogle–Nadig collection is small"**
   (`knowledge/synthesis/related-work-positioning.md:26–27`). This is **supported in direction, but needs rewording
   on two points.**
   (a) The Fig. 2 benchmark is **not zero-shot**. It is leave-one-line-out on 4 lines **with 30 % of the held-out
   line's perturbations in training** (p. 9; eqs. 4–6, p. 21). So "define the zero-shot unseen-context benchmark" is
   wrong for this paper.
   (b) The gain is **not uniformly small, it is mixed**. Against the best baseline per metric (Fig. 2, p. 8):
   discrimination −3 % (parity with perturbation-mean), Pearson-Δ +5 %, log-fold-change Spearman −12 % (context mean
   better), DE overlap −21 % (linear better), but DE AUPRC +184 % and effect-size Spearman +22 %. For Tahoe the same
   six metrics give +54 / +63 / +36 / +10 / +96 / +70 %, and for Parse +29 / +47 / +29 / +22 / +43 / +53 %. The
   authors acknowledge the contrast (p. 11, Kernthese 9).
   (c) In the **true zero-shot** task (held-out line, Tahoe pre-training, SE), the main text reports only
   effect-size ranking for Replogle–Nadig: +23 % over the best of ST+HVG / perturbation-mean (≈ 0.57 vs ≈ 0.48;
   Fig. 3G, p. 12). In Fig. S7 (p. 66) the perturbation-mean bars match or exceed ST+SE on discrimination, Pearson
   and fold change (Kernthese 13). Pre-training lowers Replogle effect-size Spearman by 14.8 % (Fig. S8A, p. 67).
   **Proposed replacement sentence** (for `related-work-positioning.md` §1):
   > STATE [@adduri-2025-state] benchmarks cross-context transfer mainly as an *underrepresented-context* task:
   > on the four-line Replogle–Nadig collection each line is held out but 30 % of its perturbations stay in training
   > (pp. 9, 21). There STATE only matches the perturbation-mean baseline on perturbation discrimination (−3 %), gains
   > 5 % on Pearson-Δ and loses on DE-gene overlap (−21 % vs linear) and fold-change Spearman (−12 % vs context mean),
   > while gaining on DE-gene AUPRC (+184 %) and effect-size ranking (+22 %) (Fig. 2D–J, p. 8; pp. 10–11). The
   > authors note that the gains over mean baselines are far larger on Tahoe-100M (50 lines) and Parse-PBMC than on
   > genetic screens «conducted in just a few cell lines» (p. 11). Its zero-shot test (unseen line, Tahoe
   > pre-training) reports only effect-size ranking in the main text (+23 %, Fig. 3G, p. 12).
2. **Context counts.** Replogle–Nadig has 4 contexts (Table 1, p. 31), Tahoe 50 with 5 held out (p. 7; Fig. S3),
   Parse 12 donors / 18 cell types with 4 donors held out (p. 9), and Jiang 30 contexts. The zero-shot task totals
   40 contexts over 5 datasets (p. 14), and the paper as a whole has "over 70 cell contexts" (p. 19). None of these
   is varied.
3. **How the destination enters.** ST has no context token. The destination is represented only by its own control
   cells, which are available in training for every context, including held-out ones (pp. 21–24; Kernthese 3). This
   matters for H3a (see the table).
4. **Baseline semantics.** STATE's "perturbation mean" is our context-free model μ + a_t anchored on the destination
   controls (x̂ = μ_t^ctrl + δ_p, p. 39). In the few-shot task its "context mean" uses the destination's 30 % support
   perturbations, so it is not available zero-shot.
5. **Model size on genetic data.** The Replogle–Nadig ST has 10 M parameters (set size 32, hidden 128, GPT-2) against
   244 M for Tahoe and Parse (Table 3, p. 33). Capacity therefore differs by dataset in the comparison of Kernthese 9.

### Overlap with H1–H4

| Hyp. | What the source shows exactly | What we add | Overlap level |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error transferable) | No ceiling, no noise correction, no split-half. Performance is only relative to the best baseline (Fig. 2). The authors say subtle metrics are «more sensitive to dataset size and quality» (p. 18). | A noise-corrected ceiling and the attained fraction. | **not addressed** |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | Tahoe has 50 lines with a fixed 5-line few-shot test (p. 7) and is the zero-shot pre-training corpus (p. 14). There is no sweep over C. The qualitative claim is that gains are larger «across dozens of baseline cellular contexts» than in «just a few cell lines» (p. 11), confounded across datasets. | A matched context-free vs context-aware ablation over C on one dataset, with a zero crossing. | **not addressed (adjacent qualitative claim)** |
| **H2b** (theory break-even within ×2) | No theory of context sample complexity. The OT theory (pp. 42–45) concerns the control → perturbed map. | — | **not addressed** |
| **H2c** (genetic C*_min > C_public) | On a 4-line genetic collection a 10 M-parameter ST only matches perturbation-mean on discrimination, even few-shot, and the zero-shot supplement shows perturbation-mean ≥ ST+SE on 3 of 6 metrics (Fig. 2, p. 8; Fig. S7, p. 66). The authors attribute the weak genetic gains to subtle effects (p. 10) and few lines (p. 11). This is consistent with H2c but gives no C*_min or C_public. | The count C*_min vs the census of public genetic contexts. | **not addressed (consistent adjacent evidence)** |
| **H3a** (transfer residual loads on the destination's top-8 control PCs) | The architecture presupposes that the destination's control population carries the context information: ST sees only control cells + a one-hot perturbation (p. 24), and destination controls are always available (pp. 21–22). There is no residual decomposition, no PC basis and no other-line/permuted null. | Evidence that the error left after transfer sits in the destination control subspace. | **not addressed (design precedent only)** |
| **H3b** (M3 beats five deflation / common-response corrections on PDS) | Its perturbation discrimination score (Manhattan, normalised inverse rank, pp. 34–35) is the PDS family we score with. No correction method is compared. | The correction comparison. | **not addressed (metric shared)** |
| **H4** (same-line genome-wide data beats the largest modelling gain) | The headline setting gives every model 30 % of target perturbations, so target data are assumed in the main benchmark. On Replogle–Nadig effect size, few-shot State (≈ 0.70, Fig. 2J) sits above zero-shot ST+SE (≈ 0.57, Fig. 3G). This is a reader's cross-figure comparison with different pipelines (pre-training, SE). There is no controlled "add target data vs improve model" contrast. | The controlled H1GW contrast on the VCC composite. | **partly shown (qualitative, inferred)** |

### Direct quotes (for the draft)

> "we implemented an underrepresented context generalization task (Section 4.2.1), in which each model had access to
> 30% of perturbations in the test context during training (Fig. 2A)." (p. 7)

> "For the Replogle-Nadig dataset, we conducted an evaluation by iteratively holding out one cell line as a test set.
> For each iteration, models were trained on the remaining three cell lines plus an additional 30% of perturbations
> randomly sampled from the test cell line." (p. 9)

> "On genetic perturbation datasets, State matched the performance of the perturbation mean baseline and
> significantly outperforms all other baseline models." (p. 10)

> "On this metric, State outperformed baselines by 63% on Tahoe-100M, 47% on Parse-PBMC and 5% on Replogle-Nadig."
> (p. 10)

> "the performance gains over mean baselines were notably larger on Tahoe-100M, which includes 100 million cells
> spanning thousands of perturbations across dozens of baseline cellular contexts, and Parse-PBMC, which includes
> 10 million cells across 12 donors and 18 cell types, as compared to the genome-scale genetic perturbation datasets
> conducted in just a few cell lines." (p. 11)

> "we held out one cell context at a time from the query datasets, thus focusing on zero-shot context-level
> generalization within a dataset rather than zero-shot dataset transfer." (p. 14)

> "Although State has been trained and evaluated on over 70 cell contexts, we have not tested its usability on
> entirely held-out datasets for which no cell contexts have been seen during training." (p. 19)

### Examples & illustrations (for later drafting)

- **One dataset, six verdicts:** Replogle–Nadig few-shot. −3 % discrimination, +5 % Pearson-Δ, +184 % AUPRC, −12 %
  fold change, −21 % DE overlap, +22 % effect size (Fig. 2, p. 8). The choice of metric decides whether State "wins".
- **Mean baseline ≥ model in zero-shot:** Replogle–Nadig Fig. S7 (p. 66). Perturbation-mean ≈ 0.25 vs ST+SE ≈ 0.18
  Pearson, while the printed label (+58.4 %) compares ST+SE with ST+HVG.
- **Pre-training can hurt:** Replogle effect size −14.8 % with Tahoe pre-training (Fig. S8A, p. 67).

### Boundary of this focus

This focus does not cover the OT theory (Kernthese 18), the SE architecture details (4) or the Tahoe drug case study
(15–16), except as background. A reader looking for a context-count curve, a ceiling or a residual analysis will not
find one in this source.

## Preprint vs published version

- **Version history** (bioRxiv API, `api.biorxiv.org/details/biorxiv/10.1101/2025.06.26.661135`): v1 2025-06-27, v2
  2025-07-10. Both are linked to doi:10.1016/j.cell.2026.07.052, and the abstracts of v1 and v2 are identical in the
  API. The file on disk is v2 (PDF creation date 2025-07-09; footer "this version posted July 10, 2025").
- **Published version** (Cell 189(19):5914–5931.e20, epub 2026-08-31; PubMed 42673963; compared only to note
  differences): 32 authors vs 27. New authors include Verma, Chalissery, and three from Myllia Biotechnology
  (Fellinger, Krejci, Bürckstümmer), which suggests added material. The abstract drops "transformer model" and
  "trained using gene expression data from over 100 million perturbed cells", adds "compared with baselines", and
  replaces the Cell-Eval clause ("highlights State's ability to detect cell type-specific perturbation responses, such
  as cell survival") with "can be used to evaluate future models".
- **Not audited:** whether the published figures, the Replogle–Nadig numbers, the split definitions or Fig. S9D differ.
  The literaturguide's note that "Fig. S9D scales perturbations, not lines" does **not** match this preprint, where
  Fig. S9D is an attention map for a mixed K562/Jurkat set (p. 68). It may refer to the Cell version.

## Mentioned entities

- Datasets: [[replogle-nadig-crispri-collection]], [[tahoe-100m]]; Parse-PBMC, Jiang, McFaline-Figueroa, Srivatsan
  (no pages).
- Concepts: [[simple-perturbation-baselines]], [[context-count-learning-curve]], [[systematic-variation]],
  [[zero-shot-transfer-ceiling]] (not addressed).

## Connections

- Contradicts (claim level): [[simple-perturbation-baselines]]. It claims to be the first model to consistently beat
  linear baselines (pp. 4–5), but on genetic data it only reaches parity with perturbation-mean (p. 10).
- Builds on: [[ahlmanneltze-2025-linear-baselines]] (linear comparator, p. 40).
- Cites / uses: [[replogle-nadig-crispri-collection]] (4 lines, p. 7, Table 1), [[tahoe-100m]] (50 lines, p. 7).
- Supports (inferred, cross-dataset only): [[context-count-learning-curve]] (p. 11).
- Evaluated by: [[molina-2026-response-decomposition]] (STATE as a zero-shot comparator in Table 17, p. 40 there).
