---
title: "Li et al. 2024/2026 — A systematic comparison of single-cell perturbation response prediction models (bioRxiv v3)"
type: source
created: 2026-11-06
updated: 2026-11-06
status: review
author: llm
bibkey: li-2024-systematic-comparison
based_on: preprint
depth: standard
tags: [perturbation-prediction, benchmark, evaluation-metrics, simple-baselines, cell-type-transfer, combinatorial-perturbation, signal-decomposition, variance-compression, unseen-perturbation]
relations:
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "«simple baselines proved highly competitive with complex deep learning models … additive Mean and Linear models demonstrated superior performance in seen1 and seen2 scenarios» (p. 9, ll. 217–221)."
  - target: ahlmanneltze-2025-linear-baselines
    type: supports
    confidence: extracted
    because: "«simple baselines proved highly competitive with complex deep learning models, which is consistent with a previous study (25)» (p. 9, ll. 217–219); ref. 25 = Ahlmann-Eltze et al. 2025 (p. 35)."
  - target: perturbation-response-decomposition
    type: supports
    confidence: extracted
    because: "Diagnostic decomposition Y_{c,p} = B_c + G_c + S_p + I_{p1,p2} + C_{c,p} + ε, with «C_{c,p} … the cell-state-specific residual of the perturbation response, corresponding to the cell-state-by-perturbation interaction term» (pp. 30–31, ll. 960–968). Not noise-corrected, not orthogonal (p. 31, ll. 987–989)."
  - target: vinastorne-2026-systema
    type: cites
    confidence: extracted
    because: "Ref. 26 (p. 35); compared in Fig. 5E as «Reference redesign + centroid accuracy» (p. 49; legend p. 44)."
  - target: systematic-variation
    type: supports
    confidence: inferred
    because: "G_c, «the systematic shift shared across perturbations within the same cell state» estimated as «the average change within each cell state» (pp. 30–31), is the same object as Systema's shared perturbed-vs-control shift; the source does not name the concept."
---

# A Systematic Comparison of Single-Cell Perturbation Response Prediction Models

> Preprint version (li-2024-systematic-comparison, bioRxiv v3 posted 2026-08-03) of li-2026-systematic-comparison; claims quoted from the preprint.

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

## Bibliographic Details

Li, Lanxiang, Yue You, Yunlin Fu, Wenyu Liao, Xueying Fan, Shihong Lu, Ye Cao, Bo Li, Wenle Ren, Jiaming Kong,
Shuangjia Zheng, Jizheng Chen, Xiaodong Liu, and Luyi Tian. "A Systematic Comparison of Single-Cell Perturbation
Response Prediction Models." bioRxiv preprint doi:10.1101/2024.12.23.630036, **version 3, posted 3 August 2026**
(v1 2024-12-23, v2 2025-09-16). CC BY-NC-ND 4.0. Guangzhou National Laboratory / Westlake University. Code:
github.com/TianGzlab/scPerturbBench; data and scripts Zenodo 10.5281/zenodo.20823660; website
luyitian.github.io/PerturbArena, which «may contain new models and show different results compared to the
manuscript» (p. 41). Published version: @li-2026-systematic-comparison (Science Advances 12(37):eaed3414, 2026; not read).
Read from `<library>/pdf/li-2024-systematic-comparison.pdf` (49 pp.). Citekey: @li-2024-systematic-comparison.

**Page numbers.** For pp. 1–44 the printed page number equals the PDF page number. The manuscript also has running
line numbers (1–1034, pp. 3–32), cited here as "ll.". Pp. 45–49 are unnumbered figure pages: Fig. 1 = p. 45,
Fig. 2 = p. 46, Fig. 3 = p. 47, Fig. 4 = p. 48, Fig. 5 = p. 49. Figure legends are on pp. 41–44 and Table 1 is on p. 44.
The supplementary materials (figs. S1–S27, Tables S1–S5, Table S2 model descriptions) are **not** in this PDF. Numbers
read off bar charts are marked "≈ (read from figure)".

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Title, affiliations | 1 |
| ✓ F1 | Abstract, teaser | 2 |
| ✓ F1 | Introduction (prior benchmarks; VCC, Cell-Eval, Tahoe-100M) | 3–4 |
| ✓ F1 | Results: benchmarking models and experimental design (25 datasets, 3 tasks, 13 models, baselines, 24 metrics) | 4–7 |
| ✓ F1 | Results: prediction of unseen single-gene perturbations (Task 1) | 7–8 |
| ✓ F1 | Results: perturbation interaction prediction (Task 2) | 8–9 |
| ✓ F1 | Results: cell state transfer (Task 3) | 9–10 |
| ✓ F1 | Discussion (signal-decomposition framework, three principles, limitations, future designs) | 10–14 |
| — | Methods: ethics; cell culture and drug treatment (Perturb KHP, Perturb cmo V1); CellPlex multiplexing | 14–15 |
| ✓ F1 | Methods: datasets and preprocessing (incl. Tahoe filtering, TRADE top/bottom-5 selection, HVGs) | 16–17 |
| ✓ F1 | Methods: baseline models (Baseline1/2/3 MLP, ContextMean, AdditiveMean, PerturbMean, Linear; eqs. 1–9) | 18–19 |
| — | Methods: foundation, GRN, VAE and other deep-learning models (settings) | 19–22 |
| ✓ F1 | Methods: evaluation metrics (eqs. 10–22) | 22–27 |
| — | Methods: ranking-based meta-analysis | 27–28 |
| — | Methods: TRADE, E-distance | 28–29 |
| ✓ F1 | Methods: statistical analysis (task-specific ANOVA; granular variance decomposition) | 29–30 |
| ✓ F1 | Methods: response-component decomposition analysis | 30–31 |
| ✓ F1 | Methods: split protocols for unseen genes, combinations, unseen cell types | 31–32 |
| — | References (63) | 33–40 |
| — | Acknowledgments, funding, contributions, competing interests | 40–41 |
| ✓ F1 | Data and code availability | 41 |
| ✓ F1 | Figure legends (Figs. 1–5); Table 1 (datasets) | 41–44 |
| ✓ F1 | Figures 1–5 | 45–49 |

## Kernthesen

### 1. ⭐ A task-spanning benchmark: 13 models, task-specific baselines, 25 datasets, 3 tasks, 24 metrics (p. 2; pp. 4–7; Table 1, p. 44)

The authors argue that earlier benchmarks were «limited in dataset diversity and task coverage» (p. 3, ll. 32–35). Their
benchmark uses 23 public and 2 in-house datasets (Fig. 1A legend, p. 41). These cover CRISPR-Cas9, CRISPRi, CRISPRa, ORF
overexpression, drugs, IFN-β, Hpoly, LPS and development, in human, mouse, rat, rabbit and cross-species data (Table 1,
p. 44). The 13 models fall into four classes (p. 6, ll. 112–118): foundation models (scGPT, scFoundation), GRN-based
models (GEARS, scELMo), VAE-based models (PerturbNet, scGen, trVAE, scVIDR, scPRAM) and other deep models (CPA,
scPreGAN, BioLord, Squidiff). The baselines are «simple statistical baselines (context-mean and perturbation-mean
predictors) and model-based baselines (linear regression and MLP-based neural network models)» (p. 6, ll. 119–121).
Metrics are reported on all genes and on the top 20/50/100 DE genes (p. 6, ll. 131–133).

### 2. Effect size varies widely between and within datasets (p. 5; Fig. 1B, p. 45)

«broad dispersion in the number of perturbations, cell numbers, and perturbation effect sizes, with substantial
heterogeneity not only between studies but also among conditions within the same study» (p. 5, ll. 73–75).
Well-powered screens such as Replogle RPE1 and the curated Virtual Cell Challenge set often show more than 500 DEGs per
condition (p. 5, ll. 76–79). Some immune subsets respond weakly to specific perturbations (p. 5, ll. 79–82).

### 3. ⭐ Absolute-expression metrics are high and barely separate models; small effects score best on them (p. 7; Fig. 2B–C, p. 46)

For unseen single-gene perturbations, Pearson correlation on expression exceeds 0.95 for most datasets, with «generally
small performance differences between models» (p. 7, ll. 163–165). Datasets and perturbations with small effects score
higher. Large-effect datasets such as Replogle and Sunshine score lower (p. 7, ll. 165–171). The authors attribute this
to small effects being «more difficult to distinguish from technical noise» (p. 7, ll. 172–174).

### 4. ⭐ Delta metrics reverse the ordering: larger effects are predicted better (p. 7; Fig. 2D)

«When moving to an analysis based on expression changes relative to the control, the relationship was reversed. Delta-
based Pearson correlation and delta direction accuracy were strongly positively correlated, and perturbations with
larger effect sizes yielded better performance» (p. 7, ll. 175–178). Fig. 2D shows this for PerturbNet on Replogle.

### 5. ⭐ Fine-tuned foundation models compress variance toward the dataset mean (pp. 7–8; Fig. 2E)

«Foundation models exhibited variance compression after fine-tuning, reducing the spread of predicted deltas toward a
dataset-wide mean» (p. 7, ll. 179–180). This shows up as higher MMD and Wasserstein distances. «scGPT converges to the
population average and has worse distribution alignment as the training epoch increases» (p. 8, ll. 181–183). The
authors' explanation is that distribution preservation is not part of the fine-tuning objective (p. 8, ll. 183–185).
The abstract generalises this to a «conservative bias» shared by all models (p. 2).

### 6. PerturbNet is the only model that recovers DE genes reliably (p. 8; p. 9; Fig. 2F, Fig. 3E)

PerturbNet «emerged as the clear top performer, achieving a mean overlap accuracy of 0.56, whereas most other methods
failed to recover DEGs reliably in several datasets» (p. 8, ll. 189–191). It also leads DE overlap in Task 2 (p. 9,
ll. 232–236). *Internal inconsistency:* Methods says PerturbNet «requires chemical structure information (SMILES
representations) as input and is therefore only applicable to drug perturbation datasets» (p. 20, ll. 627–629). Yet
Tasks 1 and 2 are genetic tasks, and Fig. 2D uses PerturbNet on the Replogle CRISPRi screen (legend, p. 42). The
preprint does not resolve this.

### 7. In Task 1 the data explain absolute and delta metrics, and the model explains DE accuracy (p. 8; Fig. 2G; pp. 29–30)

An ANOVA on metric values (value ~ C(perturbation) + C(dataset) + C(model), p. 29) shows that «dataset and perturbation
identity were primary factors for absolute and delta metrics, while DE accuracy was strongly driven by model choice,
dominated by PerturbNet» (p. 8, ll. 192–194). The granular version bins sequencing depth, effect size, KO/KD efficiency,
preprocessing, loss and latent architecture into three levels each (p. 30, ll. 932–946). Dataset attributes,
«particularly perturbation efficiency», explain numerical accuracy, and latent-space architecture explains DE accuracy and
distributional alignment (p. 8, ll. 197–201). **This "variance decomposition" is of metric scores across benchmark runs.
It is not a decomposition of expression responses** (that is Kernthese 15).

### 8. ⭐⭐ For combinations, additive-mean and linear baselines beat complex models when the singles were seen (pp. 8–9; Fig. 3B, p. 47)

Seen levels are seen0, seen1 and seen2, following GEARS (p. 8, ll. 204–207; p. 32). Most models improve from seen0 to
seen2 (p. 8, ll. 213–216). «simple baselines proved highly competitive with complex deep learning models, which is
consistent with a previous study (25). Specifically, additive Mean and Linear models demonstrated superior performance
in seen1 and seen2 scenarios, often surpassing more complex architectures in these settings where constituent
perturbations were known» (p. 9, ll. 217–222). In Fig. 3B (delta Pearson), Baseline2_AdditiveMean is near 0.1–0.2 at
seen0 but ≈ 0.8–0.9 at seen2 (read from figure). The Methods defines AdditiveMean as x_ctrl + Δ_A + Δ_B (eqs. 3–5,
pp. 18–19).

### 9. Buffering and synergistic combinations swap ranks depending on the metric (p. 9; Fig. 3C–D; p. 12)

Combinations are classed by whether the observed delta is smaller (buffering) or larger (synergistic) than the additive
expectation (p. 9, ll. 223–226). «Buffering combinations consistently achieved higher accuracy in absolute profile
correlation (PCC) compared to synergistic ones. In contrast, for Delta metrics, synergistic combinations outperformed
buffering combinations» (p. 9, ll. 228–231). The Discussion reads this as response magnitude rather than biology:
benchmarks «must separate response magnitude from interaction class before attributing performance differences to
biological non-additivity» (p. 12, ll. 342–344).

### 10. In Task 2 the model choice dominates metric variance (p. 9; Fig. 3F)

«model choice was the predominant determinant across nearly all metrics, overshadowing the contributions of seen
levels and additive interaction types» (p. 9, ll. 237–239). Fig. 3F gives the model share as 0.956 (MSE), 0.607
(Pearson), 0.701 (Pearson Δ), 0.956 (Δ agreement) and 0.922 (Top-50 DE). The seen-level shares are 0.042 / 0.264 /
0.189 / 0.038 / 0.064 (p. 47).

### 11. ⭐⭐ Cell-type transfer depends on baseline distance and on whether the response is conserved (pp. 9–10; Fig. 4B–D, p. 48)

Task 3 evaluates six methods (scPreGAN, scGen, trVAE, scVIDR, scPRAM, Squidiff) plus baselines (p. 9, ll. 245–248). On
absolute metrics all except scPreGAN score similarly high (p. 9, ll. 249–251). Datasets with smaller E-distances transfer
better (p. 10, ll. 252–258). «But E-distance is not fully predictive»: Burkhardt transfers poorly despite a small
E-distance, and Hagai transfers well despite a large domain shift, because of «conserved immune activation signatures»
(p. 10, ll. 259–263). «Task 3 performance is shaped not only by baseline domain differences between cell types, but
also by the consistency of perturbation effects across cell types» (p. 10, ll. 263–265). *Ambiguity:* E-distance is
defined as perturbed-vs-control separability (p. 10, ll. 252–254), but it is used as a source-vs-target distance. Fig. 4C
plots «the average E-distance among cell types for each dataset» (legend, p. 43).

### 12. No method generalises consistently across cell types; DE recovery is poor for all (p. 2; p. 10; Fig. 4E)

trVAE leads direction accuracy, Squidiff leads delta Pearson (p. 10, ll. 266–268). «At the DE level, all models
performed poorly» (p. 10, ll. 268–271). The abstract states «no method consistently generalized across cell types in
Task 3, where biological consistency dominated outcomes» (p. 2). The text gives **no** explicit comparison of Task 3
models against the Task 3 baselines. In Fig. 4D (RMSE, Kang/Hagai/Perturb_KHP), the Baseline3_ContextMean,
Baseline3_PerturbMean and Baseline3_Linear bars are often as low as or lower than the deep models (≈, read from figure;
for example K562 in Perturb_KHP).

### 13. In Task 3 algorithm design drives global metrics and target cell identity drives DE accuracy (p. 10; Fig. 4F)

«Global alignment metrics (MSE and Delta Pearson) were predominantly determined by algorithmic design factors,
specifically loss formulation and latent space structure. Conversely, the accuracy of differentially expressed gene
identification (Top 50 DE Accuracy) was governed by target cell identity» (p. 10, ll. 273–277). Fig. 4F shares (methods
/ perturbation / cell type / dataset) are MSE 0.57 / 0.13 / 0.23 / 0.08, Pearson 0.25 / 0.61 / 0.04 / 0.11, Pearson Δ
0.59 / 0.19 / 0.13 / 0.09, Δ agreement 0.32 / 0.45 / 0.07 / 0.16, and Top-50 DE 0.20 / 0.36 / 0.37 / 0.07 (p. 48).

### 14. ⭐ A single global score cannot summarise perturbation prediction (p. 10; p. 11)

«perturbation prediction performance cannot be summarized by a single global score; each metric must be interpreted
in terms of the biological signal it is able to measure» (pp. 10–11, ll. 285–288). Under weak perturbations, «a model can
appear accurate by reconstructing the basal expression state while failing to recover the perturbation-specific
response» (p. 11, ll. 292–294).

### 15. ⭐⭐ A diagnostic signal decomposition explains why the metric families disagree (p. 11; pp. 30–31; Fig. 5A–C, p. 49)

The observed pseudobulk profile is written as **Y_{c,p} = B_c + G_c + S_p + I_{p1,p2} + C_{c,p} + ε** (p. 30, l. 962).
The verbatim definitions (pp. 30–31, ll. 963–968) are:
B_c «the baseline control expression of cell state c»; G_c «the systematic shift shared across perturbations within the
same cell state»; S_p «the perturbation-specific response»; I_{p1,p2} «the non-additive interaction term for
combinatorial perturbations»; C_{c,p} «the cell-state-specific residual of the perturbation response, corresponding to
the cell-state-by-perturbation interaction term»; ε «the empirical noise floor».

The estimation (p. 31, ll. 969–980) goes as follows. B_c comes from control cells. G_c is «the average change within
each cell state». The rest is S_p. For combinations, I is the deviation from the additive expectation. For transfer, the
residual is split into «a perturbation-average component shared across cell states and a cell-state-specific residual
C_{c,p}». Inapplicable terms are set to zero. «The empirical noise floor ε was estimated from control cells using a
split-half procedure». Energies are sums of squares over entries and genes, reported in an uncentred "energy view" and a
"centered view" (after subtracting the gene-wise grand mean). Because the decomposition is «a component-wise diagnostic
rather than a strictly orthogonal ANOVA partition, normalized fractions were not constrained to sum to 100%» (p. 31,
ll. 987–989). The authors call it «a conceptual evaluative scaffold, not … a mechanistic linear model of gene regulation»
(p. 11, ll. 309–310).

The finding: «raw expression-level metrics are dominated by B_c, whereas centered, delta-based, DE-based, and
distributional metrics expose additional contributions from perturbation-specific, interaction-specific, and
context-dependent components» (p. 11, ll. 312–315).

Fig. 5B (≈, read from figure, p. 49). In the energy view B_c is ≈ 90–100 % in all tasks, and G_c is ≈ 9 % in Task 3.
In the centered view, Task 1 has S_p ≈ 100 % and ε ≈ 55 %. Task 2 has S_p ≈ 100 %, I ≈ 210 % and ε ≈ 5 %. Task 3 has
B_c ≈ 90 %, G_c ≈ 14 %, and thin S_p, C_{c,p} and ε layers (each under ≈ 10 %; exact values not legible).

### 16. ⭐ Transfer is a matter of "biological transportability" of the response mechanism (p. 12; p. 13)

«transferability depends jointly on basal-state similarity and conservation of the perturbation response program across
contexts … the challenge in Task 3 is not simply statistical alignment between cell types, but biological
transportability of the perturbation mechanism» (p. 12, ll. 349–355). Future cross-context benchmarks «should
deliberately sample four regimes: similar basal states with conserved responses, distant basal states with conserved
responses, similar basal states with divergent responses, and distant basal states with divergent responses» (p. 13,
ll. 391–394).

### 17. Three principles: target S_p, model I_{p1,p2}, and model whether responses are conserved (C_{c,p}) (pp. 12–13; Fig. 5D)

First, «models should be trained and evaluated against perturbation-specific effects (S_p) rather than raw profiles
alone». Second, combinatorial prediction «should explicitly target non-additive interaction residuals». Third,
«cross-context prediction should move beyond basal-state alignment and explicitly model whether perturbation response
programs are conserved or rewired across cellular states (C_{c,p})» (pp. 12–13, ll. 357–367). Fig. 5D's take-home
messages restate these (p. 49). A method-selection decision framework is in fig. S27 (not in this PDF).

### 18. Dataset properties matter as much as model choice; data design is the bottleneck (p. 13)

«Our variance analysis indicates that intrinsic dataset properties can contribute as much to performance variation as
model choice, emphasizing the need for standardized perturbation datasets with matched controls, sufficient replicates,
and carefully measured perturbation efficiency» (p. 13, ll. 374–378). Other stated limits are destructive sequencing,
so models infer population-level counterfactuals, and sparse temporal sampling (p. 13, ll. 378–383). «the path toward
virtual-cell models is not simply to train larger foundation models» (p. 13, ll. 398–399).

### 19. Two new drug-response datasets in primary immune cells (pp. 14–15; Fig. 1A)

Perturb KHP has HDFs, PBMCs (two donors) and K562, with drugs for 3–4 days (pp. 14–15; drug list Table S3). Perturb cmo
V1 has PBMCs from three donors treated with 111 repurposing drugs at 1 µM for 4 days, multiplexed with CellPlex (p. 15).
The Fig. 1A legend gives «133 drug treatments across lymphoid immune subsets» (p. 41). In Task 3 each Task 3 drug
dataset is reduced to the TRADE top-5 and bottom-5 conditions ("sub10"; p. 16, ll. 493–497).

## Figures and tables

- Fig. 1 (p. 45): data overview (A), effect size vs cell number (B), the three tasks and four metric families (C).
- Fig. 2 (p. 46): Task 1. Absolute Pearson (B), effect-size dependence (C–D), distribution compression (E), DE overlap
  (F), **metric-variance decomposition (G)**.
- Fig. 3 (p. 47): Task 2. **Seen0/1/2 by model incl. baselines (B)**, buffering vs synergy (C–D), radar (E),
  **ANOVA shares (F)**.
- Fig. 4 (p. 48): Task 3. Cosine similarity (B), **mean E-distance per dataset (C)**, RMSE per target cell type incl.
  Baseline3 models (D), radar (E), **variance partition by model / perturbation / cell type / dataset (F)**.
- Fig. 5 (p. 49): **decomposition schematic (A), component energies per task (B)**, across-model mean vs variance (C),
  take-home messages (D), **comparison with Ahlmann-Eltze, Systema and scPerturBench (E)**.
- Table 1 (p. 44): all 25 datasets with original and subset sizes.

## What this source does not address

- **No noise-corrected components and no ceiling.** ε is a descriptive "noise floor" from control-cell split halves
  (p. 31). Components are not noise-corrected, not orthogonal, and are computed on the observed profiles of the common
  evaluation subset only (p. 30). Nothing turns them into an attainable-accuracy bound. See
  [[zero-shot-transfer-ceiling]] and [[molina-2026-response-decomposition]] (cross-half inner products).
- **No context-count analysis.** No result is reported as a function of the number of training cell types or contexts.
  Tahoe enters Task 3 only as a 10-drug subset of conditions with ≥ 40 cell lines per plate (p. 16). See
  [[context-count-learning-curve]].
- **No cross-cell-line genetic transfer.** Task 3 uses stimulus, development and drug datasets (Table 1, p. 44). The
  genetic datasets (Tasks 1–2) are each single-context within-line splits.
- **No Systema-style reference or correction in the metrics.** Pearson Δ subtracts the control (p. 25, ll. 776–782).
  The "centered view" appears only inside the decomposition (p. 31) and as a recommendation (p. 12). See
  [[systematic-variation]].
- **No control-state or PC analysis** of responses or residuals. See [[control-subspace-projection]].
- **No CRISPRi multi-line resource** such as [[replogle-nadig-crispri-collection]] beyond the single RPE1 screen (Table 1).

<!-- ═══ FOKUS ═══ -->

## Focus: facts and claims the manuscript cites; exactly what overlaps with H1–H4 — 2026-11-06

> [!warning] Provenance: based on the **preprint** (bioRxiv v3, 2026-08-03), not the published original
> (@li-2026-systematic-comparison). Pagination and final claims may differ.

**Carries Kernthesen 1, 8, 11, 12, 15, 16, 17, 18** (and 5, 14 for metric framing). The hypotheses are those of
`input/ideas/zero-shot-limits-context-complexity-plan.md`.

### What this focus needs beyond the Kernthesen

1. **Verification of the wiki's summary sentence** ([[related-work-positioning]] §3, logged search summary).

   | Element in the wiki sentence | Preprint | Verdict |
   |---|---|---|
   | "decomposes responses" | Decomposes the **observed post-perturbation profile** Y_{c,p} (expression levels, B_c included), not the response δ, as a diagnostic for metrics (pp. 11, 30) | partly right: it decomposes profiles, as a metric diagnostic |
   | "baseline" | B_c «baseline control expression of cell state c» (p. 30, l. 963) | supported |
   | "shared shift" | G_c «the systematic shift shared across perturbations within the same cell state» = mean change within the cell state (pp. 30–31) | supported |
   | "perturbation-specific" | S_p «the perturbation-specific response» (p. 30) | supported |
   | "interaction" | **Two** terms: I_{p1,p2} (combinatorial non-additivity, Task 2) and C_{c,p} (cell-state × perturbation, Task 3) (pp. 30–31) | ambiguous: must name both |
   | "residual with a noise floor" | ε is one term, «the empirical noise floor», «estimated from control cells using a split-half procedure» (p. 31, ll. 968, 979–980). Elsewhere it is «residual noise or unexplained variation» (p. 11, l. 309) | supported as wording, but not a replicate-based ceiling: no formula, no scaling, no noise correction of the components |
   | (implied quantitative partition) | «not a strictly orthogonal ANOVA partition … not constrained to sum to 100%» (p. 31) | not supported: there is no variance partition |

2. **Mapping to project notation (reader's mapping, not in the source).** B_c + G_c corresponds to μ + b_c, with B_c
   the control level and G_c the within-line mean shift. Li's G_c is Molina's α_c and Systema's shared shift. S_p (in
   Task 3, the «perturbation-average component shared across cell states») corresponds to a_t. C_{c,p} corresponds to
   i_{c,t}. ε corresponds to ε, but Li estimates it from **control** cells, not from perturbed-cell split halves.
   I_{p1,p2} has no counterpart in the project.
3. **Simple baselines.** The strong statement is for combinations: AdditiveMean and Linear are best at seen1/seen2
   (p. 9). Fig. 5E marks "Simple baselines are not consistently outperformed by deep or foundation models" as covered
   (√) for this study (p. 49). For Task 1 the text gives no explicit baseline-vs-model statement. For Task 3 there is
   none either; Fig. 4D shows baselines competitive on RMSE (≈, read from figure).
4. **Unseen-context protocol.** «the training set included perturbation response data from only a subset of cell
   types, while data from the target cell type were completely excluded» (p. 32, ll. 1029–1031). The baselines still
   use target control expression, e.g. PerturbMean x_{c_target,ctrl} + Δ_p averaged over other cell types (eqs. 8–9,
   p. 19). The target's controls are therefore available. This is the same zero-shot setting as the project. *Ambiguity:*
   Baseline3_ContextMean averages «the target cell type across all perturbations observed for that cell type in the
   training data» (p. 19, ll. 587–590), which conflicts with complete exclusion.
5. **Metrics.** Pearson Δ = Pearson after «subtracting the control values from both the predicted and actual data»
   (p. 25). There is no perturbed-centroid reference, and Systema appears only as a compared framework (Fig. 5E).
   Delta direction accuracy, DE overlap top-N (N ∈ 10…500), effect-size Spearman, MMD and Wasserstein are per-gene
   aggregated (pp. 25–27).

### Overlap with H1–H4

| Hyp. | What the source shows exactly | What we add | Overlap |
|---|---|---|---|
| **H1** (pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | No ceiling and no noise correction. «no method consistently generalized across cell types in Task 3» (p. 2). In the centered Task 3 decomposition C_{c,p} is a thin layer next to B_c and G_c (Fig. 5B, ≈). | A noise-corrected transfer ceiling on CRISPRi panels and the methods placed on it (H1 untestable as registered). | low (framing only) |
| **H2a** (Tahoe learned context correction negative at C = 12; C*_emp ≥ 24) | Tahoe is used only as a Task 3 10-drug subset with ≥ 40 cell lines per plate (p. 16). There is no C manipulation. | A matched-ablation context curve. | none |
| **H2b** (theory break-even within ×2) | No theory. | Prop. 2 tested (rejected). | none |
| **H2c** (genetic C*_min > C_public) | No count of contexts and no genetic cross-line task. Its own call is to sample basal-distance × response-conservation regimes (p. 13). | C*_min and C_public. | none (adjacent motivation) |
| **H3a** (transfer residual on destination control PCs) | B_c is only the control mean. There are no control PCs and no residual analysis. Its qualitative claim is that transfer needs more than «basal-state alignment» (pp. 12–13). | Noise-controlled residual projection with other-line and permuted nulls. | none (opposite emphasis: Li puts the missing part in response conservation C_{c,p}, not in the control state) |
| **H3b** (M3 vs five corrections incl. Systema-style removal) | G_c is a Systema-like shared shift, and the "centered view" and centered metrics are recommended (pp. 12, 31). No correction is compared. | Five-way correction comparison on PDS. | low |
| **H4** (same-line data beats the largest modelling gain) | «intrinsic dataset properties can contribute as much to performance variation as model choice» (p. 13). Data-design recommendations (p. 13). No same-line data-addition experiment. | Measured H1GW gain vs modelling gains. | low (qualitative support) |

### Direct quotes (for the draft)

> "We represented each observed post-perturbation profile using the diagnostic decomposition Y_{c,p} = B_c + G_c + S_p +
> I_{p1,p2} + C_{c,p} + ε" (p. 30, ll. 960–962)

> "C_{c,p} denotes the cell-state-specific residual of the perturbation response, corresponding to the
> cell-state-by-perturbation interaction term; and ε denotes the empirical noise floor." (pp. 30–31, ll. 966–968)

> "The empirical noise floor ε was estimated from control cells using a split-half procedure." (p. 31, ll. 979–980)

> "Because this decomposition was used as a component-wise diagnostic rather than a strictly orthogonal ANOVA
> partition, normalized fractions were not constrained to sum to 100%." (p. 31, ll. 987–989)

> "additive Mean and Linear models demonstrated superior performance in seen1 and seen2 scenarios, often surpassing
> more complex architectures in these settings where constituent perturbations were known" (p. 9, ll. 219–222)

> "the challenge in Task 3 is not simply statistical alignment between cell types, but biological transportability of
> the perturbation mechanism." (p. 12, ll. 354–355)

### Examples & illustrations (for later drafting)

- Burkhardt (small E-distance, poor transfer) vs Hagai (large shift, good transfer via conserved immune activation)
  (p. 10, ll. 259–263). This is the clearest case that basal distance alone does not set transfer difficulty.
- Task 2 AdditiveMean: delta Pearson ≈ 0.1–0.2 at seen0 vs ≈ 0.8–0.9 at seen2 (Fig. 3B, ≈).
- Task 1 centered decomposition: ε ≈ 55 % of centered energy next to S_p ≈ 100 % (Fig. 5B, ≈). Even a descriptive
  noise floor is large relative to the perturbation-specific signal.

### Boundary of this focus

This focus does not cover the model-by-model settings (pp. 19–22), the ranking meta-analysis (pp. 27–28), TRADE and
E-distance details (pp. 28–29) or the wet-lab protocols (pp. 14–15). It reads no supplementary figure or table.

## Preprint vs published version

- bioRxiv API (`api.biorxiv.org/details/biorxiv/10.1101/2024.12.23.630036`, queried 2026-11-06): three versions, v1
  2024-12-23, v2 2025-09-16, v3 2026-08-03. All are CC BY-NC-ND and all are linked to the published DOI
  10.1126/sciadv.aed3414. The API returns the same abstract text for all three versions, so v1/v2 differences cannot be
  told apart from this source and are **not audited**.
- Crossref (10.1126/sciadv.aed3414): published 2026-09-11 (the .bib note says 2026-09-09). The published abstract matches
  the v3 abstract except for copy-editing: "differential expression (DE) recovery" and a capital "Expression-level"
  after the colon. No substantive difference is found in the abstract.
- The body text, figures and supplementary materials of the published version are **not audited**. In particular, it is
  not known whether the decomposition section (pp. 30–31, Fig. 5) and the PerturbNet/SMILES inconsistency (Kernthese 6)
  are unchanged.

## Mentioned entities

- Concepts: [[perturbation-response-decomposition]], [[simple-perturbation-baselines]], [[systematic-variation]],
  [[zero-shot-transfer-ceiling]] (not addressed), [[context-count-learning-curve]] (not addressed)
- Datasets: [[tahoe-100m]] (Task 3 subset), [[replogle-nadig-crispri-collection]] (RPE1 screen only, Task 1)

## Connections

- Supports: [[simple-perturbation-baselines]] (p. 9), [[ahlmanneltze-2025-linear-baselines]] (p. 9, ref. 25).
- Parallel decomposition to [[perturbation-response-decomposition]] / [[molina-2026-response-decomposition]]. Li's is
  non-orthogonal, not noise-corrected and metric-oriented. Molina's is cross-half noise-corrected (inferred contrast).
- Cites: [[vinastorne-2026-systema]] (ref. 26), Wei et al. 2026 scPerturBench (ref. 31), Wenteler PertEval-scFM
  (ref. 22), Wu PerturBench (ref. 23), Kedzierska (ref. 24), Roohani VCC (ref. 27), Adduri State (ref. 28), Zhang
  Tahoe-100M (ref. 29), Mejia (ref. 30), Wong et al. 2025 (ref. 54) (pp. 35–40).
- Cited by: [[ahlmanneltze-2025-linear-baselines]] ref. 25 cites this exact preprint DOI (10.1101/2024.12.23.630036).
- Positioning: [[related-work-positioning]] §3 (sentence to be corrected by the coordinator; see focus item 1).
