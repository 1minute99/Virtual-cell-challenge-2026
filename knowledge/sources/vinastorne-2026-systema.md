---
title: "Viñas Torné et al. 2026 — Systema: evaluating perturbation response prediction beyond systematic variation"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: vinastorne-2026-systema
depth: standard
tags: [perturbation-prediction, benchmark, evaluation-metrics, systematic-variation, perturbed-reference, centroid-accuracy, crispri, unseen-perturbation]
relations:
  - target: systematic-variation
    type: supports
    confidence: extracted
    because: "Source of the concept and its measure: «We define perturbation-specific shifts as vectors that point to the centroid of cells that underwent the same type of perturbation using the centroid of control cells as reference, while the average perturbation effect is the vector pointing to the centroid of all perturbed cells» (p. 1054)."
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "«the two simple baselines performed comparatively or outperformed state-of-the-art methods across different datasets and evaluation metrics» (p. 1051; Fig. 1c, p. 1052)."
  - target: ahlmanneltze-2025-linear-baselines
    type: cites
    confidence: extracted
    because: "Cited as ref. 4 (bioRxiv 2024 preprint) for «considerably simpler approaches, including linear models» that «attain comparable predictive performance» (pp. 1050–1051, 1058)."
---

# Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation

## Bibliographic Details

Viñas Torné, Ramon, Maciej Wiatrak, Zoe Piran, Shuyang Fan, Liangze Jiang, Sarah A. Teichmann, Mor Nitzan, and Maria
Brbić. 2026. "Systema: a framework for evaluating genetic perturbation response prediction beyond systematic
variation." *Nature Biotechnology* 44 (6): 1050–1059. doi:10.1038/s41587-025-02777-8. Received 22 Nov 2024, accepted
14 Jul 2025, published online 25 Aug 2025 (p. 1050). CC BY 4.0 (p. 1059). Code: github.com/mlbio-epfl/systema,
Zenodo 10.5281/zenodo.15746994 (p. 1058; PDF p. 13).
Read from `<library>/pdf/vinastorne-2026-systema.pdf` (16 PDF pages). Citekey: @vinastorne-2026-systema.

**Page numbering.** PDF pp. 1–10 carry printed pp. 1050–1059 (printed = PDF + 1049); these are cited by printed page.
The Methods (PDF pp. 11–13) carry no printed page number and are cited as "PDF p. N". PDF pp. 14–16 have no text layer
(p. 14 holds a single image). The Supplementary Information (Supp. Figs. 1–37, Supp. Tables 1–6, Supp. Information A
on model details) is **not** in this PDF; claims that rest on it are marked "(Supp., not read)".

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; introduction (simple baselines vs state-of-the-art; what Systema is) | 1050–1051 |
| ✓ F1 | Results: Simple baselines perform comparably to existing perturbation response prediction approaches (Fig. 1, p. 1052) | 1051 |
| ✓ F1 | Results: Systematic differences between control and perturbed cells lead to high predictive scores (GSEA/AUCell on Adamson, Norman; cell cycle and CIN in Replogle RPE1/K562; Fig. 2, p. 1053) | 1051–1054 |
| ✓ F1 | Results: Standard reference-based metrics are susceptible to systematic variation (cosine measure; correlation with PearsonΔ; RMSE; cell-cycle downsampling; Fig. 3, p. 1054) | 1054–1055 |
| ✓ F1 | Results: Evaluating the prediction of perturbation-specific effects (perturbed-centroid reference; Fig. 4, p. 1056) | 1055 |
| ✓ F1 | Results: Dissecting the predictable perturbations of perturbation response prediction methods (centroid accuracy, CIN classification; Fig. 5, p. 1057) | 1055–1058 |
| ✓ F1 | Discussion (limits of reference-sensitive metrics; heterogeneous panels; downstream utility; optical/spatial screens) | 1058 |
| — | References 1–29 | 1058–1059 |
| ✓ F1 | Methods: Nonparametric baselines (perturbed mean, matching mean); Perturbation response prediction methods; Note about CPA | PDF p. 11 |
| ✓ F1 | Methods: Evaluation (setting/splits; differential expression profiles; reference-based metrics; reference-insensitive metrics; Systema perturbed-centroid reference, eq.) | PDF pp. 11–12 |
| ✓ F1 | Methods: Limitations of reference-sensitive metrics; Centroid accuracy (eq.) | PDF p. 12 |
| — | Methods: Gene set enrichment analysis; AUCell | PDF p. 12 |
| ✓ F1 | Methods: Dataset details (data processing, control types, CIN quantification, cell-cycle annotation); Data availability | PDF pp. 12–13 |
| — | Code availability, references 30–35, acknowledgements, author contributions, competing interests | PDF p. 13 |
| — | Supplementary Information (Supp. Figs. 1–37, Supp. Tables 1–6, Supp. Information A) — not in the PDF | — |

## Kernthesen

### 1. ⭐⭐ Two nonparametric "average effect" baselines match or beat CPA, GEARS and scGPT on standard metrics (p. 1051; Fig. 1c, p. 1052)

The authors benchmark CPA, GEARS and scGPT (zero-shot and fine-tuned) on unseen one-gene perturbations in ten
datasets against two baselines that «exclusively capture average treatment effects» (p. 1051): the **perturbed mean**
μ_pert, the average expression of all perturbed training cells, used as the prediction for every perturbation; and the
**matching mean** for combinations, ½(μ(X) + μ(Y)) with an unseen single replaced by μ_pert (PDF p. 11; Fig. 1a). With
control-referenced PearsonΔ (all genes) «the perturbed mean baseline outperformed other methods across all datasets»;
with PearsonΔ20 it was higher or comparable everywhere except Frangieh, where scGPT won (p. 1051). RMSE shows the same
trend (Supp. Fig. 3, not read). The perturbed mean reaches PearsonΔ 0.70 on Adamson, 0.55 on Replogle RPE1 and 0.32 on
Replogle K562 (Fig. 1c; full table in the focus block).

### 2. The matching mean wins on unseen two-gene combinations (p. 1051; Fig. 1d, p. 1052)

On Norman's two-gene perturbations, the matching mean outperformed every other method «by a considerable margin, with
relative improvements of 11% for PearsonΔ over the best alternative method (GEARS)» (p. 1051). Performance rises with
the number of constituent singles seen in training (0-, 1-, 2-seen: n = 39, 147, 39 over three splits; Fig. 1d
legend, p. 1052).

### 3. ⭐ The methods' predicted profiles look alike and track the perturbed mean (p. 1051; Fig. 1e)

Predicted differential expression profiles «were similar to each other and correlated with those of the perturbed
mean», e.g. for FOSB+OSR2 (Fig. 1e). The authors conclude that the methods «predominantly capture systematic
differences between control and perturbed cells and may fail to uncover perturbation-specific effects» (p. 1051).

### 4. ⭐⭐ Systematic variation: consistent perturbed-vs-control differences with design, confounder or biological causes (pp. 1050–1051)

Defined in the abstract as «the consistent transcriptional differences between perturbed and control cells arising
from selection biases or confounders» (p. 1050). Three origins are named (p. 1051): a perturbation panel of
functionally related genes (selection bias); unmeasured variables such as cell-cycle phase, chromatin landscape and
target efficiency; and responses that are «biological in origin but systematic in effect, such as stress response or
cell-cycle arrest», shared by many perturbations. Its danger: it «can lead to overestimated predictive performance of
perturbation response models when they primarily capture the average perturbation effect» (p. 1051).

### 5. Focused panels (Adamson, Norman) carry pathway-level systematic shifts (pp. 1051–1052; Fig. 2a–b, p. 1053)

GSEA (GO Biological Process 2023) and AUCell show that perturbed cells in Adamson (UPR/ER-homeostasis panel) differ
from controls in response to external stimuli, chemical stress and positive regulation of cell death; in Norman (cell
cycle / growth panel) cell-death activation, down-regulated stress-response pathways, cell cycle and erythrocyte
differentiation (pp. 1051–1052; Supp. Figs. 8–13, Supp. Table 1, not read).

### 6. ⭐ Even a genome-scale screen is not free of it: RPE1 shows a cell-cycle shift, K562 much less (pp. 1052–1054; Fig. 2c–f, p. 1053)

In Replogle RPE1, 46 % of perturbed cells vs 25 % of control cells are in G1 (Jensen–Shannon 0.16, χ² = 196.35,
P = 2.3 × 10⁻⁴³; p. 1053). The authors attribute this to perturbation-induced chromosomal instability (CIN) and the
G1 arrest of p53-positive RPE1 cells with abnormal karyotypes. In the genome-wide K562 screen the cell-cycle
difference is «considerably smaller», likely because K562 is p53-negative (p. 1053). GSEA shows negative enrichment of
DNA repair and replication in RPE1 and of ribosome biogenesis in K562 (pp. 1053–1054).

### 7. ⭐⭐ Systematic variation is measured as the mean cosine between each perturbation's shift and the average perturbation effect (p. 1054; Fig. 3a–b)

«We define perturbation-specific shifts as vectors that point to the centroid of cells that underwent the same type
of perturbation using the centroid of control cells as reference, while the average perturbation effect is the vector
pointing to the centroid of all perturbed cells» (p. 1054). The dataset score is the cosine similarity averaged over
perturbations (Fig. 3a). Adamson and Norman are high, Replogle RPE1 is higher than K562, and the three Frangieh
datasets are lowest (pp. 1054–1055; Fig. 3b: roughly 0.2–0.7, read from the plot). The choice of control type
(non-targeting vs intergenic guides) does not matter much (p. 1055; Supp. Fig. 20, not read).

### 8. ⭐⭐ Control-referenced Pearson scores track systematic variation almost perfectly (p. 1055; Fig. 3c)

Across the ten datasets, mean cosine (systematic variation) vs PearsonΔ correlates at 0.91 for fine-tuned scGPT and
0.95 for GEARS (p. 1055); the figure annotates an average correlation of 0.97 (Fig. 3c). PearsonΔ and PearsonΔ20 fall
monotonically when test perturbations are subsampled by how much systematic variation they carry (Supp. Fig. 22, not
read). RMSE and RMSE20 do not show this trend, because they are reference-insensitive (p. 1055). Down-sampling RPE1
perturbed cells to match the control cell-cycle distribution gave «significantly lower» PearsonΔ and PearsonΔ20
(p. 1055; Supp. Figs. 25–26, not read). Conclusion: high reference-based scores «may reflect dataset-specific
confounding factors rather than accurate modeling of perturbation biology» (p. 1055).

### 9. ⭐⭐ Systema moves the reference from the control centroid to the "perturbed centroid", the mean of perturbation centroids (p. 1055; PDF p. 11; Fig. 4a–b, p. 1056)

«Instead of using control cells as a point of reference, Systema allows using custom references that better isolate
perturbation-specific effects. As a key example, we propose shifting the reference to the average of
perturbation-specific centroids, which we refer to as perturbed centroid» (p. 1055). Formally
O_pert = (1/|𝒫_train|) Σ_{X∈𝒫_train} O(X), with O(X) the mean of all cells carrying perturbation X (PDF p. 11). The
same reference replaces the control mean in both the true and the predicted differential profile. With the control
reference, perturbation shifts and the average shift have mean cosine 0.41 ± 0.18 across datasets (Adamson
0.76 ± 0.30, Norman 0.50 ± 0.26, Replogle K562 0.32 ± 0.16). With the perturbed reference the mean cosine is −0.06
(Adamson 0.05 ± 0.23, Norman −0.12 ± 0.45) (p. 1055; Fig. 4b). The reference change has «a weak impact on the norms
and biological interpretation» of the shifts (p. 1055; Supp. Figs. 27–32, Supp. Tables 3–6, not read).

### 10. ⭐⭐ Under the perturbed reference every method scores low, and the link to systematic variation largely disappears (p. 1055; Fig. 4c, e, p. 1056)

«Applying Systema resulted in substantially lower evaluation scores» (p. 1055): all methods attain low correlations
on unseen one-gene perturbations, with scGPT «marginally higher» than the others (Fig. 4c). The correlation between
perturbed-reference PearsonΔ and systematic variation, averaged over methods, drops to 0.21 (p. 1055; Fig. 4e). The
authors read this as showing that «predicting perturbation-specific effects for unseen genetic perturbations is a
remarkably challenging task» (p. 1055).

### 11. Under the perturbed reference only GEARS and scGPT stay non-negative on fully unseen combinations (p. 1055; Fig. 4d)

For Norman's two-gene perturbations, «GEARS and scGPT were the only methods to achieve non-negative PearsonΔ scores
in the 0-seen perturbation setting»; for 1-seen and 2-seen, GEARS, scGPT and the matching mean are best (p. 1055).
Fisher-combined paired t-tests of the top two methods: 0-seen P = 0.002, 1-seen P = 0.0003, 2-seen P = 0.6 (Fig. 4
legend, p. 1056).

### 12. ⭐ Centroid accuracy: is the predicted centroid closer to its own true centroid than to the others? (pp. 1055–1057; PDF p. 12; Fig. 5a–c)

Centroid accuracy(X) = (1/(|𝒫| − 1)) Σ_{Y∈𝒫} 𝟙[d(O_pred(X), O(X)) < d(O_pred(X), O(Y))], with d the Euclidean
distance between ground-truth and predicted perturbation centroids (PDF p. 12). It is computed for test
perturbations; 1 means the prediction recovers the expected effect (p. 1055). Plotted per perturbation against the
perturbed mean's own centroid accuracy, points on the diagonal are «uninformative» because the perturbed mean
«captures systematic variation exclusively» (Fig. 5b legend, p. 1057). Across the ten datasets average accuracies
«barely exceeded those of the perturbed mean»; fine-tuned scGPT is highest overall and wins on Adamson, Norman and
Replogle K562 (p. 1055; Fig. 5c).

### 13. Fine-tuned scGPT partly recovers coherent gene groups (translation, ribosome) and CIN class (pp. 1055–1058; Fig. 5d–e)

On genome-wide K562, scGPT is especially accurate for perturbations of core processes (gene expression, translation,
DNA replication), including ribosomal-protein genes; the authors credit the «strong and stereotyped transcriptional
programs» of translation disruption and atlas pretraining (pp. 1055–1057; Fig. 5d). The class-centroid variant of
centroid accuracy classifies unseen perturbations as low CIN (z ≤ 0) vs high CIN (z > 2): only fine-tuned scGPT is
substantially above chance (ROC AUC 0.73; CPA 0.54, GEARS 0.53, scGPT 0.56) (pp. 1057–1058; Fig. 5e; PDF p. 12).

### 14. ⭐ Reference-sensitive metrics have three limitations; reference-insensitive ones avoid them but are harder to read (PDF pp. 11–12; p. 1058)

MSE between Δ_X and Δ̂_X «does not depend on the reference choice 𝒞 (the reference cancels out)», so it equals MSE on
absolute post-perturbation profiles (PDF p. 11). Reference-sensitive metrics (cosine, Pearson on shifts): (1) are
noisy for weak perturbations, whose shifts have small norm, so «the presence or absence of a few single cells might
have a large impact on the orientation»; (2) are often invariant to perturbation strength; (3) depend strongly on the
reference, and can be undefined when a centroid coincides with the reference (PDF p. 12). The Discussion repeats that
the perturbed reference «inherits» these limitations (p. 1058).

### 15. ⭐ Evaluation protocol: GEARS splits, 25 % of perturbations held out, three seeds, within one dataset (PDF pp. 11–12)

The task is generalisation to **unseen genetic perturbations** within each dataset. Following GEARS (also used by
scGPT), data are split by perturbation, with 25 % of perturbations never seen in training; for combinations the test
set splits into 0/2, 1/2 and 2/2 seen; the split is repeated with three random seeds (PDF p. 11). Processing: total-
count normalisation to the median, log transform, top 5,000 HVGs plus perturbed non-HVGs, and removal of
perturbations whose target is not in the gene panel (needed for scGPT). Replogle follows scGPT's processing: the
subset matching the 1,973 perturbations with strong transcriptional changes in the original study, 100 cells per
perturbation, 2,500 control cells. Frangieh follows Lopez et al. (≥ 500 genes per cell, genes in ≥ 500 cells,
single-gene guides only) and is split by condition into control, co-culture and IFN-γ datasets (PDF p. 12). Controls
are non-targeting guides; Frangieh adds intergenic guides (PDF p. 12). CPA has no mechanism for unseen genetic
perturbations and is run counterfactually from random control cells (PDF p. 11).

### 16. Evaluate predictions by downstream biological utility, and build more heterogeneous panels (p. 1058)

«larger datasets with highly heterogenous perturbation gene panels will be essential to improve the robustness and
generalizability of perturbation response models», although systematic variation «may persist even in large-scale
perturbation screens, such as Replogle RPE1» (p. 1058). Models should be judged by «their biological utility (how can
inferred perturbation profiles help us answer downstream queries about relevant cellular phenotypes?)»; optical
pooled screens and spatial functional genomics are named as richer read-outs (p. 1058).

## Figures and tables

- **Fig. 1** (p. 1052): (a) perturbed-mean / matching-mean schematic; (b) dataset ring chart (perturbation counts,
  technology, cell line); (c) control-referenced PearsonΔ / Δ20 heat-table, 5 methods × 10 datasets; (d) Norman
  two-gene PearsonΔ by number of seen singles; (e) FOSB+OSR2 predicted vs measured changes.
- **Fig. 2** (p. 1053): Norman UMAP + AUCell, GSEA NES; Replogle RPE1 UMAP by cell-cycle phase and CIN z-score;
  phase densities (G1 46 % vs 25 %); volcano plot of DNA-damage / cell-cycle genes.
- **Fig. 3** (p. 1054): (a) the cosine measure of systematic variation; (b) systematic variation vs number of
  perturbations, coloured by technology; (c) PearsonΔ vs systematic variation per dataset (avg. corr. 0.97).
- **Fig. 4** (p. 1056): (a) perturbed-centroid reference; (b) cosine distributions, control vs perturbed reference;
  (c) perturbed-reference PearsonΔ / Δ20 heat-table; (d) Norman two-gene under the perturbed reference; (e) PearsonΔ
  vs systematic variation under the perturbed reference.
- **Fig. 5** (p. 1057): (a–b) centroid accuracy and the "informative vs uninformative" diagonal; (c) centroid accuracy
  per dataset; (d) scGPT-ft vs perturbed-mean accuracy on Replogle K562 (ribosomal genes); (e) CIN ROC curves.
- No numbered tables in the main text; Supp. Tables 1–6 (GSEA results, cosine statistics) are not in the PDF.

## What this source does not address

- **No cross-cell-line or unseen-context evaluation.** Every split is by perturbation *within one dataset*
  (PDF p. 11); the five cell lines are never used to train on one and test on another. Zero-shot transfer to a new
  line is the subject of [[molina-2026-response-decomposition]] (which adapts the perturbed reference to that setting),
  [[piao-2026-adaptive-contexts]] and, for drugs, [[ji-2024-universal-prediction]]; the context-count axis is
  [[dibaeinia-2026-context-scale]] and [[context-count-learning-curve]].
- **No measurement-noise model or ceiling.** There is no split-half reliability, no noise-corrected correlation and no
  upper bound on achievable scores; weak-perturbation noise is named only as a metric limitation (PDF p. 12). Ceilings
  are treated in [[zero-shot-transfer-ceiling]] and [[molina-2026-response-decomposition]].
- **No correction of predictions.** Systema changes the *evaluation reference*; it does not subtract a common response
  from predictions or deflate any subspace, and it proposes no new predictor.
- **No variance decomposition.** The "average perturbation effect" is one vector per dataset; there is no
  cell-line × perturbation interaction (see [[perturbation-response-decomposition]]).
- **No chemical perturbations** and no Tahoe-100M; no linear model is benchmarked (linear models are cited, refs. 4–6,
  p. 1050). Linear baselines are [[ahlmanneltze-2025-linear-baselines]].
- **No per-dataset cell counts** in the main text or Methods (only the Replogle 100-cells-per-perturbation subsample
  and n = 324 cells for FOSB+OSR2, Fig. 1e); presumably in the Supplementary Information, which is not in the PDF.

<!-- ═══ FOKUS ═══ -->

## Focus: Data / benchmark facts to cite — datasets, cell lines, metrics, baselines and evaluation protocol; and what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 1, 4, 7, 8, 9, 10, 12, 14, 15** (and 6 for RPE1/K562, 16 for the outlook). Hypotheses are those
of `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3 at the end).

### What this focus needs beyond the Kernthesen

1. **The ten datasets** (six sources, three technologies, five cell lines; p. 1051, Fig. 1b, p. 1052; Fig. 3/4
   legends, pp. 1054, 1056; Fig. 5 legend, p. 1057; PDF p. 12). Cell line and technology are read from the Fig. 1b
   inner and middle rings; cell counts per dataset are not given (see "does not address").

   | Dataset | Ref. | Cell line (Fig. 1b) | Technology (Fig. 1b) | Perturbations after processing | Test perts, 3 splits (Figs. 1c/4c) | Unique test perts (Fig. 5c) |
   |---|---|---|---|---|---|---|
   | Adamson | 13 (Cell 2016; UPR / ER panel) | K562 | CRISPRi | 81 | 63 | 49 |
   | Norman | 14 (Science 2019; cell cycle / growth; singles + pairs) | K562 | CRISPRi wedge (see note) | 276 | 108 singles; pairs 39 / 147 / 39 (0/1/2-seen) | 85 |
   | Replogle K562 | 15 (Cell 2022), genome-wide screen | K562 | CRISPRi | 1,813 | 1,362 | 1,043 |
   | Replogle RPE1 | 15 | RPE1 | CRISPRi | 1,410 | 1,059 | 806 |
   | Tian CRISPRa | 16 (Neuron 2019; iPSC-derived neurons) | iPSC | CRISPRa | 97 | 75 | 61 |
   | Tian CRISPRi | 16 | iPSC | CRISPRi | 181 | 138 | 107 |
   | Xu | 17 (Nat. Biotechnol. 2024) | HEK293 | CRISPRi wedge | 198 | 150 | 116 |
   | Frangieh control | 18 (Nat. Genet. 2021; Perturb-CITE-seq) | Ma-Mel (melanoma) | CRISPR KO | 167 | 126 | 99 |
   | Frangieh co-culture | 18 | Ma-Mel | CRISPR KO | 167 | 126 | 99 |
   | Frangieh IFN-γ | 18 | Ma-Mel | CRISPR KO | 167 | 126 | 99 |

   Data access: GEO GSE90546 (Adamson), GSE146194 (Norman), GSE218566 (Xu); Replogle on figshare
   10.25452/figshare.plus.20022944; Frangieh on Broad SCP1064; Tian via scPerturb (Zenodo 13350497) (PDF pp. 12–13).
   Two figure/caption inconsistencies: the Fig. 1b ring puts the **181**-perturbation wedge in the CRISPRa colour,
   while the Fig. 3/4 legends give Tian CRISPRa n = 97 and CRISPRi n = 181 (the test counts 75 / 138 fit the
   legends). Fig. 1b also draws Norman inside the CRISPRi ring, although Norman et al. 2019 is a CRISPRa screen
   (reader's knowledge, not stated in this source). Cite modality from the primary papers, not from Fig. 1b.

2. **Baseline and method scores, control reference** (Fig. 1c, p. 1052; PearsonΔ / PearsonΔ20, mean of 3 splits;
   transcribed from the heat-table).

   | Method | Adamson | Norman | Repl. K562 | Repl. RPE1 | Tian a | Tian i | Xu | Fr. ctrl | Fr. co-cult. | Fr. IFN |
   |---|---|---|---|---|---|---|---|---|---|---|
   | CPA | .02/.11 | .21/.34 | .06/.07 | .10/.21 | .10/.17 | .24/.44 | .11/.23 | .08/.08 | .08/.20 | .07/.12 |
   | GEARS | .65/.72 | .41/.59 | .22/.31 | .48/.55 | .31/.56 | .24/.45 | .13/.21 | .03/.01 | .05/.07 | .05/.07 |
   | scGPT | .69/.77 | .40/.44 | .29/.46 | .52/.66 | .28/.55 | .37/.74 | .18/.28 | .07/.39 | .06/.44 | .05/.48 |
   | scGPT (ft) | .62/.79 | .36/.47 | .27/.43 | .51/.66 | .26/.62 | .24/.70 | .22/.44 | .04/.48 | .08/.57 | .06/.48 |
   | Perturbed mean | .70/.73 | .49/.52 | .32/.41 | .55/.62 | .50/.70 | .51/.75 | .32/.49 | .15/−.14 | .15/.13 | .12/−.03 |

   **Perturbed (Systema) reference** (Fig. 4c, p. 1056):

   | Method | Adamson | Norman | Repl. K562 | Repl. RPE1 | Tian a | Tian i | Xu | Fr. ctrl | Fr. co-cult. | Fr. IFN |
   |---|---|---|---|---|---|---|---|---|---|---|
   | CPA | .01/−.31 | .15/−.17 | .05/−.21 | .08/−.32 | .09/−.17 | .00/−.10 | .02/−.04 | −.01/.06 | −.01/.20 | −.02/.11 |
   | GEARS | .13/.03 | .16/.07 | .00/.00 | .19/−.00 | −.05/.13 | −.00/.05 | .01/−.04 | .01/.08 | .05/.08 | .03/.13 |
   | scGPT | .13/.36 | .14/−.01 | .03/.34 | .08/.21 | .08/−.05 | .03/.21 | .13/.13 | .08/.48 | .11/.44 | .11/.54 |
   | scGPT (ft) | .16/.51 | .11/.10 | .06/.32 | .13/.21 | .03/.05 | .01/.25 | .12/.28 | .03/.51 | .06/.58 | .05/.52 |
   | Perturbed mean | .00/.09 | .10/−.17 | .06/−.12 | .08/−.31 | −.11/.17 | .01/.10 | .15/−.09 | .03/−.17 | .10/.12 | .03/−.06 |

   The perturbed mean is not exactly zero under the perturbed reference because μ_pert is a cell-weighted mean of
   training cells while O_pert is a perturbation-weighted mean of training centroids (PDF p. 11) — a reader's
   explanation; the source does not comment on these non-zero values.

3. **The exact Systema operation, for H3b's "Systema-style common-response removal".** Verbatim definitions (PDF
   p. 11): Δ_X = (1/|𝒫(X)|) Σ_{i∈𝒫(X)} x_i − (1/|𝒞|) Σ_{i∈𝒞} x_i, with 𝒞 «the set of reference cells (for example,
   control cells)»; for a method, the first term is replaced by the mean of its predicted cells, Δ̂_X. Systema sets
   the reference to **O_pert = (1/|𝒫_train|) Σ_{X∈𝒫_train} O(X), O(X) = (1/|𝒫(X)|) Σ_{i∈𝒫(X)} x_i**, «the average
   O_pert of all perturbation-specific centroids O(X)» (PDF p. 11; the PDF's glyph for O is a centroid symbol). Then
   PearsonΔ = Pearson(Δ_X, Δ̂_X) over all genes, and PearsonΔ20 over the top 20 DE genes «from both the ground truth
   Δ_X and inferred Δ̂_X» (PDF p. 11), averaged over held-out perturbations. Properties that a faithful
   implementation must keep:
   - (a) the reference is computed from **ground-truth training perturbations**, not from predictions or test targets;
   - (b) it is the **unweighted mean of per-perturbation centroids** (each perturbation counts once), not the mean of
     all perturbed cells (that is μ_pert, the perturbed-mean baseline);
   - (c) the **same** vector is subtracted from truth and prediction, so the operation is a change of origin inside
     the metric. It is a full (η = 1) rank-one translation, not a subspace projection, and there is no shrinkage
     parameter;
   - (d) for any reference-insensitive metric (MSE, and distances between absolute centroids such as centroid
     accuracy) it changes nothing, because «the reference cancels out» (PDF p. 11).
   Consequence for the plan: Task 17's arm («Systema = weak-half subtraction of η × target-mean (override entries
   excluded), best η per panel», log 2026-11-05) is a **prediction-side adaptation**, not Systema itself. It departs on
   (a) (no destination training perturbations exist in zero-shot; [[molina-2026-response-decomposition]] likewise
   substitutes «the observed target perturbation centroid as an evaluation-only reference», p. 15 there), on (c) (it is
   subtracted from predictions only, with η ∈ {0.25 … 1}), and possibly on (b) (if "target-mean" is cell-weighted). The
   paper should call it "a prediction-side analogue of Systema's perturbed-centroid reference", cite this source for
   the reference (PDF p. 11) and not for the correction.
4. **The separate "systematic variation" statistic uses a different mean.** Its "average perturbation effect" is the
   vector from the control centroid to «the centroid of all perturbed cells» (p. 1054; Fig. 3a legend, p. 1054), i.e.
   cell-weighted, whereas the Systema reference O_pert is perturbation-weighted. Fig. 4b's legend calls it «the
   average perturbation-specific shift» (p. 1056), which is ambiguous between the two. Cite the p. 1054 definition.
5. **C_public (H2c, definition 5; Amendment A2).** Of the five cell lines here, only **K562** (Replogle genome-wide)
   and **RPE1** (Replogle; essential-scale in the original study, not called so here) have genome- or essential-scale
   CRISPRi; both are already in the Replogle–Nadig family ([[replogle-nadig-crispri-collection]]). Adamson and Norman
   are K562 again. Tian (iPSC-derived neurons, 97/181 perturbations), Xu (HEK293, 198) and Frangieh (Ma-Mel melanoma,
   CRISPR KO, 167) are targeted panels below 300 perturbations, so they fail definition 5 and A2 Count 2's
   ≥ 300-target rule; Tian CRISPRa is excluded by modality anyway. This source therefore adds **no** new cell type to
   C_public = 11 (Count 1). Whether HEK293 or iPSC-neuron genome-scale screens exist elsewhere is not answered here.
6. **Which of the project's metrics are reference-sensitive (reader's mapping onto Kernthese 14).** H2a's residual MSE
   is reference-insensitive in Systema's sense and so is immune to the inflation shown in Fig. 3c. Pearson-type
   metrics on control-referenced deltas (H1's noise-corrected correlation) are reference-sensitive. The VCC PDS used in
   H3b/H4 is a centroid-accuracy-type discrimination score; Systema's centroid accuracy is computed on absolute
   centroids and is reference-free (PDF p. 12). A common-response subtraction applied to predictions only (Task 17)
   moves the predicted centroids relative to the fixed true ones and so *does* change such a score; a Systema-style
   change of reference applied to both sides would not.
7. **Literaturguide check.** The pre-read note («H3b H4: systematic variation; Systema-style common-response removal
   is an H3b alternative») is right that the source defines the centring H3b borrows. Two corrections: Systema
   performs no "removal" on predictions (focus item 3), and its link to H4 is only the outlook on heterogeneous panels
   and downstream utility (Kernthese 16). There is no same-line-data vs modelling comparison.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partial (framing only)** | The source shows that control-referenced Pearson is inflated by the shared perturbed-vs-control shift: a perturbation-agnostic perturbed mean scores PearsonΔ 0.32–0.70 on the CRISPRi datasets, and scores correlate 0.91–0.95 with systematic variation. This supports reporting H1 on a reference-explicit metric and setting any ceiling against what a template-only predictor already reaches. There is no ceiling, no noise correction, no r_max, no attained fraction and no unseen-line setting. | p. 1051, p. 1055, Fig. 1c p. 1052, Fig. 3c p. 1054 |
| **H2a** (Tahoe: context gain negative at C = 12; C*_emp ≥ 24) | **none** | No context axis, no chemical data, no learned context correction. The only adjacent point: residual MSE, the H2a metric, is reference-insensitive in Systema's sense, so it is not exposed to the inflation documented here. Nothing here bears on C*_emp ≈ 37.4. | PDF p. 11 (MSE invariance) |
| **H2b** (linear-Gaussian C* within ×2 of C*_emp) | **none** | No sample-complexity theory. Nothing bears on the rejected C* formula (C*_theory 230 vs 37.4). | — |
| **H2c** (genetic C*_min > C_public) | **partial (data facts only)** | Five cell lines, of which only K562 and RPE1 have genome/essential-scale CRISPRi (both already counted). The source adds no cell type to C_public = 11 (Count 1) or to the donor-line count (Count 2). Its outlook, that larger, more heterogeneous perturbation panels are needed, concerns perturbation diversity, not cell-line count. | Fig. 1b p. 1052; PDF p. 12; p. 1058 |
| **H3a** (transfer residual on destination control PCs > other lines' / permuted PCs) | **none (adjacent)** | No residual or control-PC analysis. Adjacent: part of the systematic shift in RPE1 is a cell-cycle occupancy change (G1 46 % vs 25 %) and cell-cycle-matched downsampling lowers scores. A state axis that also varies among controls thus carries common-response signal, which is consistent with control-state axes mattering. This is not tested on a residual and not across lines. | pp. 1053, 1055; Fig. 2e |
| **H3b** (M3 beats five alternatives incl. Systema-style removal on PDS) | **partial (source of one alternative)** | Defines the perturbed-centroid reference O_pert that the Systema-style arm imitates (focus item 3), and warns that perturbed-reference metrics are noisy for weak perturbations, which matters because H3b's arms act on the weak half (M2 rank). It does not compare corrections and does not apply O_pert to predictions. The arm is the project's adaptation. | p. 1055; PDF pp. 11–12; Fig. 4 p. 1056 |
| **H4** (same-line genome-wide data beat the largest modelling gain) | **partial (qualitative)** | The source's broader message is that data design dominates modelling: average-effect baselines rival deep models, and the authors call for «larger datasets with highly heterogenous perturbation gene panels». It does not test adding same-line data, has no w-sweep and no VCC composite. | p. 1051, p. 1058 |

**Positioning note.** Nothing in the source touches C*_emp ≈ 37 or the rejected linear-Gaussian C* (H2a/H2b). For
H3b it changes how the Systema-style arm must be *described*: it is a prediction-side analogue, not Systema, and its
win on the H1 panel offline (log: Systema-best H1 −0.0076) should not be reported as "Systema beats M3".

### Direct quotes (for the draft)

> «Here, we show that current methods struggle to generalize beyond systematic variation, the consistent
> transcriptional differences between perturbed and control cells arising from selection biases or confounders.»
> (p. 1050)

> «For example, the Pearson correlation between the average amount of systematic variation per dataset (mean cosine
> similarity) and PearsonΔ was 0.91 for fine-tuned scGPT and 0.95 for GEARS. This suggests that these commonly used
> metrics can be influenced by systematic variation and may lead to inflated performance» (p. 1055)

> «As a key example, we propose shifting the reference to the average of perturbation-specific centroids, which we
> refer to as perturbed centroid» (p. 1055)

> «Reference-insensitive metrics like RMSE may offer more robust evaluation in the presence of systematic variation,
> as they do not rely on the control population as a reference point.» (p. 1055)

> «First, they might be inappropriate for evaluating weak perturbations because the norm of their shifts (with
> respect to the reference) is small, leading to noisy targets.» (PDF p. 12)

### Examples & illustrations (for later drafting)

- **FOSB+OSR2** (Norman, only OSR2 seen in training): every method's predicted change resembles the perturbed mean
  (Fig. 1e, p. 1052; n = 324 cells).
- **Replogle RPE1 cell cycle**: 46 % vs 25 % of cells in G1, perturbed vs control, attributed to CIN-induced arrest
  in p53-positive RPE1. K562 (p53-negative) shows much less (p. 1053; Fig. 2e).
- **Ribosomal-protein perturbations in genome-wide K562** are the ones fine-tuned scGPT predicts better than the
  perturbed mean (Fig. 5d, p. 1057).

### Boundary of this focus

This focus takes datasets, metrics, baselines and the exact centring operation. It does not take the GSEA/AUCell
pathway details (Supp. Figs. 8–19, not in the PDF) or the per-model architecture details (Supp. Information A, not
in the PDF). No result here is a cross-cell-line result.

## Mentioned entities

- Datasets: [[replogle-nadig-crispri-collection]] (Replogle K562 genome-wide, RPE1); Adamson, Norman, Tian, Xu,
  Frangieh (no pages)
- Concepts: [[systematic-variation]], [[simple-perturbation-baselines]], [[perturbation-response-decomposition]],
  [[zero-shot-transfer-ceiling]]
- Methods: CPA, GEARS, scGPT (no pages)

## Connections

- Supports: [[systematic-variation]] (source of the concept, the cosine measure and the perturbed-centroid reference).
- Supports: [[simple-perturbation-baselines]]. Perturbed mean and matching mean match or beat GEARS, scGPT and CPA
  (Fig. 1c).
- Cites: [[ahlmanneltze-2025-linear-baselines]] (ref. 4, preprint).
- Built on by: [[molina-2026-response-decomposition]] (perturbed-reference Pearson, adapted to zero-shot with the
  target centroid) and [[piao-2026-adaptive-contexts]] (uses the Systema protocol).
- Contrast: [[control-subspace-projection]]. That removes a subspace learned from controls; Systema translates by one
  centroid inside the metric.
- Referenced in: `input/ideas/zero-shot-limits-context-complexity-plan.md` (H3b alternative; Task 17).
