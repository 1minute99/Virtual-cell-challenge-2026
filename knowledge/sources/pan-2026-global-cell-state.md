---
title: "Pan et al. 2026 — Global cell-state and gene-program representations (SCMG; genome-scale CRISPRi Perturb-seq in H1 hESCs)"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: pan-2026-global-cell-state
depth: standard
tags: [dataset, perturb-seq, crispri, hesc, h1, few-shot, mixture-of-experts, cross-context, gene-programs, conserved-response, gse295214]
relations:
  - target: h1-hesc-crispri-perturb-seq
    type: supports
    confidence: extracted
    because: "Generates and releases the dataset the project calls H1GW: «we performed a genome-scale CRISPRi Perturb-seq experiment in human embryonic stem cells (hESCs)» with 6,638 sgRNAs / 2,978 genes and «approximately 1.3 million perturbed cells» (p. 14); «accession number: GSE295214» (p. 20)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "Uses the Replogle 2022 K562 and RPE1 genome-scale CRISPRi screens as its comparison contexts: «genome-scale CRISPR interference (CRISPRi) knockdown Perturb-seq datasets of K562 and RPE1 cells57» (p. 12), and holds them out one at a time for few-shot evaluation (p. 11)."
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "Reproduces the train-mean baseline and reports that its own zero-shot model loses to it: «while zero-shot prediction performance is worse than the baseline, the few-shot prediction performance surpassed the baseline» (p. 11); the baseline is the one «previously used to shown that previously existing deep-learning approaches do not output [sic] simple linear baselines in prediction accuracy9» (p. 11)."
  - target: ahlmanneltze-2025-linear-baselines
    type: cites
    confidence: extracted
    because: "Ref. 9 is the motivation for the whole paper: «recent benchmarking studies indicate that complex deep-learning models did not outperform simple linear baselines in prediction accuracy9» (p. 2), and supplies the baseline used in Figs. 4D,E and 6G (p. 11)."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "Implements a different factorisation of the same object — a context-independent functional embedding of the perturbed gene times a cell-state-dependent linear map (mixture of experts, p. 10; Methods pp. 57–58) — i.e. a conserved factor plus context-specific rewiring, but it never quantifies variance shares."
  - target: context-count-learning-curve
    type: supports
    confidence: inferred
    because: "Supplies the within-context analogue (accuracy vs number of measured target perturbations K, plateau at K ≈ 100; pp. 11, 15) and states the context axis is blocked by data: «publicly available Perturb-seq collections currently cover a limited number of cell-type contexts» (p. 10). It runs no sweep over the number of contexts."
---

# Global cell-state and gene-program representations reveal conserved and context-specific perturbation responses of cells

## Bibliographic Details

Pan, Xingjie, Reuben Saunders, Joseph M. Replogle, Jonathan Weissman, and Xiaowei Zhuang. 2026. "Global cell-state and
gene-program representations reveal conserved and context-specific perturbation responses of cells." bioRxiv preprint,
posted 18 May 2026 (HHMI / Harvard, Whitehead Institute, MIT). doi:10.64898/2026.05.16.725005. CC BY 4.0.
Read from `<library>/pdf/pan-2026-global-cell-state.pdf` (62 pp.; the printed page footers match the PDF page numbers,
so "p. N" below is both). Citekey: @pan-2026-global-cell-state.

**Data and code.** hESC Perturb-seq raw data: NCBI GEO **GSE295214**; down-sampled representation dataset and
standardised pseudo-bulk perturbation database on HuggingFace (`xingjiepan/SCMG_data`); browser viewer
`xingjiepan.github.io/SCMG`; code `github.com/xingjiepan/SCMG` and `SCMG_scripts` (p. 20).

**Identity check.** This is the source paper for the dataset the project registry calls **H1GW**
(`families.json`: "H1GW — H1 hESC — GSE295214 (Weissman lab, 10x 5′) — genome-scale (2,747 targets passing QC)").
It is *not* the VCC 2025 H1 training data (family `H1`, 300 targets).

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; Introduction | 1–3 |
| — | Results: Construction of a global manifold of cell states (SCMG architecture, 36 datasets, Fig. 1) | 4–5 |
| ✓ F1 | Results: Understanding cell states and state transitions with the global manifold (zero-shot integration benchmark vs 5 foundation models, decoding, projection; Fig. 2) | 5–6 |
| ◐ F1 | Results: Global identification of gene-expression programs (54 programs; Fig. 3) | 6–9 |
| ✓ F1 | Results: Global identification of perturbation classes (cell-state dependence of responses; mixture-of-experts; **few-shot prediction**, Fig. 4) | 9–11 |
| ✓ F1 | Results: Conservation patterns and cell-state-specific rewiring of perturbation effects (K562 vs RPE1 vs hESC ORF; Fig. 5) | 11–13 |
| ✓ F1 | Results: **Genome-scale Perturb-seq on human embryonic stem cells** (library, cells, knockdown, effect sizes, conservation, few-shot; Fig. 6) | 13–15 |
| — | Results: Perturbation-induced differentiation of embryonic stem cells (single-cell clusters, germ-layer vs mesenchymal transitions; Fig. 7) | 15–17 |
| ◐ F1 | Discussion (incl. data-sufficiency and experimental-design argument) | 17–19 |
| ✓ F1 | Data and code availability | 20 |
| — | References | 21–30 |
| ✓ F1 | Figure legends 1–7 | 31–42 |
| — | Extended Data Figures 1–6 | 43–49 |
| ✓ F1 | Supplementary table captions (Tables 1–6) | 50 |
| ✓ F1 | Methods: SCMG, dataset collection, pairwise integration, standard gene set, benchmarking | 51–55 |
| ✓ F1 | Methods: perturbation-dataset preprocessing, human/10x decoding, gene programs, program scores | 55–57 |
| ✓ F1 | Methods: mixture-of-experts, perturbation classes, **few-shot prediction** | 57–59 |
| ✓ F1 | Methods: hESC CRISPRi library, lentivirus, balancing screen, library design, Perturb-seq experiment, sgRNA calling, pseudo-bulk shifts | 59–61 |
| — | Methods: gene-set over-representation, single-cell clustering, enrichment of perturbations in clusters | 61–62 |

## Kernthesen

### 1. ⭐ The stated problem is that data abundance has not bought prediction accuracy (pp. 2–3)

«Notwithstanding the success in transfer-learning tasks, these models provide limited improvement in predicting
gene-expression responses to genetic perturbations; notably, recent benchmarking studies indicate that complex
deep-learning models did not outperform simple linear baselines in prediction accuracy9» (p. 2; ref. 9 =
@ahlmanneltze-2025-linear-baselines). The proposed remedy is representational, not architectural: a shared coordinate
system for cell states, state transitions and perturbations (pp. 2–3).

### 2. SCMG is a contrastive autoencoder over 36 datasets and ~12 M cells (pp. 4, 52–53; Fig. 1, p. 31)

An MLP encoder maps a transcriptome to a **512-dimensional** latent vector with no dataset-identity input; the decoder
is conditioned on dataset identity; the loss is reconstruction MSE + a triplet contrastive term + L2 (pp. 51–52).
39 datasets were collected, **36 used for training and 3 held out** for benchmarking (p. 52). The compendium «includes
12 million cells from 36 datasets» (Fig. 1 legend, p. 32). The standard gene set is **18,108** human/mouse-ortholog
genes (p. 53); the visualisation "representation dataset" is **133,061 cells, 797 author cell-type labels, 31 manually
grouped major cell types** (p. 54). Species identity is treated as a batch effect, which the authors flag «may
attenuate genuine species-specific differences» (p. 4).

### 3. Training supervision comes from 125 manually annotated dataset pairs (p. 53)

Five large datasets serve as integration references; «In total, we manually annotated 125 pairs of datasets
encompassing 4,918 matched cell-type pairs», yielding **4,917,573** cross-dataset mutual-kNN (K = 10) cell pairs and
**5,585,245** pairs including intra-dataset ones (p. 53). The manual annotation cost «scales quadratically with the
total number of datasets» (p. 53) — the admitted scaling limit of the method.

### 4. ⭐ SCMG beats five single-cell foundation models at zero-shot integration (pp. 5–6; Fig. 2A–B, p. 33; Methods p. 54)

Ten dataset pairs were co-embedded, **five of them including a dataset withheld from SCMG training**. Two metrics on a
k = 30 nearest-neighbour graph: fraction of **inter-dataset** edges (batch removal) and fraction of
**same-cell-type** edges (biology preservation), compared by two-sided Welch's t-test (p. 54). Against scVI,
Geneformer, scGPT, SCimilarity and UCE, «SCMG resulted in a substantially higher inter-dataset neighbor fraction
without significantly reducing the same-cell-type fraction» (p. 5), with «performance … comparable for dataset pairs
composed entirely of training datasets and those including an unseen dataset» (p. 5).

### 5. Decoded expression tracks measured expression at cell-type level but smooths single cells (p. 6; Fig. 2C–F)

On two unseen datasets (intestine, lung), per-gene cell-type-level Pearson correlations cluster near 1 while
single-cell correlations span a broad range; «SCMG attenuated stochastic fluctuations of gene expression, yielding
smoother expression profiles», which «may also obscure some genuine fine-scale variability» (p. 6). The
human/10x-equivalent decoding correlates with measured mean expression at **Pearson r = 0.936** (p. 7; Fig. 3B).

### 6. 54 gene-expression programs are derived jointly from the manifold and a 9-dataset perturbation database (pp. 7–9, 56–57; Fig. 3)

Nine RNA-seq-before/after-perturbation datasets (single-cell and bulk, human and mouse) covering perturbations of
**7,270 genes** form the perturbation database (p. 7). Physiological filtering: drop 2,376 low-expression genes, 6,451
high-entropy genes (entropy > 6.4) → 9,281; keep genes with ≥ 5 partners at r > 0.5 → **8,544**; PCA to 50 dims.
Perturbation side: **1,613** readout genes with |pseudobulk shift| ≥ 0.5 in ≥ 5 conditions, response vectors across
**11,406 perturbation conditions** (p. 56). Leiden (resolution 10) on the 100-dim joint space gives 187 clusters,
curated down to **54 programs of 8–448 genes**, split into "essential" and "differentiation" classes (pp. 8, 57).

### 7. ⭐⭐ Perturbation responses depend strongly on the initial cell state — measured directly (p. 9; Fig. 4A, p. 37)

Cosine similarities between expression-shift vectors induced by the **same** perturbation from **different** initial
cell states are plotted against the initial states' distance in the 512-dim latent space and against effect magnitude:
«Plotting these similarity scores against both the distances between initial cell states … showed a strong dependence
of the perturbation effects on the initial cell state» (p. 9). In Fig. 4A most perturbations sit between cosine 0 and
0.8, with the low-similarity points concentrated at small effect size and large cell-state distance. «This dependency
indicates that directly comparing shift vectors across contexts conflates intrinsic gene function with context-specific
response programs» (p. 9).

### 8. ⭐⭐ The model factorises response into a context-free gene embedding and a cell-state-dependent linear map (pp. 10, 57–58; Fig. 4B)

A mixture of experts: Δŷ = Σ_k π_k(z) W_k v_g, where v_g is the perturbed gene's **context-independent functional
embedding**, W_k are K expert matrices and the gating network π(z) is conditioned on the SCMG cell-state embedding z
(p. 58). Hyperparameters: embedding dimension **M = 16**, **K = 8 experts**, **N = 18,108 readout genes** (p. 58).
«Because the functional embedding is shared across all training data and thus context independent, this architecture
supports joint learning across multiple perturbation datasets — even when different experiments perturb different
subsets of genes — while explicitly modeling context dependence through the gating network» (p. 10).

### 9. ⭐ Training uses 7 Perturb-seq knockdown datasets, only 3 of them genome-scale, in K562 and RPE1 (p. 10)

«We trained the mixture-of-experts model on seven published Perturb-seq knockdown datasets57–59,61,62, three of which
contain genome-scale perturbations (on K562 and RPE1 cells) and the rest at a smaller scale» (p. 10). The admitted
bottleneck is contexts, not perturbations: «Although the model can, in principle, incorporate arbitrarily many
datasets, publicly available Perturb-seq collections currently cover a limited number of cell-type contexts» (p. 10).

### 10. ⭐⭐ Reproducibility gate: only 608 of the perturbed genes replicate between two K562 screens at r > 0.4 (pp. 10, 58)

«many perturbations produce only weak transcriptional shifts, yielding low signal-to-noise ratios and limited
reproducibility across biological replicates. To ensure that downstream clustering is driven by reproducible
perturbation signatures …, we restricted downstream analyses to 608 genes, the perturbations of which led to
pseudo-bulk expression-shift vectors that were concordant between two independent genome-scale Perturb-seq experiments
in K562 cells (Pearson correlation > 0.4)» (p. 10). The Methods name the two experiments as **day 6 vs day 8** K562
Perturb-seq (p. 58). This is a reproducibility *filter*, not a noise-ceiling estimate: no split-half variance
decomposition and no reliability-corrected bound appear anywhere in the paper.

### 11. ⭐⭐ Few-shot beats both the zero-shot model and the train-mean baseline; zero-shot does not (p. 11; Figs. 4D–G, p. 37)

Protocol (pp. 11, 58–59): retrain the MoE with the K562 (resp. RPE1) datasets **withheld**; choose K ∈ {5, 10, 20, 50,
100, 200} perturbed genes by K-means on the functional embeddings, at the cluster centroids, to cover the embedding
space; fit **ridge regression** from functional embedding to measured pseudobulk shift vector on those K; predict the
rest; score by **Pearson correlation between predicted and measured pseudobulk expression-shift vectors**, with the K
training perturbations held out of the evaluation. Comparators: (i) the **zero-shot** MoE, (ii) a **simple baseline
that predicts the shift of a gene knockdown as the mean of the shift vectors for the same perturbation across the
training datasets** (p. 11). Result: «our few-shot models substantially outperformed both this baseline model and the
Zero-shot model — even with as few as K = 10 training perturbations», accuracy «increased with larger K and plateaued
around K ~ 100», and «Few-shot prediction performed better in RPE1 cells» (p. 11). Read off the violins (medians,
approximate): RPE1 baseline ≈ 0.15, zero-shot ≈ 0.00, K = 5 ≈ 0.41, K ≥ 10 ≈ 0.44–0.46; K562 baseline ≈ 0.10,
zero-shot ≈ 0.04, K = 5 ≈ 0.10, K = 20 ≈ 0.16, K = 200 ≈ 0.21 (Figs. 4D–E, p. 37).

### 12. ⭐ Accuracy is driven by effect size, and the zero-shot/few-shot gap is explained as context dependence (p. 11; Figs. 4F–G)

«predictions were more accurate for genes with larger perturbation effects, consistent with improved signal-to-noise
ratios» (p. 11); Figs. 4F–G plot accuracy against measured perturbation effect size for K = 100 and the baseline. The
authors' explanation of their own zero-shot failure: «potentially because the functional pathway relationship is less
context-dependent and therefore easier to learn, whereas the specific perturbation effects is context-dependent and
hence more difficult to learn from the limited existing data» (p. 11). They concede «the absolute prediction accuracy
remains imperfect given the limited size of available training data» (p. 11).

### 13. ⭐ Essential-gene knockdown moves essential programs, not differentiation programs — in all three contexts (pp. 12, 14; Figs. 5A, 6E)

Applied to K562 CRISPRi, RPE1 CRISPRi and an hESC **ORF overexpression** screen (ref. 63): «the gene-expression
programs belonging to the differentiation category … were largely insulated from knockdown of essential genes, whereas
these knockdowns substantially impacted gene-expression programs belonging to "essential" category» (p. 12), while TF
overexpression does induce differentiation. The same insulation holds in the new hESC CRISPRi screen (p. 14).

### 14. ⭐⭐ Conserved vs rewired stress response, named pathway by pathway (pp. 12–13, 14)

Conserved across K562 and RPE1: lysosome/autophagy as a general stress response to nearly all essential knockdown
classes; UPR almost exclusively from ER-protein-processing perturbations; broad repression of cell-cycle,
spliceosome, proliferation and ribosome-biogenesis programs. Context-specific: p53 responds broadly in RPE1 but is
absent in p53-deficient K562; «mitochondrial transcription and translation perturbations activate ISR in K562 cells
but not in RPE1 cells, showing context-dependent rewiring of this pathway» — and this «cannot be attributed to
variation in CRISPRi efficiency» (p. 13). Ribosomal-protein genes are regulated oppositely in the two lines after
large-subunit ribosome-biogenesis knockdown (p. 13). In hESCs the same structure recurs with its own rewiring:
mitochondrial perturbations give the *largest* stress responses, not attenuated ones (p. 14).

### 15. ⭐⭐ The new dataset: genome-scale CRISPRi Perturb-seq in H1 hESCs, ~1.3 M cells, 2,978 genes (pp. 14, 59–61; Fig. 6A–C, p. 40)

Inducible **Zim3-dCas9** CRISPRi **H1 hESCs** (StemFlex, ROCKi passaging), library of **6,638 sgRNAs targeting 2,978
genes** (two sgRNAs per gene) plus **370 non-targeting controls**, selected from CRISPRiV2/V2.1 by a prior
genome-wide **balancing screen** (growth phenotype or altered OCT4/NANOG/SOX2 by FACS), with growth-dropout sgRNAs
deliberately over-represented in proportion to their expected dropout (pp. 59–60). Transduction at > 1000× coverage,
< 20 % rate; **8 days after dCas9 induction**, GFP-sorted, **10x 5′ scRNA-seq with direct sgRNA capture** (p. 60).
Yield: «approximately 1.3 million perturbed cells» (p. 14); «Each perturbation was represented by tens to hundreds of
cells» (p. 14; Fig. 6B peaks near 10^2.3 ≈ 200 cells per perturbation); «most sgRNAs achieved effective knockdown of
their target genes, with a median expression reduction of ~70%» (p. 14; Fig. 6C). Cell filter: ≥ 500 detected genes,
exactly one sgRNA with ≥ 20 barcode counts (cellranger 8.0.1, p. 60). Pseudobulk filter: ≥ 50 cells per perturbed
gene, target expression < 90 % of non-targeting mean, one-sided Mann–Whitney p < 0.01; a further filter on
inconsistent guide pairs (shift norm > 4 and between-guide r < 0.2) «removed 259 genes» (pp. 60–61).

### 16. ⭐⭐ hESCs are measurably more perturbation-resistant than K562, and this is not a knockdown artefact (p. 14; Fig. 6D)

«perturbation-induced transcriptional effect sizes were substantially smaller in hESCs than in K562 cells, as measured
by the distribution of the expression shift vector magnitude across perturbations … This attenuation cannot be
explained solely by differences in knockdown efficiency: even when restricting analysis to genes with high knockdown
efficiency in hESCs, the magnitudes … remained significantly lower than those observed in K562» (p. 14). Fig. 6D shows
both matched-efficiency subsets with p < 10⁻⁴ (hESC medians ≈ 3 vs K562 ≈ 4–5 in the figure's effect-size units).

### 17. ⭐ Few-shot transfers to hESCs too, from as few as K = 10 (p. 15; Fig. 6G–H)

«The few-shot prediction performance on hESCs was again substantial improved as compared to the zero-shot model and
the simple baseline models, starting from as few as ten shots … and the accuracy increased as K grew» (p. 15).
Medians read off Fig. 6G (approximate): baseline ≈ 0.05, zero-shot ≈ 0.03, K = 10 ≈ 0.09, K = 50 ≈ 0.17, K = 200
≈ 0.18 — the lowest of the three contexts. The authors' reading: «perturbation prediction may be more difficult for
cells with higher multipotency or pluripotency, potentially due to their stronger differentiation potential and more
diverse differentiation outcome» (p. 15).

### 18. hESC-specific biology: stress-induced mesenchymal transitions and germ-layer exits (pp. 15–17; Fig. 7)

Single-cell Leiden clusters split into stress states that keep pluripotency, germ-layer differentiation (6 enriched
perturbed genes, incl. POU5F1/SOX2/NANOG), mesenchymal transitions (**56** enriched perturbed genes, spanning DNA
repair, transcription, mRNA processing, translation, protein degradation) and damaged cells (pp. 16, 48). eIF2B–eIF2
components are the only eIFs enriched in the mesenchymal class, tying the transition to ISR (pp. 16–17).

### 19. ⭐ The design argument: sample cell-state space broadly, then sample the functional space sparsely (pp. 18–19)

«Rather than requiring dense perturbation data in every cell type of interest, our framework suggests that strategic,
sparse sampling of perturbations across the functional space may be sufficient to enable prediction of perturbation
effects of all genes in new contexts» (p. 18). And on contexts: «Currently, the accuracy of generalizable perturbation
prediction models is limited by the low coverage of starting cell states across the cell state space … perturbation
screens can be prioritized in representative cell states spanning major identity domains, followed by denser sampling
of the transitional regions … Model performance can then serve as a quantitative criterion for data sufficiency»
(p. 19). No number of contexts is ever named.

### 20. Framework-level summary claim (pp. 15, 17–19)

«some core perturbation-response patterns are conserved across highly distinct cell types, suggesting the
generalizability of these regulatory principles broadly across the cell-state landscape. In the meantime, we also
observed distinct context-specific rewiring of some stress response pathways» (p. 15). Conservation is asserted at the
level of **program-score patterns**, never as a fraction of response variance.

## Figures and tables

- **Fig. 4** (p. 37): **A** cosine similarity of same-perturbation shifts vs initial-state distance and effect size;
  **B** mixture-of-experts schematic; **C** Leiden clusters of functional embeddings; **D, E** prediction accuracy
  violins (baseline / zero-shot / K = 5…200) for **RPE1** and **K562**; **F, G** accuracy vs effect size, K = 100 vs
  baseline.
- **Fig. 6** (p. 40): **A** hESC screen schematic (6,638 sgRNAs, 1.3 M cells); **B** cells per perturbation (log10);
  **C** on-target knockdown (residual expression); **D** hESC vs K562 effect-size boxplots on both matched-KD subsets;
  **E** program shift scores per perturbation; **F** perturbation-class × program heat map; **G** few-shot violins for
  hESC; **H** accuracy vs effect size.
- **Fig. 5** (p. 39): program shift scores for RPE1 CRISPRi, K562 CRISPRi and the hESC ORF screen (A); essential-program
  heat maps per perturbation class for RPE1 (B) and K562 (C).
- **Supplementary Tables** (captions p. 50): T2 = the scRNA-seq datasets with training flags and cell counts; T5 = the
  nine perturbation datasets with "number of perturbations passed QC"; T1 = enriched perturbed genes per hESC class.
- No table in the paper reports prediction accuracies numerically; every accuracy quoted above is read off a violin or
  scatter plot.

## What this source does not address

- **No variance decomposition and no ceiling.** Conservation is shown as shared pathway-level program responses, never
  as a share of reproducible variance, and there is no noise-corrected ceiling, no split-half estimator and no
  reliability correction. The r > 0.4 K562 concordance gate (p. 10) is the only replicate statistic in the paper.
- **No context-count curve.** Nothing varies the number of training contexts. The K sweep is over **perturbations
  measured in the target context**, with the context held fixed.
- **No cell-line descriptor ablation.** The gating network is always on; there is no context-free-vs-context-aware
  comparison, and no learned-correction-vs-no-correction contrast.
- **No evaluation on the VCC task or metrics.** Scoring is per-perturbation Pearson on pseudobulk shift vectors; no
  DE overlap, no discrimination score, no composite, no single-cell distributional metric.
- **No drug perturbations, no combinatorial perturbations** (both named as future work, p. 19), and no
  knockout (Cas9-cutting) data — all genetic data are CRISPRi knockdown or ORF overexpression.

<!-- ═══ FOKUS ═══ -->

## Focus: Data / benchmark facts to cite — datasets, cell lines, metrics, baselines and evaluation protocol; and what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 1, 4, 7–12, 14–17, 19.** The hypotheses are those of the plan at
`input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 near the top, Amendments A1–A3 at the end);
C_public is operational definition 5, counted two ways under A2 (Count 1 = 11 cell types/lines, Count 2 = 25 donor
lines).

### What this focus needs beyond the Kernthesen

1. **The dataset facts, in one place (for the methods table).** Cell line **H1 hESC** (inducible Zim3-dCas9 CRISPRi);
   modality **CRISPRi knockdown**; **2,978 target genes × 2 sgRNAs = 6,638 sgRNAs** plus **370 non-targeting**;
   **≈ 1.3 M** perturbed cells after filtering; **tens to hundreds of cells per perturbation** (Fig. 6B mode ≈ 200);
   median on-target knockdown **≈ 70 %**; harvest **8 days** after dCas9 induction; chemistry **10x 5′ with direct
   sgRNA capture**; QC ≥ 500 genes/cell, one sgRNA at ≥ 20 barcode UMIs; pseudobulk gate ≥ 50 cells, KD < 90 % of NTC,
   one-sided MWU p < 0.01, plus a guide-consistency filter that **removed 259 genes**; accession **GEO GSE295214**
   (pp. 14, 20, 59–61). No per-cell UMI depth figure is reported for the hESC screen anywhere in the paper — only
   cells per perturbation — so the project cannot cite a depth number from here.
2. **Library selection is not unbiased genome-wide.** The 2,978 genes were pre-selected by a genome-wide **balancing
   screen** for growth phenotype *or* altered OCT4/NANOG/SOX2 expression, keeping genes with ≥ 2 sgRNAs showing
   evidence, and dropout-prone sgRNAs were deliberately over-represented in the pool (pp. 59–60). So the screen is
   "genome-scale" in the sense of *screened genome-wide, profiled on the phenotypic subset* — a selection effect that
   matters when H1GW is used as a source for an unselected target panel, and that also makes its per-target effect-size
   distribution non-comparable with an unselected genome-wide screen.
3. **Contexts actually used.** Perturb-seq contexts in the paper: **K562** and **RPE1** (Replogle 2022, genome-scale
   CRISPRi, 3 of the 7 training datasets), four smaller CRISPRi screens (Adamson UPR, Frangieh melanoma Perturb-CITE,
   Tian iPSC-neurons ×2, Jiang pathway signatures), the **hESC ORF overexpression** screen (Nakatake, ref. 63) as a
   control contrast, and the **new H1 hESC CRISPRi** screen. **Overlap with the project's harness: K562, RPE1, H1
   hESC.** No HepG2, Jurkat, HCT116 or KOLF2.1J.
4. **C_public (definition 5) is unchanged by this paper.** Its only new genetic screen is H1 hESC, which the registry
   already carries as family **H1GW** (`output/data-analysis/results/families.json`: "GSE295214 (Weissman lab, 10x 5′)
   — genome-scale (2,747 targets passing QC)", `counts_for_C_public: true`). Pan et al. is the **citation** for that
   row, not an addition to it: **Count 1 stays 11, Count 2 stays 25.** Note the project's two H1 families are distinct
   screens of the same line — `H1` (VCC 2025, 300 targets, not counted) and `H1GW` (this paper, counted); the
   `line_to_families` map already merges them under "H1 hESC", so no double count.
5. **Shared vs line-specific response — what can and cannot be cited.** The paper measures context dependence
   directly (Fig. 4A: cosine similarity of the same perturbation's shift from two initial states, falling with
   cell-state distance and with small effect size) and demonstrates conserved program-level responses across K562,
   RPE1 and hESC (pp. 12–15). But it reports **no share of variance** that is conserved vs line-specific, so it
   cannot be cited for a number next to Molina & Zhang's 29.4 % conserved / 23.5 % interaction. Its contribution to
   the project's a_t story is qualitative and pathway-level, plus the negative result that a context-conditioned
   model trained on 7 datasets predicts a held-out context **worse than the train-mean baseline** (p. 11).
6. **Baselines and metric, exactly.** Metric = Pearson correlation between predicted and measured **pseudobulk
   expression-shift vectors**, per perturbation, with the K calibration perturbations excluded from evaluation
   (pp. 11, 59). Baselines = (i) train-mean over the same perturbation in the training datasets (the
   train-mean baseline of @ahlmanneltze-2025-linear-baselines, cited here as ref. 9) and (ii) the paper's own **zero-shot** MoE.
   Few-shot model = **ridge regression** on K functional embeddings, K chosen by K-means centroids. There is no
   perturbed-reference / centroid subtraction (unlike [[systematic-variation]]), no DE-overlap metric and no
   discrimination score, so none of these numbers is comparable to a VCC composite.
7. **Replicate statistics available as noise-ceiling input (weak).** Only one: 608 genes pass r > 0.4 between the two
   genome-scale K562 screens (day 6 vs day 8), chosen out of the perturbations in the training screens (p. 10, 58).
   That is a between-timepoint concordance filter, not a split-half reliability, and the paper gives no distribution
   of the correlation, so it supports the *framing* of H1's noise ceiling but supplies no estimate for it.
8. **An hESC-specific caution for the project's H1 panel.** Effect sizes in hESCs are significantly smaller than in
   K562 at matched knockdown efficiency (Fig. 6D), and prediction accuracy is lowest of the three contexts
   (Fig. 6G vs 4D–E). If effect size is smaller at comparable depth, the achievable noise-corrected ceiling on an
   hESC destination is expected to be lower than on a somatic line — a prediction the project's own r_max estimates
   on H1 vs RN panels can check.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **not addressed (two inputs only)** | No ceiling, no variance components, no split-half. Two usable inputs: the 608/r > 0.4 K562 concordance gate as evidence that a large share of perturbations is irreproducible (p. 10), and the hESC-vs-K562 effect-size gap at matched knockdown (Fig. 6D), which predicts a lower r_max on the project's hESC destination. Conservation across K562/RPE1/hESC is pathway-level and qualitative, so it backs the *existence* of a conserved component a_t without quantifying it. | pp. 10, 14, 15; Figs. 4A, 6D |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **not addressed** | No drug data, no Tahoe, no sweep over number of training contexts, no gain-vs-C curve. The gating network is always on; there is no context-free arm to difference against. | pp. 10–11, 58 |
| **H2b** (theory predicts the break-even within ×2) | **not addressed** | No sample-complexity theory, no descriptor dimension, no noise-to-interaction ratio, no break-even of any kind. | — |
| **H2c** (genetic C*_min > C_public) | **not addressed; relevant to the C_public count and strongly consistent in direction** | Adds **no** new line/type to C_public — its hESC screen is already the registry's H1GW (Count 1 = 11, Count 2 = 25 unchanged; focus item 4). Direction of evidence: with **7 Perturb-seq datasets and effectively two genome-scale contexts (K562, RPE1)**, the context-conditioned model's **zero-shot prediction is worse than the train-mean baseline** in every held-out context (p. 11, Figs. 4D–E, 6G) — exactly the "below break-even" signature H2c asserts, though no C*_min and no count is computed. The paper names the same bottleneck in words: «publicly available Perturb-seq collections currently cover a limited number of cell-type contexts» (p. 10) and «the accuracy … is limited by the low coverage of starting cell states» (p. 19). | pp. 10, 11, 15, 19; Figs. 4D–E, 6G |
| **H3a** (transfer residual loads on the destination's control PCs) | **not addressed** | No prediction residual is analysed and controls are never used as a projection basis. The adjacent fact is that the **initial cell state** governs response similarity (Fig. 4A) and gates the expert mixture (p. 10) — i.e. the basal state carries context information — but the paper's own evidence that this helps is the gating network's *training* performance, never a residual decomposition. | p. 9, p. 10, Fig. 4A |
| **H3b** (M3 beats five corrections on PDS) | **not addressed (a common-response *model*, not a correction)** | The context-free functional embedding is a common-response factorisation in the spirit of Systema-style shared-variation removal, but it is a model component, not a post-hoc correction applied to predictions, and nothing is deflated, removed or compared. No PC deflation, no permuted-PC null, no PDS. | p. 10; Fig. 4B; pp. 57–58 |
| **H4** (adding H1GW, w = 1, beats the largest modelling gain) | **partly shown (the dataset itself, plus a strong qualitative analogue)** | This paper **is** the source of H1GW (GSE295214, p. 20) and supplies every methods fact the project must cite for it (focus item 1) and the selection caveat (item 2). The analogue: K measured target-context perturbations beat both the train-mean baseline and the zero-shot model from K = 10 and plateau at K ≈ 100 (medians ≈ 0.44 vs 0.15 baseline in RPE1; ≈ 0.18 vs 0.05 in hESC) — data in the target context beats the model by a wide margin. Caveats: a different metric (shift-vector Pearson, no composite, no scorer), K shots are *measurements of the target screen itself* rather than a second same-line screen, and no w-sweep or fusion weighting exists. | pp. 11, 15, 20, 59–61; Figs. 4D–E, 6G |

**Positioning consequence for the manuscript.** By this paper's own taxonomy, using H1GW as a source for an H1
destination is not cross-context transfer at all — it is the *same cell state*, a second screen. That is why the H4
comparison must be labelled as "more data in the destination context" and set against Pan's few-shot curve, not
against its zero-shot arm. Nothing here touches C*_emp ≈ 37.4 or the rejected linear-Gaussian C* formula: the paper
contains no context-count axis and no theory to compare the formula with.

### Direct quotes (for the draft)

> "Although the model can, in principle, incorporate arbitrarily many datasets, publicly available Perturb-seq
> collections currently cover a limited number of cell-type contexts." (p. 10)

> "many perturbations produce only weak transcriptional shifts, yielding low signal-to-noise ratios and limited
> reproducibility across biological replicates. … we restricted downstream analyses to 608 genes, the perturbations of
> which led to pseudo-bulk expression-shift vectors that were concordant between two independent genome-scale
> Perturb-seq experiments in K562 cells (Pearson correlation > 0.4)." (p. 10)

> "It is interesting to note that while zero-shot prediction performance is worse than the baseline, the few-shot
> prediction performance surpassed the baseline, potentially because the functional pathway relationship is less
> context-dependent and therefore easier to learn, whereas the specific perturbation effects is context-dependent and
> hence more difficult to learn from the limited existing data." (p. 11)

> "We designed a CRISPRi library comprising 6,638 sgRNAs targeting 2,978 genes selected for putative roles in hESC
> growth or pluripotency marker expression (Methods), and profiled approximately 1.3 million perturbed cells
> (Figure 6A). Each perturbation was represented by tens to hundreds of cells (Figure 6B), and most sgRNAs achieved
> effective knockdown of their target genes, with a median expression reduction of ~70%." (p. 14)

> "This attenuation cannot be explained solely by differences in knockdown efficiency: even when restricting analysis
> to genes with high knockdown efficiency in hESCs, the magnitudes of the perturbation-induced transcriptional shift
> remained significantly lower than those observed in K562." (p. 14)

> "Currently, the accuracy of generalizable perturbation prediction models is limited by the low coverage of starting
> cell states across the cell state space. … Model performance can then serve as a quantitative criterion for data
> sufficiency." (p. 19)

> "The genome-scale hESC perturb-seq data is available at NCBI GEO data repository (accession number: GSE295214)."
> (p. 20)

### Examples & illustrations (for later drafting)

- **Zero-shot below the mean baseline, few-shot far above it** (Figs. 4D–E, p. 37; Fig. 6G, p. 40): RPE1 medians
  ≈ 0.00 (zero-shot) vs ≈ 0.15 (train-mean) vs ≈ 0.44 (K = 10); hESC ≈ 0.03 vs ≈ 0.05 vs ≈ 0.09, rising to ≈ 0.18 at
  K = 50–200. The cleanest published instance of a context-conditioned model losing to a context-free mean at the
  number of contexts currently available.
- **Same perturbation, different starting state** (Fig. 4A, p. 37): shift-vector cosine similarities spread from
  ≈ 0 to ≈ 0.8, with the lowest values at the largest cell-state distances — a direct picture of the interaction term.
- **Pathway rewiring with a named mechanism check** (p. 13): mitochondrial transcription/translation knockdown
  activates ISR in K562 but not RPE1, and the authors exclude differential CRISPRi efficiency as the cause; p53 is
  absent in K562 because the line is p53-deficient. Useful as a concrete answer to "what does line-specific response
  look like biologically".
- **hESC resistance** (Fig. 6D, p. 40): effect sizes significantly lower than K562 on both matched-knockdown subsets
  (p < 10⁻⁴) — the illustration for a per-destination ceiling.

### Boundary of this focus

This focus leaves aside the SCMG manifold's own biology (Kernthesen 2–6: the 36-dataset atlas, the 54 gene programs,
the transdifferentiation excursion) except where the integration benchmark bears on the "foundation models do not
help" framing, and leaves aside the hESC differentiation results (Kernthese 18: germ-layer vs stress-induced
mesenchymal transitions, POU5F1/SOX2/NANOG). A reader looking for a context-count curve, a noise ceiling, a
cell-line descriptor ablation or any VCC-comparable metric will not find one here.

## Mentioned entities

- Datasets: [[h1-hesc-crispri-perturb-seq]] (introduced here), [[replogle-nadig-crispri-collection]] (K562 and RPE1
  only, via ref. 57)
- Concepts: [[simple-perturbation-baselines]], [[context-count-learning-curve]],
  [[perturbation-response-decomposition]], [[systematic-variation]], [[directional-cell-quota]] (cells per
  perturbation, Fig. 6B), [[zero-shot-transfer-ceiling]] (what the source lacks)

## Connections

- Introduces / releases: [[h1-hesc-crispri-perturb-seq]] = the project's H1GW source (GSE295214; pp. 14, 20, 59–61).
- Cites and reuses as contexts: [[replogle-nadig-crispri-collection]] — K562 and RPE1 genome-scale CRISPRi (ref. 57;
  pp. 10, 12), day 6 vs day 8 K562 concordance (p. 58).
- Supports: [[simple-perturbation-baselines]] — its own zero-shot model loses to the train-mean baseline (p. 11);
  cites @ahlmanneltze-2025-linear-baselines as the reason the problem is open (p. 2).
- Supports (inferred): [[perturbation-response-decomposition]] — a context-free gene factor times a
  cell-state-conditioned map (p. 10), without variance shares; [[context-count-learning-curve]] — supplies only the
  within-context K-shot analogue and names context coverage as the bottleneck (pp. 10, 19).
- Does **not** support: [[zero-shot-transfer-ceiling]] — no ceiling, no split-half, no reliability estimate anywhere
  in the paper.
- Cites (ingested elsewhere in this wiki): @ahlmanneltze-2025-linear-baselines (ref. 9), @replogle-2022-genome-scale
  (ref. 57), @jiang-2025-pathway-signatures (ref. 62).
