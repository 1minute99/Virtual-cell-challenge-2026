---
title: "Dibaeinia et al. 2026 — Virtual Cells Need Context, Not Just Scale"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: dibaeinia-2026-context-scale
depth: standard
tags: [position-paper, context-diversity, causal-transportability, cross-context, cd4-t-cells, deg-recovery, scaling]
relations:
  - target: zhu-cd4-t-cell-perturb-seq
    type: cites
    confidence: extracted
    because: "«we analyze a genome-scale Perturb-seq dataset of primary human CD4+ T cells (Zhu et al., 2025)» (p. 4); the only dataset analysed (App. C.1, p. 13)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "Cites Nadig et al. 2025 for «The same genetic perturbation produces dramatically different phenotypes depending on cell type» (p. 1); the data are not analysed."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "Qualitative parallel only: «observing the same perturbation across diverse contexts ... allows models to disentangle context-invariant effects from context-specific modulations» (p. 7). No components are estimated."
---

# Virtual Cells Need Context, Not Just Scale

## Bibliographic Details

Dibaeinia, Payam, Sudarshan Babu, Mei Knudson, Ali ElSheikh, Yibo Wen, Han Liu, Jason Perera, and Aly A. Khan. 2026.
"Virtual Cells Need Context, Not Just Scale." bioRxiv preprint, this version posted 9 February 2026 (Biohub; University
of Chicago; Northwestern University). doi:10.64898/2026.02.04.703804. CC BY 4.0. Code: github.com/czbiohub-chi/scldm_cd4;
checkpoints: huggingface.co/biohub/scldm_cd4 (p. 9).
Read from `<library>/pdf/dibaeinia-2026-context-scale.pdf` (16 pp.; printed page numbers = PDF page numbers; main text
pp. 1–8, references pp. 9–11, appendices A–D pp. 12–16). The PDF itself does not mention ICML; the ICML 2026
position-paper status comes from the bib note only and is not verified here.
Citekey: @dibaeinia-2026-context-scale.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; 1 Introduction (position statement, three contributions) | 1 |
| ✓ F1 | 2 Related Work (scaling paradigm; empirical challenges; context-dependence; causal transportability) | 2 |
| ✓ F1 | 3.1–3.3 The standard formulation f(x, p); the context problem f_c; why implicit context inference fails | 2–3 |
| ✓ F1 | 3.4–3.5 Connection to causal transportability (Eq. 3, Def. 3.1, Prop. 3.2); contrast with distribution shift | 3–4 |
| ✓ F1 | 4.1 Simple baselines match deep learning (review of published benchmarks) | 4 |
| ✓ F1 | 4.2 A controlled test of context versus scale (CD4+ T-cell task, Fig. 2, metrics) | 4–5 |
| ✓ F1 | 4.3 Central finding: context diversity drives generalization (Figs. 3–4) | 5–6 |
| ✓ F1 | 4.4 Interpretation ("1 million cells across 100 contexts") | 7 |
| — | 5 Alternative views (four counterarguments and responses) | 7–8 |
| ✓ F1 | 6 Call to action (data generation, benchmarking, modelling) | 8 |
| ✓ F1 | 7 Discussion and conclusion (incl. Limitations) | 8 |
| — | App. A Model descriptions (scVI … scLDM) | 12 |
| — | App. B Extended context examples | 12–13 |
| ✓ F1 | App. C Full methodology (data, QC, perturbation filter, split, scLDM, baselines, evaluation; Fig. C1) | 13–15 |
| — | App. D Discussion of evaluation metrics | 15–16 |

## Kernthesen

### 1. ⭐⭐ The main failure of virtual-cell models is missing context coverage, not model capacity or data volume (p. 1)

The position, in the authors' words: «current predictive models of cellular responses fail primarily due to insufficient
coverage of biological contexts, rather than limitations in data volume or model capacity» (p. 1). The reason given is
that cellular response, unlike protein folding, is context-dependent: «The same genetic perturbation produces
dramatically different phenotypes depending on cell type, activation state, microenvironment, and epigenetic landscape»
(p. 1, citing Nadig et al. 2025 and Frangieh et al. 2021). If the mapping varies with context, «accumulating
observations from the same contexts yields diminishing returns» (p. 1). The paper calls for course correction, not
abandonment of scaling (p. 8).

### 2. ⭐ Context is an effect modifier, so the prediction target is a family of functions f_c (pp. 2–3)

The standard formulation f : (x, p) ↦ y (Eq. 1, p. 3) hides the context c. «When c varies, the function f itself
changes» to f_c (Eq. 2, p. 3). Context «determines not just initial baseline expression but how perturbations affect
the cell» (p. 3). Context covers microenvironment, tissue, developmental and epigenetic history, activation state,
temporal dynamics and genetics (p. 3; App. B, pp. 12–13). Recovering f_c from data where c is implicit «is a problem
that is not identifiable in general without additional assumptions» (p. 3).

### 3. ⭐ Implicit context inference from the transcriptome fails for three reasons (p. 3)

(1) *Lossy projection*: chromatin, protein localisation, PTMs and metabolites are not seen in RNA, so two cells with
similar transcriptomes can respond differently. (2) *Identifiability*: «most training datasets cover only a few
contexts with limited overlap of perturbations across contexts. The model cannot uniquely separate context effects
from perturbation effects» (p. 3). (3) *Causal non-invariance*: predictive signals for c in x «need not correspond to
the causal factors that determine f_c» (p. 3). The same three points are restated against the counterargument
"context could be inferred from transcriptomes" (p. 7).

### 4. ⭐⭐ Perturbation prediction across contexts is a causal transportability problem (pp. 3–4)

The context-specific effect is defined as a CATE: δ(p; c) = E[Y | do(P = p), C = c] − E[Y | do(P = ctrl), C = c]
(Eq. 3, p. 3). Contexts are domains π, π*; selection diagrams mark which mechanisms differ (Def. 3.1, p. 3, after
Pearl & Bareinboim 2011; Bareinboim & Pearl 2016). Proposition 3.2 (informal) states: «When changes in context
modulate the mechanism linking perturbation P to response Y, predictions learned in source domain π need not transfer
to target domain π* with different contexts, without additional invariance assumptions or target-domain
interventions» (p. 4). The proposition is informal and carries no proof, no rate and no count of contexts.

### 5. More data from the same contexts refines f_c on the training contexts only (p. 4)

From Prop. 3.2: more data «refines performance for the training distribution of contexts, but does not inform the
context-specific mechanism f_c* in a new target context c* ∉ supp(P_π(C))» (p. 4). This is the theoretical basis of the
paper's "diminishing returns" claim.

### 6. ⭐ This is mechanism shift, not covariate shift; reweighting cannot fix it (p. 4)

Covariate-shift and domain-adaptation methods assume a stable conditional P(Y | X) (p. 4). Here P(Y | X, P) itself varies
(«a mechanism shift»), so «no amount of reweighting or representation learning on source data alone can recover the
target mechanism» (p. 4). The two remedies named: «(1) sufficient contextual diversity in training data to observe
perturbation effects across the range of mechanisms encountered at test time, or (2) explicit modeling of how context
modulates response mechanisms» (p. 4).

### 7. Published benchmarks show simple baselines match deep models and many models collapse to averages (p. 4)

Reviewed evidence (no new analysis): on Norman combinations an additive baseline beats all deep models; for single
perturbations no deep model consistently beats a mean baseline (Ahlmann-Eltze et al. 2025; Li et al. 2025b). scGPT, UCE
and scBERT predictions vary little across perturbations ("mode collapse", "variance compression"). Systematic variation
correlates with performance at r = 0.91 (fine-tuned scGPT) and r = 0.95 (GEARS) (Viñas Torné et al. 2025, cited p. 4).
The paper's reading: «narrow context coverage may bias models toward learning context-specific correlations rather than
transportable mechanisms» (p. 4).

### 8. ⭐⭐ The controlled test uses one cell type with up to 12 donor × timepoint contexts (pp. 4–5; App. C, pp. 13–14)

Data: genome-scale CRISPRi Perturb-seq of primary human CD4+ T cells ([[zhu-cd4-t-cell-perturb-seq]]): 22 M cells,
≈12,000 knockdowns, 4 donors, 3 timepoints (rest, 8 h, 48 h) (p. 4). Context = donor × timepoint, «yielding up to 12
distinct biological contexts» (p. 5), cell type held fixed. After filtering (knockdown ratio < 0.5 in ≥ 2 donors per
timepoint, ≥ 100 cells): ≈18 M cells, 10,570 perturbations (p. 13). Split (Fig. 2A, p. 5; C.4, p. 14): donors 1–2 at all
three timepoints and donors 3–4 at rest are training contexts (8); the held-out contexts are donors 3–4 at 8 h and 48 h
(4). Within each held-out context 30 % of eligible perturbations go to validation, 40 % to test, and «the remaining
perturbations were included in training, providing some in-context supervision» (p. 14). Test perturbations must be
observed in training contexts (p. 14). Final split: 100,125 / 7,449 / 16,527 configurations (p. 14).

### 9. scLDM beats perturbation-mean baselines on Δ metrics and ranks perturbations well (pp. 5–6; Fig. 3; Fig. C1)

scLDM (latent diffusion, joint donor/timepoint/perturbation conditioning; C.5, p. 14) outperforms the mean baselines P,
P+D and P+T (C.6, p. 14) on Pearson-Δ, cosine-Δ and L2-Δ (Fig. 3A–C, p. 6). L2-Δ discrimination scores have median > 0.9
in held-out contexts (Fig. 3D; Fig. C1, p. 15). The authors conclude that «with explicit context conditioning, some
cross-context transfer is achievable in T cells» (p. 5).

### 10. ⭐ DEG recovery is poor (recall ≈ 9 %) and barely correlated with aggregate metrics (p. 6; Fig. 3E–F)

DEG precision has median ≈ 0.67, recall median ≈ 0.09: «the model recovers only ∼9% of true DEGs» (p. 6). DEG-F1
correlates weakly with Correlation-Δ (r = 0.26) and L2-Δ (r = −0.05) (Fig. 3E–F, p. 6). «Many configurations achieve
Correlation-∆ ≥ 0.7 yet DEG-F1 < 0.1» (p. 6). Baselines are excluded from DEG analysis because they output one
pseudobulk vector (p. 15).

### 11. ⭐⭐ DEG-F1 rises with the number of training contexts per perturbation, not with cell count (p. 6; Fig. 4A–D)

For each test triplet (d, t, p) the authors record (i) the number of distinct training donor-timepoint contexts in which
p was observed and (ii) the number of training cells with p (p. 6). «Mean DEG-F1 increases from < 0.1 to ∼ 0.19 when
moving from ≤ 3 to 8 contexts; Figure 4A» (p. 6); L2-Δ shows «a weaker trend» (Fig. 4B). DEG-F1 vs training cell count:
r = 0.11 (Fig. 4C); L2-Δ vs cell count: r = −0.38 (Fig. 4D), «suggesting that scale can improve standard scores even
when biologically meaningful DEG recovery does not improve» (p. 6). Exact panel contents are listed in focus item 1.

### 12. ⭐⭐ At matched cell count, 8 contexts beat < 8 contexts (p. 6; Fig. 4E–F; C.7, p. 15)

Restricting to configurations whose perturbation had 700–1700 training cells, the 8-context group beat the < 8 group:
p = 0.002 (DEG-F1) and p = 0.01 (L2-Δ), one-sided Mann–Whitney U (p. 6; test stated p. 15). «Even at controlled data
volume, context diversity improves generalization» (p. 6).

### 13. ⭐ Diverse contexts let a model separate context-invariant effects from context-specific modulations (p. 7)

Interpretation (§4.4): more cells from the same context refine that context's mechanism «but provides limited
information about held-out contexts», whereas observing a perturbation across contexts «allows models to disentangle
context-invariant effects from context-specific modulations. This is precisely the knowledge required for
cross-context generalization» (p. 7).

### 14. ⭐⭐ "1 million cells across 100 contexts can be more valuable than 10 million cells across 10 contexts" (p. 7)

The practical implication drawn (p. 7). The main barrier is «datasets that sample contexts too narrowly to support causal
transport across the biological contexts where predictions must apply» (p. 7). The 100-vs-10 comparison is a rhetorical
example. It is not measured; the paper's own data reach at most 8 training contexts per perturbation.

### 15. Scaling worked for AlphaFold and LLMs because their targets are context-invariant or context is explicit (p. 7)

Against "scale will eventually solve this" (p. 7): protein folding is largely governed by universal physics (Anfinsen
1973), and language carries explicit context cues; transcriptomes «typically lack explicit context markers for many
relevant variables» (p. 7). The authors «remain open to surprises» (p. 7).

### 16. Discriminative benchmarks do not test the counterfactual task (p. 7)

Foundation models help cell-type annotation and batch correction, which can rely on context-invariant patterns.
Perturbation prediction «is generative and counterfactual» and depends on context-specific mechanisms «that standard
discriminative benchmarks do not probe» (p. 7).

### 17. ⭐ Data generation should prioritise breadth of contexts; budget forces a breadth-depth trade-off (p. 8)

«Contextual diversity matters as much as cell count for context generalization (if not more)» (p. 8). Screens should
prioritise coverage across donors, cell types, activation states and timepoints. «Under a fixed sequencing budget, there
is an inherent tradeoff between broad context coverage and cell depth per context» (p. 8). No optimum is derived.

### 18. ⭐ Benchmarks should report new perturbations, new contexts, and both separately, with DEG metrics (p. 8; App. D)

Few benchmarks stress-test cross-context transport (p. 8). The authors propose three difficulty levels: «new
perturbations, new contexts, and new perturbations in new contexts» (p. 8), plus multi-metric evaluation including DEG
recovery (p. 8; App. D, pp. 15–16). App. D notes that MAE/MSE reconstruction error is often «on the same order as the
typical per-gene perturbation effect size» and that the 2025 Virtual Cell Challenge found most models worse than a
mean baseline on MAE (p. 15).

### 19. Models should make context explicit or be context-specialised (p. 8)

Options named: conditioning on observed context covariates, in-context learning with cell sets as prompts (Dong et al.
2026, "Stack"), causal representation learning under sparse mechanism shift (Lopez et al. 2023), or context-specialised
models: «a Virtual CD4+ T cell» «may be more actionable than a universal Virtual Cell» (p. 8).

### 20. ⭐⭐ Stated limits: transcriptomics only, discrete contexts, and no answer to "how much diversity is enough" (p. 8)

Limitations (p. 8): the empirical analysis is transcriptomic only; «Our theoretical framework assumes context is
discrete and annotatable, whereas real biological contexts often exist on a continuum»; and «We do not address how much
contextual diversity is "enough," as this likely depends on the target application» (p. 8).

## Figures and tables

- Fig. 1 (p. 2): schematic, protein folding (context-invariant) vs cellular response (context-dependent).
- Fig. 2 (p. 5): **(A) donor × timepoint split** (8 training contexts, 4 held-out with train/valid/test perturbations);
  (B) UMAP true vs generated for 25 test configurations; (C) LCP2 knockdown shift.
- Fig. 3 (p. 6): (A–C) scLDM vs mean baselines P, P+D, P+T on Pearson-Δ, cosine-Δ, L2-Δ; (D) discrimination score for
  donor CE0008678 at 48 h; (E–F) DEG-F1 vs Pearson-Δ (r = 0.26) and vs L2-Δ (r = −0.05).
- **Fig. 4 (p. 6): DEG-F1 and L2-Δ vs number of training contexts per perturbation (A–B), vs training cell count
  (C–D), and at matched cell count, < 8 vs 8 contexts (E–F).** Details in focus item 1.
- Fig. C1 (p. 15): discrimination scores for the remaining held-out contexts; DEG-F1-based discrimination is lower.
- No tables in the paper.

## What this source does not address

- **No context-count curve with a baseline, and no break-even.** Fig. 4 stratifies one trained model's test triplets by
  per-perturbation coverage; it does not vary the number of training contexts, and the mean baselines are excluded
  from DEG analysis (p. 15). The project's break-even design (context-aware vs context-free, C varied in training) is
  in the plan (`input/ideas/zero-shot-limits-context-complexity-plan.md`, H2a/H2b) and has no counterpart here.
- **No quantitative theory.** Prop. 3.2 is informal (p. 4); there is no sample-complexity bound, no effective
  dimension d, no noise term. Domain-generalisation counts (E > d) are in Rosenfeld et al. 2021 and related work
  (literature guide), not here.
- **No cell lines, no unseen donor or unseen timepoint.** Contexts are donor × timepoint within CD4+ T cells; each held-out
  context combines a donor and a timepoint that are both seen in training (Fig. 2A, p. 5), and 30 % of held-out-context
  perturbations are trained on (p. 14). Zero-shot transfer to unseen cell lines is not tested. For cell-line
  zero-shot results see [[molina-2026-response-decomposition]].
- **No noise ceiling or split-half reliability.** DEG recall of ≈ 9 % is not compared with a measurement ceiling; see
  [[zero-shot-transfer-ceiling]] and @wang-2026-reliable-perturbations (not yet ingested).
- **No count of available public contexts.** «most training datasets cover only a few contexts» (p. 3) is not
  quantified.
- **No variance decomposition** of responses into conserved vs context-specific parts; that is
  [[perturbation-response-decomposition]].

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 3, 4, 8, 11, 12, 14, 17, 20** (and 10, 18 for metrics; 19 for H4). The hypotheses are those of the
plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3 at the end).

### What this focus needs beyond the Kernthesen

1. **Fig. 4 exactly (read from the rendered figure, p. 6; bar heights are approximate visual readings).**
   - (A) y: DEG-F1 (0–0.3); x: «No. Train Configs With Same Perturbation», bins 1–8. N per bin: 1, 9, 14, 34, 55, 68,
     87, **5,331**. Bar heights ≈ 0.14, 0.13, 0.04, 0.13, 0.13, 0.15, 0.17, 0.19. The curve is not monotone: bins 1–2
     (N = 1 and 9) sit near 0.13–0.14, bin 3 is the minimum, then a monotone rise from bin 3 to 8. The text's
     «< 0.1 … from ≤ 3» holds for the pooled ≤ 3 group (N-weighted ≈ 0.08 from these readings), not for each bin.
   - (B) y: L2-Δ (0–15+); x: bins 2–8; N: 1, 6, 38, 72, 116, 316, **15,599**; heights ≈ 15, 10.7, 10.4, 10.5, 9.3,
     12.4, 9.3. Not monotone (bin 7 is the second worst).
   - (C) DEG-F1 vs «No. Train Cells With Same Perturbation» (0–≈5,500): r = 0.11, p = 1.7e-16. (D) L2-Δ vs the same:
     r = −0.38, p reported as 0.0e+00.
   - (E) DEG-F1, < 8 (N = 157) vs 8 (N = 3,628) contexts at 700–1,700 training cells: p = 0.0023. (F) L2-Δ, < 8
     (N = 125) vs 8 (N = 3,660): p = 0.0130. Medians differ modestly (E: ≈ 0.14 vs ≈ 0.18, visual).
   - **Dataset:** CD4+ T-cell CRISPRi only; **model:** scLDM only; **no baseline** on any panel; **no crossing** of
     any kind. Over 95 % of test triplets sit in the 8-context bin, so the low-context bins are thin.
   - DEG metrics use only test triplets with ≥ 100 perturbed cells (p. 15), hence N = 5,599 in (A) vs ≈ 16,150 in (B).
   - Ambiguity: the axis tops out at 8, matching the 8 fully trained contexts of Fig. 2A. The 30 % training share inside
     held-out contexts (p. 14) would allow up to 11 unless perturbations were assigned to train/valid/test jointly
     across held-out contexts. The paper does not say.
2. **Confound not discussed by the authors (reader's inference).** The context count per perturbation is observed, not
   assigned. Perturbations are kept only where knockdown ratio < 0.5 in ≥ 2 donors per timepoint and ≥ 100 cells remain
   (p. 13), so perturbations seen in fewer contexts are plausibly weaker or less efficient knockdowns. Low DEG-F1 at
   low coverage may partly reflect perturbation strength. Fig. 4E–F controls cell count, not effect size.
3. **How many contexts are needed.** The paper gives no number. The closest statements: the qualitative 100 vs 10
   contexts example (p. 7), «sufficient contextual diversity» (p. 4), and the explicit non-answer (p. 8). Its largest
   measured coverage is 8 contexts per perturbation; the project's C*_emp = 37.4 [33.3, > 40] (Tahoe, cell lines) lies
   far outside the range the paper explores.
4. **"Existing data are insufficient" (H2c).** Stated qualitatively twice: «most training datasets cover only a few
   contexts with limited overlap of perturbations across contexts» (p. 3) and «datasets that sample contexts too
   narrowly» (p. 7). No inventory of public contexts. Note for C_public: the paper treats donor × activation state as
   12 contexts within one cell type (4 donors × 3 timepoints; App. C.1, p. 13). Under the plan's operational definition 5
   (distinct cell lines/types) this dataset counts as **one** context (primary CD4+ T cells); it is relevant to the
   A2 Count-2 ambiguity (donors), not to Count 1.
5. **Identifiability argument (relation to Proposition 1).** The paper's argument is verbal: «most training datasets
   cover only a few contexts with limited overlap of perturbations across contexts. The model cannot uniquely separate
   context effects from perturbation effects; multiple context-dependent mechanisms can fit the training data equally
   well, making the inference problem underdetermined» (p. 3). This is the same qualitative claim as the project's
   Proposition 1 (rank ≤ C − 1), without a rank statement, a descriptor dimension or a noise term; Proposition 3.2 is
   labelled «(Informal)» (p. 4).

### Overlap with H1–H4

| Hyp. | Overlap | What the source shows (pages) |
|---|---|---|
| H1 (ceiling) | none | No noise ceiling, no split-half reliability; DEG recall ≈ 9 % is not compared with any ceiling (p. 6). |
| H2a (empirical break-even) | partial (qualitative precedent only) | Fig. 4 (p. 6): DEG-F1 rises with per-perturbation context coverage (≤ 8 donor × timepoint contexts) for one model, observational, no context-free baseline, no crossing. Our C*_emp = 37.4 is a different design (training C varied, matched context-free ablation, unseen cell lines). |
| H2b (theory predicts C*) | none | Prop. 3.2 is informal and non-quantitative (p. 4); no formula to compare with the rejected linear-Gaussian C*. |
| H2c (genetic data below break-even) | partial (claim without number) | «most training datasets cover only a few contexts» (p. 3); «We do not address how much contextual diversity is "enough"» (p. 8). No count of public contexts. |
| H3a (residual on destination control PCs) | none | Not addressed. |
| H3b (M3 vs alternatives) | none | Not addressed. |
| H4 (data vs modelling gain) | partial (position) | Breadth over depth: «A dataset with 1 million cells across 100 contexts can be more valuable than 10 million cells across 10 contexts» (p. 7); Fig. 4C–F (p. 6): context coverage matters more than cell count. H4 adds a same-line screen (depth for the destination), which this paper does not test. |

### Direct quotes (for the draft)

> "current predictive models of cellular responses fail primarily due to insufficient coverage of biological contexts,
> rather than limitations in data volume or model capacity" (p. 1)

> "most training datasets cover only a few contexts with limited overlap of perturbations across contexts. The model
> cannot uniquely separate context effects from perturbation effects; multiple context-dependent mechanisms can fit the
> training data equally well, making the inference problem underdetermined." (p. 3)

> "A dataset with 1 million cells across 100 contexts can be more valuable than 10 million cells across 10 contexts."
> (p. 7)

> "We do not address how much contextual diversity is "enough," as this likely depends on the target application."
> (p. 8)

> "the Virtual Cell Challenge found that most models performed worse than a mean baseline on MAE" (p. 15, App. D)

### Examples & illustrations (for later drafting)

- Primary CD4+ T cells, 4 donors × 3 activation timepoints (resting, 8 h, 48 h), ≈ 12,000 CRISPRi knockdowns, 22 M
  cells (p. 4; App. C.1, p. 13); ≈ 18 M cells and 10,570 perturbations after filtering (App. C.3, p. 13).
- LCP2 knockdown shift in a held-out context (Fig. 2C, p. 5).
- Analogy: characterising drug responses in naïve T cells and expecting transfer to exhausted T cells (p. 7).

### Boundary of this focus

The source supplies a position and a within-cell-type, observational context-coverage trend. It supplies no
cell-line zero-shot test, no baseline-referenced context curve, no break-even number and no theory to calibrate
against; for those see [[molina-2026-response-decomposition]] and the project's own Task 7/8 results.

## Mentioned entities

- Datasets: [[zhu-cd4-t-cell-perturb-seq]], [[replogle-nadig-crispri-collection]] (cited only, p. 1)
- Concepts: [[perturbation-response-decomposition]] (qualitative parallel), [[zero-shot-transfer-ceiling]] (absent here)

## Connections

- Cites: [[zhu-cd4-t-cell-perturb-seq]] (only dataset analysed), [[replogle-nadig-crispri-collection]] (p. 1).
- Supports (qualitatively, inferred): [[perturbation-response-decomposition]] — context-invariant vs context-specific
  effects (p. 7).
- Supplements: [[molina-2026-response-decomposition]] (cell-line zero-shot evidence that this paper lacks); plain
  mention, no typed edge.
- Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (H2a, H2c, H4).
