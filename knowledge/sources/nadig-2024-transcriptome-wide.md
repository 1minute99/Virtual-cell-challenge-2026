---
title: "Nadig et al. 2024 (preprint) — Transcriptome-wide characterization of genetic perturbations (TRADE)"
type: source
created: 2026-11-06
updated: 2026-11-06
status: review
author: llm
bibkey: nadig-2024-transcriptome-wide
based_on: preprint
depth: standard
tags: [perturbation-effects, differential-expression, crispri, perturb-seq, noise-correction, cross-cell-line, effect-size-distribution, preprint]
relations:
  - target: replogle-nadig-crispri-collection
    type: supports
    confidence: extracted
    because: "Primary source of the HepG2 and Jurkat screens: «we performed two additional large-scale Perturb-seq experiments targeting common essential genes in Jurkat and HepG2 cell lines» (p. 7, ll. 151–153); 262,956 / 145,473 cells retained (p. 25, ll. 660–661)."
  - target: zero-shot-transfer-ceiling
    type: supports
    confidence: inferred
    because: "Shows on the same screens that sampling noise attenuates raw correlations (replicate log2FC r 0.16 vs noise-corrected 0.90, p. 13), the premise of a noise-corrected ceiling; the source states no ceiling and no split-half estimator."
  - target: perturbation-response-decomposition
    type: supports
    confidence: inferred
    because: "Noise-corrected cross-cell-line correlations of the same perturbation are well below 1 (median 0.40–0.75 by line pair, p. 13–14), consistent with a non-trivial line × perturbation component; the source does not decompose responses."
---

# Transcriptome-wide characterization of genetic perturbations

> Preprint version (nadig-2024-transcriptome-wide, bioRxiv v1, 2024-07-03) of nadig-2025-transcriptome-wide; claims quoted from the preprint.

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

## Bibliographic Details

Nadig, Ajay, Joseph M. Replogle, Angela N. Pogson, Steven A. McCarroll, Jonathan S. Weissman, Elise B. Robinson, and
Luke J. O'Connor. 2024. "Transcriptome-wide characterization of genetic perturbations." bioRxiv preprint, v1 posted
3 July 2024 (only version). doi:10.1101/2024.07.03.601903. CC BY-NC-ND 4.0. Code: R package at
github.com/ajaynadig/TRADE (p. 34). Data: SRA BioProject PRJNA1100571; processed single-cell data GEO GSE264667
(p. 35). Read from `<library>/pdf/nadig-2024-transcriptome-wide.pdf` (40 pp.). Citekey: @nadig-2024-transcriptome-wide;
published version @nadig-2025-transcriptome-wide (Nature Genetics 57(5):1228–1237, 2025).

**Page numbers.** The PDF has no printed page numbers; it carries continuous line numbers (1–1165) in the margin.
"p." below is the PDF page, and "l." the printed line number. The PDF contains main text, Methods and references only:
**Supplementary Figures 1–14, Supplementary Tables 1–20 and Supplementary Appendices 1–4 are not in the file** and were
not read. Figures sit on separate pages (Fig. 1 p. 6, Fig. 2 p. 10, Fig. 3 p. 12, Fig. 4 p. 15, Fig. 5 p. 18,
Fig. 6 p. 20; Table 1 p. 9). Greek symbols in the text layer are mis-encoded glyphs; quotes below restore them as β.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract | 2 |
| ✓ F1 | Introduction (SNP-heritability analogy; energy distance, GSEA, RRHO, iDEA as prior art) | 3–4 |
| ✓ F1 | Results: Overview of methods (β̂ = β + ε; ash mixture; transcriptome-wide impact; simulations at 200 cells) (Fig. 1, p. 6) | 5–6 |
| ✓ F1 | Results: Transcriptome-wide impact of 9,866 genetic perturbations (datasets; GATA1 vs EIF4A3; fraction in significant genes; downsampling; π_DEG) (Table 1 p. 9; Fig. 2 p. 10) | 7–10 |
| — | Results: Two types of gene-set enrichment (perturbation impact vs response enrichment; LOEUF, expression, growth, localization) (Fig. 3, p. 12) | 11–12 |
| ✓ F1 | Results: Consistency of transcriptome-wide effects across cell types (TI correlation; 241 cell-type-specific perturbations; replicate 0.16 vs 0.90; six line-pair medians; two clusters) (Fig. 4, p. 15) | 13–15 |
| — | Results: Dosage sensitivity of transcriptome-wide impact (Jost titration, SOX9 and Polycomb dTAG) (Fig. 5, p. 18) | 16–18 |
| — | Results: Greater transcriptomic than genetic correspondence across neuropsychiatric conditions (Fig. 6, p. 20) | 19–20 |
| ✓ F1 | Discussion (limitations: pseudobulk only, predefined labels; atlas design trade-off) | 21–22 |
| ✓ F1 | Methods: Jurkat and HepG2 Perturb-seq (cell lines, libraries dJR058/dJR092, growth screens, 10x, QC filters) | 23–25 |
| ✓ F1 | Methods: Modeling transcriptome-wide responses; features (TI, π_DEG, enrichments); univariate implementation | 25–28 |
| ✓ F1 | Methods: Simulations (N = 20 / 200 / 2000 cells per condition) | 28–29 |
| ✓ F1 | Methods: Genome-wide Perturb-seq (per-GEM-group pseudobulk, DESeq2, 272 batches → 4 mega-batches) | 29 |
| — | Methods: Gene annotations | 29–30 |
| ✓ F1 | Methods: Bivariate TRADE (mash, adaptive covariance grid); cell-type-specific perturbations; GMM clustering | 30–31 |
| — | Methods: attenuated guides, dTAG, dose-response simulations, PsychENCODE, OneK1K | 31–33 |
| — | Code / Data availability, acknowledgements, interests, references | 34–40 |
| ✗ | Supplementary Figures, Tables and Appendices (not in the PDF) | — |

## Kernthesen

### 1. ⭐⭐ TRADE estimates the distribution of true differential-expression effects from noisy estimates (pp. 2, 5, 26–27)

Each gene's estimated log2FC is modelled as β̂_g = β_g + ε_g (p. 5, l. 104). TRADE takes DESeq2 point estimates and
standard errors from pseudobulk counts and fits the distribution of β with ash, a mixture of half-uniform components
with likelihood β̂_g ~ N(β_g, ŝ_g²) (p. 27, ll. 731–747). It uses a uniform rather than null-biased prior «in order to
prevent bias in our distribution estimate» (p. 28, ll. 755–756). The aim is to interpret the fitted distribution itself,
not to shrink individual genes (p. 5, ll. 114–116). The approach «distinguishes properties of a cellular response from
those of an experiment, such as sample size and sequencing depth» (p. 5, ll. 108–109).

### 2. ⭐⭐ "Transcriptome-wide impact" (TI) is the variance of the true effect distribution (pp. 5, 26)

«The transcriptome-wide impact is defined as the variance of the distribution of differential expression effect sizes,
in units of log2-fold change» (p. 26, ll. 690–691). It is measured in log2FC² units: TI = 0.25 with normally
distributed effects on all genes means a typical effect of 0.5 (p. 26, ll. 693–695). A large TI can come from large
effects on a few genes or small effects on many (p. 5, ll. 121–123). The abstract calls it «an estimator of the overall
transcriptional effect of a perturbation which is stable across sampling depths» (p. 2, ll. 34–36). Fig. 2C reports TI
«scaled by the number of measured genes» (p. 10, ll. 220–221).

### 3. TI is robust in simulation, with a small downward bias that vanishes with more cells (pp. 5, 28–29; Fig. 1)

Simulations use negative-binomial counts with means and dispersions from the first 10 control batches of K562-GW,
3 sparsity levels × 10 TI values from 0.05 to 0.5, and 100 replicates (pp. 28–29). At 200 cells per condition the
estimates «were robust, with slight downward bias. The downward bias results from inadequate power to detect very small
effects at finite sample size, even aggregating power across genes; at increased sample size, the bias disappears»
(p. 5, ll. 130–134). The inferred distributions slightly overestimate the density at zero (p. 5, ll. 129–130).

### 4. ⭐ Raw log2FC distributions conflate effect and sampling noise; cell count drives the difference (p. 7; Fig. 2A–B)

GATA1 and EIF4A3 knockdowns in K562-Essential have near-identical observed log2FC distributions. TRADE gives
TI = 0.4 (mean |log2FC| 0.63) for GATA1 and TI = 0.004 (0.06) for EIF4A3 (p. 7, ll. 164–167). «the screen sequenced only
7 cells with EIF4A3 knockdown (as opposed to 108 cells with GATA1 knockdown), likely leading to large sampling variance
that inflated the observed effect size distribution» (p. 7, ll. 168–170). Fig. 2A–B annotate the observed variance as
0.59 and 0.60 (p. 10).

### 5. ⭐⭐ Significant genes carry only a minority of the true effect (p. 7; Fig. 2C)

«In the K562-GenomeWide experiment, only 36% of transcriptome-wide impact was explained by FDR-significant effects»
(p. 7, ll. 176–177). In the essential screens the share is 18 % (K562), 35 % (RPE1), 13 % (Jurkat) and 14 % (HepG2)
(p. 7, ll. 179–180). The Methods warn that these significant-gene shares «are expected to be upwardly biased by winner's
curse» (p. 27, ll. 716–718). Non-targeting guides give minimal TI in all cell types (p. 7, ll. 180–182; Supp. Fig. 2,
not in PDF).

### 6. ⭐ TI is stable under downsampling and across repeat screens; DEG counts are not (pp. 7–8)

With 50 % of the K562-GW GEM groups, signal in significant genes drops substantially while total TI «remained
relatively consistent» (p. 7, ll. 183–186; Fig. 2C). The small drop comes from TRADE being conservative when point
estimates are non-significant (p. 8, ll. 187–188). For perturbations shared between K562-GW and K562-Essential,
«transcriptome-wide impact R2 = 59.7%; number of DEGs R2 = 28.4%» (p. 8, ll. 190–192).

### 7. A typical perturbation affects ~45 genes; a typical essential-gene perturbation > 500 (pp. 2, 8; Fig. 2D)

The effective number of DE genes is π_DEG = 3M / κ, with M the number of measured genes and κ the kurtosis of the
inferred distribution (p. 26, ll. 699–707). The median π_DEG is 45 in K562-GW (p. 8, ll. 201–203). Essential-gene
perturbations in all four cell types have median π_DEG > 500 (p. 8, ll. 203–205). Fig. 2D is restricted to
perturbations «with nominally significant transcriptome-wide impact» (p. 10, ll. 222–223). The typical effect σ is
mostly in [0.1, 1] (p. 8, ll. 208–211).

### 8. Constrained and essential genes have high impact but low response; highly expressed genes have both (p. 11; Fig. 3)

In K562-GW, constrained (LOEUF) genes are enriched ~1.57× for perturbation impact and depleted 0.40× for perturbation
response. Growth-effect genes are 4.22× and 0.71×. Highly expressed genes are 2.26× and 4.44× (p. 11, ll. 244–252).
CRISPRi knockdown efficiency does not drive impact (p. 11, ll. 240–243). The authors read this as «population-level
constraint is mirrored by regulatory robustness» (p. 11, ll. 246–248).

### 9. ⭐⭐ Noise-corrected replicate correlation is 0.90 where the raw one is 0.16 (p. 13; Fig. 4B)

For the same perturbations in K562-GW and K562-Essential: «The median correlation between log-fold-change point
estimates - not using TRADE - was only 0.16, suggesting very low replicability (Figure 4B). However, the median
transcriptome-wide impact correlation between replicates was 0.90, implying excellent replicability» (p. 13,
ll. 285–288). The reason is that «uncorrelated sampling variation causes downward bias in correlation estimates»
(p. 13, ll. 289–290). Low replicate correlations occur mostly for perturbations with very low TI (p. 13, ll. 291–293).
The bivariate method is mash with a diagonal sampling covariance, an «adaptive grid» of 21 correlation values (pp. 30–31).
The "transcriptome-wide impact correlation" is therefore a noise-corrected correlation across genes of the true
log2FC vectors of one perturbation in two experiments (p. 13, ll. 266–269).

### 10. ⭐⭐ Perturbation magnitude is conserved across the four lines (average r = 0.62); RPE1 responds most (p. 13; Fig. 4A)

«As expected, transcriptome-wide impact was correlated across cell types (average correlation = 0.62; Supplementary
Figure 8)» (p. 13, ll. 270–271). «On average, transcriptome-wide impact was larger in the RPE1 cell line than in the
other three, indicating that this cell line is more sensitive to generic perturbations than the others» (p. 13,
ll. 271–273). This 0.62 correlates TI magnitudes across perturbations, not effect vectors.

### 11. ⭐ 241 perturbations have line-specific impact, most without a known explanation (p. 13; p. 31)

Regressing each line's log TI on the median log TI across the four lines (2,050 perturbations, intercept and slope fitted)
and flagging standardized residuals > 1.64 (p ≈ 0.1) gives 241 perturbations: K562 47, RPE1 118, Jurkat 10, HepG2 66
(p. 13, ll. 274–275; p. 31, ll. 888–899). Examples are GATA1 in K562 and HMGCR in Jurkat. «most had no previously
documented explanation for their cell-type dependent effects» (p. 13, ll. 278–279). The authors call the threshold
«liberal» (p. 13, l. 274).

### 12. ⭐⭐ Noise-corrected cross-line correlation of effect vectors is 0.40–0.75 and follows biology (pp. 13–14; Fig. 4C)

On 2,053 shared essential genes, restricted to 1,660 perturbations with TI Z > 0.5 in all four lines (p. 13,
ll. 294–299), the median TRADE correlations are K562/Jurkat 0.74, HepG2/RPE1 0.75, K562/HepG2 0.64, Jurkat/HepG2
0.69, K562/RPE1 0.40 and RPE1/Jurkat 0.56 (pp. 13–14, ll. 300–307). The high pairs share biology: K562 and Jurkat are
hematopoietic, p53-mutant and in suspension; HepG2 and RPE1 are epithelial, p53 wild-type and adherent (p. 13,
ll. 302–304). The correlations «did not vary strongly with essentiality» (p. 14, ll. 308–309).

### 13. ⭐ Two classes of perturbations: line-general (56 %) and line-group-specific (44 %) (p. 14; p. 31; Fig. 4D–E)

A bivariate Gaussian mixture on within- vs outside-similar-pair correlations gives «56% of the perturbations had high
correlations across all cell types (mean correlation within similar cell type pairs: 0.75; outside similar pairs:
0.66); 44% had higher correlations across similar cell types (mean correlation within similar cell type pairs: 0.61;
outside similar pairs: 0.35)» (p. 14, ll. 314–318). DDX41 is the example of the first class and NIFK of the second
(Fig. 4E: NIFK K562/RPE1 −0.14, RPE1/HepG2 0.85; p. 15).

### 14. Larger atlases should reveal more cell-type specificity (p. 13)

«As this dataset targeted primarily common essential genes which are expected to be important for growth across most
cell types, there are expected to be many more examples of cell-type-specific effects in a larger cellular perturbation
atlas» (p. 13, ll. 279–282).

### 15. Dose changes which genes respond, not only how much (pp. 16–17; Fig. 5)

Stronger perturbations have larger TI. For SOX9 and Polycomb dTAG, the dose-response is nonlinear, and the weakest and
strongest doses correlate only r = 0.60 and 0.48 (p. 16, ll. 382–397). K562 CRISPRi titration (25 genes, Jost 2020)
shows constant (BCR), threshold (GATA1) and gradient (ATP5E) patterns (p. 16, ll. 398–405). A Hill-equation toy model
reproduces the patterns (p. 17; p. 32).

### 16. Noise-corrected transcriptomic correlations between psychiatric disorders exceed genetic correlations (pp. 19–20)

Re-analysing PsychENCODE case/control statistics, TRADE correlations exceed Spearman correlations for 9/9 trait pairs
and exceed the genetic correlations, with low correlation to IBD as a control (p. 19, ll. 444–452). Microarray vs
RNA-seq cohorts correlate 0.96 / 0.91 / 0.78 for autism, bipolar disorder and schizophrenia (p. 19, ll. 462–464).

### 17. TI could serve as an unbiased cell-type distance (p. 21)

Unlike Euclidean distance, TI «is unbiased at finite sample size» and gives a more coherent PBMC hierarchy in OneK1K
(p. 21, ll. 500–507). TRADE needs predefined labels and cannot cluster cells (p. 21, ll. 507–508).

### 18. ⭐ Limitations and the atlas design argument: pseudobulk only; shallower sampling per perturbation (p. 22)

«A key limitation of TRADE is that it currently uses only a simple readout from single cell RNA-seq experiments, the
pseudo-bulk mean RNA expression level» (p. 22, ll. 520–521). On atlas design across contexts: «TRADE shifts the balance
in this tradeoff by allowing stable quantification of highly informative metrics including transcriptome-wide impact,
correlations between perturbations, and context-dependent effects at much shallower sampling depths» (p. 22,
ll. 534–536). The authors speculate that the GRN is densely connected with low modularity (pp. 21–22, ll. 516–519).

### 19. ⭐⭐ Two new essential-gene CRISPRi Perturb-seq screens in Jurkat and HepG2 (pp. 7, 9, 23–25)

«we performed two additional large-scale Perturb-seq experiments targeting common essential genes in Jurkat and HepG2
cell lines: Jurkat-Essential (2,393 essential gene perturbations in a T-cell leukemia cell line) and HepG2-Essential
(2,393 essential gene perturbations in a hepatocellular carcinoma cell line)» (p. 7, ll. 151–155). Design details are
in focus item 2 below. The effectors differ: Jurkat Perturb-seq uses Zim3-dCas9, HepG2 uses dCas9-BFP-KRAB (KOX1)
(p. 23, ll. 546–556; p. 24, ll. 606, 622).

## Figures and tables

- Table 1 (p. 9): **the seven datasets** with cell type, gene set, number of perturbations and **median cells per
  perturbation** (see focus item 2).
- Fig. 1 (p. 6): TRADE schematic (A), simulated distributions (B), estimated vs true TI (C).
- Fig. 2 (p. 10): GATA1 vs EIF4A3 (A–B), **TI in significant vs all genes, incl. 50 % downsample (C)**, π_DEG by
  dataset (D).
- Fig. 3 (p. 12): perturbation impact vs response enrichment by gene set.
- Fig. 4 (p. 15): line TI vs median TI with outliers (A), **replicate raw vs TRADE correlation (B)**, **median
  pairwise cross-line correlation (C)**, within vs outside similar pairs (D), DDX41 / NIFK joint distributions (E).
- Fig. 5 (p. 18): dose–TI curves (A), cross-dose correlation (B), toy model (C).
- Fig. 6 (p. 20): psychiatric cross-diagnosis correlations, TRADE vs Spearman vs genetic.

**Inconsistencies in the source.** (a) Supplementary Figure 8 is cited for both the knockdown-efficiency control
(p. 11, ll. 242–243) and the cross-cell-type TI correlation (p. 13, l. 271). (b) The psychiatric results cite
«Figure 5A and 5B» (p. 19, ll. 446, 449–450), but these panels are Figure 6 (p. 20). (c) The Jurkat/HepG2 QC sentence
gives «low UMI content (Jurkat: < 14%, HepG2: <18%) and high mitochondrial RNA content (Jurkat: > 1750 UMIs, HepG2:
> 3000 UMIs)» (p. 25, ll. 658–659), so the units look swapped. (d) The shared essential set is 2,053 genes (p. 13,
l. 295), the regression uses 2,050 (p. 31, l. 892), and Table 1 lists 2,057 for K562-Essential and 2,393 elsewhere.

## What this source does not address

- **No prediction.** No model predicts responses in an unseen line, and no transfer baseline or zero-shot metric
  appears. The cross-line correlations describe measured pairs only.
- **No variance decomposition into shared / line / interaction parts.** TI is the per-perturbation variance across
  genes, and the cross-line results are pairwise correlations. There is no μ + a + b + i partition (see
  [[perturbation-response-decomposition]]).
- **No signal-vs-noise share of the response variance.** The paper reports the share of TI in significant genes
  (13–36 %) and raw vs corrected correlations, not "x % of observed variance is noise". It gives observed variance
  0.59–0.60 against TI 0.4 / 0.004 only for the two examples in Fig. 2A–B.
- **No control-cell analysis** (no control PCs, no baseline state, no technical axes), and no split-half estimator: the
  correction uses DESeq2 standard errors in ash / mash.
- **No required cells-per-perturbation number.** Sample size enters only through simulations (20 / 200 / 2000 cells)
  and 50 % downsampling.
- **No context-count argument.** Four lines (two-by-two biological grouping); no drugs.

<!-- ═══ FOKUS ═══ -->

## Focus: facts and claims the manuscript cites; exactly what overlaps with H1–H4 — 2026-11-06

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

**Carries Kernthesen 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 18, 19.** Hypotheses are those of
`input/ideas/zero-shot-limits-context-complexity-plan.md`.

### What this focus needs beyond the Kernthesen

1. **TRADE in one sentence (for Methods / Related work).** «we introduce TRanscriptome-wide Analysis of Differential
   Expression (TRADE), a statistical framework which estimates the transcriptome-wide distribution of true differential
   expression effects from noisy gene-level measurements» (p. 2, ll. 31–33). Noise is corrected through per-gene DESeq2
   standard errors (ash for one experiment, mash for two), **not** through split halves (pp. 26–31). Our split-half
   estimator is therefore a different correction for the same attenuation problem.
2. **HepG2 / Jurkat screen design and numbers (Table 1, p. 9; pp. 23–25).**

   | Dataset | Source | Perturbations | Median cells / perturbation |
   |---|---|---|---|
   | K562-GenomeWide | Replogle 2022 | 9,866 genes | 178 |
   | K562-Essential | Replogle 2022 | 2,057 | 121 |
   | RPE1-Essential | Replogle 2022 | 2,393 | 72 |
   | Jurkat-Essential | new | 2,393 | 83 |
   | HepG2-Essential | new | 2,393 | 45 |

   Cells retained after QC: **Jurkat 262,956** (7,408 removed) and **HepG2 145,473** (15,952 removed) (p. 25,
   ll. 660–661). The design: dual-sgRNA CRISPRi library dJR092 (2,688 elements), targeting DepMap 20Q1 common
   essentials plus 5 % non-targeting controls, with guide abundance re-weighted from growth screens to reduce dropout
   («a guide with roughly four-fold depletion in growth screens was four-fold overrepresented», p. 23, ll. 572–580).
   It adds some genes with interesting K562-GW phenotypes (p. 23, ll. 580–582). Target infection is ~10 % (7 %
   measured), coverage > 1,000 cells per element, harvest at day 7, 10x 3′ v3 across **56 GEM groups** per line, aiming
   for ~15,000 cells per GEM group (pp. 24–25). Guide calling keeps cells with two same-gene sgRNAs or one sgRNA
   (p. 25, ll. 655–657). Effector: Zim3-dCas9 (Jurkat) vs dCas9-BFP-KRAB KOX1 (HepG2) (pp. 23–24). Compared with K562
   and RPE1, HepG2 has the fewest median cells (45 vs 72–178), and the HepG2 / Jurkat screens target the same 2,393
   genes as RPE1-Essential.
3. **Noise vs signal (H1 / H3a context).** The paper never states a noise share of perturbation-effect variance. What
   it shows: significant genes hold 13–36 % of TI (p. 7). The raw replicate log2FC correlation is 0.16 vs 0.90 after
   correction (p. 13), so with independent noise the raw correlation is attenuated about 5.6-fold. The EIF4A3 example
   (7 cells, observed variance 0.60, TI 0.004) shows that at very low cell counts almost all observed log2FC variance
   can be noise (p. 7; Fig. 2B, p. 10). *A reader's derivation, not in the source:* if attenuation were the only
   factor, r_raw / r_true ≈ 0.16 / 0.90 = 0.18 would imply a typical per-perturbation reliability product of about
   0.03. That is far below the project's split-half measurement ceilings (e.g. RN_K562 0.636). The two numbers are not
   comparable: TRADE correlates across all measured genes of one perturbation, including the many near-null genes,
   while ours is computed on 6,009 genes after filtering, per panel. Cite only the 0.16 / 0.90 contrast.
4. **Cross-line conservation (H2a / H2c / H3a context).** Four lines; TI magnitude correlation 0.62 (p. 13). Effect
   vectors: median noise-corrected r 0.40 (K562/RPE1) to 0.75 (HepG2/RPE1) (pp. 13–14). These are restricted to 1,660
   perturbations with TI Z > 0.5 in all four lines (p. 13), so they describe responsive perturbations only. *Reader's
   derivation:* on noise-free effects, a linear copy-the-other-line predictor leaves 1 − r² = 0.44 (HepG2/RPE1) to
   0.84 (K562/RPE1) of a perturbation's true effect vector unexplained (medians). This is a pairwise analogue of the
   interaction share, not the four-way partition of [[perturbation-response-decomposition]].
5. **Cells needed per perturbation.** No threshold is stated. The facts: simulations at N = 20 / 200 / 2000 cells per
   condition, and «N = 200 is similar to the typical sample size regime for the Replogle et al (2022) dataset»
   (p. 29, ll. 801–803). The TI bias «disappears» at larger N (p. 5, l. 133). Halving the K562-GW GEM groups leaves TI
   «relatively consistent» (p. 7, l. 186). TRADE allows key metrics «at much shallower sampling depths» (p. 22, l. 536).
   Fig. 2D and the cross-line analysis exclude non-significant perturbations, so low-cell perturbations drop out.

### Overlap with H1–H4

| Hyp. | What the source shows | What we add | Overlap level |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | Sampling noise strongly attenuates raw correlations (0.16 → 0.90 between K562 repeat screens) and hides most of the true effect (13–36 % in significant genes). A noise-aware estimator (ash / mash on DESeq2 SEs) is needed. | A ceiling r_max from a four-way decomposition, a split-half estimator, an attained fraction and an error split, all on transfer predictions. | **low (motivation only)** — supports the need for noise correction; no ceiling, no prediction (pp. 7, 13) |
| **H2a** (Tahoe learned context correction negative at C = 12; C*_emp ≥ 24) | Nothing on drugs or context-count curves. | Everything. | **none** |
| **H2b** (theory predicts break-even within ×2) | No theory of context sample complexity. | Everything. | **none** |
| **H2c** (genetic C*_min > C_public) | Supplies two of the public contexts (HepG2, Jurkat). Effects are only moderately conserved across lines (median r 0.40–0.75), 44 % of perturbations are line-group-specific, and «many more» cell-type-specific effects are expected in larger atlases (pp. 13–14). | C*_min, C_public counting, d_gene, σ²_ε / Var(i). | **low (adjacent; qualitative support for a non-trivial interaction)** |
| **H3a** (transfer residual on destination control PCs) | Line-pair similarity tracks p53 status, lineage and growth mode (p. 13). No control-cell or PC analysis. | Everything. | **none** |
| **H3b** (M3 vs five corrections) | No correction compared. | Everything. | **none** |
| **H4** (same-line genome-wide data beats modelling) | Argues for shallower sampling to scale atlases across contexts (p. 22). No data-vs-model comparison. | Everything. | **none** |

### Direct quotes (for the draft)

> "we introduce TRanscriptome-wide Analysis of Differential Expression (TRADE), a statistical framework which estimates
> the transcriptome-wide distribution of true differential expression effects from noisy gene-level measurements"
> (p. 2, ll. 31–33)

> "The transcriptome-wide impact is defined as the variance of the distribution of differential expression effect
> sizes, in units of log2-fold change." (p. 26, ll. 690–691)

> "In the K562-GenomeWide experiment, only 36% of transcriptome-wide impact was explained by FDR-significant effects"
> (p. 7, ll. 176–177)

> "The median correlation between log-fold-change point estimates - not using TRADE - was only 0.16, suggesting very
> low replicability (Figure 4B). However, the median transcriptome-wide impact correlation between replicates was 0.90,
> implying excellent replicability" (p. 13, ll. 285–288)

> "The highest median correlations were for K562/Jurkat (median correlation: 0.74) and HepG2/RPE1 (median correlation:
> 0.75)." (p. 13, ll. 300–302)

> "still weaker correlations for K562/RPE1 (0.40) and RPE1/Jurkat (0.56), suggesting that RPE1 cells tend to have
> especially unique responses to perturbations" (p. 14, ll. 306–307)

> "These filters removed 7408 cells from the Jurkat-Essential experiment (262956 cells retained), and 15952 cells from
> the HepG2-Essential experiment (145473 cells retained)." (p. 25, ll. 660–661)

> "N = 200 is similar to the typical sample size regime for the Replogle et al (2022) dataset." (p. 29, ll. 802–803)

### Examples & illustrations (for later drafting)

- GATA1 vs EIF4A3 in K562-Essential: the same observed spread, TI 0.4 vs 0.004, at 108 vs 7 cells (p. 7). This is the
  clearest single illustration of why noise correction must precede any ceiling or comparison.
- K562/RPE1 at 0.40 vs K562/Jurkat at 0.74 (pp. 13–14). Which source line you transfer from matters, and lineage,
  p53 status and growth mode predict it.
- NIFK: K562/RPE1 −0.14 vs RPE1/HepG2 0.85 (Fig. 4E, p. 15). One perturbation whose response flips between line
  groups.

### Boundary of this focus

This focus leaves out gene-set enrichment (Kernthese 8), dosage (15), psychiatric case/control (16) and the cell-type
distance (17). A reader looking for a ceiling, a transfer model, a variance partition, a control-PC analysis or a
cells-per-perturbation requirement will not find one here.

## Preprint vs published version

- **Version history (bioRxiv API, `api.biorxiv.org/details/biorxiv/10.1101/2024.07.03.601903`, queried 2026-11-06):**
  a single version, v1 dated 2024-07-03, type "new results", licence cc_by_nc_nd, linked to published DOI
  10.1038/s41588-025-02169-3.
- **Published version (Crossref + PubMed PMID 40259084, metadata and abstract only; not a quote source):** "Transcriptome-wide
  analysis of differential expression in perturbation atlases", Nature Genetics 57(5):1228–1237, epub 21 Apr 2025.
  Differences visible from the metadata and abstract: (1) **title changed**; (2) **an added author, Murthy** (eight
  authors vs seven); (3) the abstract's cross-cell-type statement changes from identifying «perturbations with
  cell-type dependent effects» (preprint, p. 2, ll. 42–43) to reporting moderate consistency of effects across cell
  types; (4) the preprint's phrase «stable across sampling depths» (p. 2, l. 35) and its claim that transcriptomic
  correlations «are greater than genetic correlations» (p. 2, ll. 46–47) are not in the published abstract, which
  instead says it "clarifies the relationship". The 45 / > 500 gene figures appear in both.
- **Body text, figures, numbers and supplements of the published version: not audited.** Every number on this page is
  from the preprint and must be cited as @nadig-2024-transcriptome-wide unless checked against the published PDF.

## Mentioned entities

- Datasets: [[replogle-nadig-crispri-collection]] (K562-GW, K562-Essential, RPE1-Essential re-analysed; HepG2- and
  Jurkat-Essential introduced here). Also Jost 2020 K562 titration, Naqvi 2023 SOX9 dTAG, Weber 2021 Polycomb dTAG,
  PsychENCODE (Gandal 2018), OneK1K (not wiki pages).
- Concepts: [[perturbation-response-decomposition]], [[zero-shot-transfer-ceiling]].
- Methods (no wiki page): TRADE, ash (Stephens 2016), mash (Urbut), DESeq2.

## Connections

- Supports: [[replogle-nadig-crispri-collection]] (extracted; primary source of HepG2 / Jurkat, pp. 7, 9, 23–25).
- Supports (inferred): [[zero-shot-transfer-ceiling]] (noise attenuation of raw correlations, p. 13) and
  [[perturbation-response-decomposition]] (cross-line correlations < 1, pp. 13–14).
- Data reused by: [[molina-2026-response-decomposition]] (four-line tensor), [[huuloc-2026-cell-quota]] ("TRADE
  essential-gene screens"), [[piao-2026-adaptive-contexts]] (all four lines, without citing Nadig).
- Cited by: [[dibaeinia-2026-context-scale]] (p. 1, for context dependence of perturbation phenotypes).
- Cites (not yet ingested): @replogle-2022-genome-scale and Jost et al. 2020 (titration).
