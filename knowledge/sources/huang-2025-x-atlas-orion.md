---
title: "Huang et al. 2025 — X-Atlas/Orion: genome-wide Perturb-seq via Fix-Cryopreserve (FiCS) Perturb-seq"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: huang-2025-x-atlas-orion
depth: standard
tags: [dataset, perturb-seq, crispri, genome-wide, hct116, hek293t, knockdown-efficiency, sgrna-dosage, data-quality, batch-effects]
relations:
  - target: x-atlas-orion
    type: supports
    confidence: extracted
    because: "Introduces and releases the dataset: «we release X-Atlas: Orion edition (X-Atlas/Orion), the largest publicly available Perturb-seq atlas» (p. 1); Table 1 QC metrics (p. 7); Figshare DOI (p. 23)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "Uses Replogle K562 (ref. 17) as the benchmark: «Both genes and transcripts detected per cell compared favorably with the benchmarking dataset from Replogle et al. 2022» (p. 5); batch correlation 0.993/0.988 vs 0.967 (Fig. 2F, p. 7); discernible phenotypes 48.36/40.37 vs 54.30 % (p. 8)."
  - target: tahoe-100m
    type: cites
    confidence: extracted
    because: "Depth comparator only: HCT116 median UMIs «8.45 times as deep as the Tahoe-100M chemical perturbation dataset … (median UMIs per cell = 2298)» (p. 5; Fig. 2E)."
  - target: simple-perturbation-baselines
    type: supports
    confidence: extracted
    because: "Motivating citation (refs. 9–11, Bendidi / Li / PertEval): foundation models «often struggle … to predict the outcome of specific interventions, underperforming simpler models» (p. 2). Framing only; no models are evaluated."
  - target: huuloc-2026-cell-quota
    type: contradicts
    confidence: ambiguous
    because: "Headline detectability diverges: this source's own classifier finds 48.36 % (HCT116) / 40.37 % (HEK293T) of expressed-gene perturbations discernible (q < 0.05; Table 1 p. 7, p. 8), whereas huuloc report 14.7 % / 35.3 % clearing their angular detection floor. Criteria, denominators and cell filters differ, so whether this is a contradiction is unclear."
---

# X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models

## Bibliographic Details

Huang, Ann C., Tsung-Han S. Hsieh, Jiang Zhu, Jackson Michuda, Ashton Teng, Soohong Kim, Elizabeth M. Rumsey, Sharon K.
Lam, Ikenna Anigbogu, Philip Wright, Mohamed Ameen, Kwontae You, Christopher J. Graves, Hyunsung John Kim, Adam J.
Litterman, Rene V. Sit, Alex Blocker, and Ci Chu. 2025. "X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable
Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models." bioRxiv preprint, version posted
16 June 2025. doi:10.1101/2025.06.11.659105. CC BY-NC 4.0. Xaira Therapeutics and Foresite Labs; all authors hold
equity in Xaira (p. 23). Data: Figshare, doi:10.25452/figshare.plus.29190726 (p. 23).
Read from `<library>/pdf/huang-2025-x-atlas-orion.pdf` (27 PDF pages). **The preprint has no printed page numbers, so
every "p. N" below is the PDF page.** Citekey: @huang-2025-x-atlas-orion.

**Arithmetic notes on the source's own numbers.** "Eight million cells" (pp. 1, 6) = 3,409,169 + 4,534,299 = 7.94 M
two-sgRNA cells (Table 1, p. 7). "Approximately 30 million sequenced barcodes" (p. 6) = 9,886,825 + 21,550,681 =
31.4 M aligned cells. "Over 16,000 UMIs per cell" (p. 1) is the lower of the two medians (19,416 HCT116 / 16,557
HEK293T). The library is given as «41,780 sgRNAs targeting 18,903 genes and 1,026 non-targeting control pairs» (p. 5).
2 × 18,903 + 2 × 1,026 = 39,858, not 41,780, so the extra ~1,900 guides are not explained (possibly alternative pairs
for some genes; not stated).

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract | 1 |
| ✓ F1 | Introduction (foundation models underperform simple models; need for interventional data; dose as a continuous variable) | 2–3 |
| ✓ F1 | Results: Establishing the FiCS Perturb-seq Platform (four bottlenecks, five optimisations; Fig. 1, p. 4) | 3–4 |
| ✓ F1 | Results: FiCS Perturb-seq generates high-quality genome-wide data (Fig. 2, pp. 6–7; Table 1, p. 7) | 5–7 |
| ✓ F1 | Results: FiCS Perturb-seq detects perturbation effects and recovers known biology (Fig. 3, pp. 9–10) | 8–10 |
| ✓ F1 | Results: sgRNA abundance stratifies KD efficiency and detects dose-dependent genetic effects (Fig. 4, pp. 11–12) | 10–12 |
| ✓ F1 | Conclusions | 12–13 |
| ✓ F1 | Methods: Library design and cloning (dual-sgRNA, tRNA separators, no expression cutoff) | 13–14 |
| — | Methods: Lentivirus production | 14 |
| ✓ F1 | Methods: Cell culture and cell line generation (dCas9-mScarlet-ZIM3, PiggyBac) | 14–15 |
| ✓ F1 | Methods: Perturbation screening (MOI 0.1, coverage, puromycin, day 7) | 15 |
| ✓ F1 | Methods: Cell sorting, fixation and cryopreservation (Zombie NIR, DSP, Bambanker + RNasin) | 15 |
| — | Methods: Cell thawing and preparation for 10x scRNA-seq | 15–16 |
| ✓ F1 | Methods: Cell counting & 10x GEM-X Single Cell 5′ v3 (superloading 100k / lane) | 16 |
| — | Methods: Automation (Hamilton Vantage configuration, liquid classes; Table 2, pp. 17–18) | 16–18 |
| ✓ F1 | Methods: QC, Pooling & Sequencing (50,000 / 5,000 reads per cell) | 18 |
| ✓ F1 | Methods: Basecalling, alignment, and sgRNA calling (GMM thresholds; expressed-gene definition) | 18–19 |
| ✓ F1 | Methods: Batch-batch correlation | 19 |
| ✓ F1 | Methods: Binary classification (L1 logistic regression, Storey q) | 19 |
| — | Methods: Gene-gene correlation and STRING scores | 20 |
| ✓ F1 | Methods: Clustering of perturbations with strong phenotypes (NTC PCA, 20 PCs) | 20 |
| ✓ F1 | Methods: Stratification by sgRNA and energy distance analysis | 21 |
| ✓ F1 | Supplemental Figures 1–2 (fresh vs fixed, sorting gates, tRNA recombination, cross-line batch correlation, classifier controls) | 22–23 |
| ✓ F1 | Acknowledgements; Data Availability; Declaration of Interests | 23 |
| — | References | 24–26 |

## Kernthesen

### 1. ⭐ Genome-wide Perturb-seq is blocked by four logistics bottlenecks, which FiCS removes with five workflow changes (pp. 3–4; Fig. 1B)

The authors name the obstacles to scaling: (1) single-cell platforms want freshly harvested cells, and handling stress
confounds phenotypes, made worse by the mismatch between tens of millions of input cells and about 200,000 cells per
GEM-X batch; (2) FACS enrichment is too slow and stressful at genome scale; (3) fresh cells allow no upfront QC on a
subset; (4) hundreds of library reactions invite operator variability (p. 3). FiCS answers with (1) immediate
reversible fixation with DSP (Lomant's reagent), (2) a fixation-compatible viability dye (Zombie) plus sorting for
dCas9 + sgRNA double-positive cells, (3) long-term cryopreservation of fixed cells in Bambanker with RNase inhibitor,
(4) "superloading" up to 100k cells per GEM-X 5′ channel, and (5) library-prep automation on a Hamilton Vantage
(p. 3). Superloading targets about 70,000 recovered cells per lane (p. 16).

### 2. DSP fixation does not degrade scRNA-seq sensitivity (p. 5; Table 1, p. 7; Fig. 2A; Supp. Fig. 1A, p. 22)

«DSP fixation did not interfere with scRNAseq and yielded a highly sensitive transcriptome capture (5,387 and 5,871
median genes per cell on HCT116 and HEK293T respectively, Table 1)» (p. 5). Genes and transcripts per cell «compared
favorably with the benchmarking dataset from Replogle et al. 2022» (p. 5; Fig. 2A). Fresh-vs-fixed agreement is shown
only as a pseudobulk scatter (Supp. Fig. 1A, p. 22); no correlation value is given in the text.

### 3. Cryopreservation decouples harvest from library prep and enables a pilot QC (p. 5; Fig. 2B–C)

With RNase inhibitor, cryostorage «maintained stable RNA-seq quality for up to 140 days in storage» (p. 5; Fig. 2B:
up to 40 days for HCT116 and 140 days for HEK293T, p. 6). Cells were stored at 2–2.5 M per vial and thawed 4–8 M per
batch (p. 5). A pilot QC run used 321,706 of 9,886,825 cells (3.25 %) of the HCT116 screen and 1,084,054 of
21,550,681 (5.03 %) of the HEK293T screen before the full library preparation (p. 5; Fig. 2C).

### 4. Automated library preparation matches manual preparation (p. 5; Fig. 2D)

Pseudobulk counts (log2(cp10k + 1)) from the same cells prepared manually or on the Hamilton Vantage correlate at
Pearson 0.996 (p. 5; Fig. 2D caption, p. 6). The data «were generated using a combination of automated and manual
library preparation methods» (p. 16), so the released atlas mixes both.

### 5. ⭐⭐ Orion is deeper per cell than Replogle K562 and far deeper than Tahoe-100M (p. 5; Table 1, p. 7; Fig. 2E)

Median UMIs per cell are 19,416 (HCT116) and 16,557 (HEK293T); mean reads per cell 54,826 / 63,426; 109 / 223 GEM
groups (Table 1, p. 7). The HCT116 median is «1.68 times as deep as the Replogle K562 Perturb-seq dataset … (median UMIs
per cell = 11499), and 8.45 times as deep as the Tahoe-100M chemical perturbation dataset … (median UMIs per cell =
2298)» (p. 5). Sequencing targeted 50,000 reads per cell for gene expression and 5,000 reads per cell for the CRISPR
library (p. 18). Median CRISPR UMIs per cell differ tenfold between lines: 950 (HCT116) vs 184 (HEK293T) (Table 1).

### 6. ⭐⭐ Control pseudobulks are more reproducible across batches than in Replogle K562 (pp. 5–6; Fig. 2F, p. 7; Supp. Fig. 1E; Methods p. 19)

Pairwise Spearman correlations of pseudobulked non-targeting-control counts between batches have medians 0.993
(HCT116), 0.988 (HEK293T) and 0.967 (Replogle K562) (Fig. 2F caption, p. 7). FiCS is «significantly greater» than
K562 (Wilcoxon rank-sum, p = 5.77e-16 and 2.80e-14; p. 5). MAD of the correlation distribution is 2.87e-3, 3.58e-3 and
9.09e-3 (p. 6). As a negative control, randomly paired HCT116 × Replogle K562 batches give a median Spearman of 0.74
(p. 5; Supp. Fig. 1E). Method: raw counts of expressed genes summed over NTC cells, batches stratified by median UMI
with top and bottom quartiles randomly sampled (p. 19). This is **technical reproducibility of the control profile**,
not of perturbation effects.

### 7. ⭐⭐ About 8 M of about 31 M barcodes carry the intended guide pair, giving a median of 141 / 188 cells per perturbation (pp. 6–7; Table 1)

Barcodes lose the intended single-target state through Poisson lentiviral delivery, superloading multiplets and
sgRNA-calling limits (p. 6). «From approximately 30 million sequenced barcodes, we identified around eight million
cells harboring two sgRNAs targeting the same intended gene, making the "Orion" dataset the biggest publicly deposited
Perturb-seq atlas» (p. 6). Table 1 gives 3,409,169 (HCT116) and 4,534,299 (HEK293T) two-sgRNA cells, a **median of
141 and 188 two-sgRNA cells per perturbation**, plus 1,290,167 / 2,326,883 one-sgRNA cells (p. 7). One-sgRNA cells have
similar knockdown (73.3 % / 45.3 % vs 75.4 % / 51.5 %) but all analyses use only the two-sgRNA cells «to ensure
stringent filtering» (p. 6).

### 8. ⭐ The library targets every protein-coding gene with tRNA-separated dual sgRNAs, without an expression cutoff (pp. 5, 13; Supp. Fig. 1D)

Both lines carry dCas9-KRAB-ZIM3 (CRISPRi), chosen as «more universally applicable» with higher knockdown (p. 5),
stably integrated by PiggyBac as dCas9-mScarlet-ZIM3 and sorted to > 95 % purity (pp. 14–15). Each gene gets one
U6-driven sgRNA pair separated by one of four tRNA sequences (p. 13); guides are read by direct capture (p. 5). Unlike
earlier genome-scale screens, which used an expression threshold, the library includes «pairs of sgRNAs targeting every
gene in the human genome» to obtain empirical data on expression vs effect (p. 13). The tRNA design halves the
guide-pair recombination rate in HCT116 from 26.1 % to 12.6 % (Supp. Fig. 1D, p. 22). Screening: > 1,000-fold coverage
at MOI 0.1 (2.3e8 cells, 2.3e7 transducing units), puromycin days 3–5, harvest on day 7 post-transduction (p. 15).

### 9. ⭐⭐ About 40–48 % of expressed-gene perturbations have a discernible single-cell phenotype (p. 8; Fig. 3A; Supp. Fig. 2; Methods p. 19)

A per-perturbation L1-logistic classifier separates perturbed from a balanced set of NTC cells on held-out cells
(perturbations with ≥ 100 cells, 5,000 most variable genes, 80/20 split, binomial test, Storey q < 0.05; p. 19). «a
considerable proportion of on-target perturbations in expressed genes result in discernible single-cell
transcriptomic changes (48.36% in HCT116, 40.37% in HEK293T, compared to 54.30% in K562 …)» (p. 8). Removing the target
transcript gives similar results (Supp. Fig. 2A). Non-targeting guides as negatives come out discernible at 4.8 %
(HCT116) and 0.0 % (HEK293T) (p. 8; Supp. Fig. 2B–C), which the authors read as possible off-target guides. "Expressed"
genes are those making up 99 % of mean NTC cp10k: 11,181 in HCT116 and 12,573 in HEK293T (p. 19).

### 10. Perturbation-correlation networks recover known complexes (p. 8; Fig. 3B–C; Methods p. 20)

STRING physical-interaction scores rise with the Spearman correlation of perturbation profiles (Fig. 3B), based on
1,574 (HCT116) and 1,944 (HEK293T) strong-phenotype perturbations (≥ 25 cells, ≥ 30 % KD, significant classifier;
p. 20). Clustering the top 100 HCT116 perturbations by effect gives ribosome biogenesis (POLR1G, POLR1E, NOL12, WDR75,
URB1/2, RPL4 …), exosome, Mediator, nuclear pore, eIF2B and aminoacyl-tRNA-synthetase clusters (p. 8). Mitochondrial
proteins (TIMM44, PHB1, SAMM50, PRELID1) cluster with eIF2B, consistent with mitochondrial injury triggering an
integrated-stress-like response (p. 8). The effect score is computed in a 20-PC PCA fitted on NTC cells (p. 20).

### 11. ⭐⭐ Knockdown strength differs between the two lines despite identical vectors (p. 10; Table 1; Fig. 4A–B)

Median on-target knockdown of two-sgRNA cells is **75.4 % in HCT116 vs 51.5 % in HEK293T** (Table 1, p. 7; p. 10).
«As both cell lines received identical dCas9 and sgRNA delivery vectors, we hypothesized this variability arises from
differing expression levels of key CRISPRi components» (p. 10). Baseline sgRNA and TRIM28 levels are higher in HCT116
(fold change 2.26 and 2.61; Mann–Whitney p = 1.61e-32, 3.48e-150; Fig. 4B), suggesting these may be rate-limiting and
«could predict CRISPRi performance across cell systems» (p. 10).

### 12. ⭐⭐ Per-cell sgRNA UMI count is a proxy for per-cell knockdown (p. 10; Fig. 4C–E)

Knockdown varies strongly between cells with the same guide pair (p. 10). sgRNA transcripts are Pol III-driven and
detected at «often hundreds of UMIs per cell» (p. 10). KD efficiency correlates with sgRNA level per cell at R = 0.910
(HCT116, p = 2.54e-4) and 0.901 (HEK293T, p = 3.67e-4; Fig. 4C). The p-values imply a correlation over KD bins
(Fig. 4C plots median sgRNA UMIs per KD bin), not over individual cells. Top- vs bottom-quartile cells by
per-perturbation sgRNA UMIs differ in KD (two-sample KS p = 9.08e-39 / 1.23e-100; Fig. 4D); 15.8 % and 4.2 % of
perturbations show 100 % KD in one group, attributed to sparsity, and were dropped from that plot (p. 12). Keeping the
top 25 % of HEK293T cells by sgRNA count shifts the KD distribution upward (Fig. 4E).

### 13. ⭐ High-sgRNA cells show larger transcriptomic effects: a dose-response (p. 11; Fig. 4F; Methods p. 21)

On a subset of essential-gene perturbations, the high-sgRNA group has a higher energy distance to NTC cells than the
low-sgRNA group: median e-distance fold change 1.66 (HCT116, p = 2.87e-31) and 5.21 (HEK293T, p = 2.86e-167; p. 11).
Guide UMIs are robust-z-scored within GEM group before quartiling (p. 21); e-distance follows Replogle 2022 on 20 PCs
fitted to NTC plus high and low groups (p. 21). «sgRNA abundance is a valuable instrument for stratifying KD
efficiency, filtering out cells with weak perturbations, and enabling the detection of dose-dependent genetic effects»
(p. 11). The larger fold change in the weaker-knockdown line (HEK293T) is not discussed.

### 14. ⭐ Perturbations should be modelled as continuous strengths, which requires many cells per perturbation (pp. 3, 12)

Current analyses and «the training of many biological foundation models, treat genetic perturbations as categorical
variables (e.g., 'on' vs. 'off')» (p. 3). «A key motivation for scaling – specifically, increasing cell coverage per
perturbation – was the shift from viewing CRISPRi outcomes as binary (KD vs. non-KD) to recognizing a spectrum of KD
efficiencies across cells» (p. 12). Dose-aware models may capture haploinsufficiency and compensation (p. 12). No model
is trained or evaluated in the paper, despite the title.

### 15. The fixation principle should generalise beyond DSP (p. 12)

The batch-effect and scaling advantages are attributed to fixation itself, «suggests that similar advantages … could
be realized with other fixation-based single-cell methodologies, such as FLEX (10x Genomics)» (p. 12). This is
asserted, not tested.

### 16. Classifiers detect sparse effects that e-distance misses (p. 19)

The authors defend the classifier as the detectability test: L1 logistic regression quantifies «accuracy and
statistical significance of transcriptional signatures of a perturbation, even when only a small subset of transcripts
are affected», whereas e-distance «average[s] transcriptional effects … across the entire genome» and has «limited
sensitivity to detect sparse or subtle transcriptomic effects» (p. 19). The 48 % / 40 % discernible rates therefore
rest on a sensitive, single-cell, held-out criterion.

## Figures and tables

- Fig. 1 (p. 4): CRISPRi Perturb-seq overview (A), **FiCS workflow (B)**, computational pipeline to one AnnData (C).
- Fig. 2 (pp. 6–7): UMIs vs genes vs Replogle K562 (A), cryostorage stability (B), pilot-QC share (C), automated vs
  manual (D), **median UMIs per GEM group vs other atlases (E)**, **NTC batch-batch Spearman (F)**.
- **Table 1 (p. 7): all QC metrics** for HCT116 and HEK293T (the citable dataset numbers).
- Fig. 3 (pp. 9–10): **classifier accuracy distributions vs K562 (A)**, STRING vs correlation bins (B), cluster
  heatmap (C).
- Fig. 4 (pp. 11–12): **KD CDFs (A)**, sgRNA / TRIM28 / HP1a levels (B), **sgRNA UMIs vs KD bin (C)**, KD by sgRNA
  quartile (D), HEK293T KD after sgRNA filtering (E), **e-distance by sgRNA quartile (F)**.
- Table 2 (pp. 17–18): Hamilton liquid-class performance (automation only).
- Supp. Fig. 1 (p. 22): fresh vs DSP-fixed (A), viability and dCas9/sgRNA gates, 92 % live and 63.7 % double-positive
  (B–C), tRNA recombination 26.1 % → 12.6 % (D), HCT116 × K562 cross-batch correlation (E).
- Supp. Fig. 2 (p. 23): target-ablated classifier (A), NTC classifier controls (B–C).

## What this source does not address

- **No prediction benchmark.** There are no models, metrics, baselines, splits or held-out lines, despite "for
  Training … Foundation Models" in the title (the sharded "model training" on p. 4 is never reported). Cross-line
  prediction on this atlas is evaluated in the harness of this project and, for detectability, in
  [[huuloc-2026-cell-quota]].
- **No reproducibility of perturbation effects.** Reproducibility is shown only for NTC pseudobulks across batches
  (Kernthese 6). There are no replicate screens, no split-half correlations of perturbation responses and no noise
  ceiling. Split-half reliability is the subject of @wang-2026-reliable-perturbations and the variance partition of
  [[molina-2026-response-decomposition]].
- **No cross-line comparison of responses.** HCT116 and HEK293T are compared only on QC and knockdown, never on
  whether the same perturbation produces the same response in both lines.
- **No KO arm, no third line.** CRISPRi only, two transformed lines. Other genome-scale screens (Replogle 2022, Feng
  2024, Nourreddine 2024) are cited only as precedents (p. 2, refs. 17–19).
- **No guide-level efficacy model** beyond sgRNA and TRIM28 abundance; HP1a is measured but its result is not reported
  in the text (p. 10).

<!-- ═══ FOKUS ═══ -->

## Focus: Data / benchmark facts to cite — datasets, cell lines, metrics, baselines and evaluation protocol; and what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 5, 6, 7, 8, 9, 11, 12, 13** (and 1–3 for the Fix-Cryopreserve protocol claims). The hypotheses
are those of the plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3
at the end). The dataset itself is described on [[x-atlas-orion]].

### What this focus needs beyond the Kernthesen

1. **Citable dataset card (all from the source).** Two genome-wide CRISPRi (dCas9-KRAB-ZIM3) Perturb-seq screens,
   HCT116 (colorectal carcinoma) and HEK293T (embryonic kidney), 18,903 protein-coding genes, one dual-sgRNA pair per
   gene plus 1,026 NTC pairs (p. 5); 10x GEM-X 5′ v3 with CRISPR feature barcoding, direct capture (pp. 5, 16);
   Cell Ranger 8.0.1 on GRCh38 2024-A; sgRNA threshold per guide = max(3, lowest UMI of the high component of a
   two-component GMM) (p. 18); NovaSeq X 25B (p. 18). Released on Figshare (doi:10.25452/figshare.plus.29190726) as
   h5ad files **only for cells with two sgRNAs from the same guide pair**; all aligned cells «will be released … at a
   later date» (p. 23). All QC numbers: Table 1 (p. 7).
2. **The source's QC does not agree with the "≈ 85 % below a detection floor" reading of HCT116.** The source's own
   classifier finds 48.36 % (HCT116) and 40.37 % (HEK293T) of expressed-gene perturbations discernible at q < 0.05
   (Table 1, p. 7; p. 8), so about 52 % / 60 % are *not* discernible by that test. [[huuloc-2026-cell-quota]] report
   only 14.7 % / 35.3 % of knockdowns clearing their detection floor (85.3 % / 64.7 % below). The criteria differ: a
   single-cell held-out L1 classifier on 5,000 genes, restricted to expressed genes with ≥ 100 two-sgRNA cells
   (pp. 19), against a bias-corrected pseudobulk magnitude in a fixed embedding. The ordering also flips: the source
   ranks HCT116 *more* detectable than HEK293T, huuloc rank it less. The 4.8 % NTC false-discernible rate in HCT116
   (vs 0.0 % in HEK293T; p. 8) suggests the classifier is somewhat liberal in HCT116. Both analyses agree that a large
   share of genome-wide knockdowns produce no resolvable response. The size of that share depends on the criterion.
3. **Cells per perturbation and the project's panel.** The source's median is 141 two-sgRNA cells per perturbation in
   HCT116 (Table 1). The project's OR_HCT116 harness panel uses 263 official targets × 163 cells with about 12k NTC
   (project fact from `AGENTS.md`, not from the source). The panel is therefore at about the source's median depth,
   restricted to the VCC official targets.
4. **Knockdown comparison across public screens.** The source gives 75.4 % (HCT116) and 51.5 % (HEK293T) median
   on-target KD (Table 1). It does not report a KD for Replogle K562; the project's own measurement (85 % for K562
   GWPS, `AGENTS.md`) is not from this source. The source compares with K562 only on UMIs, genes per cell, NTC batch
   correlation and discernible fraction (pp. 5–8).
5. **Dose heterogeneity is a within-perturbation variance source.** Per-cell KD varies with sgRNA abundance
   (Kernthesen 12–13), and the top vs bottom sgRNA quartiles differ in effect size by 1.66× (HCT116) to 5.21× (HEK293T)
   in e-distance (p. 11). *Reader's inference, not in the source:* random cell halves draw from the same KD
   distribution, so a split-half ceiling is unbiased by this in expectation. However, the within-perturbation variance
   it adds lowers the reliability of each perturbation mean at fixed cell count, and line-level KD differences
   (75 % vs 52 %) put a strength factor into the cross-line response that looks like a line × perturbation
   interaction.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partial (noise-ceiling inputs only)** | Supplies the inputs that set the size of the OR_HCT116X ceiling: depth (median 19,416 UMIs), cells per perturbation (median 141), knockdown (75.4 %), the discernible share (48.36 %) and NTC batch reproducibility (Spearman 0.993). No split-half reliability of perturbation effects, no replicate screen, no r_max, no prediction. Many non-discernible targets mean a small Var_reproducible on part of the panel, and the ceiling there will be low and noisy. | Table 1 p. 7; pp. 5–8; Fig. 2F p. 7; Fig. 3A |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **none** | Tahoe-100M appears only as a depth comparator (median 2,298 UMIs per cell). | p. 5; Fig. 2E |
| **H2b** (theory predicts the break-even within ×2) | **none** | No theory, no d, σ²_ε/Var(i) or R². | — |
| **H2c** (genetic C*_min > C_public) | **partial (registry input)** | Contributes **two** of the project's 11 Count-1 cell types (HCT116, HEK293T; `output/data-analysis/results/families.json`). Both are genome-wide CRISPRi with ≥ 100 cells per perturbation at the median, so both clearly qualify. It claims to be «the largest publicly available Perturb-seq atlas» (p. 1) and cites Replogle 2022, Feng 2024 (HipSci iPSC) and Nourreddine 2024 (iPSC) as earlier genome-scale screens (p. 2, refs. 17–19). It does not count public lines and has no C*_min. The KD gap between its two lines (75 % vs 52 %) is a reason to treat CRISPRi efficiency as part of b_c / i in σ²_ε/Var(i). | pp. 1–2, 5, 7, 10 |
| **H3a** (transfer residual on destination top-8 control PCs vs other lines / permuted, split-half) | **none (method echo only)** | PCA is fitted on NTC cells (20 PCs) only to score effect size and e-distance. There is no transfer residual and no subspace comparison. | pp. 20–21 |
| **H3b** (M3 beats five PC-deflation / common-response corrections) | **none** | No correction methods. | — |
| **H4** (adding same-line genome-wide data beats the largest modelling gain) | **none (framing only)** | The introduction motivates data generation because foundation models underperform simpler models (refs. 9–11, p. 2). There is no experiment on adding same-line data and no model. Kernthese 14 ("data at scale per perturbation") is a design argument, not a test. | pp. 2–3, 12 |

**Positioning consequence.** Nothing here moves C*_emp ≈ 37 or the rejected linear-Gaussian C* formula. The source has
no context axis. Its relevance is (i) the H1 noise ceiling on OR_HCT116X, where the low discernible share and the dose
heterogeneity predict a low, uncertain r_max, and (ii) two entries in C_public.

### Direct quotes (for the draft)

> "This atlas, generated from two genome-wide FiCS Perturb-seq experiments targeting all human protein-coding genes,
> comprises eight million cells deeply sequenced to over 16,000 unique molecular identifiers (UMIs) per cell." (p. 1)

> "Moreover, FiCS Perturb-seq yielded more consistent data with significantly greater batch-to-batch correlation than
> that of Replogle K562 (Wilcoxon rank-sum test, 5.77e-16 and 2.80e-14 for HCT116 and HEK293T vs. Replogle K562
> respectively, Methods, Figure 2F)." (p. 5)

> "Our analysis indicates that a considerable proportion of on-target perturbations in expressed genes result in
> discernible single-cell transcriptomic changes (48.36% in HCT116, 40.37% in HEK293T, compared to 54.30% in K562 17,
> Figure 3A)." (p. 8)

> "Perturbed cells exhibited notable differences in median KD efficiency, the degree to which expression of a target
> gene is inhibited by the CRISPRi system, between cell lines: 75.4% in HCT116 versus 51.5% in HEK293T (Table 1,
> Figure 4A). As both cell lines received identical dCas9 and sgRNA delivery vectors, we hypothesized this variability
> arises from differing expression levels of key CRISPRi components." (p. 10)

> "Effectively interrogating this spectrum requires a substantial number of cells per perturbation to achieve adequate
> statistical power and resolution." (p. 12)

### Examples & illustrations (for later drafting)

- Same vectors, different lines: HCT116 75.4 % vs HEK293T 51.5 % median knockdown, with sgRNA and TRIM28 2.26× and
  2.61× higher in HCT116 (p. 10; Fig. 4A–B). This is a concrete case of a line-level factor that scales every
  perturbation's response.
- Dose within one perturbation: top- vs bottom-quartile sgRNA cells differ in e-distance to NTC by 1.66× (HCT116) and
  5.21× (HEK293T) on essential genes (p. 11; Fig. 4F).
- Mitochondrial perturbations (TIMM44, PHB1, SAMM50, PRELID1) cluster with eIF2B members, read as an integrated-stress
  phenotype (p. 8; Fig. 3C).

### Boundary of this focus

This focus records the dataset and QC facts and their bearing on H1–H4. It does not work up the automation hardware,
lentivirus production or thawing protocol (pp. 14–18). A reader looking for prediction metrics, baselines, a
cross-line evaluation protocol or a reproducibility ceiling for perturbation effects will not find them in this
source.

## Mentioned entities

- Datasets: [[x-atlas-orion]], [[replogle-nadig-crispri-collection]] (Replogle K562 genome-wide as comparator),
  [[tahoe-100m]] (depth comparator)
- Concepts: [[simple-perturbation-baselines]] (motivating citation only)

## Connections

- Introduces / supports: [[x-atlas-orion]] (pp. 1, 5–7, 23).
- Cites as benchmark: [[replogle-nadig-crispri-collection]] (Replogle K562; pp. 5–8) and [[tahoe-100m]] (depth only;
  p. 5).
- Supports (framing citation): [[simple-perturbation-baselines]] (p. 2, refs. 9–11).
- Ambiguous disagreement on detectability: [[huuloc-2026-cell-quota]] re-analyse Orion and find far fewer detectable
  knockdowns in HCT116 (14.7 %) than this source's classifier (48.36 %); criteria differ (focus item 2).
- Cites (not yet ingested): @replogle-2022-genome-scale (ref. 17), @zhang-2026-tahoe (ref. 31, preprint version),
  Feng 2024 and Nourreddine 2024 (refs. 18–19), Song 2025 (perturbation-response score, ref. 20).
