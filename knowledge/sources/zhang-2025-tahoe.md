---
title: "Zhang et al. 2025 — Tahoe-100M preprint (giga-scale single-cell drug perturbation atlas)"
type: source
created: 2026-11-06
updated: 2026-11-06
status: review
author: llm
bibkey: zhang-2025-tahoe
based_on: preprint
depth: standard
tags: [dataset, tahoe, chemical-perturbation, single-cell, cell-village, mosaic, demultiplexing, context-dependence, preprint]
relations:
  - target: tahoe-100m
    type: supports
    confidence: extracted
    because: "Primary description of the dataset: «Tahoe-100M comprises over 100 million single-cell transcriptomes collected from 50 cancer cell lines subjected to more than 1,100 small-molecule perturbations at varied doses» (p. 2); 47 lines, 379 drugs, 1,135 drug-dose combinations after QC (p. 7)."
  - target: replogle-nadig-crispri-collection
    type: cites
    confidence: extracted
    because: "Uses Replogle 2022 genome-wide CRISPRi as an E-distance benchmark: «when examining genome-wide CRISPRi perturbations (Replogle et al. 2022) … we found drug perturbations across this dataset to exhibit an intermediate effect» (p. 9; Fig. 4B, p. 10)."
---

# Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling

> Preprint version (zhang-2025-tahoe, bioRxiv v3, posted 2025-05-10) of zhang-2026-tahoe; claims quoted from the preprint.

> [!warning] Provenance: based on the **preprint**, not the published original.
> Pagination and final claims may differ. Re-ingest the original when acquired
> (set `based_on: original`).

## Bibliographic Details

Zhang, Jesse, Airol A. Ubas, Richard de Borja, Valentine Svensson, Nicole Thomas, et al. (32 authors; Tahoe
Therapeutics and Parse Biosciences; corresponding J. Yu and H. Goodarzi). 2025. "Tahoe-100M: A Giga-Scale Single-Cell
Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling." bioRxiv preprint
doi:10.1101/2025.02.20.639398, **version 3, posted 10 May 2025** (v1 2025-02-24, v2 2025-04-23; bioRxiv API). Licence:
«All rights reserved. No reuse allowed without permission» (header, every page), so **no figures are extracted** to
`knowledge/assets/`. Data: «Data is available on HuggingFace» (p. 23). Read from
`<library>/pdf/zhang-2025-tahoe.pdf` (27 pp.). Citekey: @zhang-2025-tahoe (preprint of @zhang-2026-tahoe, Cell
189(19):5945–5961.e8, 2026).

**Page numbers.** The preprint carries no printed page numbers. Every "p." below is the **PDF page**. Figures sit on
their own pages without extractable text, with the caption on the next page: Fig. 1 image p. 4 / caption p. 5; Fig. 2
p. 6 / p. 7; Fig. 3 p. 8 / p. 9; Fig. 4 p. 10 / p. 11; Fig. 5 p. 12 / p. 13; Fig. 6 p. 15 / p. 16. Table 1 is on p. 5.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; Introduction (motivation, Mosaic, scale claim) | 1–2 |
| ✓ F1 | Results: Generating the Tahoe-100M atlas using the Mosaic platform (Fig. 1; Table 1) | 3–5 |
| ✓ F1 | Results: Cell line and treatment diversity (Fig. 2) | 5–7 |
| ✓ F1 | Results: The global transcriptomic landscape (scVI, tSNE, plate-replicate correlation; Fig. 3) | 7–9 |
| ✓ F1 | Results: Drug-induced transcriptomic effects across cell lines and MOAs (LISI, E-distance; Fig. 4) | 9–11 |
| — | Results: gene-set-score embedding (tSNE, nnMDS; Fig. 5) | 11–13 |
| ✓ F1 | Results: Context-dependent inhibition of the RAS/RAF pathway; cell-cycle effects (Fig. 6) | 13–16 |
| ✓ F1 | Discussion (AI-model claims, heterogeneity, translational, future directions) | 16–17 |
| ✓ F1 | Methods: cell culture, compounds, dissociation, barcoding, sublibraries, sequencing, split-pipe | 18–19 |
| ✓ F1 | Methods: genotype-based demultiplexing; QC filters; replication analysis | 19–20 |
| — | Methods: drug targets and MOA (GPT-4o annotation) | 20 |
| ✓ F1 | Methods: scVI model; LISI; visualisation | 20–21 |
| — | Methods: differential gene-set scores; nnMDS | 21–22 |
| ✓ F1 | Methods: E-distance; cell-cycle phase analysis; gene-set signature analysis | 22–23 |
| — | Data availability; acknowledgments; references | 23–27 |

## Kernthesen

### 1. ⭐⭐ Tahoe-100M is the largest single-cell perturbation dataset: 100 M cells, 50 lines, > 1,100 drug perturbations (pp. 1–2, 16)

«a giga-scale single-cell atlas of 100 million transcriptomic profiles measuring how each of 1,100 small-molecule
perturbations impact cells across 50 cancer cell lines» (p. 1). The Introduction words it as «more than 1,100
small-molecule perturbations at varied doses» (p. 2), and the Discussion as «more than 1,100 small-molecule agents»
(p. 16). The "1,100" counts drug × dose combinations (1,135 after QC, p. 7), not distinct drugs (379, p. 5); the
abstract's wording blurs this (see "Inconsistencies").

### 2. ⭐⭐ The Mosaic "cell village" design removes most batch effects by pooling lines in one well (pp. 2–3, 5)

Tens to hundreds of cell models are pooled and then genetically demultiplexed, which «substantially reduces batch
effects while scaling to an increasingly larger number of perturbation sets» (p. 2). Spheroids that mix all lines in
suspension each receive one treatment, including DMSO vehicle controls, for 24 h. They are then dissociated, fixed and
profiled with the Parse GigaLab combinatorial-barcoding kit (p. 3). Perturbation identity comes from the first barcode
round; line identity comes «by the SNPs present in the cDNA insert» (Fig. 1B legend, p. 5).

### 3. ⭐ Cell-line identity is assigned by SNP demultiplexing with > 98 % accuracy (p. 19)

demuxlet is run against a curated SNP reference: exonic/UTR dbSNP variants with gnomAD AF > 0.1, genotyped per line
with mpileup, filtered to DP > 6, MQ > 30, QUAL > 50, with indels, multiallelic and all-reference SNPs removed (p. 19).
On cells from known lines, «accuracy exceeding 98%, with more than 90% of cells being labeled as singlets» (p. 19).
Only demuxlet singlets enter both filter levels (p. 19).

### 4. ⭐ Production scale: 14 plates, 1,786 sublibraries, about 1.4 trillion reads, 153 M cells (p. 3; Table 1, p. 5)

«a total of 14 96-well plates of spheroids were seeded and drug treated» (p. 3). Each plate gave 72–165 sublibraries,
1,786 in total. About 1.4 trillion raw reads covered 153 million cells (8.2–14.5 M per plate). Cells had «a mean of
2,288 transcripts (median of 1,890)» (p. 3). 100.6 M cells pass minimal filters and 95.6 M pass full filters, and
«We used the latter set for all the downstream analyses» (p. 3). Per-plate median transcripts per cell range from
1,462 (plate 3) to 2,159 (plate 6) (Table 1, p. 5).

### 5. ⭐⭐ After QC: 47 lines, 379 drugs, 1,135 drug-dose combinations, 52,886 conditions (pp. 5, 7, 19)

«Of the 50 cell lines in Tahoe-100M, 47 had sufficient representation across the conditions to be evaluated in
downstream analyses» (p. 5). The excluded three are named in Methods: «three cell lines (NCI-H661, NCI-H596,
NCI-H2122) that consistently exhibited low counts across conditions were excluded» (p. 19; red in Fig. 3E, p. 8/9).
Condition-level QC: «each cell line-drug treatment condition was required to have a minimum of 50 cells» (p. 19). The
result is «47 cell lines, 379 drugs, 1,135 drug-dose combinations and 52,886 unique cell line-drug-dose conditions,
with a median of 1,287 cells per condition» (p. 7). Full filters: ≥ 700 UMIs, < 20 % mitochondrial reads, UMI and
mito z-scores within ±3, ≥ 250 genes, demuxlet singlet (p. 19). The three doses are 0.05, 0.5 and 5 µM (Fig. 4B, p. 10;
Fig. 6 legend, p. 16: «all 3 drug concentrations»).

### 6. The panel is diverse in organ, driver mutation and drug mechanism (p. 5; Fig. 2, pp. 6–7)

The 47 lines come from 13 organs (mostly lung, bowel, pancreas, skin). TP53, KRAS and CDKN2A are altered in about half
of them (p. 5). Of the 379 drugs, 180 fall into 25 MOAs (median 5 drugs per MOA, range 3–27). 69 % are approved agents,
and they reportedly target 325 genes, 120 of them by more than one drug (p. 5). Targets, MOAs and approval status were
annotated with GPT-4o against MedChemExpress text, with manual inspection of 97 of 380 compounds (p. 20).

### 7. Tahoe is much larger than earlier screens per condition and per drug count (p. 5; Fig. 2D)

«This dataset captures 17,813 unique cell line-drug conditions, representing a 31-fold increase in the number of
drug-perturbed conditions compared to a benchmark single-cell drug screen (Srivatsan et al. 2020) and a 29-fold increase
in the number of observations per condition compared to a benchmark genetic perturbation screen (Replogle et al.
2022)» (p. 5).

### 8. ⭐⭐ Cells separate by line and cell cycle, not by plate, so batch effects are not substantive (pp. 7, 9; Figs. 3A–C, 4A)

A 10-dimensional scVI embedding of 140,000 cells from each of 47 lines (6.58 M cells) gives tSNE clusters that follow
«their genetic identity (and cell cycle phase) rather than their plates of origin … suggesting a lack of substantive
batch effects» (p. 7). The authors still control for plate: «the differential gene expression analyses described in
the following sections account for the possible batch effect by leveraging the DMSO-treated wells included in all
plates» (p. 7). LISI over 100 subsamples of 10,000 cells (30 neighbours; p. 21) ranks the factors: «cell line
identity, followed by cell cycle phase, are the most major factors in driving stratification of the data. And as
expected, drug treatment and dose have a more modest effect» (p. 9). Read off the bar chart (Fig. 4A, p. 10), median
LISI is about 1 for cell line, about 1.3 for cell cycle, and about 6 for drug and for drug × concentration.

### 9. ⭐⭐ Plate 14 is a biological replicate of plate 6; matched pseudobulks correlate at r ≈ 0.975 (pp. 7, 19–20; Fig. 3D)

«To highlight the reproducibility of the Mosaic platform, we included Plate 14 as a biological replicate of Plate 6»
(p. 7). Pseudobulks per (line, treatment) gave 4,800 profiles on plate 6 and 4,796 on plate 14, restricted to the 4,796
common pairs (p. 19). Pearson r on log(UMI + 1) for 4,000 sampled pairs per class: matched treatment and identity q25
0.97 / median 0.975 / q75 0.98; unmatched q25 0.892 / median 0.915 / q75 0.929 (p. 7; method p. 20). This is the only
designed replicate in the atlas. Plate 14 is also held out from scVI training «for validation and model criticism»
(p. 21).

### 10. ⭐ Drug effect size grows with dose and sits between non-essential and essential CRISPRi (p. 9; Fig. 4B)

E-distance (computed in the scVI latent space against same-plate DMSO_TF cells, with the median taken over plates and
then over lines; p. 22) is larger at the highest dose. Its distribution resembles sci-Plex3. Against genome-wide
CRISPRi split by essentiality, «drug perturbations across this dataset … exhibit an intermediate effect» (p. 9).

### 11. Some MOAs and outlier drugs dominate effect size (p. 9; Fig. 4C)

Harringtonine and homoharringtonine (direct ribosome-elongation inhibitors) and dinaciclib (CDK1/2/5/9) are outliers.
Proteasome, HDAC and PI3K/AKT inhibitors have larger effects across doses. Drugs active in epithelial solid tumours
(MAPK, PI3K/AKT, MTOR and CDK inhibitors, microtubule agents, DNA synthesis/repair inhibitors) exceed retinoid,
adrenoceptor, COX inhibitors and unclassified compounds (p. 9).

### 12. A gene-set-score embedding shows an E-distance gradient and genotype-dependent separation for RAS inhibitors (p. 11; Fig. 5)

Each (line, drug, dose) pseudobulk is scored on 7,390 MSigDB gene sets (p. 21) and embedded by tSNE and nnMDS. Both
embeddings show «a clear E-distance gradient from the center to the periphery» (p. 11). High-effect MOAs overlap but
partly separate, and adagrasib in KRAS G12C lines separates from adagrasib in non-G12C lines, alongside RMC-6236 (p. 11).

### 13. ⭐⭐ RAS/RAF inhibitor responses depend on the line's driver genotype (p. 13; Fig. 6A–B)

On a 600-gene RAS/RAF signature (KRAS.600.LUNG.BREAST.V1; p. 23), dabrafenib (BRAF V600E inhibitor) «elicited the
most pronounced effect in cell lines harboring BRAF mutations, as evidenced by significantly lower scores at both the
lowest (P = 0.0013, Welch's t-test) and highest (P = 0.00023, Welch's t-test) doses, whereas KRAS-mutated lines showed
minimal change» (p. 13). The pan-RAS inhibitor RMC-6236 had its largest effect in KRAS-mutated lines, at p = 0.065
(lowest dose) and p = 0.063 (intermediate dose) against RAS-unmutated lines (p. 13). Those two p-values do not reach
0.05, but the text still reads them as underscoring «the KRAS-specific action of RMC-6236» (p. 13).

### 14. Cell-cycle shifts are drug-class and drug specific (pp. 13–14; Fig. 6C–E)

The per-line, per-plate log2 odds ratio of phase fractions against DMSO_TF (p. 23) gives these results: palbociclib
enriches G1 and dinaciclib enriches G2/M, while abemaciclib, AT5719 and ribociclib show little change. HDAC and
microtubule inhibitors give the largest G2/M arrest (belinostat, panobinostat). Carbamazepine shows little change and
tucidinostat negative enrichment, explained by weak HDAC activity and CDK4 down-regulation respectively. Paclitaxel and
tubulin inhibitor 6 are exceptions among the microtubule agents (pp. 13–14).

### 15. A single scVI model represents the whole atlas (pp. 20–21)

scVI with a negative-binomial likelihood models all 62,710 genes in 10 dimensions. It was trained on 89,423,257
full-filter cells from plates 1–13 (90/10 train/validation split, 10 epochs, «over 800 million examples»), with plate 14
excluded (pp. 20–21). The minified model covers 95,624,334 cells in 41 GB (p. 21). This is the only model the preprint
trains. It is used for embedding, E-distance and gene-set scores, **not** for response prediction.

### 16. ⭐ The authors claim the atlas will power context-aware AI models, by analogy to scaling in NLP and vision (pp. 2, 16–17)

The goal is «to catalyze the development of robust, context-aware AI models of cellular behavior» (p. 2). Broad,
balanced conditions mean that «AI models see complex, context-dependent gene regulatory modulations rather than uniform
or highly specialized responses» (p. 16). The «cell village» framework lets AI systems «discern how genetic and
epigenetic contexts modulate drug efficacy and cell-state transitions» (p. 16). The scaling expectation is stated as an
analogy, not tested: «The success of AI frameworks in fields like natural language processing and image recognition
highlights how performance can scale dramatically with increased data size and diversity. We anticipate a similar
trajectory for single-cell biology» (p. 16).

### 17. Drug response is heterogeneous within and between lines, and models must capture distributions (pp. 16–17)

The atlas «exposes the profound heterogeneity in how cells respond to identical pharmacological perturbations» (p. 16).
This covers subpopulation substates within a line, driver-mutation-dependent pathway effects, and transcriptional
correlates usable as biomarkers (pp. 16–17). AI models of drug response «must account for the full distribution of
potential outcomes rather than an "average" or bulk response» (p. 17).

### 18. Every differential analysis is referenced to same-plate vehicle controls (pp. 7, 18, 21–23)

Wells H11 and H12 of each drug plate are «control media with no compounds» (p. 18). Differential Vision scores compare
«drug-treated cells compared against vehicle (DMSO)-treated controls within the same plate», per plate and line,
dropping comparisons with < 50 cells on either side (62,927 kept; p. 21). E-distance (p. 22) and cell-cycle odds
ratios («To address potential plate batch effects, the score was computed on a per cell line, per plate basis», p. 23)
use the same per-plate reference.

## Figures and tables

- Fig. 1 (image p. 4, caption p. 5): Mosaic workflow; read structure (perturbation barcode + SNP line identity); QC
  pipeline.
- Table 1 (p. 5): per-plate sublibraries, total cells, full-filter cells, median and mean transcripts per cell.
- Fig. 2 (p. 6 / 7): oncoplot of 47 lines; MOA and target-gene groupings of the 379 drugs; **scale comparison with other
  datasets (D)**.
- Fig. 3 (p. 8 / 9): tSNE by line, plate and cell cycle; **plate 6 vs plate 14 replicate correlation (D)**; **coverage
  heatmap and per-line cell counts, three low-abundance lines in red (E)**.
- Fig. 4 (p. 10 / 11): **LISI by factor (A)**; E-distance by dose against sci-Plex3 and Replogle (B); by MOA (C).
- Fig. 5 (p. 12 / 13): gene-set-score tSNE and nnMDS; RAS inhibitors by KRAS G12C status.
- Fig. 6 (p. 15 / 16): **dabrafenib and RMC-6236 by genotype (A–B)**; cell-cycle log2 odds ratios (C–E).

**Inconsistencies in the source.** (i) "1,100 small-molecule perturbations" (p. 1) vs 379 drugs and 1,135 drug-dose
combinations (pp. 5, 7). (ii) 379 drugs (p. 5) vs «97 of the 380 compounds» (p. 20). (iii) 61,927 pseudobulk groups
(Fig. 5 legend, p. 13) vs 62,927 comparisons (pp. 21–22). (iv) The controls are called «DMSO vehicle controls» (p. 3)
and «DMSO_TF» (p. 22), but the plate layout calls H11/H12 «control media with no compounds» (p. 18). (v) The Fig. 3D
panel title reads "Mosaic Assay Correlation from Plate Controls" (p. 8), while text and Methods describe all
(line, treatment) pseudobulks (pp. 7, 19–20). (vi) Three different condition counts appear: 17,813 line–drug (p. 5),
52,886 line–drug–dose (p. 7) and 62,927 plate–line–drug–dose comparisons (p. 21).

## What this source does not address

- **No prediction model and no benchmark.** The preprint trains only an scVI representation (Kernthese 15). It reports
  no response-prediction task, no held-out-line evaluation, no baselines and no metric for predictive accuracy. The
  claim that data scale and diversity will improve models (p. 16) is an expectation, not a result. Benchmarks on Tahoe
  are in [[molina-2026-response-decomposition]] (5 lines) and in the project's own curve ([[context-count-learning-curve]]).
- **No variance decomposition and no noise model.** There is no partition into drug, line and line × drug interaction
  components, no measurement-noise estimate, and no split-half or replicate-based noise variance. The only
  reproducibility number is the plate 6 vs plate 14 pseudobulk correlation (Kernthese 9). Noise regimes are treated in
  [[huuloc-2026-cell-quota]]; the decomposition in [[perturbation-response-decomposition]].
- **No context-count analysis.** The number of lines is never varied, and nothing is said about how many contexts a
  model needs.
- **No replicate structure beyond one duplicated plate**, and no per-condition plate counts. It says only that «For cell
  line–drug conditions with multiple plates, a median was taken across all plates» (p. 22).
- **Context dependence is shown qualitatively**: two RAS/RAF inhibitors, one KRAS G12C contrast, and cell-cycle odds
  ratios. There is no atlas-wide test of line × drug interaction.
- **No genetic perturbations.** CRISPRi appears only as an E-distance comparator (p. 9).

<!-- ═══ FOKUS ═══ -->

## Focus: facts and claims the manuscript cites; exactly what overlaps with H1–H4 — 2026-11-06

**Carries Kernthesen 1, 2, 3, 4, 5, 8, 9, 13, 15, 16, 18** (and 10, 17 for framing). The hypotheses are those of
`input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 plus Amendments A1–A3). Verdicts are from
[[results-summary]].

### What this focus needs beyond the Kernthesen

1. **The 48-line drug data (H2a/H2b) against the preprint's 47.** The preprint keeps 47 of 50 lines and excludes
   NCI-H661, NCI-H596 and NCI-H2122 for low counts (p. 19). The project's 48 Tahoe contexts
   (`vcc2026_data/data/processed/drug_context_pretrain/contexts.npz`, `names`) contain **NCI-H661** (`tahoe:NCIH661`)
   and contain neither NCI-H596 nor NCI-H2122. So the project's 48 = the preprint's 47 + NCI-H661. This is a reader's
   derivation from the project's own file, not a statement of the preprint, and it assumes `HEPG2_FAMILY` is the
   atlas's HepG2/C3A. The manuscript should say "48 of the 50 Tahoe-100M lines; the original analysis retained 47 and
   excluded NCI-H661 for low counts [@zhang-2025-tahoe, p. 19], which we keep" (or drop it, for a sensitivity check).
2. **Cite-able design facts for the Methods.** 50 lines pooled per well (p. 3; Fig. 1A, p. 5); SNP demultiplexing with
   demuxlet, > 98 % accuracy (p. 19); 24 h treatment (pp. 3, 18); 14 96-well plates (p. 3); three doses 0.05 / 0.5 /
   5 µM (Fig. 4B, p. 10); 379 drugs, 1,135 drug-dose combinations (p. 7); 95.6 M full-filter cells (p. 3); ≥ 50 cells
   per condition in the original QC (p. 19); same-plate vehicle controls, two wells per plate (H11, H12; p. 18). The
   project's 5 µM subset (185 drugs, 4.43 M cells, 13,011 conditions) is a project construction. The preprint gives no
   per-dose drug count to check it against.
3. **Noise / interaction ratio (Tahoe 2.87).** **Not in the preprint.** It has no Var(i), no σ²_ε and no ratio, so the
   2.87 must be cited as the project's own estimate. Two preprint facts give qualitative context only:
   - Drug and dose stratify cells far less than line identity and cell cycle (LISI, p. 9; Fig. 4A).
   - Unmatched pseudobulks from the two replicate plates still correlate at a median of 0.915, against 0.975 for matched
     ones (p. 7).

   Reader's inference, not the source's: most pseudobulk log-expression variance is shared baseline, and the
   drug-specific signal is a small residual. Neither number can be turned into σ²_ε or Var(i). The correlation is on
   raw log-UMI pseudobulks, not on responses.
4. **Replicate structure.** The only designed replicate is plate 14 = plate 6 (pp. 7, 19, 21). This fits the project's
   choices: cell-level split halves for σ²_ε, and the use of Tahoe "primary vs repeat" in [[appx-drug-gene-bound]].
   The descriptor pools DMSO_TF from plates 3/6/14 ([[rq2_tahoe_curve-assumptions]]). The preprint does not say why
   plate 3 carries those lines' controls.
5. **Context dependence of drug responses.** The preprint's evidence is qualitative and genotype-specific: dabrafenib
   in BRAF V600E vs KRAS lines (P = 0.0013 / 0.00023), RMC-6236 in KRAS lines (p ≈ 0.06), adagrasib by KRAS G12C, and
   CDK inhibitor cell-cycle profiles (pp. 11, 13–14). This supports the *existence* of a line × drug interaction i,
   which motivates a learned context correction. It says nothing about its size relative to noise.
6. **Modelling benchmarks.** **None reported** (Kernthese 15). The citable modelling claim is the scaling expectation
   (p. 16), which the project's H2a result (context-pathway gain negative at C = 12; C*_emp = 37.4 > 24) puts in
   perspective. The gain from more contexts is real but small, and it only appears after about 37 lines. That
   contrast is the project's argument; the preprint makes no quantitative claim to refute.

### Overlap with H1–H4

| Hyp. | What the source shows | What we add | Overlap |
|---|---|---|---|
| **H1** (pipeline ≥ 0.7 of r_max) | Nothing on ceilings or noise-corrected accuracy. The replicate correlation (0.975, p. 7) is a raw-expression reproducibility figure, not a response ceiling. | Variance components, r_max, attained fraction (genetic panels; untestable as registered). | **none** |
| **H2a** (Tahoe: context gain < 0 at C = 12; C*_emp ≥ 24) | Supplies the data (50 lines, 24 h, 3 doses, within-plate controls; pp. 3, 5, 7, 18) and the qualitative claim that line context modulates drug response (p. 13). It never varies the number of lines and benchmarks no model. | The seven-point context curve with matched constant-descriptor ablation, cluster bootstrap, C*_emp = 37.4 (confirmed). | **none (data provider only)** |
| **H2b** (Prop. 2 predicts C*_emp within ×2) | No d, σ²_ε, Var(i) or R². No noise model; LISI (p. 9) is qualitative. | Blind parameter estimates; C*_theory 230 vs 37.4 (rejected). | **none** |
| **H2c** (genetic C*_min > C_public) | Not genetic. CRISPRi appears only as an E-distance comparator (p. 9). | Genetic C*_min and C_public counts. | **none** |
| **H3a** (transfer residual on destination control PCs) | Line identity dominates expression variation (tSNE, LISI; pp. 7, 9). This is consistent with control-state structure being line-specific, but there is no PCA of controls and no residual. | The whole residual analysis (confirmed on the genetic H1 panel). | **none (background only)** |
| **H3b** (M3 vs five corrections) | Uses within-plate DMSO referencing only as an analysis choice (pp. 21–23); no correction is compared. | The correction comparison (inconclusive). | **none** |
| **H4** (same-line genome-wide data beats modelling) | Not addressed. The Discussion's call for more, more diverse data (p. 16) is rhetorically aligned but is about scale, not same-line data. | H1GW vs modelling gains on the real scorer (confirmed). | **none** |

### Direct quotes (for the draft)

> "Tahoe-100M comprises over 100 million single-cell transcriptomes collected from 50 cancer cell lines subjected to
> more than 1,100 small-molecule perturbations at varied doses." (p. 2)

> "Of the 50 cell lines in Tahoe-100M, 47 had sufficient representation across the conditions to be evaluated in
> downstream analyses." (p. 5)

> "three cell lines (NCI-H661, NCI-H596, NCI-H2122) that consistently exhibited low counts across conditions were
> excluded from subsequent analyses." (p. 19)

> "Taken together, after applying our QC filtering, we obtained a dataset comprising 47 cell lines, 379 drugs, 1,135
> drug-dose combinations and 52,886 unique cell line-drug-dose conditions, with a median of 1,287 cells per
> condition." (p. 7)

> "The performance evaluation demonstrated an accuracy exceeding 98%, with more than 90% of cells being labeled as
> singlets." (p. 19)

> "cell line identity, followed by cell cycle phase, are the most major factors in driving stratification of the data.
> And as expected, drug treatment and dose have a more modest effect relative to these factors" (p. 9)

> "Our results show that targeted inhibitors, such as those directed against RAS or RAF, elicit markedly different
> transcriptomic shifts in a driver mutation-dependent manner." (p. 17)

> "The success of AI frameworks in fields like natural language processing and image recognition highlights how
> performance can scale dramatically with increased data size and diversity. We anticipate a similar trajectory for
> single-cell biology as Tahoe-100M provides a critical leap in the availability of meaningful training data." (p. 16)

### Examples & illustrations (for later drafting)

- Dabrafenib vs genotype: BRAF V600E lines respond and KRAS lines barely change (P = 0.0013 lowest dose, 0.00023 highest;
  p. 13). This is the cleanest published example of a line × drug interaction in Tahoe.
- Plate 6 vs plate 14: matched median r 0.975 vs unmatched 0.915 (p. 7). Use it to illustrate that shared baseline
  dominates raw pseudobulks, which is why responses (differences from controls) are noise-dominated.
- LISI: line about 1, drug about 6 (Fig. 4A, p. 10). Drug identity barely structures the embedding.

### Boundary of this focus

This focus does not work up the gene-set embedding (Fig. 5), the GPT-4o annotation pipeline or the per-MOA biology,
beyond their Kernthesen. A reader looking for a modelling benchmark, a noise estimate or a context-count result will
not find one in this preprint.

## Preprint vs published version

- **Version history (bioRxiv API, 2026-11-06):** v1 2025-02-24, v2 2025-04-23, v3 2025-05-10. All three have the same
  title, category genomics and licence `cc_no`. Only v3 was read. **Differences between v1, v2 and v3: not audited.**
- **Published:** Cell 189(19):5945–5961.e8 (17 Sep 2026; PMID 42753697), doi:10.1016/j.cell.2026.08.035, retitled
  "Tahoe-100M: Mapping drug-induced molecular phenotypes at single-cell resolution". It has 40 authors (preprint 32),
  with Arc Institute added as an affiliation. These differences come from PubMed and Crossref metadata and are noted
  here, never quoted.
- **Abstract (PubMed, compared only to note differences):**
  - The count is reworded from "1,100 small-molecule perturbations" to 1,100 drug-dose conditions.
  - The published abstract adds systematic phenotype readouts (proliferation, cytotoxicity, lineage-specific
    vulnerabilities, cell-cycle changes), population heterogeneity (divergent vs convergent fates), pathway-based MOA
    classification, off-target activities and adaptive stress responses linked to resistance.
  - The published version therefore has substantial analyses that are absent from this preprint. Whether it adds a
    variance decomposition, a modelling benchmark or a changed line count is **not audited**.
- **Consequence:** the wiki's earlier "50 lines × 1,100 drug-dose conditions" wording (literaturguide) follows the
  published abstract. When citing @zhang-2025-tahoe, use the preprint's numbers.

## Mentioned entities

- Datasets: [[tahoe-100m]], [[replogle-nadig-crispri-collection]] (E-distance comparator)
- Concepts (absent from this source; linked for orientation): [[perturbation-response-decomposition]],
  [[context-count-learning-curve]]

## Connections

- Supports / is the primary description of: [[tahoe-100m]] (pp. 1–7, 18–19).
- Cites as comparator: [[replogle-nadig-crispri-collection]] (p. 9).
- Used as data by: [[molina-2026-response-decomposition]] (5 lines), [[huuloc-2026-cell-quota]] (noise regime),
  [[huang-2025-x-atlas-orion]] (depth comparator), and the project's [[rq2_params-assumptions]] / [[rq2_h2a-assumptions]].
- Published version: @zhang-2026-tahoe (not ingested; paywalled).
