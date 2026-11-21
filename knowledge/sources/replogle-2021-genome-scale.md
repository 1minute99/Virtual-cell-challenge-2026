---
title: "Replogle, Saunders et al. 2021/22 (preprint) — Genome-scale Perturb-seq"
type: source
created: 2026-11-06
updated: 2026-11-06
status: review
author: llm
bibkey: replogle-2021-genome-scale
based_on: preprint
depth: standard
tags: [perturb-seq, crispri, genome-scale, k562, rpe1, energy-distance, knockdown, transcriptional-phenotype, preprint]
relations:
  - target: replogle-nadig-crispri-collection
    type: supports
    confidence: extracted
    because: "Primary source of the K562 (genome-wide day 8, essential day 6) and RPE1 (essential, day 7) screens of the collection: «In total, we obtained >2.5 million high-quality cells with a median coverage of >100 cells per perturbation (fig. S1 E-G; table S4-S6).» (p. 6; design pp. 5–6)."
  - target: systematic-variation
    type: supports
    confidence: inferred
    because: "Reports response components shared by many perturbations (growth-gene downregulation after essential knockdowns, p. 12; phenotype strength tracking growth defect, p. 7; a homogeneous ISR-type nuclear response to mitochondrial knockdowns, p. 16; CIN-linked G1/G0 arrest in RPE1, p. 15) without naming them systematic variation."
---

# Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq

> Preprint version (replogle-2021-genome-scale, bioRxiv v3 posted 2022-05-28) of replogle-2022-genome-scale; claims quoted from the preprint.

> [!warning] Provenance: substitute edition (user-approved)
> This page is based on the **open-access bioRxiv preprint** (v3, 28 May 2022), not on the published Cell article
> (@replogle-2022-genome-scale, Cell 185(14):2559–2575.e28, 2022). The user approved preprint substitution on
> 2026-11-06 (option (b), `knowledge/_meta/log.md`). Every claim, number and quote below comes from the preprint PDF
> only. Page numbers refer to the preprint. Before citing a number as "Replogle et al. 2022", check it against the
> published version or cite the preprint key.

## Bibliographic Details

Replogle, Joseph M., Reuben A. Saunders, Angela N. Pogson, Jeffrey A. Hussmann, Alexander Lenail, Alina Guna, Lauren
Mascibroda, Eric J. Wagner, Karen Adelman, Gila Lithwick-Yanai, Nika Iremadze, Florian Oberstrass, Doron Lipson,
Jessica L. Bonnar, Marco Jost, Thomas M. Norman, and Jonathan S. Weissman. "Mapping information-rich
genotype-phenotype landscapes with genome-scale Perturb-seq." bioRxiv preprint doi:10.1101/2021.12.16.473013, v3
posted 28 May 2022 (v1 17 Dec 2021, v2 19 Dec 2021). CC BY-NC-ND 4.0. Co-first authors Replogle and Saunders;
corresponding authors Weissman and Norman (p. 1).
Read from `<library>/pdf/replogle-2021-genome-scale.pdf` (44 pp.). Citekey: @replogle-2021-genome-scale.

**Page numbers.** Printed page numbers equal PDF page numbers. The main text also carries running line numbers
(ll. 1–894). Figures 1–7 are images on pp. 32, 33, 35, 37, 39, 41 and 43, with legends on pp. 32–44.

**Not in this PDF.** The Materials and Methods, Figs. S1–S13 and Tables S1–S9 are in a separate supplement (p. 21,
p. 31) that was not read. The energy-distance test, the Anderson–Darling test, the z-normalisation, the sgRNA
design and the leverage scores are therefore described here only as far as the main text describes them. The data
statement on p. 31 gives no accession: raw data «will be deposited into SRA» and a data browser «will be made
available upon publication».

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract | 2 |
| ✓ F1 | Main Text (introduction; why genome-scale; open questions incl. cell-type rewiring) | 3–4 |
| ✓ F1 | A multiplexed CRISPRi strategy for genome-scale Perturb-seq (modality, library, three screens, knockdown; Fig. 1) | 5–6 |
| ✓ F1 | A robust computational framework to detect transcriptional phenotypes (z-normalisation, energy distance, AD test, growth, lower bound, neighbour knockdown) | 6–8 |
| ✓ | Annotating gene function from transcriptional phenotypes (1,973 strong perturbations, CORUM, STRING, clusters, ribosome validation; Fig. 2) | 8–10 |
| ✓ | Delineating functional modules of the Integrator complex (C7orf26 / INTS15; Fig. 3) | 10–11 |
| ✓ F1 | Data-driven definition of transcriptional programs (expression-program map, ISR vs UPR, differentiation; Fig. 4A–E) | 12–13 |
| ✓ F1 | Hypothesis-driven study of composite phenotypes (TE fraction, total RNA; Fig. 4F–H) | 13–14 |
| ✓ F1 | Exploring genetic drivers and consequences of aneuploidy in single cells (leverage scores, inferCNV, CIN vs cell cycle, CIN scores; Fig. 5) | 14–15 |
| ✓ F1 | Discovery of stress-specific regulation of the mitochondrial genome (Figs. 6–7) | 16–18 |
| ✓ F1 | Discussion (interpretability, composite phenotypes, CIN, mitochondria, cost, Ultima sequencing) | 18–20 |
| — | Materials and Methods pointer; References | 21–29 |
| ✓ F1 | Acknowledgements, author contributions, funding, competing interests, data availability | 30–31 |
| ✓ F1 | Figure 1–7 legends and figure images | 32–44 |
| — | Supplement (Methods, Figs. S1–S13, Tables S1–S9) | not in PDF |

## Kernthesen

### 1. ⭐ Genome-scale Perturb-seq is feasible, and scale has an argument of its own (pp. 2, 4)

The paper presents what it calls the first genome-scale Perturb-seq screens: CRISPRi against all expressed genes with
single-cell RNA-seq readout across more than 2.5 million human cells (pp. 2, 4). The case for scale is that cost grows
linearly in the number of perturbations while the number of pairwise comparisons — and hence the value for
unsupervised classification of gene function — grows quadratically, and that well-characterised genes serve as
built-in anchors for interpretation (p. 4). Two questions named as only answerable at this scale are what fraction of
genetic changes elicit global transcriptional phenotypes and how transcriptional programs are rewired between cell
types (p. 4).

### 2. CRISPRi was chosen over knockout for measurability and homogeneity (p. 5)

Three reasons are given: knockdown efficacy can be read directly from the scRNA-seq data, so each gene needs only one
library element and unperturbed genes can be excluded empirically; CRISPRi is more homogeneous than nuclease knockout
(which leaves cells with active in-frame indels), limiting selection for unperturbed cells among essential genes; and
CRISPRi does not activate the DNA-damage response (p. 5).

### 3. Library design: two sgRNAs per construct, with toxic constructs over-represented (p. 5)

Each construct carries two distinct sgRNAs against the same gene (tables S1–S3, not in this PDF). Preliminary growth
screens were run first, and constructs causing strong growth defects were over-represented during library synthesis so
that essential-gene guides would not be lost to dropout (p. 5; fig. S1A–D, not in this PDF).

### 4. ⭐⭐ Three screens, two cell lines, three time points (pp. 5–6; Fig. 1A, p. 32)

- **K562 genome-wide**, all expressed genes, day 8 after transduction: n = 9,866 genes, n = 10,673 perturbations (some
  genes have several independent transcripts), in K562 expressing dCas9-KRAB (p. 5).
- **K562 essential**, common essential genes, day 6: n = 2,057 genes, n = 2,176 perturbations (pp. 5–6).
- **RPE1 essential**, common essential genes plus a subset of non-essential genes that gave phenotypes in K562, day 7:
  n = 2,393 genes, n = 2,549 perturbations (p. 6).

RPE1 expresses dCas9 fused to a ZIM3-derived KRAB domain, reported as a stronger repressor than the KOX1-derived
domain used previously, and is described as non-cancerous, hTERT-immortalised, near-euploid, adherent and p53-positive,
in contrast to K562 (pp. 5–6). Note that the Fig. 1A schematic (p. 32) gives slightly different gene counts for the
same screens — n = 9,867 genome-wide and n = 2,285 essential-wide — than the running text.

### 5. ⭐⭐ Scale of the data: >2.5 M cells, median >100 cells per perturbation, median ~10,000 UMIs per cell (p. 6; Fig. 1A, p. 32)

After sgRNA identification and removal of cells carrying guides against different genes (a product of lentiviral
recombination or droplet doublets), the three screens together yield «>2.5 million high-quality cells with a median
coverage of >100 cells per perturbation» (p. 6). The Fig. 1A panel adds the sequencing depth: «median ~10,000 UMIs per
cell» (p. 32). The screens use 10x Genomics droplet-based 3′ scRNA-seq with direct sgRNA capture (p. 6). Per-screen
cell counts are in tables S4–S6, which are not in this PDF.

### 6. ⭐⭐ Median on-target knockdown is 85.5 % in K562 and 91.6 % in RPE1 (p. 6; Fig. 1B, p. 32)

Knockdown is measured from the scRNA-seq data itself. The difference between the lines is attributed to the
ZIM3-derived KRAB domain used in RPE1, and the authors take it as a sign that «future efforts would benefit from
improved CRISPRi efficacy» (p. 6). The Fig. 1B cumulative-density plot is over n = 9,464 target genes in K562 and
n = 2,333 in RPE1 (p. 32) — fewer genes than were targeted, so some targets are not represented in the knockdown
statistic.

### 7. ⭐ Non-targeting controls allow internal z-normalisation, which also corrects batch effects (pp. 6–7)

The design includes many control cells carrying diverse non-targeting sgRNAs. Expression is z-normalised internally
against these, and the authors report that this corrected batch effects arising from parallelised scRNA-seq and
sequencing (p. 7; fig. S2, not in this PDF). Because each cell carries its own perturbation label, statistical tests
treat cells as independent samples, and the tests chosen are deliberately conservative and non-parametric (pp. 6–7).

### 8. ⭐⭐ By a permuted energy-distance test, 31.1 % of K562 perturbations have a global transcriptional phenotype, against a 1.9 % control rate (p. 7)

Cells carrying each perturbation are compared with non-targeting control cells at the level of principal components,
and the energy distance is tested against a permuted null — the question being whether the perturbed cells could have
been drawn from the control population. The result: «2,987 of 9,608 genetic perturbations targeting a primary
transcript (31.1%) compared to 11 of 585 non-targeting controls (1.9%) caused a significant transcriptional phenotype
in K562 cells» (p. 7). The non-targeting rate is the paper's own noise floor for this test. The corresponding figures
for the RPE1 and K562-essential screens are in tables S4–S6, which are not in this PDF.

### 9. ⭐⭐ An Anderson–Darling DE test gives 30.5 %, concordant with the energy distance at 78.7 % Jaccard (p. 7)

The energy-distance test measures a global shift without naming genes, so a per-gene Anderson–Darling test is run in
addition; it is sensitive to changes in a subset of cells, i.e. to incomplete penetrance. «2,935 of 9,608 genetic
perturbations targeting a primary transcript (30.5%) compared to 12 of 585 non-targeting controls (2.1%) caused >10
differentially expressed genes in K562 cells» (p. 7). The authors report these results as «well-correlated between
time points and cell types» (fig. S3A,B and tables S4–S6, not in this PDF) and «concordant with the energy distance
test (78.7% concordance by Jaccard index)» (p. 7). The two tests therefore give the same order of magnitude by two
different definitions of "phenotype".

### 10. ⭐ Transcriptional phenotype tracks growth defect, but the two are not the same screen (p. 7; Fig. 1C, p. 32)

Phenotype strength correlates with growth defect (Spearman's rho = −0.51), and «86.6% of essential genetic
perturbations (gamma < –0.1) [lead] to a significant transcriptional response in K562 cells» (p. 7). Yet n = 771
perturbations produce a transcriptional phenotype with a negligible growth phenotype, «indicating that many genetic
perturbations influence cell state but not growth or survival» (p. 7). Highly expressed genes are also more likely to
give a phenotype (rho = 0.42) (p. 7).

### 11. ⭐⭐ The detectable fraction is explicitly a lower bound (pp. 7–8)

Because some perturbations did not achieve strong on-target knockdown, «our estimate of the fraction of genetic
perturbations that cause a transcriptional phenotype is likely to be a lower bound» (p. 7). On the other side,
about 7.5 % of perturbations significantly knock down a neighbouring gene in K562, but such neighbour knockdown is not
enriched among the perturbations with a phenotype and no growth defect (p. 8; fig. S4, not in this PDF). The summary
offered is «a coherent picture where knockdown of a significant fraction of expressed genes causes a transcriptional
response» (p. 8).

### 12. 1,973 perturbations are carried into the function analysis as the "strong phenotype" subset (p. 8; Figs. 2A, 4A)

For gene-function inference the authors restrict to 1,973 perturbations with more than 50 differentially expressed
genes by the AD test, compared and clustered on the correlation of mean expression profiles over 2,319 highly
variable genes (p. 8; Fig. 2A legend, p. 33; Fig. 4A legend, p. 37). Correlation is used because related perturbations
may differ in effect magnitude, so a scale-invariant similarity is wanted (p. 8). This 1,973-perturbation subset is
the one most downstream benchmarks of this dataset inherit.

### 13. ⭐ Expression-profile correlation recovers curated biology (pp. 8–9; Fig. 2B–D)

Among the 1,973 genes there are 327 CORUM 3.0 complexes with at least two thirds of their members present
(14,165 confirmed protein–protein interactions). Their expression-profile correlations have median 0.61 against a
background median of 0.10 over all gene pairs (p. 9). High STRING scores likewise go with high expression correlation
over the 243,558 scored gene pairs (p. 9). An unbiased clustering yields 64 perturbation clusters, visualised in a
minimum-distortion embedding (pp. 9–10; Figs. 2D, 4A).

### 14. Correlation-based function prediction validates experimentally in 9 of 10 cases (pp. 9–10; Fig. 2E)

Ten poorly annotated genes whose perturbation profile correlated (r > 0.6) with ribosomal subunits or biogenesis
factors were retargeted by CRISPRi and assayed by 28S/18S rRNA ratio. Nine showed ribosome-biogenesis defects (all but
ABCF1), and in each case the affected subunit matched the Perturb-seq cluster across two independent sgRNAs (p. 10).
The authors caution that some such relationships «might be explained by indirect or off-target effects» (p. 10).

### 15. C7orf26 is a fourth subunit of an Integrator module, proposed as INTS15 (pp. 10–11; Fig. 3)

Integrator subunits split into three transcriptional modules (INTS1-2-5-7-8; INTS3-4-9-11; INTS10-13-14 with the
uncharacterised C7orf26), mirroring the structural architecture. C7orf26 depletion destabilises INTS10; His-INTS10
pulldown recovers C7orf26 with INTS13/14; and the four proteins co-purify by size-exclusion chromatography (pp. 10–11).
Only the cleavage and shoulder/backbone modules cause splicing defects and snRNA readthrough by PRO-seq, not the
INTS10-13-14-C7orf26 module (p. 11). The renaming to INTS15 is proposed on p. 11.

### 16. ⭐ A program × perturbation-cluster map recovers known regulatory logic, including a growth program (p. 12; Fig. 4B–C, p. 37)

Genes are clustered into expression programs by co-regulation across perturbations, perturbations are clustered by
profile, and the mean activity of each program in each perturbation cluster is tabulated (p. 12; 37 gene-expression
clusters × 64 perturbation clusters, Fig. 4A–B, p. 37). Recovered relationships include proteasome upregulation after
proteasome dysfunction, NFkB activation on ESCRT loss, cholesterol-biosynthesis upregulation on trafficking defects,
and «downregulation of growth-related genes in response to essential genetic perturbations» (p. 12). The map also
separates the two canonical stress branches, UPR and ISR, as independently regulated axes: ISR by loss of
mitochondrial proteins, aminoacyl-tRNA synthetases and translation-initiation factors, UPR by loss of ER chaperones
and translocation machinery (p. 12; Fig. 4C).

### 17. Many perturbations drive erythroid or myeloid differentiation in K562, including common essential genes (pp. 12–13; Fig. 4D–E)

Loss of erythropoiesis regulators (GATA1, LDB1, LMO2, KDM1A) drives myeloid differentiation; knockdown of BCR-ABL and
GAB2 drives erythroid differentiation. Unexpectedly, loss of common essential genes also pushes cells towards myeloid
(e.g. Integrator subcomplex) or erythroid (NuRD, DNA replication) markers (p. 12). PTPN1, selectively essential in
K562, drives myeloid differentiation, and PTPN1 + KDM1A knockdown together gives more differentiation and a larger
growth defect than either alone (p. 13; Fig. 4E).

### 18. Composite phenotypes are reproducible across replicates and cell types (pp. 13–14; Fig. 4F–H)

Phenotypes computed across the transcriptome — transposable-element RNA fraction, mitochondrial RNA fraction, total
RNA content (UMIs) — are «under strong genetic control, with highly reproducible effects across screen replicates and
cell types» (p. 13; Fig. 4F). Between the two K562 time points, TE fraction correlates r = 0.79 and total RNA content
r = 0.74 (Fig. 4G–H legend, p. 38). Exosome loss and CPSF/Integrator loss both raise TE RNA; S-phase and mitosis
perturbations raise total RNA content, while loss of general transcription machinery lowers it (pp. 13–14).

### 19. ⭐ Single-cell leverage scores quantify phenotype magnitude and heterogeneity (p. 14; Fig. 5A–B)

Leverage scores from an SVD measure how outlying each perturbed cell is relative to non-targeting controls, without
assuming a single axis of variation. Mean leverage per perturbation correlates with DEG count (rho = 0.71) and
reproduces across the day-6 and day-8 K562 screens (rho = 0.79) (p. 14; fig. S10A–B, not in this PDF). The standard
deviation of leverage scores is then used as a heterogeneity statistic, and chromosome-segregation genes (TTK, SPC25,
DSN1) are among the most heterogeneous (pp. 14–15).

### 20. ⭐ Chromosomal instability can be read per cell, and its consequences differ between RPE1 and K562 (pp. 14–15; Fig. 5C–G)

inferCNV on the expression data detects karyotype changes in single cells: on TTK knockdown, 61/80 (76 %) of RPE1
cells show karyotypic changes against 274/13,140 (2 %) of unperturbed cells (pp. 14–15). Cell-cycle consequences
differ by line: RPE1 cells with abnormal karyotypes arrest in G1/G0 (fraction 0.68 vs 0.44 for stable karyotypes),
while «K562 cells with altered karyotypes had less significant shifts in cell cycle occupancy» (p. 15) — the
difference being framed by euploid, p53-positive RPE1 versus constitutively aneuploid, p53-deficient K562 (p. 15).
Cells with both gains and losses upregulate the ISR (p. 15). Across perturbations, the CIN score correlates
r = 0.69 between the two cell lines (p. 15; Fig. 5G).

### 21. ⭐⭐ The nuclear response to mitochondrial perturbations is largely generic; the mitochondrial transcriptome is specific (pp. 16–17; Figs. 6–7)

For 268 mitochondrial perturbations with strong phenotypes, nuclear transcriptional responses are relatively
homogeneous and well-correlated, and «nuclear transcriptional responses generally failed to discriminate genetic
perturbations by function» (p. 16), beyond ISR magnitude. The 13 mitochondrially encoded genes, by contrast, separate
Complex I, III, IV, V, the two mito-ribosome subunits, chaperone/import machinery and RNA-processing factors; a random
forest classifier reaches accuracy 0.64 on the mitochondrial transcriptome against 0.25 on the nuclear one (p. 17;
fig. S12F, not in this PDF). Bulk RNA-seq without poly-A selection reproduces the Perturb-seq pattern (cophenetic
r = 0.79, p. 17). The model proposed is «a general nuclear response … layered over perturbation-specific changes in
the expression level of mitochondrially encoded genes» (p. 19; Fig. 7F). TMEM242 is assigned to ATP synthase on this
signature in both K562 and RPE1 and validated by Seahorse respirometry (pp. 17–18).

### 22. ⭐ The stated limitation of single-cell CRISPR screens is cost, not biology (pp. 19–20)

«At present, the major limitation of single-cell CRISPR screens is cost» (p. 20). The remedies named are multiplexed
libraries, compressed sensing, advances in single-cell phenotyping, and cheaper sequencing; the authors sequenced the
genome-scale K562 libraries on the Ultima Genomics platform and report results «equivalent to those sequenced on
Illumina instruments» (p. 20; fig. S13, not in this PDF). They also argue these screens need far fewer cells than
other approaches, making them suited to iPSC-derived and in vivo material (pp. 19–20).

## Figures and tables

- **Fig. 1 (p. 32):** (A) screen schematic with the three screens and the scale annotations (>2.5 M cells after
  filtering, median >100 cells per perturbation, median ~10,000 UMIs per cell); (B) cumulative on-target knockdown
  for 9,464 K562 and 2,333 RPE1 target genes; (C) growth phenotype (gamma) vs log10 DEG count per guide pair in K562,
  points coloured by energy-distance significance (rho = −0.51).
- **Fig. 2 (p. 33):** 1,973 strong perturbations over 2,319 HVGs; CORUM and STRING correlation distributions;
  minimum-distortion embedding; 28S/18S validation.
- **Fig. 3 (pp. 35–36):** Integrator modules in K562 and RPE1; biochemistry; splicing and PRO-seq; structure.
- **Fig. 4 (p. 37):** 37 gene-expression clusters × 64 perturbation clusters heatmap; ISR vs UPR; erythroid vs myeloid;
  composite-phenotype correlations across time points and cell types (TE r = 0.79, total RNA r = 0.74).
- **Fig. 5 (pp. 39–40):** leverage-score heterogeneity; inferCNV heatmap; cell-cycle occupancy by karyotype status
  (RPE1 240,768 stable / 5,522 loss / 1,987 gain / 904 both); CIN score correlation r = 0.69 across lines.
- **Figs. 6–7 (pp. 41–44):** mitochondrial perturbations (n = 268) by nuclear vs mitochondrial transcriptome;
  mitochondrial genome read positions; Perturb-seq vs bulk RNA-seq (r = 0.79); TMEM242.
- **Tables S1–S9 and Figs. S1–S13 are not in this PDF** — including the per-screen cell counts and per-perturbation
  statistics (tables S4–S6) that later re-users of this dataset cite.

## What this source does not address

- **No cross-cell-line prediction and no variance decomposition.** The two lines are compared qualitatively
  (Integrator modules, CIN scores r = 0.69, TMEM242) but there is no model that predicts one line's responses from the
  other, no train/test split over perturbations, and no decomposition of responses into shared, perturbation-specific
  and line-specific parts. The decomposition applied to these data is @molina-2026-response-decomposition's, not this
  paper's.
- **No ceiling, no noise model.** There is no split-half reliability estimate, no noise-corrected correlation, and no
  statement of an achievable prediction accuracy. The only "floor" is the non-targeting significance rate of the two
  tests (1.9 % / 2.1 %, p. 7).
- **No context-count analysis.** Two cell lines, three screens. The paper names cross-cell-type rewiring as a
  motivation (p. 4) but does not quantify it.
- **No control-subspace or residual geometry.** Principal components are used inside the energy-distance test and the
  leverage scores (pp. 7, 14), never as a basis onto which responses or prediction residuals are projected.
- **No machine-learning benchmark.** No perturbation-prediction model of any kind is trained or evaluated.
- **Methods, supplementary figures and all tables are outside this PDF** (p. 21, p. 31), so the energy-distance
  implementation, the exact per-screen cell counts, the sgRNA sequences and the knockdown distributions cannot be read
  here. The preprint also gives **no data accession**: «Raw sequencing data will be deposited into SRA» (p. 31).

<!-- ═══ FOKUS ═══ -->

## Focus: facts and claims the manuscript cites; exactly what overlaps with H1–H4 — 2026-11-06

**Carries Kernthesen 4, 5, 6, 8, 9, 10, 11, 12, 16, 20, 21.** The hypotheses are those of
`input/ideas/zero-shot-limits-context-complexity-plan.md`. This is the primary screen paper behind the K562 and RPE1
halves of [[replogle-nadig-crispri-collection]]; the project cites it for screen-level facts, not for any modelling
result.

### What this focus needs beyond the Kernthesen

1. **Screen design, as the manuscript must state it.** K562 genome-wide (all expressed genes, day 8, 9,866 genes /
   10,673 perturbations), K562 essential (day 6, 2,057 / 2,176), RPE1 essential plus K562-phenotype non-essentials
   (day 7, 2,393 / 2,549); CRISPRi with two sgRNAs per construct; dCas9-KRAB in K562, ZIM3-KRAB in RPE1; 10x 3′ with
   direct guide capture (pp. 5–6). Caution: Fig. 1A (p. 32) prints n = 9,867 and n = 2,285 for the same screens, so
   the counts in the running text and in the figure disagree slightly — cite the text and note the figure.
2. **Depth and coverage.** >2.5 M cells after filtering; median >100 cells per perturbation (p. 6); median ~10,000
   UMIs per cell (Fig. 1A, p. 32). The project's own comparison numbers — e.g. the 11,499 median UMIs that
   [[huang-2025-x-atlas-orion]] quote for Replogle K562 — are **not** in this preprint; only the ~10,000 figure is.
3. **Knockdown.** Median on-target 85.5 % (K562) and 91.6 % (RPE1), measured from the scRNA-seq itself (p. 6), over
   9,464 and 2,333 target genes (Fig. 1B, p. 32). This confirms the 85 % figure for K562 GWPS carried in the project's
   `AGENTS.md` and noted as unsourced in [[huang-2025-x-atlas-orion]] (p. 288 of that page's numbering).
4. **Detectable-phenotype fraction and its noise floor.** Energy distance: 2,987 / 9,608 = 31.1 % significant vs
   11 / 585 = 1.9 % of non-targeting controls. AD test: 2,935 / 9,608 = 30.5 % with >10 DEGs vs 12 / 585 = 2.1 %.
   Concordance 78.7 % by Jaccard (p. 7). All of these are K562 genome-wide, day 8; the equivalent RPE1 numbers are in
   tables S4–S6 and are not in this PDF. The 31 % is declared a **lower bound** because of incomplete knockdown
   (p. 7), and 86.6 % of essential perturbations (gamma < −0.1) are significant (p. 7) — so the fraction is strongly
   conditional on the gene panel, which is why essential-enriched screens look so different from genome-wide ones.
5. **K562 vs RPE1.** The comparison here is qualitative and is about *mechanistic* rewiring, not about prediction
   transfer: the CIN score correlates r = 0.69 across lines (p. 15; Fig. 5G); the Integrator module structure is
   reproduced in both (Fig. 3B); TMEM242's ATP-synthase signature appears in both (p. 17). The sharpest documented
   difference is the aneuploidy response — RPE1 (near-euploid, p53-positive) arrests in G1/G0 with abnormal karyotypes
   (0.68 vs 0.44) while K562 (aneuploid, p53-deficient) does not shift much (p. 15). That is the mechanism behind the
   cell-cycle difference [[vinastorne-2026-systema]] measure as [[systematic-variation]] in the same two lines.
6. **Shared / generic response component.** The paper has no decomposition, but it names three generic components:
   downregulation of growth-related genes after essential knockdowns (p. 12), phenotype strength tracking growth
   defect (rho = −0.51, p. 7), and a mitochondrial-stress case where «a general nuclear response is layered over
   perturbation-specific changes in the expression level of mitochondrially encoded genes» (p. 19). The layered
   framing is the closest this source comes to the project's δ = μ + a + b + i notation, and it is only stated for
   mitochondrial perturbations. Ribosome biology appears as a *cluster* (ribosome biogenesis, 40S; p. 10) and as the
   validation assay, not as a generic response axis.

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **not addressed (data only)** | No ceiling, no noise correction, no prediction model. It supplies the measurement substrate for the RN_K562 panel and its depth/knockdown characteristics, which is what r_max is estimated on. The only "floor" stated is the non-targeting significance rate (1.9 % / 2.1 %) of two phenotype tests, which is not a noise ceiling on correlation. | pp. 6–7 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **not addressed** | Genetic CRISPRi only, two cell lines, no chemical screen, no context curve. | pp. 5–6 |
| **H2b** (theory predicts the break-even within ×2) | **not addressed** | No theory of context sample complexity. | — |
| **H2c** (genetic C*_min > C_public) | **partial (registry input)** | Contributes **two** public contexts, K562 and RPE1, both genome- or essential-scale CRISPRi with median >100 cells per perturbation and strong knockdown — i.e. exactly the usability criteria C_public applies. It gives no count of available cell lines and no C*_min. | pp. 5–6 |
| **H3a** (transfer residual loads on the destination's top-8 control PCs) | **not addressed (indirect support for the premise)** | PCs appear only inside the energy-distance test and the leverage scores, never as a projection basis, and there is no transfer residual. Indirectly it documents why a destination's control subspace should be informative: the two lines differ in ploidy, p53 status and cell-cycle response to CIN, i.e. in control-state structure, not only in perturbation response. | pp. 7, 14, 15 |
| **H3b** (M3 beats five alternative corrections on PDS) | **not addressed** | No correction method, no PDS, no deflation. Internal z-normalisation against non-targeting controls (p. 7) is a normalisation of the data, not a correction of a prediction. | p. 7 |
| **H4** (same-line genome-wide data beats the largest modelling gain) | **partial (the data object itself)** | H1GW's counterpart on the RN panels is this K562 genome-wide screen; the source establishes what "same-line genome-wide data" is and what it costs (>2.5 M cells; «the major limitation … is cost», p. 20). It runs no model and therefore makes no data-vs-model comparison. | pp. 5–6, 20 |

### Direct quotes (for the draft)

> "we found that 2,987 of 9,608 genetic perturbations targeting a primary transcript (31.1%) compared to 11 of 585
> non-targeting controls (1.9%) caused a significant transcriptional phenotype in K562 cells" (p. 7)

> "Considering that some of our genetic perturbations did not yield strong on-target knockdown, our estimate of the
> fraction of genetic perturbations that cause a transcriptional phenotype is likely to be a lower bound." (p. 7)

> "We observed a median target knockdown of 85.5% in K562 cells and 91.6% in RPE1 cells" (p. 6)

> "In total, we obtained >2.5 million high-quality cells with a median coverage of >100 cells per perturbation" (p. 6)

> "A substantial number of genetic perturbations that cause a transcriptional phenotype nonetheless have a negligible
> growth phenotype (n=771 …), indicating that many genetic perturbations influence cell state but not growth or
> survival." (p. 7)

> "a general nuclear response is layered over perturbation-specific changes in the expression level of
> mitochondrially encoded genes" (p. 19)

### Examples & illustrations (for later drafting)

- **The 31 % figure is panel-dependent.** 31.1 % of all expressed genes give a phenotype, but 86.6 % of growth-defect
  perturbations do (p. 7). Any project statement about "how many perturbations do something" must say which panel.
- **Knockdown is a design variable, not a constant.** 85.5 % (KOX1-KRAB, K562) vs 91.6 % (ZIM3-KRAB, RPE1) in the same
  study (p. 6) — the line comparison is confounded with the effector.
- **Line difference with a named mechanism.** Abnormal-karyotype cells: G1/G0 fraction 0.68 vs 0.44 in RPE1, little
  shift in K562, attributed to p53 status (p. 15). Useful whenever the manuscript explains why RPE1 carries more
  systematic variation than K562.

### Boundary of this focus

This focus takes only screen-level facts and the two cell lines' qualitative comparison. It does not cover the
Integrator / INTS15 biochemistry (Kernthese 15), the differentiation and composite-phenotype results (17–18), or the
mitochondrial-genome mechanism (21) beyond the single layered-response statement. Nothing here can be cited for
ceilings, transfer or context counts.

## Preprint vs published version

**Version history (bioRxiv API, 2026-11-06):** v1 2021-12-17, v2 2021-12-19, v3 2022-05-28 — all CC BY-NC-ND, all
typed "new results", all already linked to the published DOI 10.1016/j.cell.2022.05.013. The PDF on disk is v3,
posted after Cell acceptance (p. 1 header). Published version: Cell 185(14):2559–2575.e28, July 2022
(@replogle-2022-genome-scale; Crossref).

**Differences: largely not audited.** Only the abstracts were compared (Crossref has none; the published abstract was
read from OpenAlex purely for this comparison and is not quoted anywhere above). They are near-identical; the
published abstract drops «pooled» before "CRISPR-based screens", drops the preprint's clause «and present a framework
to power biological discovery with the resulting genotype-phenotype map», writes "poorly characterized" without the
hyphen, and ends «gene and cellular function» instead of «gene function and cellular behavior». No numerical claim
differs between the two abstracts.

Everything else — the three screens' gene and perturbation counts, >2.5 M cells, the knockdown medians, the 31.1 % /
30.5 % phenotype fractions and their control rates, the 1,973-perturbation subset, figure numbering, and all
supplementary tables — is **not audited** against the published article. The published version has 28 supplementary
items (".e28") against the preprint's Figs. S1–S13 / Tables S1–S9 (p. 31), so the supplement was clearly reorganised.
The preprint's data statement is also obsolete: it promises SRA deposition and a browser «upon publication» (p. 31),
whereas the released dataset (figshare / GEO, as cited by later papers in this wiki) postdates it. **Cite screen
numbers from this page as preprint numbers**, or verify them against the Cell version before the manuscript goes out.

## Mentioned entities

- Datasets: [[replogle-nadig-crispri-collection]] (the K562 genome-wide, K562 essential and RPE1 essential screens
  introduced here)
- Concepts: [[systematic-variation]] (not named here; the p53 / cell-cycle mechanism on pp. 15 is what
  [[vinastorne-2026-systema]] later measure), [[perturbation-response-decomposition]] (applied to these data by
  [[molina-2026-response-decomposition]], not by this paper)

## Connections

- **Introduces the data for:** [[replogle-nadig-crispri-collection]] — K562 and RPE1 halves (pp. 5–6).
- **Underlies (as data, without stating the result):** [[perturbation-response-decomposition]] and the variance
  shares reported by [[molina-2026-response-decomposition]]; the cell-cycle / CIN contrast measured as
  [[systematic-variation]] by [[vinastorne-2026-systema]]; the depth and batch-reproducibility comparisons in
  [[huang-2025-x-atlas-orion]]; the unseen-perturbation benchmarks of [[ahlmanneltze-2025-linear-baselines]] and
  [[piao-2026-adaptive-contexts]]; the held-out contexts of [[pan-2026-global-cell-state]] and
  [[adduri-2025-state]].
- **Published counterpart:** @replogle-2022-genome-scale (Cell 185, 2022) — same work, not ingested.
- **Cites (not ingested):** Adamson et al. 2016 (ref. 9), Dixit et al. 2016 (ref. 10), Alerasool et al. 2020
  (ZIM3-KRAB, ref. 24), Rosenbluh et al. 2017 (neighbour knockdown, ref. 25), CORUM (ref. 26), STRING (ref. 27),
  inferCNV / Patel et al. (ref. 43).
