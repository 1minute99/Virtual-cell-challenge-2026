---
title: "H1 hESC genome-scale CRISPRi Perturb-seq (GSE295214; the project's H1GW)"
type: entity
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [dataset, crispri, perturb-seq, hesc, h1, genome-scale, gse295214, h1gw]
sources: ["[[pan-2026-global-cell-state]]"]
relations:
  - target: pan-2026-global-cell-state
    type: cites
    confidence: extracted
    because: "The dataset is generated and released by this paper: «we performed a genome-scale CRISPRi Perturb-seq experiment in human embryonic stem cells (hESCs)» (p. 14); «accession number: GSE295214» (p. 20)."
---

# H1 hESC genome-scale CRISPRi Perturb-seq (GSE295214)

## Location and authority data

NCBI GEO accession **GSE295214** (@pan-2026-global-cell-state, p. 20). Generated in the Weissman lab (Saunders and
Replogle designed and performed the screen; author contributions, p. 20) and published with the SCMG framework paper
of Pan, Saunders, Replogle, Weissman and Zhuang (bioRxiv, 18 May 2026, doi:10.64898/2026.05.16.725005). Processed
pseudo-bulk perturbation data: HuggingFace `xingjiepan/SCMG_data`; code `github.com/xingjiepan/SCMG` (p. 20).

In this project's registry (`output/data-analysis/results/families.json`) the dataset is the family **`H1GW`** —
"H1 hESC · GSE295214 (Weissman lab, 10x 5′) · genome-scale (2,747 targets passing QC)" — distinct from the family
`H1`, which is the VCC 2025 H1 training data (300 targets, Arc Institute, 10x Flex). Both map to the one cell line
"H1 hESC" in `line_to_families`, so the pair counts **once** towards C_public.

## What one sees

A CRISPR-interference Perturb-seq screen in **H1 human embryonic stem cells** carrying inducible **Zim3-dCas9**,
cultured in StemFlex and passaged with Accutase under ROCKi (@pan-2026-global-cell-state, pp. 59–60).

- **Library:** 6,638 sgRNAs targeting **2,978 genes** (two sgRNAs each) plus **370 non-targeting controls**, drawn
  from CRISPRiV2 / CRISPRiV2.1 (pp. 14, 60).
- **Target selection:** not unbiased genome-wide. A prior genome-wide CRISPRi **balancing screen** identified genes
  with a growth phenotype or altered OCT4 / NANOG / SOX2 expression (day-11 FACS on triple-positive vs
  single-negative cells); genes with ≥ 2 sgRNAs showing evidence were kept, and dropout-prone sgRNAs were
  over-represented in the pool in proportion to their expected dropout (pp. 59–60).
- **Protocol:** transduction > 1000× coverage at < 20 % rate; puromycin selection; dCas9 induced with 25 nM
  doxycycline; harvest **8 days** after induction; GFP sort; **10x 5′ scRNA-seq with a spike-in primer against the
  sgRNA constant region** for direct capture (p. 60).
- **Scale:** **≈ 1.3 M** perturbed cells with assigned perturbations; **tens to hundreds of cells per perturbation**
  (Fig. 6B, p. 40; mode ≈ 200); median on-target knockdown **≈ 70 %** residual-expression reduction (Fig. 6C).
- **QC:** cells with ≥ 500 detected genes and exactly one sgRNA at ≥ 20 barcode UMIs (cellranger 8.0.1);
  perturbations with ≥ 50 cells, target expression < 90 % of the non-targeting mean, one-sided Mann–Whitney
  p < 0.01; a guide-consistency filter (shift norm > 4 with between-guide r < 0.2) **removed 259 genes** (pp. 60–61).
- **Not reported:** per-cell UMI or gene depth, and any cross-batch reproducibility statistic.

## History of research

Released in May 2026 as the validation context for the SCMG framework: the authors chose hESCs because they are
«developmentally and transcriptionally distinct from K562 and RPE1 cells — the two lines with previously published
genome-scale CRISPRi Perturb-seq datasets» (@pan-2026-global-cell-state, p. 14), and to separate "is the difference
the cell state or the perturbation type?" from an earlier hESC **ORF overexpression** screen (Nakatake et al.,
ref. 63 there).

## The evidence

- **Attenuated responses.** Perturbation-induced shift magnitudes are significantly smaller in hESCs than in K562,
  on both matched-knockdown-efficiency subsets (p < 10⁻⁴; Fig. 6D, p. 40) — read as genuine stress resistance of the
  hESC state, not a knockdown artefact (p. 14).
- **Conserved structure.** Essential-gene knockdown moves essential programs and largely spares differentiation
  programs, as in K562 and RPE1; lysosome/autophagy and p53 respond broadly, ISR and UPR selectively (p. 14).
- **Its own rewiring.** Mitochondrial transcription/translation perturbations produce the *largest* stress responses
  here, unlike in RPE1 (p. 14); stress-induced mesenchymal transitions (56 enriched perturbed genes) and germ-layer
  exits after POU5F1 / SOX2 / NANOG knockdown appear at single-cell resolution (pp. 16–17).
- **Hardest of the three contexts to predict:** few-shot accuracy medians ≈ 0.09 (K = 10) to ≈ 0.18 (K = 200),
  against ≈ 0.44 in RPE1 (Fig. 6G vs Figs. 4D–E).

## Contested

- Whether "genome-scale" is the right label: the profiled set is the phenotypic subset (2,978 genes) of a
  genome-wide balancing screen, not the genome. The project's registry counts it as genome-scale
  (2,747 targets passing its own QC); an unselected genome-wide screen would have a different effect-size
  distribution.
- Whether the lower effect sizes reflect biology (stress resistance) or an hESC-specific CRISPRi/assay effect. The
  source rules out knockdown efficiency alone (p. 14) but no second hESC screen exists to confirm it.

## Relevance to this project

This is the **H1GW** source of hypothesis **H4** (does adding same-line genome-scale data at w = 1 gain more composite
than the largest single modelling gain?), and it is the counted H1-hESC row of **C_public** for H2c. Because it is a
second screen of the *same* cell line as the H1 destination panel, it is "more data in the destination context",
not cross-context transfer — the distinction the H4 write-up has to make explicitly.

## Multiple roles

Source dataset in the M1 registry (weight 0 for somatic destinations, same-line fallback for stem destinations);
the validation context of @pan-2026-global-cell-state's conservation claims; and, through its attenuated effect
sizes, an argument for a per-destination noise ceiling.

## Open / to verify

- Whether the project's 2,747 QC-passing targets correspond to the source's post-filter set (2,978 targeted minus
  259 guide-inconsistent = 2,719 — close but not equal; the difference is unexplained and should be checked against
  the GEO processed files).
- Per-cell sequencing depth, absent from the paper.

## Sources

- [[pan-2026-global-cell-state]] — pp. 14, 20, 40 (Fig. 6), 59–61.
