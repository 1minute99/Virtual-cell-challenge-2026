---
title: "X-Atlas/Orion (genome-wide CRISPRi Perturb-seq, HCT116 and HEK293T)"
type: entity
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
tags: [dataset, crispri, perturb-seq, genome-wide, hct116, hek293t]
sources: ["[[huang-2025-x-atlas-orion]]", "[[huuloc-2026-cell-quota]]"]
---

# X-Atlas/Orion

## Location and authority data

This is a dataset, so no ORCID, GND or Wikidata identifier applies. Primary publication: [[huang-2025-x-atlas-orion]]
(@huang-2025-x-atlas-orion; bioRxiv, June 2025; Xaira Therapeutics). Data: Figshare,
doi:10.25452/figshare.plus.29190726. The h5ad files there contain only cells with two sgRNAs from the same guide pair.
The release of all aligned cells was announced «at a later date» (p. 23).

## What one sees

Two genome-wide CRISPRi Perturb-seq screens (dCas9-KRAB-ZIM3; one tRNA-separated dual-sgRNA pair per gene; 18,903
protein-coding genes plus 1,026 non-targeting pairs). They were produced with the Fix-Cryopreserve ("FiCS") workflow:
DSP fixation, viability and dCas9/sgRNA sorting, cryostorage, superloaded 10x GEM-X 5′ lanes and automated library
preparation ([[huang-2025-x-atlas-orion]], pp. 3, 5, 13). QC (Table 1, p. 7):

| | HCT116 | HEK293T |
|---|---|---|
| GEM groups | 109 | 223 |
| Aligned cells | 9,886,825 | 21,550,681 |
| Mean reads per cell | 54,826 | 63,426 |
| Median UMIs per cell | 19,416 | 16,557 |
| Median genes per cell | 5,387 | 5,871 |
| Median CRISPR UMIs per cell | 950 | 184 |
| Cells with the expected two sgRNAs | 3,409,169 | 4,534,299 |
| Median two-sgRNA cells per perturbation | 141 | 188 |
| Median on-target KD (two-sgRNA cells) | 75.4 % | 51.5 % |
| Cells with one sgRNA (KD) | 1,290,167 (73.3 %) | 2,326,883 (45.3 %) |
| Expressed-gene perturbations discernible (q < 0.05) | 48.36 % | 40.37 % |
| Expressed genes (99 % of NTC cp10k) | 11,181 | 12,573 |

The last row is from p. 19. The source calls it «the largest publicly available Perturb-seq atlas» (p. 1), with about
8 M two-sgRNA cells in total.

## History of research

- [[huang-2025-x-atlas-orion]] introduce it as a platform demonstration and a training resource for dose-aware
  foundation models. They benchmark it against Replogle K562 (from the [[replogle-nadig-crispri-collection]] source
  paper) and [[tahoe-100m]] on depth and NTC batch reproducibility (pp. 5–8). They train no models.
- [[huuloc-2026-cell-quota]] re-analyse both lines for per-perturbation cell quotas: per-cell noise σ² = 0.913
  (HCT116) and 1.033 (HEK293T), in the same band as Tahoe. Only 14.7 % / 35.3 % of knockdowns clear their detection
  floor (median bias-corrected magnitude m 0.14 / 0.35 at 155–204 cells), so they call genome-wide CRISPRi
  detection-limited.

## The evidence

- NTC pseudobulk batch-to-batch Spearman: median 0.993 (HCT116), 0.988 (HEK293T) vs 0.967 (Replogle K562). Random
  HCT116 × K562 batch pairs give 0.74 (pp. 5–7; Fig. 2F; Supp. Fig. 1E).
- Depth: HCT116 median UMIs 1.68× Replogle K562 (11,499) and 8.45× Tahoe-100M (2,298) (p. 5).
- Knockdown varies within a perturbation with sgRNA abundance (R = 0.91 over KD bins), and high-sgRNA cells have
  1.66× (HCT116) / 5.21× (HEK293T) larger e-distance on essential genes (pp. 10–11; Fig. 4).

## Contested

- **How much of the atlas is signal.** The source's classifier calls 48.36 % (HCT116) of expressed-gene
  perturbations discernible, while [[huuloc-2026-cell-quota]] find 14.7 % above their floor. The criteria (single-cell
  classifier vs pseudobulk magnitude), the denominators (expressed genes with ≥ 100 cells vs all knockdowns) and the
  line ordering differ. Neither number is a split-half reproducibility estimate.
- **Line comparability.** The same vectors give 75 % vs 52 % knockdown (p. 10), so HCT116 vs HEK293T response
  differences mix biology with CRISPRi efficiency.

## Relevance to this project

- **H1 panel.** The OR_HCT116 / OR_HCT116X harness panels are built from the HCT116 screen: 263 official targets ×
  163 cells, about 12k NTC; OR_HCT116X also excludes HEK293T from the sources (`AGENTS.md`). The panel's r_max is
  expected to be low and noisy because many genome-wide knockdowns are weak.
- **H2c.** HCT116 and HEK293T are two of the 11 Count-1 cell types in C_public
  (`output/data-analysis/results/families.json`).
- **Source registry.** Both lines are M1 source families (HCT116, HEK293T).

## Multiple roles

The atlas is a destination (OR_HCT116X), a source family for other destinations, and a test case for detectability and
cell-quota arguments.

## Open / to verify

- Whether the project's OR_HCT116 cells are restricted to two-sgRNA same-pair cells (the Figshare release is).
- Whether the 41,780-guide library figure (p. 5) implies alternative guide pairs for some genes.
- Which detectability criterion the project's per-target noise estimates align with.

## Sources

- [[huang-2025-x-atlas-orion]] — pp. 1, 3, 5–8, 10–13, 18–19, 23; Table 1, Figs. 2–4.
- [[huuloc-2026-cell-quota]] — reanalysis of detectability and cell quotas (its Tables 1 and 4).
