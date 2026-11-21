---
title: "RQ2 Task 9 — Proposition-2 parameters on the public genetic screens and H2c: assumptions and results"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ2 Task 9: genetic C*_min and H2c

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (H2c, operational definitions 3–5, Task 9, Amendments
A2/A3). Pre-specifications: `knowledge/_meta/log.md`, 2026-11-05 (d_gene over the Count-1 contexts, capped at C − 1;
H2c against both counts; secondary Count-2 contexts; secondary S1–S3). Estimator library: `output/data-analysis/rq2_params.py`
(Task 8, unchanged; see [[rq2_params-assumptions]]).

Scripts (run from `output/data-analysis/`, env `ct_transfer`, nice 10, 4 threads, detached):

| Step | Script | Output |
|---|---|---|
| 1. Count 2 (A2), before any C*_min | `rq2_gene_count2.py` | `results/rq2_params_gene.json["count2"]` |
| 2. Cell-level split-half extraction | `rq2_gene_extract.py plan`, `extract <SRC>` | `results/rq2_gene_cache/*.npz` (4.8 GB, derived) |
| 3. Estimates, sensitivities, verdict | `rq2_gene_params.py` | `results/rq2_params_gene.json`, `results/rq2_h2c.json` |

Runtime: plan 5 min; extraction ≈ 61 min wall time (4 parallel jobs; the Feng 4.5 GB gzip CSV took 61 min, the others
0.5–14 min); analysis 55.5 min. Peak RSS was about 14 GB during the analysis, above the ~12 GB guideline, inside the
job slice limits.

## Step 1: C_public, Count 2 (A2)

The Feng/HipSci targeted screen has 19 lines × 445 targets. Seventeen lines have ≥ 300 targets with ≥ 20 cells.
The two that do not are fiaj_3 (164 targets) and pipw_5 (198). Cell counts come from the RNA-matched
`feng_targetedscreen_target_line` store, and a direct metadata recount agrees. kolf_2 and kolf_3 qualify but are
merged with KOLF2.1J, which is already counted. **Count 2 = 11 − 1 (pooled HipSci family) + 15 = 25.**

- The genome-wide Feng arms qualify no line (at most 17 targets).
- Sensitivity, outside the A2 wording: counting clones of one donor once gives 19.

## Assumptions

- **G1 contexts.** There is one context per Count-1 line:
  - K562 (GWPS, the registry's K562 family)
  - RPE1, HepG2, Jurkat (Replogle/Nadig essential screens)
  - HCT116, HEK293T (Orion; all called cells)
  - CD4 (three states pooled)
  - KOLF2.1J (raw counts layer)
  - the HipSci pool (Feng targeted screen)
  - DLD-1 (singlets, control label `NTC`)
  - H1 = H1GW (guide call ≥ 5 UMIs, runner-up ≤ 20 %, as in `ml_h1gw_build.py`)

  VCC2025-H1 enters only the sensitivities.
- **G1b HipSci pool.** For each target, the pool's response is the 1/(1/n_t + 1/n_c)-weighted mean of the
  within-line responses. Each line is compared with its own control halves, so differences in line mix between
  target cells and control cells cannot create a response. The registry-style summed pseudobulk is a sensitivity
  and gives identical results to 3 digits.
- **E1/G3 gene axis.** 6,109 official-axis genes measured in every source. The intersection is limited by the
  ~8k-gene scPerturb matrices.
- **E2 library.** The native total UMI per cell.
- **E3 target universe.** 2,366 targets with ≥ 20 cells in at least one of the RPE1/HepG2/Jurkat essential screens.
  A target enters a line with ≥ 20 cells, and the analysis keeps targets present in ≥ 3 lines (filter_min_lines).
  Lines per target: median 8–9; 36 targets are in all 11 lines.
- **E4/G2 split halves.** Five split seeds. Within every group, including the control pool, cells are permuted with
  `default_rng([20261108, seed, src, group])` and floor(n/2) go to half A. **A and B have separate control
  halves**, with no shared baseline.
- **CD4 deviation.** No cell-level CD4 data are on disk; there are only run × donor × guide × condition pseudobulk
  units. CD4 halves are therefore unit-level, balanced greedily by cell count, and targets with
  |n − 2nA|/n > 0.25 in any seed are dropped (239 targets). The sensitivity that drops CD4 is reported below.
- **G4 descriptor.** Full-depth control pseudobulk in log1p CP10K on the axis. d is the participation ratio, capped
  at C − 1. The cap does not bind: d = 5.08 against a cap of 10. A participation ratio over 11 contexts is still a
  finite-sample quantity (on Tahoe, 12 contexts gave d 8.6 against 20.6 for all 48), so it is likely biased low.
  That direction is conservative for H2c, as pre-specified.
- **G5 depth rescaling.** σ²_ε is rescaled from half to full depth by × ½ (rq2_params A1). This is exact for
  cell-level floor/ceil halves; for CD4 the error is ≤ ~7 % at 25 % imbalance.
- **G6 R² and C*_theory.** They are reported because 1,801 targets have n_t − k − 1 ≥ 1 at k = ⌈d⌉ = 6. With
  ≤ 4 residual d.o.f. per target they are descriptive only.
- **Bootstrap.** B = 2,000 over targets × split seeds (seed 20261109); d and the additive fit are held fixed.
  S3: B = 1,000 line resamples (seed 20261110).

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Count 2 written before any C*_min | yes | written and logged first | pass |
| Half balance (cell-level sources) | ‖n − 2nA‖ ≤ 1 | max 0.048·n (n = 21) | pass |
| Half balance (CD4 units) | ≤ 25 % after filter | 239 / 1,878 targets dropped | pass (deviation) |
| Half B = full − A ≥ 0 | ≥ 0 | min 0 | pass |
| Own-gene knockdown (log effect, median) | < 0 | −0.16 (H1GW, KOLF) … −0.80 (VCC-H1); all 12 sources negative | pass |
| Extra stats vs halves_stats (Var(i), σ²_ε) | identical | rel. diff < 1e-13 | pass |
| Split-seed stability of C*_min | small | 16.39–16.51 | pass |
| d cap C − 1 | not binding | 5.08 < 10 | pass |
| Pooling artefact (H1, summed counts across platforms) | none | first run's summed H1GW + VCC pool inflated Var(i) to 1.60e-3; replaced by within-dataset responses (1.147e-3) | fixed before final run |

## Coverage per line (primary, after the ≥ 3-line filter)

| Line | Targets used | Median cells/target | Control cells |
|---|---|---|---|
| K562 (GWPS) | 2,344 | 175 | 75,328 |
| hTERT-RPE1 | 2,204 | 77 | 11,485 |
| HepG2 | 2,160 | 49 | 4,976 |
| Jurkat | 2,270 | 86 | 12,013 |
| HCT116 | 2,176 | 122 | 165,948 |
| HEK293T | 2,182 | 160 | 219,049 |
| CD4 (units) | 1,639 | 620 | 940,539 |
| KOLF2.1J | 1,899 | 148 | 147,058 |
| HipSci pool | 300 | 1,332 | 8,241 |
| DLD-1 | 580 | 109 | 16,751 |
| H1 (H1GW) | 1,717 | 216 | 81,617 |

The primary analysis has 2,366 targets, 19,471 (line, target) rows and 6,109 genes.

## Results (primary, C = 11)

| Quantity | Estimate | 95 % CI |
|---|---|---|
| d_gene (cap 10; LOO 4.46–6.63; jackknife SE 1.85) | 5.08 | — (fixed) |
| Var(i) | 1.145e-3 | 1.103e-3 – 1.187e-3 |
| σ²_ε (full depth) | 2.337e-3 | 2.296e-3 – 2.378e-3 |
| σ²_ε / Var(i) | 2.04 | 1.97 – 2.11 |
| **C*_min** | **16.45** | **16.09 – 16.82** |
| R²_adj (k = 6; descriptive) | 0.337 | 0.317 – 0.356 |
| C*_theory (descriptive) | 46.8 | 43.9 – 50.2 |

## H2c verdict

The verdict is **"holds for distinct cell types, not for donor lines"**:

- **Count 1 (C_public = 11): confirmed.** The lower bound 16.09 is above 11.
- **Count 2 (C_public = 25): rejected.** The point estimate 16.45 is ≤ 25.

## Sensitivities (C*_min, 95 % CI over targets × seeds)

| Variant | C | d | C*_min |
|---|---|---|---|
| Descriptor projected on the top 8 pooled-control-cell PCs | 11 | 4.31 | 14.12 [13.82, 14.43] |
| Top 24 PCs | 11 | 5.05 | 16.35 [15.99, 16.71] |
| Top 50 PCs | 11 | 5.05 | 16.36 [16.00, 16.72] |
| KOLF excluded | 10 | 4.46 | 14.13 [13.80, 14.47] |
| H1GW + VCC2025-H1 as one line | 11 | 5.16 | 16.67 [16.31, 17.04] |
| VCC2025-H1 as a 12th context | 12 | 4.12 | 13.48 [13.20, 13.77] |
| HipSci pool as one summed pseudobulk | 11 | 5.09 | 16.46 [16.11, 16.83] |
| CD4 excluded | 10 | 6.63 | 18.18 [17.89, 18.48] |
| **Secondary: Count-2 contexts (10 lines + 15 HipSci lines)** | 25 | 5.11 | 17.79 [17.40, 18.21] (< 25) |

Every variant has a lower bound above 11 and a point estimate below 25, so none changes the verdict.

About the "control PCs" sensitivity:

- "Top-k control PCs" was read as the PCs of the pooled control **cells** (≤ 2,000 per line). Under the other
  reading, the cross-context PCs of the descriptor, every k ≥ C − 1 = 10 equals the primary.
- In the secondary Count-2 analysis, R²_adj is 0.48 and C*_theory is 35.8.

## Secondary S1–S3 (pre-specified; they do not change the verdict)

- **S1 (target b + i).** Var(b) = 3.7e-4 (R²_b 0.70). C*_min(b + i) = 13.9 (13.9 with the empirical b noise) and
  C*_theory = 31.3.
- **S2 (line component).** C*_b,min = 6.1 and C*_b,theory = 8.3. The line component alone is below both C_public
  counts.
- **b noise.** The modelled b noise σ²_ε/n̄ (1.3e-6) is about 7.6 times below the empirical half-difference
  (1.0e-5): the shared control half does not average out over targets. The library reports both.
- **S3 (line bootstrap that recomputes d).**

  | Quantity | 95 % CI | Bootstrap median |
  |---|---|---|
  | d | [1.97, 5.23] | 3.84 |
  | C*_min | [7.95, 17.07] | 12.5 |
  | S1 C*_min | [6.5, 15.5] | — |
  | S2 C*_b,min | [3.0, 6.3] | — |

  With 11 lines a resample has about 7 distinct lines, so the duplicated rows bias the participation ratio
  downward, and this percentile interval is shifted below the point estimate. A jackknife normal interval on d
  (5.08 ± 3.6) instead gives C*_min ≈ 3.04·d + 1 ∈ [5.5, 27.5]. **Once d's uncertainty is propagated, the Count-1
  confirmation is not robust.** H2c's registered CI (targets × split halves) holds d fixed, so it does not show
  this.

## Evidence base

Cell-level sources:

- `vcc2026_data/raw/scperturb/{ReplogleWeissman2022_K562_gwps,ReplogleWeissman2022_rpe1,NadigOConner2024_hepg2,NadigOConner2024_jurkat}.h5ad`
- `raw/orion/data/*.parquet`
- `raw/kolf/KOLF_Pan_Genome_QC_Filtered.h5ad` (`layers/counts`)
- `raw/cd4/GWCD4i.pseudobulk_merged.h5ad` (units)
- `raw/feng/TargetedScreen_*`
- `raw/dld1_gse337988/*Low_assays.h5`
- `raw/h1_gw_gse295214/*.h5`
- `raw/vcc2025/adata_Training.h5ad`

Registry: `output/data-analysis/results/families.json`. All inputs were read-only, and no existing file was modified.

## What stays open

- **CD4 has no cell-level halves.** Its unit-level noise includes donor and run variation; dropping CD4 raises
  C*_min to 18.2.
- **d rests on 11 points.** Its uncertainty (jackknife SE 1.85) dominates every other source of uncertainty and is
  not in the registered CI.
- **The gene axis is about 6.1k genes**, limited by the scPerturb matrices.
- **The universe is essential-screen targets.** The genome-wide lines contribute only these targets.
- **The HipSci pool and DLD-1 cover fewer targets** (300 and 580).
- **Linear-Gaussian theory** underlies all of the estimates.
