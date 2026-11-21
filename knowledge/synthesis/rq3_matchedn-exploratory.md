---
title: "RQ3 matched-n sensitivity for H3a and H3b (exploratory, matched-n)"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
analysis_status: exploratory   # schema enum allows only draft|review|stable; exploratory label carried here + tags
author: llm
tags: [exploratory, matched-n, rq3, h3a, h3b]
---

# RQ3 matched-n sensitivity (exploratory, matched-n)

**Exploratory. Not pre-registered.** This page does not change the registered verdicts in
`output/data-analysis/results/rq3_h3a.json` (H3a CONFIRMED) and `rq3_h3b.json`, and those files were not touched.

- **Plan:** `input/ideas/zero-shot-limits-context-complexity-plan.md` (H3a, H3b, operational definition 7, Tasks 16–17).
- **Trigger:** the T15+T16 review in `knowledge/_meta/log.md` (2026-11-05). In the registered runs, U_j is built from
  **all** control cells of panel j (800 to 12,000 cells), while U_dest comes from a single third, C_pc (267 to 4,000 cells).
- **Registered notes:** [[rq3_h3a-assumptions]], [[rq3_h3b-assumptions]].

## What was run

**H3a.** Script: `output/data-analysis/rq3_matchedn_h3a.py`, which imports `rq3_energy.py` unchanged. Outputs:
`results/rq3_h3a_matchedn.{json,parquet,png}` and per-destination files in `results/rq3_matchedn/`.

Everything except the subspaces is identical to the registered run: residuals, target halves, C_a/C_b, the bootstrap
weights W (B = 2,000) and the two-level resampling. The variants:

- **reg:** the registered design, reproduced. It matches `rq3_energy.parquet` to ≤ 1.1e-16 on all 7 panels
  (point, lo and hi).
- **mn (matched n per comparison):** n′ = min(|C_pc|, n_j).
  - U_j comes from an n′-subsample of j's controls.
  - U_dest comes from C_pc, or from an n′-subsample of C_pc when j is smaller. That happens for H1, OR and KOLF
    against the four RN lines (800 cells).
  - 5 reps per (destination, seed, j).
- **c267 / c800:** every subspace is built from exactly 267 or 800 cells, with 5 reps.
  - c800 cannot be computed for the RN destinations, because C_pc has only 267 cells. It is reported as NA.
  - For the RN destinations, mn and c267 are the same design with independent draws.

**H3b.** Script: `output/data-analysis/rq3_matchedn_h3b.py`, which imports `rq3_corrections.py` unchanged. Output:
`results/rq3_h3b_matchedn.json`.

- Offline pack evaluator only; no cell-eval2.
- M3_dest and M3_other_j (6 lines) are rebuilt with every control-PC basis taken from n cells:
  - n = 800, which is both the smallest control set (the RN lines) and the requested 800;
  - n = 267 as an extra.
- 5 reps per line, using the builder recipe. Rep r of the destination is paired with rep r of each other line.
- The registered arms are reproduced exactly (Δ = 0 on 7/7). The builder-recipe destination basis at 6,000 cells equals
  the frozen ctrl_V to within 1.4e-5 PDS.

## H3a: D_other, registered vs matched n (k = 8, target-bootstrap 95 % CI)

| Panel | reg | mn | c267 | c800 |
|---|---|---|---|---|
| **H1** | +0.0182 [+0.0092, +0.0267] | **+0.0093 [+0.0003, +0.0176]** | +0.0011 [−0.0067, +0.0076] | +0.0085 [−0.0003, +0.0164] |
| OR_HCT116X | +0.0171 [+0.0100, +0.0244] | +0.0134 [+0.0065, +0.0206] | +0.0122 [+0.0046, +0.0204] | +0.0133 [+0.0063, +0.0208] |
| KOLF | +0.1059 [+0.0949, +0.1182] | +0.0746 [+0.0660, +0.0835] | +0.0295 [+0.0255, +0.0337] | +0.0618 [+0.0545, +0.0696] |
| RN_K562 | −0.0303 [−0.0356, −0.0250] | **+0.0236 [+0.0183, +0.0294]** | +0.0237 [+0.0184, +0.0295] | NA |
| RN_HepG2 | +0.0402 [+0.0229, +0.0559] | +0.0897 [+0.0780, +0.1027] | +0.0898 [+0.0779, +0.1028] | NA |
| RN_Jurkat | −0.0458 [−0.0541, −0.0366] | **+0.0101 [+0.0028, +0.0183]** | +0.0098 [+0.0026, +0.0181] | NA |
| RN_RPE1 | +0.1698 [+0.1565, +0.1827] | +0.2544 [+0.2370, +0.2705] | +0.2529 [+0.2357, +0.2689] | NA |

- **Two-level CIs** nearly coincide with these. For H1 mn the two-level CI is [+0.0001, +0.0177].
- **Rep-to-rep SD** of the point estimate is small: 0.0002 for H1 mn and ≤ 0.0015 everywhere.
- **k = 24 / 50 on H1:**
  - reg: 0.133 / 0.135
  - mn: 0.062 / 0.069
  - c800: 0.034 / 0.043
  - c267: 0.009 / 0.017

  All of these CIs are above 0.
- **D_other > 0 with the CI above 0:**
  - reg: 5/7
  - mn: **7/7**
  - c267: 6/7 (all except H1)
  - c800: 2/3 attainable panels (H1 lower bound −0.0003)

**Reading.**
1. **On H1 the mismatch inflated D_other by about 2×.** In the registered design, U_dest (4,000 cells) was compared
   with RN U_j built from 800 cells.
   - For the 4 RN lines, the per-line D drops from about 0.02 to about 0.002–0.009 once U_dest is cut to 800.
   - Against OR it is negative under both designs (−0.007 / −0.005), so OR's controls capture H1's residual slightly
     better at equal n.
   - At the H1 decision point (k = 8), the matched-n CI only just excludes 0. At common n = 267 or 800 it does not.
   - The registered H3a verdict is unaffected (registered rule, registered design). The H1-specific D_other margin is
     **fragile to PC-estimation n** and should be stated that way.
2. **The negative RN_K562 and RN_Jurkat D_other was an artefact of the mismatch.** In the registered design, U_dest came
   from 267 cells and U_j from 800 to 12,000. At matched n both are positive with CIs above 0, so the "other lines' PCs
   carry more residual energy" statement in [[rq3_h3a-assumptions]] does not survive at equal n.
3. **Control-PC capture depends steeply on n.** For example, F_dest on H1 is 0.084 at 4,000 cells and 0.044 at 267.
   Any F comparison across subspaces therefore needs matched n.

## H3b: ΔPDS (M3_dest − M3_other), 7-panel mean, paired target bootstrap

| Design | ΔPDS | 95 % CI | Panels > 0 | Rep SD |
|---|---|---|---|---|
| Registered (6,000-cell recipe, so 800 for RN) | +0.0058 | [+0.0039, +0.0076] | 7/7 | — |
| **Matched n = 800** | **+0.0058** | **[+0.0042, +0.0075]** | **7/7** | 0.0002 |
| Matched n = 267 | +0.0051 | [+0.0036, +0.0065] | 6/7 (H1 −0.00001) | 0.0005 |

M3_dest − none at matched n:
- n = 800: +0.0128 [+0.0101, +0.0154] (registered +0.0131)
- n = 267: +0.0105 [+0.0084, +0.0126]

Per panel at n = 800, H1 has the smallest margin: +0.0005. Single other lines still beat M3_dest on a few panels:
- H1 dest with KOLF PCs: 0.9109 vs 0.9103
- RN_Jurkat dest with OR PCs (0.9299) and RN_RPE1 PCs (0.9295) vs 0.9287

**Reading.** The H3b M3-vs-M3_other margin does **not** depend on the cell-count mismatch. It is unchanged at n = 800
and shrinks by about 12 % at n = 267.

## Assumption checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| H3a 'reg' reproduces `rq3_energy.parquet` (F_dest, F_other_mean, D_other; both bootstraps; 3 k) | ≤ 1e-8 | ≤ 1.1e-16 on 7/7 | yes |
| H3b registered M3_dest / M3_other / none reproduce `rq3_h3b.json` | ≤ 1e-9 | 0.0 on 7/7 | yes |
| Builder-recipe dest basis (6,000 cells) vs frozen ctrl_V | ≈ 0 | ≤ 1.4e-5 PDS | yes |
| Denominator > 0 (all variants, point) | yes | den_sign_ok on all rows | yes |
| Subsampling variability (rep SD) small vs CI width | yes | H3a ≤ 0.0015, H3b ≤ 0.0005 | yes |

## Runtime
Detached units, nice 10, 4 threads, on a loaded machine.

| Unit | Content | Runtime |
|---|---|---|
| `paper-mn-h3a-rn` | 3 RN panels (K562, HepG2, Jurkat) | 13 min |
| `paper-mn-h3a-rpe1` | RN_RPE1 re-run | 2 min |
| `paper-mn-h3a-big` | H1, OR_HCT116X, KOLF | 41 min |
| `paper-mn-h3b` | H3b, 7 panels | 51 min |

- `paper-mn-h3a-rn` was stopped by systemd-oomd (machine-wide memory pressure) while working on RN_RPE1, so RN_RPE1 was
  re-run in its own unit.
- Wall clock: about 53 min (15:19–16:12).
- Logs: `results/logs/paper-mn-*.log`.

## Evidence base
- `output/data-analysis/results/rq3_h3a_matchedn.{json,parquet,png}`
- `output/data-analysis/results/rq3_h3b_matchedn.json`
- The registered `rq3_energy.parquet` and `rq3_h3b.json`, used for the reproduction checks.

## What stays open
- **Bootstrap CIs** are conditional on the drawn subsamples. Subsampling variability is reported as the rep SD and is not
  folded into the CI. PC-estimation noise from C_pc is still in neither CI.
- **c800** is not attainable for the RN destinations without breaking the C_pc / C_a / C_b disjointness.
- **The manuscript should report** the H1 D_other as positive but n-sensitive: about +0.009 at matched n, versus +0.018
  registered.
- **The RN_K562 / RN_Jurkat "negative D_other" claim** should be retracted as a mismatch artefact, with a pointer to this
  page.
