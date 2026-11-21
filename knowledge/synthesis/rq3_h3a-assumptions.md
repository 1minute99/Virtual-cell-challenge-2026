---
title: "RQ3 Task 16 — E1a / H3a run (residual energy on control PCs): assumptions, checks, verdict"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ3 Task 16: E1a / H3a

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` — H3a, its falsification criterion, operational
definition 7, Task 16. No amendment (A1–A3) changes RQ3. The pre-data clarification in `knowledge/_meta/log.md`
("Pre-data operationalisation for Task 16") fixes two things:
- the H3a decision uses the registered target bootstrap on H1; every panel also reports a two-level bootstrap;
- the negative-control rule uses the numerator form.

Upstream notes: [[rq3_energy-assumptions]] (Task 15 estimator) and [[rq3_build_preds-assumptions]] (Task 14 predictions).

**Code**
- `output/data-analysis/rq3_energy.py`: the estimator, unchanged, plus the Task 16 additions T16-1 to T16-4 in its header.
- `output/data-analysis/rq3_h3a.py`: the verdict and the summary.
- Tests: `output/data-analysis/tests/test_rq3_energy.py`, 8/8 pass. Six are Task 15 tests. Test (vi) checks that the
  new summariser reproduces the old one exactly and that its two-level row 0 equals the point estimate. Test (vii)
  checks that the two-level bootstrap is wider on 100-cell control thirds.

**Outputs**
- `output/data-analysis/results/rq3_energy.parquet`
- `output/data-analysis/results/rq3_negctrl.parquet`
- `output/data-analysis/results/rq3_h3a.json`
- Per-run files in `results/rq3/`: `<mode>_<panel>_<pred>[_full].{parquet,json}` plus `_parts.npz` (per-target
  projections, for re-analysis without re-reading cells).
- Logs: `results/logs/paper-t16{,-full}.log`, `results/logs/rq3_t16_calib_{300,800,3000}.log`.

## Assumptions stated

1. **Primary input.** `results/preds/<panel>_prem3_ctrlAB_seed{0..4}.npz`. These are the pre-M3 predictions that
   Task 14 rebuilt from C_a ∪ C_b only, so the prediction is independent of C_pc and of U_dest. The splits are
   `<panel>_ctrl_split_seed{0..4}.npz`, reused and validated.
2. **Gene space.** G∩ = 6,014 genes. k ∈ {8, 24, 50}, nested.
3. **Subspaces.**
   - U_dest: from C_pc.
   - U_j: from all control cells of each of the other 6 panels.
   - U_perm: 20 gene-wise permutations of C_pc.
   - U_haar: reported only.
4. **Target bootstrap (registered).** B = 2,000 over targets. In each draw the same resampled target set is used for
   every split seed; the result is the seed mean of the 5 split seeds. 95 % percentile CI. Seed 20261105.
5. **Two-level bootstrap (T16-1).**
   - Each draw uses the same target draw as in item 4. Independently for each split seed, it also resamples the cells
     of C_a and of C_b with replacement (seed `[20261105, split_seed, 17]`).
   - The perturbation of each baseline enters exactly in both the numerator and the denominator.
   - C_pc is not resampled, so the noise from estimating the PCs is in neither CI. That noise varies across the
     averaged split seeds.
6. **Negative control.**
   - 300 pseudo-targets per split seed, drawn from C_a ∪ C_b, with Δpred = 0. Sizes are drawn from the panel's real
     per-target cell counts.
   - Each pseudo-target has leave-out baselines C_a∖S and C_b∖S.
   - Numerator form: D_neg = N_neg / (observed denominator per target). It is compared with 0.5 × D_obs, separately
     for D_other and D_perm.
   - Only the target bootstrap is used here (T16-4).
7. **H3a rule.** H3a is rejected if either of these holds:
   - the target-bootstrap CI of D_other or of D_perm at k = 8 on H1 is not entirely above 0;
   - D_neg ≥ 0.5 × D_obs for D_other or for D_perm (point estimates).

   Otherwise H3a is confirmed. The plan names no "inconclusive" branch for H3a.
8. **Strata.** Quartiles of ‖Δpred‖ on G∩. The amplitude comes from the main `<panel>_prem3.npz` (Q1 = smallest);
   each quartile holds 65–75 targets.
9. **Secondary analysis (T16-2): destination only, full gene space.**
   - Gene axis: the destination's full truth var_names (7,678–18,533 genes). Task 14 writes Δpred on that axis.
   - Subspaces: U_dest, U_perm and U_haar only, so this analysis has D_perm but no D_other.
10. **Loader fix (T16-3).** Task 14 stores `genes` as an object array, so `load_pred` now uses `allow_pickle=True`
    (the files are the project's own).

## Assumption checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Unit tests (i)–(v), smoke, (vi), (vii) | pass | 8/8 | yes |
| summarize2(two_level=False) vs Task 15 `summarize` | identical | rtol 1e-8 on point / lo / hi / se | yes |
| Two-level point vs target point (all rows) | identical | asserted in `rq3_h3a.py` | yes |
| F_haar at k = 8 vs k/G = 0.00133 | ≈ | 0.0012–0.0014 on all 7 panels | yes |
| Denominator > 0 in every seed (point) | yes | `den_sign_ok` true on all real runs | yes |
| Two-level calibration, synthetic, true D = 0, 2 seeds, 20 replicates (`rq3_t16_calib.py`) | ≈ 95 % coverage | see the table below | partly |

Two-level calibration. Each cell gives the CI half-width ÷ 1.96, then the coverage of the true D = 0.

| Control thirds | sd(point) | Target bootstrap | Two-level bootstrap |
|---|---|---|---|
| 100 cells | D 0.158, N 0.106 | 0.031, cover 0.50 | D 0.43, N 0.17; cover 1.00 |
| 266 cells (RN geometry) | D 0.061 | 0.034, cover 0.60 | 0.091 (≈ 1.5×), cover 1.00 |
| 1,000 cells | D 0.035 | 0.029, cover 0.90 | 0.038, cover 0.95 |

- With small control thirds, the target bootstrap is anti-conservative. The two-level bootstrap is conservative,
  because the bilinear baseline term is degenerate: it picks up the realized ε × δ cross-terms.
- The truth lies between the two CIs. On the real panels the two CIs nearly coincide, even on RN (e.g. RN_HepG2
  D_perm: target [0.168, 0.219], two-level [0.164, 0.220]). Between-target variance dominates there, so the
  calibration issue does not matter in practice.

## Results (k = 8 unless stated; full table in `rq3_h3a.json`)

**H3a on H1: CONFIRMED.**

| Quantity | Point | Target 95 % CI | Two-level 95 % CI |
|---|---|---|---|
| F(U_dest) | 0.0842 | [0.0711, 0.0983] | |
| mean_j F(U_j) | 0.0659 | [0.0531, 0.0826] | |
| F(U_perm) | 0.0060 | | |
| F(U_haar) | 0.0014 | | |
| D_other | 0.0182 | [0.0092, 0.0267] | [0.0090, 0.0269] |
| D_perm | 0.0782 | [0.0653, 0.0922] | [0.0651, 0.0925] |

Negative control on H1 (numerator form):

| | D_neg | 0.5 × D_obs | D_neg as % of D_obs |
|---|---|---|---|
| D_other | 0.00034 | 0.0091 | 1.8 % |
| D_perm | 0.00053 | 0.0391 | 0.7 % |

The ratio form is unstable, as Task 15 predicted: 0.115 / 0.221, with CIs of about ±3.

At k = 24 and k = 50, D_other on H1 grows to 0.133 / 0.135, while F(U_j) saturates at about 0.09–0.10.

**Cross-panel** (two-level CIs; target CIs give the same sets at every k):
- D_perm > 0 on **7/7** panels at every k.
- D_other > 0 on 5/7: H1, OR_HCT116X, KOLF, RN_HepG2, RN_RPE1.
- **D_other < 0, with CIs excluding 0,** on RN_K562 (−0.030 [−0.037, −0.024]) and RN_Jurkat (−0.046 [−0.054, −0.036]).
  On those two panels the other lines' control PCs carry more residual energy than their own PCs. The large
  contributions come from OR_HCT116X and the other RN lines (F 0.18–0.25 vs F_dest 0.13–0.14); H1 and KOLF
  contribute ≈ 0.05–0.08.
- **Both D > 0: 5/7 at k = 8, 24 and 50.**

**Negative control on the other panels.** Everywhere D_obs > 0, D_neg < 0.5 × D_obs.
- Largest: OR_HCT116X, at 31 % of D_other (CI 0.0018–0.0087 vs D_obs 0.0171) and 21 % of D_perm.
- All other panels: ≤ 3 %.
- On RN_K562 and RN_Jurkat the rule is vacuous for D_other, because D_obs < 0 (flagged `rule_applicable: false`).
- The noise-only numerator N_perm is small but has a CI above 0 on several panels. The most likely cause is the
  O(1/m) log-pseudobulk bias, which is shared by the two halves of a small target (A6 is not exactly negligible).
  The same bias is present in the observed D, so on OR the bias could account for up to about a third of D_other.

**Strata** (Q1 = smallest predicted amplitude). F_dest, D_perm and D_other mostly fall with amplitude:

| Panel | D_perm, Q1 → Q4 |
|---|---|
| KOLF | 0.41 → 0.11 |
| RN_HepG2 | 0.40 → 0.16 |
| RN_RPE1 | 0.67 → 0.40 |
| OR_HCT116X | 0.11 → 0.04 |

- On H1 the strata are flat: D_perm 0.07–0.10, D_other 0.031 → 0.011, and the Q4 two-level lower bound is ≈ 0.
- On RN_K562 and RN_Jurkat, D_other is most negative in Q4.
- So the destination-control share of the residual is largest for weakly predicted targets.

**Secondary analysis: destination only, full gene space.** D_perm > 0 on **7/7** panels at every k, with two-level
CIs excluding 0. At k = 8 the magnitudes are close to those on G∩:

| Panel | D_perm (full space) | Two-level 95 % CI | D_perm on G∩ |
|---|---|---|---|
| H1 | 0.092 | [0.076, 0.110] | 0.078 |
| OR_HCT116X | 0.045 | | 0.050 |
| KOLF | 0.192 | | 0.169 |
| RN_K562 | 0.124 | | 0.136 |
| RN_HepG2 | 0.185 | | 0.192 |
| RN_Jurkat | 0.111 | | 0.125 |
| RN_RPE1 | 0.438 | | 0.482 |

The finding is therefore not an artefact of restricting to G∩.

## Runtime

Two detached units ran in parallel, with nice 10 and 4 threads, on a loaded machine (load 15–26):
- `aivc-paper-t16` (G∩, real + negctrl, 7 panels): about 33 min. Per panel: H1, OR and KOLF 4–6 min each per mode;
  RN panels under 1 min.
- `aivc-paper-t16-full` (full gene space, real only): 38 min, 14:17–14:55. H1, OR and KOLF took 9–14 min each.
- Tests: about 16 s.

Peak memory of the G∩ unit was 5.7 GB (systemd MemoryPeak).

## Evidence base

- Real-data runs above.
- Synthetic calibration: `output/data-analysis/rq3_t16_calib.py`, with logs `results/logs/rq3_t16_calib_{300,800,3000}.log`.
- Task 14 inputs: `rq3_build_preds_summary.json`.

## What stays open

- **OR_HCT116X negative control (31 % of D_other).** The log-bias component could be removed with a bias-corrected
  pseudobulk or equal half sizes for target and baseline. That would be exploratory.
- **The two-level bootstrap is conservative with small control thirds.** A calibrated alternative, such as subtracting
  the ε × δ cross-term variance, was not pre-specified.
- **PC-estimation noise (C_pc)** is in neither CI.
- **D_other < 0 on RN_K562 and RN_Jurkat.** The control PCs of other lines measured with the same technology (Replogle
  and Nadig 10x, Orion) capture shared technical axes better than the destination's own 266-cell C_pc. This is
  relevant for H3b's M3_other arm (Task 17).
