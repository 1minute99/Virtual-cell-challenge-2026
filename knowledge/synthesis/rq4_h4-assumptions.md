---
title: "RQ4 Task 18 — same-line data value vs modelling gains (H4): assumptions and checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ4 Task 18: same-line data value vs modelling gains (H4)

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, H4, its falsification criterion, and Task 18.
Builds on [[rq3_build_preds-assumptions]] (Task 14 frozen recipe).

- **Script:** `output/data-analysis/rq4_sameline.py`. Subcommands: `build`, `real ARM…`, `verdict`.
- **Queue:** `output/data-analysis/rq4_queue.sh`, unit `aivc-paper-t18-queue`.
- **Outputs:**
  - packages: `results/rq4/pkg/<arm>/{pred.h5ad, inputs.npz, meta.json, result.json}`
  - the DE delta: `results/rq4/de_delta_H1.npz`
  - `results/rq4/build_report.json`
  - `results/rq4_sameline.parquet`
  - `results/rq4_h4.json`

**Pre-registered H4.** On the H1 panel with the real scorer, adding H1GW (w = 1) gains more composite than the largest single
modelling gain (M2, M3, or the best of the 17 architectures relative to its transfer baseline). It is *rejected* if
composite(base + H1GW, w = 1) − composite(base) ≤ that gain.

## Arms (H1 panel: 300 targets × 300 cells, 12,000 controls)

| arm | sums | DE side |
|---|---|---|
| `w0` = base | Task 14 recipe: harness M1 + F2 sum side, incl. M2 and M3 | F2 DE side (below) |
| `w0.25`, `w0.5`, `w0.75`, `w1` | base + H1GW same-line replacement at weight w | same replacement, then drc |
| `noM2` | base without the `--sum-tamp 0.5` amplitude (M3 kept) | = base |
| `noM3` | base without M3 (= Task 14 pre-M3) | = base |

- Verdict arms: `w0`, `w1`, `noM2`, `noM3`.
- The scorer runs them first, then `w0.5`, `w0.25` and `w0.75`.

## Assumptions (stated in the script header)

1. **Sums = Task 14 verbatim.**
   - The helpers are copied from `rq3_build_preds.py`.
   - The budget Bt is the ERA-pack Bt. Same-line data does not move it: in the builder, the budget reference D and the overrides OV are untouched.
2. **DE side = the F2 build** (`tv2_submit_build`, dual_mult).
   - *Delta:* Dde is the same M1 recipe, built on the **denoised** source library.
     - F2 sets `--de-denoise intended`. The harness mirror is `denoise_src/dnK100r0.25` (token `dsrcdnK100r0.25`).
     - The recipe: S4 Yorig + an 8-member min_cells ensemble + a 0.3 neighbour blend, with H1 excluded.
   - *Scale:* the DE side sits at γ 1; the budget scale applies to the sums only.
   - *Transforms, in order:* the same-line replacement (the builder applies it to Dover too), then drc 0.8 / 1.1.
   - *Emission* (`sch_emit.emit`, arm `em_mult_share_k3000_du0.5_ag_cpc_ci5_lsb`, seed 0):
     - K 3000, depth_up 0.5, cpc, ci 5, lsb;
     - strength order = |lfc_de| × consensus agreement, using the builder formula with the own family left out (`agc`).
   - *Overrides:* OV is the main recipe's, by the builder convention. The denoised recipe's OV differs in 2 of 389 entries.
   - *Left out:* F2 M1 items that are not reproducible on the harness are omitted on both sides, as in Task 14: cens, DLD1, promoter-off, knn_noK, resid, S4noFK.
3. **Same-line formula** = builder `same_line`, verbatim, on the official axis:
   - d ← (1 − w) d + w v ‖d‖/‖v‖;
   - applied to targets with ≥ 20 source cells, on the genes the source measures, excluding override entries and the target's own gene;
   - source: H1GW = `transfer_v2_2026-09-23/sources_new/H1GW.npz`, sha256 `417c861c…`, in the frozen manifest.
4. **The M2 rank and the M3 weak half do not change with w.** Both are source-side.
5. **Composite** = `tv2_score.composite` (vcc2026 constants), one run per arm.
   - Emission-seed noise reference: the earlier harness H1GW w1 arm scored 0.3772 / 0.3769 / 0.3768 at seeds 0 / 1 / 2, a range of 0.0004.
6. **Best architecture** (`arch_screen_2026-09-22/final_ranking.json`, `<line>_c`). Fixed before any Task 18 real result.
   - *Primary:* the largest line-mean over K562 and HepG2. This is the **additive shift, +0.0036**.
   - *Sensitivity:* the largest single cell over 34 line × architecture cells. This is the **LoRA hypernetwork on HepG2, +0.0169**.
   - *Setup differences:*
     - lines K562 / HepG2, not H1;
     - base = the refit_v2 transfer model, not the fusion pipeline;
     - a different emission;
     - seed 0, α 1;
     - same scorer preset.
   - Task 10 re-checks this number.
7. **Per-context value:** `null`. It comes from Task 11 (LOFO converted to composite on H1), which has not run yet. The H1GW / per-context ratio is therefore also `null`.

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Sums recomputed vs Task 14 `H1_m3` / `H1_prem3` lfc | rel < 1e-5 | 4.0e-8 / 4.6e-8 | yes |
| Harness `s4Yorig_ensP1` delta vs ERA `base_total` | identical | max abs 8.8e-8, cos 1.0, Od identical | yes |
| Denoised DE delta: cos to the sum recipe | 0.9723 (earlier harness dsrc arm) | 0.9723 | yes |
| Denoised vs main override pattern | small | 2 / 389 entries differ; main OV used | noted |
| Emitted per-target column sums = round(fm e^L n) | exact | asserted on 3 targets per arm | see build log |
| Harness controls = Task 14 controls (fm) | equal | equal | yes |

## Operations

- **The first build was killed by systemd-oomd.** User-session memory pressure reached 51.9 %, alongside the T17 scorer (22 GB) and T9 (13 GB).
- **Fix:** the build was made lean (one control load, one arm in memory at a time).
- **It now runs inside the queue, under the scorer lock, after `aivc-paper-t17-real` exits.** The packages are therefore built just before scoring, not "now".
- **Every `cell-eval2` call goes through** `flock results/.celleval.lock` with `--preset vcc2026 --set de.backend=pdex`.
- **Real-moments cache:** `results/rq3_cache/cache_real_H1`, the one written by Task 17.

## Evidence base

- `results/rq4/build_report.json`, the per-arm `meta.json` and `result.json`, and `results/rq4_h4.json`
- The Task 14 preds and their note [[rq3_build_preds-assumptions]]
- The ERA packs `score_H1` / `inputs_H1`
- `arch_screen_2026-09-22/final_ranking.json` and its FINDINGS
- `sc_hybrid_2026-09-25` FINDINGS §20.20: earlier harness same-line numbers, base 0.3029 → 0.3772 at w 1

## What stays open

- **The H4 verdict is PROVISIONAL** until the 4 verdict arms are scored. `rq4_h4.json` is rewritten after each run.
- **Single runs.** There are no seeds per arm. The harness emission-seed range is about 0.0004.
- **The best-architecture gain comes from a different panel, base and emission.** Task 10 updates it.
- **Related-line cases were skipped**, as the plan's cut-first item: KOLF ← H1GW and OR_HCT116 ← DLD-1 at w ∈ {0, 0.5, 1}. The reason is scorer-queue time, and KOLF would also need a kolf-free denoised DE library.
- **Per-context value / ratio:** waiting on Task 11.
