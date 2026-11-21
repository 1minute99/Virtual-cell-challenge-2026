---
title: "Task 19 — official VCC record (E4): trajectory, rank, single-change gains, shrinkage"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# Task 19: official record

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 19 and Amendment A3 (38 official submissions =
36 validation + 2 final). VCC rules: `input/notes/external-gates.md` §1. This note reports only our own scores and
statistics derived from them. Leaderboard snapshots appear only as aggregate quantiles plus our rank. No other team
is named or identified.

**Code / outputs**
- `output/data-analysis/e4_official.py` (deterministic, no seed needed)
- `output/data-analysis/results/e4_official.csv`: one row per submission, with its official parent, the number of
  levers changed, the evidence, the lineage confidence, the attribution class, Δ composite and Δ per component
- `output/data-analysis/results/e4_official.json`: summary, snapshot aggregates, checks
- `output/data-analysis/results/e4_official_trajectory.png`: left panel, score per submission with best-so-far;
  right panel, our rank as a % of teams for each snapshot
- Log: `output/data-analysis/results/logs/e4_official.log`

## Assumptions stated
E1–E7 are in the script header. The load-bearing ones:
- **E4.** A snapshot shows each team's **latest** entry, not its best.
- **E6.** A gain counts toward a component only when the child differs from its nearest **officially scored** parent
  by one lever. Lever counts come from a recursive diff of the two build specs (`<label>.pred_build.json`) where
  both exist. Otherwise they come from the submission notes, AGENTS.md or PROJECT_HISTORY.md, and are marked as
  lower confidence.

## Assumption checks
| Check | Expected | Got | Pass |
|---|---|---|---|
| E1 submissions | 38 = 36 val + 2 final | 38 / 36 / 2; no entry id for B10, B11, B14 | ✓ |
| E2 composite = mean of the 6 components | 0 | max abs diff 8e-17 | ✓ |
| E4 one row per team in each snapshot | yes | yes, all 10 snapshots | ✓ |
| Final rank vs AGENTS.md line 25 ("S18 0.2545 ended rank 129; top-10 ≥ 0.319") | 129 / 1,635 | 129 / 1,635, displayed entry = S18, 0.25445, top-10 threshold 0.31933 | ✓ |
| Our submission count on the final board | 36 | 36 | ✓ |
| S4 − S1 vs submission_2026-09-22/FINDINGS.md | +0.001310 | +0.001310 | ✓ |
| S3 − S1 vs AGENTS.md | −0.000501 | −0.000501 | ✓ |
| b18 (byte-identical to B16) | Δ = 0 | Δ = −1e-16 (the scorer is deterministic) | ✓ |

## Results

**Trajectory.** The first entry scored −0.0228 (prior A). The best validation entry is **S18 = 0.25445**. The
best-so-far composite rose in steps:
- B8 0.122 (delta ×3)
- b19 0.136
- S4 0.139
- S10b 0.192 and S10c 0.217 (dual-moment emission and the cpc pin)
- S15 0.246
- S18 0.254

**Standing.** At the end of the validation phase (snapshot 2026-10-27) we were **rank 129 of 1,635 teams**, i.e.
the top 7.9 % (we scored above 92.1 % of teams). The aggregate scores on that board were:

| statistic | score |
|---|---|
| top-10 threshold | 0.3193 |
| 99th percentile | 0.3106 |
| 90th percentile | 0.2392 |
| median | 0.0807 |

Our gap to the top-10 threshold was 0.0649. The 129 / 1,635 is the **final validation standing, not our best rank**.
Our best rank fraction was 34 / 1,244 (2.7 %, S14fk, snapshot 09-30b). Rank fell to 252 / 1,434 on 10-10 because
the board showed our latest entry, the S17u budget probe, rather than our best.

| snapshot | teams | our rank | % | top-10 threshold | median |
|---|---|---|---|---|---|
| 09-18 | 1,036 | 205 | 19.8 | 0.229 | 0.048 |
| 09-22 | 1,104 | 266 | 24.1 | 0.242 | 0.056 |
| 09-27 | 1,171 | 52 | 4.4 | 0.253 | 0.068 |
| 10-01 | 1,250 | 35 | 2.8 | 0.273 | 0.078 |
| 10-10 | 1,434 | 252 | 17.6 | 0.294 | 0.090 |
| 10-27 | 1,635 | 129 | 7.9 | 0.319 | 0.081 |

**Single-change official gains (Δ composite vs the official parent).** All rows below have high confidence (the
lever was checked by a build-spec diff).

| pair | change | Δ |
|---|---|---|
| **S14fk2 → S15 (M2)** | agreement amplitude | **+0.0106** (MSE component +0.0559, PDS +0.0046, reach +0.0034, NMAE −0.0008) |
| S10b → S10c | cpc composition pin | +0.0247 |
| S4mix → S8 | neighbour blend 0.2 | +0.0041 |
| S14fk → S14fk2 | DE-source bug fix | +0.0027 |
| S15 → S16 | learned per-target amplitude | +0.0014 |
| S11 → S12 | DE side from the sum override | +0.0008 |
| S10c → S10e | S4 half swapped | −0.0046 |
| S4mix → S7a | fusion share 0.75 | −0.0063 |
| S16 → S17d / S17u | sum-budget probes | −0.0073 / −0.0462 |

Lower-confidence single-change pairs (lineage taken from the notes only):

| pair | change | Δ |
|---|---|---|
| B6 → B8 | delta ×3 | +0.0447 |
| B1 → B2 | 3-source mean | +0.0084 |
| b18 → b19 | context-forcing residual | +0.0073 |
| B8 → B16 (parent inferred) | | +0.0075 |
| B8 → B15 (parent inferred) | | −0.0150 |
| B6 → B7 | | −0.0994 |

S1 → S2 (−0.0020) and S1 → S3 (−0.0005) are single changes but **invalid**: they used scrambled context vectors
(see `vcc2026_data/reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`).

**The M3 bundle is not attributable to M3.**
- S17 was built (`submission_2026-09-23/S17.pred_build.json`) but **never uploaded**, so the official S17 → S18
  pair does not exist.
- S17 → S18 differs at the spec level only by `sum.ctrl_defl` (k 8, η 0.75, weaker half).
- The scored pair is **S16 → S18 = +0.0070** (PDS +0.0342, NMAE +0.0189, MSE −0.0053, reach −0.0055). It bundles
  four changes:
  1. four source families (HEPG2, RNJURKAT, RPE1, DLD1)
  2. the LFC-pseudocount ensemble
  3. DE calibration drc 0.8 / 1.1
  4. M3 deflation
- The S17 bundle alone was worth +0.0055 R21 on the local harness, a different scale from the official score.

Every other multi-change pair is flagged `MULTI_CHANGE_not_attributable` in the CSV:
B2→B6, B8→B10, B8→B11, B8→B14, S1→S4, S4→S4gp, S4gp→S4mix, S8→S10b, S10c→S10f, S10f→S10k, S10k→S11, S12→S14fk
and S16→S18.

**Local → official shrinkage.** Both contrasts are against the refit S1.

| case | local Δ | official Δ | ratio | valid? |
|---|---|---|---|---|
| S4 (refit_v2 + FiLM, five-line local composite, 15/15 seeds) | +0.0247 | +0.00131 | **0.053** (~19× smaller; sign held) | yes. Note the bundle: refit_v2 alone was +0.0060 locally |
| S3 (scratch FiLM) | +0.0180 | −0.00050 | −0.028 | **no**. The local number comes from the drug-context-pretraining experiment and both arms used scrambled gene-context vectors |

**Final phase.** F1 was uploaded 2026-10-27 and superseded. F2 (2026-10-28) is the counting final entry. Its score
is **not yet published**: the rules say the final leaderboard is released with the winners in late November 2026.

## Evidence base
- `output/data-analysis/results/official_record.csv` (Task 4; see [[task04-freeze-inputs-assumptions]])
- `vcc2026_data/reports/leaderboard_*.json` (10 snapshots, reduced to aggregates)
- `vcc2026_data/reports/submission_2026-09-23/*.pred_build.json`
- `vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`
- `vcc2026_data/reports/PROJECT_HISTORY.md` (table "Official submissions, all 14")
- `AGENTS.md`: lines 25, 736–745 and 811–816, plus the BUG section at about line 780

## What stays open
- B-series lineage comes from one-line notes. B15 and B16 have an inferred parent (B8).
- Official deltas are single uploads with no seed replicates. The noise floor is only bounded by b18 ≡ B16
  (the scorer is deterministic), not by model-seed variance.
- An M3-only official contrast does not exist. The offline M3 effect is RQ3/Task 17's job (see
  [[rq3_h3b-assumptions]]).
- The final-phase score is pending publication by the organisers.
