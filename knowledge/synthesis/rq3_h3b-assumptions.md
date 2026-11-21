---
title: "RQ3 Task 17 — E1b / H3b correction comparison: assumptions and checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ3 Task 17: M3 vs alternative corrections (H3b)

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (H3b + falsification criterion, Task 17).
Script: `output/data-analysis/rq3_corrections.py` (`offline` → `real ...` → `verdict`). Inputs: [[rq3_build_preds-assumptions]].
Outputs: `output/data-analysis/results/rq3_h3b.json`, `results/rq3_corrections.parquet`, `results/rq3_corr/`.

**STATUS: PROVISIONAL.** The offline part (all 7 panels, bootstrap, LOPO) is complete. The H1 real-scorer runs
(detached unit `aivc-paper-t17-real`, 6 arms, serial) were still running at the time of writing. Re-run `verdict` once they finish.

## Assumptions
1. **Chain position.** X = Task 14 pre-M3 lfc / its budget scale; each arm acts on the weak half (M2 rank ≤ median), k 8, η 0.75, then budget match to the same Bt.
2. **Offline scorer.** The ERA-pack evaluator (= `sch_sumeval`). It is exact for budget-invariant arms; every arm is budget-matched. Checked on H1 against the real scorer (tolerance 0.001), pending.
3. **Alternatives:**
   - *M3_other:* the builder recipe on each other line's controls, mapped to the destination axis and re-orthonormalised (98.7–96.8 % of PC energy survives the mapping on H1). Primary = mean over the 6 lines.
   - *M3_perm:* 20 gene-wise permutations of the destination's 6,000-cell control matrix.
   - *predPC:* top-8 right singular vectors of X (override entries zeroed, centred).
   - *Systema:* X_w − η·mean_t X, η ∈ {0.25, 0.5, 0.75, 1}, best per panel on the scored data (favours the alternative).
4. **Bootstrap.** Per-target PDS; paired target bootstrap within panel, B 2,000, panels weighted equally. The CI conditions on the fixed truth cells.
5. **Real-scorer emission.** The pinned sums round(fm·e^L·n) are spread over 300 random controls (share emission). Only PDS is meaningful; DE metrics and the composite are not comparable with F2.
6. **H1 real-scorer matrices.** M3_perm = permutation 0 (pre-specified). M3_other = the line with the best H1 PDS (KOLF PCs).

## Checks
| check | expected | got | pass |
|---|---|---|---|
| none / M3_dest PDS vs Task 14 | \|Δ\| < 2e-4 | 0.0 on 7/7 | yes |
| recomputed dest control PCs vs frozen `ctrl_V` | abs dot > 0.99 | ≥ 0.99999 | yes |
| frozen V gram deviation | small | 1.9e-4 (H1) | yes |
| budget invariance on H1 (M3_dest, Systema-best) | \|PDS_real − PDS_off\| ≤ 0.001 | pending | — |

## Result (offline, 7 panels)
| M3 − alternative | mean ΔPDS | 95 % CI | panels > 0 |
|---|---|---|---|
| none | +0.0131 | [+0.0101, +0.0160] | 7/7 |
| M3_other (mean of 6) | +0.0058 | [+0.0039, +0.0076] | 7/7 |
| M3_perm (mean of 20) | +0.0131 | [+0.0101, +0.0159] | 7/7 |
| predPC | +0.0065 | [+0.0036, +0.0095] | 7/7 |
| Systema-best | +0.0145 | [+0.0108, +0.0179] | 6/7 (H1 −0.0076, η = 1) |

- **All five means are > 0 with CIs excluding 0, so H3b is not rejected.**
- **Systema on H1.** On H1, Systema with η = 1 beats M3 offline (0.9180 vs 0.9104). If the real scorer confirms this, the "same sign on H1 against all five" condition fails, and the verdict is **inconclusive**.
- **Single other lines.** Some individual other-line PCs beat M3 on 3 panels (H1 −0.0013 with KOLF PCs, OR_HCT116X −0.0005, RN_Jurkat −0.0021).
- **LOPO.** The (k, η) choice is k 24, η 1 on 6/7 folds; it is k 8, η 1 for OR_HCT116X.
  - Mean held-out gain over none: +0.0147, vs +0.0131 for the frozen k 8, η 0.75. The held-out gain is positive on 6/7 panels; on H1 it is −0.0038, i.e. −0.0050 vs frozen.

## Evidence base
- `results/rq3_h3b.json` and `results/rq3_corr/offline_summary.json`.
- Logs: `results/logs/paper-t17-{offline,real}.log`.

## What stays open
- The H1 real-scorer signs and the budget-invariance check (pending).
- The selection effect in Systema-best is not reflected in its CI.
- The CIs ignore resampling of the truth cells.
