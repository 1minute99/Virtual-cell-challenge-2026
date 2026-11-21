---
title: "RQ2 Task 10 — architecture seed replication: assumption checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-06
status: draft
author: llm
---

# RQ2 Task 10: seed replication of the architecture screen, assumption checks

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 10 (a descriptive consistency check for H2c; no
verdict attached). Seed-0 source: `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md` (and `_HepG2/`).
Scripts: `output/data-analysis/rq2_arch_replicate.sh` (queue), `rq2_arch_seed_screen.py` (seed wrapper around the
unchanged `vcc2026_data/scripts/arch_screen.py`), `rq2_arch_preset.py` (wrapper around the unchanged
`ext_preset_score.py`), `rq2_arch_seeds.py` (collector). Outputs: `output/data-analysis/results/rq2_arch_seeds.csv`,
`rq2_arch_seeds_summary.json`, `rq2_arch_seeds.png`, run folders `results/rq2_arch/s{1,2}_{K562,HepG2}/`.

**Status 2026-11-06: BLOCKED, no seed-1/2 result yet.** Run 1 failed on the deleted base-table `y.npy`, now fixed with the stub (A9). Run 2 then stalled in the seed-0 refit inside the training scope: anon 9.1 GB is over the 8 GB cap (MemoryHigh 8G), so the job was throttled to 6 s of CPU in 4 min, and it was stopped. Raising the cap, or patching the refit in a bit-identical way, is a main-session decision. See the log, 2026-11-06.

~~Status at writing: QUEUED, no seed-1/2 result yet.~~ The unit `aivc-paper-t10-queue` waits for
`aivc-paper-t17-real` and `aivc-paper-t18-queue` to exit. The tables below fill in from the summary JSON.

## Selection (recorded before any run: `results/rq2_arch/selection.json`)

Top 3 by seed-0 two-line-mean composite: additive shift (+0.0036), mixture of experts (+0.0034), bilinear (+0.0010).
Added because the plan names them: full hypernetwork (141,888 params; −0.0011), token Transformer (86,496; −0.0034),
cross-attention (−0.0254). The LoRA hypernetwork (seed-0 single-line best, HepG2 +0.0169) is **not** among the six.

## Assumptions

| # | Assumption | How it is checked |
|---|---|---|
| A1 | Composite = `avg_score` (from_replicate) of `cell-eval2 score` against the frozen real bundle, as FINDINGS ranks | collector reproduces FINDINGS seed-0 deltas (additive K562 +0.0154, HepG2 −0.0083, mean +0.0036): **pass** |
| A2 | "Same protocol" = the seed-0 code with only the seed changed. The seed sets training init, minibatch RNG and the cell-reconstruction seed together, as in seed 0 | wrapper changes only `arch_screen.SEED` and the output folder |
| A3 | Paired baseline: refit_v2 is re-scored at reconstruction seeds 1 and 2 (4 extra runs, 28 in total) so each Δ is paired | in the CSV, `base_composite` is per (line, seed) |
| A4 | Seeds are independent replicates; the lines are not exchangeable (K562 has no context specificity) | design; the two-line mean is primary, as in FINDINGS / Task 18 |
| A5 | 95 % CI = Student t, n = 3, df = 2; normality is assumed | **untestable at n = 3**, so a sign count is reported as well. The exact sign test cannot reach p < 0.05 with n = 3 |
| A6 | Protocol reproduction: the K562 real panel `b19_sweep2_2026-09-18/full_metrics/K562/real_disjoint.h5ad` is **missing on disk**. It is regenerated with the same code into `results/rq2_arch/real/K562/` | SHA vs b19 receipt (7c7a777b…). HepG2 is regenerated as a determinism control (4d37d551…). **Gate:** the refit_v2 seed-0 K562 composite re-scored on the regenerated panel must equal 0.3480609 within 1e-6, or the queue stops |
| A7 | Current torch environment (ct_transfer) ≈ the September environment | additive + moe are retrained at seed 0 on K562 and compared with the stored harness PDS, best step and composite. Reported, not a gate |
| A8 | No multiplicity correction (descriptive check) | — |
| A9 | (added 2026-11-06, execution fix) The base table `simple_residual_2026-09-16/data/table/{y,controls}.npy` was deleted by the 2026-10-11 cleanup; the first queue run failed 3 s into seed-0 training (`FileNotFoundError …/table/y.npy`). The documented NaN stub `reports/final_2026-10-27/base_table_stub` is used via `TV2_BASE_TABLE` | code read: `arch_screen.py` takes from the base table only `features` (shape), `feature_map`, `targets`; the training and refit rows come from `ext_contexts_2026-09-21/table_v1` (intact). Load smoke test under the 8 G scope passed. A NaN leak would show up in the A7 seed-0 additive/moe comparison |

## Checks (filled from `rq2_arch_seeds_summary.json` → `seed0_reproduction`)

| Check | Expected | Got | Pass |
|---|---|---|---|
| Seed-0 deltas reproduced from the stored scored.csv | FINDINGS table | identical (6 archs × 2 lines) | pass |
| Regenerated real panels byte-identical | SHA equal | HepG2 4d37d551…, K562 7c7a777b… = reference | pass |
| refit_v2 seed-0 composite on the regenerated K562 panel | 0.3480609 ± 1e-6 | pending | — |
| additive/moe seed-0 retrain matches stored | harness PDS 0.8931 / 0.8740 | pending | — |
| PDS parity, preset vs harness (built into ext_preset_score) | ≤ 1e-9 | pending | — |
| Memory: training tree ≤ 8 GB (MemoryHigh 8G / MemoryMax 10G scope; `time -v` max RSS in screen.log) | ≤ 8 GB | anon 9.1 GB in the seed-0 K562 refit (`ext_film.fold`), throttled + 1 GB swap full | **fail** (run stopped) |
| Memory under the main-session 12G exception (training step only: MemoryHigh 12G / MemoryMax 14G) | ≤ 12 GB | anon 13.3 GB + 1 GB swap, throttled; estimated peak ≈ 15–16 GB | **fail** (run stopped 2026-11-06) |

## Evidence base

FINDINGS.md seed-0 table (17 architectures, 104 runs). [[rq4_h4-assumptions]] uses the best-architecture gain
(seed 0: additive +0.0036 primary; LoRA HepG2 +0.0169 sensitivity). Task 10 replaces the primary value with the
seed mean over the six replicated architectures.

## What stays open

- n = 3 seeds gives very wide t intervals (×4.30 SE), so "CI excludes 0" is a high bar. A null here is weak evidence.
- LoRA (seed-0 single-line best) is not replicated. Task 18's sensitivity value therefore mixes a seed-0-only figure
  with seed means.
- Scorer setup differs from Task 18 (K562/HepG2 panels, ext_preset emission, not H1 / F2).
