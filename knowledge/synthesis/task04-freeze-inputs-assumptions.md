---
title: "Task 4 — frozen inputs, source registry and leave-destination-out audit: assumption checks"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# Task 4 — frozen inputs, source registry, leave-destination-out audit

Scripts: `output/data-analysis/task04_freeze_inputs.py`, `output/data-analysis/common.py`,
`output/code/build_inputs_manifest.py`. Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` Task 4, operational definition 5.
Results: `output/data-analysis/results/{panels,families,task04_checks,g_cap}.json`, `official_record.csv`; manifest
`output/code/INPUTS.md`, `inputs.sha256`.

## Assumptions and checks

| # | assumption | check | result |
|---|---|---|---|
| S1 | ERA-pack truth = `final_harness_2026-09-24/panels.json` paths | pack target set == truth target set, 7/7 panels | pass |
| S2 | "excluded sources" = union over M1 first part, fusion, neighbour library, S4; parsed by `ast` from `fh_harness.EXCL`, `ml_h1fuse.OUTER/S4KEY`, `sch_emit.S4_LINE`, `ext_baselines.CONFLICT_GENE` | fusion list in packs == FAMS minus OUTER, 7/7 | pass |
| S3 | leak = destination's own line (incl. parental line / same donor) or the truth's own dataset not excluded | per-panel table in `panels.json` | **KOLF: LEAK** (FENG pools HipSci kolf_2 = KOLF2.1J parental line + kolf_3; 11.0–17.6 % of FENG cells); 6/7 OK |
| S4 | C_public = distinct lines with ≥ essential-scale CRISPRi/KO Perturb-seq among active F2 families + H1GW | `families.json` | 11 (sensitivity 9–11) |
| S5 | submission count = distinct ids + history rows without id | `task04_checks.json` | 36 validation (= leaderboard `submission_count` 36 on 2026-10-27) + 2 final = 38 |
| A3 | G∩ on the official axis via `vcclib.GeneMap.canon`, ≥ 1 control UMI per panel | `g_cap.json` | 6,014 genes (bounded by the RN panels' 7.7–9.0 k gene axes) |

## Caveats
- The `ext_baselines.CONFLICT_GENE['KOLF']` labels (`FENGLKOLF2/3`) are v1 names absent from `table_v2`; there the per-line
  kolf_2/kolf_3 Feng contexts are filed under biology KOLF (so S4 does exclude them), but pooled Feng rows (biology FENG) keep kolf cells.
- `rho_hat` (M2 in the packs) is "leave-this-line-out" per the ERA task docstring; not re-derived here.
- RN panel truth is a subsample (800 controls, ≤ 200 cells/target); RN_K562 comes from GWPS, not K562 essential (the plan's path).
- Tahoe cell-level data and the original `tahoe/table.npz` are not on disk (deleted 2026-10-11) → affects Task 8.
- Lab/scale annotations in `families.json` combine the source JSON metadata with hand annotation; citations are resolved in Task 2.
