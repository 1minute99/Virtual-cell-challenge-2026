# INPUTS — frozen input manifest for plan `zero-shot-limits-context-complexity` (Task 4, 2026-11-05)

- **Hashes:** `output/code/inputs.sha256` (440 files, 54.68 GB; paths relative to `/home/wonmin-kim/AIVC`).
  Per-file group, size and hash: `output/code/inputs_filelist.tsv`. Builder: `output/code/build_inputs_manifest.py`
  (import closure: `output/code/import_closure.py`, static AST over `vcc2026_data/scripts` + `scripts/*.py|sh` references).
- **Backup:** `/media/wonmin-kim/Wonmin/aivc_backup_2026-10-27/paper_inputs/<same relative path>` (rsync -a, ionice -c3),
  verified with `cd .../paper_inputs && sha256sum -c /home/wonmin-kim/AIVC/output/code/inputs.sha256`.
- Verify locally: `cd /home/wonmin-kim/AIVC && sha256sum -c output/code/inputs.sha256`.

## Frozen pipeline (S18 / F2)

| item | value |
|---|---|
| `vcc2026_data/scripts/tv2_final.py` sha256 | `113b3902810362e5bd86e496ff9ca95a0885e929a814dfae80caec4d843b6842` |
| `tv2_submit_build.py` | `e5165893ba830ffa865dfa5637b53333634fe0502efae56cef10886a37efa0d2` |
| `tv2_model.py` (source registry `FAMILY`) | `4c673dfc789de28bae4e76944758e75a41993477331ace035418706a208e8e3b` |
| `sch_sumeval.py` | `ddcb39fd1d8fd9b8bbdb52bcf5e57ab5896465cd25bd94636b7d13dd5be632fd` |
| `ml_ctx_curve.py` / `ml_h1gw_value.py` | `a381fcdd…494d` / `8933b181…b036` |
| F2 spec `reports/final_2026-10-27/F1/F2_spec.json` | `895a3d86650cd73a93ef29b77caf79dcef3fdef026152da3daf924d8104beebc` |
| F1 spec `…/F1/F1_spec.json` | `3a45b1272b03e272d77c261948ba345b6c7775926d188cd0d0c76199102a8740` |
| F2 package `reports/submission_2026-09-23/F2.vcc` | sha256 `0c2a01803ee94da365cd8a39f74120b6cd64fe3bde66d212fe8c069e8a8046c8` (entry oYO6WAS8x4FgQt9ypgmq, last final entry) |

**F2 build (exact).** Spec step (`reports/final_2026-10-27/run_F1_spec.sh`, from `vcc2026_data`, env `vcc2026`):

```bash
export TV2_BASE_TABLE=/home/wonmin-kim/AIVC/vcc2026_data/reports/final_2026-10-27/base_table_stub
~/miniconda3/envs/vcc2026/bin/python -P scripts/tv2_final.py --controls-dir data/final_2026-10-22/controls --contexts D,E,F \
  --out reports/final_2026-10-27/F1 --tag F1 --rn-sources --fam-weights DLD1=1 --s4 --blend 0.3 \
  --dual 1.0,3000,0.5,ag,cpc,ci5,lsb --kink 1.17,0.8 --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g --sum-policy auto \
  --sum-ens 5,10,15,30,50,100,200 --de-follow-sum --sum-resid 0.25 --s4-excl KOLF,FENG --de-denoise intended \
  --sum-tamp 0.5 --same-line 0.5 --sum-cens 5,10,20,100 --de-drc 0.8,1.1 --sum-defl 8,0.75,0.5
```

F2_spec.json = F1_spec.json with exactly one change, `sum.ew_adaptive: true -> false` (diffed 2026-11-05). Build:
`bash scripts/tv2_submit_chain.sh reports/final_2026-10-27/F1/F2_spec.json F2` (`run_F2_build.sh`).
**M3** = `--sum-defl 8,0.75,0.5` (spec `sum.ctrl_defl` k 8, eta 0.75, weaker 50 %). **M2 in F2** = `--sum-tamp 0.5`
(cross-family agreement rank^0.5, `tv2_tamp_prep.py`).

**S18 (validation, entry etQ6SGWT1uu0wk24xuvQ, 0.2545, rank 74 at upload / 129 at close).** No spec file on disk; the full
spec is embedded in `reports/submission_2026-09-23/S18.pred_build.json["spec"]`: S17 sum side + `ctrl_defl {k 8, eta 0.75,
weak 0.5, beta 0.75, ncells 6000, seed 0}`, `tamp {tlrn_official.npz, beta 0.75}` (learned S16 amplitude), fusion weights
K562 8 / HCT116 1 / HEK293T 0.5 / CD4 1 / FENG 0.25 / H1 2 / HEPG2 1 / RNJURKAT 1 / RPE1 1 / DLD1 1, `ew_adaptive true`.
Uploaded by `reports/submission_2026-10-11/upload_s17d_s18.sh` ("S17 + control-cell PC deflation of the weaker half (k8, eta 0.75)").
**Note for Task 14/18:** the plan's "M2 = sum.tamp beta 0.5" is F2's M2; S18's M2 is the learned `tlrn` amplitude at beta 0.75,
and the ERA packs' S16 reference (`aux_<P>.npz rho_hat`) also uses beta 0.75. Use F2 settings as "frozen" unless the plan says S18.

## Groups (file lists in `inputs_filelist.tsv`)

| group (plan tasks) | files | GB |
|---|---|---|
| T6–T8 Tahoe context curve: `data/processed/drug_context_pretrain/*` (tahoe.npz sha `1976c53f…525f`), `context_splits.json` (`7d20381e…6e19`), `scaling_table.json` (`3330549d…1ec5`), FINDINGS, `dcp_{data,model,prepare,select}.py` + closure | 14 | 1.04 |
| T9, T12–T18 harness truth (cell level), harness `panels.json` / `neighbours.json`, real-moments caches, RN mock panel300, raw scPerturb parents (K562 GWPS + essential, RPE1, HepG2, Jurkat), KOLF/Jurkat harness files, Orion real_thin | 44 | 25.23 |
| T12–T17 ERA packs `reports/era_2026-10-07/packs/*` (8 panels × inputs/score/ctrlpc/ctrlvar/aux/post/sunny) + `sch_sumeval.py` closure + ERA `evaluator.py`/`prep.py`/`task.py` | 73 | 11.15 |
| T9, T11, T14, T18 M1 source registry: `transfer_v2_2026-09-23/sources/*`, `sources_new/*` (incl. H1GW.npz `417c861c…94d`), harness `final_harness_2026-09-24/sources/*`, pack-builder inputs `sc_hybrid_2026-09-25/h1fuse/*_data.npz`, S4 parts per panel, Feng cell metadata, `ml_h1fuse` / `fh_harness` / `sch_emit` / `s4i_refit_*` closure | 136 | 17.25 |
| T10 arch screen: `arch_zoo.py` closure, `arch_screen_2026-09-22{,_HepG2}/*.json`, FINDINGS | 45 | <0.01 |
| T11, T18 `ml_ctx_curve.py`, `ml_h1gw_value.py`, H1GW, `ctx_curve.json`, `h1gw_value.json` | 14 | 0.24 |
| T14, T17, T18 frozen pipeline: `tv2_final.py` + `tv2_submit_build.py` + `tv2_submit_chain.sh` + `tv2_tlrn_prep.py` closure (49 scripts), F1/F2 specs, README, run scripts, F2/S18 receipts | 67 | <0.01 |
| T13, T19, T20, T22+ reports: arch_screen / ext_contexts / sc_hybrid FINDINGS, PROJECT_HISTORY.md, `drug_gene_concordance_2026-09-16/*` (RESULTS.md) | 20 | <0.01 |
| T19 official record: `official_scores.csv`, both `submission/official_scores_new.csv`, S1/S2 status logs, status/receipt JSONs, leaderboard snapshots | 79 | 0.02 |

(A file can belong to several groups; the sha file lists each once.)

## Missing / deviating inputs (logged in `knowledge/_meta/log.md`)

1. **`reports/simple_residual_2026-09-16/data/tahoe/table.npz` — MISSING** (deleted by `reports/cleanup_2026-10-11`, line 1287 of
   its `delete_list.txt`; not on the external drive). Only the processed condition-level `tahoe.npz` (13,011 conditions × 18,150
   genes, 48 lines, `is_repeat`, `plate`, `weight`) survives. **No cell-level Tahoe data is on disk**, so plan Task 8's cell-level
   split halves for σ²_ε / Var(i) cannot be computed as written (raw Tahoe would have to be re-downloaded, or replicate
   conditions `is_repeat` used instead — a method change that needs a logged deviation).
2. **`reports/transfer_v2_2026-09-23/cache_real_H1/pseudobulk_moments_*.moments.npz` — MISSING**, so `sch_sumeval.real_moments("H1")`
   fails. The same moments are inside `packs/score_H1.npz` (`mr`, `cr`, `m0`, `c0`); use the ERA evaluator, or regenerate the
   cache with one real-scorer run.
3. Plan truth paths that differ from what the ERA packs were built from (common.py uses the harness paths; `which="full"` gives
   the parents): KOLF `checkpoints/state/kolf_real.h5ad` and RN Jurkat `checkpoints/state/jurkat_real.h5ad` do not exist;
   **RN_K562 truth is subsampled from Replogle K562 GWPS, not K562 essential**; RN panels are 800-control / ≤200-cell subsamples
   (`simple_residual_2026-09-16/mock/<line>/panel300/real.h5ad`); OR_HCT116X truth is `real_panel.h5ad` (263 targets), not `real_thin.h5ad`.
4. Not backed up (recomputable or too large, not read directly by Tasks 6–21): `final_harness_2026-09-24/cache` (52 GB delta caches),
   the S4 training tables `ext_contexts_2026-09-21/table_v2` (31 GB dir), `final_2026-10-27/F1/{de_dnK562,knn*,s4}` (2.8 GB build intermediates).
