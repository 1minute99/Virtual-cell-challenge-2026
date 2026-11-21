---
title: "zero-shot-transfer-destination-state — research plan"
type: plan
created: 2026-10-31
updated: 2026-10-31
status: superseded   # design replaced by [[zero-shot-limits-context-complexity-design]] (2026-11-05); never pre-registered
methodology: quantitative
author: mixed
design: "[[zero-shot-transfer-destination-state-design]]"
---

# 처음 보는 세포주로의 교란 반응 전이: 오차는 목적지 control 상태 공간에 있는가 — 연구 계획

> **For agentic workers:** REQUIRED SUB-SKILL: Use `research-superpowers:executing-research-plan` to run this plan task by
> task. Tasks use checkbox (`- [ ]`) syntax for TodoWrite integration. 경로는 모두 프로젝트 루트 `~/AIVC/` 기준.
> 오래 걸리는 계산은 반드시 `scripts/run_detached.sh NAME LOG_DIR /abs/program args...`로, `cell-eval2`는 **한 번에 하나만**
> (회당 약 25분, 약 15 GB). `vcc2026_data/reports/simple_residual_2026-09-16/data/table`, `disc_residual/fam_*`,
> `denoise_src/full`은 **절대 지우지 말 것** (AGENTS.md).

**Research Question:** 공개 Perturb-seq 스크린을 다른 세포주에서 옮겨 와 처음 보는 세포주의 CRISPRi 반응을 예측할 때,
전이가 놓친 오차(정답 − 예측)는 측정 잡음을 통제한 뒤에도 목적지 세포 **자신의** control 세포 상위 주성분 공간에
집중되는가? 그리고 그 공간을 이용한 보정(M3)은 대안 보정보다 판별력(PDS)을 더 높이는가?

**Hypotheses (pre-registered, immutable after sign-off):**

- **H1 (primary, E1).** *(a)* On the H1 panel, after split-half control of measurement noise, the cross-line transfer
  residual (truth − pre-M3 prediction) places a larger fraction of its energy in the top-8 control PCs of the destination
  than in the top-8 control PCs of the other harness lines (mean over the six) and than in the top-8 PCs of gene-wise
  permuted destination controls; *(b)* averaged over the seven harness panels, deflating the weaker half of target
  predictions along the destination's control PCs (M3) raises PDS more than the same deflation along other lines'
  control PCs, along permuted-control PCs, along the top PCs of the prediction matrix itself, and more than a
  Systema-style common-response removal.
- **H2 (secondary, E2).** The PDS gain from source contexts has not saturated at the full set: on the H1 panel the
  leave-one-family-out loss at k = K is positive.
- **H3 (secondary, E3).** On the H1 panel, adding one other-lab same-line screen (H1GW, w = 1) is worth more composite
  score than the largest single modelling gain (M2, M3, or the best of the 17 architecture comparisons), measured on
  the same base, panel, and scorer.
- **H4 (secondary, E5).** Applying M2 + M3 post hoc to the zero-shot pseudobulk deltas of a deep model never trained on
  the destination line (PIE `replogle_xdataset`) raises its PDS, averaged over the honest panels (RN ×4 + KOLF).

**Falsification Criteria:**

- **H1a rejected** if, on H1 at k = 8, the target-bootstrap 95 % CI (B = 2,000, averaged over 5 random splits) of
  either D_other = F(U_dest) − mean_j F(U_j) or D_perm = F(U_dest) − F(U_perm) includes 0 or lies below 0 (both must
  exclude 0 above zero — intersection-union, no multiplicity correction needed). Also rejected if the control-only
  negative control (Task 8) shows D ≥ 50 % of the observed D (i.e. the estimator's structural bias explains the effect).
- **H1b rejected** if the 7-panel mean ΔPDS (M3_dest − alternative) is ≤ 0 against **any** of the five alternatives.
  **Confirmed** if it is > 0 against all five with a 95 % CI (target bootstrap within panel, panels equally weighted)
  excluding 0, and the H1 real-scorer run has the same sign against all five. Anything between → "inconclusive".
- **H2 rejected** if the mean leave-one-family-out ΔPDS on H1 at k = K has a 95 % CI that includes 0 or lies below it, or
  if a saturating model (a − c·e^(−k/τ)) beats the log-linear model (a + b·ln k) by ΔAIC ≥ 2.
- **H3 rejected** if composite(base + H1GW, w = 1) − composite(base) ≤ the largest single modelling gain, both measured
  with the real `cell-eval2 vcc2026` scorer on the H1 panel. (The ≥ 10× ratio from the design is reported
  descriptively, not tested.)
- **H4 rejected** if the mean ΔPDS (PIE + M2/M3 − PIE) over the five honest panels is ≤ 0. Confirmed if > 0 with a 95 % CI
  excluding 0. Otherwise → "inconclusive".

**Stop Criterion:** Analysis freeze on **2026-11-14**. Data collection is complete when: E1 covers all 7 panels × 5
splits × k ∈ {8, 24, 50}; E2 has R = 20 random subsets per k (or every combination if there are fewer than 20; minimum
R = 10 if it runs late — the only reduction allowed); E3 has w ∈ {0, 0.25, 0.5, 0.75, 1}; E5 covers the 5 honest panels.
No panel, k, or split is added after the freeze. Any further data counts as exploratory and goes in the Deviations Log.

**Method:** A frozen interpretable fusion pipeline (M1 multi-source transfer + M2 agreement amplitude + M3
destination-control deflation, settings of S18/F1) predicts each harness panel leave-destination-out. E1 measures residual
energy on candidate subspaces with a split-half cross-product estimator that removes measurement noise. Its controls use
a three-way split: PCs, truth-A baseline, and truth-B baseline each get their own control cells. E1 then compares M3
against four alternative corrections plus no correction. E2–E5 measure, on the same metric, the scaling of source
contexts, the value of same-line data, the official validation record, and the transfer of M2/M3 to a deep model.

**Output target:** article (arXiv preprint, cs.LG primary + q-bio.GN cross-list; bioRxiv fallback), 8–10 pages main
+ appendix, code on GitHub. Submission window 2026-11-28 … 12-01.
**Language:** en (manuscript); plan/wiki notes ko + en.

---

## Design → plan decisions (resolved before pre-registration; flag at sign-off)

1. **E5 model: PIE primary, STATE conditional.** The design names STATE first. But `checkpoints/state/ST-SE-Replogle` has
   **no weights** locally, uses a one-hot perturbation vocabulary (RN panels only), and was trained on the Replogle–Nadig
   lines. So it is not zero-shot on the RN panels. PIE `arcinstitute/PIE_replogle_xdataset` never saw Replogle or KOLF
   data, which makes it honest on RN ×4 + KOLF. STATE enters only if Task 15 finds a public checkpoint with the
   destination line held out.
2. **Three-way control split for noise control.** If truth halves A/B shared one control baseline, its noise
   (≈ Σ_ctrl/n, which lies on the control PCs) would appear in the cross-product and bias the test toward H1. So the
   control cells split into C_pc (PC estimation), C_a (baseline for half A), and C_b (baseline for half B).
3. **Cross-line comparison in a shared gene space.** The other lines' PCs are only defined on genes measured in their
   own panel. The primary E1 analysis therefore uses the intersection gene set G∩ shared by the 7 panels, and recomputes
   every subspace there. A destination-only analysis on the full scored gene space is secondary.
4. **Stronger alternatives than the design.** H1b also tests **prediction-PC deflation**: deflating along the top-k PCs
   of the predicted response matrix (target × gene). This is a destination-agnostic null for "is it just the dominant
   axes of the predictions?". The Systema-style arm gets its own η grid {0.25, 0.5, 0.75, 1.0} and keeps the best value
   (favours the alternative, so the test is conservative).
5. **M3 hyperparameters are not out of sample.** k = 8, η = 0.75, and the weak 50 % were picked on these harness panels
   before this plan. H1b is therefore not a held-out test for M3's settings. Secondary analysis: leave-one-panel-out
   (choose (k, η) on 6 panels, evaluate on the 7th). E4's official scores are the only fully out-of-sample evidence,
   and the paper says so.
6. **Strength strata use test-time amplitude.** Strata are quartiles of predicted amplitude, as in the 2026-10-29
   strong screen. Stratifying by truth amplitude would select on the outcome.

## Data Sources

**Literature** (provisional bibkeys; Task 2 confirms or renames them and logs every rename):

- [[adapert-2026]] AdaPert (ICML 2026; Appendix Table 13: learned models vs cross-line mean transfer) — *acquire pending*
- [[pan-2026]] Pan et al. 2026 (genome-scale H1 CRISPRi = H1GW, GSE295214) — *acquire pending*
- [[state-2025]] STATE (Arc Institute, 2025) — *acquire pending*
- [[molina-zhang-2026]] Molina & Zhang 2026 (cell-line patterns and residuals on control PCs) — *acquire pending*
- [[cipher-2026]] CIPHER (2025/2026; control covariance does not transfer across cell types) — *acquire pending*
- [[dibaeinia-2026]] Dibaeinia et al. 2026 (context diversity over model size) — *acquire pending*
- [[systema-2025]] Systema (systematic variation in perturbation-response evaluation) — *acquire pending*
- [[replogle-2022]] Replogle et al. 2022 K562 GWPS / RPE1 — *acquire pending*
- [[nadig-2025]] Nadig et al. 2025 HepG2 / Jurkat — *acquire pending*
- [[orion-2025]] X-Atlas/Orion (HCT116, HEK293T) — *acquire pending*
- [[vcc-2025]] Virtual Cell Challenge description + cell-eval — *acquire pending*
- [[pie-2026]] PIE (Arc; model card / preprint) — *acquire pending*
- [[ahlmann-eltze-2025]] Linear baselines for perturbation prediction — *acquire pending*
- Dataset citations to resolve in Task 2: KOLF2.1J iPSC CRISPRi screen, HipSci (Feng), CD4 T-cell genome-scale CRISPRi,
  DLD-1 (GSE337988), DepMap 23Q4.

**Data on disk** (all *on disk*; frozen in Task 4):

- Harness truth (cell level): `vcc2026_data/reports/transfer_v2_2026-09-23/h1/real_thin.h5ad` (H1),
  `vcc2026_data/reports/sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad` (OR_HCT116X),
  `checkpoints/state/kolf_real.h5ad` (KOLF), `checkpoints/state/jurkat_real.h5ad` (RN Jurkat),
  `vcc2026_data/raw/scperturb/ReplogleWeissman2022_K562_essential.h5ad` (RN K562),
  `vcc2026_data/raw/scperturb/NadigOConner2024_hepg2.h5ad` (RN HepG2),
  `vcc2026_data/raw/scperturb/ReplogleWeissman2022_rpe1.h5ad` (RN RPE1)
- Offline evaluator packs: `vcc2026_data/reports/era_2026-10-07/packs/{inputs,aux,post,score,ctrlpc,ctrlvar}_<panel>.npz`;
  evaluator `vcc2026_data/scripts/sch_sumeval.py` (exact only for budget-invariant arms)
- Pipeline: `vcc2026_data/scripts/tv2_final.py` (`--sum-defl 8,0.75,0.5` = M3; `sum.tamp` β 0.5 = M2),
  `vcc2026_data/scripts/harness_pred.py`, `vcc2026_data/scripts/ml_ctx_curve.py` (E2 seed),
  `vcc2026_data/scripts/ml_h1gw_value.py` (E3 seed), `sources_new/H1GW.npz`
- Official record: `vcc2026_data/reports/official_scores.csv` (29 rows) + `submission/official_scores_new.csv` — the
  design says 36 submissions; Task 4 reconciles the count. Leaderboard snapshots `vcc2026_data/reports/leaderboard_*.json`.
- Architecture comparison: `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md` (17 architectures)
- PIE: `arcinstitute/PIE_replogle_xdataset` (public HF) + PIE_sources (~44 GB, external drive
  `/media/wonmin-kim/Wonmin/`) — *acquire pending*

**Reference values to reproduce** (the paper reports only the frozen-pipeline values; disagreements go in the Deviations
Log): H1 PDS ≈ 0.796 + 0.048 ln k; H1GW +0.074 composite (w = 1); strong-target truth control-PC energy H1 PC24 0.188 vs
prediction 0.070; M3 positive on 7/7 harness panels; official M2 +0.0106, M3 bundle +0.0070.

## Tasks

### Task 0: Scaffold the research wiki inside AIVC (week 1, 11-01)
**Files:** `CLAUDE.md`, `knowledge/{sources,entities,concepts,synthesis,_meta}/`, `knowledge/_meta/log.md`,
`knowledge/_meta/gate-overrides.log`, `output/{data-analysis,article,bibtex,code}/`, `scripts/lint-wiki.py`

- [ ] Copy `knowledge/`, `output/`, `schema/`, `scripts/lint-wiki.py`, and `scripts/vendor/` from the plugin template
      `research-project-template` without touching existing `scripts/*`. Use `cp -n` (never overwrite).
- [ ] Write `CLAUDE.md` with frontmatter `methodology: quantitative`, `discipline: "ML for computational biology"`,
      `languages: [en, ko]`. Its body points to `AGENTS.md` for the compute rules.
- [ ] Add a line to `knowledge/_meta/log.md`: "2026-11-01 plan pre-registered: [[zero-shot-transfer-destination-state-plan]]"
- [ ] Check `quarto --version` and `xelatex --version`. If either is missing, record it in the log and plan Task 23's
      render with plain LaTeX (`output/article/main.tex`).
- [ ] `python scripts/lint-wiki.py` exits 0 on the empty wiki.

### Task 1: External gates — competition rules and arXiv endorsement (week 1)
**Files:** `knowledge/_meta/log.md`, `input/notes/external-gates.md`

- [ ] Read the VCC 2026 terms on publishing results and data use. Record the exact clause (quote + URL) in
      `input/notes/external-gates.md`, and note whether it permits publishing before the final results.
- [ ] Ask for arXiv cs.LG endorsement: find the endorsement code under the user's arXiv account and list 2–3 candidate
      endorsers (authors of A-graded cs.LG papers from Task 2). The user sends the request; this task only drafts the
      email text.
- [ ] If there is no endorsement by 2026-11-21, switch to bioRxiv and log the decision.

### Task 2: Literature search (week 1)
**Files:** `input/bibliography/literaturguide.md`, `output/bibtex/references.bib`

- [ ] Run `literature-review` on: zero-shot / unseen-cell-line perturbation prediction; mean-transfer baselines; control-
      state / control-PC structure of responses; Systema-type systematic-variation corrections; context-count scaling.
      Seed set = the Data Sources literature above.
- [ ] Grade A/B/C. A = the 13 named works plus any work that already runs the E1-type comparison (residual on own vs other
      control axes). **If such a work exists, stop and raise it with the user before Task 7** (risk to novelty).
- [ ] Fix the bibkeys in `references.bib` and log any rename of the provisional keys above.

### Task 3: Acquire source PDFs (week 1)
**Files:** `<library>/pdf/*.pdf`, `input/bibliography/acquisition-todo.md`

- [ ] Run `acquire-sources` on the A+B set from `literaturguide.md`.
- [ ] The user manually downloads paywalled originals (ICML proceedings are open; check CIPHER/Molina & Zhang).
- [ ] Re-run `acquire-sources` until all A sources are on disk or explicitly deferred.

### Task 4: Freeze inputs and the pipeline (week 1)
**Files:** `output/code/INPUTS.md`, `output/code/inputs.sha256`, `output/data-analysis/common.py`,
`output/data-analysis/results/panels.json`, `output/data-analysis/results/official_record.csv`

- [ ] List every input file used by Tasks 6–17 with its sha256 → `inputs.sha256`. Back up the set to
      `/media/wonmin-kim/Wonmin/aivc_backup_2026-10-27/paper_inputs/` (prevents a repeat of the 10/11 deletion).
- [ ] `common.py`: load panel truth (cell level), controls, packs, and G∩ (gene intersection over the 7 panels; write
      |G∩| and per-panel gene counts to `panels.json`).
- [ ] Confirm leave-destination-out per panel (no K562 source when RN K562 is the destination, etc.). Write each panel's
      excluded source list to `panels.json`. **If any leak is found, stop and record it in the Deviations Log.**
- [ ] Reconcile official submission count (29 vs 36): merge both CSVs → `official_record.csv`; record the true count.
- [ ] Freeze the pipeline: record `tv2_final.py` sha + the exact S18 arm string in `INPUTS.md`.

### Task 5: Ingest A-graded sources (week 1–2, parallel to Tasks 6–9)
**Files:** `knowledge/sources/*.md`, `knowledge/entities/*.md`, `output/bibtex/references.bib`

- [ ] Ingest via `ingest-source`, focus "what is shown about zero-shot transfer error, control-state structure, and
      corrections": adapert-2026, molina-zhang-2026, cipher-2026, dibaeinia-2026, systema-2025, pan-2026, state-2025.
- [ ] Ingest (focus "data/benchmark facts to cite"): replogle-2022, nadig-2025, orion-2025, vcc-2025, pie-2026,
      ahlmann-eltze-2025.
- [ ] For each of the first seven, record in its source page exactly what it tests that overlaps with H1–H4.

### Task 6: Build frozen pre-M3 and M3 predictions for all 7 panels (week 1)
**Files:** `output/data-analysis/e1_build_preds.py`, `output/data-analysis/results/preds/<panel>_{prem3,m3}.npz`

- [ ] Script: for each panel, run the frozen S18 recipe with `--sum-defl` off (pre-M3) and on (M3), leave-destination-out.
      Save the pseudobulk predicted deltas (target × gene, log1p CP10K) and M2 confidence ranks.
- [ ] Sanity check: M3 − pre-M3 PDS from `sch_sumeval` is positive on 7/7 panels (reproduces the reference value).
      Record agreement or disagreement.
- [ ] Check that the pre-M3 prediction depends on destination controls only through the baseline profile. If it depends
      on them more, rebuild with C_a ∪ C_b as the control input, so the PCs (from C_pc) stay independent.

### Task 7: Split-half residual-energy estimator (week 1)
**Files:** `output/data-analysis/e1_energy.py`, `output/data-analysis/tests/test_e1_energy.py`

- [ ] Implement: control three-way split (C_pc / C_a / C_b, equal sizes, seeds 0–4); per-target truth halves A/B;
      Δtruth_A = pb(A) − pb(C_a), Δtruth_B = pb(B) − pb(C_b); r_A, r_B = Δtruth − Δpred;
      F(U) = Σ_t ⟨P_U r_A,t, P_U r_B,t⟩ / Σ_t ⟨r_A,t, r_B,t⟩.
- [ ] Subspaces on G∩: U_dest (top-k PCs of C_pc, log1p CP10K, centred); U_j (top-k control PCs of each other panel j,
      all its control cells); U_perm (top-k PCs of C_pc after independent per-gene permutation across cells; 20
      permutations, averaged); U_haar (uniform random, reported only).
- [ ] Unit tests on synthetic data: (i) pure noise with residual ≡ 0 → F ≈ 0 numerator in expectation; (ii) a residual
      planted on U_dest → F(U_dest) ≫ F(U_perm); (iii) shared-baseline variant shows the bias that the three-way split
      removes.
- [ ] Target bootstrap (B = 2,000) of D_other and D_perm, averaged over the 5 splits.

### Task 8: Run E1a — residual concentration (week 1)
**Files:** `output/data-analysis/results/e1_energy.parquet`, `output/data-analysis/results/e1_negctrl.parquet`

- [ ] Run all 7 panels × k ∈ {8, 24, 50} × 5 splits (detached job).
- [ ] Negative control: within the destination controls, draw 300 pseudo-targets with the same cell counts as the real
      targets, set Δpred = 0, and compute D_other and D_perm. This measures the estimator's leftover structural bias.
- [ ] Strata: quartiles of predicted amplitude (pre-M3 ‖Δpred‖). Report D per quartile.
- [ ] Secondary: destination-only analysis on the full scored gene space (no U_j).

### Task 9: Run E1b — M3 against alternative corrections (week 1–2)
**Files:** `output/data-analysis/e1_corrections.py`, `output/data-analysis/results/e1_corrections.parquet`

- [ ] Arms applied to the same weaker-half targets (by M2 rank), k = 8, η = 0.75: M3_dest, M3_other (each of 6 lines;
      primary = mean), M3_perm, predPC (top-8 PCs of the prediction matrix), Systema-style (subtract η × mean predicted
      response across targets; η ∈ {0.25, 0.5, 0.75, 1.0}, keep the best), none.
- [ ] Budget-invariance check for each arm: compare the `sch_sumeval` output against the real scorer on H1 for M3_dest
      and Systema-style. If PDS differs by > 0.001, that arm uses the real scorer on every panel (queued one at a time).
- [ ] Per-target PDS (rank of the true target) → bootstrap within panel; panels equally weighted.
- [ ] Real `cell-eval2 vcc2026` on H1 for: none, M3_dest, the best M3_other line, M3_perm, predPC, Systema-best (6 runs
      ≈ 2.5 h, serial, detached).
- [ ] Secondary: leave-one-panel-out choice of (k ∈ {8, 24}, η ∈ {0.5, 0.75, 1.0}) for M3_dest, evaluated on the
      held-out panel.

### Task 10: E1 assumption-check note (week 2)
**Files:** `knowledge/synthesis/e1-transfer-residual.md`

- [ ] Report H1a/H1b verdicts exactly against the pre-registered criteria (confirmed / rejected / inconclusive).
- [ ] Assumptions: halves independent; PC estimation stable (compare the C_pc PCs across seeds by principal angles);
      G∩ coverage per panel; depth of shallow panels (OR 163 cells, RN 53–96 cells) and the width of their CIs.
- [ ] If H1a is rejected, write the pre-planned reading: "the apparent concentration is sampling noise; M3 acts as
      noise/hub suppression".

### Task 11: E2 — controlled scaling of source contexts (week 2)
**Files:** `output/data-analysis/e2_scaling.py`, `output/data-analysis/results/families.json`,
`output/data-analysis/results/e2_curve.parquet`, `output/data-analysis/results/e2_fit.json`

- [ ] Enumerate source families (destination and same-line family excluded per panel) → `families.json` with K per
      panel.
- [ ] Extend `ml_ctx_curve.py`: R = 20 random subsets per k = 1…K (all combinations if fewer); the k = K−1 subsets are
      exactly leave-one-family-out (LOFO). Rebuild M1 + M2 (no M3) and score PDS and MSE.
- [ ] Exactness spot-check: 3 real-scorer runs on H1 (k = 1, ⌈K/2⌉, K). If PDS deviates by > 0.002, report curves as
      relative and record it in the Deviations Log.
- [ ] Fit log-linear vs saturating models (least squares on subset means, AIC). Bootstrap the LOFO mean over subsets ×
      targets (B = 2,000). H2 verdict per the criteria; other panels secondary.

### Task 12: E3 — value of same-line data vs modelling (week 2)
**Files:** `output/data-analysis/e3_sameline.py`, `output/data-analysis/results/e3_sameline.parquet`

- [ ] H1 panel, real scorer: base (frozen S18) with H1GW mix w ∈ {0, 0.25, 0.5, 0.75, 1} (5 runs).
- [ ] Same scorer and base: base − M2, base − M3 (2 runs). Best architecture gain = max over the 17 arch-screen results
      relative to its transfer baseline (from `arch_screen_2026-09-22/FINDINGS.md`; note if they used a different
      scorer).
- [ ] Per-context marginal value from Task 11 (LOFO, converted to composite with the same scorer on H1).
- [ ] Related-line cases (cut first if late): KOLF ← H1GW, OR_HCT116 ← DLD-1, w ∈ {0, 0.5, 1}.
- [ ] H3 verdict; report the ratio H1GW gain / per-context value descriptively.

### Task 13: E4 — official validation record (week 2)
**Files:** `output/data-analysis/e4_official.py`, `output/data-analysis/results/e4_official.csv`

- [ ] From `official_record.csv`: score trajectory, rank 129 / 1,635, per-component gains attributable to single-change
      submissions only (M2: S14fk2 → S15 +0.0106; M3 bundle: S17 → S18). Flag any submission that changed more than
      one thing.
- [ ] Leaderboard snapshots: aggregate only (quantiles), never other teams' names.
- [ ] Add final-phase results only if they are published before 2026-11-25.

### Task 14: Acquire PIE (week 2; needs a GPU window)
**Files:** `checkpoints/pie/`, `input/notes/e5-models.md`

- [ ] Download `arcinstitute/PIE_replogle_xdataset` and the PIE_sources (~44 GB) to the external drive. Record licence and
      version.
- [ ] Confirm from the model card which lines/datasets PIE was trained on. Fix the honest panel set (expected RN ×4 +
      KOLF) in `e5-models.md` **before** running inference.

### Task 15: STATE availability check (week 2, ≤ 30 min)
**Files:** `input/notes/e5-models.md`

- [ ] Check whether public ST checkpoints exist with the destination line held out, and with weights. If yes, add STATE
      to E5 on those panels only. If no, record "STATE excluded: no zero-shot-honest checkpoint" and stop.

### Task 16: E5 — plug-in test on a deep model (week 2)
**Files:** `output/data-analysis/e5_plugin.py`, `output/data-analysis/results/e5_plugin.parquet`

- [ ] Run PIE zero-shot inference on the honest panels (short GPU window; CPU fallback if slow) → pseudobulk deltas in
      the same space as Task 6.
- [ ] Apply M2 (agreement ranks from the source screens) and M3 (destination control PCs, same k/η/weak half) post hoc.
      Score PIE, PIE + M2, PIE + M3, PIE + M2 + M3.
- [ ] H4 verdict (on PIE + M2 + M3 vs PIE). Also report the E1a energy fractions of the PIE residual.

### Task 17: Figures (week 2)
**Files:** `output/data-analysis/figures.py`, `output/article/figures/fig{1..5}.pdf`

- [ ] Use the `dataviz` skill before plotting.
- [ ] Fig 1: pipeline + estimator schematic. Fig 2: F(U) by subspace × panel (k = 8), plus the negative control.
      Fig 3: ΔPDS of corrections. Fig 4: E2 curves + LOFO. Fig 5: value comparison (H1GW vs per-context vs modelling).

### Task 18: Synthesis pages (week 3)
**Files:** `knowledge/synthesis/related-work-positioning.md`, `knowledge/synthesis/results-summary.md`

- [ ] Positioning: for each of the 7 focus works, state what it shows and what this paper adds (claims must cite
      [[bibkey]]). Required: at least 10 cited sources.
- [ ] Results summary: H1–H4 verdicts, effect sizes with CIs, deviations. User review → `status: stable`.
- [ ] `python scripts/lint-wiki.py` exits 0.

### Task 19: Draft — Introduction + Related Work (week 3)
**Files:** `output/article/sections/01-intro.qmd`, `output/article/sections/02-related.qmd`

- [ ] Via `drafting-manuscript`: Introduction 1,000–1,400 words; Related Work 800–1,100 words. Citations `[@bibkey]`.

### Task 20: Draft — Setup and Methods (week 3)
**Files:** `output/article/sections/03-methods.qmd`

- [ ] 1,600–2,200 words: problem setup, M1/M2/M3, estimator with the three-way split, alternatives, statistics, and
      the pre-registration statement (incl. decisions 1–6 above).

### Task 21: Draft — Results (week 3)
**Files:** `output/article/sections/04-results.qmd`

- [ ] 1,800–2,500 words: E1 (main), E2, E3, E4, E5 in that order. Every number traces to a file in
      `output/data-analysis/results/`.

### Task 22: Draft — Discussion, Limitations, Appendix (week 3)
**Files:** `output/article/sections/05-discussion.qmd`, `output/article/sections/99-appendix.qmd`

- [ ] Discussion + limitations 900–1,300 words (M3 hyperparameters not held out; shallow panels; one competition;
      evaluation caveats kept short).
- [ ] Appendix 2,500–4,000 words: per-panel tables, the leave-one-panel-out analysis, estimator derivation and unit
      tests, input manifest, competition-rule note.

### Task 23: Render, code release, and review (week 4)
**Files:** `output/article/main.qmd`, `output/article/main.pdf`, `output/code/README.md`

- [ ] Render the PDF (quarto, or the LaTeX fallback from Task 0); main text ≤ 10 pages.
- [ ] Code release: `output/code/` contains Tasks 6–17 scripts, `INPUTS.md`, `inputs.sha256`, and an `environment.yml`.
      Re-run Task 8 on H1 from a clean checkout and confirm identical numbers.
- [ ] `requesting-peer-review` (constructive, then adversarial). Fix, then re-render.
- [ ] Submit to arXiv (or bioRxiv) 2026-11-28 … 12-01, with the user's go-ahead.

## Verification

The project is complete when:

- [ ] `python scripts/lint-wiki.py` exits 0
- [ ] ≥ 13 sources ingested (the A set), each with `status: stable` or `review`
- [ ] `knowledge/synthesis/e1-transfer-residual.md` and `results-summary.md` have `status: stable`
- [ ] Every result file named in Tasks 6–16 exists, and each number in `04-results.qmd` matches it (checked by a grep
      script)
- [ ] The Task 7 unit tests pass (`pytest output/data-analysis/tests`)
- [ ] `main.pdf` renders without errors; main text ≤ 10 pages
- [ ] A clean re-run of Task 8 on H1 reproduces F(U) to 4 decimals
- [ ] Peer review round completed via `requesting-peer-review`
- [ ] Each of H1a, H1b, H2, H3, H4 is explicitly marked confirmed / rejected / inconclusive per its pre-registered
      criterion
- [ ] Preprint submitted (arXiv ID or bioRxiv DOI recorded in `knowledge/_meta/log.md`)

**Fallback cuts if late (in order; each one is logged):** E3 related-line cases → STATE in E5 → E2 down to R = 10 →
PIE in E5 (H4 then reported as "not run"). H1a/H1b are never cut.

## Deviations Log

Append to `knowledge/_meta/log.md` whenever scope, method, or data sources change after pre-registration. Reference the
plan slug `zero-shot-transfer-destination-state-plan`. A result reached after a deviation is labelled **exploratory**
in the manuscript.
