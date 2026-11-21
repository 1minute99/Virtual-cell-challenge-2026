---
title: "zero-shot-limits-context-complexity — research plan"
type: plan
created: 2026-11-05
updated: 2026-11-05
status: pre-registered   # confirmed by the user 2026-11-05 ("진행"); hypotheses + falsification criteria frozen
methodology: quantitative
author: mixed
design: "[[zero-shot-limits-context-complexity-design]]"
---

# 처음 보는 세포주로의 교란 반응 예측은 데이터에 막혀 있는가: 상한과 맥락 표본 복잡도 — 연구 계획

> **For agentic workers:** REQUIRED SUB-SKILL: Use `research-superpowers:executing-research-plan` to run this plan task by
> task. Tasks use checkbox (`- [ ]`) syntax for TodoWrite integration. 경로는 모두 프로젝트 루트 `~/AIVC/` 기준.
> 오래 걸리는 계산은 `scripts/run_detached.sh NAME LOG_DIR /abs/program args...`로, `cell-eval2`는 **한 번에 하나만**
> (회당 약 25분, 약 15 GB). `vcc2026_data/reports/simple_residual_2026-09-16/data/table`, `disc_residual/fam_*`,
> `denoise_src/full`은 **절대 지우지 말 것** (AGENTS.md).

**Research Question:** 공개 스크린으로 처음 보는 세포주의 CRISPRi 반응을 예측할 때, 달성 가능한 정확도는 모델 용량이
아니라 (i) 측정 잡음 천장, (ii) 다른 세포주에서 옮겨 올 수 없는 세포주×교란 상호작용, (iii) 상호작용을 배우기에
부족한 맥락(세포주) 수에 막혀 있는가?

**Notation.** δ_{c,t} = μ + a_t + b_c + i_{c,t} + ε_{c,t} (shared response, conserved perturbation component,
cell-line component, line × perturbation interaction, measurement noise; decomposition of Molina & Zhang 2026). A
learned context model predicts i from a control descriptor x_c with effective dimension d (participation ratio across
contexts). Proposition 2 (linear-Gaussian): break-even context count C* = d + 1 + d·σ²_res / (R²·Var(i)), with the
model-free lower bound **C*_min = d + 1 + d·σ²_ε / Var(i)** (R² = 1).

**What has already been seen (honesty statement).** Seed-0 results already exist for the Tahoe context curve
(C = 6/12/24/40: residual-MSE gain −0.00046 / −0.00039 / −0.00010 / +0.00002), the 17-architecture screen, and the
learned 6→35 context curve. The hypotheses below are therefore **not** blind to those numbers. What is pre-registered:
(1) replications with new seeds and new C points, (2) a theory prediction from parameters estimated **without** the gain
curve, and (3) every analysis not yet run (ceilings, C*_min for genetic data, RQ3, RQ4). The existing seed-0 curve does
not fit Proposition 2's functional form (the increments grow from 12 → 24), so no curve-shape fit is pre-registered.

**Hypotheses (pre-registered, immutable after sign-off):**

- **H2a (RQ2, primary — empirical break-even, Tahoe).** With 3 new seeds and C ∈ {6, 8, 12, 16, 24, 32, 40} training
  contexts, the learned context correction's held-out residual-MSE gain over the context-free model is negative at
  C = 12, and its break-even C*_emp (zero crossing) is ≥ 24.
- **H2b (RQ2, primary — theory predicts the break-even).** The Proposition-2 break-even C*_theory, computed from d, σ²_ε,
  Var(i) and R² estimated on Tahoe **without using the gain curve**, lies within a factor of 2 of C*_emp.
- **H2c (RQ2, primary — genetic data are below break-even).** On public genetic screens, the model-free lower bound
  C*_min = d + 1 + d·σ²_ε/Var(i) exceeds C_public, the number of distinct public cell lines with a usable CRISPRi/KO
  Perturb-seq response in the M1 source registry.
- **H1 (RQ1 — ceiling).** On the H1 panel, the frozen pipeline (M1 + M2 + M3) attains ≥ 0.7 of the zero-shot transfer
  ceiling r_max = √[(Var(a) + Var(b)) / Var_reproducible] in noise-corrected correlation, and < 50 % of its remaining
  noise-corrected squared error lies in the transferable components (a_t + b_c).
- **H3a (RQ3).** On the H1 panel, after split-half noise control, the transfer residual (truth − pre-M3 prediction) puts
  a larger fraction of its energy on the destination's top-8 control PCs than on the other harness lines' top-8 control
  PCs (mean over six) and than on the top-8 PCs of gene-wise permuted destination controls.
- **H3b (RQ3).** Averaged over the 7 harness panels, M3 raises PDS more than each of five alternatives: other-line PC
  deflation, permuted-PC deflation, prediction-matrix PC deflation, Systema-style common-response removal, no correction.
- **H4 (RQ4).** On the H1 panel with the real scorer, adding H1GW (w = 1) gains more composite than the largest single
  modelling gain (M2, M3, or the best of the 17 architectures relative to its transfer baseline).

**Falsification Criteria** (each verdict is confirmed / rejected / inconclusive):

- **H2a.** *Confirmed* if the gain at C = 12 has a 95 % CI (bootstrap over held-out contexts, B = 2,000) entirely below
  0 **and** the lower 95 % bound of C*_emp is ≥ 24 (if the gain never crosses 0 by C = 40, C*_emp > 40 counts as ≥ 24).
  *Rejected* if the point estimate of the gain at C = 12 is ≥ 0 or the point estimate of C*_emp is < 24. Otherwise
  inconclusive.
- **H2b.** *Rejected* if C*_theory / C*_emp is outside [0.5, 2] (if C*_emp is censored at > 40, rejected only if
  C*_theory < 20). On rejection the paper keeps Proposition 1 (identifiability) and the empirical curves, and drops the
  quantitative form of Proposition 2, as the design specifies.
- **H2c.** *Rejected* if the point estimate of C*_min ≤ C_public. *Confirmed* if the lower 95 % bound (bootstrap over
  targets and split halves) of C*_min is > C_public. Otherwise inconclusive.
- **H1.** *Rejected* ("model-limited", reported as such) if the attained fraction is < 0.7 **and** ≥ 50 % of the
  remaining error is in the transferable components. If any method's noise-corrected r exceeds r_max with a 95 % CI
  excluding r_max, the ceiling estimator is declared invalid and H1 is reported as untestable. Confirmed if both parts
  of H1 hold. Otherwise inconclusive.
- **H3a.** *Rejected* if the 95 % target-bootstrap CI (B = 2,000; mean over 5 split seeds) of D_other or D_perm at k = 8
  is not entirely above 0, or if the control-only negative control yields D ≥ 50 % of the observed D.
- **H3b.** *Rejected* if the 7-panel mean ΔPDS (M3 − alternative) is ≤ 0 against any alternative. *Confirmed* if > 0
  against all five with 95 % CIs excluding 0 **and** the H1 real-scorer run has the same sign against all five.
  Otherwise inconclusive.
- **H4.** *Rejected* if composite(base + H1GW, w = 1) − composite(base) ≤ the largest single modelling gain, both
  measured with the real `cell-eval2 vcc2026` scorer on the H1 panel.

**Stop Criterion:** Analysis freeze on **2026-11-18**. Data collection is complete when:
- Tahoe: 7 C values × 6 folds × 3 seeds.
- Genetic parameters: computed on every line in `families.json`.
- Ceilings: computed on all 7 panels.
- RQ3: 7 panels × 5 splits × k ∈ {8, 24, 50}.
- Architecture replication: 6 architectures × 2 lines × seeds {1, 2}.
- RQ4: 5 w values.
No panel, C, seed or k is added after the freeze. Anything later is exploratory and goes in the Deviations Log.

**Method:** The frozen interpretable fusion pipeline (S18/F2 settings) and the existing learned-model results are put
under one theory. Tahoe (48 lines) tests whether the theory predicts the empirical break-even from independently
estimated parameters. The theory's model-free lower bound then says whether public genetic data are below break-even.
Split-half variance components give measurement and transfer ceilings. A noise-controlled residual analysis and a
correction comparison show where the remaining error lives and what a principled fix gains.

**Output target:** article (arXiv preprint, cs.LG primary + q-bio.GN cross-list; bioRxiv fallback), 8–10 pages main +
appendix, code on GitHub. Submission window 2026-11-28 … 12-01.
**Language:** en (manuscript); plan/wiki notes ko + en.

---

## Operational definitions (frozen with the hypotheses)

1. **Gain (Tahoe).** Mean over held-out contexts of [residual MSE without context − residual MSE with context], using
   the existing `dcp_*` protocol: FiLM r32, weight decay 1e-3, 6 folds with 8 held-out contexts each, identical budget
   across C. Total-response cosine gain is reported as a secondary metric. A seed changes both the training-subset draw
   and the initialization.
2. **C*_emp.** Zero crossing of the seed-and-fold-mean gain, linearly interpolated in ln C.
3. **d.** Participation ratio (Σλ)² / Σλ² of the descriptor matrix (contexts × descriptor features), centred across
   contexts. The descriptor is the exact input the learned model receives (Tahoe: the `dcp_model.py` context vector;
   genetic: control pseudobulk in log1p CP10K on the shared gene set).
4. **σ²_ε and Var(i).** These are estimated from cell-level split halves (A/B) on targets measured in ≥ 3 contexts:
   - Var(i): the cross-product ⟨i_A, i_B⟩, where i is the residual after removing μ, a_t and b_c.
   - σ²_ε: ½·E‖i_A − i_B‖², scaled to full depth.
   - R² (H2b only): adjusted R² of an OLS regression of i on the top-d descriptor PCs across all 48 Tahoe contexts, with
     the classical degrees-of-freedom correction.
5. **C_public.** The number of distinct cell lines/types with a public genome-scale or essential-scale CRISPRi/KO
   Perturb-seq response in the M1 registry. Same-line screens from different labs count once.
6. **r_max and attained fraction.** r_max is the noise-corrected correlation achievable by an oracle that knows a_t and
   b_c exactly but not i_{c,t}. The attained fraction is r_ours / r_max. The remaining-error share in the transferable
   components is ‖(μ + â_t + b̂_c) − δ̂‖² / (‖(μ + â_t + b̂_c) − δ̂‖² + Var(i)), noise-corrected.
7. **RQ3 estimator.** The control cells are split three ways: C_pc for PC estimation, C_a as the half-A baseline, and
   C_b as the half-B baseline. F(U) = Σ_t ⟨P_U r_A,t, P_U r_B,t⟩ / Σ_t ⟨r_A,t, r_B,t⟩, computed on the gene set G∩
   shared by the 7 panels. The weak half is defined by the M2 rank. The strength strata are quartiles of predicted
   amplitude.

## Data Sources

**Literature** (provisional bibkeys; Task 2 confirms or renames them and logs every rename):
- [[molina-zhang-2026]] perturbation response decomposition — *acquire pending*. **Gate:** full text read in Task 2.
- [[cell-quota-2026]] control-aware cell quota — *acquire pending*
- [[sciadv-benchmark-2026]] systematic comparison of perturbation models (Science Advances) — *acquire pending*
- [[dibaeinia-2026]] Position: virtual cells need context, not just scale (ICML 2026) — *acquire pending*
- [[ahlmann-eltze-2025]] linear baselines — *acquire pending*
- [[adapert-2026]] AdaPert (ICML 2026) — *acquire pending*
- [[state-2025]] STATE (bioRxiv 2025 / Cell 2026) — *acquire pending*
- [[tahoe-2025]] Tahoe-100M (bioRxiv 2025 / Cell 2026) — *acquire pending*
- [[systema-2025]] Systema — *acquire pending*
- [[vcc-2026]] VCC 2026 commentary (Cell) — *acquire pending*
- [[replogle-2022]], [[nadig-2025]], [[orion-2025]], [[pan-2026]] — *acquire pending*
- Dataset citations to resolve in Task 2: KOLF2.1J CRISPRi, HipSci (Feng), CD4 T-cell genome-scale CRISPRi,
  DLD-1 (GSE337988), GSE306429, Jiang CRISPRi, LINCS shRNA.

**Data on disk** (frozen in Task 4):
- **Tahoe subset:** 4.43 M cells, 48 lines, 13,011 conditions.
  - `vcc2026_data/reports/simple_residual_2026-09-16/data/tahoe/table.npz`
  - `vcc2026_data/data/processed/drug_context_pretrain/tahoe.npz`
  - Protocol and seed-0 results: `vcc2026_data/reports/drug_context_pretrain_2026-09-20/` (`scaling_table.json`,
    `context_splits.json`)
  - Code: `vcc2026_data/scripts/dcp_{data,model,prepare,select}.py`
- **Harness truth (cell level):**
  - H1: `vcc2026_data/reports/transfer_v2_2026-09-23/h1/real_thin.h5ad`
  - OR_HCT116X: `vcc2026_data/reports/sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad`
  - KOLF: `checkpoints/state/kolf_real.h5ad`
  - RN Jurkat: `checkpoints/state/jurkat_real.h5ad`
  - RN K562 / RN RPE1 / RN HepG2:
    `vcc2026_data/raw/scperturb/{ReplogleWeissman2022_K562_essential,ReplogleWeissman2022_rpe1,NadigOConner2024_hepg2}.h5ad`
- **Evaluator:** offline packs `vcc2026_data/reports/era_2026-10-07/packs/*_<panel>.npz`; evaluator
  `vcc2026_data/scripts/sch_sumeval.py` (exact only for budget-invariant arms).
- **Pipeline:**
  - `vcc2026_data/scripts/tv2_final.py`: M3 = `--sum-defl 8,0.75,0.5`; M2 = `sum.tamp` β 0.5
  - `vcc2026_data/scripts/ml_ctx_curve.py` (source-count curve)
  - `vcc2026_data/scripts/ml_h1gw_value.py` and `sources_new/H1GW.npz`
- **Existing results:**
  - `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md` (17 architectures, 104 real-scorer runs, seed 0)
  - `vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md` (learned 6→35 curve)
  - `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (ceilings)
  - `vcc2026_data/reports/PROJECT_HISTORY.md` (negative-result record)
  - `vcc2026_data/reports/drug_gene_concordance_2026-09-16/RESULTS.md` (appendix)
- **Official record:** `vcc2026_data/reports/official_scores.csv` (29 rows) + `submission/official_scores_new.csv`. Task 4
  reconciles the submission count (the design says 36).

## Tasks

### Task 0: Scaffold the research wiki inside AIVC (11-05)
**Files:** `CLAUDE.md`, `knowledge/{sources,entities,concepts,synthesis,_meta}/`, `knowledge/_meta/log.md`,
`knowledge/_meta/gate-overrides.log`, `output/{data-analysis,article,bibtex,code}/`, `scripts/lint-wiki.py`

- [ ] Copy `knowledge/`, `output/`, `schema/`, `scripts/lint-wiki.py`, and `scripts/vendor/` from the plugin's
      `research-project-template` with `cp -rn` (never overwrite existing `scripts/*`).
- [ ] Write `CLAUDE.md` with frontmatter `methodology: quantitative`, `discipline: "ML for computational biology"`,
      `languages: [en, ko]`. Its body points to `AGENTS.md` for the compute rules.
- [ ] Log in `knowledge/_meta/log.md`: the plan is pre-registered, and the earlier destination-state design and plan
      are superseded (never pre-registered).
- [ ] Check `quarto --version` and `xelatex --version`; record the result. If either is missing, Task 27 uses plain
      LaTeX at `output/article/main.tex`.
- [ ] `python scripts/lint-wiki.py` exits 0.

### Task 1: External gates (11-05 … 11-07)
**Files:** `input/notes/external-gates.md`

- [ ] VCC 2026 terms on publishing results and data use: quote the clause and give its URL. Record whether publishing
      before the final results is allowed.
- [ ] arXiv cs.LG endorsement: draft the request text and list 2–3 candidate endorsers (authors of A-graded cs.LG
      papers). The user sends it. If there is no endorsement by 11-21, switch to bioRxiv and log the switch.

### Task 2: Literature search + overlap gate (11-05 … 11-07)
**Files:** `input/bibliography/literaturguide.md`, `output/bibtex/references.bib`, `input/notes/overlap-molina-zhang.md`

- [ ] Run `literature-review` on: zero-shot / unseen-cell-line perturbation prediction; response decomposition and
      noise ceilings; sample complexity in the number of contexts / environments (incl. domain generalization and
      multi-environment learning theory); learned models vs baselines; drug ↔ genetic perturbation concordance.
- [ ] Grade A/B/C. The A set = the literature above plus any work that reports a context-count break-even or a
      context sample-complexity bound. **If such a work exists, stop and raise it with the user before Task 8.**
- [ ] Read the Molina & Zhang full text. In `overlap-molina-zhang.md`, list exactly what it already shows for RQ1 and
      RQ3. If it already compares transfer residuals against other lines' control axes, shrink RQ3 to the correction
      comparison (H3b) and log the change.
- [ ] Fix the bibkeys in `references.bib` and log any rename of the provisional keys.

### Task 3: Acquire source PDFs (11-06 … 11-08)
**Files:** `<library>/pdf/*.pdf`, `input/bibliography/acquisition-todo.md`

- [ ] Run `acquire-sources` on the A+B set.
- [ ] The user manually downloads paywalled originals. Re-run `acquire-sources` until every A source is on disk or
      explicitly deferred.

### Task 4: Freeze inputs, pipeline and source registry (11-05 … 11-06)
**Files:** `output/code/INPUTS.md`, `output/code/inputs.sha256`, `output/data-analysis/common.py`,
`output/data-analysis/results/panels.json`, `output/data-analysis/results/families.json`,
`output/data-analysis/results/official_record.csv`

- [ ] Write sha256 for every input of Tasks 6–21 to `inputs.sha256`. Back them up to
      `/media/wonmin-kim/Wonmin/aivc_backup_2026-10-27/paper_inputs/`.
- [ ] `common.py` loaders: panel truth (cell level), controls, packs, the Tahoe table, and G∩. Write the gene counts to
      `panels.json`.
- [ ] `families.json`: the M1 source registry, giving each family its cell line, lab, scale, and target count, plus
      **C_public** (definition 5).
- [ ] Check leave-destination-out for each panel and write each panel's excluded source list to `panels.json`.
      **If there is any leak, stop and log it.**
- [ ] Merge the two official-score CSVs into `official_record.csv` and record the true submission count.
- [ ] Record the `tv2_final.py` sha and the exact S18/F2 arm string in `INPUTS.md`.

### Task 5: Ingest A-graded sources (11-07 … 11-12, parallel)
**Files:** `knowledge/sources/*.md`, `knowledge/entities/*.md`, `output/bibtex/references.bib`

- [ ] Ingest via `ingest-source` with focus "limits / ceilings / context count / what is already shown": molina-zhang-2026,
      cell-quota-2026, sciadv-benchmark-2026, dibaeinia-2026, ahlmann-eltze-2025, adapert-2026, state-2025, tahoe-2025.
- [ ] Ingest with focus "data / benchmark facts": systema-2025, vcc-2026, replogle-2022, nadig-2025, orion-2025, pan-2026.
- [ ] On each source page, state what overlaps with H1–H4.

### Task 6: Tahoe context-curve runner (RQ2; 11-06 … 11-07)
**Files:** `output/data-analysis/rq2_tahoe_curve.py`, `output/data-analysis/tests/test_rq2_tahoe_curve.py`

- [ ] Wrap `dcp_model.py` / `dcp_data.py` without changing them. The arguments are C, a seed, and a fold from
      `context_splits.json`. The wrapper draws the C training contexts from each fold's training pool with that seed.
- [ ] Output per (C, seed, fold, held-out context): residual MSE with and without context, and the cosines.
- [ ] Test: C = 40, seed 0 reproduces the existing `scaling_table.json` row (MSE gain +1.9e-05) to within 10 %.
- [ ] Record the descriptor (the exact context input vector and its dimension) in the test docstring. Task 8 uses it.

### Task 7: Run the Tahoe curve and H2a (11-07 … 11-09)
**Files:** `output/data-analysis/results/rq2_tahoe_curve.parquet`, `output/data-analysis/results/rq2_h2a.json`

- [ ] Run C ∈ {6, 8, 12, 16, 24, 32, 40} × folds 0–5 × seeds {1, 2, 3} as a detached CPU job (seed 0 is kept but
      reported separately).
- [ ] Compute the gain and C*_emp with a bootstrap over held-out contexts (B = 2,000). Write the H2a verdict to
      `rq2_h2a.json`.

### Task 8: Independent parameter estimates on Tahoe and H2b (11-08 … 11-10)
**Files:** `output/data-analysis/rq2_params.py`, `output/data-analysis/results/rq2_params_tahoe.json`,
`output/data-analysis/results/rq2_h2b.json`, `output/data-analysis/tests/test_rq2_params.py`

- [ ] Estimate d from the Task 6 descriptor (definition 3).
- [ ] Split Tahoe conditions into cell halves. Estimate σ²_ε and Var(i) (definition 4) on drugs measured in ≥ 3 lines.
- [ ] Compute the adjusted R² of i on the top-d descriptor PCs across the 48 contexts.
- [ ] Compute C*_theory and C*_min, and write the H2b verdict. **These scripts must not read `rq2_tahoe_curve.parquet`**
      (enforced by a test that greps the imports).

### Task 9: Genetic parameters and H2c (11-10 … 11-12)
**Files:** `output/data-analysis/results/rq2_params_gene.json`, `output/data-analysis/results/rq2_h2c.json`

- [ ] Use the same `rq2_params.py` on the genetic lines in `families.json`, restricted to targets shared by ≥ 3 lines
      (primarily the essential overlap of RN K562 / RPE1 / HepG2 / Jurkat, plus H1, H1GW, KOLF and Orion where they
      overlap).
- [ ] Compute d (from control pseudobulk across lines), σ²_ε, Var(i), and C*_min with a bootstrap CI. Write the H2c
      verdict against C_public.
- [ ] Sensitivity analysis: C*_min with the descriptor restricted to the top 8 / 24 / 50 control PCs.

### Task 10: Replicate the learned-model results with seeds (11-09 … 11-14)
**Files:** `output/data-analysis/rq2_arch_replicate.sh`, `output/data-analysis/results/rq2_arch_seeds.csv`

- [ ] Choose 6 of the 17 architectures from `arch_screen_2026-09-22/FINDINGS.md`: the top 3 by seed-0 composite, plus
      the hypernetwork (141 k), the Transformer (86 k), and cross-attention.
- [ ] Rerun them on K562 + HepG2 with seeds 1 and 2 through the same real-scorer protocol: 24 runs ≈ 10 h, serial and
      detached.
- [ ] Report the mean ± sd over 3 seeds against the refit/fusion baseline. State whether any learned model beats the
      baseline with a CI excluding 0. This is a descriptive consistency check for H2c.

### Task 11: Transfer-pipeline source-count curve (RQ2 contrast + RQ4 input; 11-10 … 11-13)
**Files:** `output/data-analysis/rq2_source_curve.py`, `output/data-analysis/results/rq2_source_curve.parquet`

- [ ] Extend `ml_ctx_curve.py`: for k = 1 … K, use R = 20 random family subsets (all of them if there are fewer). The
      k = K−1 subsets are exactly the leave-one-family-out (LOFO) set. Rebuild M1 + M2 and score PDS and MSE.
- [ ] Exactness spot check: run the real scorer on H1 at k = 1, ⌈K/2⌉ and K. If PDS deviates by more than 0.002, report
      the curve in relative terms and log it.
- [ ] Report descriptively: the log-linear vs saturating fits (AIC) and the LOFO marginal value with a bootstrap CI.
      There is no hypothesis attached to this task.

### Task 12: Variance components on the harness panels (RQ1; 11-11 … 11-13)
**Files:** `output/data-analysis/rq1_components.py`, `output/data-analysis/results/rq1_components.parquet`, `output/data-analysis/tests/test_rq1_components.py`

- [ ] For each panel, estimate μ, a_t (from the other lines' truth), b_c, Var(i), σ²_ε and Var_reproducible by split
      halves, on targets shared with ≥ 2 other lines.
- [ ] Compute the measurement ceiling (split-half reliability) and r_max (definition 6) for each panel.
- [ ] Unit test on synthetic data with known components: the components are recovered to within 5 %.

### Task 13: Methods vs ceilings and H1 (11-13 … 11-14)
**Files:** `output/data-analysis/rq1_ceiling.py`, `output/data-analysis/results/rq1_ceiling.csv`

- [ ] Compute the noise-corrected r for: the frozen pipeline (Task 14 predictions), plain cross-line mean transfer, and
      the 17 architectures on the panels where their predictions exist.
- [ ] Compute the attained fraction and the remaining-error share in the transferable components (definition 6).
      Write the H1 verdict.
- [ ] Real-scorer PDS for H1: the pipeline vs the oracle-transfer predictor (μ + â_t + b̂_c), 2 runs.

### Task 14: Frozen pre-M3 and M3 predictions (RQ3; 11-06 … 11-08)
**Files:** `output/data-analysis/rq3_build_preds.py`, `output/data-analysis/results/preds/<panel>_{prem3,m3}.npz`

- [ ] For each panel, leave-destination-out, run the frozen recipe with `--sum-defl` off and on. Save the pseudobulk
      deltas (target × gene, log1p CP10K) and the M2 ranks.
- [ ] Sanity check: M3 − pre-M3 PDS from `sch_sumeval` is > 0 on 7/7 panels. Record whether this agrees.
- [ ] Check that the pre-M3 prediction depends on destination controls only through the baseline profile. If it
      depends on them more, rebuild with C_a ∪ C_b only.

### Task 15: Split-half residual-energy estimator (RQ3; 11-08 … 11-09)
**Files:** `output/data-analysis/rq3_energy.py`, `output/data-analysis/tests/test_rq3_energy.py`

- [ ] Implement definition 7, with these subspaces: U_dest; U_j (each of the 6 other lines); U_perm (20 permutations);
      U_haar (reported only).
- [ ] Synthetic tests:
  - [ ] zero residual → E[F numerator] ≈ 0;
  - [ ] a residual planted on U_dest → F(U_dest) ≫ F(U_perm);
  - [ ] a shared-baseline variant shows the bias that the three-way split removes.
- [ ] Target bootstrap (B = 2,000) of D_other and D_perm.

### Task 16: Run E1a / H3a (11-09 … 11-11)
**Files:** `output/data-analysis/results/rq3_energy.parquet`, `output/data-analysis/results/rq3_negctrl.parquet`

- [ ] Run 7 panels × k ∈ {8, 24, 50} × 5 splits (detached).
- [ ] Negative control: 300 pseudo-targets drawn from the destination controls with Δpred = 0.
- [ ] Report the predicted-amplitude quartile strata. Write the H3a verdict.

### Task 17: Run E1b / H3b (11-11 … 11-15)
**Files:** `output/data-analysis/rq3_corrections.py`, `output/data-analysis/results/rq3_corrections.parquet`

- [ ] Run the arms on the weak half (M2 rank) with k = 8 and η = 0.75: M3_dest, M3_other (6 lines; primary = their
      mean), M3_perm, predPC, Systema-style (η ∈ {0.25, 0.5, 0.75, 1.0}, keep the best), none.
- [ ] Budget-invariance check on H1 for M3_dest and Systema-style: `sch_sumeval` vs the real scorer. If PDS differs by
      more than 0.001, that arm uses the real scorer on every panel.
- [ ] Per-target PDS → bootstrap within each panel, with panels weighted equally.
- [ ] Real scorer on H1: none, M3_dest, the best M3_other, M3_perm, predPC, Systema-best (6 runs ≈ 2.5 h).
- [ ] Secondary: leave-one-panel-out choice of (k ∈ {8, 24}, η ∈ {0.5, 0.75, 1.0}) for M3. Write the H3b verdict.

### Task 18: Same-line data value and H4 (11-13 … 11-16)
**Files:** `output/data-analysis/rq4_sameline.py`, `output/data-analysis/results/rq4_sameline.parquet`

- [ ] On H1 with the real scorer: base (frozen S18/F2) with H1GW at w ∈ {0, 0.25, 0.5, 0.75, 1} (5 runs), plus
      base − M2 and base − M3 (2 runs).
- [ ] Best-architecture gain: from `arch_screen_2026-09-22/FINDINGS.md` updated by Task 10, relative to its own transfer
      baseline (note if the scorer setup differs).
- [ ] Per-context value: the LOFO value from Task 11, converted to composite on H1.
- [ ] Related-line cases (cut first): KOLF ← H1GW and OR_HCT116 ← DLD-1, w ∈ {0, 0.5, 1}.
- [ ] Write the H4 verdict. Report the H1GW / per-context ratio descriptively.

### Task 19: Official record (11-14)
**Files:** `output/data-analysis/e4_official.py`, `output/data-analysis/results/e4_official.csv`

- [ ] From `official_record.csv`: the score trajectory and rank 129 / 1,635. Attribute gains only from single-change
      submissions: M2 is S14fk2 → S15 (+0.0106); the M3 bundle is S17 → S18.
- [ ] Local → official shrinkage cases: FiLM S3 (+0.018 local → −0.0005 official) and S4 (+0.0247 → +0.0013).
- [ ] Use leaderboard snapshots only in aggregate (quantiles), never with other teams' names. Add final-phase results
      only if they are published before 11-25.

### Task 20: Appendix — drug→gene bound (11-14 … 11-15)
**Files:** `output/data-analysis/appx_drug_gene.py`, `knowledge/synthesis/appx-drug-gene-bound.md`

- [ ] From `drug_gene_concordance_2026-09-16`: derive the correlation upper bound for target-identity-mediated transfer
      (via Cauchy–Schwarz on the noise-corrected concordance), with a bootstrap CI.
- [ ] State explicitly that the earlier gene-side drug-pretraining experiment is invalid (scrambled context vectors) and
      is neither cited nor used.

### Task 21: Figures (11-15 … 11-17)
**Files:** `output/data-analysis/figures.py`, `output/article/figures/fig{1..6}.pdf`

- [ ] Load the `dataviz` skill first.
- [ ] Figures:
  - [ ] Fig 1: the decomposition and setup schematic.
  - [ ] Fig 2: the Tahoe gain vs C with C*_emp and C*_theory.
  - [ ] Fig 3: genetic C*_min vs C_public, with the learned models vs baseline as an inset.
  - [ ] Fig 4: methods on the ceiling scale.
  - [ ] Fig 5: RQ3 energy fractions and the correction ΔPDS.
  - [ ] Fig 6: value comparison (H1GW vs per-context vs modelling).

### Task 22: Synthesis pages (11-18 … 11-19)
**Files:** `knowledge/synthesis/theory-context-complexity.md`, `knowledge/synthesis/related-work-positioning.md`,
`knowledge/synthesis/results-summary.md`

- [ ] Theory page: Proposition 1 (identifiability: rank ≤ C−1), Proposition 2 (break-even and its lower bound), and the
      nonparametric version (rate with sample size = C). Give full derivations and assumptions.
- [ ] Positioning page: for each focus work, what it shows vs what this paper adds. Cite at least 12 sources as
      [[bibkey]].
- [ ] Results summary: every verdict, with effect sizes, CIs, and deviations. Set `status: stable` after the user
      reviews it. `lint-wiki.py` exits 0.

### Task 23: Draft — Introduction + Related Work (11-19 … 11-20)
**Files:** `output/article/sections/01-intro.qmd`, `output/article/sections/02-related.qmd`

- [ ] Via `drafting-manuscript`: the Introduction in 1,000–1,400 words and Related Work in 900–1,200 words, cited as
      `[@bibkey]`.

### Task 24: Draft — Theory + Setup (11-20 … 11-21)
**Files:** `output/article/sections/03-theory.qmd`, `output/article/sections/04-setup.qmd`

- [ ] Theory section: 1,200–1,700 words (propositions, the lower bound, assumptions; proofs go to the appendix).
- [ ] Setup section: 1,000–1,400 words (data, pipeline M1–M3, estimators, ceilings, the pre-registration statement and
      honesty statement).

### Task 25: Draft — Results (11-21 … 11-23)
**Files:** `output/article/sections/05-results.qmd`

- [ ] 2,000–2,800 words, in order RQ2 → RQ1 → RQ3 → RQ4 → official record. Every number must trace to a file in
      `output/data-analysis/results/`.

### Task 26: Draft — Discussion, Limitations, Appendix (11-23 … 11-25)
**Files:** `output/article/sections/06-discussion.qmd`, `output/article/sections/99-appendix.qmd`

- [ ] Discussion + limitations: 900–1,300 words. Cover: the linear-Gaussian theory; the Tahoe subset; one competition;
      seed-0 history; M3 hyperparameters not held out; shallow panels.
- [ ] Appendix: 3,000–4,500 words. Cover: proofs, per-panel tables, the leave-one-panel-out analysis, estimator tests,
      the drug→gene bound, the input manifest, and the competition-rule note.

### Task 27: Render, code release, review, submit (11-26 … 12-01)
**Files:** `output/article/main.qmd`, `output/article/main.pdf`, `output/code/README.md`, `output/code/environment.yml`

- [ ] Render the PDF (quarto, or the LaTeX fallback); the main text is ≤ 10 pages.
- [ ] Code release: the scripts from Tasks 6–21, `INPUTS.md`, `inputs.sha256` and `environment.yml`. From a clean
      checkout, rerun Task 8 and the H1 part of Task 16, and confirm the numbers match to 4 decimals.
- [ ] Run `requesting-peer-review` (constructive, then adversarial). Revise and re-render.
- [ ] Submit to arXiv (or bioRxiv), 11-28 … 12-01, with the user's go-ahead.

## Verification

The project is complete when:

- [ ] `python scripts/lint-wiki.py` exits 0
- [ ] ≥ 14 A-graded sources are ingested, each with `status: stable` or `review`
- [ ] `theory-context-complexity.md`, `related-work-positioning.md` and `results-summary.md` have `status: stable`
- [ ] Every result file named in Tasks 7–20 exists, and every number in `05-results.qmd` matches it (checked by a grep
      script)
- [ ] `pytest output/data-analysis/tests` passes (Tasks 6, 8, 12, 15)
- [ ] `main.pdf` renders without errors, with ≤ 10 main-text pages
- [ ] A clean re-run reproduces `rq2_params_tahoe.json` and the H1 rows of `rq3_energy.parquet` to 4 decimals
- [ ] Peer review round completed via `requesting-peer-review`
- [ ] H2a, H2b, H2c, H1, H3a, H3b and H4 are each marked confirmed / rejected / inconclusive per their criteria
- [ ] The preprint is submitted, and its arXiv ID or bioRxiv DOI is recorded in `knowledge/_meta/log.md`

**Fallback cuts if late (in order; each one is logged):**
1. RQ4 related-line cases
2. Task 10 seed replication (report seed 0 with the caveat)
3. RQ3 reduced to H1 + 3 panels
4. Task 11 reduced to R = 10

RQ2 (Tasks 6–9) and RQ1 (Tasks 12–13) are never cut.

## Deviations Log

Append to `knowledge/_meta/log.md` whenever scope, method, or data sources change after pre-registration. Reference
the plan slug `zero-shot-limits-context-complexity-plan`. Any result reached after a deviation is labelled
**exploratory** in the manuscript.

## Amendments (pre-data; approved by the user 2026-11-05; the text above is left unchanged)

**A1 — H2a/H2b quantity and CI** (reason: an independent review of Task 6 showed the frozen `dcp` protocol partly
guarantees a negative gain at small C — checkpoint selection can never choose the zero-residual model, and the
"without context" arm is the analytic zero residual, not a network trained under the same budget; logged
2026-11-05 before any Task 7 data were collected, except one disclosed smoke run):
- Primary quantity for H2a and for C*_emp in H2b: **context-pathway gain** = residual MSE of the constant-descriptor
  ablated network (identical data, budget, early stopping, selection, seed) − residual MSE of the context network,
  unweighted mean over held-out contexts. Thresholds unchanged (gain at C = 12 entirely below 0; C*_emp ≥ 24; H2b
  factor 2).
- CI: cluster bootstrap — resample the 18 (seed, fold) runs, then held-out contexts within each, paired across the
  nested C values (B = 2,000). Also report per-seed C*_emp and the between-seed SD of the gain at C = 12.
- If C*_emp is censored (> 40), H2b is reported as **"untestable (censored)"**, because C*_theory ≥ d + 1 ≈ 21.6 by
  construction and the registered rule could then never reject.
- The registered quantity (gain vs the analytic zero residual under the frozen protocol) is still computed and
  reported, with the confound stated. An exploratory arm with step 0 allowed in checkpoint selection is labelled
  exploratory.

**A2 — C_public counted two ways for H2c** (reason: the registry's HipSci/Feng family pools donor iPSC lines; the
definition is ambiguous for donor lines):
- Count 1 (cell types/lines): as computed in Task 4 (`families.json`: C_public = 11).
- Count 2 (donor lines): every distinct cell line with a usable per-line response — ≥ 300 targets with ≥ 20 cells
  each — counted separately; kolf_2/kolf_3 are merged with KOLF2.1J (same line); computed in Task 9 before any
  C*_min estimate is looked at.
- H2c is **confirmed** only if it holds under both counts; if it holds only under Count 1 it is reported as
  "holds for distinct cell types, not for donor lines"; otherwise rejected/inconclusive per the original rule.

**A3 — corrections forced by Task 4** (no change to hypotheses):
- Truth paths are those in `output/data-analysis/results/panels.json` (RN_K562 truth was built from K562 GWPS;
  KOLF/Jurkat truth paths in the Data Sources section above were wrong).
- KOLF leak: the FENG (HipSci) family contains kolf_2/kolf_3 (parental/same-donor of KOLF2.1J; 11–18 % of FENG
  cells) and was not excluded for the KOLF panel. All analyses with KOLF as destination exclude kolf_2/kolf_3 cells
  from every FENG arm (or drop FENG for KOLF if a cell-level rebuild is infeasible; whichever is used is logged).
- Official submissions: 38 (36 validation + 2 final), `official_record.csv`.
- Framing: Prophet (ji-2024-universal-prediction) is prior empirical evidence; Proposition 1 is cited from
  domain-generalisation theory (rosenfeld-2021-risks-irm); the claimed contribution is the explicit break-even C*
  with its noise lower bound, its measurement, and the de-confounded learning curve.
