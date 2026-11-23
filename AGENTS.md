# AIVC / Virtual Cell Challenge 2026 — handoff

## 2026-10-31 PAPER (preprint for a PhD application) — design APPROVED, plan document still to write
- **2026-11-05 PIVOT (user chose option B, "integrate"):** new design `input/ideas/zero-shot-limits-context-complexity-design.md` (status APPROVED 2026-11-05); plan `input/ideas/zero-shot-limits-context-complexity-plan.md` written (status PRE-REGISTERED 2026-11-05; analysis freeze 11-18; executing) = "data-limited, not model-limited": RQ2 context sample complexity (main; Tahoe pre-registered prediction at C=32/40), RQ1 ceilings, RQ3 = old E1/M3, RQ4 = old E3; PIE/STATE plug-in dropped. The old destination-state design and plan are marked superseded. VCC final: F2 remains the last entry (nothing after 10/28).

- User goal: an arXiv preprint (bioRxiv fallback) supporting PhD applications in ML methods for virtual-cell
  perturbation prediction, THIS cycle (December deadlines) -> submit 2026-11-28..12-01; single author; workshop optional
  (v2 candidate: an ICLR 2027 workshop, suggested deadline 2027-02-01). Compute: CPU only except short GPU windows.
- Approved design: `input/ideas/zero-shot-transfer-destination-state-design.md` (status approved). Next step: run
  `research-superpowers:writing-research-plan` on it -> `input/ideas/zero-shot-transfer-destination-state-plan.md`.
  The plan file currently on disk is an INCOMPLETE stub (24 lines, the writing session was interrupted) -- overwrite it.
- Literature check (2026-10-31): "learned models do not beat cross-line mean transfer in unseen contexts" is ALREADY
  SHOWN (AdaPert ICML 2026 Table 13; Pan 2026; STATE; Molina & Zhang 2026) -> not a claim of ours. Open: whether the
  transfer residual concentrates on the destination's OWN control axes vs other lines' axes (noise-controlled), the
  deflation vs Systema-style corrections, a controlled source-count scaling curve, same-line data value vs modelling.
- Facts for the plan: cell-level truth (needed for split halves) exists for H1 (`reports/transfer_v2_2026-09-23/h1/real_thin.h5ad`),
  OR_HCT116 (`reports/sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad`), KOLF (`harness/kolf_real.h5ad`), RN Jurkat
  (`harness/jurkat_real.h5ad`); exact offline packs `reports/era_2026-10-07/packs` (7 panels, incl. ctrlpc_*). STATE
  `checkpoints/state/ST-SE-Replogle` has NO model weights locally and uses a one-hot perturbation vocabulary (RN panels
  only); plug-in check (E5) therefore prefers PIE `arcinstitute/PIE_replogle_xdataset` (never trained on Replogle or
  KOLF -> honest on RN x4 + KOLF; needs PIE_sources ~44 GB -> external drive).

## 2026-10-27 FINAL PHASE BUILD F1 (user: "진행해줘"; `vcc2026_data/reports/final_2026-10-27/README.md`)

- Validation closed 10/22 (S18 0.2545 ended rank 129; top-10 >= 0.319). Final set D/E/F downloaded with vcc-cli **0.3.0**
  (`data/final_2026-10-22/controls`); deadline **2026-11-05 23:59 UTC**. Identity: D AML myeloid (MOLM-13), E U-87MG, F
  hTERT fibroblast; all somatic; no public same-line screen (ELEM-1200RS bulk KO is gated, needs the user to accept terms).
- **F1 UPLOADED 2026-10-27 15:45Z (user go-ahead), entry HPPRo7pAvkuTyzfik2Xl, model b24-460545, status published, final score hidden until winners are announced; receipt `reports/submission_2026-09-23/F1.upload_receipt.json`**: `reports/submission_2026-09-23/F1.vcc`, sha256 4605456755d4626f71a98c9a6e00d28afadd598244b45ad5c840fb16047c90fd,
  vcc prep passed. = 10/22 command + `--rn-sources --fam-weights DLD1=1 --sum-defl 8,0.75,0.5 --sum-cens 5,10,20,100 --de-drc 0.8,1.1`.
- **F2 UPLOADED 2026-10-28 21:05Z = CURRENT FINAL ENTRY (last entry counts), entry oYO6WAS8x4FgQt9ypgmq, model b24-0c2a01, published; receipt reports/submission_2026-09-23/F2.upload_receipt.json. F1 superseded.** Was: (F1 + ew kept; `reports/submission_2026-09-23/F2.vcc`, sha 0c2a0180...; upload: `bash ../scripts/run_detached.sh aivc-final-F2-upload reports/final_2026-10-27 /usr/bin/bash reports/final_2026-10-27/upload_F2.sh` from vcc2026_data, only with the user's go-ahead). PIE rejected (cross-dataset discrimination 0.53).
- **2026-10-28 big-lever search (compute-careful, CPU-only; `reports/final_2026-10-27/README.md`):** the 0.31 tier's edge is MSE
  (+0.025) = strong-target DIRECTION (global cosine ~0.4 vs our ~0.2); amplitude is at the kink for everyone. Exact offline screens
  of every untested public lever ([another team]'s K562 response-space kNN, learned per-gene transfer weights, VIPerturb-seq family,
  ew re-test with budget match, stacks with the deflation) all land at <= +0.005 composite; ew (+gw) adds OR_X PDS +0.012 over F1
  but H1 +0.0007. VIP as a source is rejected. **Final rule: only the LAST final entry counts (2/day).** PIE (Arc, public) is the
  one DE-side candidate with external evidence; needs a GPU window.
- The 10/11 cleanup had deleted three S4/residual/denoise inputs; restored as described in that README (stub table via
  `TV2_BASE_TABLE`, regenerated fam_* and full/K562). **Never delete `simple_residual_2026-09-16/data/table`,
  `disc_residual/fam_*`, `denoise_src/full` again.**

## 2026-10-03 PDS/MSE ROUND (user: "PDS, MSE를 높이기위해서는 어디를 건드려야해" / "mean transfer을 개선할수는없어?"; FINDINGS §25)

- **SCORER FACT, from source (`streaming_bulk.py:176`, `prep.bulk_lognorm_means`, `metrics/delta.py:1078`):** the
  `bulk_lognorm` comparator is `log1p(TS * S_g / sum_g S_g)` -- a function of the 300 group SUMS alone. "fixed raw group
  totals AUTOMATICALLY fix the pseudobulk". So PDS and MSE cannot be touched by the within-cell emission, and there is
  NO `cpc`-style functional mismatch on the sum side (hypothesis falsified from source; do not re-test).
- **PDS and MSE lose in DISJOINT target populations** (H1, S16, scorer space; rho median 0.085): MSE numerator is 97 %
  in the strongest 3 deciles (rho 0.32-0.18); PDS is 0.99 there and 0.86 / 0.74 / 0.67 in the weakest 3 (rho ~0.01),
  126/300 targets already rank-1. **No single lever raises both** -- that premise was wrong.
- **NEW RULE -- `J = 2.23 dPDS - 1.05 dMSE` IS BLIND TO NMAE.** It was derived for arms that leave the DE side alone; a
  per-target AMPLITUDE arm moves NMAE through the `cpc` pin. OR_X `tzb0.25`: J +0.0398, real composite **-0.0108**, the
  whole gap being NMAE +0.0373. On H1 (NMAE -0.0001) J/6 = +0.0020 vs real +0.0022. **Never decide an amplitude arm on
  J.** Every past amplitude screen used it (incl. `tamp`/S15, which survived only because official NMAE behaves like H1's).
- **NEW RULE -- the offline evaluator is EXACT iff the arm is budget-INVARIANT.** `tzb0.25` and `tzb0.6` on H1 both
  reproduced PDS and MSE to four decimals; the budget-CHANGING `tz0.25x0.9` sign-flipped MSE (offline -0.0211, real
  +0.0055). Cause: `sch_sumeval.evaluate` stands in the REAL panel's jackknife for `C_pred_hat` while the #247 cap binds
  on every row. So this axis is now screenable in seconds instead of 25 min/panel -- but only when budget-invariant.
- **Everything tried this round closed; net score 0.** per-target amplitude budget half (unreachable) / shape half
  (`tzb`, PDS +3/3 but R21 **-0.0042**, DE-following repair made NMAE WORSE +0.0530); better rho-hat (v2 +0.0003, v3
  **null**: mean Spearman .6722 -> .6736, the new agreement features are redundant); strength-gated hubness `disw`
  +0.0007; neighbour-cut ensemble `nbe` +0.0005 (**no member diversity** -- lists are <= 20 long so top-20/40/80 are
  identical); continuous DE calibration `drp` -0.0013. **Ensemble diversity axes are exhausted**: `min_cells` shipped,
  `pcr-lambda`/centring provable no-ops, neighbour cut flat, ROW bootstrap INFEASIBLE (`replogle_k562_gwps` has exactly
  1 row/target). New tokens `tz`, `tzb`, `nbe`, `drp`, `drc`, `disw` are all default-off; S16 arm parses bit-identically.
- **Correction to §25.4 made in-round:** "amplitude is closed" was first reached with a GREEDY oracle that collapses
  `B_pred` (r 0.89 -> 0.287). The budget-held oracle is +0.0139 in J units on H1 -- but J is NMAE-blind and it is an
  oracle, so it is inflated twice. The axis is closed, now for the right reasons.
- **ONE POSITIVE (§25.13): the LFC PSEUDOCOUNT ensemble `cens`.** The transfer's COMBINATION side is saturated, but its
  ESTIMATOR side had an untried axis: the fusion part's pseudocount is fixed at c=50, and `_C3` (c 10/20/100) sat
  DEFINED in `ENS["mcwc"]` with zero scored arms. Built as `delta_s4(..., c_list=)`, token `cens5.10.20.100`, unioned
  with the `min_cells` members. Real scorer: H1 +0.0003, OR_HCT116X +0.0042, RN_Jurkat +0.0063 -- **3/3 positive**,
  with PDS / NMAE / reach positive on 3/3 (reach +0.0033..+0.0042) and only MSE panel-dependent. **R21 +0.0024 < the
  +0.005 bar -> RIDE-ALONG only**, same class as `drc 0.8/1.1` (H1 +0.0018, 3/3). Valid c values are only
  `tv2_sources.CS = (1,5,10,20,50,100)`; a parse-time guard enforces it. `cens` absent => bit-identical to before.
  A DE-side-only variant (sums pinned, PDS/MSE invariant) estimates to R21 ~+0.0033 -- still under the bar and it
  gives up the 3/3 PDS gain, so it was not built.
- **RIDE-ALONG LIST, now PORTED and audited** (never a slot on their own; both default off):
  `tv2_final.py --de-drc 0.8,1.1 --sum-cens 5,10,20,100`. They were harness-only until 2026-10-03, i.e. the list was
  unusable on 10/22. `cens` needed NO builder change -- `parts[].cfg_ens` already takes arbitrary cfg overrides, so the
  pseudocount members join the same union as `min_cells` (at `sum.spec_override.parts[1].cfg_ens`); invalid c is
  rejected at parse time (only `tv2_sources.CS = 1/5/10/20/50/100` have precomputed `center_c<c>`). `drc` is a new
  builder block (`de.drc`) placed after the DE `pcr` step. AUDITED on spec_S16 (`audit_only`, A/B/C):
  drc formula bit-identical to the harness, `pred_dsum_{A,B,C}` BIT-IDENTICAL with and without it (PDS/MSE exactly
  invariant), override columns untouched; cens builds on all three contexts with 11 members, sums cosine 0.9938 vs
  base and gamma_sum 1.399 -> 1.28 (the budget match absorbing the ensemble's reduced spread, as documented for S11).
- **Standing conclusion:** PDS/MSE improvement is NOT finished (S15 moved MSE 0.749 -> 0.696 for +0.0106 on
  2026-10-02); what is finished is cheap RESHAPING of the delta we already have. Remaining: (1) change how the delta is
  BUILT -- 17 architectures already lose 0.04-0.22 via PDS, low prior; (2) **data** -- +0.003 per source context,
  +0.02..+0.04 for a same-line screen, which opens on 10/22. Spend the remaining weeks on (2).

## 2026-10-03 12-HOUR LOOP (user: "12시간동안 개선 루프 돌리면서 실행해줘"; FINDINGS §24) — pipeline VALIDATED, scoring saturated

- **FINAL-PHASE PIPELINE VALIDATED (§24.1, 24.7).** `tv2_final.py` ran the 10/22 command on the H1 panel as a new context with
  H1 + H1GW removed from sources / neighbour library / library PCs / S4, then scored with the real scorer. It first read -0.030
  vs the harness arm; attribution showed the rebuilt neighbour library and library PCs are fine, the pipeline's own sums are
  slightly BETTER (J +0.016) and its budget choice (gamma_sum 1.294) is the measured MSE optimum. The whole gap was the
  rehearsal feeding a DIFFERENT 12,000-cell H1 control draw than the panel's real data (profile cosine 0.99997, worth -0.023 PDS
  / +0.130 MSE). 10/22 rule: emit from the controls the organisers ship for that context; when rehearsing on a harness panel use
  that panel's own controls and `TV2_CELLS` to match its cells/target. New env hooks: `TV2_SRC_DIR`, `TV2_HOLDOUT_FAMS`, `TV2_CELLS`.
- **Closed this loop:** transform-side ensembling (pcr-lambda and centring averages are provably no-ops; the rest <= +0.002);
  the hubness operator again (helps Orion, -0.04..-0.14 MSE on RN); the kink FACTOR (floor-free panels disagree violently:
  RN_Jurkat +0.0075, OR_X **-0.181** at the same +7 % budget -> margin 0.8 stays); `nt`, `le`, plain-|d| DE ranking, K 2000.
  **Never price a budget / amplitude-scale arm offline** -- the probe said +0.077 where the scorer said -0.181 (§24.6).
- **Survivor:** `drc 0.8 / 1.1` (new token; per target the top 5 % of genes by |DE lfc| x 0.8 and the bottom 70 % x 1.1, bulk
  untouched so PDS / MSE are invariant): H1 +0.0018, OR_X +0.0004, RN_Jurkat +0.0005, 3/3 positive, H1 seed noise +0.0002.
  Worth ~+0.001 officially -> ride-along only. `K 5000` gains on H1 (+0.0022) but loses on OR_X (-0.0055): not adopted.
- **Data search exhausted (§24.10).** No public transcriptome-wide low-MOI screen exists for HeLa, CAL-33, HaCaT, BJ-hTERT,
  RPTEC/TERT1, HUVEC/TERT2 or any new stem line of use; the two new stem screens (WTC-11 GSE176038, KOLF2.2J GSE343369) cover
  0-7 of 300 official targets, i.e. chance. **X-Atlas/Pisces now has a public HF repo (`Xaira-Therapeutics/X-Atlas-Pisces`,
  CC-BY-NC-SA) with README + LICENSE only, unchanged since 2026-03-25** -- its Jurkat release would be a validation-phase
  (context A) asset, not a final-phase one. Watcher stays on.

## 2026-10-03: NEW BEST OFFICIAL — S16 0.24742 (entry lDWUfeHo4pEpR42pEc17, b24-92d6ca; rank 33 at retrieval)

- 3-hour multi-angle loop (FINDINGS §23.12): offline screens closed gene-level sign-agreement shrink, per-target neighbour blend,
  per-target S4 share, a higher budget margin (cliff), better agreement estimates and the DE-side amplitude (real scorer -0.006).
  The winner was a LEARNED amplitude signal: GBM on 8 harness panels predicting each target's true cosine from 8 test-time features
  (agreement and |delta| dominate); real scorer vs S15 positive on 8/8 runs (H1 +0.0014, OR_X +0.0093 / +0.0070, KOLF +0.0022,
  RN +0.0007..+0.0049).
- S16 = S15 with `sum.tamp` fed by `scripts/tv2_tlrn_prep.py` predictions (`reports/submission_2026-10-01/tlrn_official.npz`), beta
  0.75; spec `reports/submission_2026-10-01/spec_S16.json`. Official: +0.0014 over S15 (PDS +0.0016, MSE -0.004) = the H1 harness
  forecast. Builder now dumps `*_dbase_<ctx>.npy` (sum override before resid). Final-phase command unchanged (`--sum-tamp 0.5`).

## 2026-10-02: NEW BEST OFFICIAL — S15 0.24603 (rank 32; entry jvcHmRtUHN9TpQyD2rHb, b24-4b990c)

- S15 = S14fk2 + per-target SUM amplitude rank(cross-family agreement)^0.5 before the budget match (`sum.tamp`). +0.0106 over
  S14fk2: MSE 0.749 -> 0.696 (-0.053), PDS +0.002, reach +0.003, DE members otherwise flat. Forecast +0.005 (H1 harness +0.0065).
- **Recommended 10/22 command now adds `--sum-tamp 0.5`** (rehearsal 9 PASS): `tv2_final.py --controls-dir <D/E/F controls>
  --contexts D,E,F --out <OUT> --tag <TAG> --rebuild-sources --s4 --blend 0.3 --dual 1.0,3000,0.5,ag,cpc,ci5,lsb --kink 1.17,0.8
  --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g --sum-policy auto --sum-ens 5,10,15,30,50,100,200 --de-follow-sum
  --sum-resid 0.25 --s4-excl KOLF,FENG --de-denoise intended --sum-tamp 0.5 [--same-line 0.5]
  [--sum-cens 5,10,20,100] [--de-drc 0.8,1.1] --build`  (the two bracketed ride-alongs are new 2026-10-03, see FINDINGS 25.13); `find <OUT> -xtype l` empty.
- Follow-ups CLOSED (FINDINGS §23.11, user "다음 후보 진행해줘"): beta x centring x K562-weight grid (offline, 8 panels) -> S15 is the
  best cell on the non-Orion mean; real scorer vs S15 for the H1-favoured cells: ct1.0+ta0.5 mean(H1, OR_X) +0.0012, ct1.0+ta1.0
  -0.0001, ta0.75 +0.0026, ta1.0 -0.0010 (all RN negative) -> none reaches +0.005; S15 stays the final configuration.

## 2026-10-01 SUM-SIDE AGREEMENT AMPLITUDE (user "Sum쪽을 개선할수없어?"; FINDINGS §23.7-23.8)

- Leaderboard inference (§23.7): the top tier is pseudobulk transfer (described deep models max 0.178); ~46 % of their gains
  since 09-18 are MSE credit, with PDS-flat MSE jumps of -0.08..-0.16.
- New lever: per-target SUM amplitude a_t = rank(cross-family agreement_t)^0.5 before the budget match (PDS-invariant).
  Real scorer, S14fk structure: H1 +0.0065, OR_X +0.0215 / +0.0206 (2 seeds), RN_Jurkat +0.0045, RN_HepG2 -0.0045, RN_RPE1
  +0.0025, RN_K562 +0.0047, KOLF +0.0013. R21 composite +0.0054 (>= 0.005) but the pre-registered per-panel floor (each RN >=
  -0.003) fails on HepG2 -> gate FAIL by one criterion. Ported default-off: `scripts/tv2_tamp_prep.py`, builder `sum.tamp`,
  `tv2_final.py --sum-tamp 0.5`; official agreement `reports/submission_2026-10-01/tamp_official.npz`.
- **S15 BUILT, NOT UPLOADED** (user: "S15 패키지 빌드만, 업로드는 나중에"): spec `reports/submission_2026-10-01/spec_S15.json` =
  S14fk2 + `sum.tamp` beta 0.5. Package `reports/submission_2026-09-23/S15.vcc`, sha256 4b990cc8ebf56951149845ee8bd351e729c62916068a222e44c1ed9081d6b1c2,
  vcc prep passed. Audit vs S14fk2: dsum = S14fk2 x a_t to 1.7e-7, DE / main delta / overrides bit-identical, ew dropped on A/B/C as
  in S14fk2 (a first build had ew flipped back on by the shrunk sums -- builder fixed: the ew decision now uses the pre-tamp sums).
  Expected official ~ +0.005. Upload only with the user's go-ahead (from vcc2026_data): `bash ../scripts/run_detached.sh
  aivc-tv2-upload-s15 reports/submission_2026-09-23 /usr/bin/bash scripts/tv2_upload.sh S15 "S14fk2 + agreement amplitude"`.

## 2026-10-01 SPIKE "크게 개선할 수 있는 방향을 탐색해서 테스트해줘" (FINDINGS §23; Codex R24) — no big lever left

- **Leaderboard 10-01:** #10 = 0.2744, we are #35 at 0.2354; gap by member credit MSE +0.019 (50 %), PDS +0.012, reach +0.006,
  NMAE +0.004, fid -0.002. #5 [another team] has PDS 0.841 (= ours) with MSE 0.592: the top wins on strong-target DIRECTION. No team in the
  top 20 publishes a method description any more. Fact check of §1-22 (subagent): ensemble axes beyond min_cells dilute; untested =
  PC-removal-k ensemble axis, cell bootstrap, TF-regulon priors (75/300 targets are TFs, 45 of them ZNFs) -- all small.
- **P1 KILLED (§23.1):** GSE249595 (Jurkat, high-MOI genome-scale CRISPRi, 298/300 official targets, 827k cells) re-estimated by
  ridge regression over each cell's called guides (`scripts/ml_jurkat_lr.py`, pre-registered gate in `jurkat_lr/PREREG.md`):
  split-half per-target cosine +0.003 (official targets, 412 cells / target) and +0.008 (RN_Jurkat essential targets); LR =
  marginal estimator (cosine 0.98); barcode / feature alignment audited OK. Cells are shallow (median 2,900 UMI) and the knockdown
  weak (~40-60 %). Context A has NO usable same-line data; never retry this dataset (memory vcc2026-jurkat-gse249595-is-noise).
- **P2: dual token `lcp0.25` measured on the real scorer** (DE-side re-allocation weakened on high-expression genes; sums pinned):
  OR_X +0.0108 / +0.0093 (2 seeds), H1 +0.0015, RN_Jurkat +0.0055, RN_HepG2 +0.0073, RN_RPE1 -0.0027, KOLF +0.0008; gate
  composite +0.0026 < 0.005; lcp0.5 costs MSE +0.049 on OR_X (scatter below the #247 cap). Mechanism (§23.5): the cpc pin is
  unattainable for ~10 hyper-expressed MT genes when a down-gene is concentrated into the top-10 % cells; OR_X's NMAE is scored
  on 43 targets with ~24 real-significant genes each, so SLIRP alone (-0.775) makes the panel gain. Official budget +0.001..+0.002.
  **Ported as a DORMANT option** (`tv2_emit.emit_dual_iter(lcp=)`, spec `dual.lcp`, `tv2_final.py --dual ...,lcp<eta>`; parity:
  lcp 0 = old module bit-exact; rehearsal 8 `scripts/tv2_final_sim8.sh` PASS). **OFF by default; the 10/22 command is unchanged.**
- **Harness lesson:** on the shallow panels (OR 163 cells, RN 53-96) the NMAE member rests on a handful of targets with 10-30
  real-significant genes; one target can move the panel NMAE by -0.018 and the emission seed alone by -0.012. Only H1's NMAE
  behaves like the official one. Discount OR_X / RN NMAE deltas.
- **H1-TRUTH-TRAINED FUSION KILLED (user "가보자", FINDINGS §23.6, design `docs/superpowers/specs/2026-10-01-h1-truth-fusion-design.md`,
  Codex R25):** the leaderboard hint ([another team]: "donor-attention net fit on all 300 H1-2025 targets") was tested as the sixth
  source-weighting learner, now with the official-pipeline H1 truth as the training label (`scripts/ml_h1fuse.py` assembly +
  exact per-fold full-panel OOF evaluation, `scripts/ml_h1fuse_train.py` torch mixer, masked softmax over 9 donors). OOF full
  recipe: variant A (MSE + own cosine) PDS -0.0098 / MSE -0.0346, J +0.0145 KILL; variant D (MSE + in-batch InfoNCE) PDS -0.0085 /
  MSE -0.0299, J +0.0125 KILL; per-target NNLS oracle J +0.104. Own-truth cosine rose 0.096 -> 0.149 while PDS fell: the learned
  mix captures the truths' shared component (hubness 0.0034 -> 0.0067, truth 0.037). Same PDS-vs-MSE frontier as all five before.
  Untested part of [another team]'s recipe: the call-set-dependent EMISSION, not the fusion.
- **Codex R24 verdict (`debate_loop/codex_round24.md`):** no in-repo modelling direction with a defensible >= +0.01 remains;
  effort -> final-phase same-line data acquisition with a frozen S14fk2 fallback; each new source = independent same-target
  reproducibility first, then one held-out same-line full-recipe comparison at replacement weight 0.5.

## 2026-10-01: NEW BEST OFFICIAL — S14fk2 0.23543 (rank 35; entry gJHB2FvYxsWSnoTCp9zk, b24-54c07a)

- S14fk with the de.src links fixed (the INTENDED DE side): +0.0027 over S14fk (NMAE -0.0014, fid +0.0014, reach +0.0057,
  jac +0.0010; PDS / MSE identical). Confirms `--de-denoise intended` for the final phase on the official scorer.
- DE calibration (FINDINGS §22.9-22.10): no single harness panel predicts official DE members (NMAE sign: OR_X 2/3,
  majority panels 1/3); the R21 composite 0.5 x (0.75 OR + 0.25 RN) predicted S14fk (+0.0065 vs +0.0060) and S14fk2
  (+0.003 vs +0.0027). Decide on that composite, not on one panel's NMAE.
- Disk: freed 55 -> 181 GB (closed caches, redundant pred.h5ad). `src_ablation/` (57 GB) still to delete by the user (the
  safety check blocked my rm because it was the shell cwd); `_trash_2026-09-21/` (290 GB) kept (user did not choose it).

## LOOP 2026-09-30 06:00-18:00Z (user: "B부터 진행해줘 12시간동안 자동으로 돌리면서 판단하고 분석하고 테스트해서 알려줘"; FINDINGS §21)

- **NEW BEST OFFICIAL — S14fk: 0.23273, rank 34** (entry ckxrFX2RigxKeIvjbuUL, b24-5cbcbf; user go-ahead "업로드 승인할게",
  uploaded 2026-09-30 18:52Z; +0.00597 over S12, forecast +0.0057). Composition (composite): NMAE +0.0037 (0.9179 -> 0.9047),
  PDS +0.0008, MSE +0.0006, reach +0.0006, fid +0.0003, jac -0.0001. The sum-side S4 change transferred at ~1/4 of the OR
  harness; the DE side over-delivered (first time). Final phase: the S14fk structure is now officially supported -- mirror it
  (`--s4-excl KOLF,FENG` exists; the K562-denoised DE side still has to be ported into tv2_final.py).
- **BUG FOUND IN THE UPLOADED S14fk (FINDINGS §22.6):** `denoise_src/official_dnK562` had all 26 symlinks broken (relative
  paths), so the uploaded DE-side fusion saw denoised K562 ONLY. The harness then scored that "shipped" DE side against the
  intended one (§22.8): intended better on 5/6 panels, by +0.0057 on average (OR +0.0103, RN_HepG2 +0.0111, H1 +0.0070). The
  official 0.23273 was reached WITH the weaker DE side. Fixes: `ml_denoise_apply.py` resolves paths and asserts no broken
  links; the builder aborts on a broken de.src link. Corrected package **S14fk2** (`spec_S14fk2.json`, de.src ->
  `official_dnK562_fixed`) is being built (unit aivc-build-s14fk2-20260930) and is NOT uploaded. Expected gain over S14fk is
  ~+0.002..+0.003, below the +0.005 bar alone, but it is the only official check of which DE side the final phase carries
  -> the user decided NOT to upload it (2026-09-30 ~23:00Z); the package is kept: `reports/submission_2026-09-23/S14fk2.vcc`,
  sha256 54c07ad2846f7d5bcbeb632741beaba13fc95038221b8e100bd9cc97d7182259, vcc prep passed, dumps bit-identical to its audit.
- **FINAL-PHASE COMMAND (10/22), S14fk structure, ported and rehearsed (rehearsal 7 PASS, §22.7):**
  `tv2_final.py --controls-dir <D/E/F controls> --contexts D,E,F --out <OUT> --tag <TAG> --rebuild-sources --s4 --blend 0.3
  --dual 1.0,3000,0.5,ag,cpc,ci5,lsb --kink 1.17,0.8 --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g --sum-policy auto
  --sum-ens 5,10,15,30,50,100,200 --de-follow-sum --sum-resid 0.25 --s4-excl KOLF,FENG --de-denoise intended
  [--same-line 0.5] --build`. The new options default off; the old path is bit-exact.
- **Closed after S14fk (§22):** ② removing K562-specific programs (RN -0.015..-0.064); ① learned source-mixing gate pilot
  (PDS down on 6/6, J < 0 on 5/6 -- the 5th failed source-weighting learner; per-target mixing headroom is not learnable from
  test-time features). Validation top-10 judged out of reach by modelling; effort moved to the final phase (user's choice).

- **NEW SUM-SIDE LEVER: S4 without the FENG (HipSci iPSC) biology for somatic contexts** (§21.9; memory
  vcc2026-s4-feng-exclusion). S4 averages biologies at EQUAL weight, so a stem biology weighs far more there than in the
  fusion (FENG 0.25). Real scorer (agc, S12 structure): OR_HCT116 +0.0092 / +0.0117 (2 seeds; PDS +0.0073, MSE -0.030),
  RN_HepG2 +0.0018, RN_RPE1 +0.0031, RN_Jurkat -0.0008, RN_K562 +0.0003. Codex R21 pre-registered gate:
  0.5 x (0.75 OR + 0.25 RN) = +0.0041 < 0.005 -> S14f ALONE FAILS narrowly. Inputs: `scripts/tv2_s4_excl_final.py` ->
  `reports/submission_2026-09-30/S4noKF_ctx{A,B,C}`; specs `spec_S14f.json`, `spec_S14fk.json` (+ K562-denoised DE side),
  both audited vs S12 (main delta / overrides bit-identical). `tv2_final.py --s4-excl KOLF,FENG` (default KOLF) is dormant; rehearsal 6 verified it routes S4noFK to somatic contexts only.
- **S13k (K562 program-denoised DE side) is 6/6 positive on the builder-exact agc check** (H1 +0.0061, OR +0.0067,
  RN_Jurkat +0.0007, RN_HepG2 +0.0057, RN_RPE1 +0.0036, KOLF +0.0028; mean +0.0043 harness, ~+0.0015 official at the
  1/3 DE-side transfer) -> ride-along only.
- **BUNDLE S14fk (S14f + K562-denoised DE side) PASSES both gates** (FINDINGS §21.11): OR +0.0144 / +0.0178 (2 seeds),
  RN_HepG2 +0.0092, RN_RPE1 +0.0076, RN_Jurkat -0.0013, RN_K562 +0.0003; R21 discounted +0.0065, R22 forecast (DE increment
  at 1/3) +0.0057. Budget ~+0.005..+0.006 official (Codex range -0.003..+0.013). Package `reports/submission_2026-09-23/S14fk.vcc`
  BUILT 14:28Z, sha256 5cbcbf8230a58331a8fb7e3fe041a75e4010deeccf3c66f428a549f0cd9d1f78, vcc prep passed, dumps bit-identical to
  the audit. NOT uploaded: needs the user's go-ahead (from vcc2026_data: `bash ../scripts/run_detached.sh aivc-tv2-upload-s14fk
  reports/submission_2026-09-23 /usr/bin/bash scripts/tv2_upload.sh S14fk "S12 + S4 no-FENG + K562-denoised DE side"`;
  only the opaque b24-<sha6> name is sent).
  Robustness (§21.12): OR_HCT116X (HCT116 + HEK293T held out) base 0.3530 -> S14fk 0.3634 (+0.0104; PDS +0.0101).
  Caveat: RN_Jurkat is reproducibly -0.0014 (2 seeds) -- context A is where it is least likely to help.
  Do not change the final-phase default (no-KOLF) on this validation upload alone (Codex R22).
- **Closed today (measured):** B (destination-feature strength prediction; even the oracle re-rank is negative on 3/7
  panels); program-space denoising on the sum side (hubness) and on the DE side for all sources (dnALL); DE-side library
  projection (dpj: H1 -0.006, OR -0.015); Cfg.xw per-gene source-expression weighting (7-panel mean J 0.000); cis
  modelling for weak targets (their reproducible signal is ~88 % trans); FENG also out of the fusion (worse); CD4 out of
  S4 (OR offline MSE +0.08, real scorer in §21.11).
- **Leaderboard 09-30:** #5 = 0.274 (+0.047 over us); top teams were already 0.22-0.33 on 09-18 and then climbed at our rate;
  the 0.26 cluster looks like one group on several accounts; no public method >= 0.25 exists (subagent survey §21.8).
- C (single-cell learned model) NOT started: Codex R20 and the evidence say it does not address PDS.

## LOOP 2026-09-29 17:00-21:00Z (user: "가장 점수 상승폭이 높을 것으로 예상되는 방향 순서대로 테스트 ... 4시간 루프")

- **NEW BEST OFFICIAL — S12: 0.22676, rank 38** (entry 6NHNPxARhFMzvsG6RC8a, b24-61dd00; +0.0008 over S11). PDS / MSE identical
  to S11 (sums unchanged); DE side reach +0.0053 / fid +0.0023 / jac +0.0013 but NMAE -0.0042 (scaled). Harness said +0.0024
  -> a third transferred (DE-side harness gains over-predict again). The job sat 3 h in `launching` before scoring normally.
- **SAME-LINE DATA IS THE LEVER, AND 70 % OF IT IS ON THE DE SIDE** (FINDINGS §20.20; real scorer, H1 panel, H1GW =
  GSE295214 other-lab H1 screen covering 131/300 targets; harness tokens `sl<SRC>x<w>` / `sls<SRC>x<w>` = builder `same_line`):
  base 0.3029 -> sums only 0.3249 (+0.022) -> **sums + DE side 0.3772 (+0.074) at w 1, 0.3798 (+0.077) at w 0.5**
  (reach +0.13, NMAE -0.05, fid +0.02). Same-line STRENGTH is no better than ours -> keep the per-target norm match.
  The builder-exact agreement formula (`agc`) keeps it: +0.0735. RELATED lines (§20.23): KOLF <- H1GW +0.021 (w 0.5) /
  +0.024 (w 1); OR_HCT116 <- DLD1 +0.005 (w 0.5) / **-0.035 (w 1)**. Policy (Codex R17): `--same-line 0.5` recommended, every
  non-exact source capped at 0.5; budget +0.02..+0.04 per comparable destination officially (DE-side gains transfer uncertainly).
  => On 10/22: identify D/E/F, find a LOW-MOI same-line screen for each (audit provenance / state / multiplicity), ingest
  (`ml_newsrc_build.py`, now drops unknown-nperts cells), register in `tv2_model.FAMILY` + `tv2_final.SAME_LINE_SRC`.
  Validation lines: NO public same-line screen exists for Jurkat / HeLa / CAL-33 (§20.21); X-Atlas/Pisces (Jurkat, 10,872
  targets) is still unreleased -- watcher on; its release would be the one event worth a validation-phase rebuild.
- **FINAL PIPELINE = S12 STRUCTURE** (FINDINGS §20.19, memory vcc2026-final-pipeline-s12-port): `tv2_final.py ... --sum-policy auto
  --sum-ens 5,10,15,30,50,100,200 --de-follow-sum --sum-resid 0.25 --dual ...,ci5,lsb [--same-line W]`; residual / agree / PCs /
  S4 bit-exact on A/B/C; per-panel KNN differs by construction; CD4 x2 rule removed; x3 identity boost opt-in only.
  Builder audit (`audit_only`) with the shipped KNN: deltas / overrides / DE side bit-identical to S12, sums 1e-7, kink = S12.
  Rehearsal 5 PASS (§20.24): new panel, D/E/F somatic + G = H1 controls -> stem, H1GW fallback at 0.5, residual exclusions
  per context -- every branch ran. Still running (units aivc-loopq4 / aivc-loopq5): confirmation seeds @s1/@s2 for the H1
  same-line arms and the KOLF <- VCC2025-H1 source-emphasis control (read with `scripts/ml_sameline_compare.py`).

## MULTI-ANGLE BOTTLENECK REVIEW, 2026-09-29 (user: "크게 개선할수있는 부분을 찾고 파고드는게 중요해"; sc_hybrid FINDINGS §20)

- **Gap:** leader [top team] 0.397 (MSE 0.399, reach 0.42) vs us 0.226 (rank 40); MSE is 36 % of the gap, the four DE
  members half. The leader very likely has a proprietary same-lineage atlas (Billion Cell Atlas: cervix / lymphocyte /
  esophagus); benchmark against #2-#10 (0.26-0.33). No public method above 0.24 is written up anywhere.
- **Ceiling, scorer space (`ml_ceiling_v2.py`):** our sum side is at 80-91 % of the best per-target nonnegative mix of our
  sources; that mix reaches only ~45-60 % of the truth's noise ceiling on the strong responders that decide the MSE ratio.
  (A first version mixed spaces; Codex R12 caught it -- compare candidates in ONE space.)
- **What same-line data is worth (`ml_h1gw_value.py`):** GSE295214 (Weissman genome-scale H1 CRISPRi, built as
  `sources_new/H1GW.npz`, 2,747 targets) on the H1 panel: strong-target cosine 0.49 vs ours 0.34; blended into the 44 % of
  targets it covers, J +0.118 (~+0.02 composite). An order of magnitude above any within-data lever. FINAL-PHASE PRIORITY:
  once D/E/F are named, find a LOW-MOI same-line screen for each; the `tv2_final.py` same-line boost (x3 fusion weight) is too
  weak for that -- use near-full per-target replacement (design in FINDINGS 20.13). H1GW: stem-like final lines only; never a
  source on the H1 harness.
- **Closed today:** boundary ensembles; agreement-ranked DE side; LINCS L1000 shRNA; CD4 as Jurkat source; same-line
  other-target transform; strength-adaptive S4 share (H1-only); GSE249595 Jurkat as a same-line source (high-MOI noise:
  STIM vs UNSTIM cosine 0.001); emission sign loss (none); silent-target zeroing (~+0.001).
- **S12 CANDIDATE (not built, not uploaded):** the official DE side still used the pre-S10k main recipe (KOLF in, no
  ensemble). DE side = S11 sum-side recipe: real scorer, sum side identical, **6/6 panels positive**; the exact builder
  equivalent (ranking also from the new delta) mean +0.0024 (OR_HCT116 +0.0046, H1 +0.0014); Codex R14 audit: no defect;
  S12 BUILT `reports/submission_2026-09-23/S12.vcc` sha256 61dd0011...0161 (vcc prep passed), NOT uploaded; it also explains the S11 harness->official gap (~0.0014). Builder flag
  `de.from_sum_override` (regression: S11 main delta reproduced to 1.2e-7); `reports/submission_2026-09-27/spec_S12.json`.
  Expected official ~+0.002..+0.005 -- at the edge of the no-marginal-submission rule: needs the user's decision.

## SUBMISSION POLICY (user ruling 2026-09-27 night): do NOT spend slots on marginal gains

User: "아주 미미한 개선은 제출하지말자 신뢰성도없고, 의미도 크게 없어". Rule from now on:
- **Do not submit anything whose expected official gain is below ~+0.005**, standalone. Our harness->official transfer is
  unreliable at that scale (S10f: harness +0.006 -> official -0.0015; S10e -0.0046) and a slot spent there buys no
  information either.
- Verified CORRECTNESS fixes that cost nothing and cannot hurt (PDS/MSE provably invariant) may RIDE ALONG in the next
  substantial build - they just never justify a slot by themselves. `ci 5` (composition fixed point, H1 +0.0022 /
  KOLF +0.0004, NMAE better on both, PDS/MSE exactly unchanged) is in this category: include it, do not submit it alone.
- Effort goes where the measured stake is large: the final phase's coverage risk is ~0.015 composite per 20 % uncovered
  targets, i.e. an order of magnitude more than any remaining validation-side knob.

## Codex debate loop + final-phase coverage risk, 2026-09-27 night (details sc_hybrid FINDINGS §12-13, debate_loop/)

- **Loop**: `~/.local/bin/codex exec -s read-only ... resume <session>` (gpt-6-astra, medium). 7 rounds. Codex critiques /
  proposes, Claude implements / measures. State `reports/sc_hybrid_2026-09-25/debate_loop/{STATE.md,ROUNDS.md,codex_round*.md}`.
- **Produced one real gain**: the composition pin solved against the ORIGINAL control row totals while the scorer uses the
  FINAL ones. The existing but never-scored `ci` fixed point fixes it: H1 +0.0022, KOLF +0.0004 (NMAE better on both,
  PDS/MSE exactly unchanged). Residual error after ci5 is integer rounding (|log err| tracks 1/sqrt(counts/cell) over two
  decades) -> at its floor. **Candidate for the next build: add `ci 5` to the dual emission.**
- **Closed by measurement**: scored-geometry residual training; distal cis prior (oracle-bounded); override-consistent
  target (sd 6e-5); logL vs log(Lcur) (0/15,000 non-monotone); MODZ per-target source weighting; LPO-learned source weights.
- **Source transfer is saturated** (three independent confirmations): weak responders' cross-line transferable fraction
  0.12-0.23; our recipe beats perfect same-lab transfer from three lines; the per-target source-weighting ORACLE has
  +0.048/+0.130 PDS headroom but is unpredictable from test-time features. The remaining PDS gap to the top-10 is most
  consistent with more source CONTEXTS, not better use of ours.
- **FINAL-PHASE PRIORITY (new): coverage** (figure CORRECTED, see FINDINGS §13). The underlying stores cover 17,978 of
  18,533 genes (97 %); the earlier "13.9 %" was the size of the source files we had built, not the data. Over the whole
  axis: 74.2 % well covered (>= 3 families), 22.7 % limited, 3.1 % unmeasured, versus 97.7-100 % well covered on A/B/C
 . 20 % fully uncovered costs ~0.015 composite (`ml_uncovered_sim.py`); the LIMITED bucket's penalty is being measured.
  Status: (1) the
  final-phase fallback ALREADY uses concat DepMap+STRING (SA 0.25) for uncovered targets - the sweep's `depmap` baseline
  was the A/B/C path, where coverage is 100 %; the only open point is that Orion dislikes the STRING component;
  (2) DONE 2026-09-27: RN HepG2/Jurkat/RPE1 (896 targets, highest KD) can now enter the FUSION via
  `tv2_final.py --rn-sources [--rn-weight W]` (they were only in the neighbour library); S10k rebuilds bit-identically,
  so A/B/C is untouched; (3) still to build: Codex's sealed pseudo-final rehearsal (`debate_loop/codex_round7.md`).

## FINAL-PHASE RISK REPRICED — coverage is NOT the dominant risk; destination IDENTITY is (2026-09-29)

**Sealed rehearsal** (`scripts/ml_final_rehearsal.py`, Codex round-7 spec; `tv2_model.DROP_TARGETS` deletes a target's
measurements from EVERY source and the neighbour library, not merely from the fusion output):
  natural coverage after the destination's own family is excluded: H1 100 % / Orion 98 % / KOLF 99 % well covered (>= 3 families)
  delete 5 % of targets everywhere:  H1 -0.00086 | Orion +0.00019 | KOLF -0.00056 composite
  delete 20 % of targets everywhere: H1 -0.00171 | Orion -0.00004 | KOLF -0.00182 composite
That is **an order of magnitude smaller than the -0.015 recorded earlier**. The earlier figure came from
`ml_uncovered_sim.py`, which measured the DELETED TARGETS' own PDS (0.91 -> 0.71) rather than the panel composite, and
the panel-level loss is small because part 0 (S4) still predicts those targets from its 18,711-target training table --
which is the REAL final-phase condition, since a new 300-target panel's genes will also be in that table.
=> Stop buying coverage insurance. Residual caveat: S4 is not rebuilt under the exclusion in this rehearsal, so these
are lower bounds for a target absent from the table as well; that case is rare (97 % of expressed genes are covered).

**What IS large: the stem-vs-somatic source decision, per destination line.** Keeping KOLF on the SUM side is
+0.0079 PDS on H1 (stem) and -0.0126 on OR_HCT116 (somatic); shipping that decision for somatic destinations (S10k)
scored **+0.0061 officially**. The final phase has THREE new lines, so the decision is worth about that much three
times -- if it can be made from D/E/F's controls alone. `ml_identity_policy.py` validates a frozen control-based stem
score (pinned pluripotency vs somatic marker panels) against the measured dPDS of dropping each stem-like family.

**Which panel to trust** (FINDINGS §18.3): for SOURCE-FAMILY composition, use OR_HCT116 + the RN panels -- H1 is itself
a stem line and gives the wrong sign (the shipped no-KOLF change reads H1 J -0.011 / Orion J +0.037). Keep H1 primary
for emission and metric-profile questions. Orion's offline MSE is unreliable: read dPDS there.

## NEW BEST OFFICIAL — S11, 2026-09-28 23:14Z: 0.22596 (model b24-78b629, rank 41) — +0.00296 over S10k

S11 (spec `reports/submission_2026-09-27/spec_S11.json`) = S10k + three changes, of which only the FIRST paid:
1. **HYPERPARAMETER ENSEMBLE of the sum-side fusion** (`parts[1].cfg_ens`, min_cells 5/10/15/30/50/100/200 averaged with
   the shipped 20). Official: MSE 0.7647 -> 0.7529 (-0.0118), PDS 0.8205 -> 0.8217. This is the whole gain.
2. `ci5` composition fixed point (ported to `tv2_emit.emit_dual_iter`, `dual.ci`): official NMAE moved -0.0000165. NIL.
3. `lsb` depth-balanced emitted block (`dual.lsb`): official fid -0.0001, jac +0.0001, reach +0.0021. NIL/marginal.
=> H1 predicted +0.0044 (ens+ci5, 3 paired seeds sd 0.0003) and +0.0053..+0.0064 with lsb; official +0.00296. The
SUM-side part transferred; both DE-side fixes did NOT. Harness DE metrics are a poor guide to official DE metrics.

**Offline evaluator FIXED (important for every future screen).** Its MSE leg was biased by a CONSTANT +0.0103 for any
arm changing the delta's across-target spread, because the offline screen matches to the BASE budget while the real
scorer matches to the candidate's own (H1 `kink_base` 0.901 base -> 0.939 ensemble). Token `bmb` in `sch_emit.py`
matches to the non-ensembled recipe and the offline prediction then reproduces the real scorer exactly. Offline PDS was
accurate all along. The OFFICIAL builder already uses the main recipe delta as the budget reference, so it gets this for
free. Details FINDINGS §17.0.

## SOURCE CONTEXTS ARE STILL THE BINDING CONSTRAINT AND THEY ARE PRICED — 2026-09-29 (`ml_ctx_curve.py`)

PDS as a function of the NUMBER of source families (random subsets, both recipe parts and the neighbour library rebuilt):
  H1    k=1 0.7958 | 3 0.8576 | 5 0.8803 | 6 0.8884   marginal for the 6th family **+0.0081 PDS ~ +0.0030 composite**
  Orion k=1 0.6025 | 3 0.6939 | 5 0.7180 | 6 0.7272   marginal **+0.0092 PDS ~ +0.0034 composite**
  KOLF  k=1 0.5757 | 3 0.5896 | 5 0.6067 | 6 0.6048   marginal -0.0019 (KOLF is the stem outlier, as everywhere)
The curve is STILL RISING at k=6 and the log slope is nearly constant (0.053 at k=1->2, 0.045 at k=5->6), i.e.
PDS ~ 0.796 + 0.048 ln(k) on H1. **Every additional source context is worth ~ +0.003 composite.**
CONSEQUENCE, and it reverses a §15 decision: the seven public screens dismissed there as "small and not tissue-matched"
(THP-1 Yao, IGVF hPSC / cardiomyocyte / NPC, Tian-Kampmann iPSC neurons, Frangieh melanoma, Papalexi THP-1, Schmidt
primary T) are each worth ~+0.003 if the curve is about COUNT rather than VOLUME - `ml_fam_value.py` separates those two
by leave-one-family-out vs family size. X-Atlas/Pisces (16 contexts, Jurkat = context A) extrapolates to >= +0.02.

## NEW BEST OFFICIAL — S10k, 2026-09-28 01:18Z: 0.22300 (model b24-b73891, rank 50) — +0.0061 over S10c

S10k (spec `reports/submission_2026-09-27/spec_S10k.json`) = S10i sum side (centring 0.5 + guarded PC removal + budget
match) built from S4 / fusion / neighbour library WITHOUT the KOLF biology + no-embedding learned residual alpha 0.25; DE side
identical to S10c. Raw vs S10c: PDS 0.810 -> 0.8205, MSE 0.778 -> 0.7647, NMAE 0.916 -> 0.9157, fid 0.5234 -> 0.5235,
reach 0.278 -> 0.2748, jac 0.029 -> 0.029. Pre-upload estimate was +0.006 (S10i part +0.002-0.003 + S10k part +0.003-0.007):
matched. S10i itself was never uploaded (stopped by the user). One 09-28 UTC slot used; the second is free.
Receipt `reports/submission_2026-09-23/S10k.upload_receipt.json`. Debate loop with Codex started (`reports/sc_hybrid_2026-09-25/debate_loop/`;
Codex runs through `~/.local/bin/codex exec -s read-only` because the plugin's Node launcher is missing).

## S10k BUILD, 2026-09-27 16:20 — S10i upload STOPPED by the user (details sc_hybrid FINDINGS §10)

- User 16:13 PDT: "s10i제출 일단 정지하자. 우리가 했던걸로 개선한 모델 만들수있잖아". Unit `aivc-tv2-upload-sched-20260928`
  stopped before its 17:03 fire time (no receipt); `upload_queue_20260928.json` emptied (original saved as `.held_S10i_20260927`).
- **S10k** (`reports/submission_2026-09-27/spec_S10k.json`, package `reports/submission_2026-09-23/S10k.vcc`) = S10i with the
  SUM side from S4 / fusion / neighbour library without the KOLF biology + no-embedding residual alpha 0.25; DE side
  bit-identical to S10i. Harness vs S10i structure: somatic panels 6/6 positive (J +0.023..+0.097). NOT uploaded; needs the
  user's go-ahead. S10i stays built (`S10i.vcc`) as the fallback.
- BUILT 16:51 PDT: vcc prep passed (360,000 cells, 18,533 genes, targets verified, counts preserved); sha256
  `b73891cc86871bb609ce1224eec8bb014f47e3a4f1fa131c0b5f81a60c454382` (3.19 GB); pred.h5ad deleted. gamma_sum A/B/C
  0.931 / 0.933 / 0.928 (S10i 1.019 / 1.020 / 1.019); ew dropped on all three; pcr lambda 0.75; DE delta and overrides
  bit-identical to S10i; sums cos 0.934 vs S10i. To upload: put S10k in `upload_queue_20260928.json` and relaunch the
  scheduler (see the S10i entry below) or run `scripts/tv2_upload.sh S10k` - only after the user approves.

## Data-quality audit, 2026-09-27 15:00 (details sc_hybrid FINDINGS §9)

- No source family is harmful across panels (leave-one-dataset-out, real S4 kept). **KOLF inside S4** is the exception
  for somatic contexts: S4 refit without the KOLF biology (rows + empirical, `s4i_refit_excl.py`) is J > 0 on all 6
  somatic panels (RN +0.002..+0.015, Orion +0.015/+0.028), stem panels need it. Not built for A/B/C yet (needs a
  `tv2_s4_rebuild.py` exclusion option); candidate for one leaderboard A/B with the user's go-ahead.
- Weak-knockdown row filtering (`Cfg.min_kd`) hurts; keep all rows. RN essential screens are not in the official sources.

## Multi-angle bias review, 2026-09-27 afternoon (details sc_hybrid FINDINGS §8; S10i upload scheduled 17:03 PDT)

- **New harness panel OR_HCT116** (X-Atlas/Orion HCT116, somatic, 263 OFFICIAL targets x 163 cells, 12k NTC; OR_HCT116X also
  drops HEK293T). Official metric profile matches H1 far better than Orion -> H1 stays primary. S10i sum side at the official
  budget on Orion: +0.0023 vs S10c (ct +0.0038, pcr -0.0015: PDS +0.044 / MSE +0.080) -> no downside for tonight's S10i.
- **Closed with evidence**: gene-level sign prior (<= +0.001), control-cell covariance (fluctuation-response, chance),
  destination-line DepMap (< our norm), per-target strength models (oracle H1 +0.015 composite but true strength agrees
  across lines only at Spearman 0.29-0.60; our norm is at that ceiling; the old `so` oracle arm was invalid), NMAE amplitude.
- **fid** = sign accuracy of ALL our significant calls (n_pred 2,820 vs 145 real DE); low PDS + high fid = strong emission,
  weak biology. Top-3 teams' MSE implies ~1.8x our global cosine (MSE ~ 1 - 2 kappa rho sqrt(B/claim)).
- **S4 + discriminative (InfoNCE) residual** (user/ChatGPT proposal, `ml_dr_*.py`): with ESM-2 it HURTS H1 (lookup key);
  without embedding (target-agnostic 64-PC reshaping, target-disjoint training too) J > 0 on 7/7 panels x 3 seeds at alpha
  0.05-0.25, MSE only, ~+0.0015-0.003 composite on H1. Not wired into any build yet.
- **Open candidates**: stem-source down-weighting for somatic contexts (Orion PDS +0.007/+0.029; RN mixed; official mirror
  needs S4 rebuilt without stem biologies), the no-embedding residual, sum-side power rule (non-essential panels only).

## S4 improvement round — 2026-09-27 midday (user: "일단 S4모델 개선부터 시도해보자"); details sc_hybrid FINDINGS §7

- **Harness with the real S4.** The harness now takes S4 itself as part 0 (`sch_emit.delta_s4`):
  - saved leave-one-biology-out S4 for KOLF and the RN panels;
  - a leave-H1-out refit for H1.
  It reproduces the official S4 effect (H1 PDS +0.018, MSE −0.014).
- **Rejected:**
  - row z-shrinkage of the training table (H1 −0.002 to −0.025, RN negative);
  - knockdown-quality biology weights (no effect);
  - refit rank 16/64/128 (±0.001 on H1);
  - FiLM α 0.5/1 (KOLF −0.010/−0.032).
- **S4 share.** The sum-side S4 share of 0.5 is near-optimal. **S10j (0.75) is expected neutral to negative, so keep it
  unsubmitted.**
- **S10i expectation revised.** With the real S4, S10i's sum-side gain on H1 is +0.003 (vs +0.010 with the stand-in),
  because the S4 refit already centres its means. S10i is still expected positive, about +0.002–0.003 officially.

## 12-hour loop, 2026-09-26 22:46 → 2026-09-27 ~08:45 — COMPLETE (details sc_hybrid FINDINGS §6, LOOP_SUMMARY_2026-09-27.md)

**Next submissions (NOT uploaded; the user must approve each slot, 09-28 UTC opens 17:00 PDT).**
**S10i BUILT 2026-09-27 05:56**: `reports/submission_2026-09-23/S10i.vcc`, sha256
`3353638a2ec0398ea82f7ba417d45dc3a7c9af2e14c7bf87e1db495d9007cd6b`.
- vcc prep passed.
- Build decisions for A/B/C: guard f 0.037, lambda 0.75; ew dropped (kink with ew 0.82–0.87); kink_base 0.89/0.94/0.96;
  gamma_sum 1.019/1.020/1.019.
- The audit dump matches the formula to 1e-7.
**S10j BUILT 06:21**: `S10j.vcc`, sha256 `ca341c31ac0da3e32ced757fa818c0a111ba3d42227b942f3a2ad032bc01ce10`.
- vcc prep passed.
- gamma_sum 0.994 / 0.989 / 0.989.
- Its sums have cosine 0.971 with S10i's; the DE side is bit-identical, so S10i vs S10j isolates the sum-side S4 share.

Final-pipeline rehearsal 4 with the new flags passed (vcc prep rejected only the fake D/E/F labels, as before):
`tv2_final.py ... --kink 1.17,0.8 --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g`
(`scripts/tv2_final_sim4.sh`).
**S10i UPLOAD SCHEDULED (user go-ahead 2026-09-27 09:10 PDT, S10i only).**
- Unit `aivc-tv2-upload-sched-20260928` sleeps until 2026-09-28T00:03:00Z, then runs `tv2_upload.sh S10i`.
- Queue `upload_queue_20260928.json` now holds S10i only. Log: `reports/submission_2026-09-23/tv2-upload-sched-20260928.log`.
- Afterwards, check `S10i.upload_receipt.json` and the score.
- S10j is built but NOT queued (the user deferred it).

Earlier plan (superseded) — upload queue with both tags: `reports/submission_2026-09-23/upload_queue_20260928.json` (S10i, then S10j).
After the user's go-ahead:
`bash scripts/run_detached.sh aivc-tv2-upload-sched-20260928 reports/submission_2026-09-23 /usr/bin/bash
scripts/tv2_scheduled_upload.sh reports/submission_2026-09-23/upload_queue_20260928.json 2026-09-28T00:03:00Z`
(run from `vcc2026_data`).
- **S10i** (`spec_S10i_sumct_pcr_kink.json`) = S10c + sum-side centring 0.5 + GUARDED top-3 library-PC removal (lambda 0.75)
  + budget-matched kink.
  - The guard scales lambda by clip((0.25 − f)/0.10, 0, 1), where f is our own sum-side energy on those axes: official
    0.037, H1 0.098, KOLF 0.072, RN 0.32–0.40, so RN gets lambda 0.
  - Budget matching rescales transformed sums to S10c's own B_pred at max(1, kink); official s = 1.02.
  - The ew check drops ew on A/B/C (kink with ew 0.82–0.87 < 1).
  - Offline budget-matched J (≈ 6 × d composite): H1 +0.036, KOLF +0.062.
  - **Real scorer, H1: 0.2842 vs S10c 0.2741 (+0.0100; PDS +0.030, MSE +0.0045, reach −0.002, DE flat).**
    Seeds 1 and 2 give +0.0097 and +0.0100 (3/3).
    Centring + kink alone was +0.0048; without budget matching +0.0078.
  - **Real scorer, KOLF: 0.0580 vs 0.0423 (+0.0157; PDS +0.021, MSE −0.091, reach +0.012).**
    Seeds 1 and 2 give +0.0115 and +0.0128. Together with H1, that is 6/6 paired runs positive.
  - **Exact S10i configuration on the four RN panels:** K562 +0.018, HepG2 +0.084, Jurkat +0.079, RPE1 +0.044. With the
    guard at lambda 0 there, the configuration beats S10c on all six panels.
  - Official gain is expected smaller, because official sums carry the least generic-axis energy of any panel.
- **S10j** (`spec_S10j_sumct_pcr_kink_s4sum75.json`) = S10i with the SUM side mixed 0.75 S4 / 0.25 fusion (DE side
  unchanged). This tests the one lever the harness cannot see: S4 is worth +0.015 PDS officially (S10e). Official s = 0.989.
  - Harness proxy with the equal-weight stand-in for S4: 0.75 share gives H1 J −0.005, KOLF −0.012.
  - So S10j is an information experiment with an uncertain official sign; S10i is the expected gain.
- Real-scorer checks before upload: 028b (H1, lambda 0.75 guard ± budget matching), 028c (KOLF), 028d (DE-side
  `dpr3x0.75g`). Builder audit dump `S10i.pred_dsum_{A,B,C}.npy`.

**Closed this loop (do not retry):**
- Emitting more per-cell scatter: C_pred/C_real is 1.47 H1 and 1.18 KOLF, so the cap binds. 400 cells per target is fixed.
- NMAE magnitude profiles (NMAE is sign accuracy).
- Reach re-ordering: ~+0.002.
- Cis prior: already calibrated (≤1 kb-neighbour targets PDS 0.93).
- Sum-only blend 0.5: official PDS/MSE trade, net ~0.
- Agreement-weighted target amplitude: mixed.
- Per-source z-shrinkage:
  - Fusion-only: PDS down 5/6.
  - Library-only: negative.
  - Part-0 only: PDS +0.003..+0.016 on 6/6 but MSE worse, and officially part 0 is S4. Future idea: shrink S4's training
    rows.
- Knockdown-efficiency weighting: HEK293T median KD 0.50 but weight 0.5; small.
- No public recipe reaches PDS ≥ 0.84 (`vcc2026_data/reports/sc_hybrid_2026-09-25/public_repo_survey_2026-09-27.md`).

New code:
- `tv2_pc_prep.py`
- `tv2_emit.budget_B/budget_match`
- Builder options: `sum.pcr{guard}`, `dual.kink.budget_match`, `de.pcr`, dump `_dsum_`.
- `tv2_final` flags: `--sum-pcr k,l[g]`, `--de-pcr`; budget matching is automatic with sum-side options.
- `tv2_model.Cfg.zshrink`.
- Harness tokens: `pr<k>x<l>g`, `bm`, `dpr…`, `zs/zf/zl/zlp<k>`, `wm<FAM>x<m>`.

Screens: `sch_zs_screen.py`, `sch_zl_screen.py`, `sch_pr_sweep.py`, `sch_budget_match.py`, `sch_zp_floor.py`,
`sch_pc_energy.py`, `sch_wm_screen.py`, `sch_nmae_shape.py`, `sch_reach_order.py`, `sch_cis_diag.py`, `sch_tgt_amp.py`.

## Scorer-structure improvement loop, 2026-09-27 (harness; S10c remains the best official; details sc_hybrid FINDINGS §4)

Test-set criterion: follows from the scorer definition, test-time computable, wins on all six panels.
- **MSE budget kink** (`tv2_emit.kink_gamma`, builder `dual.kink` factor 1.17 margin 0.8): sums at max(1, 0.8 x kink);
  RN +0.022 (12/12), H1 0, KOLF -0.001 (its MSE ~1 is clamped, DE drift costs). Official A/B/C kink ~1.0-1.1 -> inert there.
- **Sum-side centring 0.5** (+kink): H1 +0.0045 (3/3, PDS +0.013, DE unchanged), RN +0.015..+0.069 (12/12), KOLF -0.003
  (0/3, DE drift at gamma 2.7). **Expression down-weighting** (1+cpm/c0)^-0.25 on the sums: offline LOPO best (6/6).
- DE drift at large sum amplitude: the per-cell composition pin (cpc) becomes unattainable; test `gdf` (DE follows gamma^a).
- Rejected: DE call counts per target (power-dependent), matched filter, gene reproducibility, per-gene calibration,
  low-rank, coverage energy, adaptive blend, SNR/agreement call-set membership. Held: top-PC removal (H1/KOLF up, RN MSE down).
- Tools: `sch_sumeval.py` exact offline MSE/PDS; `sch_lopo_screen*.py`; serial queue `sch_queue_runner.sh`
  (reports/sc_hybrid_2026-09-25/queue/*.job). Builder/tv2_final options: `--sum-center --sum-ew --sum-snr --kink`.

## S10f (09-27 slot 2): 0.21540, rank 57 — small harness knob gains do NOT transfer

S10c + DepMap blend 0.5 + nt 0.2 + K 5000. Harness predicted H1 +0.006 (3/3 seeds) -> official -0.0015: PDS 0.810 -> 0.800,
MSE 0.778 -> 0.763, jac 0.029 -> 0.031, NMAE 0.916 -> 0.919. Big structural changes transferred (S10b +0.035, S10c +0.025
matched the harness); sub-0.006 knob tuning did not. **S10c stays the best and the final-phase base** (keep --s4, blend 0.3,
dual 1.0,3000,0.5,ag,cpc). Adopt a knob only with a harness gain well above +-0.005 AND a mechanistic reason.

## S10e (09-27 slot 1): 0.21233, rank 57 — S4 is worth keeping (-0.0046 without it)

S10c with the S4 half (refit+FiLM deltas x3) replaced by the harness first part (equal-weight fold changes c=100 x3).
PDS raw 0.810 -> 0.795 (-0.015), MSE 0.778 -> 0.771, NMAE 0.916 -> 0.914, reach 0.278 -> 0.271. The learned S4 part carries
target-discrimination signal even under dual+cpc emission -> keep `--s4` in tv2_final.py. Harness cannot measure this
(S4 is trained on the panels' own lines). Iterative tuning rounds 1-12 and prior tests: `reports/sc_hybrid_2026-09-25/FINDINGS.md`.

## NEW BEST OFFICIAL — S10c, 2026-09-26: 0.21687 (entry 3L6CahZV3KUSYZc37Q4t, model b24-aa6a4c, rank 49)

S10b + **cpc** (spec `reports/submission_2026-09-23/spec_S10c_dual_cpc.json`, `dual.cpc: true`; `tv2_emit.emit_dual_iter(cpc=True)`,
bit-identical to harness `sch_emit.py ..._cpc`). After the dual re-allocation every expressed non-override gene column is
depth-tilted (x_c ~ w_c L_c^s, bisection) so the pinned group SUM is kept while the mean per-cell composition -- the
functional the scorer's DE LFC / NMAE uses -- equals the cellshift emission's expectation over all controls.
**+0.0247 over S10b, rank 107 -> 49.** NMAE raw 1.017 -> 0.916 (credit -0.025 -> +0.138); PDS/MSE identical to S10b
(sums pinned, verified exact); reach 0.292 -> 0.278, fid 0.5228 -> 0.5234. Harness: +0.020..+0.038 on 6/6 panels, H1 +0.023/+0.026
(2 seeds) -> official gain ~ harness gain. Why: pinned multiplicative sums make weak-effect genes' DE LFC sign-random (sum vs
per-cell-CPM functional mismatch); the cellshift target keeps a sign-consistent amplified shift on the DE side only
(`sch_nmae_diag.py`: H1 sign agreement by |pred| quintile 0.51..0.80 dual vs 0.53..0.75 cellshift). cpm (honest c0 x FC
target) does not help. Both 09-26 UTC slots used (S10b, S10c).

## NEW BEST OFFICIAL — S10b, 2026-09-26: 0.19215 (entry wRIB5uPULTJTaVJM0qZd, model b24-66bcd2, rank 107)

S9 recipe (S6 + DepMap neighbour blend 0.3, energy-matched) with DUAL-MOMENT emission (dual_mult gamma 1, K=3000,
du0.5 + ag). **+0.0349 over S8 (0.15732), rank 225 -> 107.** User ruled 2026-09-25: rules do not constrain emission
(memory vcc2026-emission-rules-no-constraint); S10b replaced S9 in the 09-26 00:03Z slot, S9 never uploaded.
Official vs S8: MSE ratio 6.73 -> 0.778 (first MSE credit, 0.221), PDS raw 0.764 -> 0.810, reach 0.271 -> 0.292, NMAE raw
0.923 -> 1.017 (credit +0.127 -> -0.025, the one loss). Local +0.064 over blend30_e -> ~55 % transferred; the harness
metric pattern (MSE/PDS up, NMAE down) held officially. Next lever: recover NMAE while keeping pinned sums.
Receipt `reports/submission_2026-09-23/S10b.upload_receipt.json`. Second 09-26 slot unused.
Single-cell hybrid study 2026-09-25 (`reports/sc_hybrid_2026-09-25/`, scripts `sch_*.py`): cell-cycle composition
correction and phase-stratified emission are ~0 on H1 (composition = 0.1-0.5 % of H1/KOLF response energy; oracle
positive 6/6, transferred version mixed); knockdown-responder Orion sources = 0 on H1 (escape 2 % on H1). Emission
sweep (`sch_emit.py`): pinning sums to the log1p-cellshift expectation loses (H1 0.202, KOLF -0.102); geometric
blend of cellshift amplitude into dual loses on H1 (b0.25 0.226 vs 0.251).

## NEW BEST OFFICIAL — S8, 2026-09-25: 0.15732 (entry aCQnYadgrqQiYYOuwcnX, model b24-9bbfa9, rank 225)

S6 recipe + DepMap co-essentiality neighbour blend (blend20_e: 0.8 x delta + 0.2 x neighbour delta, energy-matched;
`scripts/tv2_blend_prep.py` -> `knn_official.npz`, builder spec key `blend`). **+0.0041 over S6 (0.15321)**; local
harness predicted +0.0065 (18/18 paired) -> 63 % transferred, the best local->official ratio so far. Receipts
`reports/submission_2026-09-23/S8.*`. Both 2026-09-25 UTC slots used (S7a, S8).
10-hour sweep (FINDINGS §'10-hour sweep'): every learned change model (GBM/RF stacking, gain, learned family weights,
GBM/MLP residuals) loses or is nil on H1/KOLF; best is blend30_e (18/18 vs base, H1 +0.008, ~+0.003 over S8 on H1/KOLF).
combo (sim_soft + blend20) is +0.016 but its edge is RN same-lab. Next candidate submission: S6 + blend30_e.
DUAL-MOMENT EMISSION (fh_dual.py; builder mode dual_mult): pinned group sums at gamma 1.0 + top-3000 within-column
re-allocation for Wilcoxon calls: +0.061 vs blend30_e on 6/6 panels (H1 +0.029, KOLF +0.034). S10 package built
(reports/submission_2026-09-23/S10.vcc), NOT uploaded: rule status unconfirmed (#348). S9 scheduled 2026-09-26T00:03Z.
Refined: du0.5+ag (+0.003, 6/6 on H1/KOLF seeds) -> S10b (spec_S10b_dual_duag.json), NOT uploaded. tv2_final.py now has
--blend 0.3 (default; tv2_final_knn.py + uncovered-target fallback) and --dual 1.0,3000,0.5,ag; rehearsal built on final_sim.

## Final-phase harness + ML arms — COMPLETE 2026-09-24

`vcc2026_data/reports/final_harness_2026-09-24/FINDINGS.md`; code `scripts/fh_prep.py`, `fh_sources.py` (2,581-target
library incl. RN HepG2/Jurkat/RPE1), `fh_harness.py` (6 held-out panels H1/KOLF/RN x4, official preset, own family
excluded everywhere; idempotent). Base (final recipe, no S4) mean composite 0.226. Per-target amplitude -0.091 (0/6),
amp_shrink -0.047, expression multiplier -0.010 (3/6), low-rank lr64 -0.150, lr128_half -0.043: keep the recipe.
Unseen targets: DepMap co-essentiality KNN gives PDS 0.59-0.70 vs 0.51 for mean/zero (+0.20 composite, 6/6) -> use it
as the fallback for uncovered final-panel targets (not yet wired into tv2_final.py).
Follow-up: blend20_e (0.8 measured + 0.2 DepMap-neighbour delta, rescaled to base energy) beats base 18/18 paired
(3 seeds x 6 panels), mean +0.0065, min +0.004 (fid/reach/NMAE/Jaccard up, PDS -0.002). Not yet in tv2_final.py or any
submission. Catalog of public perturbation data: `reports/data_catalog_2026-09-24/CATALOG.md`.

## DE-set classifier — COMPLETE 2026-09-24 (user: "predict how many genes move and which")

`vcc2026_data/reports/de_classifier_2026-09-24/FINDINGS.md`; code `scripts/dec_sources.py` (source library with raw
counts for z-scores, 846 targets x 14 screens) + `scripts/dec_eval.py` (leave-one-context-out over six real DE tables:
H1, KOLF, RN K562/HepG2/Jurkat/RPE1; 6.5 min). With the true count, top-n precision is H1 0.207 (target-shuffled
0.186, expression-only 0.173), RN 0.22-0.34, KOLF 0.128 (expression-only 0.311: power, not biology). Count per target
is weakly predictable (Spearman H1 0.39, KOLF 0.16, RN 0.5-0.66). H1 sig_jaccard: oracle-k 0.151 vs S6 0.066 vs
automatic k 0.02-0.05. Not wired into any submission.

## Follow-up round 2026-09-24 (user: "run steps 1-4") — report §6-7 of transfer_v2 FINDINGS

**1. Scheduled uploads (09-25 UTC slots)**
- The detached unit `aivc-tv2-upload-sched-20260925` sleeps until 2026-09-25T00:03Z, then uploads the queue in
  `reports/submission_2026-09-23/upload_queue_20260925.json`. **The user dropped S7b on 2026-09-24; only S7a is queued**
  (S7a scored below the S6 recipe on H1, 0.1928 vs 0.2003; it is an explicit test of H1's under-prediction):
  - S7a = 0.75 fusion share;
  - S7b = S6 recipe with CD4 ×2 in context A only.
- Edit the queue JSON to swap packages before then. Check the receipts (`S7*.upload_receipt.json`) afterwards.
- **RESULT (uploaded 2026-09-25T00:03Z, entry djlCozt4S9QCdludHqWp, model b24-678ccd): S7a 0.14692 vs S6 0.15321 (-0.0063)**;
  PDS 0.756 vs 0.758, NMAE worse. H1 predicted the sign (0.1928 vs 0.2003). Keep fusion share 0.5 (tv2_final default).

**2. New data — neither helps, both excluded (kept as sources)**
- GSE249595 Jurkat CRISPRi is high-MOI: ~12 guides per cell; the marginal estimate gives on-target KD −0.26 to −0.46
  and transfers at chance.
- VIPerturb-seq K562 Flex (CC-BY) has KD −1.0 and ~47 cells per target; it is worse than Replogle K562 on the same
  targets.
- Data quality beats platform match.
- `_trash_2026-09-21` partially purged (432 GB; see its README); disk went from 99% to 75%.

**3. Final-phase pipeline `scripts/tv2_final.py`**
- Steps: HPA identity (`tv2_identity.py`: A/B/C → Jurkat/HeLa/CAL-33), a KD-quality-gated same-line boost, new-panel
  source rebuild, the S4 rebuild (`tv2_s4_rebuild.py`, bit-exact to S4 on A/B/C), then build/pack.
- The end-to-end simulation (`final_sim`, D/E/F) passed everything except `vcc prep`, which rejected the fake
  labels D/E/F. Use the real label names on 10/22.

**4. Emission decoupling (measured only)**
- Pinned group totals on H1: PDS 0.878 and MSE 0.955, but the DE members collapse (composite 0.098).
- Not submitted. Organiser-inquiry draft: `reports/submission_2026-09-23/ORGANIZER_INQUIRY_DRAFT.md`.

## NEW BEST OFFICIAL — S6, 2026-09-24: 0.15321 (entry 9ply7uUQW54qb6QZPeSt, model b24-004fe0)

Report: `vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md` (Korean).
Receipts: `reports/submission_2026-09-23/{S4gp,S4mix}.upload_receipt.json`. Both 2026-09-24 UTC slots used.

**What was submitted**
- **S5 0.14549** = S4 deltas ×3 with gated cellshift (genes <5 CPM in the context's controls are multiplicative, so no
  counts are invented) and a MULTIPLICATIVE CRISPRi promoter-neighbour prior (GENCODE v47 TSS ≤5 kb).
- **S6 0.15321** = 0.5×(S4 ×3) + 0.5×transfer_v2 fusion ×1.6 (true fold changes from raw sums, c=50 CPM, per-source
  genome-wide mean removed, weights K562 8 / HCT116 1 / HEK293T 0.5 / CD4 1 / KOLF 0.25 / FENG 0.25 / H1 2), same
  emission and prior. Its gain came from NMAE 0.949→0.926 and reach 0.259→0.265; PDS stayed ~0.761.

**Code** (run Python with `-P`; import cell_eval2 before adding scripts/ to sys.path, because `scripts/coverage.py`
shadows numba's coverage import)
- `scripts/tv2_sources.py`, `tv2_model.py`, `tv2_emit.py`: model and emission.
- `tv2_score.py`: H1 harness.
- `tv2_submit_build.py`, `tv2_submit_chain.sh`, `tv2_upload.sh`: build, pack, upload.

**New harness** = VCC2025 H1 (10x Flex, the official chemistry)
- 300 targets × 300 cells, 12k NTC, thinned to 20k UMI, official preset.
- It predicted S5's official gain exactly (+0.006 vs +0.0065) and under-predicted S6.

**Measured on H1**
- The new fusion lifts PDS 0.777→0.865 but loses reach/fid/NMAE.
- The log1p cell-shift transform amplifies 5–50 CPM genes, which is what gives the DE metrics their calls.
- Independent per-target draws kill PDS at small amplitude; the shared block is better.
- Never put a knockdown through the log1p shift (NMAE blows up).

**Not done / next**
- Next slot: raise the fusion share (0.75) and the A-context CD4 weight.
- Emission decoupling (dual-moment / pinned totals = #348 rebate residual) was deliberately not implemented; it
  needs a user decision and an organiser ruling.

## Leaderboard gap analysis — COMPLETE 2026-09-23 (read before any new modelling)

`vcc2026_data/reports/gap_analysis_2026-09-23/FINDINGS.md` (Korean; code + results in that folder).

**What top teams do.** The leaderboard API `description` field holds method notes from 428 teams
(`results/leaderboard_descriptions_2026-09-22.txt`).
- **Public recipe.** Most top-200 entries are pseudobulk transfer too. The open-source AtlasShift recipe
  (github.com/kaipengm2/Virtual-Cell-Challenge-2026) gets official **PDS 0.808 / MSE 0.943** from the same public
  sources we use; we are at PDS 0.750 / MSE 14.68. Its ingredients:
  - source weights K562 GWPS 2, H1 2, CD4 0.5;
  - per-source whole-screen mean removed;
  - fold changes multiplied onto the target's own controls;
  - amplitude 0.3 (bulk) / 0.6–1.0 (log2FC);
  - promoter-neighbour cap;
  - a renderer that pins group sums.
- **Context identities.** A/B/C are publicly identified as **Jurkat / HeLa / CAL-33**. The official data is
  10x Flex; the closest platforms are CD4 and VCC2025 H1.

**Our defects (measured).**
- **Reconstruction.** The per-cell log1p(CP10K) additive shift ×3 fabricates counts in genes the context does
  not express. 70% of S4's scorer-space energy sits in the 12,941 non-table genes (control median CP10K
  0.002–0.005). Those genes carry the raw multi-source mean ×3.
- **Energy.** S4's predicted energy is ~14× the true per-target signal (≈14, about equal to 400-cell noise).
- **Source weighting.** Equal/precision weights dilute K562 GWPS, the best single source on every proxy
  (e.g. HCT116: 0.686 vs mean 0.606). Its knockdown is also strongest: 85% vs HCT116 75% vs HEK293T 52%.
- **Promoter neighbours.** Not modelled. 60 targets have a neighbour TSS within 1 kb; in CD4 it moves by
  −0.45 to −0.56; adding the prior gives +0.01 to +0.02 PDS on proxies.
- **Shared 400-cell control block.** It adds a common offset and forfeits the #348 correction at low scale.

**RN harness is not the official task.**
- RN lines: essential-gene panels (0/300 official targets), ~8k genes, shared response 35–66% of signal.
- Official A/B/C: shared response ~1% (organiser baseline MSE 0.986–0.992). This explains the local→official
  reversals.
- Transfer matrix: same cell type, other state 0.97–0.99 PDS; same lab, other cell line 0.63; stem-cell sources
  transfer worst.

**Rules and next steps.**
- **Rules:** revised 2026-09-16 (published results may train models but must not be submitted or used to revise
  predictions). Ask the organisers whether direct mean transfer qualifies.
- **Final phase:** new lines and a new 300-target panel, no leaderboard; D/E/F controls arrive 10/22, deadline 11/5.
- **Next steps:** replace the harness with H1 + CD4 + Orion + KOLF (official targets, full axis, official preset),
  then fix reconstruction, weights and the promoter prior. Realistic target 0.17–0.20; top-10 needs new signal.

## Pseudobulk construction + within-condition uncertainty — NEGATIVE 2026-09-23 (stopped by user)

`vcc2026_data/reports/pb_uncertainty_2026-09-23/FINDINGS.md`. Cell-level robust aggregation (median,
winsor 99/99.5, trim 1/5, unweighted mean) on Feng-targeted/KOLF/Orion-eval rows (27 % of eval-target
rows; CD4/RN have no cells) and bounded gene-wise inverse-SE weighting within biologies: every arm
loses on all 5 lines locally (weight -0.0006 corr, robust -0.002, median -0.029, median zeroes
34-84 % of effects); officially weight -0.0007, nested robust -0.0004 (seed 0). Within-experiment z
adds only +0.0014 AUC to sign-confidence ranking. Direction closed. `recon_sweep.py` needs
`--threads 2` on this box. Caches in that folder (cells/, cd4/, arms2/) can be reused.

## NEW BEST OFFICIAL — S4, 2026-09-22: 0.138975 (entry df2PmtmQLtThacyTKu00, model b23-4b7c04)

Extended-table refit (CD4 donor x state + Feng targeted lines + Feng fitness arm nested inside
CD4/FENG, 200,163 rows) + 0.25 x FiLM r32 residual trained on all 38 gene contexts with
table-derived control profiles. **+0.001310 over the refit (0.137665)**; PDS +0.0061, NMAE +0.0019,
reach +0.0018, fidelity −0.0013, Jaccard −0.0006. Local five-line composite predicted +0.0247
(15/15 seeds): direction held, magnitude ~19x smaller. Report
`vcc2026_data/reports/submission_2026-09-22/` (FINDINGS.md, OOD_NOTE.md, receipts); build code
`scripts/ext_submit_build.py` + `ext_submit_chain.sh` + `ext_submit_upload.sh`. One 2026-09-22
slot used. **Leaderboard API** `https://virtualcellchallenge.org/api/leaderboard` gives every
entry's raw six metrics (`reports/leaderboard_2026-09-22.json`): our raw expr-MSE ratio 14.7 vs
top-30 0.58-0.79 (credit below 1.0), PDS 0.75 vs 0.85, gap decomposition MSE .072 / PDS .048 /
NMAE .040 / reach .033. `reports/recon_sweep_2026-09-22/FINDINGS.md`: delta scale is amplitude
restoration for the direction metrics, not calibration (least-squares slope 0.5, variance-matching
2, composite optimum 2.5 vs 3 = +0.004, line-dependent 2-3); cell-dispersion shrink is harmful;
S4 residual worth +0.014-0.019 over refit at every scale. Do NOT spend slots on scale.
Architecture screen COMPLETE (`reports/arch_screen_2026-09-22/FINDINGS.md`, 17 architectures x
K562 + HepG2, 104 official runs): every JOINT/fully-conditional model (CPA, biolord, joint MLP,
cond. flow, flow matching) loses 0.04-0.22 composite through PDS -> keep B(g) + residual. Residual
models differ by thousandths; candidates for multi-seed/multi-line validation: **bilinear** and
**MoE** (full FiLM safe third). Cross-attention worst residual on both lines; capacity buys nothing.

## Feng + CD4 gene-context experiment — COMPLETE 2026-09-22 (best local composite so far)

Report `vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md` (+ `CONTROLS_COMMON_BUG.md`,
`CRASH_RECOVERY_NOTES.md`, `coverage.json`, `final_table.json`). Added 12 CD4 donor x state
contexts (all targets, from the 42 GB author pseudobulk; stores `cd4_target_condition_donor_all_*`)
and 19 HipSci line contexts from Feng's TARGETED screen (`feng_targetedscreen_target_line`, 445
targets, 23 cells/pair); the broad Feng screens have 1-2 cells per line-target and stay pooled
(fitness arm added to FENG). Extended table 213k rows: `reports/ext_contexts_2026-09-21/table_v2`
(contexts nested in CD4/FENG) and `table_v1` (each context its own biology).
Official five-line composite vs refit 0.2755: refit_v2 **+0.0060** (12/15); refit_v1 -0.057; mean
transfer -0.04 / -0.12; **film_v2 = FiLM r32 residual (35 gene contexts, table-derived control
profiles) on refit_v2 at nested alpha: +0.0247, 15/15 paired seeds, 5/5 lines** (K562 +0.021,
RPE1 +0.036, HepG2 +0.039, Jurkat +0.022, KOLF +0.005); correct beats wrong context +0.052 and
zero +0.090 officially; gains are NMAE/reach/fid, PDS -0.014. Context-count scaling 6->35 is FLAT
(discrimination and PDS), so more contexts is not the lever; correct control profiles + finer data
in the refit are. NOT uploaded. Next: one leaderboard slot for film_v2 (build via
`dcp_submit_build.py`-style full-data fit: refit_v2 on all biologies + FiLM residual with A/B/C
profiles on the TABLE gene axis), then decide. Code: `scripts/ext_*.py`, `ingest_cd4_donor_all.py`,
`ingest_feng_targeted.py`; chain `ext_chain.sh` + `ext_followup.sh`. Never poll `is-active` on
run_detached units (RemainAfterExit); use SubState. `_trash_2026-09-21/` (718 GB, moved by
someone else at 08:03 09-21) holds S1-S3 packages and the real panels I restored; do not delete
without the user.

## BUG FOUND 2026-09-21: scrambled gene-context control profiles (read before citing dcp/cellfm/S2/S3)

`reports/cellfm_context_2026-09-20/controls{,_common}/*.h5ad` carry sorted gene names over
source-order counts for HEPG2, JURKAT, K562, RPE1 and KOLF (`cellfm_controls.py::write` reorders
columns only when duplicate symbols exist). Consumers: the cell-FM context experiment, the gene
side of the drug-context pretraining experiment (its gene context vectors were random keys, so
"gene contexts sit outside the drug manifold" and the gate/transfer conclusions are NOT
established), and the S2/S3 submission builds. Drug-side results, the refit, mean transfer and all
table-`controls`-based experiments are unaffected. Full account:
`vcc2026_data/reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`. Do not reuse those files
or the gene rows of `data/processed/drug_context_pretrain/contexts.npz`.


## Disk cleanup — 2026-09-21

The repo was 1.5 TB: 412 GB raw sources, ~1 TB of our own outputs, of which only ~26 GB were results.
Everything regenerable was moved (not deleted) to `_trash_2026-09-21/` (~770 GB) with `README.md`,
`MANIFEST.txt`, `PURGE.sh`, `RESTORE.sh`. Run `bash _trash_2026-09-21/PURGE.sh` to reclaim the space.
Staged: every submission's 17 GB predicted-cell `.h5ad` (all `.vcc`, SHA256 receipts and deltas kept);
per-line scorer inputs (`prediction.h5ad`, `real*.h5ad`, `baseline.h5ad`); per-fold dumps
(`predictions.npz`, `anchors.npy`, `training_anchors.npy`, `prepared_*.npz`, `targetfold*.npz`); the
data-audit `table_fullaxis*` tables (rebuild: `fullaxis_assemble.py`, `fullaxis_masked_assemble.py`);
old Sep-7 pipeline `vcc2026/{data/vcc_infer,data/st_input,work/st_main}` and `work/cellline_cache`
(rebuild: `scripts/build_cache.py`). KEPT because current scripts read them:
`reports/simple_residual_2026-09-16/data/` (the canonical 5,592-gene training table + Tahoe table),
`reports/context_model_2026-09-19/ctx_table/`, `reports/context_signal_2026-09-15/tahoe/` (downloaded
Tahoe shards), `cellfm_context_2026-09-20/controls*`, every `final/` and `final_decision/` dir, all
models, deltas, scores, receipts and FINDINGS. Rule going forward: delete scorer h5ads as soon as
their score CSV is written; never copy the training table into a report folder.

## S3 official submission of the scratch FiLM residual — LOST 2026-09-21

User authorised ("submit it"). Entry `5BTNATrQRIsC3XvRQayZ`, model b22-cbe57f, **0.137164 vs refit
0.137665 (−0.000501)**, PDS −0.0014, rank 241. Report `vcc2026_data/reports/submission_2026-09-21/`
(FINDINGS.md, build_receipt.json, upload_receipt.json). The +0.0180 local composite gain did not
transfer, the fourth such reversal. **Refit (entry 3CFK0NSTG9Fw6SfjD3IC, 0.137665) remains the best
official entry.** One daily slot used on 2026-09-21; check `upload_receipt.json` before any upload.
Build code `scripts/dcp_submit_build.py` reproduces the stored fold to 2.9e-5 and can build any
`dcp_gene_transfer` arm on all seven biologies for A/B/C.

## Drug-context pretraining experiment — COMPLETE 2026-09-21 (negative for context, positive for calibration)

Report `vcc2026_data/reports/drug_context_pretrain_2026-09-20/FINDINGS.md` (read §5-7 and the
recommendation). 192 official `cell-eval2 --preset vcc2026` runs, 3 seeds x 5 lines, PDS parity
1e-9 on every run. Five-line official composite: refit 0.2755; **scratch FiLM residual +0.0180**
(10/15 seeds, 3/5 lines); Tahoe+GSE-pretrained frozen trunk +0.0093 (14/15 seeds, 5/5 lines);
Tahoe-only frozen -0.0092; joint +0.0006; gated ensemble = refit. **RPE1 carries 80 % of every
gain** (refit's own composite there is 0.055 because its NMAE anchor is -1.14); without RPE1
scratch is +0.0043 and tgse_frozen +0.0025. Gains are NMAE/reach/fidelity, PDS flat or worse
(tgse_frozen -0.04). Correct context beats a real wrong context on the official composite 5/5
lines, 15/15 seeds (+0.026) - a first - but at alpha 1 it is below refit on 4/5 lines.
Drug side: discrimination grows 4x from 6 to 40 Tahoe contexts and correction is harmful below
~24, break-even at 40; gene contexts sit 2.3x outside the drug context manifold. Pretrained gene
heads all early-stop at step 100. Conclusion: context diversity IS the bottleneck for
discrimination, but drug contexts cannot substitute for gene contexts, and the composite gains
are magnitude calibration, not context. The top local arm was then submitted and lost officially (see above). PDS-only ranking reversed on the
full preset for the third time (tahoe_frozen +0.006 PDS -> -0.009 composite). No upload made or
authorised. Training curves are in `curves.json`, NOT W&B (protocol deviation).

Operations: scoring was killed three times by **systemd-oomd** (memory pressure > 50 % on the user
session): five terminal-parented scorers at 55 GiB took VS Code's pty host with them at 16:39 PDT
on 09-20; three and then two scorers under systemd hit the 40 GiB slice cap. One preset scorer
peaks ~16 GiB. **Run cell-eval2 preset scorers ONE at a time**, always via
`scripts/run_detached.sh`; pass parameters as arguments (systemd-run drops the caller env). After
any unexplained death run `journalctl -u systemd-oomd | grep Killed` first. Details in
`CRASH_RECOVERY_2026-09-20.md`; runner `scripts/dcp_preset_resume.sh N`, idempotent.

## Gene-embedding + projection-width arms — RUNNING 2026-09-19

Report `vcc2026_data/reports/esm_arms_2026-09-19/FINDINGS.md`. Units
`aivc-esm-arms-a-20260919` (done) and `aivc-esm-arms-b-20260919` (running the rest).
Code: `esm_arm_eval.py`, `esmc_embed.py`, `esm_arms_rest.sh`.

READ FIRST, it saves a wasted experiment: **gene embeddings cannot affect the refit baseline.**
`SourceBaseline.predict` overwrites the feature ridge with each target's empirical mean, and all
300 official / 279 KOLF evaluation targets are seen in every fold. Permuting or zeroing the
features changes the refit prediction by exactly 0.000. Embeddings only reach the b19 residual's
gene branch and the anchor ridge.

RESULT SO FAR — projection width is a clean null. ESM-2 15B five-line b19 mean: 64D 0.7958
(reproduces frozen b19 exactly), 128D 0.7940 (-0.0017, 8/15 paired), 256D 0.7946 (-0.0012, 6/15).
So the 64D bottleneck is NOT what limits the gene features.

WHY, measured: the adapter preserves a large shared direction. Median cosine between two
DIFFERENT genes is +0.65 for ESM-2 raw and +0.76 for ESM-C raw, and stays +0.60-0.67 after
projection at every width — a random projection preserves geometry so it cannot remove it.
Centring the embeddings first drops it to +0.02 (ESM-2) and +0.001 (ESM-C). The adapter
L2-normalises AFTER projecting, which is nonlinear, so the shared part is not an offset the
network can absorb in a bias. Centred arms were added for this reason; if ESM-C also comes back
flat, the conclusion is "this adapter cannot distinguish protein language models", not "the models
are equivalent".

ESM-C assets built: public ungated EvolutionaryScale checkpoints, isolated venv `.venv-esmc` (do
NOT install esm into ct_transfer, it would move torch). Sequences from UniProt human reviewed,
19,084/19,790 symbols matched. `processed/esmc/esmc_300m.pt` (960D) and `esmc_600m.pt` (1152D),
mean-pooled over residues, truncated at 2,046. Caveat: Arc's ESM-2 map pooling is undocumented, so
part of any ESM-2 vs ESM-C gap could be pooling.

## Data diagnosis — complete 2026-09-19

`vcc2026_data/reports/data_diagnosis_2026-09-19/FINDINGS.md`. The data is correctly processed and
not short of rows; it is short of transferable signal. Only 2.4% of rows have a true replicate
(Jiang only — KOLF's three batches screen DISJOINT genes, they are not replicates). Where
measurable, replicate reproducibility is r=0.209. The same CD4 cells 8 hours apart agree at only
0.270; different cell lines 0.03-0.17. Sources are not on a common scale (RN reports 1.6x larger
effects than Orion on identical targets).

## Gene subsets — complete 2026-09-19, negative

`reports/gene_subset_2026-09-19/FINDINGS.md`. Ranking genes by cross-biology ICC x between-target
variance and keeping the top 500/1k/2k/3k/5k NEVER beats all 5,592 (mean -0.029 to -0.001, 1-2/5
lines). The ranking IS real — it beats size-matched random by up to +0.063 — but cannot pay for
the discarded genes, and it does NOT help on KOLF.

## Wide gene axis — complete 2026-09-19, negative

`reports/data_audit_2026-09-18/FINDINGS.md`. The training axis had been cut to 5,592 of 18,533 by
intersecting with EVALUATION panels and the unused Tahoe drug cache. Rebuilt at 17,282 genes with
all biologies (masked): five-line refit 0.8016 vs frozen 0.8005, but only 2/5 lines and 4/15
paired seeds — FAILS the gate, and KOLF loses. Biologies matter ~10x more than genes.

## Data audit and wide-gene-axis rebuild — RUNNING 2026-09-18

User asked whether the data is insufficient or preprocessed wrong, then asked to fix it.
Report `vcc2026_data/reports/data_audit_2026-09-18/FINDINGS.md`. Read it before acting.

Two first-pass claims are RETRACTED there and must not be repeated: CP10K-vs-CP50K is NOT a
preprocessing bug (the reconstruction is the exact inverse of log1p(CP10K)), and cutting the x3
scale does NOT recover expression accuracy (our delta scores 0 on that metric at either scale).

The real defect: `simple_residual_assemble.common_axis` intersected the training gene axis with
the four RN EVALUATION panels, the Jiang panels and the unused Tahoe drug cache, cutting it from
17,847 to 5,592 of 18,533. Mean transfer was never affected. But the control arm shows biologies
matter ~10x more than genes: dropping the RN/Jiang rows costs -0.062 refit / -0.068 b19, and the
loss falls on the BASELINE not the neural part. KOLF, the line most like the unseen contexts, is
the one line where extra biologies gain nothing (+0.002), so the five-line mean flatters.

New scripts under `vcc2026_data/scripts`: `fullaxis_assemble.py`, `fullaxis_verify.py`,
`fullaxis_masked_assemble.py`, `fullaxis_masked_data.py`, `fullaxis_eval.py`,
`fullaxis_masked_baseline.py`, `score_delta_full.py`, `fullaxis_sequence.sh`.
Tables: `table_fullaxis` (84,752 x 17,282, verified bit-identical to frozen on shared genes) and
`table_fullaxis_masked` (103,295 x 17,282, all biologies, 9.1% NaN, needs mask-aware prepare).

Detached unit `aivc-fullaxis-sequence-20260918.service` waits for sweep 2, then runs the masked
baseline, then the wide-arm control, strictly serialized and gated on 24 GiB free. MEMORY IS THE
BINDING CONSTRAINT: the baseline fit peaks ~18 GiB and was OOM-killed once at 21.9 GiB; do not run
it beside `cell-eval2`. The narrow-arm control is complete and preserved.

No new baseline beats the frozen five-line refit mean of 0.8005 yet; that and b19 (0.7958) stand.
No official upload was made or authorized.

## B19 sweep 2 (weight decay x batch, seed ensembles, six-metric diagnostic) — RUNNING 2026-09-18

User approved arms A+B+C after the LR sweep came back flat, and requested launch.
Report `vcc2026_data/reports/b19_sweep2_2026-09-18/` (read `LAUNCH.md`, `status.json`,
`b19-sweep2-20260918.log`); plan `docs/superpowers/plans/2026-09-18-b19-sweep2-proposal.md`.
Detached unit `aivc-b19-sweep2-20260918.service` is the sole runner; never duplicate it.
Code: `vcc2026_data/scripts/b19_sweep2.py`, `b19_sweep2_full_metrics.py`, `b19_sweep2_final.py`,
`test_b19_sweep2.py`. Run tests as `ct_transfer/bin/pytest vcc2026_data/scripts/test_*.py`
from the repo root: from inside `scripts/` a local `coverage.py` shadows the coverage
package and breaks numba import under pytest. W&B group `b19-sweep2-2026-09-18`.
Not a re-run of LR/epochs/cosine (complete and negative, see below). Gates and the
conditional full-data refit are frozen in `protocol.json`; no official upload authorized.

## B19 LR / duration sweep — COMPLETE 2026-09-18 (gate failed)

Job exited 0 at 06:26 PDT; `status.json` stage=complete, pending=[]. Selected
0.794087 vs b19 reference 0.795777 (-0.001689), 1/5 line wins, 5/15 paired wins,
refit baseline 0.800538. All 30 LR x budget x schedule configs within noise; longer
training monotonically worse; cosine = constant. No promotion, no full-data refit,
no upload. Follow-up proposal (unapproved): `docs/superpowers/plans/2026-09-18-b19-sweep2-proposal.md`.
Original launch note follows for provenance.

User approved the proposed30-configuration sweep and requested launch. New report:
`vcc2026_data/reports/b19_lr_sweep_2026-09-18/`; plan:
`docs/superpowers/plans/2026-09-18-b19-lr-sweep.md`.
Detached systemd unit `aivc-b19-lr-sweep-20260918.service` is the sole runner.
Check its state, `status.json` and log before restarting; never duplicate it.
W&B project `wonminkim-uc-san-diego-health/aivc-vcc2026`, group
`b19-lr-sweep-2026-09-18`. First live run56b30bdd858c.
Grid: LR1e-4/3e-4/1e-3/3e-3/1e-2; sampled updates2500/10000/40000;
constant or cosine to1% base LR at each run's own budget. Fixed b19 architecture,
penalty.5,AdamW weight decay.03,batch256,alpha1. Source-only nested selection,
3-seed finalists/outer checks, original b19 and refit-baseline controls.
200 search trajectories/600 inner scores; <=120 repeat fits and<=40 outer/control
fits before deduplication.19 tests passed before launch and15 prepared scopes /
25 baseline caches verified. Historical scientific sources/artifacts remain read-only.
Full-six diagnostics and new full-data local refit only if frozen promotion gate
passes. No official upload and no100k extension authorized. Work is running;
do not label study complete until actual closeout status and artifacts agree.

## W&B logging required for future training — 2026-09-18

The user requests training logs in W&B from now on. Defaults are in
`wandb_config.json`; setup/use instructions are `scripts/WANDB.md`.
For future residual runners use `vcc2026_data/scripts/run_with_wandb.py` through
`scripts/run_detached.sh`, or directly use `context_forcing_wandb_model.fit` /
`b19_wandb_model.fit`. These new trainers log actual live losses and preserve
bit-exact trained weights in CPU regression checks. Other trainers must be
explicitly instrumented with `aivc_wandb.TrainingLog` before future training.
Log real validation results where the runner evaluates them; never label a
training probe as validation. Name/group each fit by experiment, fold and seed.
Keep original frozen trainer sources and historical results read-only. No
scientific retraining is authorized merely by this logging setup request.
Credentials are in the user's W&B login store, never in source, chat or logs.

## Crash recovery and future run protection — 2026-09-18

VS Code main process crashed on September 17 at00:54 PDT with a V8 ~4GB heap
limit error. The underlying allocation source is unproven; there was no matching
Linux OOM kill. Detached jobs finished at01:01. See
`vcc2026_data/reports/crash_recovery_2026-09-18/RECOVERY.md`.
For future long compute jobs use `scripts/run_detached.sh NAME LOG_DIR /absolute/program [args...]`
from the intended working directory. It uses user systemd, file logs, duplicate-name
protection and the shared `aivc-jobs.slice` (40GiB soft /48GiB hard RAM,2GiB swap,
8 CPU-equivalents total). Inspect status before any restart; no automatic retries.
Use absolute environment executables. Limits do not constrain GPU memory; run
GPU-heavy jobs serially. This protects jobs from editor closure, not reboot/logout.

## B19 convergence and separate-branch optimization — complete 2026-09-18

The user approved the proposed convergence/branch/capacity study. Work is in
`vcc2026_data/reports/b19_optimization_2026-09-16/`; plan is
`docs/superpowers/plans/2026-09-16-b19-optimization.md`. New scripts are
`b19_optimization_model.py`, `b19_optimization_data.py`, and
`b19_optimization_run.py`. Historical code/model/submission files are read-only.
All training and120 outer evaluations /20 context controls are complete. The
source-selected pipeline scores0.796162 versus fixed b19 recipe0.795777:
+0.000385,2/5 line wins,8/15 paired-seed wins. The frozen promotion gate fails;
`final/decision.json` retains b19, with no new full-data fit. Rank64 and the gated
hidden128 comparison lose; rank128 was not triggered. Longer training lowers
source MSE but hurts mean inner PDS. All10 full-six diagnostics and the final
report finished on September17; crash recovery on September18 reconciled the
stale status and handoff. Fresh verification passed:1292 checked files,20
preserved references,169 fits and209 checkpoints. Read `FINAL_REPORT.md`,
`DELIVERABLES.md`, `CLOSEOUT_VERIFICATION.json` and `research_status.json` in
the report directory. No fit, score job, refit, upload or poll remains pending.
This is a completed local study; no official upload was authorized or made.
B19 remains the preserved official best.

## Stronger residual sweep and official b19 — complete 2026-09-16

Read [the final report](vcc2026_data/reports/context_forcing_strength_2026-09-16/FINAL_REPORT.md),
[deliverable index](vcc2026_data/reports/context_forcing_strength_2026-09-16/DELIVERABLES.md),
and [status](vcc2026_data/reports/context_forcing_strength_2026-09-16/research_status.json).

- **b19 is published**, ID **pamfw9b1glTBvlXvEyQA**, official overall
  **0.13646991123621474**, rank193 at retrieval. This is the project's best
  recorded overall: +0.0072768684 / +5.63% versus preserved b16.
- Calibrated official PDS is **0.5431564640579106**, below b16's
  **0.554670505212238**. The overall gain is mainly NMAE. Same panel/anchor.
  This does not prove incremental context learning; baseline rank/data also differ.
- Alpha1/1.5/2/3/4 sweep finished: 120 source-inner checks, 90 three-seed outer
  outcomes and20 context controls. All five inner selections choose alpha1.
  Local five-line raw PDS: alpha0 0.800538, alpha1 0.794715, alpha1.5 0.790091.
  Five stronger outputs are scored but count-ineligible, three unscorable.
  All selected-alpha1/control outputs valid. Stronger corrections help KOLF
  but do not improve consistently across lines.
- New full-data checkpoint: standard small MLP,24,192 parameters, rank32 residual,
  context16, hidden64, alpha1, seed0,103,295 genetic rows/18,711 targets;
  exhaustive coverage then2,500 balanced batches. No drug rows. A/B/C residuals
  are nonzero and distinct. This is not the b16/b18 package and is not newly
  trained historical b16. Final model lives in the new report's `final/`.
- Package `submission/b19.vcc`, SHA256
  `ca8d60b1cd8f5f39b215cfde2cd4973f2fc2047967463384760ff0c125f92f47`.
  Official validator, physical reconstruction review, all-grid audit and eight
  targeted scientific tests pass. Thirteen prior references are hash-identical.
- One new entry uploaded and actual score retrieved. No public-score tuning,
  pending fit, upload or poll. Read the new `submission/receipt.json` before
  any future submission to avoid duplication. Original b16/b18 records and
  checkpoints remain preserved; new scores are in `submission/official_scores_new.csv`.

## Context forcing and training without empirical mean — completed 2026-09-16

The user's requested follow-up is complete. Read the
[final report](vcc2026_data/reports/context_forcing_2026-09-16/FINAL_REPORT.md),
[deliverable index](vcc2026_data/reports/context_forcing_2026-09-16/DELIVERABLES.md),
and [status](vcc2026_data/reports/context_forcing_2026-09-16/research_status.json).

- **158 distinct fits**: 110 source-inner pilot fits and 48 full fits, plus two
  explicitly labelled gene-only control aliases. All 210 full count-score records
  are complete. Five whole-line holdouts, three primary seeds, genetic data only.
- Tested weaker penalties, baseline-input removal/dropout, strict context
  interactions and contrast losses; also direct models with no empirical mean
  in their input, target offset or output. Residual alpha is fixed at one.
- Exact biology-balanced source MSE improves in all 30 primary full models:
  residual reductions 11.73–42.20%, direct 11.87–23.38%. Nonzero corrections do
  not establish transferable context. Previous b18 equalled b16 because its
  package deliberately reused b16; the prior models did not all learn zero.
- With identical predicted-effect gene support, five-line raw PDS is refit
  alpha-zero **0.786234**, residual **0.783225**, direct **0.574450**. Direct loses
  all 15 paired seeds. Full-axis residual **0.794715** loses on average to refit
  **0.800538**, winning 6/15 paired seeds. Full-axis KOLF improves slightly
  (**0.568947 vs 0.565788**) but remains below archived b16 (**0.596134**).
- Correct context beats separately trained no-context in only 2/5 residual
  lines and 2/4 active-context direct lines. K562 direct selected gene-only;
  its context-control ties are aliases. KOLF's residual gain does not survive
  attribution controls. No stable new context benefit or model promotion.
- Forty scientific tests, independent artifact/arithmetic review and exact-count
  reconstruction checks pass. All historical references and the 82 pre-acceleration
  receipts are preserved. Final report review passes with no material corrections.

These are local reconstructed PDS results, not new official scores or proof that
context learning is impossible. All full fits use every permitted genetic row,
excluding the destination biology across studies/stimuli; there is no new
all-training deployment model. No fit, scoring job, upload or score poll remains.
No official submission was made in this follow-up. The b16/b18 result below and
the historical best checkpoint remain unchanged.

## Simple pseudobulk residual workflow completed through official scoring — 2026-09-16

The resumed genetic-only versus genetic+Tahoe request is complete through actual
official validation and score retrieval. Read the
[completed report](vcc2026_data/reports/simple_residual_2026-09-16/independent_verification_1937/FINAL_REPORT.md),
[deliverable index](vcc2026_data/reports/simple_residual_2026-09-16/DELIVERABLES_COMPLETE.md),
and [current status](vcc2026_data/reports/simple_residual_2026-09-16/research_status.json).
This newer experiment supersedes the older statements that drug-transfer
comparisons or official reconstruction remain unexecuted.

- Linear, bilinear and small-MLP residuals, genetic-only/joint/pretraining arms,
  sharing ablations, matched drug controls, whole-line selection and count-based
  evaluation finished. No stable incremental context or drug benefit was found.
- The four-line primary mock selected a refitted rank-32 baseline at alpha zero:
  raw PDS 0.902890 versus archived b16 0.854647. The closer-panel KOLF check
  (279 official targets, 17,971 genes) instead favored b16 0.596134 versus
  refit 0.565788 in all three paired seeds. The separately frozen official
  fallback preserves the original primary decision and all results.
- The expanded-data challenger was fully fitted on 103,295 genetic rows over
  18,711 targets and archived in `final_decision/final/`, relative to the new
  report directory. It was not submitted. Do not describe historical b16 as
  newly trained on that expanded table. Wholly unseen-target evaluation remains
  a separately reported pilot, with performance near chance.
- The one new entry **b18**, ID **Ax90uDSJhyMEFBacq5WB**, is **published** at
  **0.12919304282260788**. Its package is byte-identical to b16 and its score
  equals b16 to numerical precision. This is baseline revalidation, not an
  improved model. The selected predictor is historical rank-128 b16, new
  context-residual alpha zero, scale 3, with no drug contribution.
- Full-six, random-null, disjoint-control and Jiang-repeat diagnostics finished;
  official validators and independent artifact/preservation checks passed.
  The b16 checkpoint, package and original submission ledger are preserved.

No fit, upload or score poll is pending. Before any future submission, inspect
`vcc2026_data/reports/simple_residual_2026-09-16/submission/receipt.json` to avoid
duplicating b18. The current best predictor remains preserved b16.

## Resumed original research prompt: Stage 10 completed — 2026-09-16

The bounded model runs below had finished, but the requested quantitative recovery
of repeatable cell-specific variation was missing. That analysis is now completed
in a separate supplement; all earlier reports and checkpoints remain unchanged.
Read [the current research summary](vcc2026_data/reports/context_signal_2026-09-15/resume_2026-09-16/context_signal_summary.md),
[deliverable index](vcc2026_data/reports/context_signal_2026-09-15/resume_2026-09-16/DELIVERABLES.md),
and [restart status](vcc2026_data/reports/context_signal_2026-09-15/resume_2026-09-16/research_status.json).

- **Tahoe has a small exploratory RNA context signal on the new repeat-based
  diagnostic.** Saved full-GP predictions improve latent source-to-destination
  contrast recovery by **1.445%** over matched no-context, descriptive interval
  **0.540–2.495%**, with **7 improving lines / 9 zero-strength ties**. This uses
  exactly 91 repeated conditions and replicate-specific source means. It is not
  total expression variance explained or an official score. No new Tahoe fit.
- **The original Tahoe PDS promotion gate remains failed.** Its line threshold
  and independent-repeat PDS evidence have not changed. Full-GP wrong/training-
  shuffle vectors were not archived for the new diagnostic. Do not equate the
  small new signal with robust context transfer or proof that context is absent.
- **Jiang: 92 separate-replicate small-kernel folds.** All fitting, feature inputs
  and inner residual-error selection use one biological replicate; the other is
  evaluation-only, with directions reversed. PCA2/HVG256 recovery is **−11.227% /
  −0.930%**, both worse than no correction in 6/6 line averages. Historical Jiang
  mean/denoiser/PDS-selected kernel scores remain unchanged.
- **LINCS: 64 separate-half folds.** Technical and disjoint-reagent recovery
  means are negative. These endpoints do not identify biological reproducibility.
  DepMap and the four-line Replogle/Nadig panel explicitly lack a complete
  biological-repeat noise decomposition. The latter has a 1,479-target K562
  screen-agreement check, not another independent context.
- **26 scientific tests and 16,128 historical Tahoe metric reconstructions pass**;
  actual response/control poisoning, installed scorer parity and independent
  reviews pass. All 13 original reference hashes and 21 key prior report hashes
  match; metadata for all 11,659 prior research files is unchanged.
- Drug-to-gene A–E, the mechanistic bridge, larger architectures and challenge
  reconstruction remain conditional and unexecuted. **b16 remains the recorded
  official best, 0.12919304.** No new download or official submission occurred.

All jobs in this continuation have finished. The older unconditional completion
wording below described model runs and omitted Stage 10; this supplement closes
that gap without replacing the historical evidence. The next eligible research
direction is independent confirmation of the small Tahoe residual signal with a
separately frozen analysis, before drug-to-gene pretraining or larger models.

## Context-signal workflow recovered and completed — 2026-09-15

The interrupted session had finished Tahoe preparation and leakage checks, but
had not started its main model evaluation. The requested recovery **completed**
that evaluation and the already prepared auxiliary DepMap fitness experiment.
Read [the combined research summary](vcc2026_data/reports/context_signal_2026-09-15/context_signal_summary.md)
and [the recovery record](vcc2026_data/reports/context_signal_2026-09-15/recovery_2026-09-15/RECOVERY.md).

- Tahoe: **16 lines, 150 drugs, 128 model-folds / 3,600 fits**; all outputs verified.
  Seven rank-16 primary families plus a separately frozen full-response GP.
  **No family passes the predefined context criteria.** Mean PDS 0.99895414,
  denoiser 0.99935123, full GP 0.99899049. Full GP wins only 5/16 lines and loses
  on the independent biological repeat on average. These are native pseudobulk
  scores with little discrimination headroom, not official scores or proof that
  context cannot be learned. Drug-to-gene A–E training remains deferred by its gate.
- The dose metadata has a documented exception: **149 drugs at 5 µM and Adagrasib
  at 0.05 µM**, consistently matched across lines. Neural penalties were fixed
  despite duplicate nominal candidates. Seed-matched shuffle diagnostics and
  reporting qualifications are saved in `recovery_2026-09-15/`.
- DepMap: **1,021 model rows / 984 biological groups, 512 genes, 15 model-folds**.
  Ridge reduces native fitness error by 4.52% and gene-specific error by 4.51%,
  with 5/5 outer-fold wins and passed predictor-level controls/poisoning checks.
  This is processed cancer knockout fitness; it does not establish RNA-response,
  unseen-tissue, or drug-to-gene transfer. Historical preparation files are retained.
- Earlier Jiang full-response kernel follow-ups show a small exploratory context
  lead; the inner-selected full-assay version scores 0.742486 versus mean 0.736738,
  still below denoiser 0.779015. The independently validated no-context MoE loses
  on broader genetic panels and was not promoted. See the combined summary.
- **b16 remains the best recorded official submission, 0.12919304.** All 13 frozen
  reference hashes match. No official reconstruction or upload occurred. The
  original DeltaNet objective mismatch remains a separate historical item.

These runs are finished; no training process or user input remains pending for
this bounded workflow. The completed report and recovery receipts supersede older
status/preparation statements below and in the preserved preliminary reports.

Written 2026-09-14. Everything below was measured in this repo unless marked as an assumption.
Numbers with no hedge are arithmetic on files you can re-open; anything uncertain says so.

## Eight-architecture context test — 2026-09-15

**Subsequent diagnostic review:** read
[the ML training diagnosis](vcc2026_data/reports/ml_training_diagnosis_2026-09-15.md).
Pairwise DeltaNet's main b16-residual objective and auxiliary raw source-to-target
difference objective have a confirmed mismatch of `b16 - mean_source_response`.
The discrepancy is nonzero on actual permitted-source IFNB training episodes.
Original scores and artifacts are unchanged; a corrected DeltaNet implementation
and rerun remain outstanding. This does not explain the other seven families.
The small shared-feature sweep is not a verdict that ML or these architecture
families cannot work. No new training or official upload occurred in this review.

The user's requested eight residual architectures are **implemented and fully evaluated**.
Read [the completed comparison](vcc2026_data/reports/context_architectures_2026-09-15/RESULTS.md)
and [full table](vcc2026_data/reports/context_architectures_2026-09-15/model_comparison.csv).
The official b16 checkpoint and original Jiang preprocessing/denoiser/scorer are unchanged.

- One interface in `scripts/context_architectures.py` and one nested runner in
  `scripts/run_context_architectures.py`; Pairwise DeltaNet, hypernetwork, conditional
  DeepSets, MoE, tensor, CNP, learned metric, GP all completed in the requested order.
- **184 outer model-folds / 3,300 training fits**, about 445 seconds for the main run.
  Same 23 Jiang folds, response ranks 8/16, training-only control PCA rank 2, inner-selected
  residual strengths 0/.25/.5/1, separate replicate labels, seeds 0/1/2. Every architecture
  has no-context, wrong/zero query, training-label/context-shuffle and 32 target-shuffle
  controls. Seeds for GP are identical by design, not independent robustness evidence.
- Fixed permitted-source denoiser PDS **0.77901524**, mean **0.73673816**. Correct-context
  PDS: Pairwise **0.77928236**, metric **0.77925908**, GP **0.77923649**, hypernetwork
  **0.77895653**, CNP **0.77865385**, MoE **0.77826601**, tensor **0.77820737**,
  DeepSets **0.77751312**. These are local CP10K pseudobulk scores, not official scores.
- **No reliable transferable context-learning success.** Pairwise loses to matched
  no-context (0.78006935); its three-seed mean 0.77889754 also falls below b16. Metric
  wins only 3/6 lines and its zero-query control is better. GP wins 4/6 lines but both
  training-label and training-context shuffles outperform it. Do not promote any of them.
- Limited next context experiments: **learned metric first, shared-kernel GP second**,
  explicitly as exploratory leads. Test perturbation-independent/fixed-kernel controls,
  more training permutations and stronger inner-fold regularization/uncertainty calibration
  before considering broader biological validation or official reconstruction. The no-context
  MoE's 0.78206266 is a separate denoising lead, not context-learning evidence.
  Its three experts are linear and its context-free routing weights are constant within
  a fitted model, so the residual module is exactly equivalent to a single linear map.
  Do not describe this result as evidence for nonlinear neural modeling or adaptive
  expert specialization; compare a single linear residual with matched tuning next.
  The independently selected no-context MoE wins pooled PDS in 5/6 Jiang line
  averages and improves both separate-replicate scores in 5/6 lines. Its exact
  selected settings still need seed repeats and broader-biological validation;
  original repeated no-context runs used the context-selected settings instead.
- Core scientific tests and an independent nested-selection poisoning test pass. Fresh
  full-length IFNB/A549 fits for all eight families preserve all 80 inner scores, selections,
  and unshrunk predictions exactly under held-out response poisoning. GP/CNP latent
  uncertainty is saved; it does not include b16 or omitted-gene uncertainty.

Paths in these bullets are relative to `vcc2026_data`. b16 remains the best recorded
official submission at **0.12919304**; no new reconstruction/upload was performed.
All runs above have finished; no user input or pending leaderboard result is required.

## Implemented follow-up — 2026-09-15

**b16 is the best recorded official submission: 0.12919304**, versus b8's preserved
0.12170879. b15 lost at 0.10671379. The b15 → b16 submission sequence is **complete**;
there is no pending login or b16 result, and b17 remains deferred. b8's old entry returns
`not_found`; its comparison score is the preserved official record.

The next context-aware model is now implemented and evaluated, not merely proposed.
See [the completed experimental report](vcc2026_data/reports/context_interaction_2026-09-15/RESULTS.md).

- Target-response components x untreated-expression components, signed ridge correction,
  with inner whole-line selection and source-only baseline predictions for training examples.
  All **23 Jiang outer folds completed**. Context model PDS **0.778350**, refitted denoiser
  **0.779015**, no-context residual **0.780934**. Context assignments shuffled: **0.778289**.
  Correct context information has not established an incremental benefit. Do not promote it.
- A 1.5x correction improves Jiang slightly and KOLF PDS across all four paired seeds
  (**0.257194 vs current 0.230554**), but loses in **all four broader held-out biologies**
  (**0.569491 vs current 0.610668**). The current correction remains the reference.
  No new global scale sweep or leaderboard upload was performed.
- `scripts/harness_pred.py` now exposes its original reconstruction as `reconstructed_blocks`.
  `scripts/check_correction_strength.py` sums those exact reconstructed counts and calls the
  installed bulk normalizer and PDS scorer. Zero bulk error and identical per-target PDS
  against a retained complete prediction; all four historical current-model seeds reproduced.
  This helper is **PDS-only**, not a full six-member score. Forty paired/null runs completed.
- Four context-model scientific tests passed; real IFNB held-out response poisoning leaves
  all 17 candidates bit-identical. A separate test proves the refactored full builder matches
  the saved builder cell for cell for both rounding modes and unchanged genes.
- A public, versioned DepMap 23Q4 prior is acquired in `processed/depmap_23q4/`: 18,442 genes
  x 1,100 cell models, publisher checksums verified. Canonical coverage is **300/300** official
  targets (at least 324 measurements each) and 1,656/1,661 b16 training targets (1,644 have
  at least 100 measurements). Preserve missing values. It measures cancer-cell knockout
  fitness, not CRISPRi RNA responses. No model has been fit with this new prior yet.

All paths in these bullets are relative to `vcc2026_data`. The next supervised experiment can
test this external gene prior with evaluated targets excluded from response-label fitting and
whole destination biologies held out. Compare against b16-style denoising and shuffled
feature/response correspondence. No user input is missing for this work. The sections below
retain historical plans and numbers; the completed reports and receipts take precedence.

## Reassessment update — later on 2026-09-14

Read [the new measured report](vcc2026_data/reports/strategy_reassessment_2026-09-14.md)
before interpreting the historical conclusions below. The recorded official best is 0.1217.
Authentication was recovered on 2026-09-14 from the existing Linux Secret Service keyring
(`vcc`, profile `default`); a live `vcc whoami` confirmed `can_submit=true`. No new token was
needed. The earlier restricted session's “Not logged in” result must not be taken as proof that
the saved credential is missing. Check keyring access before asking the user to regenerate it.
Never print the token or copy it into the repository. The b15 upload and the bounded
`scripts/finish_b15_sequence.py` runner were started. It waits for b15, then submits exactly
one prescribed follow-up (b17 if b15 beats b8's 0.12170879332843752 overall, otherwise b16),
using normal CLI allowance checks. It stops on failure or missing published scores.
**Inspect `vcc2026_data/reports/submission_receipts/confirmation_sequence.json` before another
upload**, to avoid duplicating the running sequence. Published scores are appended to
`vcc2026_data/reports/official_scores.csv`; receipts are in `submission_receipts/`.

- **A learned signed correction now beats plain mean transfer on KOLF.** Full `cell-eval2`,
  scale 2, seed 0: mean 0.15703 PDS; correction on mean 0.23586; correction on raw precision
  0.28562 (that input alone is 0.19565). Candidate-specific shuffled PDS is 0.02998 and 0.00289.
  Across four paired seeds, mean PDS is 0.16549 for plain mean, 0.23055 for its correction,
  and 0.25818 for precision plus correction: both win 4/4. Compare
  `reports/kolf_reassessment_2026-09-14/paired_seeds.csv` (seeds 1–3 use exact PDS-only checks).
- **Do not claim this proves a learned context function.** A training-label-shuffled residual
  operator scores 0.24491 on KOLF. The symmetric pooled regression learns source-disagreement
  covariance, consistent with denoising. A separate direct ranking-loss model gives only 0.16345,
  versus 0.16988 with shuffled training labels. Neither establishes useful supervised transfer
  beyond the denoising effect. These negative controls matter as much as the positive scores.
- **Jiang now has a nested whole-cell-line benchmark.** On official-gene intersections, raw
  pseudobulk PDS is mean 0.73674, centered mean 0.76485, learned 0.77825; 6/6 cell-line average
  wins. Model selection excludes the outer line; a held-out-response poisoning check passes
  for all 39 candidates. These CP10K pseudobulk ranks are not the CP50K reconstructed-cell
  competition score, and Jiang success alone does not prove challenge transfer.
- **CD4 also has biological replication:** four donors are already available in the small
  `pseudobulk/cd4_target_condition_donor` store. Matched donor effects and reproducibility
  diagnostics are now in `processed/cd4_donor_effects`. The older “only Jiang has replicates”
  statement below is incorrect. CD4 donor folds remain within one biology.
- **“Precision” is currently effective cell count**, not estimated gene-wise sampling variance.
  Shared controls and donor/guide heterogeneity are not represented by the current weights.
- **Ready at ×3:** `submissions/b16.vcc` (mean plus correction) and `submissions/b17.vcc`
  (precision plus correction), with validated receipts and SHA256s. Frozen operator:
  `processed/cross_source_residual/crsres128/model.npz`. Preserve the existing b15 clean weighting
  test: after login, b15 then b17 if weighting wins, or b16 if it loses. Use opaque names only.
- The new `scripts/pds_prebuilt.py` calls the installed bulk preparation and discrimination
  functions and matches a completed full run exactly. It permits fast paired-seed PDS checks;
  it does not produce or substitute for a full six-member score.
- Official rules distinguish the three live validation contexts from three additional final
  test contexts. Marker evidence supports A being T-cell-like, B epithelial, C basal epithelial;
  exact identities remain unverified.

All paths in the bullets, unless otherwise stated, are relative to `vcc2026_data/`.
The rest of this handoff is preserved as the historical state that motivated the reassessment.

---

## 1. The task, in one paragraph

Predict the transcriptomic response (18,533 genes) to CRISPRi knockdown of **300 target genes** in
**3 cell contexts (A, B, C) that appear in no training data**. The challenge gives only control
cells for those contexts plus the target list — there is no 2026 training set. All perturbation
data is public data we acquired ourselves. Submissions are scored by `cell-eval2` with the
`vcc2026` preset, **2 submissions per day**.

**Standing best: 0.1217, rank 214.** That is a "mean-transfer" baseline: for each target, average
its measured pseudobulk delta across our source datasets and submit the same delta to all three
contexts, scaled x3. Leader is ~0.28.

### Scorer facts that decide strategy

`Overall = plain mean of 6 members` (verified to 4 decimals). Our 0.1156 submission decomposed as
`(pds 0.5075 + mse 0 + nmae 0.0050 + fid 0.0165 + reach 0.1621 + jac 0.0023)/6`.
**pds is 73% of our score, reach 23% — the other four are ~4% combined.**

- `pds_cosine` deletes **all 300 panel target genes** from both operands first, so the on-target
  knockdown contributes exactly zero. It is a **discrimination rank** (is this target's predicted
  delta closer to its own truth than to the other 299?), not a raw cosine.
- `fid = k / max(n_pred, n_conf)` — over-calling costs as much as being wrong. Never prune/gate,
  only re-weight.
- `reach` = depth of a per-target prefix ranked by predicted |log2FC| holding 0.9 direction purity.
- `mse` has been integer 0 on every submission ever scored here.

---

## 2. Leaderboard history (the only ground truth)

| submission | what it was | Overall | rank |
|---|---|---|---|
| x1.0 | plain mean-transfer | 0.0771 | 342 |
| **x3.0 (b8)** | **plain mean-transfer** | **0.1217** | **214** |
| x4.0 (b10) | plain mean-transfer | 0.0868 | 328 |
| b11 | expression-gated, per-context, x2.25 | 0.1156 | 234 |
| b14 | precision-weighted sources, x3.5 | 0.0888 | 345 |

Per-metric for b14 vs the x3 best: pds **0.462 vs 0.504**, nmae **−0.158 vs 0.036**, fid 0.035 vs
0.016, reach 0.196 vs 0.180. nmae ≈ the known x4 value (−0.188) ⇒ **b14 was over-scaled**, and it
changed *two* things at once (weighting + scale), so it does not cleanly convict the weighting.

---

## 3. THE MOST IMPORTANT THING: the offline harness is unreliable

`harness/` holds a KOLF2.1J (iPSC) leave-one-cell-line-out harness scored with the real
`cell-eval2`. **It has mispredicted the sign of both submissions it was used to justify.**

### Why (measured, 2026-09-14)

Take the mean-transfer KOLF delta, **permute its 279 target rows** so none is in place, score it
normally (`bash scripts/xl_harness_x.sh permut 2.0 0`):

| metric | correct | target-SHUFFLED | reachable with NO target info |
|---|---|---|---|
| **pds_cosine** | 0.1570 | **0.0044** | **3%** |
| lfc_nmae | 0.1029 | 0.0632 | 61% |
| reach | 0.0703 | 0.0470 | 67% |
| sig_jaccard | 0.0878 | 0.0761 | 87% |
| fidelity | 0.7653 | 0.7438 | **97%** |
| **avg_score** | 0.1972 | **0.1558** | **79%** |

**Only pds discriminates targets.** 79% of the composite is reachable by a prediction that has the
targets completely wrong. Expression gating "won" the harness composite via `sig_jaccard` (87%
target-blind) and lost the leaderboard.

### Mandatory protocol from here

1. **Judge on `pds_cosine` only.** Never gate on `avg_score`.
2. **Always score a target-shuffled copy of the same delta** and report
   `(pds − shuffled)/(baseline − shuffled)`. Costs one extra 5-minute run.
3. Seed noise: 4 seeds of an **identical** delta span composite 0.1972–0.2085 (sd 0.0076). Compare
   only paired on `--seed`; treat <0.01 composite as undecided.
4. Harness pds (0.157) and leaderboard pds (0.507) are on very different scales. Do not port
   thresholds between them.

---

## 4. What has been tried

### 4a. Neural networks — 7 variants, none transfers

All predict a pseudobulk delta; all evaluated on held-out KOLF.

| # | model | objective | result on the held-out line |
|---|---|---|---|
| v1 | FC 56.5M, raw delta target | MSE | collapsed to ~0 output |
| v2 | FC, residual on mean-transfer, standardised | MSE | flat; residual is 8% of true size |
| v2cos | + cosine aux term | MSE+cos | identical to v2 |
| v3 | + top-2000 |y| weighting | weighted MSE | **worse from epoch 1** (weights leak the row's own noise) |
| v4multi | same, **58 contexts** instead of 8 | MSE | flat |
| gate_v1 | learned multiplicative gate, 4.8k params | MSE + L2 | multiplier came out flat (0.97–1.01) |
| setpred | anchor-free DeepSets over other contexts' observations, 11k params | MSE, then cosine | learns in-distribution, degrades on held-out |

**The setpred/cosine run is the instructive one.** It doubled cosine on 4 clean bio-group folds
(+0.0208, t=2.69, 4/4) and looked like the first real win. Then: subtract the component it predicts
identically for every target and its entire advantage goes **negative** (row cosine 0.0421 → −0.0020).
42.5% of its output energy was a shared all-target component vs 2.6% for mean-transfer. It had
learned the **generic perturbation response**, which inflates row cosine and destroys pds. Harness
composite 0.1333 vs mean-transfer 0.1972, pds 0.0152 vs 0.1570.

**Do not** re-justify a cosine objective with "pds is cosine" — that was my error. pds is a
discrimination rank; raw cosine's global optimum is the generic response.

### 4b. Preprocessing / denoising — one real effect

All precision-weighted, magnitude-matched, same code path, x2.0, judged on pds (shuffle floor 0.0044):

| change | pds | vs baseline |
|---|---|---|
| mean-transfer, equal-weight (current pipeline) | 0.1570 | — |
| **cross-source PRECISION weighting** | **0.2062** | **+0.0492** |
| + pseudocount c=0.2 instead of implicit c=1 | 0.2098 | +0.0036 more |
| per-source amplitude calibration | 0.1539 | −0.003, and it makes precision weighting worse |
| true logFC (remove the log1p compression) | 0.1323 | −0.025 |

**Precision weighting** is the whole story: `build_baseline_delta.source_profile()` already
inverse-variance weights rows *within* a source and returns the accumulated precision as `den`,
but `harness_pred.py` combined sources with a plain mean and discarded it (`prof, *_ =`).
Resulting weights: cd4 20%→49%, vcc2025 20%→33%, orion 20%→11.5%, replogle 20%→5.6%,
**feng 20%→0.8%** (feng has 23 cells/target vs cd4's 518).

**Pseudocount**: the pipeline's delta is `log1p(CP10K_pert) − log1p(CP10K_ctrl)`, i.e. c=1, which is
a compression not a fold change — 84–88% of genes sit below CP10K 1 where it captures ≤44% of the
true logFC. But removing the compression is *worse* (it amplifies noise on the low-abundance tail).
The whole sweep c=0.1…3.0 spans only **0.0071 pds**, shallow peak near c≈0.2. Real but tenth-order.

### 4c. Measurement corrections that invalidated earlier conclusions

- **Leave-one-CONTEXT-out inflates transfer ~6x.** cd4 Rest/Stim8hr/Stim48hr are one biology; holding
  one out leaves two near-replicates in training. Honest cross-*biology* base correlation is
  **0.0218**, not 0.1332, and expression gating's gain is **+0.0028**, not +0.0202. Always hold out
  whole bio groups: {HCT116} {HEK293T} {cd4 ×3} {feng} {vcc2025} {replogle}.
- **Measurement noise attenuates cross-line agreement 8–12x.** Two independent Replogle K562 screens
  (gwps vs essential, 2,053 shared targets) agree at cosine 0.28 ⇒ one measurement is ~30% signal.
  On Jiang (which has true replicates) noise-corrected cross-cell-line cosine is 0.066–0.173 across
  4 conditions vs 0.006–0.015 observed. **There is ~10x more real shared biology than raw numbers show.**

---

## 5. Current hypotheses (ranked, with what would falsify each)

**H1 — Precision weighting is real and b14 failed on scale, not weighting.** *Status: the only live
lead.* Falsify: submit precision weighting at x3.0 (file already built, §7). If it loses to 0.1217,
the weighting does not transfer and the harness is 0-for-3.

**H2 — The binding constraint is measurement noise, not absent biology.** Supported: 30% signal per
measurement; noise-corrected cross-line agreement ~0.1 vs ~0.01 observed. Implication: methods that
raise effective sample size (better weighting, replicate averaging) should beat methods that add
model capacity. Consistent with everything observed so far. Falsify: a denoising method that
demonstrably raises effective n and still fails to move pds.

**H3 — Direction is the wall, and it is not fixable by re-weighting.** On genes that really changed
in KOLF, mean-transfer's direction accuracy is **0.505** overall — 0.994 in the top 1% of |delta|,
0.592 at top 25%. A strictly positive multiplier cannot flip a sign, so every gate/weight scheme
inherits 0.505. Cross-source sign agreement is useless as a confidence signal (0.477–0.538 = chance).
**|mean-transfer| is already the best confidence ranking available.**

**H4 — "Which genes change" is largely answered by control expression, and it does not pay.**
Ranking genes by the target context's own control expression retrieves real DE genes far better than
our delta (AUC 0.826 vs 0.619, P@10 0.569 vs 0.132), and a classifier trained on real DE labels from
*another* cell line (Jurkat→KOLF) matches expression alone and adds nothing (AUC 0.828). But
expression-selected genes have **coin-flip direction** (0.525), and only `sig_jaccard` (1/6 of the
score, and 87% target-blind) rewards set membership without direction.

**H5 — Cell-level realism is untested.** Predictions are built by adding a *constant* delta to
sampled control cells, so within-group variance is just the control's. The DE metrics run Wilcoxon
on those cells. Whether realistic per-cell dispersion matters has never been measured. Cheap test
exists (perturb the reconstruction, same delta).

**Assumption, not measured:** contexts A/B/C are T-cell and keratinocyte (from an audit note, not
verified against the data). If true, CD4 is our closest source and Jiang's cancer lines are distant.

---

## 6. Data on disk

`/home/wonmin-kim/AIVC/vcc2026_data/`

| dir | what |
|---|---|
| `effects/<store>/` | `delta.npy` (n_pairs × 18533, NaN = unmeasured), `pairs.csv`, `measured_mask.npy` — the per-source pseudobulk deltas everything is built from |
| `pseudobulk/<store>/` | `sum_counts.npz` + `obs.csv` — raw summed counts, the correct primitive |
| `processed/xl/` | training matrices for the NN experiments (`contexts.csv`, `contexts.npz`, `pairs.parquet`) |
| `processed/xl_multi/` | same, cut into 58 contexts |
| `processed/jiang/` | **NEW, underused.** `{TGFB,IFNB,TNFA,INS}.{bin,obs.csv,var.csv}` — pseudobulk per (cell_line, perturbation, replicate) |
| `processed/baseline_deltas/` | official-300 deltas, incl. `b4_sevencontext_consensus_delta.npy` (the 0.1217 submission) and `b14_precision_delta.npy` |
| `harness/` | KOLF harness: `kolf_real.h5ad`, `de_real_kolf.parquet`, all `xl_delta_*_kolf.npy` candidates, all `scored_out_*.csv` |
| `raw/` | 409 GB single-cell source data, mostly unused (see below) |

**Source coverage of the official 300:** orion 299, cd4 294, kolf 282, replogle 272, feng 151,
vcc2025 25. Union = 300/300.

**Replogle measures only 7,681 of 18,533 features**; every other source measures ~18,000.

**`raw/jiang/` (19 GB) — the best unused asset.** 6 cell lines (A549, BXPC3, HAP1, HT29, K562, MCF7)
× 4 usable stimuli × ~50 shared perturbations × **2 biological replicates**, 1.38M cells. This is the
only properly factorial cross-cell-line design we have and the only source with replicates, which is
what makes noise-corrected measurement possible. Only 5 of its 106 genes are official targets, so it
is for **learning/measuring the context function**, not for the base table.
`Seurat_object_IFNG_Perturb_seq.rds` is **corrupt** (9 of 44 download blocks missing) — re-fetch from
Zenodo 14518762 if wanted.

**Single-cell available:** orion (118 GB parquet), kolf (177 GB h5ad), vcc2025 (33 GB), feng, scperturb.
**Pseudobulk only, cannot go to cells:** cd4, replogle.

**STATE assets:** `checkpoints/state/SE-600M/` (full `se600m_epoch16.ckpt` + protein embeddings),
`checkpoints/state/ST-SE-Replogle/`. `state` and `cell_load` are installed in the `ct_transfer` env.

### Environments
- `vcc2026` — anndata, cell-eval2, scoring, all delta building
- `ct_transfer` — torch (RTX 5090, 32 GB), all model training
- `rbase` — minimal R 4.6.1 + Matrix, for reading Seurat `.rds` (no Seurat pkg needed; use
  `attributes(o)`, not `slotNames`, which requires the class definition)

---

## 7. Key files

### Scoring (use these, do not re-invent)
- `scripts/xl_harness_x.sh <tag> <scale> [seed]` — score `harness/xl_delta_<tag>_kolf.npy` end to end.
- `scripts/harness_pred.py` — delta → predicted cells. `--delta-file` to score a pre-made delta,
  `--save-delta[-only]` to dump an aggregate, `--agg precision` for the new weighting,
  `--drop-store kolf_target_batch` is **required** for the KOLF harness or it reads the answer.
- `scripts/build_submission.py` — official submission. `--delta-per-context A=..,B=..,C=..` supports
  a different delta per context.

### Delta builders
- `scripts/build_baseline_delta.py` — `source_profile()` is the workhorse; it returns
  `(prof, den, rows, cells, match, prov, pairs)` and **`den` is the precision you want**.
- `scripts/build_precision_delta.py` — official-300 precision-weighted delta.
- `scripts/pseudocount_delta.py` — arbitrary pseudocount, asserts its transform pair is invertible.
- `scripts/xl_gate_delta.py`, `scripts/amp_calib.py`, `scripts/logfc_delta.py` — the tried-and-rejected variants.

### Analysis
- `scripts/jiang_transfer.py` — noise-corrected cross-cell-line transfer with a zero-signal null.
- `reports/xl_model_audit_2026-09-13.md` — 75-agent audit of the NN pipeline.
- `reports/xl_next_model_spec_2026-09-13.md` — 111-agent spec; §7 lists ~30 measured dead ends.

### Ready to submit
**`submissions/b15_precision_scale3p0.vcc`** (3.47 GB, packaged, validated, sha256 in
`reports/submission_receipts/b15.sha256`). Precision-weighted sources, scale 3.0 — **exactly one
variable changed** from the 0.1217 submission. This is the clean test of H1.
```
vcc submit submissions/b15_precision_scale3p0.vcc -m "b15"
```
Use an opaque model name and no description; the user does not want the method visible on the
leaderboard. Note that earlier submissions (b8, b10) already went up with descriptive names.

---

## 8. Traps that have already cost time

1. **Target ordering.** `harness_pred.py` indexes rows by `sorted(obs.target.categories)`;
   `harness/kolf_real.targets.csv` is in a *different* order — **1 of 279 positions matches**. Reading
   the CSV order silently permutes everything, no error, plausible-looking score. This invalidated
   three experiments in one session. Detector: two deltas correlating ~1.0 must score within noise.
2. **Non-invertible transform pairs.** `log((cp+c)/(cp0+c))` paired with `cp0*exp(lfc)+c` is not an
   inverse; it produced a non-monotonic scan. Assert the round trip.
3. **`pgrep -f` wait loops.** Failed 4 times in one session, both ways: a self-matching pattern waits
   forever, an over-specific one (`-x python3.12` when conda reports `python`) exits instantly and
   launches jobs concurrently. **Just run jobs sequentially in one script.**
4. **`cell_load` with `num_workers>0`** leaks memory and OOM-kills STATE training. Use 0.
5. **cell-eval2 `--set de.backend=pdex`** must be explicit; `auto` raises when a CUDA device is
   visible. pdex is ~18x faster than the scanpy backend.
6. **`ScaleLimitError`** if reconstructed per-cell totals exceed 1e6 — happens when an L1-renormalised
   delta is scaled up.
7. **`vcc cancel <entry> --yes`** frees a stuck upload and does **not** count against the daily limit.

---

## 9. Next experiments, in order

**1. Submit `b15` (H1).** Zero work, already built. Decides whether precision weighting transfers.
If it wins, rebuild the per-context gated and pseudocount variants on top of it. If it loses, the
KOLF harness is 0-for-3 and should be abandoned as a gate.

**2. Build a Jiang-based evaluation harness.** This is the highest-value infrastructure work. Jiang
gives 6 cell lines × the same ~50 perturbations × 2 replicates, so you can (a) hold out a whole cell
line, (b) noise-correct both sides, (c) build a proper zero-signal null. Everything the KOLF harness
cannot do. `scripts/jiang_transfer.py` already does the measurement half.

**3. The decisive model experiment, on Jiang.** Train on 5 cell lines, predict the 6th, compare
against the mean of the same 5, all noise-corrected. Our existing data cannot pose this question
(different targets per context, no replicates). If a learned context function beats the mean here it
transfers to the real task; if it does not, the direction is closed with real evidence rather than
the seven inconclusive attempts in §4a.

**4. Untried from the literature review** (agent output is in this session's transcript, not saved
to a file — regenerate if needed): `ash`/`apeglm` shrinkage of each source's LFC before combining
(TRADE reports replicate LFC correlation 0.16→0.90 after modelling sampling variance on this exact
Replogle data); limma-style **variance moderation of the precision weights** before inverting, since
noisy weights can make inverse-variance worse than an unweighted mean; and LINCS **MODZ** weighting
(weight by pairwise agreement between sources rather than by sampling precision).

**5. H5, cheap:** give the count reconstruction realistic per-cell dispersion instead of a constant
shift and re-score the same delta. Tells you whether cell-level realism matters at all before anyone
invests in a cell-level model.

### Do not repeat
Global scale tuning (leaderboard optimum is x3, x1/x4 both measured and much worse); pruning or
thresholding the delta (fid punishes it structurally); sign-agreement weighting; empirical-Bayes
shrinkage as a *ranking* (top-1% direction 0.800 vs 0.970 for plain |mt|); adding the other on-disk
contexts to the base table (0–5/300 panel coverage); per-source amplitude calibration.
