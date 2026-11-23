# Single-cell hybrid + emission optimisation — 2026-09-25/26

Official result: **S10c 0.21687, rank 49** (S10b 0.19215 rank 107, S8 0.15732 rank 225).

## 1. Single-cell information for the response (negative for the official task)
Harness = six held-out panels (H1 Flex, KOLF, RN K562/HepG2/Jurkat/RPE1), official preset, base blend30_e, paired seeds.
- Cell-cycle composition (`sch_phase.py`): share of response energy explained by composition H1 0.5 %, KOLF 0.1 %,
  K562 1.6 %, Jurkat 2.7 %, HepG2 9 %, RPE1 20 %. Phase shifts transfer across RN lines (Spearman 0.28-0.60).
- `sch_arms.py` (tally.csv/json): oracle composition (true pi_t) 6/6 wins (H1 +0.015, KOLF +0.041, RPE1 +0.33);
  transferred pi (other panels) mean +0.019 but H1 +0.0005 (only 25-33 H1/KOLF targets transferable), K562/Jurkat -0.02;
  phase-stratified shared draw 13/18, mean +0.004 (KOLF +0.016 3/3, H1 +0.0004). Half shift (h50) and oracle seeds scored.
- Knockdown escape (`sch_escape.py`, `sch_orion_rsel.py`): median escape H1 2 %, KOLF 31 %, RN 7-15 %, Orion HCT116 10 % /
  HEK293T 16 %; responder-weighted Orion sources change LFC by cos 0.96-0.97; rsel H1 0.2156 vs 0.2157 (null).
  Parity of the rebuild vs the tv2 sources: 1e-5 (denominator = counts on mapped genes, not the parquet library).

## 2. Emission (positive) — `sch_emit.py`, parity with fh_dual / tv2_emit exact
- Pinning sums to the cellshift expectation (em_cs_cs_k0): H1 0.202, KOLF -0.102 -> the MSE penalty is the cellshift
  amplitude itself, not block noise. Geometric blend of cellshift amplitude into dual sums: H1 0.226 (loses).
  Stratified draw on dual: nil. Multiplicative sums with cellshift shape and no re-allocation: H1 0.146 (DE collapses).
- NMAE anatomy on the scorer's own DE tables (`--cache-pred`, `sch_nmae_diag.py`, H1): dual under-states and randomises
  the sign of weak-effect genes' DE LFC (sign agreement by |pred| quintile 0.51/0.55/0.60/0.69/0.80 vs cellshift
  0.53/0.59/0.65/0.70/0.75); amplitude scan optimum c~0.75 (tiny), boosting hurts.
- **cpc**: after re-allocation, depth-tilt every expressed gene column (x_c ~ w_c L_c^s) so the SUM stays pinned and the mean
  per-cell composition equals the cellshift expectation over all controls. 6/6 panels at seed 0 (+0.020..+0.038), H1 +0.023/
  +0.026/+0.025 and KOLF +0.020/+0.022/+0.015 over 3 seeds. cpm (target c0 x FC) does not help.
- Official S10c vs S10b: NMAE raw 1.017 -> 0.916 (credit -0.025 -> +0.138), PDS/MSE identical, reach 0.292 -> 0.278.

## Code / artefacts
`scripts/sch_phase.py sch_arms.py sch_escape.py sch_orion_rsel.py sch_tally.py sch_emit.py sch_lfcdiag.py sch_nmae_diag.py`;
`tv2_emit.emit_dual_iter(cpc=True)`, builder spec `dual.cpc`, `tv2_final.py --dual 1.0,3000,0.5,ag,cpc`.
Not run: cc_oracle_emit/add and cc_xfer_emit/add decomposition arms (killed by oomd 13:34, low priority); rsel on the other
five panels (H1 null).

## 3. Iterative tuning on top of S10c (2026-09-26, rounds 1-7, H1/KOLF paired seeds; `sch_emit_tally.py`)
Reference = S10c settings (`em_mult_share_k3000_du0.5_ag_cpc`).
- Amplitude split: DE target gd 0.75/1.25 and sum gamma 0.8/1.2 all lose (current 1/1 is optimal).
- Recipe under dual+cpc: DepMap blend 0.4 +0.003, **0.5 +0.004**, 0.6 +0.002 (2/2); fusion share 0.75 PDS +0.014 but MSE/NMAE
  worse (H1 -0.010), per-target energy match (pt) does not rescue it; sim_soft / combo / no blend lose.
- Re-allocation: concentration share nt 0.2 +0.001 (2/2), 0.05 loses; strength floor le 0 / 0.6 nil; K 2000 loses;
  **K 5000 +0.004 on H1 (3/3) but -0.002 on KOLF (0/3)**; expression gate 2 CPM nil; cellshift pseudocount 0.5 / 2 nil.
- **bw50 + nt0.2**: H1 +0.0019/+0.0014/+0.0023, KOLF +0.009/+0.009/+0.010 (6/6), RN 1/4 (-0.004 mean).
- **bw50 + nt0.2 + K5000** vs S10c: H1 ~+0.006 (0.2756 -> 0.2813), KOLF ~+0.007. Parameter space is at a flat optimum.
Packages built (not uploaded): S10d (bw50 nt0.2), S10f (bw50 nt0.2 K5000), S10e (S4 half -> equal-weight fold changes,
not measurable on the harness because S4 is a learned model trained on the panels' own lines). S10e scheduled 09-27 00:03Z.

## 4. Scorer-structure levers that should transfer to the test set (2026-09-27)
Criterion: a change must follow from the scorer's definition (identical on the test set), be test-time computable from the
test context's controls + our own prediction (no labels), and win on all six harness panels.
- Official S10e (S4 half -> equal-weight fold changes) 0.21233 (-0.0046, PDS -0.015): keep S4. Official S10f (blend 0.5,
  nt 0.2, K 5000; H1 +0.006) 0.21540 (-0.0015): small knob gains do not transfer.
- **DE call counts** (fid = k / max(n_pred, n_conf), `sch_calls_diag.py`): oracle per-target n_pred = n_conf gives H1 fid
  0.539 -> 0.590, but n_conf's LEVEL is a property of the test data's power (median 145 H1, 27 KOLF, 16-444 RN) and a
  leave-one-panel-out predictor under-calls catastrophically (H1 fid 0.23-0.34). Over-calling (~3000) is the robust choice.
- **Reach ordering** (`sch_reach_sim.py`, `tmp ties`): the scorer ranks real-DE genes by our p-values; ideal |d| order gives
  H1 reach 0.288 vs realised 0.263 (S10c). The agreement ordering (ag) is worse than |d| on 5/6 panels offline. Not yet exploited.
- **MSE rebate kink** (`sch_mse_kink.py`, metrics/delta.py mse_unbiased_capped): num = ||mu_pred - mu_real||^2 - r min(C_pred, C_real)
  - C_real, r = min(1, B_pred / claim). Below the kink (B_pred < claim) noise energy is free and MSE falls linearly with
  <pred, truth>; above it MSE rises like a cliff. Our plug-in error is WORSE than predicting the control on H1 (D 25410 vs
  24503); the MSE credit is the rebate. The kink is set by the TEST DATA's noise (cells, depth): H1 ~1.0, official A/B/C
  ~1.0-1.1, K562 1.1, Jurkat 1.6, RPE1 1.8, HepG2 2.0, KOLF 2.8. Sums moved to the kink (DE side fixed): RN +0.013..+0.046,
  11/12 paired; H1 neutral. Estimator = controls' delete-1 jackknife x 1.17 (perturbed/control noise), margin 0.8 because
  overshoot is a cliff (KOLF gamma 2.8 MSE 0.998, 3.0 1.279). Implemented: `tv2_emit.kink_gamma`, builder `dual.kink`.
- **Offline exact sum-side evaluator** (`sch_sumeval.py`): MSE and PDS are functions of the pinned sums only; reproduces the
  scorer (H1 MSE 0.7613 = 0.7613, PDS 0.8581 vs 0.8576) in seconds, no emission noise. Used for the screens below.
- **Sum-side centring 0.5** (subtract half the across-target mean; override entries excluded): PDS up on 6/6 panels (H1
  +0.013, RN +0.02..+0.06; mean MSE -0.029) -- a hubness effect of the shared component on the cosine rank. **Source-family
  SNR weighting** (mu^2/(mu^2+sd^2/n)): PDS up 6/6. Matched filter on destination noise, gene reproducibility, per-gene
  transfer calibration, low-rank denoising: rejected. Real scoring RN seed 0: centring + kink = +0.016/+0.073/+0.068/+0.038.
- **ScaleLimit guard**: large sum amplitude + unscaled DE composition target makes cpc pile counts into the deepest cells
  (a HepG2 cell hit 1.85M > 1e6 -> scorer rejects the whole file). Both emitters now relax the tilt until cells < 8e5.
- **Concatenated DepMap + STRING (0.25) neighbours for UNMEASURED targets**: fallback H1 0.136 -> 0.156, KOLF 0.028 -> 0.041
  (list-construction control +0.004 only); offline neighbour quality +9-18 %. In `tv2_final_knn.py` for uncovered targets.
- Literature (PerturbMap, sampling-depth entanglement, class-overlap PDS, DeMixPert, VCBench, VCC2025 wrap-up, TransPert):
  no zero-shot method expected to beat the recipe; learned fusion loses to averaging there too.

## 5. 12-hour improvement loop (2026-09-27, serial queue `sch_queue_runner.sh`, no OOM)
Real scoring vs S10c settings (paired seeds):
- Round 19: ct0.5 + kink(mg0.8): H1 +0.0045 (3/3, PDS +0.013, DE unchanged), KOLF -0.003 (0/3; gamma 2.7 -> NMAE +0.007,
  reach -0.009; KOLF MSE ~1 is clamped so its MSE gain -0.073 earns nothing), RN +0.015..+0.069 (12/12); overall 15/18
  +0.032. snr+ct0.5+kink 15/18 +0.028. kink only 12/14 +0.022.
- RN (12 runs each): **ew50k0.25 + ct0.5 + kink = +0.060 (12/12)** and it removes the DE drift (NMAE +0.0075 -> -0.0004,
  reach +0.005). DE side following the sum amplitude (gdf 0.5 / 1): no help / NMAE +0.045. DE-side centring: -0.021 (0/4,
  reach -0.041). DE-side SNR weighting: -0.017 (0/4). Reach-ordering strength by expression (lcp): +0.002 / -0.001.
- Offline (exact evaluator, deployable kink, LOPO): expression-weight shape is a plateau (c0 20-200, kappa 0.15-0.35,
  mean J 0.11-0.12); sum-only neighbour blend 0.5 helps H1/KOLF; target-adaptive blend +0.013 J (5/6, too small);
  promoter prior x2 +0.001 PDS (keep); matched filter, coverage energy, per-gene calibration rejected; 30 % K562 coverage
  loss keeps the transform gains (H1 PDS +0.011).
- Scorer DE LFC = log2 of ARITHMETIC means of CPM-normalised cells (vcc2026 `mean_calc: arithmetic`); the drift at large
  sum amplitude is the composition pin becoming unattainable, not a geometric-mean effect.
- Bug fixed: option-key clash (`pc` pseudocount vs PC removal -> `pr<k>x<lam>`).

## 6. Loop continued (2026-09-27 03:00-, serial queue)
- **MSE claim is capped by the TEST data, not by our scatter** (scorer: rebate = min(B_pred, sum_p min(C_pred,p, 1.0 x C_real,p))).
  Cached pred moments of the S10c emission: C_pred/C_real = 1.47 (H1, 94 % of rows above) and 1.18 (KOLF, 85 %). The #247
  cap binds, so emitting more per-cell scatter cannot raise the claim; cells per target are fixed at exactly 400
  (`vcc prep` EXPECTED_CELLS_PER_PERT). The kink (B_pred = sum C_real) is the whole MSE budget. Closed.
- **NMAE magnitude profile is not the lever** (`sch_nmae_shape.py`, cached S10c tables): per-|pred|-quintile optimal scales
  gain 0.008 (H1) / 0.010 (KOLF) in-sample and are WORSE cross-panel (H1->KOLF +0.001, KOLF->H1 +0.008). NMAE is set by the
  sign accuracy of real-DE genes (H1 0.53..0.78 by quintile, KOLF 0.50..0.60), i.e. by delta quality. Closed.
- **Top-PC removal pr3x0.25** (on ew200 + ct0.5 + kink, seed 0): PDS up on 6/6 (+0.008..+0.012 vs ew50 arm); composite H1
  0.2830 (ct-only 0.2789, ew50 0.2751), KOLF 0.0489 (0.0422 / 0.0455); RN pays MSE (true signal on the shared axes of
  essential screens), net -0.017..-0.001. Official library top-1/3/5/10 axes carry 10 / 18 / 22 / 28 % of library energy
  (`tv2_pc_prep.py` -> `reports/submission_2026-09-23/lib_pcs_official.npz`); builder `sum.pcr`, `tv2_final --sum-pcr k,l`.
  Clean no-ew test on H1/KOLF x3 seeds queued (031).
- Knockdown efficiency per source (target's own fp/fc, >=10 CPM): K562 median 0.85-0.87 (2 % < 0.3), CD4 0.79-0.85,
  HCT116 0.70, KOLF 0.61, HEK293T 0.50 (48 % < 0.5). Per-target KD weighting would mostly touch HEK293T/KOLF on the ~28
  non-K562 targets: small expected effect, not run.
- Leaderboard 09-27 (`reports/leaderboard_2026-09-27.json`, latest entry per team): top-30 PDS 0.83-0.89, MSE 0.57-0.79,
  NMAE 0.81-0.95, reach 0.22-0.43. A cluster of teams at PDS ~0.86 / MSE ~0.75 / NMAE ~0.90 / fid 0.528 / reach 0.29 shares
  identical DE metrics (a common public template?). [another team] (#42, PDS 0.842) uses only K562 GWPS + Orion + CD4 with a
  cross-source sign-agreement filter and a jackknife-credit kink; [another team] (#103) uses a truth-free knee at g 2.4-2.7.
- **Per-source per-coordinate z-shrinkage before fusion** (`tv2_model.Cfg.zshrink`: v * z^2/(z^2+k), z = v / robust
  across-target SD of that gene in that source x sqrt(n_med/n); harness `zs<k>` both parts, `zf<k>` fusion only;
  `sch_zs_screen.py`, per-target energy held at the reference, S10g sum side). Fusion only (the only form the official
  recipe can take, since part 0 is the learned S4): PDS down on 5/6 panels, J -0.062..+0.004 -> rejected. Part 0 only
  (zp4, the equal-weight S4 stand-in): PDS up 6/6 (+0.003..+0.016) but MSE worse 5/6 (H1 +0.015, HepG2 +0.044) at the
  gamma>=1 floor -> budget test running (`sch_zp_floor.py`).
- Official S4 half vs fusion half (context A/B/C): per-target cosine median 0.57 (q10 0.46, q90 0.69), equal energy;
  S4 puts 25-27 % of its energy in genes < 5 CPM (fusion 17-20 %). Two partly independent views -> averaging helps.
- zp4's H1 MSE loss is NOT the gamma floor (B/claim at gamma 1 is 0.838 < 1): it is the amplitude profile itself. Kink rule
  accuracy (`sch_zp_floor.py`): at the rule's gamma the realised B/claim is 0.58 (K562) .. 0.87, i.e. estimated/actual claim
  0.72..1.07; margin 0.8 keeps every panel below the cliff (margin 0.9 would too, max 0.96; not worth the official risk).
- **Budget matching** (`sch_budget_match.py`: a sum-side transform is a direction change, rescaled so our own B_pred equals the
  base's at gamma_b = max(1, kink)). J = 2.23 dPDS - 1.05 dMSE (~ 6 x d composite), ct0.5 base:
  ew50  H1 +0.001  KOLF +0.004  K562 +0.013  HepG2 +0.048  Jurkat +0.060  RPE1 +0.011  (neutral at the floor = adaptive rule)
  pr3   H1 +0.017  KOLF +0.018  K562 -0.014  HepG2 -0.035  Jurkat -0.006  RPE1 -0.007  (RN pays MSE on shared axes)
  ew+pr H1 +0.019  KOLF +0.022  K562 -0.002  HepG2 +0.022  Jurkat +0.054  RPE1 +0.005  (5/6)
  Real scorer (027 arm ew200+ct0.5+pr3 vs ct0.5, seed 0): H1 +0.004, KOLF +0.007, K562 -0.0005, HepG2 +0.005, Jurkat +0.010,
  RPE1 +0.004 (5/6). Caveat: official S10c sums carry far less shared/PC structure (shared 0.009, top-3 PC 3.75 %) than H1
  (0.022 / 10 %), so the official gain of ct/pr is expected to be small (~+0.001-0.003).
- Future direction (not started): z-shrink the S4 refit's TRAINING rows (per biology, same MAD noise model) -- part-0
  shrinkage raises PDS 6/6 on the harness and S4 is the official part 0; needs the ext table pipeline + a leaderboard slot.
- **Guarded PC removal** (sum side, k 3, lambda 0.75 x clip((0.25 - f)/0.10, 0, 1), f = our own sum-side energy on the top-3
  library axes: RN 0.32-0.40, H1 0.098, KOLF 0.072, official A/B/C 0.037). Offline budget-matched sweep (`sch_pr_sweep.py`):
  H1 J +0.017 / +0.031 / +0.036 at lambda 0.25 / 0.5 / 0.75 (lib10x0.75 +0.040), KOLF +0.018 / +0.040 / +0.062; RN
  strongly negative beyond 0.25 (true essential-gene signal lies on those axes) -> guard sets lambda 0 there. Panel-own PCs
  (pan3x0.5) behave like library PCs. Builder `sum.pcr.guard`, harness `pr<k>x<l>g`, `tv2_final --sum-pcr 3,0.75g`.
- **Budget matching wired** (`tv2_emit.budget_B/budget_match`, builder `dual.kink.budget_match`, harness `bm`): transformed
  sums are rescaled to the UNtransformed recipe's B_pred at max(1, kink). Official preview: base kink A 0.893 / B 0.937 /
  C 0.960 (gamma 1); ew is dropped in all three (kink with ew 0.82-0.87); S10i s = 1.02; S10j (sum-side S4 share 0.75,
  |Dj|/|D| 1.045) s = 0.989. Library-only z-shrinkage (neighbour estimates) rejected: H1 J -0.003..-0.007, KOLF -0.014..-0.027.
- Literature for the PC-removal mechanism: PDS is a cosine nearest-neighbour rank, and a few real deltas aligned with the
  dominant generic programmes act as hubs for many queries. Removing dominant principal directions before cosine ranking is
  the "all-but-the-top" post-processing of word embeddings (Mu & Viswanath, arXiv 1702.01417, cited from memory);
  hubness reduction measurably improves cosine semantic retrieval ("Hubness Reduction Improves Sentence-BERT Semantic
  Spaces", https://www.alphaxiv.org/abs/2311.18364). Our guard is the case those papers warn about: when the task signal itself lies
  on the dominant axes (essential-gene screens), removal destroys it.
- **Source families** (`sch_wm_screen.py`, drop / re-weight on both parts + library, per-target energy held, ct0.5 + kink):
  dropping ALL stem-cell sources (KOLF, FENG, H1) hurts the somatic RN panels (K562 J -0.006, HepG2 -0.012, Jurkat 0.000):
  FENG and H1 carry real signal for somatic lines (Jurkat noFENG -0.011, noH1 -0.013). KOLF alone is the exception
  (K562 +0.003, HepG2 +0.001, Jurkat +0.024). **CD4 x2 hurts everywhere, including the lineage-matched RN_Jurkat
  (-0.038; K562 -0.043, HepG2 -0.017)** -> no lineage-matched source boost (the dropped S7b idea stays dropped).
  RPE1: noKOLF +0.007, noFENG +0.009, noH1 -0.004, noStem +0.015, CD4x2 -0.028. Across the four somatic RN panels: KOLF
  removal 4/4 positive, H1 (the Flex-chemistry source) removal 0/4. Most of the harness KOLF effect sits in part 0 (equal
  weights, KOLF ~10 %); officially part 0 is S4, whose biology-balanced refit includes KOLF rows -> candidate for a later
  slot: S4 rebuilt without the KOLF biology (`tv2_s4_rebuild.py` already has a leave-one-biology-out path). Official
  mirror (fusion 0.25 + neighbour library only) is being screened (`sch_kolf_mirror.py`).
  H1 (stem-cell panel): noKOLF -0.043, noFENG -0.023, noStem -0.074, CD4x2 -0.013 -> stem-cell sources are strongly
  lineage-matched signal for H1 and mildly harmful (KOLF only) for somatic lines. A lineage-aware KOLF weight (drop for
  somatic final contexts, keep for stem-like ones; tv2_identity gives the lineage) is the deployable form if the official
  mirror (fusion + library, S4 untouched) is positive.
  **Official mirror** (`sch_kolf_mirror.py`: KOLF removed from the fusion + both library parts, part 0 untouched):
  K562 +0.0020, HepG2 +0.0008, Jurkat +0.0050, RPE1 +0.0018 (somatic 4/4 but ~+0.0004 composite), H1 -0.0099.
  Dropping KOLF+FENG: +0.0001 / +0.0004 / +0.0053 / +0.0056; all three stem sources: -0.0012 / -0.0007 / -0.0012 / +0.0032.
  Far below the adoption bar -> keep the source weights. (The part-0 effect would need an S4 rebuild; low priority.)
- **Real scorer, H1 seed 0** (028a/028b): S10c 0.2741; ct0.5+kink 0.2789 (+0.0048); +pr3x0.25 0.2807 (+0.0066);
  **+guarded pr3x0.75 0.2819 (+0.0078; PDS 0.8576 -> 0.8881, MSE 0.7613 -> 0.7787 at gamma 1, NMAE/fid/jac flat, reach
  -0.001)**; guard stat 0.0977 -> lambda 0.75 as designed. Budget-matched arm (restores B_pred, i.e. the MSE loss) pending.
  **Budget-matched (the exact S10i configuration) 0.2842 = +0.0100 over S10c** (PDS +0.0304, MSE +0.0045 with gamma_sum
  1.058 instead of +0.0174 at gamma 1, reach -0.0019, NMAE/fid/jac flat). Budget matching is worth +0.0023 here.
- **Real scorer, KOLF seed 0: exact S10i configuration (ew kept, kink 2.46) 0.0580 vs S10c 0.0423 = +0.0157** (PDS +0.021,
  MSE -0.091 -> ratio below 1, NMAE -0.002, reach +0.012); ct+kink alone -0.0002, ew50+ct+kink +0.0032, ew200+ct+pr0.25
  +0.0066. S10i configuration: H1 +0.0100, KOLF +0.0157.
  Same without ew on KOLF: 0.0533 (+0.0110; gamma 2.82): ew adds +0.0047 where the kink leaves budget (the adaptive rule
  keeps it there and drops it on A/B/C, where gamma is at the floor).
- **DE-side guarded PC removal** (dpr3x0.75g on the S10i configuration): H1 0.2820 vs 0.2842 (-0.0022; NMAE +0.0014,
  fid -0.0026, jac -0.0007; sums identical) -> rejected; the generic axes on the DE side carry real sign information.
  KOLF: 0.0528 vs 0.0533 (-0.0005; reach -0.0036). DE-side removal rejected on both official-like panels.
- **Per-target guard** (`sch_pc_target.py`, budget matched; f_t = target t's own energy fraction on the top-3 axes):
  official A/B/C f_t 10/50/90/99 % = 0.004/0.024/0.072/0.147 (0.3 % > 0.25), so per-target and panel guards coincide
  officially. H1: panel guard J +0.036, per-target(0.15-0.25) +0.027; KOLF +0.062 vs +0.070; RN (panel guard 0):
  per-target(0.15-0.25) +0.004 / +0.032 / 0.000 / +0.005. A max(panel, per-target) rule would keep H1/KOLF and add the RN
  gains -> only relevant if a final panel were essential-gene-like; noted, not implemented.
- **Centring level with guarded PC removal** (`sch_ct_with_pr.py`, budget matched, vs S10i's ct0.5): H1 ct0.75 J +0.0075,
  ct1.0 +0.0132; KOLF +0.0063 / +0.0085; RN ct0.5 is optimal (ct0.75 -0.012..-0.036, ct1.0 -0.06..-0.10). Same guard
  logic applies (ct 0.5 + 0.5 x guard) -> final-phase refinement candidate; official gain small (shared share 0.009).
  Real-scorer check queued (048: ct1.0 on H1, and with ew on KOLF).
- **Rehearsal 4** (`scripts/tv2_final_sim4.sh`: tv2_final --s4 --blend 0.3 --dual 1.0,3000,0.5,ag,cpc --kink 1.17,0.8
  --sum-center 0.5 --sum-ew 50,0.25 --sum-pcr 3,0.75g on the simulated D/E/F panel): spec translated correctly (budget
  matching on, guard [0.15, 0.25]); per context guard f 0.065/0.066/0.065 -> lambda 0.75, ew dropped (0.837/0.820/0.855),
  kink_base 0.875/0.884/0.907 -> gamma_sum 1.040/1.038/1.038; build 1452 s; `vcc prep` rejected only the fake labels D/E/F
  (known; use the real labels on 10/22).
- **S10i configuration vs S10c, H1 paired seeds 0/1/2: +0.0100 / +0.0097 / +0.0100 (mean +0.0099, 3/3)**.
- **S10i configuration vs S10c, KOLF paired seeds 0/1/2: +0.0157 / +0.0115 / +0.0128 (mean +0.0133, 3/3)**. With H1
  (+0.0099, 3/3): 6/6 paired runs positive on the two official-like panels.
  Real scorer (048): ct1.0 instead of ct0.5 in the S10i configuration: H1 +0.0024 (PDS +0.0072, MSE +0.0028), KOLF +0.0017
  (MSE -0.0096, reach -0.003). Consistent but ~+0.002 (half of that officially, shared share 0.009) -> not rebuilt;
  final-phase option: ct = 0.5 + 0.5 x guard(f).
- **Exact S10i configuration vs S10c on all six panels (real scorer, seed 0; 049 for RN):** H1 +0.0100 (3 seeds +0.0099),
  KOLF +0.0157 (3 seeds +0.0133), RN_K562 +0.0180, RN_HepG2 +0.0838, RN_Jurkat +0.0793, RN_RPE1 +0.0436 -> **6/6**. The guard
  sets lambda 0 on every RN panel (f 0.32-0.40) and 0.75 on H1/KOLF; ew kept where the kink allows (KOLF, RN), dropped at the
  floor (H1, official). On K562 budget matching holds the sums at the margin target (gamma 0.881) where the older floor
  rule overshot it (0.4531 there, because K562's claim estimate is 28 % low) -- the conservative choice is deliberate.
- Masking override entries / panel-target columns in the PC projection: |J| <= 0.0004 (H1, KOLF) -> irrelevant.
- Sum-side neighbour blend with PC removal (`sch_blend_with_pr.py`, budget matched): H1 bw0.5 J +0.0081 with pcr vs +0.0048
  without (PDS cost -0.0013 vs -0.0020); KOLF +0.019 vs +0.043 (pcr already takes part of the gain). Official expectation
  ~+0.001 -> not pursued.
- **Fidelity call-set anatomy** (`sch_fid_calls.py`, replay reproduces the scorer: H1 0.5405, KOLF 0.5122): H1 calls/target
  median 2820 (81 % inside our top-3000 |d|); sign accuracy inside 0.563 vs outside 0.551. Suppressing the outside calls
  LOWERS fid (0.5405 -> 0.5291; top-1500 only 0.4844: recall-limited targets), KOLF neutral -> over-calling stays. Closed.
- DE-side generic-axis AMPLIFICATION (dpr3x-0.5g): H1 -0.0004 (NMAE +0.0028, fid/reach +0.0005). With removal at -0.0022,
  the DE side is at its optimum on the generic axes. Closed.
- DE-side amplification KOLF +0.0005 (reach +0.009, NMAE +0.0019) -> net ~0 with H1 -0.0004. Closed.
- **S10j proxy** (`sch_s10j_proxy.py`; part 0 = the equal-weight stand-in, weaker than S4 officially): sum-side part-0 share
  0.75 vs 0.50 under the S10i sum side: H1 J -0.005 (PDS -0.010, MSE -0.016), KOLF -0.012; share 0.25: H1 -0.013, KOLF
  +0.010. The 50/50 mix is near-optimal on the harness; S10j's official sign is uncertain (S4 carried +0.015 PDS vs the
  stand-in in S10e while costing MSE 0.007) -> S10j is an information experiment, S10i the expected gain.
- **Gene-dependent call set from a power model** (`sch_callcount_power.py`; user question 2026-09-27 "learn the between-cell
  distribution per gene"): gene g of target t predicted significant if |kappa lfc_pred| > z SE_g, SE_g from the TEST
  context's own control CPM (CV, 400-cell / control counts); (kappa, z) cross-fitted H1 <-> KOLF. Replay (fid, Jaccard):
  H1 current (0.5405, 0.0695), oracle top-n_conf (0.5852, 0.0849), power top-n_hat (0.5417, 0.0715), predicted set
  (0.5434, 0.0720); KOLF current (0.5122, 0.0158), oracle (0.5408, 0.0243), power (0.5142, 0.0197), predicted set
  (0.5005, 0.0194). The power model recovers ~5-10 % of the oracle headroom (~+0.0015 H1 / +0.003 KOLF composite at best,
  before any emission cost): n_hat vs n_conf Spearman 0.27 (H1) / 0.13 (KOLF) -- the number of genes that truly move per
  target is set by each knockdown's true strength, which our deltas rank poorly. NMAE is shape-independent (arithmetic means,
  already pinned by cpc). Per-gene distribution learning is not a lever until per-target effect strength is predictable.
- **Target's own expression does NOT predict its knockdown strength** (`sch_target_expr.py`): H1 Spearman(n_conf, own
  control CPM) +0.056, (effect norm, own CPM) +0.054 -- the panel targets are expressed (287/300 >= 20 CPM) and the few
  low-CPM targets are among the STRONGEST; KOLF +0.157 / +0.070. Our predicted energy is the better per-target strength
  signal (H1 +0.46 for n_conf and effect norm; KOLF +0.12 / +0.13). Target-expression gating is not a lever.

## 7. S4 improvement round (2026-09-27, user: "일단 S4모델 개선부터 시도해보자"; `reports/s4_improve_2026-09-27/`)
- **The harness can now carry the REAL S4 as part 0** (`sch_emit.delta_s4`, token `s4<key>`): the five lines' saved
  leave-one-biology-out predictions (`reports/ext_contexts_2026-09-21/<line>/film_v2_a0.25.npy`, target lists identical
  to the KOLF / RN panels) and a leave-H1-out refit for H1 (`s4i_refit_h1.py`), mixed exactly as the builder
  (0.5 x 3 x S4 + 0.5 x 1.6 x fusion, blend 0.3, global energy). vs the equal-weight stand-in: PDS up on 6/6 (H1 +0.018,
  RN +0.006..+0.014, KOLF +0.006), MSE worse on KOLF/RN but BETTER on H1 (-0.014); H1 J +0.055 -- the official S10e sign.
- **A. Row-level z-shrinkage of the S4 training table** (`s4i_table.py` k 4 / 1, per source group MAD noise; `s4i_refit.py`,
  reproduction of the saved refit_v2 to 1.5e-6): vs the original S4, energy-matched / raw: H1 -0.002..-0.025,
  RN -0.004..-0.045, KOLF +0.002..+0.010 -> REJECTED (S4's refit already denoises; shrinkage removes weak real signal).
- **B. Knockdown-quality weights across biologies** (`s4i_refit_q.py`, pair weight clip(q / median_q(bio), 0.1, 1)^kappa,
  kappa 1/2; 3-5 % of pairs < 0.5): |dJ| <= 0.0005 on all six panels -> no effect.
- **Refit rank 16/64/128**: H1 within +-0.001, mixed elsewhere. **FiLM residual alpha 0.5/1.0**: RN better
  (Jurkat +0.021/+0.044) but KOLF worse (-0.010/-0.032), H1 unmeasurable -> keep 0.25.
- **S4 share of the sum-side mix** (`sch_s4_share.py`, real S4): 0.75 (S10j) H1 -0.0003, KOLF -0.011, RN -0.016..-0.077;
  0.25 H1 -0.030 -> 0.5 is near-optimal; **S10j is expected neutral-to-negative -> keep it unsubmitted**.
- **S10i vs S10c with the real S4 part 0** (`sch_s10i_with_s4.py`, sum side only): H1 +0.003 composite (ct +0.0025, pcr
  +0.0004; S4's refit already centres its means), KOLF +0.020, RN_K562 -0.004, other RN +0.02..+0.05 (kink). Official
  S10i expectation revised to ~+0.002-0.003 (was +0.004).
- Verdict: the S4 refit is at a local optimum for input denoising, quality weighting, rank and share. What remains for
  S4 is new signal (data), not reprocessing.
- **Flex-platform per-gene calibration (ML idea, `ml_flex_calib_diag.py`)**: per-gene slope of observed vs fusion-predicted
  (other families) fold changes per destination. Spearman g(CD4 Flex) vs g(H1 Flex) +0.371, but g(K562 3') vs g(H1) +0.575,
  HCT116 +0.306, HEK293T +0.280, K562 vs HCT116 +0.653 -> the shared per-gene structure is generic (not Flex-specific), and
  generic per-gene calibration was already rejected (LOPO). No platform model to learn from the two Flex datasets.

## 8. Multi-angle bias review + Orion somatic panel + discriminative residual (2026-09-27 afternoon)

User: "이 외에는 개선방향이 없는지 다시 다각도로 확인해줘 증명안된 편견은 조금 줄이고", then (mid-turn) "why is our PDS low but fid
high?" and the S4 + contrastive-residual proposal. All numbers measured today; J = 2.23 dPDS - 1.05 dMSE (~6 x composite).

**Assumptions re-tested**
- Weak-gene signs are a coin flip (`ml_sign_prior.py`): H1 lowest-confidence quintile our sign 0.530 vs a gene-level library
  direction prior 0.539 (KOLF 0.52 vs 0.61); replacing low-confidence signs moves NMAE/fid/reach <= +0.001 composite. Holds.
- Test-time context signal from the destination's own controls (fluctuation-response, `ml_fdt_context.py`): -corr(target, g)
  in control cells discriminates at chance (0.49-0.52) on H1/KOLF/RN K562/Jurkat; blending never helps. Closed.
- Destination-line DepMap essentiality (`ml_depmap_line.py`; Jurkat = A and CAL-33 = C are in DepMap CRISPR, HeLa is not):
  Spearman with true per-target energy 0.16-0.36 < our own delta norm 0.46-0.61; LPO ridge gains nothing.
- MSE already optimal at the kink: only the TOTAL budget is. Per-target sum allocation at a fixed budget (PDS-invariant):
  strength oracle H1 -0.084 raw MSE (~+0.015 composite), RN -0.07..-0.12, KOLF +0.018. The 09-26 `so` oracle arm was invalid
  (raw LS slope, median 0.06, never renormalised). Reachable part small: TRUE strength of one RN line predicts another only at
  Spearman 0.29-0.60 (300 shared targets, `ml_strength_ceiling.py`) - our norm is already at that ceiling; power rule
  a = norm^k helps the three non-essential panels (H1 -0.016/-0.024 MSE at k 0.5/1) and hurts all four RN panels.
- NMAE amplitude (`ml_nmae_amp.py`): global multiplier already optimal (x1.00); per-target oracle only H1 0.874 -> 0.853.
- MSE formula (`ml_mse_formula.py`): MSE ~ 1 - 2 kappa rho sqrt(B/claim) (H1 0.748 predicted vs 0.766) -> official rho ~0.13,
  top-3 teams (MSE 0.57-0.59 at PDS similar to #4-#10 with 0.70-0.76) ~0.21-0.24.
- Leaderboard gap to the top 10 (09-27, best entry per team): PDS 34 %, MSE 24 %, NMAE 23 %, reach 11 %, fid 4 %, jac 4 %.
- fid (`direction.py _components`): scored set = our significant calls with a DEFINED real LFC (not real-DE genes); we call a
  median 2,820/target vs 145 real DE, so fid = sign accuracy of the call set (0.523; peers 0.513-0.517; top 0.547). Low PDS +
  high fid = emission strong, biology weak. Sum side is not over-generic (shared share 0.026 pred vs 0.057 truth on H1).

**Orion HCT116 somatic panel (new harness panel OR_HCT116 / OR_HCT116X)** - `orion_panel_build.py`, 263 official targets x
median 163 cells, 12k NTC, 12.1k UMI/cell on the official axis (70 % of Orion UMI; batches 10k-31k), no thinning.
- Metric profile (S10c): PDS 0.690 / MSE 0.531 / NMAE 0.745 / reach 0.707; official 0.810 / 0.778 / 0.916 / 0.278 is far
  closer to H1 (0.858 / 0.761 / 0.874 / 0.263) -> H1 stays the primary harness (data-generation regime > cell type).
- Family transfer (disc rank, `ml_family_transfer.py`): K562 alone 0.267 beats the recipe 0.274 on its targets; CD4 0.360,
  HEK293T (same lab) 0.376, FENG 0.413, KOLF 0.461.
- Real scorer: S10i configuration (kink 1.85 -> gamma 1.93) +0.0786 over S10c (MSE 0.531 -> 0.100; PDS +0.036; NMAE +0.039).
  At the official budget (mg0.2 -> gamma 1): ct0.5 alone +0.0038 (MSE -0.024); + guarded pcr 3x0.75 -0.0015 (PDS +0.044,
  MSE +0.080, reach -0.015) -> S10i sum side +0.0023 vs S10c on a somatic panel (H1 +0.010): no downside found. The guard read
  f = 0.064 (full lambda) although the truth carries 2x H1's energy on those axes (raw 0.053 vs 0.027, `ml_pc_truth.py`).
- Offline (budget-matched; PDS exact, MSE unreliable on Orion: negative at the kink, C_pred < C_real): dropping stem-cell
  sources PDS +0.007 (OR) / +0.029 (OR-X, J +0.033) - H1/KOLF, the panels that kept them, are stem cells themselves; K562
  up-weighting raises PDS but costs MSE on every panel (H1 x2 J +0.003, KOLF +0.004, OR -0.06).

**S4 + discriminative residual (user/ChatGPT proposal)** - `ml_dr_prep.py` (vcc2026) -> `ml_dr_train.py` (ct_transfer GPU) ->
`ml_dr_eval.py`; outputs `reports/sc_hybrid_2026-09-25/disc_residual/`. Pseudo-panels = each allowed source family as truth,
input = fusion-weighted mean of the other allowed families (panel family never used); r = U c, U = top-64 PCs of training
truths, c = 2-layer MLP([U^T X, ESM-2 15B PCA-64, log|X|, #families]); loss InfoNCE (tau 0.1, pull to own truth / push from
other targets' truths) + own-cosine + |r|^2/|X|^2; best step on 10 % held-out rows. Applied on the SUM side only (DE side
untouched), D' = D + alpha r |D|/|X| on the S4-faithful recipe (real LOBO S4 on H1/KOLF/RN) and on the S10i sum side,
budget-matched, exact offline PDS/MSE.
- Within the library it learns: val disc 0.238 -> 0.193 (H1 run); held-out WHOLE family improves at alpha 0.1-0.5 on all 9
  (K562 0.175 -> 0.154, RPE1 0.142 -> 0.110); shuffled-label control stays at chance.
- On the panels, WITH ESM-2 (as proposed): H1 PDS -0.003 / -0.029 / -0.062 at alpha 0.05 / 0.25 / 0.5 (J -0.012 at 0.1 on S10i
  sums) while KOLF/RN are positive -> the embedding is a lookup key for how a target looks in the (same-lab) source
  families; H1 (other lab, Flex) is hurt. Random embedding also hurts; own-cosine only hurts; shuffled hurts everywhere.
- WITHOUT embedding (target-agnostic reshaping of the transferred mean in the 64-PC basis): J > 0 on 7/7 panels x 3 seeds at
  alpha 0.05-0.25, also target-disjoint (panel targets removed from training): at alpha 0.1 H1 +0.008..+0.013, KOLF
  +0.004..+0.010, OR +0.025..+0.041, RN +0.007..+0.026 (42/42 > 0). Gains are MSE (global cosine); H1 PDS +-0.001.
  ~+0.0015-0.003 composite on H1 -> the "PDS up" goal is not met; a small robust MSE denoiser remains (candidate for the
  final-phase build at alpha 0.1-0.25, or one controlled leaderboard slot).

## 9. Data-quality audit: does every dataset help? (2026-09-27, user: "우리 데이터들은 전부 도움이 돼? ...")

Code: `ml_src_quality.py` (per-source stats), `ml_family_transfer.py` (each family alone), `sch_src_ablate.py` (leave-one-
dataset-out; mode b = real LOBO S4 kept, fusion + neighbour library ablated; mode a = everywhere incl. the equal-weight
part 0), `sch_kd_filter.py` (+ `tv2_model.Cfg.min_kd`, default 0 = unchanged), `s4i_refit_excl.py` / `s4i_refit_or.py`
(S4 refit without a biology, via `ext_baselines.CONFLICT_GENE`, rows AND empirical means), `sch_s4_nokolf_eval.py`.
Logs `reports/sc_hybrid_2026-09-25/src_ablation/`. All evaluations: S10i sum side, budget-matched, exact offline PDS/MSE
(Orion: PDS exact, MSE not), J = 2.23 dPDS - 1.05 dMSE.
- Quality: median on-target KD K562 0.80, CD4 0.73-0.78, HCT116 0.66, KOLF 0.56, HEK293T 0.47 (lowest); cells/target
  FENG_GWN 21 (lowest), FENG_GWF 41; H1 1,071 (25 official targets). Alone, KOLF transfers worst on every somatic panel
  (disc 0.37-0.41 RN, 0.46 Orion; K562 best 0.11-0.27).
- Mode b (what a submission can change without an S4 rebuild): no family is harmful across panels; every family is
  needed somewhere (K562 -0.014 H1 / -0.140 RN Jurkat; HCT116 -0.006..-0.025 on 5/6; H1 -0.002..-0.021 on all).
  The one large positive (KOLF panel, drop K562 +0.050) is a stem-lineage effect (H1 then dominates the fusion).
- Mode a on the somatic official-target panel: dropping KOLF PDS +0.018 (Orion) / +0.025 (Orion-X); stem panels need it
  (H1 -0.052); RN K562 +0.003, HepG2 -0.011, Jurkat +0.028, RPE1 +0.009 -> 5/6 somatic panels positive. CD4 out: RN
  Jurkat/RPE1 +0.013/+0.011 but Orion (official targets) -0.218 -> keep CD4.
- **S4 without the KOLF biology** (refit rows + empirical means; S4 refit cos 0.86-0.88 vs the original): J +0.002 /
  +0.015 / +0.014 / +0.006 on RN K562 / HepG2 / Jurkat / RPE1, Orion +0.015, Orion-X +0.028 (6/6 somatic > 0); plus KOLF out
  of fusion + library: -0.001 / +0.011 / +0.015 / +0.007, Orion +0.026 (PDS +0.007), Orion-X +0.055 (PDS +0.010).
  Mechanism: weakest-KD, lineage-mismatched source holding ~20 % of S4's equal-biology weight (5 biologies per official
  target) vs ~2 % in the fusion. Candidate for one leaderboard A/B (needs a tv2_s4_rebuild exclusion option).
- Row-level knockdown filter (drop source rows with on-target KD < 0.3/0.5/0.7): H1 -0.017/-0.046/-0.115, Orion +0.005/
  -0.094/-0.235, RN mostly negative, only KOLF positive -> weak-KD rows still carry response signal. Closed.
- Sub-sources: CD4 states each useful (small); FENG_GWN (21 cells) mixed (KOLF +0.020, Orion +0.017, H1 -0.013) -> keep.
- RN HepG2/Jurkat/RPE1 essential screens (KD 0.85+) are harness-only, not in the official source set (0/300 A/B/C targets);
  include them in the final-phase sources if the new panel contains essential genes.

## 10. S10k = S10i + no-KOLF sum side + no-embedding residual (2026-09-27 16:00-, user: "s10i제출 일단 정지하자. 우리가 했던걸로 개선한 모델 만들수있잖아")

S10i's scheduled upload was STOPPED at 16:13 PDT (unit aivc-tv2-upload-sched-20260928 stopped, no receipt, queue emptied;
original queue kept as `upload_queue_20260928.json.held_S10i_20260927`).
- **Design** (all changes on the SUM side; DE side bit-identical to S10i): sums from S4 without the KOLF biology (refit rows +
  empirical means; saved FiLM residual), fusion without KOLF, DepMap-neighbour library without KOLF; + no-embedding learned
  residual (target-disjoint `noesm_strict`, KOLF excluded from its data, mean of 3 seeds) at alpha 0.25; then the S10i sum
  transforms (centring 0.5, guarded top-3 PC removal) and budget matching to S10i's untransformed budget.
- **Harness** (`sch_s10k_eval.py`, builder budget rule, exact offline sum-side metrics; J vs the S10i structure; worst seed):
  RN K562 +0.027 (+0.015), HepG2 +0.032 (+0.031), Jurkat +0.040 (+0.031), RPE1 +0.023 (+0.018), Orion +0.073 (+0.063),
  Orion-X +0.097 (+0.091); no-KOLF alone -0.001 / +0.011 / +0.015 / +0.007 / +0.026 / +0.055. Stem panels (H1/KOLF) prefer
  KOLF and were not used for this decision (A/B/C are somatic); the residual alone was positive on H1 (+0.008..+0.024).
- **Official inputs** (`reports/submission_2026-09-27/`): `S4noK_ctx{A,B,C}.delta.npy` (`tv2_s4_nokolf.py`; saved S4 refit
  reproduced exactly, max abs 0.0; S4noK cos vs S4 0.857, norm x1.09), `knn_official_noK.npz` (`tv2_blend_prep_excl.py`;
  reproduces knn_official to 6e-8 without exclusion), `resid_S10k.npz` (`ml_dr_official_prep.py` + `ml_dr_train.py`;
  context genes >= 10 CPM 8,843 / 8,399 / 8,854), spec `spec_S10k.json`.
- **Builder** (`tv2_submit_build.py`): new sum-side options `spec_override` (replace parts / blend for the sums only) and
  `resid` (file, alpha). Context A: pcr guard 0.037 -> lambda 0.75, ew dropped (0.754 < 1), kink_base 0.893, gamma_sum 0.931
  (S10i 1.019); DE delta identical to S10i; sums cos 0.934 vs S10i.

## 11. Why is our PDS low? (2026-09-27 evening, `ml_pds_anatomy.py`, `ml_weak_targets.py`)

- Split-half ceiling (pred = bulk of half the target's real cells, truth = the other half) vs our S10i-structure sums:
  H1 0.963 vs 0.895 (vs full truth 0.907). By true response strength (unbiased energy quartiles): weakest 0.866 vs 0.686
  (66 % of the gap), 25-50 % 0.990 vs 0.935 (20 %), 50-75 % 0.997 vs 0.969 (10 %), strongest 0.999 vs 0.991 (3 %).
  Orion's split-half ceiling (0.671) is BELOW ours (0.727): 81 cells/half at 12k UMI is too noisy to rank; strongest
  quartile 0.931 vs ours 0.799.
- So the PDS loss is weak responders: their small, target-specific true responses are ranked by the truth itself but not
  by our transferred means (source noise at 100-200 cells/target and line-specific responses). Official PDS 0.810 << H1 0.907
  suggests the official non-essential panel has more weak responders.
- Test-time rules keyed on our predicted norm (predicted-weakest 25 / 50 %): full centring + PC removal mixed (RN negative);
  shrinking their trans part x0.5 (cis/knockdown overrides dominate) J +0.013..+0.026 on H1/KOLF/RN K562/Jurkat but
  PDS only +0.0006..+0.004 -> an MSE reallocation, not weak-target discrimination.
- **Official (2026-09-28 01:18Z): S10k 0.22300 = new best (+0.0061 over S10c)**: PDS raw 0.810 -> 0.8205, MSE 0.778 -> 0.7647,
  DE metrics unchanged (NMAE 0.9157, fid 0.5235, jac 0.029), reach 0.278 -> 0.2748. Predicted +0.006 -> transferred as estimated.

## 12. Claude <-> Codex debate loop, 2026-09-27 evening (`reports/sc_hybrid_2026-09-25/debate_loop/`)

Codex runs as `~/.local/bin/codex exec -s read-only -C <repo> --skip-git-repo-check -o OUT resume <session> -` (model
gpt-6-astra, reasoning medium; the plugin's Node launcher is absent, so the CLI is called directly). Shared state
`STATE.md`, round log `ROUNDS.md`, per-round prompts/replies `codex_round*.md`. Session 01a0e5b1-0831-7031-a675-4985734cf20f.

**Round 1 - protocol critique (accepted).** Aggregate by LABORATORY (RN / Orion / H1 / KOLF), not "6/6 panels"; evaluate the
exact submitted construction (parity: un-matched no-KOLF S4, seed-MEAN residual); internal residual-validation leaks (PCA
basis fitted before the split; `lofo_diag` leaves the held-out family inside other families' inputs) - external panels
unaffected, so S10k's verdict stands. Parity harness: `sch_geo_eval.py`.
**Round 1 experiment - scored-geometry residual: REJECTED.** Training the residual under the scorer's own weighting
(diagonal w = a/(1+a), a = 5e4 x control fraction; `ml_dr_famctrl.py` + `ml_dr_train.py` variant `*_geo`) gives J
current -> geo: H1 +0.018 -> +0.013, RN K562 +0.031 -> +0.027, HepG2 +0.017 -> +0.008, Jurkat +0.041 -> +0.037,
RPE1 +0.025 -> +0.026, Orion +0.075 -> +0.057 (RN lab mean 0.029 -> 0.025). Diagnostic: cos(exact scored delta, LFC)
0.86-0.96 and cos(exact, diagonal) 0.94-0.996 -> the geometry mismatch has little headroom (Codex's own falsification branch).
**Round 2-3 - distal cis prior: real but bounded small.** Mean signed z by TSS distance (`ml_cis_distance.py`; H1 weak
targets / our prediction): 0-5k -7.8 / -3.5, 5-20k -1.08 / -0.59, 20-100k -0.36 / -0.12, 100-500k -0.06 / -0.04 -> real
CRISPRi repression beyond the 5 kb prior that we under-predict, strongest on WEAK targets. But the oracle upper bound
(`ml_cis_oracle.py`: replace our LFC with the real one at those pairs) is H1 5-100k PDS +0.0022 / J +0.0053 (5-500k +0.0044 /
+0.0111), Orion +0.0019 / +0.0083 (5-500k +0.0097 / +0.0421), RN Jurkat +0.0005 / +0.0021; the 0-5k control leaves 1-6 pairs
and ~0 effect (existing override verified). Coverage is 2-17 pairs per target. Codex's caveats recorded: the script replaces
the c=50 LFC (not the scored delta) on the ORIGINAL S4, so it is an optimistic diagnostic rather than a strict bound.
**Round 4 - the transfer ceiling (the decisive result).**
- `ml_weak_reach.py` (RN K562/Jurkat/HepG2, 300 shared targets, genes >= 10 CPM): per target, WITHIN-line split-half cosine
  vs CROSS-line cosine, by true-strength quartile: weakest 25 % 0.108 vs 0.013-0.025 (transferable fraction 0.12-0.23),
  25-50 % 0.15-0.18 vs 0.06-0.07 (0.33-0.41), 50-75 % 0.24-0.26 vs 0.12-0.14 (0.44-0.58), strongest 0.36-0.42 vs
  0.21-0.27 (0.49-0.74). **Weak responders' real responses are measurable within a line but almost not shared between lines.**
- `ml_transfer_ceiling.py` (predict RN_Jurkat with other lines' REAL responses for the same targets, through our own sum-side
  pipeline; same lab = optimistic): ours PDS 0.9078 / MSE 0.3187 (weak quartile 0.818); oracle transfer K562 0.8863 / 0.7392
  (0.726), HepG2 0.8707 / 0.6905, RPE1 0.8194 / 0.7106, mean of the three 0.8974 / 0.4075 (0.788); energy-matched to ours
  0.8924 / 0.8741 / 0.8273 / 0.9004. **Our recipe already beats perfect same-lab transfer from one or three other lines,
  overall and on the weak quartile.** A 50/50 blend of ours with the energy-matched 3-line oracle is 0.9134 / 0.2992
  (J +0.033) - i.e. the only remaining transfer gain is ensemble diversity against lines we already use as sources.
**Conclusion: better source transfer / weak-responder modelling is closed by the data, not by method.** The weak-quartile PDS
gap (H1 0.686 vs split-half 0.866) is line-specific biology. Remaining levers are scorer-structure / emission and the DE side.

## 13. Final-phase rehearsal, step 1: COVERAGE is the dominant final-phase risk (2026-09-27 night)

Codex round 7 named "false confidence on poorly covered new targets" as the failure mode A/B/C never tested. Measured:
- **Our validation panels are 97.7-100 % well covered** (>= 3 independent source families after the own-line exclusion;
  `ml_coverage_audit.py`): H1 300/300, RN_K562 300/300, Orion 257/263, KOLF 274/279. The limited (1-2 family) targets that
  do exist score much worse (Orion 0.581 vs 0.745, KOLF 0.558 vs 0.612), on n = 5-6.
- **A mixed panel with 20 % of targets deleted from EVERY source and from the neighbour library** (`ml_uncovered_sim.py`,
  production DepMap fallback): panel PDS -0.0399 (H1) / -0.0314 (Orion); the hidden targets themselves fall 0.910 -> 0.710
  and 0.712 -> 0.563; retained targets are unaffected (+0.000 / -0.002). **~ -0.007 composite per 10 % uncovered.**
- **CORRECTED 2026-09-27 night.** The earlier "the library covers 2,581 of 18,533 genes (13.9 %)" was WRONG: 2,581 is the
  size of the source FILES we built for our own target lists, not the data's coverage. Counting the underlying pseudobulk
  stores (>= 20 cells per target) on the official axis: Orion HCT116 17,163 / HEK293T 17,355, CD4 11,921, KOLF 10,956,
  Replogle K562 9,383, FENG 4,603, H1 300, RN HepG2/Jurkat/RPE1 ~1,900 each; **union 17,978 of 18,533 = 97 %**.
  Over the whole axis (what a randomly drawn new panel looks like): **74.2 % well covered (>= 3 families), 22.7 % limited
  (1-2), 3.1 % unmeasured**, versus 97.7-100 % well covered on our validation panels. Adding the RN families moves this
  only 74.2 % -> 74.5 %, because their targets are essential genes Orion already covers - so the RN wiring is a QUALITY
  change (highest on-target KD), not a coverage change.
- **Fallback sweep** (8 neighbour sources x K x weight power, 6 panels, hidden-target PDS; `ml_uncovered_sim.py`).
  **CORRECTION (2026-09-27, after reading `tv2_final_knn.py`): the FINAL-PHASE path already uses concatenated
  DepMap+STRING neighbours (SA = 0.25) for UNCOVERED targets and already builds the RN sources for the neighbour
  library.** The `depmap` baseline of this sweep is the A/B/C path (`tv2_blend_prep.py` -> `H.NB`), where coverage is
  100 % and the fallback never fires. So the sweep does NOT show an unexploited fallback gain; what it adds is the
  Orion caveat below. Aggregated by laboratory (H1 / Orion / KOLF / RN mean) with
  Codex's rule "no laboratory may lose more than 0.005":
     concat0 K10 pw2  +0.022 / +0.027 / +0.009 / +0.009  lab mean +0.0168  PASS   <- selected
     concat0 K10 pw1  +0.019 / +0.034 / +0.006 / -0.001  lab mean +0.0146  PASS
     concat0.25 K10 pw2 +0.059 / -0.022 / +0.024 / +0.032  lab mean +0.0233  FAIL (Orion -0.022)
     concat0.5 K10 pw2  +0.048 / -0.043 / +0.049 / +0.021  lab mean +0.0188  FAIL (Orion -0.043)
  Relative to what the final phase ALREADY uses (concat SA=0.25), concat0 is better only on Orion (+0.049) and worse on
  H1 (-0.037), KOLF (-0.015) and RN (-0.023), so the shipped choice stands; the open question is only why Orion dislikes
  the STRING component. Selecting per panel would have picked concat0.25 on H1 and concat0 on Orion - the selection bias
  Codex warned about. Either way this is insurance, not a validation gain (A/B/C are fully covered).
- **The RN screens (HEPG2 / RNJURKAT / RPE1) are built for the final-phase NEIGHBOUR library (`tv2_final_knn.py` jobs)
  but carry NO weight in the FUSION** (`tv2_model.FAMILY` lacked them; `Cfg.weights` has no such keys; the official
  source dir has no RN npz). They hold 896 targets at the highest on-target KD we have (0.85+). WIRED 2026-09-27:
  `tv2_model.FAMILY` now registers them (inert without a source file AND a positive weight), `tv2_sources.py` accepts
  `RN_HepG2 / RN_Jurkat / RN_RPE1` by name, and `tv2_final.py --rn-sources [--rn-weight W]` builds them for the new panel
  and adds W to both parts' weights. Regression: the S10k context-A delta and override rebuild identically (max abs
  1.2e-7 / 5.8e-8 = float32 round-trip), so A/B/C is untouched. Harness value where already used: dropping the three
  costs J -0.009 (H1) / -0.019 (Orion).
- ci5 (the round-5 composition fix) confirmed on a second panel: H1 +0.0022, KOLF +0.0004, NMAE better on both
  (0.8743 -> 0.8697, 0.9630 -> 0.9593); PDS/MSE exactly unchanged. Second H1 seed still running.

## 14. Alternative formulations, 2026-09-27 night (user: "challenge our assumptions, prioritise empirically testable ideas")

Four non-obvious reformulations, all tested offline (budget-matched, exact sum-side PDS/MSE; J = 2.23 dPDS - 1.05 dMSE).

**Diagnostic first (`ml_gram_geometry.py`) - PDS as a MATCHING problem over 300 vectors, not 300 independent predictions.**
Gram of our predictions vs Gram of the truths (scorer space, panel genes removed):
- H1: off-diagonal cosine pred +0.001 (sd 0.049) vs truth +0.035 (0.080); effective rank 245 vs 193. Orion: 241 vs 224.
  **On non-essential panels our predictions are LESS correlated / MORE spread than the truth**, yet removing the top library
  axes (which decorrelates further) still gains +0.004 PDS -> "match the truth's Gram" is the WRONG objective; PDS is a
  nearest-neighbour rank and what matters is hub removal, not correlation level. Hypothesis killed before implementation.
- RN_Jurkat is the mirror image: effective rank 86 vs truth 171, top eigenvalue 0.236 of trace vs 0.122. **This is the
  mechanism behind the S10i guard**: essential-gene panels are over-concentrated, non-essential ones are not.
- Our pairwise over-similarity grows with the predicted similarity itself (spearman +0.46 H1 / +0.46 Orion / +0.96 RN),
  and on Orion the top-1 % most similar predicted pairs have truth cosine 0.019 vs our 0.115 (spurious), while on H1 they
  are genuine (0.314 vs 0.293).

**1. Spectral flattening** (`ml_spectral_shape.py`), s_i -> s_i (s_i/s_1)^-rho, the natural generalisation of the guarded
PC removal: never beats it. H1 J -0.010..-0.107 at rho 0.1..1; RN_Jurkat gains PDS (+0.0095) but loses MSE (+0.09).
A TARGETED operator (three known generic axes) beats an untargeted one, because flattening also destroys real structure.

**2. Destination control-coherence filtering** (`ml_ctrl_coherence.py`), new operator, test-time computable: project the
prediction on the top-k PCs of the destination's own control cells (library size regressed out). Only 2-8 % of our energy
lies there, and "keep top-k only" is catastrophic (PDS -0.10..-0.21) -> our signal is NOT in the controls' coherent
subspace. Shrinking those directions raises PDS and costs more MSE; AMPLIFYING them lowers PDS and gains MSE
(H1 k=200: PDS -0.008, MSE -0.027, J +0.010; Orion J +0.040 but its MSE is unreliable).

**3. Rank-space fusion** (`ml_rank_fusion.py`): combine the two recipe parts as normal scores instead of a linear average.
Negative where it matters (H1 J -0.045, Orion -0.213; RN_Jurkat +0.009). "Rank average then restore |linear|" is 0.0000,
i.e. the two agree on SIGNS almost everywhere - so what rank-space discards is the MAGNITUDE, which is our confidence
ranking. Winsorising the parts (a milder version) also loses. Magnitudes carry real information.

**META-CONCLUSION: the sum side sits on a PDS <-> MSE efficient frontier and the shipped operator is on it.** Every
transformation tried today (spectral, destination-coherence, rank-space, plus the earlier scored-geometry residual and
MODZ) moves along that frontier rather than beyond it, and the shipped ct0.5 + guarded top-3 PC removal + budget matching
is the best point found. Combined with section 12 (source transfer saturated) the validation-side search is finished;
remaining effort belongs to final-phase coverage and robustness (section 13).

## 15. Is the DATA optimal? exclude / process / acquire (2026-09-27 night, user question)

**Coverage, corrected and then put in context.** The stores cover 17,978 of the 18,533 official genes (97 %). Mean number
of INDEPENDENT source families (>= 20 cells) per gene, which is what predicts per-target PDS:
  whole gene axis 3.87 (74.2 % at >= 3 families) | genes >= 1 CPM in context A 4.70 (94.6 %) | >= 10 CPM 4.92 (97.1 %) |
  >= 50 CPM 5.13 (98.1 %) | **the actual 300 validation targets 5.27 (99.3 %)**.
Panels of targets are drawn from EXPRESSED genes (an unexpressed gene cannot show a knockdown effect), and published
screens target expressed genes too, so a realistic new panel sits near 4.9-5.1 families, not 3.87. Measured slope of
per-target PDS per extra family (`coverage_audit.json`, correlational so an upper bound): H1 +0.033, Orion +0.050,
KOLF +0.043, RN_K562 +0.007. **Realistic final-phase coverage penalty ~ -0.35 families x 0.04 = -0.014 PDS ~ -0.005
composite**, not the -0.015..-0.03 implied by the earlier (wrong) 13.9 % figure.

**Exclude**: done in section 9 - no family is harmful across panels; the one real case was KOLF inside S4 for somatic
destinations, already shipped in S10k. Weak-knockdown row filtering hurts.

**Process** (`ml_proc_knobs.py`, the least-explored axis; sum side S10i, budget-matched):
- fusion pseudocount c (shipped 50; the scorer's bulk_lognorm at TS 5e4 has an effective pseudocount near 20, so this
  had a mechanistic reason to re-test): c=20 gives H1 J +0.006 but Orion -0.018; c=10 H1 +0.006 / Orion -0.041;
  c=100 H1 -0.007 / Orion +0.000. **Inconsistent -> keep c=50.**
- `part0 c=50` looked like a winner (H1 J +0.011, Orion +0.016) but is a HARNESS-ONLY artifact: in the official build
  part 0 is the S4 delta file, not a fusion with a pseudocount. Not actionable.
- centring fraction (0.75 / 1.25) and the per-(source, target) cell floor (10 / 50) are inconsistent or negative.
=> No actionable processing knob; the pipeline's preprocessing is already at a local optimum.

**Acquire**: coverage is not the binding constraint (97 % of expressed genes), so new data buys CONTEXTS, which section 12
showed is exactly what we lack. Unused public sets in `reports/data_catalog_2026-09-24/CATALOG.md`: THP-1 (Yao),
IGVF hPSC / cardiomyocyte / NPC, Tian-Kampmann iPSC neurons, Frangieh melanoma, Papalexi THP-1, Schmidt primary T.
None is squamous or lymphoid-matched to A/B/C, and all are small. The one high-value item is the **Ker-CT / 3D epidermal
keratinocyte CRISPRi set (10x Flex, the challenge chemistry, closest biology to CAL-33)** - access is "on request only",
so it is a request to make, not an experiment to run. K562-only additions (Norman 2019, Adamson, Dixit) would add no
coverage: Replogle already gives 9,383 K562 targets.

## 16. Top-team gap re-decomposition + HYPERPARAMETER ENSEMBLING (2026-09-28, user: "현재 모델과 상위권 모델들을 비교했을때, 가장 개선 가능성이 높은 실험 루트")

**16.1 The gap, from all 1,188 leaderboard entries** (`reports/leaderboard_2026-09-28.json`; S10k 0.22300, team rank 48).
Gap to the top-10 median, 0.052 total, by scaled member / 6:
MSE **+0.0218** (42 %) | NMAE **+0.0161** (31 %) | PDS **+0.0148** (28 %) | reach +0.0057 | jac +0.0035 | fid +0.0017.
Raw -> composite slopes (local fit at our own value): PDS +2.234, MSE -1.051, NMAE -1.654, fid +3.343, reach +1.127,
jac +2.724 -- i.e. +0.005 composite needs dPDS +0.0134 / dMSE -0.0285 / dNMAE -0.0181 / dfid +0.0090 / dreach +0.0266 /
djac +0.0110.

**The structural correction to section 14.** Binned over all entries, the best PDS reachable RISES as MSE falls
(MSE 0.95-1.05 -> best PDS 0.835; 0.55-0.65 -> 0.875), and **30 entries beat us on PDS, MSE and NMAE simultaneously**.
So the PDS<->MSE frontier of section 14 is a property of RE-WEIGHTING ONE FIXED PREDICTION, not a bound on the task:
a change that lowers estimator ERROR moves off it. 158 entries beat our PDS and 156 our NMAE, so neither sits at a
public ceiling. Useful reference points: `[another team]` has our MSE (0.749 vs 0.765) and near our PDS (0.844) but NMAE 0.835
vs our 0.9157; `[another team]`/`[another team]`/`[another team]` are a three-way identical cluster (PDS 0.860 / MSE 0.760), i.e. a shared
public recipe that already beats our PDS by 0.04.

**16.2 DE call-set pruning: CLOSED** (`ml_callset.py`, `ml_order.py`, `ml_order2.py`, on the cached H1 scorer tables;
`fid_raw = k / max(n_pred, n_conf)` and `sig_jaccard` reproduced exactly).
- We call a mean of **2,884 genes per target**; the real side's n_conf is mean 899, median 145, p10 2, p90 2,851.
- ORACLE budget M = n_conf: fid 0.5405 -> 0.5662 (p_adj order) / **0.5801 (|lfc_pred| order)**, jac 0.0695 -> 0.0835,
  composite +0.028. Oracle ORDERING reaches fid 0.9285, so the ordering headroom is large.
- **Every test-time budget rule loses**: quantile-mapping our rank on n_p1e6 / n_p001 / sum|lfc| / n_lfc05 / n_pred onto
  the n_conf distribution gives -0.089..-0.121 composite at multipliers 1.0-2.0; a FIXED M (50..1500) gives -0.026..-0.179.
- Why: n_conf is unpredictable from our side (Spearman 0.256-0.326) and the penalty is ASYMMETRIC -- under-calling keeps
  n_conf in the denominator while k collapses (M = 0.5 n_conf -> fid 0.307); over-calling only dilutes slowly.
=> keep calling everything. The original "never prune or gate" rule stands, now with the quantitative reason.

**16.3 NMAE is sign-limited, confirmed.** `mean|lfc_pred - lfc_real| / mean|lfc_real|` over real-significant genes. For a
prediction a x (sign-correct fraction s) the optimum is a = 1 whenever s > 0.5 and the value is 2(1 - s), so magnitude
calibration cannot pay unless the signs are right: the ORACLE per-quintile rescale (scales 3.64 / 2.20 / 1.55 / 1.12 /
1.08, i.e. we under-predict by 3-20x) moves H1 NMAE only 0.8471 -> 0.8395. NMAE follows delta quality, like MSE and PDS.

**16.4 HYPERPARAMETER ENSEMBLING -- the one lever that moves both metrics** (`ml_hp_ensemble{,2,3}.py`; shipped S10i sum
side, budget-matched, exact offline PDS/MSE; J = 2.23 dPDS - 1.05 dMSE ~ 6 x composite).
Section 15 judged each processing knob alone and kept c=50 because c=20 helped H1 and hurt Orion while c=100 did the
reverse. That is the signature of estimators with OPPOSITE biases, and their MEAN need not lie between them.

| arm (members) | H1 | OR_HCT116 | KOLF | mean J | ~composite |
|---|---|---|---|---|---|
| ens mc4 = base + min_cells 5/10/50/100 | +0.0146 | +0.0566 | +0.0440 | +0.0384 | +0.0064 |
| **ens mc7 = base + min_cells 5/10/15/30/50/100/200** | **+0.0159** | **+0.0762** | **+0.0416** | **+0.0446** | **+0.0074** |
| ens mc4+w (+6 family-weight multipliers) | +0.0170 | +0.0689 | +0.0193 | +0.0351 | +0.0059 |
| ens nb (8 neighbour-similarity sources) | -0.0088 | -0.0050 | -0.0133 | -0.0090 | negative |

- **PDS and MSE improve TOGETHER on 3/3 panels** -- the first arm all year to do so; everything in section 14 traded them.
- The mechanism is measurable: every one of the 6 family-weight members is individually J < 0 on H1 (-0.005..-0.045) and
  their mean is J +0.0061; min_cells 200 alone is J -0.184 on H1 yet mc7 > mc4 there. The gain tracks member DIVERSITY --
  `center_frac` members have cosine 1.0000 with the base and gain exactly 0.0000, `min_cells` members 0.54-0.98 and gain
  the most.
- **Neighbour-source diversity does NOT work**: all seven alternatives (STRING, DepMap+STRING, concat, ESM-2, GenePT,
  GO-BP, STRING-link) are individually negative on H1/KOLF and biased the SAME way, so their mean loses and adding them
  to `ens mc` dilutes it. DepMap co-essentiality is genuinely the best single neighbour source.
- **Diagnosis of WHY min_cells is the live axis** (`tv2_model.predict_lfc`): within a family every (source, target) row is
  combined with EQUAL weight (`fnum += v; fden += 1`) behind a hard `min_cells` cut -- a 20-cell measurement counts as
  much as a 500-cell one. Averaging over thresholds is a step-function approximation to weighting rows by an increasing
  function of their cell count. `tv2_model.Cfg.ncw` (new, default 0.0 = inert, so every existing build is bit-identical)
  does it directly: row weight n / (n + ncw). Screen `ml_ncw.py`.
- Wiring: `sch_emit.py` arm token `ens<pool>` (pools `mc`, `mc7`, `mcw`, `mcwc`) applies the ensemble to BOTH the sum and
  the DE delta, so the real scorer can show whether NMAE/reach move too -- the offline evaluator only sees the sums.

## 17. 12-hour exploration loop, 2026-09-28 (user: explore untried directions, use papers/posts and Codex)

**17.0 The offline sum-side evaluator had a systematic MSE bias, and it is now FIXED and VALIDATED.**
`ml_hp_ensemble*.py` computes the budget target `Bt` ONCE from the base recipe and matches every candidate to it;
`sch_emit.score` instead sets `Lb0_` from the CANDIDATE and matches to its own budget. Averaging shrinks the
across-target spread, so the ensemble's own budget is smaller (H1 `kink_base` 0.901 base -> 0.939 ensemble) and matching
to it lands at a lower amplitude and costs MSE. Diagnosed by Codex (round 9) from the saved `info` blocks; confirmed by
the real scorer, where the offline-vs-real MSE discrepancy was a CONSTANT +0.0100..+0.0105 across three arms.
New arm token `bmb` budget-matches to the NON-ensembled recipe (fully test-time computable: build both deltas, use the
base's budget). With it the offline prediction reproduces exactly: ens MSE -0.0081 -> -0.0081, ens+dis100x0.5
+0.0059 -> +0.0059. **Offline PDS was accurate all along** (predicted +0.0042 / +0.0062, measured +0.0042 / +0.0062).

**17.1 Real scorer, H1, official `delta_s4` structure (base composite 0.2976)**
| arm | composite | delta | PDS | MSE | NMAE | fid | reach | jac |
|---|---|---|---|---|---|---|---|---|
| ensP1 (min_cells ensemble on part 1) | 0.2995 | +0.0019 | +0.0015 | +0.0024 | -0.0025 | +0.0009 | +0.0024 | +0.0002 |
| **ensP1 + bmb** | **0.3011** | **+0.0035** | +0.0012 | **-0.0081** | -0.0023 | +0.0010 | +0.0016 | +0.0003 |
| ensP1 + dis100x0.25 | 0.2993 | +0.0017 | +0.0042 | +0.0080 | -0.0025 | +0.0008 | +0.0020 | +0.0001 |
| ensP1 + dis100x0.5 | 0.2977 | +0.0001 | +0.0062 | +0.0159 | -0.0026 | +0.0006 | -0.0022 | 0.0000 |
| ensP1 + dis100x0.5 + bmb | 0.3000 | +0.0024 | +0.0062 | +0.0059 | -0.0023 | +0.0007 | -0.0001 | +0.0002 |

**17.2 NEW OPERATOR `dis<k>x<lam>` (per-target hubness reduction), and it PASSES its safety test.**
arXiv 2511.16954 attributes PDS's difficulty to HUBNESS; the retrieval literature's fixes (local scaling, DisSim,
mutual proximity) are PER-POINT while our `ct` + guarded PC removal is GLOBAL. `dis` removes lam x each target's
projection onto the centroid of its own k most similar OTHER predictions. It buys PDS reliably (real scorer H1
+0.0042 at lam 0.25, +0.0062 at 0.5) and costs MSE, i.e. it is a frontier move; on KOLF/Orion the offline J is much
larger (+0.0305 / +0.0418 combined with the ensemble vs +0.0193 / +0.0396 for the ensemble alone).
SAFETY (`ml_dis_safety.py`, Codex's design): the feared failure is a final panel enriched for one pathway, where the
local centroid IS genuine shared biology. Built CLUSTERED vs DISPERSED 120-target pseudo-panels from DepMap
co-essentiality alone (no outcome information; within-set co-essentiality 0.0122 vs 0.0001 on H1). **The opposite of
the fear is true**: `dis` is MORE positive on clustered panels (J +0.0064/+0.0073 H1, +0.0044/+0.0057 KOLF,
+0.0078/+0.0144 Orion) and near zero on dispersed H1. Clustering creates the hubness `dis` removes. The removed
direction's cosine with the real response is +0.077 (H1 clustered), i.e. small.

**17.3 Mechanism of the ensemble, narrowed by three rejections**
- CORRECTED cell-count row weighting (`Cfg.ncw`, which keeps every row's `has` flag) is INERT: |J| <= 0.003 on 3/3
  panels. The earlier large "gains" were the denominator bug Codex found (`maximum(fden, 1)` clamping a fractional
  weight sum), i.e. pure amplitude shrinkage. 120 poisoned cache files quarantined in `src_ablation/skipped_ncw_buggy/`.
- COVERAGE-PATTERN offsets (`ml_strat_center.py`): group the targets by which families measure them (36/17/16 distinct
  signatures) and subtract each group's residual mean. Nil on 3/3 (J +0.0002 / +0.0035 / +0.0001), and adding it to the
  ensemble does not improve it. So the ensemble is not fixing a coverage artefact -- which is what would NOT have
  generalised to the final phase's differently covered panel.
- FAMILY-WEIGHT hedges (`ml_lofo_ens.py`): leave-one-family-out (H1 J +0.0022) and Dirichlet-perturbed weights
  (-0.0004) are far weaker than min_cells (+0.0112), and mixing them in DILUTES it.
=> What min_cells varies and the others do not is whether a (source, target) ROW is in or out at all. Codex's reading:
stabilisation of a discontinuous SELECTION estimator. Nothing "strictly better" is justified yet.

**17.4 CLOSED this loop (measured, do not retry)**
- Weak-target directional decorrelation (leaderboard rank 29's self-description): removing lam x the projection onto the
  across-target mean direction for the weakest q of targets is J -0.0014..-0.0126 on H1 for every (q, lam); our fixed
  `ct 0.5` is better and dropping `ct` costs -0.0231.
- Divergent-promoter strand prior (a mid-leaderboard entry's self-description): rebuilt the TSS table WITH strand from
  the GENCODE v47 GTF and measured the REAL neighbour log-FC on five panels by distance x orientation. Orientation does
  not separate (0-500 bp divergent medians -0.180/-0.112/-0.133/-0.136/-0.304 vs convergent -0.078/-0.117/-0.052/-0.066),
  and essentially every TSS pair within 5 kb is on opposite strands, so there is no tandem class to down-weight.
- Neighbour-similarity-source diversity (round 3 above).
- Literature checked, nothing new: STATE (zero-shot across datasets is its stated limitation), Stable-Shift (low-rank
  basis + GO/STRING graph convolution; our lr64/lr128 and neighbour-source tests already cover it), forrestsheldon's
  linear-response baseline (control covariance is chance -- we measured the same in section 8), TransPert (pseudobulk +
  Wilcoxon + PDS-targeted global scaling; 2025's PDS was L1, 2026's is cosine and scale-free per target).

## 18. Data acquisition round + destination-conditioned weighting (2026-09-29)

**18.1 The context-count curve is real but it is priced per COVERED TARGET, not per context.**
`ml_ctx_curve.py` gave +0.008-0.009 PDS for the 6th source family (H1 +0.0081, Orion +0.0092; KOLF -0.0019) with the
curve still rising, and `ml_fam_value.py` found family SIZE does not predict marginal value (Spearman -0.406 / +0.257 /
+0.429). That read was used to reverse §15 and acquire the dismissed public screens. It was HALF right: a source only
moves the targets it MEASURED, and the curve's families covered 200-800 of each panel.
Acquired and converted (`ml_newsrc_build.py`, `raw/scperturb_ctx/`, sources in `sources_new/`):
FrangiehIzar2021 (melanoma, 217 targets / **5 official**), SunshineHein2023 (lung, 133 / 5), XuCao2023 (HEK, 185 / 5),
TianKampmann 2019 x2 + 2021 (iPSC, neuron; 26-180 / 1-6), PapalexiSatija2021 + WesselsSatija2023 (THP-1, 25-27 / 1),
ShifrutMarson2018 (primary T, 21 / 1); ZhaoSims2021 and ChangYe2021 unusable; LiangWang2023 and CuiHacohen2023 mouse.
Measured at weight 1.0 in the full recipe (`ml_newfam_value.py`): every one of them is within +-0.005 J of zero on
H1/Orion/KOLF. **The small public screens are worthless for a 300-target panel** because they overlap it at chance
(~3.5 genes), and they will overlap the final phase's new panel at chance too.

**18.2 DLD-1 (GSE337988) is built and is a DESTINATION-SPECIFIC asset.**
`ml_dld1_pseudobulk.py` streams the 1,196,592 x 36,601 count matrix (metadata from the companion .rds: gene symbols in
rowRanges' elementMetadata, per-cell target in colData `assignment_v2`, 770,717 singlets, NTC controls 16,751) into
**4,714 targets, 67 of them official**, median 137 cells/target. Marginal value: **OR_HCT116 J +0.0427** (both
colorectal), H1 -0.0007, KOLF -0.0044. On the official `delta_s4` structure, adding it to the S11 sum side gives
Orion dPDS +0.0018 (+0.0036 at weight 2; +0.0026 on a no-KOLF base) and H1 ~0 -> roughly +0.001 composite for A/B/C,
but it is a real asset if a final-phase line is GI/epithelial.

**18.3 WHICH PANEL TO TRUST for a source-family change: the somatic one, not H1.**
Re-measuring the S10k no-KOLF sum side on the official structure: H1 J **-0.0114**, Orion J **+0.0368**. That change
scored **+0.0061 officially**. So for stem-vs-somatic source composition, H1 (a stem line) gives the WRONG sign and the
somatic panels give the right one. A/B/C are Jurkat / HeLa / CAL-33, all cancer lines. Use Orion + RN for source
composition; keep H1 primary for emission and metric-profile questions.

**18.4 CLOSED: destination-conditioned source weighting from control expression.**
`ml_dest_weights.py` built the full (family x panel) marginal-value table over seven panels and asked whether the
cosine between the destination's control profile and each source line's control profile ranks it. Pooled within-panel
Spearman **+0.189** (per panel -0.140..+0.392). Too weak to set weights from. Section 12's oracle stands unreached, and
the destination side is now closed as well as the source side.

**18.5 THREE VERIFICATION BUGS found this round, all of the silent kind.**
(a) `tv2_model` clamped the family-mean denominator at 1, so fractional `ncw` weights became amplitude shrinkage
(Codex); (b) the offline evaluator matched candidates to the BASE budget while the real scorer used the candidate's own,
biasing MSE by a constant +0.0103 (Codex); (c) `delta_s4` keys its cache on the panel alone, so weight overrides
silently returned the cached base delta and all eight family arms read EXACTLY 0.0000. Rule: a result that is exactly
zero, or exactly constant across arms, is a bug report until proven otherwise.

## 19. Per-target source weighting: the fourth and strongest attempt, and it fails (2026-09-29)

`ml_learn_weights.py`. Section 12's per-target NNLS oracle (+0.048 PDS H1 / +0.130 Orion) had defeated three predictors
(MODZ agreement, an LPO-GBM on test-time features at -17 % of the oracle, and the destination control cosine at
Spearman +0.189). All three shared two weaknesses, and this attempt removed both: the LABEL is the oracle WEIGHT VECTOR
itself rather than a performance score, and the model is trained on the SEVEN PANELS' OWN TRUTHS (everything we had ever
fitted -- S4's refit, the FiLM residual, the discriminative residual -- was trained on the source table instead).
~12,000 (panel, target, family) rows; features all test-time (measured flag, the family's delta norm for that target,
its agreement with the other families, how many families measure it, the shipped weight, the target's expression in the
DESTINATION's controls, gene count, family one-hot); HistGradientBoosting; leave-one-panel-out.

**On the fusion alone it works**, and it reverses the sign of the previous best predictor: +13..20 % of the oracle
recovered at lam 0.5 on 5/7 panels (H1 +0.0145 PDS, RN_RPE1 +0.0116, RN_K562 +0.0100, RN_HepG2 +0.0096), versus
-17 % for the LPO-GBM. `--strict` (dropping training rows whose target appears in the held-out panel, so the model
cannot memorise "for gene X use K562") does NOT hurt -- H1 goes +0.0132 -> +0.0145.

**It dies when the reweighted fusion is put back into the shipped recipe.** Full recipe (0.5 x 3 x S4 + 0.5 x 1.6 x
fusion + blend 0.3), strict, lam 0.5: H1 **-0.0065**, Orion **-0.0022**, RN_Jurkat -0.0092, RN_RPE1 -0.0030,
RN_K562 -0.0001, KOLF +0.0027, RN_HepG2 +0.0072 -- mean **-0.0016**. The pre-declared bar (H1 and Orion both positive,
>= +0.003 composite) fails.
**Why:** S4 is the other half of the sum and already carries most of what the reweighting adds -- the ORACLE itself
halves, from +0.0765 on the fusion alone to +0.0372 in the full recipe on H1 -- so reweighting mainly destroys the way
the fusion COMPLEMENTS S4. A fusion-only measurement of any source-side change is therefore not evidence for a build.

=> Per-target source weighting is closed. The oracle headroom is real and still unreached (full recipe: H1 +0.0372,
Orion +0.1132); four predictor families have now failed on it.

**AMENDMENT 2026-09-29 (Codex round 11, `debate_loop/codex_round11.md`, defect A): the full-recipe numbers above are
INVALID.** `score_W(full=True)` injected the reweighted fusion as `L0 + 0.8 (F' - F)` after rescaling both to the full
recipe's per-target norm, but `delta_s4` applies the 0.3 neighbour blend and the global energy match AFTER mixing the
parts, so a part cannot be swapped post hoc; it also zero-filled unmeasured genes instead of renormalising per gene and
omitted the RN families. The "closed" verdict is therefore WITHDRAWN as unproven. It is not reopened as a priority
either: section 20's ceiling analysis bounds the whole per-target weighting headroom on H1's strong responders at
+0.06..+0.12 cosine (oracle NNLS mix, half-gene held out), of which the learned weights recovered 13-20 % on the fusion,
i.e. an expected ~+0.001..+0.003 composite. A faithful rerun must first reproduce `delta_s4` with shipped weights.

## 20. Where is the bottleneck, from every angle? (2026-09-29, user: "다각도로 병목이 어딘지 ... 크게 개선할수있는 부분")

**20.1 The gap.** Leaderboard 2026-09-29 (`reports/leaderboard_2026-09-29.json`): [top team] 0.3974 (PDS 0.890, MSE 0.399,
NMAE 0.829, fid 0.560, reach 0.423, jac 0.053), then 0.328 / 0.306 / 0.291 / 0.280; we are 0.2260, rank 40. Leader minus us
= 0.171 composite: MSE 0.062 (36 %), reach 0.028, PDS 0.026, NMAE 0.024, fid 0.021, jac 0.011 -- the four DE members are
HALF of it. [top team] already had 0.330 / reach 0.40 on 09-18; its recent jump is MSE 0.62 -> 0.40 in two days. Implied
replicate (score 1) raw values: MSE 0.04, PDS 0.95, reach 0.96, fid 0.80, NMAE 0.40, jac 0.40; baseline (score 0): MSE
0.987, PDS 0.505, reach 0.08, fid 0.513, NMAE 0.999, jac 0.031 -- we sit AT the baseline on fid and jac.
No shared public recipe above AtlasShift (0.155) is visible at 5-decimal precision.

**20.2 The MSE leg is decided by the strong responders** (`ml_err_decomp.py`, S11 sum side, offline evaluator reproduces
`SE.evaluate` exactly). H1 (MSE 0.7543, the official value is 0.7529): the 60 strongest targets hold 75 % of the true energy
(top-10 28 %, top-30 54 %); on them our cosine is 0.29 and our norm 0.34x the truth; the weakest quintile holds 0.4 %.
What-ifs: direction oracle on the top 30 (our norms) -0.243 MSE, truth on the top 30 -0.594, truth on the whole weak
half -0.388 MSE / +0.083 PDS; per-target MSE-optimal amplitude on our own directions is WORSE (+0.027) because it
shrinks B_pred and with it the #348 rebate. KOLF and Orion truths are noise-dominated (jackknife noise >= energy for
most targets; global cosine 0.02 / 0.05), so neither can price this.

**20.3 The ceiling of source-based prediction for strong responders** (`ml_strong_ceiling.py`, H1): top-10 / 30 / 60
strongest targets -- ours 0.363 / 0.309 / 0.280, equal family mean 0.364 / 0.296 / 0.257, best single family 0.426 / 0.324 /
0.284, per-target NNLS mix fitted on half the genes and scored on the other half 0.481 / 0.386 / 0.344, noise ceiling of
the truth 0.904 / 0.831 / 0.801. So even perfect per-target weighting of every source we own stops at ~0.35-0.48 against a
truth that allows 0.8-0.9: **we are at ~80 % of what our public data can express on H1**, and the leader's MSE implies about
twice our direction quality. Per-family, stem sources are best for this stem panel (FENG 0.336, KOLF 0.294 on the top 10).
Skipping ct/pcr for the top-k PREDICTED-strong targets: H1 J +0.0036 / +0.0042 / +0.0024 (k 10/30/60), Orion dPDS +0.0007 /
+0.0004 (MSE unreliable) -- small, a ride-along at most.

**20.4 The DE side is lost on strong targets** (`ml_de_strata.py`, cached H1 tables, replay reproduces fid 0.5405 / jac
0.0695): targets with <= 20 real DE genes (93 of 300) reach 0.51-0.69 with 0.80 |lfc|-weighted sign accuracy; targets with
> 100 (161) reach 0.08-0.12 at ~0.65 and hold ~70 % of the reach shortfall, although our top-1 gene is right 94 % of the time.
**Re-ranking by source agreement is negative everywhere** (`ml_conf_rank.py`, 7 panels): |fused| x agree^2 H1 -0.0014,
agree^4 -0.0069, t-stat -0.0078, count of agreeing families -0.0135; RN_K562 -0.009..-0.028; Orion -0.001..-0.007. The fused
magnitude is already the best confidence ranking.

**20.5 Closed today.**
- Other discontinuous boundaries (`ml_boundary_ens.py`): clip {2,3,5,10} exactly 0; NB_TOPK {5,10,20,40} H1 J +0.0005,
  Orion +0.0047, KOLF -0.0001.
- LINCS L1000 shRNA as extra contexts (`ml_lincs_value.py`, GSE106127 CGS, 9 lines, 98-116 official targets each): cosine
  with the truth on landmark genes H1 0.007-0.026 (ours 0.035-0.113), strong targets 0.005-0.065 (ours 0.11-0.24); LINCS HEPG2
  is not even the best LINCS line on RN_HepG2 -> generic / off-target signal, not target biology.
- Lineage source for Jurkat: CD4 is NEGATIVE on RN_Jurkat (drop-CD4 dPDS +0.0068); Jurkat's controls resemble K562 (0.847)
  more than CD4 (0.519).
- Same-line OTHER-target transform (`ml_line_transform{,2,3}.py`): a K562 -> line map learned on the destination's own other
  targets lifts RN strong-target cosine 0.30 -> 0.48 (Jurkat), 0.28 -> 0.71 (HepG2), 0.15 -> 0.68 (RPE1) -- but that is the
  shared ESSENTIAL program: learned on the K562 essential screen and applied to the OFFICIAL targets as measured in K562 GWPS
  it LOWERS PDS for every source (HCT116 0.767 -> gw 0.754 / mix30 0.700; HEK293T 0.762 -> 0.738 / 0.707; CD4 0.721 -> 0.659 /
  0.619). Not a lever for a regulator-heavy panel.
- Target silent in one context (`target_expr_ABC`): only 4 targets are silent in one context (S100A11 A, VIM C, KIF21B B/C,
  PSMB9 B), 18 differ > 10x -> ~+0.001.
- `pert_counts.csv` carries target names only (every target has exactly 400 real cells) -> no fitness signal.

**20.6 Codex round 11 audit verdicts** (`debate_loop/codex_round11.md`): A DEFECT (section 19 amended above); B DEFECT --
the final-phase coverage repricing (-0.002) holds only for "fusion missing, S4 measured" targets: S4, the source centres and
the library PCs still carry the deleted targets, the limited bucket was never built, and the reported "composite" was
dPDS x 2.234/6, not the six-member score; C no defect in the K sweep; D the K562-ONLY arm keeps S4 and the RN families, and
x2/x4 also move K562 in the equal-weight library component (16 vs 1) -- confounded; E sampled H5 row alignment passes (14
rows, totals match) but DLD-1 is normalised by MAPPED-gene totals, unlike every other source -- audit before relying on it.

=> The pipeline is near the ceiling of its data on the one panel that predicts the official MSE. The remaining large lever
is data that is closer to the destinations (same-lab or same-line), not a better use of what we have.

**20.7 CORRECTION (Codex round 12, `debate_loop/codex_round12.md`) -- the "~80 % of the data ceiling" claim in 20.3 is
WITHDRAWN until re-measured.** (1a) The half-gene per-target NNLS is a RESTRICTED-MODEL oracle (a nonnegative cone over
family means, missing coordinates zero-filled, families already averaged and cut at 20 cells), not an information ceiling:
it excludes signed, program-specific and cross-target combinations. (1b) The comparison mixed spaces: family / NNLS
cosines were in the c=50 natural-log LFC space of `r["R"]`, "ours" in scorer space (bulk_lognorm at TS 5e4 = c 20 CPM),
and the "noise ceiling" subtracted SCORER-space jackknife energy from c=50 energy -- invalid. (1c) The transform closure is
weaker than stated: genes were selected by essential-screen variance, one ridge penalty, no intercept, and the PDS mask
removed only each row's own gene. `ml_de_strata.py` ranks calls by magnitude while the scorer ranks by significance,
p_adj, p, magnitude, feature -- its fid/jac replay is exact but the reach attribution is NOT validated. The err_decomp
truth substitutions are algebraic sensitivities (they keep the jackknife rebate), not attainable gains.
Codex's alternative explanation for #2-#5 (a HYPOTHESIS): information is lost by fusing sources BEFORE learning -- the
discriminative residual sees only the fused response on 64 PCs, so a coherent source-specific program and several
contradictory ones with the same weighted mean are indistinguishable. Its proposals: a residual that keeps per-family
signatures + availability masks (expected ~+0.003, range 0..+0.010) and program-wise family reliability (+0.001..+0.004).
Re-measurement in one space: `ml_ceiling_v2.py` (running).

**20.8 The ceiling, re-measured in ONE space** (`ml_ceiling_v2.py`, `err_decomp/ceiling_v2.json`): every candidate is
applied to the destination controls exactly as the sums are built and pushed through `bulk_lognorm_means` (TS 5e4); same
genes (all panel target genes removed), same targets, noise from the same space's jackknife. Strata by the scorer's
noise-corrected distance. Mean cosine with the truth, top-10 / top-30 / top-60:

| panel | noise ceiling | per-target NNLS families+S4, half genes (in-sample all genes) | S4 alone | OURS S11 sum side |
|---|---|---|---|---|
| H1 | 0.954 / 0.948 / 0.923 | 0.439 / 0.408 / 0.364 (0.443 / 0.406 / 0.361) | 0.394 / 0.349 / 0.300 | 0.356 / 0.326 / 0.290 |
| RN_Jurkat | 0.786 / 0.763 / 0.722 | 0.466 / 0.423 / 0.378 (0.470 / 0.429 / 0.382) | 0.384 / 0.318 / 0.265 | 0.430 / 0.368 / 0.318 |
| RN_HepG2 | 0.943 / 0.906 / 0.866 | 0.550 / 0.504 / 0.465 (0.554 / 0.511 / 0.473) | 0.458 / 0.400 / 0.349 | 0.483 / 0.433 / 0.377 |

Reading: (i) we sit at 80-91 % of the best per-target nonnegative mix of our own sources, and that bound is tight (the
honest half-gene and the in-sample all-gene fits agree to 0.01); (ii) that cone reaches only ~45-60 % of the noise
ceiling on the strong targets, so most of the strong responses are NOT expressible as any per-target weighting of our
per-target source measurements -- a restricted-model statement, not an information ceiling (Codex R12 1a stands);
(iii) on H1 S4 ALONE beats the shipped mix on strong targets (0.394 vs 0.356) -- the fusion / blend / operators dilute
it; on RN_HepG2 it does not (0.458 vs 0.483). Tested next: a strength-adaptive S4 share (`ml_s4share.py`).

**20.9 Closed after the re-measurement.**
- **Strength-adaptive S4 share** (`ml_s4share.py`, parity with `delta_s4` to float32 rounding): a larger S4 share for the
  top-k PREDICTED-strong targets is H1-only. J for top-30/60 at share 0.75 / 1.0: H1 +0.0084 / +0.0124, +0.0102 / +0.0128;
  RN_Jurkat -0.029 / -0.072, -0.033 / -0.077; RN_HepG2 -0.020 / -0.047, -0.014 / -0.047; Orion +0.0017 / -0.0038, +0.0045 /
  -0.0093. Global shares 0.6 / 0.75 lose everywhere. S4 beating the mix on strong targets is a property of the STEM panel.
- **GSE249595 as a SAME-LINE source for context A (Jurkat)** (`ml_jurkat_signal.py`): it covers 299 / 300 official targets
  (STIM, median 499 cells) and the on-target knockdown is visible (median -0.31, < -0.2 in 62 % of targets), but its
  per-target responses carry no reproducible signal: STIM vs UNSTIM (independent cells) cosine +0.001, vs every other line
  +0.001..+0.004, where K562 vs HCT116 is +0.031 and CD4 Stim8hr vs Rest +0.324 on the same targets and genes. The high-MOI
  marginal estimator (~12 guides per cell) buries the target effect under the other guides. Its value for Jurkat had never
  been measured (the 09-24 test was transfer to OTHER lines); it is now measured, and it is noise.

**20.10 New data: GSE295214, genome-scale CRISPRi Perturb-seq in H1 hESC** (Weissman/Zhuang lab, public 2026-05-05, CC-BY
preprint doi 10.64898/2026.05.16.725005; found by the recon agent). `raw/h1_gw_gse295214/` (85 Cell Ranger h5, 7.8 GB),
built by `scripts/ml_h1gw_build.py` into `reports/transfer_v2_2026-09-23/sources_new/H1GW.npz` (+ `.json`): 1,521,921 cells,
1,025,672 uniquely assigned (top target >= 5 guide UMIs, runner-up <= 20 % of it; dual Replogle guides collapsed per
target), 81,617 non-targeting controls, 2,747 targets with >= 20 cells (median 284), on-target knockdown median 62 %
(> 50 % in 63 % of targets), 66 of the 300 official targets. Normalised by FULL-library totals (unlike DLD-1, Codex R11 E).
10x 5' chemistry, sequenced on Ultima UG100 (the organisers' 2026 sequencer). It is a STEM line: expected <= 0 for the
somatic A/B/C (same finding as KOLF), a first-choice source for a stem-like final line, and -- on the 131 VCC2025-H1 targets
it shares -- a same-line, other-lab measurement of Arc's H1 truth. **Never a source in the H1 harness** (leak).

**20.11 What same-lineage data is worth** (`ml_h1gw_value.py`, H1 panel, scorer space; H1GW covers 131 / 300 panel
targets). Cosine with Arc's H1 truth over the COVERED targets, top-10 / top-30 / top-60 by strength: noise ceiling 0.953 /
0.932 / 0.905; H1GW ALONE 0.482 / 0.492 / 0.412; best per-target mix of our sources + H1GW 0.547 / 0.553 / 0.478; best mix of
our sources only 0.445 / 0.419 / 0.350; S4 0.386 / 0.360 / 0.290; ours (S11 sum side) 0.347 / 0.335 / 0.276.
Blending H1GW into the covered targets' sum-side delta (energy-matched, offline evaluator, shipped operators and budget):
w 0.25 dPDS +0.0076 / dMSE -0.0505 / J +0.070; w 0.5 +0.0099 / -0.0873 / J +0.114; **w 1.0 +0.0060 / -0.0993 / J +0.118
(~ +0.020 composite from 44 % coverage, sum side only)**.
Reading: ONE same-line screen from another lab and chemistry (10x 5' vs Flex) lifts strong-target direction ~40 % over
everything we own, i.e. an order of magnitude more than any within-data lever measured this week (+0.001..+0.003), yet it
still stops at ~58 % of the noise ceiling (lab / chemistry / protocol). This is the quantitative basis for the final-phase
strategy: once D/E/F are named, a public LOW-MOI same-line screen for any of them is worth up to several hundredths.
(Diagnostic only: H1GW never becomes a source on the H1 harness.)

**20.12 Codex round 13 audit of the DE-side experiment** (`debate_loop/codex_round13.md`): OLD vs NEW differ only in
the DE delta (sum side, budget reference, overrides and ranking score shared) -- NO DEFECT. OLD is not an exact official
mirror (harness agreement is relative to the recipe's own sign, the official file to the unweighted family-mean sign;
overrides from the new fusion; no `sum.resid`) -- shared by both arms. Scope: the harness freezes the within-top-K
STRENGTH ORDER to the old recipe's |d| x agreement, while the builder would re-rank by |new DE delta| x agree file; top-K
membership and signs follow the new delta in both. A positive NEW arm supports "new DE LFCs on the old strength order"; the
builder change is: snapshot the override recipe before `sum.resid` for `Dde`, keep the old `D` as the budget reference and
the old `OV`, and either accept the re-ranking or add a frozen ranking-score input.

**20.13 Final-phase consequence: the same-line policy in `tv2_final.py` is too weak (design, not yet implemented).**
Today `context_weights` multiplies a same-line family's FUSION weight by `SAME_LINE_BOOST = 3` when the HPA identity rho is
>= 0.85 and the knockdown gate passes. With the base weights (K562 8 ... HCT116 1, HEK293T 0.5) that gives a same-line
family ~20-25 % of the fusion and nothing in S4 or the neighbour part; a zero base weight (JURKAT) stays zero; and the
sources added since (H1GW, DLD1) are not in `SOURCE_LINES`. The one valid measurement (20.11, H1GW on H1, other lab and
chemistry) found the best use to be near-FULL replacement of the covered targets' delta, energy-matched per target
(w 1.0 J +0.118 > w 0.5 +0.114 > w 0.25 +0.070). Proposed change for 10/22: per context with a same-line LOW-MOI source,
L_t <- (1 - w) L_t + w * S_t * |L_t| / |S_t| on the targets that source covers (>= 20 cells), w ~ 0.75-1.0, applied to BOTH
the main (DE) and the sum-side recipes before the sum operators; keep the KD gate and add a split-half reliability gate
(GSE249595 shows why: KD visible, responses noise). Add H1GW (stem), DLD1 (colorectal) and the Orion lines to SOURCE_LINES.
A second valid calibration point needs a same-line, other-EXPERIMENT pair: the RN_K562 panel cannot provide it (it is
sampled from the K562 GWPS itself, so the harness K562 store is the same experiment -- a leak; that test was stopped).

**20.14 The DE side was left on the pre-S10k recipe -- aligning it is the one positive, consistent lever found today.**
Official S11 builds its DE side (NMAE / fid / reach / jac, 4 of 6 members) from the MAIN recipe: S4 with KOLF, fusion
with KOLF 0.25, neighbour library with KOLF, no min_cells ensemble. Only the SUM side got the validated no-KOLF sources and
the ensemble. Real scorer, harness options `dk0` / `deo` / `bro` in `sch_emit.py` (`ml_deside_queue.sh`; Codex R13: the arms
differ only in the DE delta), sum side and budget reference identical by construction:

| panel | OLD (official S11 mirror) | NEW (DE side = S11 sum-side recipe) | d composite | d nmae | d fid | d reach | d jac |
|---|---|---|---|---|---|---|---|
| OR_HCT116 | 0.4151 | 0.4205 | **+0.0054** | -0.0065 | +0.0009 | +0.0164 | +0.0001 |
| RN_Jurkat | 0.4819 | 0.4839 | +0.0020 | -0.0010 | +0.0022 | +0.0015 | +0.0003 |
| RN_HepG2 | 0.4565 | 0.4583 | +0.0018 | +0.0018 | +0.0002 | +0.0109 | +0.0002 |
| RN_RPE1 | 0.3613 | 0.3646 | +0.0033 | -0.0003 | +0.0020 | +0.0086 | +0.0010 |
| RN_K562 | 0.4200 | 0.4225 | +0.0025 | -0.0023 | +0.0017 | +0.0041 | +0.0003 |

5 / 5 somatic panels positive, mean **+0.0030**; reach, fid and jac up on all five, NMAE better on four. H1 (stem; only the
ensemble half applies, KOLF kept): the existing ensembled-DE arms read +0.0014 (paired seed) / +0.0022 (seed 1) over the
OLD mirror, mostly NMAE (-0.003..-0.005).
The queued H1 pair itself: OLD 0.3015 -> NEW 0.3029, **+0.0014** (NMAE -0.0035, fid +0.0008, reach -0.0001, jac +0.0001)
-- 6 / 6 panels positive, mean +0.0027. **This also accounts for the S11 transfer gap:** harness S11 (ens + ci5) predicted
+0.0044 on H1 while official gave +0.0030; the harness DE side already carried the ensemble, the official one did not, and
the difference (~0.0014) is exactly the DE-side ensemble effect measured here.
Builder: `tv2_submit_build.py` option `de.from_sum_override` (flag-guarded; snapshot of the override recipe BEFORE
`sum.resid`; budget reference, overrides and agree file stay the main recipe's). Regression: the patched builder
reproduces S11's shipped main delta A to 1.2e-7. `spec_S12.json` = S11 + that flag. The new DE delta vs the old: per-target
cosine median 0.940 (p10 0.873), norm ratio 1.04, sign agreement on each target's old top-3000 genes 0.997 -- a modest
refinement, not a new regime. Follow-up arms with the within-top-K ranking score also taken from the new delta (`brs`, the
exact builder equivalent) in `ml_deside_queue2.sh`.
**20.14b The exact builder equivalent (`brs`: within-top-K ranking score ALSO from the new delta, as S12's builder does)**
holds the gain: OR_HCT116 +0.0046, RN_Jurkat +0.0015, RN_HepG2 +0.0016, RN_RPE1 +0.0033, RN_K562 +0.0023, H1 +0.0014 ->
**6 / 6, mean +0.0024** (somatic mean +0.0027, vs +0.0030 for the frozen-ranking `bro` arms). **S12 built** (`reports/submission_2026-09-23/S12.vcc`, sha256
61dd0011cbc193095e01504e9b05e3549378f7fe61c743f9270cc9d5eec20161, vcc prep passed: 360,000 cells, 18,533 genes, targets
verified, counts preserved); its main delta, sum-side delta and overrides are bit-identical to S11's (context A checked),
only the DE delta differs. NOT uploaded (needs the user's go-ahead; expected official ~+0.002..+0.005).
**20.14c Codex round 14 pre-upload audit of the S12 builder patch** (`debate_loop/codex_round14.md`, diff
`tv2_submit_build_S12_patch.diff`): NO DEFECT on the snapshot (pre-`sum.resid`, reset per context; `spec_override` still
supersedes `parts_w`), on `lfc_de = g0 * Dde` being reached, and on flag-off builds (algebraically unchanged). Pinned sums --
hence PDS -- are provably unchanged; the cpc target, top-K membership, signs and agreement-weighted order intentionally follow
the new DE delta. Open point: MSE's rebate uses the emitted cells' jackknife, which the DE-side re-allocation could move unless
the cap binds -- empirically it does not move: |dMSE| < 5e-5 in all 6 real-scorer OLD/NEW pairs (same emitter).

**20.15 S4 variant ENSEMBLE: closed** (`ml_s4ens.py`; the one operator that transferred officially is an ensemble, so the
S4 half was ensembled over its saved LOBO variants, each energy-matched per target, budget fixed to the shipped recipe).
J by pool (rank r16/64/128 | alpha a0.5/a1 | qz q1/q2/z1/z4 | rank+alpha | all):
H1 +0.0008 | n/a | +0.0026 | +0.0008 | +0.0019;  RN_Jurkat +0.0012 | +0.0260 | -0.0054 | +0.0148 | +0.0069;
RN_HepG2 -0.0016 | -0.0059 | -0.0087 | -0.0020 | -0.0051;  RN_RPE1 -0.0002 | +0.0064 | -0.0044 | +0.0031 | -0.0002;
RN_K562 +0.0006 | +0.0086 | -0.0017 | +0.0060 | +0.0040;  KOLF +0.0011 | -0.0146 | -0.0026 | -0.0054 | -0.0047.
The rank pool is ~0 everywhere (the variants are too similar to diversify); the alpha pool is really a SHIFT to alpha ~0.58
(FiLM adds linearly) and repeats the 09-27 alpha verdict (helps RN, loses KOLF / HepG2). Nothing is positive on every panel.

**20.16 Per-target source weighting, FAITHFUL rerun (loop #1, `ml_learn_weights2.py --strict`): closed with a valid
measurement.** `tv2_model.predict_lfc` now accepts per-target family-weight vectors (scalar weights bit-identical); the
learned weights are injected INTO the fusion and the recipe is rebuilt through `delta_s4`'s own path (parity with
`delta_s4` 5e-9..1e-8 on all 7 panels). The full-recipe ORACLE is large (J H1 +0.125, RN +0.09..+0.15, KOLF +0.42), but the
LOPO-learned weights (lam 0.25 / 0.5 / 1.0): H1 +0.004 / +0.005 / -0.006 (PDS -0.003..-0.014), OR_HCT116 +0.027 / +0.046 /
+0.038 (MSE only, PDS ~0), KOLF +0.003 / +0.010 / +0.023, RN_K562 -0.009 / -0.022 / -0.066, RN_Jurkat -0.011 / -0.032 /
-0.098, RN_HepG2 -0.009 / -0.025 / -0.097, RN_RPE1 -0.003 / -0.007 / -0.035. Negative on all four RN panels and never a PDS
gain on H1: the oracle is real and unpredictable. Section 19's closure stands, now on a valid measurement.

**20.17 Family-signature residual (loop #2, Codex R12 idea 1; `ml_dr_train_fam.py`, `ml_dr_eval_fam.py`, `ml_dr_eval_mix.py`;
offline, S10i sums, alpha scan; outputs `disc_residual/famsig_20260929/`).** Identity-aware per-family inputs OVERFIT: the
correction grows to 1.2-1.6x |X| and PDS falls (H1 -0.019, Orion -0.019 at alpha 0.25). Shuffling the family slots per row
helps H1 / Orion -- the set of source responses matters, their identity does not -- so a permutation-invariant encoding
(mean + std over available families + count; `famset`) was built: paired vs the shipped fused-input residual at alpha
0.25 (3 seeds) it is mixed, H1 +0.006, OR_HCT116 +0.001 (PDS +0.004), KOLF -0.005, RN_Jurkat -0.007, RN_HepG2 +0.011,
RN_RPE1 0.000, RN_K562 0.000. The AVERAGE of the two residuals is >= the shipped one on 5 / 7 panels (15 / 21 seeds): H1
+0.004, Orion +0.010, KOLF -0.000, RN_Jurkat -0.001, RN_HepG2 +0.007, RN_RPE1 +0.001, RN_K562 +0.002 -> mean +0.003 J
(~ +0.0005 composite). Real but far below the +0.005 bar: a ride-along only if the residual is ever retrained.

**20.18 Strength-aware operator skip (loop #4): closed.** Skipping ct 0.5 / guarded pcr for the top-30 PREDICTED-strong
targets: H1 J +0.004 (PDS -0.002), Orion +0.056 (MSE-driven, PDS +0.000), RN_Jurkat +0.006, RN_HepG2 +0.012, RN_RPE1
+0.007, RN_K562 +0.002 -- but KOLF -0.011 (operators are worth 0.08 there). ~+0.001 composite and not all-panel positive.

**20.19 Final-phase port of the S12 structure (task #23; Codex R15 / R16 fixes; 2026-09-29 17:00-18:10Z).** `tv2_final.py`
now emits, per context, the S10k/S11/S12 structure: `--sum-policy auto` (stem score from pinned pluripotency vs somatic
marker panels, threshold 1.5; somatic -> no-KOLF S4 via `tv2_s4_nokolf_final.py` + no-KOLF neighbour library via
`tv2_final_knn.py ... KOLF`), `--sum-ens` (min_cells ensemble), `--de-follow-sum` (S12), `--sum-resid 0.25` (new
`tv2_final_resid.py`: prep + `ml_dr_train.py noesm_strict` x 3 seeds + seed mean; KOLF excluded exactly where the policy
is somatic), and `--dual ...,ci5,lsb`. Parity on A/B/C:
- `tv2_final_resid.py` reproduces `resid_S10k.npz` BIT-EXACTLY (input X max diff 0, residual max diff 0 for every seed).
- The generated spec equals `spec_S12.json` except the per-panel KNN file; agree file and library PCs rebuild bit-identically;
  S4 / no-KOLF S4 are bit-exact (earlier). The KNN rebuilt for a panel is cosine 0.90 (median) to the shipped one because the
  shipped KNN drew incidental neighbours from the 2,581-target harness library (5.45 of the top-10 neighbours present on
  average vs 5.02 in the per-panel library); a top-10 library (Codex R15) would be a different, unvalidated estimator, so the
  top-3 construction is kept and documented.
- Builder audit (`audit_only` spec key: deltas / overrides / sums / DE side / kink without emission, in-memory HDF5) with the
  shipped KNN files swapped in: **PASS** (19:15Z) -- `_delta_`, `_override_`, `_dde_` bit-identical to S12 on A/B/C, `_dsum_`
  max |diff| 1.2e-7 / 6e-8 / 3e-8 (float32), kink_base / gamma_sum 0.893 / 0.984, 0.937 / 0.987, 0.960 / 0.982 = S12's log.
Defects fixed on Codex review: legacy CD4 x2 rule for T-lymphoid identities removed (S12 has no boosts; dropping CD4 on
RN_Jurkat was +0.0068 PDS); the x3 same-line fusion boost is now opt-in (`--identity-boost`; never validated in the S12
structure); same-line replacement (builder key `same_line`, flag `--same-line W`) keeps baseline coordinates the source does
not measure, excludes override / own-gene coordinates from the norm match, clips w to [0,1], changes only covered rows of the
DE delta (`Dsl`; was: the whole override recipe), requires an identity margin >= 0.03 (or a runner-up mapping to the same
source) and a MEASURED passing knockdown (DLD1 / H1GW now go through `TM._src_path`); HEK293 -> HEK293T and HCT 15 -> DLD1
are cross-label substitutions at half weight; stem scores are computed whenever `--same-line` is used; a reused
`resid_final.npz` must match the panel's targets and per-context exclusions or it is rebuilt. Rehearsal 5
(`tv2_final_sim5.sh`: new 300-target panel, D/E/F = A/B/C controls relabelled, G = H1 non-targeting cells on the official axis
= a STEM destination) is queued behind the audit to exercise every branch (stem policy, H1GW fallback, residual with KOLF).

**20.20 SAME-LINE DATA IS WORTH ~+0.074 COMPOSITE ON H1, 70 % OF IT ON THE DE SIDE (loop #6; real scorer; H1 panel =
Arc VCC2025 H1 truth; source H1GW = GSE295214, other lab, 10x 5', covers 131 / 300 targets with >= 20 cells; DIAGNOSTIC --
H1GW is never a source for H1 decisions).** Harness tokens `sl<SRC>x<w>` (builder `same_line` exactly: covered rows only,
genes the source measures, override / own gene excluded from the norm match, applied to the recipe delta so sums AND DE side
change; budget references untouched) and `sls<SRC>x<w>` (identical sums, DE delta AND agreement ranking kept at the
pre-replacement delta). S12-structure base (`..._s4Yorig_ensP1_bmb_ci5_lsb_brsYorig`), seed 0:
| arm | composite | PDS | MSE | NMAE | fid | reach | jac |
|---|---|---|---|---|---|---|---|
| base | 0.3029 | 0.9076 | 0.7438 | 0.8622 | 0.5427 | 0.2732 | 0.0711 |
| sls (sums only), w 1 | 0.3249 (+0.022) | 0.9137 | 0.6405 | 0.8578 | 0.5419 | 0.2752 | 0.0717 |
| sl (sums + DE side), w 1 | **0.3772 (+0.074)** | 0.9137 | 0.6405 | 0.8148 | 0.5647 | **0.4035** | 0.0796 |
The sums-only arm reproduces the offline estimate (J +0.122 real vs +0.117 offline) and its PDS / MSE equal the full arm's to
the last digit, so the DE-side share is exactly +0.052 (reach +0.128, NMAE -0.043, fid +0.023, jac +0.008). Not leakage:
H1GW's strong-target cosine with the truth is 0.49 against a noise ceiling of 0.93, i.e. an independent experiment.
Offline amplitude rules (`ml_h1gw_strength.py`, sums only): same-line STRENGTH ranks the truth no better than ours (Spearman
+0.520 vs +0.517 on covered targets; top-30 strongest +0.08 vs +0.22), and taking H1GW's relative amplitudes instead of our
per-target norm adds only J +0.004..+0.008 (w 0.5: pt +0.113 / geo +0.115 / glob +0.117; w 1: +0.117 / +0.123 / +0.125) ->
keep the per-target norm match; the value is the source's DIRECTION, and it pays most in the DE metrics.
CONSEQUENCES: (1) the final pipeline's same-line replacement must act on BOTH sides (it does: `Dsl` / `Dover`); (2) the
final-phase data search -- a low-MOI same-line screen for each of D/E/F -- is priced at ~+0.07 composite per line at 44 %
coverage, 3-10x any other lever ever measured here; (3) queued: w 0.5 on H1, DLD1 -> OR_HCT116 (cross-label related line),
H1GW -> KOLF (the stem fallback on a different pluripotent line), confirmation seeds.
Half weight (w 0.5, both sides, seed 0): 0.3798 (+0.077; PDS +0.010, MSE -0.091, reach +0.127) -- w 0.5..1 is a plateau
(PDS favours 0.5, MSE 1.0; the DE-side gain is saturated at 0.5).

**20.21 Same-line data search for the VALIDATION lines (2026-09-29 19:10Z, background agent, verified accessions): NOTHING
NEW.** Jurkat: X-Atlas/Pisces (Xaira; Jurkat resting + CD3/CD28, CRISPRi dual guide, 10,872 targets, ~260 cells/target,
median KD 79 % / 71 %) is still UNRELEASED (HF repo = README + LICENSE, last change 2026-03-25; no GEO / Figshare / Zenodo
copy) -- the one external event that would move context A (watcher `pisces_watch.sh` stays on); other Jurkat screens are
KO with ~10-31 targets. HeLa: no public transcriptomic perturbation screen (imaging only). CAL-33 / HNSCC: no single-cell
CRISPR screen in any HNSCC line (bulk fitness only); Ker-CT keratinocytes still on request. Broad/DepMap multi-line
Perturb-seq (figshare 10.6084/m9.figshare.33273600; cervix SKGII / C4I, H&N T3M5, esophageal KYSE140/450, ...) has 100
targets and 0 / 300 official -> a lineage-matched validation panel at most. Illumina's Billion Cell Atlas: members only.
Genome-wide primary CD4 Perturb-seq (EMBL, bioRxiv 2026.03.09) has no accession.

**20.22 S12 OFFICIAL: 0.22676, rank 38 (entry 6NHNPxARhFMzvsG6RC8a, b24-61dd00) -- +0.0008 over S11.** The job sat in
`launching` for 3 h 18 min (16:42 -> ~20:00Z) while later entries from other teams were scored, then scored normally.
Scaled members vs S11: PDS 0.71346 and MSE 0.24739 IDENTICAL (sum side unchanged, as designed); NMAE 0.13454 vs 0.13871
(-0.0042), fid 0.03877 vs 0.03644 (+0.0023), reach 0.22871 vs 0.22343 (+0.0053), jac -0.00234 vs -0.00366 (+0.0013); raw
reach 0.2818, NMAE 0.9179. Harness predicted +0.0024 (6/6 panels; H1 +0.0014): about a third transferred -- the DE side
moved in the predicted direction on three of four members but NMAE went the other way. Consistent with S11's lesson
(harness DE metrics over-predict official DE gains); the sum-side mechanisms (S10k, S11 ensemble) remain the ones that
transfer at full size.
CAVEAT for 20.20: S11 (ci5 / lsb) and S12 show DE-side-only harness gains transfer at ~0-1/3 officially. The same-line
DE-side share (+0.052 on H1) may therefore shrink officially; the sum-side share (+0.022; sum-side mechanisms have
transferred at full size) is the reliable floor. Even at 1/3 transfer of the DE share, same-line data is worth ~+0.04
composite per line at 44 % coverage -- still the dominant final-phase lever.

**20.23 RELATED-LINE data (real scorer, S12-structure bases, seed 0) -> the final pipeline's same-line weight policy.**
| destination panel <- source (relation) | covered | w 0.5 | w 1.0 |
|---|---|---|---|
| H1 <- H1GW (same line, other lab) | 131/300 | +0.077 | +0.074 |
| KOLF <- H1GW (another pluripotent line) | 63/279 | +0.021 (PDS +0.008, MSE -0.051, reach +0.019) | +0.024 (MSE -0.078) |
| OR_HCT116 <- DLD1 (colorectal, different donor) | 58/263 | +0.005 (PDS +0.006, MSE +0.011, NMAE -0.014) | **-0.035** (MSE +0.149, reach -0.040) |
Per unit of coverage: same line ~0.17, pluripotent-to-pluripotent ~0.09, somatic related line ~0.02 at w 0.5 and strongly
negative at w 1 (a wrong direction on strong responders costs MSE once our own direction is fully replaced).
Policy now in `tv2_final.py`: exact identity -> W; stem fallback (H1GW for a destination with stem score >= 1.5) -> W (was
W/2); cross-label somatic maps (HEK293 -> HEK293T, HCT 15 -> DLD1; note HCT-15 and DLD-1 share a donor, so that map is closer
to same-line than this test) -> min(W/2, 0.5). Default W = 1.0 (the H1 plateau is 0.5-1).
Codex R17 revision of the policy (adopted): the H1 plateau shows no advantage of w 1 over 0.5 and the related-line curve
collapses at 1, so RECOMMENDED W = 0.5 and every non-exact source is capped at 0.5 (stem fallback min(W, 0.5), cross-label
min(W/2, 0.5)). Budget +0.02..+0.04 composite per comparable destination officially (single seeds; DE-side transfer
uncertain; the H1 pair may share guide designs, which would make it unusually favourable -- establish independence by
provenance, not by cosine). Before using any new screen: audit provenance, state (stimulation / differentiation / harvest)
and multiplicity; `ml_newsrc_build.py` now drops cells with unknown `nperts` (NaN > 1 was False) and warns when the
multiplicity is unknown altogether. Open control (Codex R17 item 3): replace the same covered KOLF targets with the existing
VCC2025-H1 source to separate "more pluripotent data" from H1GW-specific information.
Builder-exact agreement check (token `agc`: DE ranking = |DE delta| x agreement with the FAMILY-MEAN sign, frozen, as
tv2_agree_prep / the builder, instead of the harness's agreement with the candidate's sign), H1, seed 0: base 0.3029
(identical to the harness-native base), same-line H1GW w 1 **0.3764 (+0.0735** vs +0.0744 harness-native; reach +0.126,
NMAE -0.048, fid +0.022). The builder's same-line implementation keeps essentially the whole DE-side gain.

**20.24 Rehearsal 5 PASS (2026-09-29 ~21:40Z; `tv2_final_sim5.sh`, new 300-target panel, audit build).** D/E/F (A/B/C
controls relabelled) -> somatic (stem scores -1.33 / -5.11 / -2.43): no-KOLF S4, no-KOLF neighbour library, residual with
KOLF excluded, DE side from the override, pcr guard lambda 0.75, ew dropped, gamma_sum 0.970 / 0.972 / 0.948. G (H1
non-targeting cells on the official axis) -> STEM (score +2.255): S4 with KOLF, full library, identity NTERA-2 (rho 0.769,
below the exact-identity bar) -> stem fallback H1GW same-line at w 0.5 (capped) on 131 targets, residual with KOLF kept
(excl_G = []), gamma_sum 1.085. No errors, no partial prediction file (in-memory audit). The D Jurkat identity is refused
for the fusion boost by the KD gate (GSE249595 KD -0.30), as intended. Every final-phase branch has now run end to end.
Seeds @s1/@s2 + the KOLF control finish ~00:30Z; full table: reports/sc_hybrid_2026-09-25/loop6_sameline_final.txt

**20.25 Loop-6 confirmations (2026-09-30 ~00:10Z).** H1 same-line H1GW w 1, 3 seeds: +0.0744 / +0.0732 / +0.0739 (mean +0.0738);
sums only +0.0216 (3/3); DE-side share +0.0522 (3/3; reach +0.128, NMAE -0.043, fid +0.023). KOLF source-emphasis control
(Codex R17 item 3): replacing KOLF's covered targets with the ALREADY-FUSED VCC2025-H1 source at w 0.5 gives +0.0046 on 25
targets (0.00018 / target) vs H1GW +0.0205 on 63 (0.00033 / target) -> part of the stem-fallback gain is replacement /
pluripotent emphasis, about half is H1GW-specific information; both positive, so the stem fallback (capped 0.5) stands.
Full table: `loop6_sameline_final.txt`.

## 21. 12-hour loop B -> A (2026-09-30; user: "B부터 진행해줘 12시간동안 자동으로 돌리면서 판단하고 분석하고 테스트해서 알려줘")

**21.1 B CLOSED: destination-feature strength prediction does not pay.** `ml_strength_dest.py` (7 held-out panels,
scorer-space noise-corrected true strength): Spearman(ours, true) H1 .501 / KOLF .211 / OR_HCT116 .214 / RN_K562 .551 /
RN_HepG2 .740 / RN_Jurkat .577 / RN_RPE1 .644; a leave-one-panel-out GBM on destination features (p53-status-matched
source strengths, DepMap essentiality overall and in the destination, destination expression, family count / dispersion)
lifts it by ~+0.04 on average (.546/.216/.257/.606/.742/.616/.735). DepMap MDM2 / TP53 dependency gives clean functional
p53 status (functional: HepG2, HCT116, RPE1, A549; inactive: K562, Jurkat, CAL-33, DLD-1, C-33 A; HeLa by E6). Re-ranking
our amplitudes by the model (`ml_strength_realloc.py`, offline exact evaluator, S12-structure sums): composite H1 +0.0004,
KOLF 0.0000, OR -0.0050, RN_K562 -0.0023, RN_HepG2 +0.0058, RN_Jurkat -0.0038, RN_RPE1 -0.0006; even the ORACLE re-rank by
the true strength is negative on KOLF (-0.006), OR (-0.009), RPE1 (-0.001). Reason: the MSE-optimal amplitude is true
strength x our direction's accuracy (the least-squares slope), not true strength -- a strong target with a wrong direction
must NOT be amplified, and direction accuracy per target is the quantity already shown to be unpredictable (§19).
TP53 direction policy (`ml_p53_signature.py`): the functional-minus-inactive RP-knockdown axis (MDM2 / TP53I3 up; MKI67 /
TOP2A / BIRC5 down) carries 0.1 % (sums) / 0.3 % (DE) of S12's A/B/C energy, no target > 5 % -> nothing to remove on
A/B/C; kept for a p53-functional final line.

**21.2 A-1: program-space Wiener denoising of a source screen -- real reproducibility gain, sum-side use CLOSED.**
`ml_program_denoise.py`: K562 GWPS (9,729 targets) whitened (genes by a noise scale, rows by sqrt(n)), PCA, per-target
Wiener weights w = n tau^2 / (n tau^2 + 1). Cosine with the independent K562 essential screen (1,938 targets): raw 0.296
(terciles 0.156 / 0.293 / 0.439) -> K100 rho 0.25 0.370 (0.189 / 0.366 / 0.556; 86 % of targets improved). Codex R18: this is
same-lab reproducibility (shared guides / lab effects count as signal), the comparison mixes centred c50 with uncentred
log1p-CP10K spaces, the Wiener variance model is uncalibrated (scale includes biology; basis and power from the same rows)
-- a reproducibility result, not a transfer result. Harness gate (`ml_denoise_harness.py`, denoised K562 swapped into the
source library, S12-structure sums, offline): PDS falls on EVERY panel -- H1 -0.0068, KOLF -0.0036, OR_HCT116 -0.0136,
RN_HepG2 -0.0046, RN_Jurkat -0.0040 -- while MSE falls (H1 -0.004 ... OR -0.060); the MSE leg is not trustworthy because each
arm was budget-matched to its OWN recipe (Codex R18: the old `bmb` trap; S12 fixes the budget reference on the main
recipe). Mechanism: the Wiener weights keep the dominant shared programs and shrink target-specific low-variance
directions, so targets look more alike (hubness) -- a move ALONG the PDS-MSE frontier. Open: the DE side only (sums,
budget, overrides from the original sources; DE delta from the denoised library; harness token `dsrc<tag>`), real scorer.
Cross-lab / cross-chemistry check (`ml_denoise_crosslab.py`; Codex R18 item 1): denoised GWPS vs VIPerturb-seq K562 (another
lab, 10x Flex, 47 cells / target; 3,337 shared targets, same c=50 centred space, target gene masked, terciles by the SOURCE
norm): raw mean cosine 0.0369 (0.014 / 0.026 / 0.071) -> K100 rho 0.25 0.0490 (0.022 / 0.037 / 0.088), 65 % improved, i.e.
+33 % relative -- the gain is transferable signal, not only a same-lab artefact. Wiener weights vs plain subspace filtering
(wk = 1): identical (0.0490 vs 0.0489) -> the operator is "project onto the source's top-K response programs, keep 25 % of
the residual". It raises per-target direction agreement but makes targets more alike, which is why sum-side PDS falls; the
DE members (per-target direction) are the place it can pay.

**21.3 A-2: DE-SIDE-ONLY use of the denoised sources is POSITIVE on the real scorer.** Harness token `dsrc<tag>`: sums,
budget reference, overrides and PDS / MSE from the original sources; the DE delta (and its ranking score) from the SAME
recipe built on the denoised library. K562 only (dnK100r0.25; DE delta cosine to the original 0.97): H1 **+0.0057**
(3-seed base 0.3032 -> 0.3088; NMAE -0.0043, fid +0.0028, reach +0.0117, jac +0.0016; ~14 seed-sd), OR_HCT116 **+0.0062**
(NMAE -0.0079, fid +0.0042, reach +0.0081), RN_Jurkat +0.0005 (NMAE +0.0072 worse, fid +0.0028, reach +0.0008, jac +0.0017;
the DE delta changed less there, cosine 0.985). All-source version (dnALL: K562, HCT116, HEK293T, CD4 x3; the Orion lines
change little -- per-target cosine to raw 0.99 -- because they carry hundreds of cells per target; K562 / CD4 change ~0.83)
queued on six panels. Official path ready: builder key `de.src` {dir, knn} (DE-side snapshot from a denoised source
directory; `ml_denoise_apply.py` built `denoise_src/official_dnALL`, `ml_denoise_knn.py` the no-KOLF neighbour file);
S13 audit vs S12: main delta, overrides and sums BIT-IDENTICAL on A/B/C, kink / gamma identical, DE delta cosine 0.963.
Caveat: DE-side harness gains have transferred at ~1/3 officially (S12).
dnALL (K562 + HCT116 + HEK293T + CD4 x3 denoised) on the DE side, 6 panels: H1 +0.0045, OR_HCT116 -0.0002, RN_Jurkat
-0.0004, RN_HepG2 -0.0012, RN_RPE1 -0.0039, KOLF +0.0022 -> mean +0.0002, NMAE WORSE on 5/6 (+0.004..+0.013) -> CLOSED;
denoising CD4 (primary T, 2-3 % of values clipped at fp' = 0) and the Orion lines undoes K562's gain. K562-only so far:
H1 +0.0057, OR +0.0062, RN_Jurkat +0.0005 (NMAE worse on RN_Jurkat, as for dnALL on every RN panel). Official inputs for the
K562-only candidate built (`denoise_src/official_dnK562`, `knn_official_dnK562_noK.npz`); spec `spec_S13k.json` audited:
main delta / overrides / sums bit-identical to S12, DE delta cosine 0.952 (K562 is 63 % of the official fusion vs 45 % in
the harness, whose fusion carries three RN families). Pending: builder-exact (agc) paired check on six panels.

**21.4 Leaderboard 2026-09-30 (`reports/leaderboard_2026-09-30.json`) -- where the gap to the top 5 sits, and when it
opened.** Us 0.22676 (rank 38): PDS 0.822 / MSE 0.753 / NMAE 0.918 / fid 0.524 / reach 0.282 / jac 0.030. #5 [another team]
0.2737 = +0.047 composite: MSE +0.019, PDS +0.012, NMAE +0.008, reach +0.007. A cluster at 0.253-0.263 ([another team],
[another team], [another team], [another team] -- identical submission counts 8 / 16-17 / 26 / 27-28 / 30 / 32 at every snapshot, similar
scores: very likely ONE group on several accounts) sits at OUR MSE (0.72-0.75) with PDS 0.857-0.862: its +0.03 gap is
PDS +0.014, NMAE +0.010. Trajectories (best score per snapshot): the top teams were ALREADY at 0.22-0.33 on 09-18 (we were
0.136); since then they gained +0.03..+0.09 over ~20-25 submissions each; we gained +0.091 over 15 -- the gap is the
starting point, not the rate. Reading of "PDS up at equal MSE" (PDS and MSE use the same group sums, bulk_lognorm
log1p(5e4 S/T); PDS is per-target scale-free, ranks against all other real deltas, weights weak targets equally; MSE is
dominated by strong targets and carries the #348 rebate / budget machinery): better per-target DIRECTIONS (weak targets
included) with less-optimised MSE mechanics than ours -- i.e. an estimator gap, not an emission trick.

**21.5 dpj (DE-side projection onto the top-100 LIBRARY programs + 0.25 residual; energy kept 36 % by the pure
projection): H1 -0.0061** (real scorer, PDS / MSE unchanged by construction) -> closed. Codex R20 caveat: the harness
recomputes the candidate-sign agreement after the projection (`sch_emit.py:598-600`, it would overwrite an `agc` ranking),
so this rejects the combined intervention, not every program-space DE estimator.

**21.6 Weak responders' reproducible signal is TRANS, not cis** (`ml_weak_cis2.py`, H1; Codex R20 fixes: cells AND
controls split into independent thirds, strength bins from the third part, unannotated genes separate). Share of the
split-half consistent energy, weakest quartile: <=1 kb 7.6 %, 1-5 kb 0.7 % (both covered by the promoter prior),
5-20 kb 0.6 %, 20-100 kb 0.3 %, 0.1-1 Mb 2.0 %, >1 Mb same chromosome 4.6 % (= the uniform expectation), other
chromosomes 83.9 %, no TSS 0.4 %. Q2-Q4: cis <= 1.5 %. => cis beyond the 5 kb prior carries ~3 % of the weak targets'
reproducible signal; not a lever. (v1 with SHARED controls read 95 % trans -- common control noise inflated it.)

**21.7 Cfg.xw -- per-gene source weighting by the gene's expression in each source's controls** (tv2_model.py; e =
(cpm/(cpm+c))^xw; a c=50 fold change of a gene the source barely expresses is ~0 but entered at the full family weight,
K562 = 63 % of the official fusion). Fusion-only offline proxy (`ml_xw_diag.py`): per-target cosine with the truth
H1 +0.003 / OR +0.001, discrimination H1 +0.001/-0.003, OR -0.008/-0.013 (xw 0.5 / 1). Codex R20: the proxy used
row-specific masks and not the scorer space -> not a valid closure; full S12-structure offline check (`ml_xw_harness.py`,
budget fixed from the orig arm, fresh cache; frozen arm xw 0.5; kill unless mean J >= 0.03 with PDS > 0 on H1 and OR and
no RN regression) -> **CLOSED**: J (xw 0.5 / 1) H1 +0.011 / +0.011 (MSE -0.009 / -0.013, PDS ~0), OR_HCT116 -0.008 / -0.000
(PDS -0.0075 / -0.0072), RN_K562 +0.003 / +0.004, RN_HepG2 +0.003 / +0.005, RN_Jurkat -0.001 / -0.003, RN_RPE1 +0.003 / +0.005,
KOLF -0.012 / -0.008 (PDS -0.009) -> mean -0.0001 / +0.0019, 4/7 positive: an MSE-for-PDS frontier move again.

**21.8 Public-method survey 2026-09-30 (subagent; GitHub / HF / bioRxiv 08-01..09-30 / alphaXiv / GEO / Figshare / Zenodo).**
No public code, write-up or preprint reaches >= 0.25, PDS >= 0.85 or MSE <= 0.70. The only information at that level is
leaderboard `description` text: [another team] (#8, PDS 0.857 / MSE 0.629) earlier described "donor-attention response net, final
fit on all 300 H1-2025 targets + 8 evidence slots" -- a learned fusion supervised on the VCC2025 H1 truth (our analogue,
LOPO-learned per-target family weights, is closed on a valid measurement, §20.16); [another team] "DonorMix retrained with public H1";
[another team] (0.240) gene-program transfer weights learned on H1 Flex, mean-direction removal for the weakest 80 % of targets
(our §17.4 test of that: J < 0 on H1), per-context bulk amplitude, and a mean-preserving-spread emission (fid / jac only);
[another team] a learned emission-dispersion head. New data: GSE345058 (A549 Cas9-KO CROP-seq, 1,000 genes, 21 / 300 official
targets) -- low value; X-Atlas/Pisces still unreleased (HF repo empty since 2026-03-25); no 2026 Perturb-seq in Jurkat /
HeLa / CAL-33. Several identical-metric clusters of accounts exist beyond the 0.26 group.

**21.9 A-5: S4 WITHOUT THE FENG BIOLOGY for somatic destinations -- the first sizable sum-side candidate since S12.**
Mechanism: S4's refit takes each target's empirical mean over table_v2 biologies at EQUAL biology weight; S10k removed the
stem biology KOLF from it for somatic contexts (+0.0061 officially as part of S10k). FENG (HipSci iPSC; 10,276 rows;
covers 94 / 300 official targets in table_v2) and H1 (218 rows, 17 official) are the other stem biologies, never tested
inside S4 (fusion weight 0.25 only). Refits `s4i_refit_excl.py noKOLFFENG KOLF,FENG` (+ `_or` for HCT116; 20-35 s /
line, 12 GB). Offline exact S12-structure sums (`ml_s4x_harness.py`; min_cells ensemble, dk0 fusion, ct 0.5, guarded
pcr, budget + kink FIXED from the base key; only the S4 key changes; base = the official-equivalent rows+empirical
exclusion `noKOLFall` / `orRnoKOLF`):
| panel | dPDS | dMSE | J |
|---|---|---|---|
| OR_HCT116 (official targets) | **+0.0058** | -0.0263 | **+0.0405** |
| RN_Jurkat | -0.0014 | +0.0038 | -0.0071 |
| RN_HepG2 | -0.0012 | -0.0071 | +0.0048 |
| RN_RPE1 | +0.0006 | -0.0054 | +0.0070 |
| RN_K562 | -0.0017 | +0.0001 | -0.0039 |
Also removing H1 (noStem) is slightly worse than noKOLFFENG on every panel. The RN panels are essential-gene panels, where
FENG's genome-wide fitness arm measures the targets densely; OR carries the official target set. Official inputs built:
`tv2_s4_excl_final.py` (parametrised tv2_s4_nokolf_final.py; with KOLF only it reproduces the shipped S4noK to 6e-7) ->
`reports/submission_2026-09-30/S4noKF_ctx{A,B,C}.delta.npy` (vs S4noK: global cosine 0.818, 202 / 300 targets changed --
FENG's rows also move the refit's shared components and the outside-table empirical means); spec `spec_S14f.json` = S12
with only the override's part-0 S4 replaced. Builder audit: main delta and overrides BIT-IDENTICAL to S12, kink_base
identical (0.893 / 0.937 / 0.960), sums cosine 0.941 and DE side 0.940 vs S12 (S10k's no-KOLF sums were 0.934 vs S10i),
gamma_sum 0.978 / 0.989 / 0.983. Real-scorer check (sum AND DE side, agc, OR two seeds, RN x3) queued in
`aivc-agcbs4x-20260930`.
PRE-REGISTERED DECISION RULE (Codex R21, frozen 2026-09-30 09:31Z BEFORE any real-scorer S4noKF result was read):
upload S14f only if (1) both paired OR_HCT116 seeds are positive, (2) the RN laboratory-mean composite delta is >= 0 with no
RN panel below -0.003 (RN_K562 included), and (3) 0.5 x [0.75 x mean dC_OR + 0.25 x mean dC_RN] >= 0.005 (0.5 = explicit
transfer discount), i.e. OR ~ +0.013 composite if RN ~ 0. Codex's prior budget: ~+0.002 officially. Final phase: keep the
no-KOLF rule as the default; adopt no-FENG only if the RN gate passes too. Offline caveats: the offline budget reference is
the ensembled no-KOLF override (not S12's old main recipe; the real-scorer `brs` arms fix that), `sum.resid` is omitted,
OR's S4 is refit-only. Also measured offline: removing CD4 as well (noKFC) -> OR PDS +0.0049 but MSE +0.082 (OR offline MSE
is unreliable), RN +0.0008 / +0.0088 / +0.0123 / +0.0132 (Jurkat / HepG2 / RPE1 / K562) -- held for the real scorer.
Codex R21 check of "202 changed vs 94 covered": the 94 FENG-covered official targets carry the change (table-gene cosine
median 0.753, outside-table 0.652 vs S4noK); the other 206 move only at the 1e-3 level (0.9986 / 1.0000) through the
refit's shared components -- the count at a 0.999 threshold was inflated by those, not an anomaly.
Also closed offline (same harness): FENG ALSO zeroed in the fusion + library (`ml_s4x_harness_f0.py`, "+f0") on top of
S4noKF: J vs S4noKF alone OR -0.010, RN_HepG2 -0.0025, RN_K562 -0.003, RN_RPE1 +0.003, RN_Jurkat +0.0005 -> FENG's fusion
weight 0.25 stays. `tv2_final.py --s4-excl` (default KOLF = the bit-exact S4noK path; KOLF,FENG -> tv2_s4_excl_final.py)
added for the final phase, dormant.

**21.10 S13k builder-exact (agc) paired check, real scorer** (K562-program-denoised DE side; sums identical):
H1 0.3029 -> 0.3089 (**+0.0060**; NMAE -0.0047, fid +0.0026, reach +0.0141, jac +0.0016), OR_HCT116 0.4193 -> 0.4260
(**+0.0067**; NMAE -0.0094, fid +0.0042, reach +0.0084), RN_Jurkat 0.4831 -> 0.4838 (+0.0007; NMAE +0.0071 worse, fid +0.0031), RN_HepG2 +0.0057 (NMAE -0.0011, fid +0.0048, reach +0.0106), RN_RPE1
+0.0036 (NMAE +0.0044 worse, fid +0.0027, reach +0.0114) -> **5 / 5 positive, mean +0.0045**; KOLF 0.0586 -> 0.0614 (+0.0028) makes it **6 / 6, mean +0.0043**; at the ~1/3
DE-side transfer seen for S12 that is ~+0.0015 officially -- a ride-along, not a slot by itself.
Bundle spec `reports/submission_2026-09-30/spec_S14fk.json` (S14f + de.src K562) audited: main delta / overrides identical
to S12, sums IDENTICAL to S14f, DE side cosine 0.957 vs S14f (0.891 vs S12).
REAL SCORER, S14f vs base (agc, sum AND DE side), first read 11:05Z: OR_HCT116 s0 **+0.0092** (PDS +0.0073, MSE -0.0302 --
the offline OR estimates +0.0058 / -0.0263 were close, NMAE +0.0074, reach +0.0166); RN_HepG2 +0.0018 (NMAE -0.0029),
RN_RPE1 +0.0031 (NMAE -0.0045), RN_Jurkat -0.0008 (NMAE -0.0023, MSE +0.0038) -> RN mean +0.0014, none below -0.003; gate
criterion 3 so far 0.5 x (0.75 x 0.0092 + 0.25 x 0.0014) = +0.0036 < 0.005. Pending: OR seed 1, RN_K562, bundle S14fk.
Precedent (context, not part of the rule): S10k's Orion offline J / 6 = +0.0061 matched its official +0.0061.
OR seed 1 (11:38Z): base 0.4204 -> S14f 0.4321 (**+0.0117**; PDS / MSE identical to s0 -- the sums are deterministic --
NMAE +0.0055, reach +0.0268). Gate so far: OR both positive (mean +0.0105) YES, RN (3 panels) mean +0.0014, min -0.0008 YES,
discounted estimate 0.5 x (0.75 x 0.0105 + 0.25 x 0.0014) = **+0.0041 < 0.005 -> S14f ALONE FAILS criterion 3** (RN_K562
still pending). The bundle S14fk (S14f + K562-denoised DE side) is scored next under the same gate.
RN_K562 (12:18Z): 0.4204 -> 0.4207 (+0.0003; NMAE +0.0060, reach +0.0126). FINAL S14f gate: OR +0.0092 / +0.0117 (both > 0),
RN mean +0.0011 over 4 panels (min -0.0008), discounted 0.5 x (0.75 x 0.0105 + 0.25 x 0.0011) = **+0.0041 -> FAILS the
+0.005 bar narrowly**; S14f alone is a ride-along, not a slot.
CD4 ALSO out of S4 (noKFC), real scorer OR_HCT116: 0.4193 -> 0.4101 (**-0.0092**; MSE 0.071 -> 0.150) -> closed (the
offline OR MSE +0.082 was right); CD4 carries the official targets' core signal.

**21.11 BUNDLE S14fk = S14f + K562-denoised DE side (spec_S14fk.json), real scorer, agc, same R21 gate.**
OR_HCT116 s0 0.4193 -> 0.4337 (**+0.0144**; PDS +0.0073, MSE -0.0302, NMAE -0.0005, fid +0.0054, reach +0.0182) -- the two
parts add (S14f +0.0092, S13k +0.0067); RN_HepG2 **+0.0092** (NMAE -0.0061, fid +0.0049, reach +0.0169), RN_RPE1 **+0.0076**
(reach +0.0161), RN_Jurkat -0.0013 (NMAE +0.0036, reach -0.0039), RN_K562 = S14f's +0.0003 (K562 excluded there, the denoised
DE side is a no-op). RN mean +0.0039, min -0.0013; discounted 0.5 x (0.75 x 0.0144 + 0.25 x 0.0039) = **+0.0059 >= 0.005**;
criterion 1 (both OR seeds > 0) waits for OR s1 (13:35Z).
Codex R22 (final review): the bundle test is legitimate but carries selection optimism (built after S14f missed). Official
FORECAST with the incremental DE gain discounted by 1/3 (S12's DE transfer): E = 0.5 W(dS14f) + 1/3 W(dS14fk - dS14f),
W(x) = 0.75 x_OR + 0.25 mean x_RN -> needs OR bundle s1 >= ~+0.0122 to reach +0.005 (the plain R21 gate needs >= +0.0096).
Budget ~+0.005..+0.006, judgemental range -0.003..+0.013. Propose ONE S14fk upload only if s1 >= ~+0.0122; do not change the
final-phase default on this validation upload alone. What an upload would teach: whether Orion-led source-exclusion gains
transfer to A/B/C at all (official PDS / MSE isolate the S4 change, since S14fk's sums equal S14f's).
OR s1 (13:15Z): base 0.4204 -> S14fk 0.4382 (**+0.0178**). FINAL S14fk gate: OR +0.0144 / +0.0178 (both > 0), RN mean +0.0039
(min -0.0013), R21 discounted 0.5 x (0.75 x 0.0161 + 0.25 x 0.0039) = **+0.0065 -> PASS**; R22 forecast
0.5 x W(dS14f) + 1/3 x W(increment) = 0.0041 + 0.0017 = **+0.0057 >= 0.005 -> PASS** (s1 >= +0.0122 met).
=> S14fk is the first candidate since S12 that clears the no-marginal-submission bar; budget ~+0.005..+0.006 official
(Codex range -0.003..+0.013). **BUILT 14:28Z**: `reports/submission_2026-09-23/S14fk.vcc`, sha256
5cbcbf8230a58331a8fb7e3fe041a75e4010deeccf3c66f428a549f0cd9d1f78 (3.12 GB; opaque name would be b24-5cbcbf); vcc prep passed
(360,000 cells, 18,533 genes, A/B/C 120k each, targets verified, counts preserved); build dumps (delta / override / sums /
DE side, A/B/C) BIT-IDENTICAL to the S14fk audit; pred.h5ad deleted. NOT uploaded: needs the user's go-ahead.
S13k alone on OR s1: 0.4204 -> 0.4288 (+0.0084; s0 +0.0067). S14fk build started 13:32Z (unit aivc-build-s14fk-20260930).
Next-candidate probe, closed: FENG out of the DepMap-neighbour LIBRARY only (its part 0 is the equal-weight recipe, FENG 1.0)
on top of S4noKF (`ml_s4x_harness_lib.py`, "+lf0"): J OR -0.0023, RN_Jurkat +0.0015, RN_HepG2 -0.0010, RN_RPE1 +0.0036,
RN_K562 -0.0022 -> mean -0.0001 (the library moves the sums by cos 0.999 only). The S4 exclusion is where FENG mattered.
RN second seed for the bundle (14:40Z; base noKOLFall agc@s1 vs S14fk@s1): RN_Jurkat -0.0014 (s0 -0.0013), RN_HepG2 +0.0057
(s0 +0.0092), RN_RPE1 +0.0079 (s0 +0.0076) -> two-seed means -0.0014 / +0.0074 / +0.0077, RN mean with RN_K562 +0.0035,
R21 discounted with two seeds everywhere **+0.0065** (unchanged -> PASS). Caveat to carry: the Jurkat-lineage panel is
reproducibly slightly negative (MSE +0.004, PDS -0.0014 from the S4 change; the DE denoise is ~0 there) -- context A is the
official context this candidate is least likely to help.
Final-phase path verified, rehearsal 6 (`scripts/tv2_final_sim6.sh`, rehearsal 5's new panel D/E/F somatic + G stem, reuses its
artefacts): `tv2_final.py ... --s4-excl KOLF,FENG` ran tv2_s4_excl_final.py on the new panel (refit_v2 reproduced max_abs 0.0,
rows 200,163 -> 177,743), routed S4noFK_FENG-KOLF to D/E/F and kept S4_ctxG for the stem context G; builder audit vs
rehearsal 5: D/E/F main delta + overrides identical, sums cos 0.954 / 0.954 / 0.955, DE side 0.952; G bit-identical.

**21.12 Robustness before the upload decision: OR_HCT116X (HCT116 AND HEK293T held out -- no close epithelial source left,
nearer to HeLa / CAL-33's situation).** S4 refit with both held out (orX / orXnoKOLF / orXnoKOLFFENG); real scorer, agc:
base 0.3530 -> S14f 0.3592 (**+0.0062**: PDS +0.0101, MSE -0.0047, NMAE +0.0110, reach +0.0235) -> S14fk 0.3634 (**+0.0104**:
PDS +0.0101, MSE -0.0047, NMAE -0.0015, fid +0.0031, reach +0.0193). The gain holds (PDS rises MORE than on OR) and the
K562-denoised DE side cancels S14f's NMAE loss -- the two parts are complementary.

**21.13 OFFICIAL: S14fk = 0.23273 (rank 34, entry ckxrFX2RigxKeIvjbuUL, model b24-5cbcbf; uploaded 2026-09-30 18:52Z with the
user's go-ahead) -- NEW BEST, +0.00597 over S12 (0.22676).** Forecasts: R22 +0.0057, R21 discounted +0.0065 -> the TOTAL was
predicted well, the COMPOSITION was not. Member contributions (composite): NMAE 0.9179 -> 0.9047 **+0.0037** (the bulk),
PDS 0.8217 -> 0.8238 +0.0008, MSE 0.7529 -> 0.7493 +0.0006, reach 0.2818 -> 0.2850 +0.0006, fid +0.0003, jac -0.0001.
Reading: (1) the sum-side S4 change transferred at ~1/4 of the OR harness (PDS +0.0073 -> +0.0021, MSE -0.030 -> -0.0036);
(2) the DE side gained far MORE officially (NMAE -0.013) than any harness panel showed (OR -0.0005, OR_HCT116X -0.0015,
RN -0.006..+0.004) -- after S11 / S12, where DE-side harness gains transferred at 0..1/3, this is the first DE-side change
that over-delivered. The bundle cannot say which part moved NMAE (the DE side carries both the new S4 and the K562 denoise).

## 22. After S14fk: what would it take to reach the top 10? (2026-09-30 19:30Z-, user: "상위 10위안에 들어가기위해서 어떤걸 다시만들거나 개선해야해? 이런식으로 조금씩 개선하는거는 가망이없어보여"; chosen plan "②로 시작, ① 병행")

**22.1 Gap to ranks 6-10 (leaderboard 2026-09-30 19:30Z, `leaderboard_2026-09-30b.json`):** us 0.2327 (rank 34) vs mean
0.2692 = +0.0365: MSE +0.0175 (48 %), PDS +0.0106, NMAE +0.0062, jac +0.0028, reach +0.0021, fid -0.0027 (we lead).
Their raw means PDS 0.852 / MSE 0.649 / NMAE 0.882 / reach 0.296 / jac 0.035. The rank-10 line moves ~+0.003/day
(0.2286 on 09-18 -> 0.2646 on 09-30). Reading: at optimal amplitude MSE ~ 1 - rho^2, so ranks 6-10 have ~1.18x our
energy-weighted cosine -- exactly the 1.10-1.25x that §20's ceiling puts the BEST PER-TARGET MIX of our own sources above
us. The top-10 gap is consistent with better per-target source mixing, not with more data.

**22.2 ② CLOSED: removing each source's SOURCE-SPECIFIC programs before transfer** (`ml_srcspec.py`, `ml_srcspec_harness.py`;
K562 only; shared subspace = top-50 PCs of HEK293T + CD4 x3 (HCT116 left out), specific programs = top-M PCs of K562's
residual). Diagnostic confirms the premise: the top-5 K562-specific programs ARE erythroid / myeloid (HBZ, HBG1, HBD, ALAS2,
BLVRB, SLC25A37, GYPA, HEMGN, PRG2, MS4A3, GP1BB) and carry 8.2 % of K562's (sqrt-cells-weighted, centred) response energy,
as much as all 50 shared programs (8.2 %). But removing them does not pay (offline exact sums, budget fixed from orig):
M = 5 / 20: H1 +0.0016 / -0.0015, OR -0.0001 / -0.0120, KOLF -0.0003 / -0.0038, RN_HepG2 -0.040 / -0.039, RN_Jurkat
-0.040 / -0.064, RN_RPE1 -0.015 / -0.017. The erythroid genes are gated out in non-erythroid destinations anyway, and the
same programs carry the integrated-stress genes (DDIT3, ATF5, BTG1/2, IER3, TXNIP) that essential-gene knockdowns move in
every line -- removing them costs MSE. => the mixing gap is per-target source choice, not K562 idiosyncrasy.
**22.3 ① design written**: `/home/wonmin-kim/AIVC/docs/superpowers/specs/2026-09-30-learned-source-mixer-design.md` (learned
per-(target, gene, source) mixing weights trained on leave-one-family-out pseudo-destinations; the 7 real panels only for
selection; full-recipe gate >= 30 % of the oracle headroom on OR / OR_X with RN >= 0). Awaiting the user's review.
**22.4 Codex R23 on the ① design:** the MSE->cosine argument only shows numerical compatibility (with the budget-matched
capped-jackknife MSE, MSE ~ 1 - K rho gives 1.40x, not 1.18x; the ceiling reports per-target cosines by stratum), so run a
<= 1-day PILOT before the 3-5-day build: a low-dimensional gate w = w_base * exp(theta^T phi) (~10 standardised features,
L2 toward theta = 0, direct response-prediction objective with whole-target batches), OUTER exclusion of the evaluated
biology from the whole training graph (H1 + H1GW together; HCT116 + HEK293T for OR_X; CD4 states grouped), theta = 0 must
reproduce the incumbent, injection before the neighbour blend inside the exact S14fk sum construction with its original
budget, DE branch fixed; kill if < 30 % of the same-construction oracle on OR or OR_X, RN mean J < 0, or PDS not positive on
both Orion variants. Spec updated (§6b).
**22.5 ① PILOT RESULT: KILL** (`scripts/ml_gate_pilot.py`, spec §6b; user "하루 파일럿 진행"). Gate w = w_base exp(theta.phi) on
the fusion's family weights, 9 features (source expression, source-destination expression gap, |z|, sign agreement with
the other sources, deviation from them, cells, knockdown, response norm, control-profile similarity), fitted by L-BFGS
(expression-weighted MSE + per-target cosine + L2, lambda by 3-fold held-out pseudo-destination families) on
leave-one-family-out PSEUDO-destinations with the evaluated panel's biology removed from the whole training graph; injected
into the exact S14f sum construction (theta = 0 parity with delta_s4 5e-8..1.4e-7 on all 6 panels), budget fixed from theta = 0.
Learned theta (OR): |z| -1.68, agree +0.51, cells +0.67, response norm +0.46, control similarity +0.66, knockdown -0.20,
expression terms ~0. Results J (dPDS / dMSE; recovered fraction of the per-target NNLS oracle through the same construction):
OR_HCT116 +0.029 (-0.006 / -0.041; 6 %), OR_HCT116X -0.045 (-0.016 / +0.010), RN_HepG2 -0.044, RN_Jurkat -0.046, RN_RPE1 -0.005,
RN_K562 -0.001; PDS falls on 6 / 6. Oracle J +0.07..+0.47. => Source-to-source weighting rules learned on pseudo-destinations
do not hold for real destinations (Codex R23's main risk). Fifth failed attempt at per-target / per-gene source weighting
(MODZ, LPO-GBM, control cosine, §19/§20.16 learned per-target weights, this gate): the oracle headroom is real and NOT learnable
from test-time features with these data. Every test-time reweighting measured today (program denoising, xw, this gate) moves
ALONG the PDS <-> MSE frontier. The 3-5-day full build is cancelled per the pre-registered kill rule.

**22.6 BUG in the official S13 / S13k / S14fk DE-side inputs (found 2026-09-30 ~21:40Z while porting S14fk to the final
pipeline).** `ml_denoise_apply.py` was called with a RELATIVE source dir, so every symlink it wrote into
`denoise_src/official_dnK562/` (26) and `official_dnALL/` (16) is relative to the wrong directory and BROKEN; the builder's
`de.src` therefore saw ONLY the real file K562.npz (every other family `_src_path(...).exists()` False -> skipped). The
UPLOADED S14fk DE side = S4noKF + a K562-ONLY (denoised) fusion + the knn from the full denoised harness library (built by
ml_denoise_knn.py with absolute paths -- unaffected). Re-audit with fixed links (`official_dnK562_fixed`): intended vs
shipped DE side cosine 0.953 (per target median 0.955, min 0.924); sums identical. The harness dsrc arms (dnK100r0.25 /
dnALL, absolute links) are NOT affected -- only the two official dirs had broken links (scan of every denoise_src subdir).
=> S14fk's official +0.0060 (NMAE -0.013) belongs to the SHIPPED variant, which was never harness-tested; the harness tested
the INTENDED one (6/6 positive). Fix: `ml_denoise_apply.py` now resolves paths and asserts no broken links. The shipped
variant is being scored on the harness (token `dsl<libtag>`: DE-side fusion from `dnK100r0.25K562only`, library from
`dnK100r0.25`; queue `aivc-shipde-20260930`) to decide which DE side the final pipeline should carry.
**22.7 Final-phase port of the S14fk structure DONE (dormant options, default paths bit-exact).** Builder: `de.src.knn` may be
a {context: path} dict and a de.src dir with a broken symlink now aborts the build. `tv2_final_knn.py`: env
`TV2_KNN_DENOISE=K562:100:0.25` denoises that library copy before the neighbour estimates. `tv2_final.py --de-denoise
{off,intended,shipped}` builds `de_dnK562` (ml_denoise_apply on the panel's sources) or the K562-only `de_dnK562only`
(shipped), the per-context denoised knn (`knn_noK_dnK562` somatic / `knn_dnK562` stem, library files reused by absolute
links) and sets spec de.src. Rehearsal 7 (`scripts/tv2_final_sim7.sh`, new panel D/E/F somatic + G stem, `--s4-excl
KOLF,FENG --de-denoise shipped|intended`, audit): sums bit-identical to rehearsal 6 in every context; DE side cosine vs
rehearsal 6 shipped 0.947 (D/E/F) / 0.958 (G), intended 0.977 / 0.980 -- the same deviation sizes as the official S14fk audit
(0.957 / 0.969); knn routed per context. Which variant the final command carries waits for the harness comparison (§22.6).
**22.8 Shipped vs intended S14fk DE side on the harness (real scorer, agc; composite delta vs the official-equivalent base):**
| panel | S14f (S4 only) | S14fk INTENDED (all families, K562 denoised) | S14fk SHIPPED (K562-only DE fusion) |
|---|---|---|---|
| OR_HCT116 | +0.0092 | **+0.0144** | +0.0041 |
| OR_HCT116X | +0.0062 | **+0.0104** | +0.0090 |
| RN_HepG2 | +0.0018 | **+0.0092** | -0.0019 |
| RN_RPE1 | +0.0031 | **+0.0076** | +0.0025 |
| RN_Jurkat | -0.0008 | -0.0013 | -0.0004 |
| H1 (stem: full S4, DE side only) | -- | **+0.0061** | -0.0009 |
Intended beats shipped on 5 / 6 (mean +0.0057); the shipped DE side is even below S14f alone on OR and RN_HepG2. The official
S14fk +0.0060 was therefore reached WITH the weaker DE side. Final pipeline carries `--de-denoise intended`. Corrected package
spec `reports/submission_2026-09-30/spec_S14fk2.json` (de.src -> official_dnK562_fixed; audit: DE side cos 0.953 vs S14fk, sums
identical); expected official gain over S14fk with the R21 discount ~+0.003 (0.5 x (0.75 x 0.0059 + 0.25 x 0.0038)), with the
1/3 DE-side prior ~+0.002 -- below the +0.005 bar by itself, but the only official check of which DE side the final phase
should carry.
User decision 2026-09-30 ~23:00Z on S14fk2: "올리지 않음" -- keep the +0.005 rule; the final phase carries the harness-preferred
`--de-denoise intended` on the 5/6 evidence. The S14fk2 package is built and kept (not uploaded).
S14fk2 BUILT 23:52Z: `reports/submission_2026-09-23/S14fk2.vcc`, sha256 54c07ad2846f7d5bcbeb632741beaba13fc95038221b8e100bd9cc97d7182259
(3.12 GB), vcc prep passed, build dumps bit-identical to the fixed-link audit; pred.h5ad deleted; NOT uploaded (user).
Default-path regression after today's tv2_final.py edits: the rehearsal 5 command without the new flags reproduces
sim5_spec.json with 0 differences.

**22.9 Harness -> official DE calibration (2026-10-01, `scripts/ml_de_calibration.py`; user "1+2 함께").** Per-member raw deltas,
official vs each harness panel, for every official transition with matching arms:
- S10b->S10c (cpc, a LARGE change): official NMAE -0.101; harness H1 -0.072, KOLF -0.068, RN_HepG2 -0.097 -> all right.
- S11->S12 (DE side follows the sums): official NMAE **+0.0022 (worse)**, fid +0.0007, reach +0.0049; harness NMAE H1 -0.0034,
  OR -0.0044, RN_Jurkat -0.0011, RN_K562 -0.0023, RN_RPE1 -0.0001 (wrong sign), RN_HepG2 +0.0021 (right).
- S12->S14fk (shipped): official NMAE **-0.0131 (better)**; harness OR +0.0218, RN_HepG2 +0.0110, RN_RPE1 +0.0118, RN_Jurkat
  +0.0109, H1 +0.0053 (all wrong sign) -- only **OR_HCT116X -0.0096** (right sign, ratio 1.37).
- fid: harness always over-predicts (official / harness 0.04..1.2, typically 0.1-0.4); reach inconsistent.
=> For SMALL DE-side changes the harness NMAE has predicted the wrong official sign in 2 of 3 transitions -- the root of the
erratic DE transfer (S11 nil, S12 1/3, S14fk over-delivered). Hypothesis: the RN panels share lab / protocol with K562 GWPS
(Replogle / Nadig) and OR_HCT116 keeps its same-lab HEK293T source, so they over-reward K562-like DE signal; OR_HCT116X
(HCT116 + HEK293T out: no same-lab, no close-lineage source) is the only panel built like the official A/B/C situation and
the only one that got S14fk's NMAE right. Test queued: S11->S12 on OR_HCT116X (unit aivc-orxcal-20261001); S14fk2's official
result (uploaded 2026-10-01 07:18Z, entry gJHB2FvYxsWSnoTCp9zk) is the next calibration point (harness intended - shipped:
composite OR +0.0103, OR_X +0.0014, RN_HepG2 +0.0111, RN_RPE1 +0.0051, RN_Jurkat -0.0009, H1 +0.0070).

**22.10 OFFICIAL S14fk2 = 0.23543 (entry gJHB2FvYxsWSnoTCp9zk, b24-54c07a; user "1+2 함께") -- NEW BEST, +0.0027 over S14fk.**
PDS / MSE identical (sums unchanged); NMAE 0.9047 -> 0.9033 (-0.0014), fid +0.0014, reach +0.0057, jac +0.0010. => The INTENDED
DE side (all families, K562 program-denoised) beats the accidentally shipped K562-only one OFFICIALLY too: the final phase
carries `--de-denoise intended` (now confirmed on the official scorer). Forecast check: R21-style composite 0.5 x (0.75 OR +
0.25 RN) said ~+0.003 -> official +0.0027 (and S14fk: +0.0065 -> +0.0060).
Calibration closed with three small DE-side transitions: NMAE sign right -- OR_HCT116X 2/3 (S11->S12 harness +0.0179 vs official
+0.0022; S12->S14fk right; S14fk->S14fk2 wrong, +0.0081 vs -0.0014), majority panels 1/3; reach / fid sign mostly right but
magnitudes ~1/3. => No single panel predicts official DE members; the discounted multi-panel COMPOSITE does (+-0.002-0.003).
Rule from now on: decide on the R21 composite (OR + OR_X + RN, x 0.5), never on a single panel's NMAE.
Leaderboard 2026-10-01 08:00Z: us rank 35 (0.23543); rank 10 = 0.2729, rank 5 = 0.2937 (the line keeps moving ~+0.004..0.008/day).

## 23. Brainstorm spike 2026-10-01 (user: "크게 개선할 수 있는 방향을 탐색해서 테스트해줘"; chosen "P1 Jurkat 재분석 + P2 병행")

**Leaderboard 10-01** (best entry per team): #1 [top team] 0.4289 (PDS .897, MSE .349, NMAE .796, reach .456); #2-5 0.29-0.33;
#10 = 0.2744; us #35 = 0.2354. Gap to #10 by member credit: MSE +0.019 (50 %), PDS +0.012, reach +0.006, NMAE +0.004,
fid -0.002 (ours is better), jac 0. Profile of #5 [another team]: PDS 0.841 (= ours) but MSE 0.592 -> the top teams win on the
STRONG targets' direction, not on discrimination. No team in the top 20 has a method description on the API any more.
Fact check (subagent over FINDINGS §1-22): ensemble axes beyond min_cells dilute (c, cf, fs, bw, LOFO, Dirichlet); untested
axes = PC-removal k and cell bootstrap; `lcp` (p-value order follows |d|, reach) was scored on RN only (+0.0018 mean, mostly
NMAE), never on H1 / OR; TF-regulon priors never tried (75/300 official targets are TFs by Lambert 2018, 45 of them ZNFs).

**23.1 P1 KILLED -- GSE249595 (Jurkat, high-MOI genome-scale CRISPRi) has no DETECTABLE per-target signal at its depth, with
either estimator (Codex R24: say "not detectable / not usable", not "proven noise").**
Pre-registration `jurkat_lr/PREREG.md`; code `scripts/ml_jurkat_lr.py` (ridge regression over the called guides of every cell,
covariates channel / log UMI / n_guides, 16 channels, 827k cells; guide call = UMI >= 3 and >= 0.25 x the cell's max guide UMI
instead of the plain >= 3 of §20.9; STIM split into random halves A / B, UNSTIM separate; the marginal estimator from the
SAME calls as the control). Run time 10 min. Results (`jurkat_lr/full/validation.json`), official targets, 299 covered,
412 cells / target (209 per half):
  split-half cosine LR_A vs LR_B: mean +0.003, median +0.002 (strong-20 % by K562 norm: +0.005 / +0.003)  -> GATE (>= 0.10) FAILS
  split-half marginal MG_A vs MG_B: +0.003 / +0.002;   LR vs MG cosine 0.98 (the regression changes nothing)
  STIM vs UNSTIM (281 targets, 31 cells in UNSTIM): +0.001;   LR vs K562 / HCT116 / CD4: +0.002..+0.007
  RN_Jurkat targets (ESSENTIAL genes, 294 covered, 387 cells): split-half LR +0.008 / +0.007 -- noise even where responses must be strong
  on-target KD (fp/fc - 1): LR -0.42, MG -0.44, old JURKAT_STIM (>= 3 UMI) -0.31 vs K562 -0.87, CD4 -0.86, HCT116 -0.70
Reading: stricter guide calls do sharpen the on-target knockdown (-0.31 -> -0.44, the channel-1 probe reached -0.74 at
0.5 x max) but the downstream response is irreproducible at ~200 cells per half, for regulators AND for essential genes;
co-perturbation contamination is not the cause (LR = MG). Likely causes: shallow cells (median 2,900 UMI, 523 genes),
weak knockdown (~40-60 %), and/or an early harvest. Alignment audited after Codex R24 (`jurkat_align_audit.py`): the
transcriptome / guide / hashing barcode lists are identical in all 16 channels and the feature lists identical across channels,
so the estimate is not an indexing artefact. Codex's power illustration: 500 affected genes at LFC sd 0.3 give an expected
split-half cosine 0.84 at 100 CPM, 0.34 at 10 CPM, 0.05 at 1 CPM for 580k UMI per half -- so a weak, low-expression response
could sit below our gate; but a source that cannot show split-half reproducibility above the cross-line agreement (K562 vs
HCT116 +0.03) cannot add information to the prediction either way. §20.9's closure stands for the USE of this dataset.
Context A has no usable same-line data; do not retry this dataset. The probe sources are kept in `jurkat_lr/sources/` (JURKAT_LR / JURKAT_MG,
both halves and the smoke run deleted; 16.6 GB of normal-equation stats deleted); nothing was added to the production
source directory.

**23.2 P2 -- `lcp` (DE-side re-allocation weakened on high-expression genes) is an NMAE lever, first real-scorer read on a
somatic official-target panel.** Token `lcp<eta>` (sch_emit.py): inside the dual emission's top-K within-column re-allocation,
lam <- lam x clip((median CPM / CPM)^eta, 0.2, 1), i.e. high-expression genes keep more of the plain control-shaped column;
sums stay pinned (PDS / MSE bit-identical). It had only been scored on the RN panels at the OLD structure (§5: mean +0.0018,
all four NMAE -0.006..-0.012 raw, fid mixed, RPE1 -0.0026); the H1 / KOLF job (026) was skipped. Today, S14fk structure,
real scorer, OR_HCT116X: base 0.3634 -> lcp0.25 **0.3742 (+0.0108)**: NMAE raw 0.830 -> 0.789 (-0.042; scaled +0.0685),
reach -0.005 raw, fid +0.0005, jac 0, PDS / MSE unchanged to 1e-14. Queue (aivc-sch-queue-20261001, serial): H1 lcp0.25,
OR_X base@s1 + lcp0.25@s1 (seed), OR_X / H1 lcp0.5 (dose), RN_Jurkat / RN_HepG2 / RN_RPE1 lcp0.25 at the current structure
(gate), KOLF lcp0.25. Decision rule (pre-registered here): H1 is the primary DE-side panel; adopt only if H1 > 0 and the
R21-style composite 0.5 x (0.75 x mean(OR_X, H1) + 0.25 x RN mean) >= 0.005 after the 1/3 DE-side discount of R22 -- i.e.
the harness mean has to be >= ~+0.015, OR a smaller gain rides along in the next substantial build (correctness-fix rule).
Builder port if adopted: `tv2_emit.emit_dual_iter` line ~198 (`lam = lam_rank[:len(top)]`) gets the same expression factor
from `full_mean`; spec key `dual.lcp`; `tv2_final.py --dual ...,lcp0.25`; default off, bit-exact old path.
H1 (primary DE-side panel), real scorer: base 0.3089 -> lcp0.25 **0.3104 (+0.0015)**: NMAE raw 0.8576 -> 0.8568 (-0.0008 only),
reach +0.0074 raw (+0.0084 scaled), fid 0, jac -0.0002; PDS / MSE unchanged. So the OR_X gain is an NMAE effect of that panel's
depth / expression profile; on the official-chemistry proxy the lever is small and comes through reach. Pre-registered rule
with the numbers so far: 0.5 x (0.75 x mean(+0.0108, +0.0015) + 0.25 x RN +0.0018 (old structure)) = +0.0025 < 0.005 ->
NOT a slot by itself. (Codex R24: use ONE discount -- this R21 composite is the gate; the separate "1/3 DE-side transfer"
reading of the undiscounted harness mean +0.0051 gives the same budget, +0.0017 official, with a possible reversal. Correction:
the official NMAE 0.903 is WORSE than both H1 0.858 and OR_X 0.830, not between them, so neither panel pins the official
emission regime.) Ride-along status is NOT automatic: PDS / MSE invariance does not make it a correctness fix (reach worsens
on OR_X); it stays OFF by default until the current-structure robustness queue and the pin-error diagnostic (23.4) are read. Remaining queue (seed, dose, RN at the current structure, KOLF) refines robustness only.

**23.3 `lcp` ported as a dormant option (user choice 2026-10-01: "빌더·최종 파이프라인에 휴면 옵션으로 포팅").**
`tv2_emit.emit_dual_iter(..., lcp=0.0)` multiplies the top-K re-allocation strength by clip((median CPM / CPM)^lcp, 0.2, 1)
exactly as sch_emit's token (after the agree ordering, before the up / down split); builder spec key `dual.lcp`
(`tv2_submit_build.py` passes `float(dual.get("lcp", 0.0))`); `tv2_final.py --dual ...,lcp<eta>`. Parity
(`/home/wonmin-kim/.claude/jobs/43119e26/tmp/lcp_parity.py`, context-A official controls, 6 S14fk2 dsum rows, cpc + ci5 + lsb,
400 cells): patched module with lcp 0 vs the pre-patch module |diff| = 0.0 (bit-exact old path); lcp 0.25 changes 9,198
columns with every column sum equal (PDS / MSE invariant); the factor on a target's top-K is median 1.0 (genes below the
median CPM untouched), min 0.33, none at the 0.2 floor. Backups: `.claude/jobs/43119e26/tmp/{tv2_emit,tv2_submit_build,
tv2_final}.py.bak_lcp`. Rehearsal 8 (`scripts/tv2_final_sim8.sh` = sim7 intended + lcp0.25, audit mode) running as
aivc-sim8-20261001. NOT in any package yet; the final-command decision waits for the robustness queue + a Codex round.

**23.4 Codex round 24 (`debate_loop/codex_round24.md`) and rehearsal 8.** Q1: kill the production candidate but do not call
the dataset proven noise; audit barcode / feature alignment (DONE, identical everywhere); the last cheap test would be
injecting a 500-gene response into the real design and recovering it (power test) -- not run: even a sub-gate real signal
cannot help context A. Q2: the pin-attainability mechanism is plausible but unproven; diagnostic = per-gene composition
error e_g = log(mean_c(W_cg / L_c) / CPC_tg) after tilt / relaxation / rounding, lcp off vs on, by expression stratum and on
the real-significant genes; lcp stays off by default; budget +0.001..+0.002 official. Q3: no in-repo modelling direction with
a defensible >= +0.01; effort -> final-phase same-line data acquisition with a frozen S14fk2 fallback; each new source:
independent same-target reproducibility, then ONE held-out same-line full-recipe comparison at replacement weight 0.5.
Rehearsal 8 (`tv2_final_sim8.sh`, audit mode) PASS: `--dual 1.0,3000,0.5,ag,cpc,ci5,lsb,lcp0.25` -> spec dual.lcp 0.25, the
builder accepts it, DE-side deltas bit-identical to rehearsal 7 (lcp touches emission only).
Pin-error diagnostic on the OFFICIAL context A (`ml_lcp_pin_diag.py`, 300 targets, 400 cells, S14fk2 sums / DE deltas,
lsb block, cpc + ci5): the UNROUNDED composition error |e_g| of the top-K genes is 0.000 in every control-CPM decile with
lcp 0 AND 0.25 (median bracket hits 8 of ~9,000 gated genes per target, no 8e5 relaxation on any target); the rounded error
is 0.045-0.083 in the four lowest deciles and 0.003-0.006 in the top three, identical with and without lcp. => on a
400-cell panel the pin IS attained for high-expression genes and lcp does not change the scorer's functional at all; its
OR_X / RN NMAE gains must come from the small-n panels (OR 163, RN 53-96, H1 300 cells per target -> n_top 5-30 cells, so
the down-gene concentration restricts the column support and the pin fails). Confirmation run on the OR_HCT116X and H1
panels themselves (aivc-lcpdiag2-20261001). If confirmed: lcp is a small-n artefact, worth ~0 officially -> stays OFF.
Panel diagnostics (`aivc-lcpdiag2-20261001`, 150 targets each, S14fk2 deltas mapped onto the panel axis): OR_HCT116X (163
cells / target) unrounded top-K |e_g| 0.000 in every decile with lcp 0 and 0.25, median bracket hits 61 (top-K 11) vs 60 (9),
no relaxation; H1 (300 cells) 0.000-0.001 everywhere, bracket hits 130 (top-K 92) vs 127 (87). => the "small-n pin failure"
explanation is REJECTED as well: the composition pin is attained on all three panels and lcp does not change the pinned
functional, yet the real scorer moved OR_X NMAE by -0.042 raw. The scorer's NMAE = mean_g |lfc_pred - lfc_real| / mean_g
|lfc_real| over the REAL-significant genes, divided by a split-half reference (lfc_nmae_ref.py), and its LFC is log2 of
arithmetic per-cell CPM means -- exactly the pinned quantity. The mechanism is therefore UNEXPLAINED; next diagnostic: keep
the two OR_X runs' DE tables (SCH_KEEP) and diff the predicted LFCs of the real-significant genes by expression stratum.
Until then lcp stays OFF by default and its official effect is unpredictable from the harness panels.
Robustness so far (real scorer, OR_HCT116X): seed 1 base 0.3641 -> lcp0.25 0.3734 (+0.0093; NMAE -0.031, reach +0.004) --
the s0 gain (+0.0108) replicates. Dose: lcp0.5 +0.0015 only, because **MSE raw +0.049**: the flatter allocation lowers the
predicted cells' within-group scatter tr Sigma_pred / n_pred below the cap k tr Sigma_real / n_real on some targets, so the
#247/#348 rebate shrinks. lcp0.25 leaves MSE unchanged to 1e-14 on OR_X (both seeds) and H1, but the official contexts'
margin above the cap (H1 1.47x, KOLF 1.18x measured earlier; A/B/C unknown) is a reason for caution before any default-on.
Emission-seed noise on OR_X: base s0 -> s1 moves NMAE -0.012 and reach -0.015 for the same recipe, i.e. the lcp0.25 NMAE
effect (-0.042 / -0.031) is ~3x the seed noise.

**23.5 lcp VERDICT (robustness queue complete, mechanism found): a real but tiny effect that lives in a dozen hyper-expressed
genes; stays OFF.** Real scorer, S14fk structure, lcp0.25 vs base: OR_X +0.0108 (s0) / +0.0093 (s1), H1 +0.0015, RN_Jurkat
+0.0055, RN_HepG2 +0.0073, RN_RPE1 -0.0027 (fid -0.0047, jac -0.0068), KOLF +0.0008; PDS / MSE unchanged to 1e-14 on every
panel at 0.25; lcp0.5 costs MSE +0.049 on OR_X (scatter falls below the #247 cap) and gives nothing on H1 (+0.0012).
Gate 0.5 x (0.75 x mean(OR_X +0.0101, H1 +0.0015) + 0.25 x RN +0.0034) = **+0.0026 < 0.005**; undiscounted mean +0.0050 ->
1/3 transfer +0.0017.
Mechanism (per-target results of the two kept OR_X runs + the cached real DE table): OR_X's NMAE member is computed on only
43 of 263 targets (the de_lfc_nmae gate drops targets with <= 10 real-significant genes; median 24 genes on the scored ones)
and the whole -0.042 comes from a few targets -- SLIRP -0.775 (its 12 real-significant genes are 11 MT-encoded transcripts
plus itself), HMGXB3 -0.19, RBM26 -0.17, MED25 -0.11; median target -0.004. The cpc pin IS unattainable for the hyper-expressed
genes (MT-, CPM >= 10^4): a down-gene with lambda ~ 1 is concentrated into the top-10 % cells, the tilt can only re-weight
inside that support, and the bisection hits its bracket (these are the ~10 "bracket hits in top-K" per target that the decile
medians of 23.4 hid). lcp (factor 0.2-0.3 for those genes) restores the pin there. Officially (400 cells, hundreds of
real-significant genes per target) the same mechanism touches ~10 genes of a few mitochondrial-regulator targets ->
budget +0.001..+0.002 (reach-led, as on H1), and MSE invariance holds only while the predicted scatter stays above the cap
-> NOT a provable ride-along under the user's rule. DECISION: `lcp` stays a dormant option (ported, parity-tested, rehearsal
8 PASS); the 10/22 command is unchanged. A per-gene attainability rule (Codex R24 Q2c) would act on the same dozen genes --
not worth a day.
**Harness lesson (important for every DE-side reading):** on the shallow panels (OR 163 cells, RN 53-96) the NMAE member
rests on a handful of targets with ~10-30 real-significant genes each, so a single target can move the panel NMAE by -0.018
and the emission seed alone moves it by -0.012 (OR_X base s0 vs s1). H1 (300 cells, median 145 real-significant genes) is the
only panel whose NMAE behaves like the official one; discount OR_X / RN NMAE deltas accordingly (this also softens §22.8's
per-panel NMAE arguments, although the S14fk2 official result confirmed that decision independently).

**23.6 H1-truth-trained learned fusion (user "가보자"; design `docs/superpowers/specs/2026-10-01-h1-truth-fusion-design.md`,
Codex R25 corrections: per-fold full-panel construction, masked-softmax weights-only mixer, training-fold standardisation).**
Code `scripts/ml_h1fuse.py` (assembly + exact per-target OOF evaluation on ml_gate_pilot's sum path; base parity vs
delta_s4 8.8e-8; offline base PDS 0.9081 vs the real-scorer base 0.9076) and `scripts/ml_h1fuse_train.py` (torch).
Variant A (pre-registered): logits = log(base weight) + MLP_2x64(9 pilot features), masked softmax over available donors
(K562, HCT116, HEK293T, CD4, KOLF, FENG, HEPG2, RNJURKAT, RPE1), loss = expression-weighted MSE + 1.0 x (1 - own cosine),
5 target-disjoint folds, early stop on 10 % of the training targets (4/5 folds ran to the 60-epoch cap, inner-val loss
-15..-20 %). OOF, full recipe: **PDS 0.8983 (base 0.9081, -0.0098), MSE 0.7197 (base 0.7543, -0.0346), J +0.0145 -> KILL**
(< 0.03). Fusion-only: PDS -0.0144 / MSE -0.093 (J +0.065); per-target NNLS oracle J +0.104. Same signature as every
learned fusion before it: MSE down, PDS down -- the H1 truth as label does not change that.
Variant D (exploratory, pre-registered in the design doc §6 before running): loss = expression-weighted MSE + in-batch InfoNCE
(rows = predictions, cols = truths, tau 0.1, 30 targets per batch) instead of the own cosine. OOF full recipe: **PDS 0.8997
(-0.0085), MSE 0.7244 (-0.0299), J +0.0125 -> KILL**; fusion-only J +0.060. The discrimination loss did not rescue PDS.
=> The H1-truth-trained weights-only mixer is the SIXTH failed source-weighting learner (memory
vcc2026-source-transfer-saturated): training on the official-pipeline truth changes nothing about the PDS-vs-MSE signature.
What remains untested of [another team]'s recipe is the call-set-dependent EMISSION ("called-set amplitude calibration + uncalled-set
pseudobulk pin"), not the fusion.
Frontier check (exploratory, `/home/wonmin-kim/.claude/jobs/43119e26/tmp/h1fuse_ct_sweep.py`): at FULL sum-side centring (ct 1.0)
the offline base reads PDS 0.9134 / MSE 0.7560 (J +0.0098 vs the shipped ct 0.5) and the learned fusions PDS 0.9051 / MSE 0.7218
(A) and 0.9055 / 0.7264 (D): J vs the same-ct base +0.0175 / +0.0135 -- still below 0.03, and the PDS deficit (-0.008) persists
at every centring. The learned mix's PDS loss is not a centring artefact. Diagnostics: own-truth cosine 0.096 -> 0.149 (A),
hubness 0.0034 -> 0.0067 (truth 0.037), fused norm 4.06 -> 1.90. Step 2 / 3 not run. Artefacts: `h1fuse/H1_data.npz`,
`H1_oof_{A,D}.npz`, `H1_eval_{A,D}.json`, `H1_train_{A,D}.json`.

**23.7 Do the top teams all use pseudobulk transfer? (user question 2026-10-01; 8 leaderboard snapshots 09-18 -> 10-01)**
- Described entries (1,206 distinct team x description pairs): pure deep / single-cell models (STATE, scGPT, CellOT, flow
  matching, transformers; n 87) max 0.178, none >= 0.20; pseudobulk-transfer descriptions max 0.242 ([another team]); [another team] 0.272 =
  attention net over donor pseudobulk effects + pseudobulk pin. Every top-20 team with a description is pseudobulk-based.
- Fingerprint: all top-20 have MSE < 0.71 (credit), fid 0.51-0.57, jac 0.03-0.05 = pinned-sum pseudobulk emission with
  over-calling; deep models with their own renderers sit at MSE 5-10. #7/#8/#9/#11/#12 share fid / jac to 4 dp (one pipeline).
- Trajectory: on 09-18 the current top 16 were already at PDS 0.82-0.88 (we are 0.824 today); since then their gains are
  ~46 % MSE credit (+0.030 of +0.065 composite on average), PDS +0.011, DE members +0.024. Several single-step jumps are large
  MSE drops at constant PDS: [another team] 10-01 MSE -0.16 (PDS -0.003), [another team] 09-27 -0.121, [another team] 09-29 -0.104 (PDS -0.0002),
  [another team] 10-01 -0.076. Descriptions name emission-scatter mechanisms ([another team] "dispersed emission template ... jackknife
  credit of the bulk moment rises", [another team] "dispersion head 1.50", [another team] "gamma-Poisson generator, depth-aware amplitude").
- Our check of the cheapest PDS-invariant MSE lever, the budget margin (offline exact, S14fk structure, amplitude scan):
  H1 MSE 0.754 at the shipped 0.8 x kink (r 0.82) -> 0.730 at 1.1 x (r 0.99) -> 0.774 at 1.2 x; OR_X 0.239 -> 0.159 at 1.1 x
  -> 2.02 at 1.2 x (cliff). Moving to the kink is worth ~-0.025 MSE on H1 (~+0.004 composite) with a cliff right above it --
  not the -0.1..-0.16 jumps. Their jumps are either emission-scatter fixes for under-dispersed renderers (ours already
  over-claims: C_pred / C_real 1.47 on H1, the cap binds) or better STRONG-target means, which PDS cannot see (H1 strong-
  tercile PDS is already 0.991; the per-target NNLS oracle gives MSE -0.045 at PDS +0.026).

**23.8 Sum-side per-target AMPLITUDE by cross-family agreement (user "Sum쪽을 개선할수없어?").** PDS is invariant to a
target's own amplitude, so per-target scaling is a pure MSE lever -- the signature of the top teams' PDS-flat MSE jumps.
Offline probe (`scripts/ml_amp_probe.py`, exact evaluator, S14fk structure, budget re-matched): oracle per-target slope is 0
for half of the H1 targets; the best test-time predictor of the oracle slope is the CROSS-FAMILY AGREEMENT (mean pairwise
cosine among donor families measuring the target; Spearman with the oracle slope H1 0.60, RN 0.45-0.58, KOLF 0.19, OR 0.34-0.42;
#families 0.56, |delta| 0.50, S4-vs-fusion cosine 0.25). Rule a_k = rank(agreement_k)^beta, J at beta 0.5 / 1.0:
H1 +0.039 / +0.050, KOLF +0.023 / +0.033, RN_Jurkat +0.030 / 0.000, RN_HepG2 -0.009 / -0.049, RN_RPE1 +0.015 / +0.008,
RN_K562 +0.030 / +0.015; OR / OR_X +0.21 / +0.17 at 0.5 (offline MSE is unreliable on Orion). PDS moves <= +0.005.
Harness token `ta<beta>` (sums only, before the budget match; agreement from `h1fuse/{panel}_data.npz`). PRE-REGISTERED
real-scorer test, beta 0.5 (chosen as the most robust exponent on the non-Orion panels; beta 1 is reported as secondary):
gate = both OR_X seeds > 0, H1 > 0, RN mean >= 0 and every RN >= -0.003, 0.5 x (0.75 x mean(OR_X, H1) + 0.25 x RN mean)
>= 0.005. Panels: H1, OR_HCT116X (s0, s1), RN_Jurkat, RN_HepG2, RN_RPE1, RN_K562, KOLF.
REAL SCORER (S14fk structure, ta0.5 vs base; RN_K562 on the arm without the K562-denoised DE source, which is a no-op there):
  H1 0.3089 -> 0.3155 **+0.0065** (PDS +0.0008, MSE -0.0357; offline predicted J +0.039 = +0.0065, exact)
  OR_X s0 / s1 **+0.0215 / +0.0206** (MSE -0.156, NMAE +0.027..+0.030 worse, PDS +0.0012)
  RN_Jurkat +0.0045, RN_HepG2 **-0.0045** (MSE +0.014), RN_RPE1 +0.0025, RN_K562 +0.0047 -> RN mean +0.0018; KOLF +0.0013
  secondary beta 1.0: H1 +0.0087 (MSE -0.049), OR_X +0.0174
Gate: OR_X both > 0 yes, H1 > 0 yes, RN mean >= 0 yes, R21 composite 0.5 x (0.75 x mean(OR_X +0.0211, H1 +0.0065) + 0.25 x
+0.0018) = **+0.0054 >= 0.005 yes**, every RN >= -0.003 **NO (HepG2 -0.0045)** -> pre-registered verdict FAIL on the per-panel
floor alone. 7 of 8 panels positive; MSE down on 7 of 8. The R21 composite predicted the last two official results within
0.0005 (S14fk +0.0065 vs +0.0060, S14fk2 +0.003 vs +0.0027) -> expected official ~ +0.005. HepG2 is an essential-gene panel
(0/300 official targets are essential); the official agreement distribution (median 0.0275) matches H1's (0.0287).
Ported, default off: `tv2_tamp_prep.py` (parity with the harness token on OR_X targets 8.8e-8), builder `sum.tamp`
{file, beta}, `tv2_final.py --sum-tamp BETA`; official agreement file `reports/submission_2026-10-01/tamp_official.npz`.
**S15 BUILT (user choice 2026-10-02 "S15 패키지 빌드만, 업로드는 나중에"; NOT uploaded):** `reports/submission_2026-10-01/spec_S15.json`
= spec_S14fk2 + `sum.tamp {file: tamp_official.npz, beta 0.5}`. First build had a defect caught by the audit: the builder took
the adaptive-ew decision on the tamp-SHRUNK sums, which flipped ew back ON (kink with ew 1.08-1.14 >= 1), i.e. a configuration
never tested; fixed in `tv2_submit_build.py` (the ew decision uses the pre-tamp sums; the no-tamp path is unchanged), wrong
artefacts deleted, rebuilt. Rebuild audit: ew dropped on A/B/C exactly as S14fk2 (0.795 / 0.814 / 0.843 < 1); dsum = S14fk2
dsum x a_t to 1.7e-7; DE delta, main delta and overrides bit-identical to S14fk2; amplitude median 0.708; budget re-match
gamma_sum 1.325 / 1.338 / 1.328 (S14fk2 0.978 / 0.989 / 0.983). vcc prep passed (360,000 cells, targets verified, counts
preserved). Package `reports/submission_2026-09-23/S15.vcc`, sha256 4b990cc8ebf56951149845ee8bd351e729c62916068a222e44c1ed9081d6b1c2,
pred.h5ad deleted. Expected official ~ +0.005 (R21 composite +0.0054); HepG2 floor caveat (§23.8).
Rehearsal 9 (`scripts/tv2_final_sim9.sh`, sim7 intended + `--sum-tamp 0.5`, audit mode) PASS: per-context agreement prep ran
inside tv2_final (somatic D/E/F median 0.045 without KOLF; stem G 0.040 with KOLF), ew dropped on all four contexts with the
pre-tamp decision, DE-side deltas bit-identical to rehearsal 7. `--sum-tamp` stays OFF in the recommended 10/22 command until
S15's official result (or the user) decides.

**23.9 Where our PDS is lost (user "PDS가 차이나는이유가뭐야"; H1 harness, S14fk structure, per-target PDS of the offline base).**
Strength terciles: weak 0.769 / mid 0.964 / strong 0.991 -> the weak third carries 84 % of the rank loss; the 30 worst targets
(PDS 0.44) carry 61 % of it, 97 % of them below median strength, cross-family agreement 0.007 (all targets 0.029). By source
coverage: targets measured by 4 families PDS 0.775, 5 families 0.877, 8-9 families 0.986 (confounded: the 8-9-family targets are
essential genes, i.e. strong). The OFFICIAL targets are measured by 2 / 3 / 4 / 5 / 6 of the shipped somatic families for
7 / 25 / 145 / 108 / 15 targets -> the 4-5-family regime, whose H1 PDS (~0.82-0.84 weighted) matches our official 0.824.
Unused coverage left: DLD1 67 official targets, H1GW 66, VIP 118 (worse than Replogle K562). Our official PDS history: dual
emission +0.046 (S8 -> S10b), S4 +0.015 (S10e), no-KOLF sum side +0.011 (S10k), ensemble +0.001, FENG-out S4 +0.002.

**23.10 S15 OFFICIAL: 0.24603, rank 32 (entry jvcHmRtUHN9TpQyD2rHb, b24-4b990c; user go-ahead "업로드해봐", 2026-10-02 09:44Z) --
+0.01060 over S14fk2, the largest single gain since S10c.** Raw: PDS 0.8238 -> 0.8258 (+0.0020), MSE 0.7493 -> 0.6962 (-0.0531),
NMAE +0.0006, fid +0.0001, reach +0.0030, jac +0.0001. Composite by member: MSE +0.0093, PDS +0.0008, reach +0.0006, NMAE -0.0001.
Forecast was ~+0.005 (R21 composite +0.0054; H1 harness +0.0065 with MSE -0.036): the official MSE moved 1.5x the H1 harness.
=> The PDS-invariant per-target amplitude lever transfers; the pre-registered HepG2 floor was too conservative for this lever
(essential-gene panel). Final phase: `--sum-tamp 0.5` becomes part of the recommended 10/22 command (rehearsal 9 PASS).
Open follow-ups: beta (1.0 was better on H1, +0.0087, worse on OR_X and offline RN), combining tamp with the PDS-raising but
MSE-costly levers (ct 1.0, K562 x2) that were rejected for their MSE cost.

**23.11 Next candidates after S15 (user "다음 후보 진행해줘"): beta x centring x K562 weight grid** (`scripts/ml_combo_probe.py`,
offline exact, S14fk structure, budget fixed; J vs S15). Mean over the non-Orion panels (H1, KOLF, RN x4): S15 is the BEST cell
of 24. Split by panel type: non-essential panels prefer more shrink / centring (H1 vs S15: beta 0.75 +0.006, beta 1.0 +0.012,
ct 1.0 + beta 1.0 +0.020; KOLF +0.008 / +0.010 / +0.013), the essential RN panels prefer less (beta 1.0 RN mean -0.023, ct 1.0
-0.06). K562 x2 loses on 5/6 (H1 -0.006). Real scorer already in hand: beta 1.0 vs beta 0.5 = H1 +0.0022, OR_X -0.0041 -> ~0
(offline OR_X had claimed +0.077: Orion's offline MSE is unreliable, again). The official S15 gain tracked H1 / OR_X, not RN,
so the non-essential panels decide for the validation set; the final panel's composition is unknown, so RN is the risk read.
PRE-REGISTERED real-scorer check (H1 + OR_X, + RN_Jurkat / RN_HepG2 as the risk read), arms vs the S15 arm (ta0.5):
ct1.0+ta0.5, ct1.0+ta1.0, ta0.75. Adopt only if mean(H1, OR_X) >= +0.005 composite vs S15 (the slot rule) and neither is negative;
otherwise S15 stays the final configuration.
REAL SCORER vs the S15 arm (ta0.5): ct1.0+ta0.5 H1 +0.0016 (PDS +0.0052, MSE +0.0014), OR_X +0.0008, RN_Jurkat -0.0100, RN_HepG2
-0.0072; ct1.0+ta1.0 H1 +0.0037, OR_X -0.0040 (MSE -0.071 but NMAE +0.047), RN -0.012 / -0.014; ta0.75 H1 +0.0016, OR_X +0.0037
(NMAE +0.016), RN -0.0012 / -0.0023; ta1.0 H1 +0.0021, OR_X -0.0041. mean(H1, OR_X): +0.0012 / -0.0001 / +0.0026 / -0.0010 ->
**none reaches the pre-registered +0.005; S15 (beta 0.5, ct 0.5, K562 x1) stays the final configuration.** Mechanism of the
ceiling: stronger shrink keeps lowering MSE but the shrunk targets' emitted column totals start to contradict the DE side
(OR_X NMAE +0.016 at beta 0.75, +0.047 at beta 1.0); full centring buys H1 PDS (+0.005) but costs RN PDS (-0.006..-0.009).

**23.12 3-hour multi-angle loop (user "다각도 개선시도 테스트 루프 3시간 돌려줘", 2026-10-02 10:00-13:00 PDT).** Every arm relative to
S15. Offline exact screens (`scripts/ml_multi_probe.py`, `scripts/ml_agree_probe.py`, 8 panels; J on H1 / KOLF / RN x4, PDS only on Orion):
- A amplitude signal: composite rank(rank(agree) + rank(#families) + rank(|delta|))^0.75 mean J +0.0095 (H1 +0.013, KOLF +0.017,
  RN_Jurkat -0.0075, RN_HepG2 +0.026; OR_X / OR dPDS +0.0064 / +0.0027) -> real-scorer candidate (token `tcomp`); same at beta 0.5
  +0.0078 (H1 +0.001); agree x nfam^0.25 +0.0071 (RN -0.02); value-based -0.002; agree + norm +0.004.
- B gene-level shrink by cross-family sign agreement: -0.006 / -0.020 / -0.059 (gamma 0.5 / 1 / 2) -> closed.
- D per-target neighbour-blend weight by agreement: -0.014 / -0.039 -> closed.
- E per-target S4 share by agreement: low agreement -> S4 +0.011 mean but H1 -0.018 (RN-driven); the opposite -0.020 -> closed.
- G budget margin with tamp: 0.9 x kink H1 +0.021 but KOLF -0.134, RN_HepG2 -0.098 (cliff); 1.0 x kink -0.29 -> closed.
- J/K agreement estimate: + KOLF voter -0.001, + KOLF / DLD1 / H1GW / VIP voters -0.002, cell-weighted pairs -0.0005, min_cells 50
  -0.006, median of pairs -0.003 -> the S15 agreement is already the best estimate; closed.
Real scorer queued: C (DE-side amplitude `tad0.5` on top of S15) and A (`ta0.75_tcomp`) on H1, OR_X, RN_Jurkat, RN_HepG2.
- H learned amplitude signal (`scripts/ml_amp_learn.py`): 8 test-time per-target features (agreement, #families, |delta|, S4-vs-fusion
  cosine, family cells, max family KD, target CPM in the destination, #promoter neighbours -- the last is constant 0, a harmless
  feature bug); label = rank of the per-target TRUE cosine (scorer space) under S15; leave-one-PANEL-out GBM (200 trees, depth 3).
  Spearman(pred, true cos) beats the agreement alone on 8/8 panels (H1 0.654 vs 0.631, OR_X 0.563 vs 0.385, RN +0.06..+0.13);
  importances agree 0.41, |delta| 0.31. Offline J vs S15 at beta 0.75: positive on 8/8 (H1 +0.006, KOLF +0.010, RN +0.003..+0.024).
  REAL SCORER vs the S15 arm (`ta0.75_tlrn`, LOPO predictions per panel): **H1 +0.0014 (PDS +0.0012), OR_X +0.0093 (PDS +0.0031,
  MSE -0.061), RN_Jurkat +0.0007, RN_HepG2 +0.0049 (PDS +0.0041)** -> 4/4 positive, PDS up on all four, mean(H1, OR_X) +0.0053 >=
  0.005 -> ADOPT per the pre-registered rule. Robustness runs (OR_X seed 1, KOLF, RN_RPE1, RN_K562) queued.
- C DE-side amplitude (`ta0.5_tad0.5`): H1 -0.0066, OR_X -0.0057, RN -0.0070 / -0.0096 (NMAE +0.009..+0.032, reach down) -> closed.
Robustness (real scorer vs S15): OR_X seed 1 +0.0070, KOLF +0.0022, RN_RPE1 +0.0015, RN_K562 +0.0012 -> the learned amplitude beats
S15 on 8 / 8 runs, PDS up on all 8. (Composite rank `tcomp` H1 +0.0028; not pursued -- dominated offline.)
Known flaw (kept for consistency with the tested runs): per-panel rank normalisation breaks ties by target-list position, so the
constant #promoter-neighbours column and ties in #families become a list-position signal (importance 0.048); fix with average ranks
if the model is ever retrained.

**23.13 S16 OFFICIAL: 0.24742 (entry lDWUfeHo4pEpR42pEc17, b24-92d6ca, rank 33 at retrieval; user "개선이 있어보이는건 제출해줘 하나",
uploaded 2026-10-03 00:02Z) -- +0.00139 over S15.** S16 = S15 with the amplitude ranked by the learned true-cosine prediction
(`scripts/tv2_tlrn_prep.py`: GBM trained on all 8 harness panels; official features from the context's family means, controls and
the builder's new `*_dbase_<ctx>.npy` dump of the pre-resid sum override; predictions written in the `sum.tamp` format) at beta 0.75.
Audit: ew dropped on A/B/C as S14fk2, dsum = S14fk2 x rank(pred)^0.75 to 2.2e-7, DE / main delta / overrides bit-identical; sha256
92d6ca5bc9b593b7563aa9c4a5b60fb349d6871135da4c37fc469139965ad75e. Official members: PDS 0.8258 -> 0.8275 (+0.0016), MSE 0.6962 -> 0.6922
(-0.0040), reach +0.0010, NMAE +0.0004. The H1 harness predicted exactly this (+0.0014, PDS +0.0012); OR_X (+0.008) over-predicted.
Final phase: the 10/22 command keeps `--sum-tamp 0.5` (agreement rank); the learned variant needs a two-pass build (audit pass for
the dbase dump, then tv2_tlrn_prep.py) for +0.0014 -- optional, not wired into tv2_final.py.

**23.14 Where to improve next, for score AND final-phase generalisation (user "어느방향으로 더 개선해야할지 ... 최종 일반화에 도움이될지 탐색").**
Current recipe (S16-equivalent arm) per harness panel, composite: H1 0.317, OR_X 0.394, RN 0.377-0.489, KOLF 0.065; official 0.247.
The official member profile is H1's (MSE 0.69 vs 0.70, reach 0.295 vs 0.288) except PDS (0.827 vs 0.910, -0.030 composite), NMAE
(0.904 vs 0.858, -0.013), jac (0.031 vs 0.073, -0.019) and fid (-0.011). KOLF (iPSC, low signal) is the panel type we serve worst.
- Amplitude lever, remaining offline headroom after S16 (true-cosine rule beta 1 minus S15 minus S16, J): H1 +0.052 (~+0.009
  composite), RN_RPE1 +0.031, RN_K562 +0.063, RN_Jurkat / HepG2 +0.10, KOLF +0.11. It is the only lever family that transferred twice
  (S15 +0.0106, S16 +0.0014 = the H1 forecast) and it is learned across 8 panels, i.e. line-agnostic -> best score x generality bet.
  Untried predictors of a target's direction accuracy: DepMap essentiality of the target (line-independent), the real #promoter
  neighbours (the S16 feature was constant 0), S4-vs-fusion agreement on the top genes, source knockdown per family; plus the
  average-rank fix.
- TF-regulon prior: NOT a lever -- TF targets are not the weak ones (H1 PDS TF 0.908 = non-TF 0.908; ZNF-like 0.855 but 11 targets);
  official 75 TFs, 35 of them ZNF-like with unknown regulons.
- Same-line / same-lineage data for D/E/F: largest single lever when it exists (H1GW on H1 +0.074) but conditional on the lines.
- Final-pipeline insurance: no SCORED end-to-end pseudo-final has been run (rehearsals 5-9 were audit-only); a silent pipeline error
  is the largest downside risk of the final phase.

**23.15 Arc's 2026-10-01 data post vs our pipeline (user pasted it: "우리 이거 사용하고있는거야?").** Already used: 400 cells per
perturbation (emission), 10x Flex (H1-2025 Flex harness is primary), ~20k UMI / cell (H1 harness thinned to 20k), >= 80 % on-target KD
curation (official targets are non-essential, KD override), NTC baseline (official controls). The post cites the Replogle 2022
dual-sgRNA library -> the same guides as K562 GWPS, consistent with K562 being our best / heaviest source. New information: the six
lines span tissue, lineage and DISEASE STATUS with stem / immortalized (hTERT) / cancer as the illustrated classes, none is H1, and
A/B/C (Jurkat, HeLa, CAL-33) are all cancer lines -> D/E/F likely include immortalized and possibly stem (non-H1) lines. Checks:
HPA (identity reference) contains hTERT-RPE1, BJ-hTERT, HaCaT, RPTEC/TERT1, HUVEC/TERT2 but NO iPSC / ESC line; SAME_LINE_SRC covers
K562, HCT116, HEK293T, DLD1 (HCT-15), HepG2, RPE1 -- KOLF2.1J (our genome-scale KOLF screen, 282/300 official targets) is not
registered. Matching destination controls against our sources' control profiles is platform-confounded (official A/B/C match the
Flex sources CD4 / H1 / VIP first: A CD4 0.64, B H1 0.55, C H1 0.59) -> not usable naively. ProPer-seq (Swinderman 2026, Flex
Perturb-seq): no public GEO deposit found. The other cited datasets (Allen IL-6, cytokine dictionary, Tahoe) are not genetic knockdowns.

**23.16 Plan "1->2->3 전부" (user, 2026-10-03).**
(1) Amplitude predictor v2 (`scripts/ml_amp_learn2.py`): + DepMap essentiality (mean gene effect; Spearman with the true cosine
-0.34..-0.51 on every panel), the real #promoter neighbours, S4-vs-fusion cosine on the top-300 genes, mean family KD; average-rank
ties. LOPO Spearman with the true cosine up on 6/8 (H1 0.697 vs 0.654) but J vs S16 at beta 0.75: H1 +0.005, KOLF +0.004, RN_Jurkat
-0.003, RN_HepG2 -0.002, RN_RPE1 +0.008, RN_K562 -0.006, OR_X +0.001, OR +0.007 (mean ~+0.002 J, ~+0.0003 composite) -> the
amplitude lever is saturated for test-time predictors; NO S17.
(3) KOLF2.1J identity for the final phase: expression correlation across platforms fails (H1 Flex controls match the HipSci pool
0.767 before H1GW 0.702), but a presence / absence fingerprint works: KOLF2.1J has DDX3Y / UTY / USP9Y (224-528 CPM) and NO KDM5D /
EIF1AY (0.0), unlike H1 / H1GW (KDM5D 8-40, EIF1AY 60-186), official A / C (male: KDM5D 52 / 39) and B (female, all Y = 0).
`tv2_final.kolf_signature()` (mean on-genes >= 20 CPM, both off-genes < 1 CPM) -> in the stem branch the same-line source becomes
our genome-scale KOLF screen at the full --same-line weight (KD gate -0.93 <= -0.8 passes) instead of the H1GW fallback.
Unit test: official A / B / C False, H1 controls False, KOLF controls True. Gap noted: same-line RN sources (hTERT-RPE1, HepG2)
exist only with --rn-sources (essential-gene screens; small value on an official-like panel).
(2) Scored holdout pseudo-final (H1 as a new context G, H1 + H1GW removed from sources (TV2_SRC_DIR), neighbour library,
library PCs (TV2_HOLDOUT_FAMS) and S4 (leave-H1-out refit); same-line off) running as aivc-hold-h1-20261003; caveat: the
no-embedding residual model was trained with H1 among its pseudo-panels (alpha 0.25, small).

## 24. 12-hour improvement loop, 2026-10-03 03:50- (user: "12시간동안 개선 루프 돌리면서 실행해줘"); plan LOOP_PLAN_2026-10-03.md

**24.0 Operational.** Starting the holdout scorer and an offline probe together tripped systemd-oomd at 03:54 (user session
73 % > 50 %): both were killed, and Discord with them. The whole loop now runs as ONE strictly serial chain
(`aivc-loopchain-20261003`), and the scorer queue waits for that unit to exit. The holdout scorer was also rewritten to
stream the build's cells onto the panel axis in 20,000-cell chunks and to score in a separate process
(`hold_score2.py prep|score`), so the 6 GB prediction and the 16 GB scorer never coexist.

**24.1 SCORED holdout pseudo-final (the insurance item, first time end to end): -0.0298 vs the harness arm, ALL of it on the
SUM side.** `tv2_final.py` ran the recommended 10/22 command on the H1 panel as a new context G, with H1 and H1GW removed
from the sources (`TV2_SRC_DIR`), the neighbour library, the library PCs (`TV2_HOLDOUT_FAMS`, new env in `tv2_pc_prep.py`)
and S4 (leave-H1-out refit); same-line off. Real scorer vs the S15 harness arm on the same panel:
  composite 0.2857 vs 0.3155 | PDS 0.8843 vs 0.9084 (-0.0241) | MSE 0.8361 vs 0.7080 (+0.1281)
  NMAE 0.8663 vs 0.8584 | fid 0.5502 vs 0.5454 (+) | reach 0.2891 vs 0.2872 (+) | jac 0.0748 vs 0.0728 (+)
The DE side is fine; the rebuilt SUM-side inputs are not. The pipeline's own kink estimate also read 0.756 vs the harness's
~0.90. Attribution running (`scripts/ml_hold_ablate.py`: swap the holdout's neighbour KNN / library PCs into the harness
recipe, and score the pipeline's own pre-gamma sums, all budget-matched to one reference).
Caveat: part of the gap is by construction (the pipeline builds a top-3 neighbour library for the new panel from ~1,000
library targets, the harness uses its 2,581-target library) and the no-embedding residual saw H1 among its pseudo-panels.

**24.2 Offline screen 1 (angle: levers that cost MSE were priced BEFORE the amplitude rule existed).** `ml_loop1_probe.py`,
J vs the S16 arm, 8 panels.
- **kink factor 1.17 -> 1.35 is positive on 8/8 panels**: core mean J +0.0238 (H1 +0.0041, KOLF +0.0042, RN_Jurkat +0.0494,
  RN_HepG2 +0.0337, RN_RPE1 +0.0214, RN_K562 +0.0298), OR_X +0.0768, OR +0.0915; dPDS ~ -0.0002 everywhere, so it is a pure
  MSE/budget move. Factor 1.0 is negative on 6/8. -> real scorer (queued).
- hubness `dis`: buys PDS on every panel (H1 +0.001..+0.004, OR_X +0.004..+0.012) but the MSE cost is large on RN
  (HepG2 -0.04..-0.14) -> core mean negative at every (k, lambda). Closed again, now with tamp in place.
- pcr k / lambda: mixed and panel-type-specific (RN sees no change: the guard zeroes lambda there). blend 0.4 +0.0028 mean,
  blend 0.2 / 0.5 negative. Not worth a slot.

**24.3 Offline screen 2 (angle: new ensemble axes on the TRANSFORM side) -- CLOSED.** `ml_loop2_probe.py`: averaging the
transformed sums over pcr lambda or over the centring fraction is EXACTLY a no-op (both are linear in their parameter, so
the average equals the mean parameter: measured 0.0000 on 8/8, a useful check that the probe is right). The nonlinear axes
are no better: pcr k {1,3,5} core -0.0096, amplitude beta {0.5,0.75,1} +0.0017, blend {0.2,0.3,0.4} +0.0008, all four
together -0.0054. The shipped min_cells ensemble stays the only ensemble that pays.

**24.4 New harness token `drc<a>x<b>`** (angle 3, from leaderboard entry v17a's public description, team [another team] 0.2423):
within each target the top 5 % of genes by |DE lfc| are scaled by a and the bottom 70 % by b, the BULK (sum) target
unchanged -> PDS / MSE invariant by construction, NMAE / reach / fid / jac move. Queued at [another team]'s own 0.64 / 1.18.

**24.5 kink factor on H1: EXACTLY zero, and the offline probe was mis-referenced.** Real scorer, `gk1.35` vs `gk` (1.17) on the
S16 arm, H1: every member identical to 4 dp (composite 0.3169 both). Reason: the arm reports `kink_base 0.974` at factor 1.35
and `gamma_base = max(1.0, kink_base) = 1.0`, i.e. the floor binds at both factors, so the budget reference
`budget_B(g0 * D)(1.0)` is bit-identical. The offline probe said +0.0041 because it recomputed the kink from the ENSEMBLED
total while the arm uses the `bmb` reference (`budget_ref: 'old'`, the non-ensembled main recipe) -- the same reference trap
as FINDINGS 17.0. Lesson for every future budget arm: price it against the arm's own budget reference, not the total.
This does NOT make the factor irrelevant for the official build: A/B/C report kink_base 0.893 / 0.937 / 0.960 at factor 1.17,
so factor 1.35 would lift gamma_base to ~1.03-1.11 there (the floor stops binding). The panels where the floor does not bind
(OR_X, RN) are therefore the relevant evidence; those runs are queued.
**24.6 kink factor 1.35 KILLED on the real scorer -- it is the cliff.** OR_HCT116X (the floor does not bind there:
kink_base 1.646 -> gamma_sum 2.302 at factor 1.17): raising the factor to 1.35 moves gamma_sum to 2.471 (+7 %) and
**MSE 0.0299 -> 1.2108** (composite 0.3942 -> 0.2130, -0.181); PDS -0.0012, NMAE +0.049. The offline probe had claimed
+0.0768 on this very panel. Two independent reasons it was wrong: the probe's budget reference (the ensembled total, not
the arm's `bmb` 'old' reference) and Orion's offline MSE being unreliable in general. The shipped margin 0.8 is what keeps
the sums below this cliff; the official contexts (kink_base 0.89-0.96, floor binding at gamma_base 1.0) sit safely below it.
RULE: never price a budget / amplitude-scale arm offline again -- go straight to the real scorer, and only downward.
RN_Jurkat (also floor-free, gamma_sum 1.794 -> 1.908) instead GAINS: composite 0.4892 -> 0.4968 (+0.0075, MSE 0.2887 ->
0.2398). So the cliff is panel-specific: the same +7 % budget is free headroom on RN_Jurkat and a -0.18 catastrophe on OR_X.
With the official contexts sitting right at the floor (gamma_base 1.0 at factor 1.17, ~1.03-1.11 at 1.35) there is no way to
know which side they are on without a slot -> NOT pursued. The margin 0.8 stays.

**24.7 HOLDOUT PSEUDO-FINAL PASSES: the -0.027 was the rehearsal's control cells, not the pipeline.** `ml_hold_pin.py`
evaluates the pipeline's own sums at its own gamma (1.294) under four input combinations:
  harness override + harness controls   pds 0.9082  mse 0.6929
  PIPELINE override + harness controls  pds 0.9082  mse 0.6929   (override is irrelevant; both have 1 override gene / target)
  harness override + PIPELINE controls  pds 0.8852  mse 0.8233
  PIPELINE override + PIPELINE controls pds 0.8852  mse 0.8233   == the actual run (0.8848 / 0.8245)
The rehearsal fed `tv2_final.py` the sim5 context_G control block, a DIFFERENT 12,000-cell draw of H1 NTCs than the one the
H1 panel's real perturbed data is paired with (profile cosine 0.99997, mean UMI 21,394 vs 21,390 -- yet worth -0.023 PDS and
+0.130 MSE). The emission pins every gene's group sum to `control_mean x exp(lfc) x n`, so the control draw's sampling noise
enters every prediction. In the real final phase D/E/F's own controls come with the panel, so this cannot happen.
**With the correct controls the final pipeline MATCHES and slightly beats the harness arm on a held-out panel**
(PDS 0.9082 vs 0.9084, MSE 0.6929 vs 0.7080, i.e. J +0.016), its own budget choice (gamma_sum 1.294) is the measured MSE
optimum, and its rebuilt neighbour library / library PCs are fine (24.1). The insurance item is delivered: no silent defect.
Residual caveats: the no-embedding residual saw H1 among its pseudo-panels, and `--same-line` was off.
Operational note for 10/22: emit with the CONTROLS THE ORGANISERS SHIP for that context, never a resampled block; and when
rehearsing on a harness panel, pass that panel's own control cells (`TV2_CELLS` now also lets the rehearsal match n).

**24.8 `drc` (rank-class DE amplitude calibration, v17a's published rule 0.64 / 1.18) -- real scorer.** PDS and MSE are
bit-identical on both panels, as designed. H1 **+0.0013** (NMAE 0.8583 -> 0.8555, fid +0.0009, jac +0.0002, reach -0.0008);
OR_HCT116X **-0.0041** (NMAE +0.0137, reach -0.0029). mean(H1, OR_X) -0.0014, below the +0.005 bar. OR_X's NMAE rests on
43 targets x ~24 real-significant genes (23.5), so its +0.0137 is inside that panel's noise band, and H1 is the panel whose
NMAE behaves like the official one -- but +0.0013 alone is a ride-along at best. Strength sweep queued
(0.5/1.3, 0.8/1.1) plus the two single-sided controls (top-only 0.64, bottom-only 1.18) to see which half carries it.

**24.9 Emission-knob sweep on H1 (real scorer, all PDS / MSE invariant by construction; d vs the S16 arm 0.3169).**
  drc 0.8 / 1.1      **+0.0018**  (NMAE -0.0042, reach +0.0014, fid +0.0005)   <- best drc strength
  drc 1.0 / 1.18     +0.0013  (bottom-70 % inflation ALONE)
  drc 0.64 / 1.18    +0.0013  (v17a's published pair)
  drc 0.64 / 1.0     -0.0006  (top-5 % shrink ALONE)  -> the gain is the BOTTOM inflation, not the top shrink
  drc 0.5 / 1.3      -0.0015  (too strong)
  **K 5000           +0.0022**  (fid +0.0053, jac +0.0057, NMAE -0.0011, reach -0.0194)
  K 2000             -0.0084
  no `agc` (rank by plain |d|)  -0.0002   -> the agreement-weighted ranking is right; closed
  `gdf` 0.5 / 1.0    exactly 0.0000 -> inert: with `gk`+`bm`+a separate DE delta the builder takes the branch that leaves
                     LdDE unscaled, so the token never fires. Not a result about the idea; it was never applied.
  nt 0.05 / 0.2      -0.0004 / -0.0006 ; le 0.1 / 0.5  -0.0011 / -0.0013  -> the shipped 0.1 / 0.3 are at their optimum.
Seed context: OR_X base vs base@s1 is -0.0024, so single-panel deltas below ~0.003 are noise-sized. H1 seed pair queued.
Cross-panel checks queued for the two positives (drc 0.8/1.1 and K 5000) plus their combination.

**24.10 Exhaustive public same-line data search (2026-10-03), after the user judged the knob gains too small.** Subagent
swept the GEO DataSets API (~750 series), Europe PMC 2024-2026, ArrayExpress/BioStudies, scPerturb, CELLxGENE Discover,
the IGVF portal and HuggingFace, for transcriptome-wide low-MOI single-cell CRISPR screens in stem / hTERT-immortalized /
our cancer lines. Verified by me afterwards:
- **X-Atlas/Pisces status upgraded from "unreleased" to "announced, upload pending":** the HF repo
  `Xaira-Therapeutics/X-Atlas-Pisces` now EXISTS (public, CC-BY-NC-SA 4.0) but holds only README + LICENSE, unchanged since
  2026-03-25. The README promises all Jurkat Resting + Jurkat Active cells (10,872 / 10,878 perturbations, 2.8M cells each,
  clean dual-guide assignment) and HepG2 / iPSC restricted to a 200-gene validation set. Jurkat = validation context A, so a
  release helps the LEADERBOARD, not the final ranking; the iPSC subset's 200 genes are useless on a 300-target panel.
- **Two new downloadable stem screens, both useless as fusion sources:** GSE176038 (WTC-11 iPSC CRISPRi, 252 targets, 63k
  cells, mtx on FTP) covers **0-3** of our official 300; GSE343369 (KOLF2.2J CRISPRi, 230 targets, h5ad, explicitly low MOI)
  covers **7**. Chance-level, exactly as the 09-29 round predicted for small screens. Only GENOME-SCALE screens can serve a
  300-target panel -- which is why our KOLF2.1J / H1GW / FENG / Replogle assets matter and these do not.
- **Confirmed absent, search now exhaustive:** HeLa (zero single-cell CRISPR screens in GEO), CAL-33 / any HNSCC (zero),
  HaCaT / keratinocyte (zero), BJ-hTERT / RPTEC-TERT1 / HUVEC-TERT2 / hTERT-HPNE (zero), hTERT-RPE1 (nothing beyond the
  Nadig essential screen we hold), H9 / WA01 / HipSci (nothing beyond GSE295214 + Feng). ENCODE's many 2026 Jurkat and
  WTC11 screens are Flow-FISH or growth-only; IGVF's 329 Perturb-seq sets are targeted cCRE panels, not gene knockdown.
=> There is no new public data lever for either phase today. The one asset that would matter (Pisces Jurkat) is not
downloadable; watcher kept on the HF repo.

**24.11 Cross-panel close-out of the two positives (real scorer).** H1 seed repeat of the base arm: **+0.0002**, so H1
deltas of ~0.002 are real, not emission noise (OR_X's seed pair is -0.0024, so that panel needs ~0.005 to be readable).
  `drc 0.8 / 1.1`         H1 +0.0018 | OR_X +0.0004 | RN_Jurkat +0.0005   -> **positive on 3/3**, PDS / MSE invariant
  `K 5000`                H1 +0.0022 | OR_X -0.0055                        -> panel-split, NOT adopted
  `K 5000 + drc 0.8/1.1`  H1 +0.0034 | OR_X -0.0048                        -> the combination inherits K's OR_X loss
=> The loop's only surviving scoring lever is `drc 0.8 / 1.1`, worth ~+0.001 officially: a ride-along for the next
substantial build, never a slot of its own (the +0.005 rule). Everything else tested this loop is closed.

**24.12 LOOP VERDICT (12 h, 2026-10-03 03:50-15:50).** The scoring side of the validation task is saturated at the
+0.002 level: 20 real-scorer arms and 2 eight-panel offline screens produced one +0.001 ride-along. The loop's real
deliverable is 24.7: the final-phase pipeline is validated end to end against a held-out panel's real scorer, its rebuilt
inputs are sound, its budget choice is the measured optimum, and the one thing that looked like a -0.030 defect was the
rehearsal feeding it the wrong control block. Second deliverable: 24.10, an exhaustive public-data search that closes the
"maybe there is a same-line screen we missed" question for HeLa, CAL-33, every hTERT line and every stem line -- there is
not, and the only asset that would matter (X-Atlas/Pisces Jurkat) is announced but still unpublished.

## 25. Aggressive directions, 2026-10-03 (user: "공격적인 시도를 해볼수있는곳은" -> "ㄱㄱ")

**25.1 Emission DISPERSION is already saturated -- option killed in 10 minutes.** The idea ([another team]'s public v16a/v17a
description: a mean-preserving per-cell spread so "the jackknife credit of the bulk moment rises") requires our claimed
correction to sit BELOW the #247 cap. Measured on a kept OR_HCT116X run (`cell_eval2/results.csv`, 263 targets):
`expr_mse_unbiased_capped` > `expr_mse_unbiased` on **263/263 targets**, i.e. the cap clips our claim on every row --
our emitted per-cell dispersion already exceeds k x trSigma_real/n_real everywhere. Algebra of the member:
total rebate = min(B_pred, SUM_p claim_p) with claim_p = min(trSigma_pred/n_pred, k trSigma_real,p/n_real). Our claim is
at its per-row maximum, so the binding term is B_pred -- the sum amplitude, which the kink already optimises (r = 0.914 at
the measured MSE optimum, 23.x/24.x). Raising dispersion changes nothing; lowering it loses. CLOSED.

**25.2 What the emission can be worth AT ALL (leaderboard 2026-10-03, #10 = 0.2823, us 0.24742 rank 34, gap 0.0349).**
Gap by member in composite units: PDS +0.0133, MSE +0.0107, NMAE +0.0071, reach +0.0020, jac +0.0011, fid +0.0006.
**The four DE members -- everything the emission can touch with the sums pinned -- total +0.0109; PDS + MSE (the delta
itself) is +0.0240.** So a PERFECT emission reaches ~0.258 and ~rank 21. Within that, NMAE is two thirds of the prize.
=> The aggressive play inside the emission is a magnitude calibration of the DE side aimed at NMAE; everything beyond
+0.011 requires a better delta, which is measured-saturated without new data (24.10).

**25.3 New token `drp<gamma>`: continuous DE rank calibration.** Generalises `drc`'s two rank classes to a smooth monotone
curve: within each target, multiply the DE lfc by (descending-|lfc| rank / N)^gamma, normalised to median 1. The bulk
target is untouched, so PDS / MSE stay invariant (verified for drc, 24.8). gamma 0 = shipped; gamma > 0 inflates the weak
tail, which is the half of `drc` that carried its gain (bottom-70 % inflation +0.0013 alone, top-5 % shrink -0.0006).
Sweep gamma 0.10 / 0.20 / 0.35 / 0.50 queued on H1 (the only panel whose NMAE behaves like the official one).

### 25.4 PDS and MSE lose in DISJOINT target populations (H1, 2026-10-03)

Answering the user's "PDS, MSE를 높이려면 어디를 건드려야 해".

**Structural fact, from the scorer source (not measured).** `streaming_bulk.py:176` + `prep.bulk_lognorm_means`:
`bulk_lognorm` is DERIVED from the counts accumulator and is `log1p(TS * S_g / sum_g S_g)` -- a function of the group
SUM alone, with no per-cell divide and no per-cell log1p. `metrics/delta.py:1078` states it outright: "fixed raw group
totals AUTOMATICALLY fix the pseudobulk". So PDS and MSE are a DETERMINISTIC function of the 300 group count totals;
there is no Jensen gap and no `cpc`-style functional mismatch on the sum side (this was the hypothesis being tested --
it is falsified from source). The within-column emission cannot move PDS or MSE at all.

**Decomposition** (exact offline evaluator, S16 structure, scorer space, panel genes excluded):
rho (per-target cosine) median 0.085 / mean 0.124, `||pred||/||real||` median 0.461, PDS 0.9098, MSE 0.704.

```
decile(real strength)  rho    |P|/|R|   %MSE numer   %MSE denom    PDS
 1 (strongest)        0.318    0.434       58.8%        52.7%     0.996
 2                    0.258    0.441       25.4%        22.3%     0.986
 3                    0.177    0.441       13.0%        11.6%     0.985
 4..7              0.168-0.042  ~0.5        10.8%         11.6%   0.939-0.986
 8                    0.025    0.358       -2.5%         0.5%     0.861
 9                    0.013    0.543       -2.6%         0.3%     0.738
10 (weakest)          0.008    0.419       -3.0%         0.1%     0.672
```
MSE is decided by the top 3 deciles (97 % of the numerator); PDS is decided by the bottom 3 (the top 4 are saturated at
0.99, 126/300 targets are already rank-1). **They do not overlap**, so no single lever raises both.

**The amplitude axis is CLOSED, with a number.** A greedy per-target amplitude ORACLE on top of `ta0.75` gives
mse 0.787 / pds 0.8314 -- WORSE than the shipped 0.704 / 0.9098 -- because zeroing 125/300 targets collapses
`B_pred` and with it the rebate (`r` 0.89 -> 0.287). Per-target amplitude is coupled to the global budget; tamp
beta 0.5-0.75 is the real optimum and there is no headroom left on this axis.

**Strength-gated hubness (`disw`), new and now closed.** The §24 rejection of `dis` was a GLOBAL application, which
costs MSE on the strong responders. Restricting it to the weakest targets by our OWN norm (test-time computable) and
making it norm-preserving (pure rotation, budget untouched) separates the two populations cleanly. H1 grid
(k 5/10/20 x lam 0.2/0.4/0.6 x frac 0.2/0.3/0.5/1.0): the mechanism is exactly as predicted --
at frac 1.0, dMSE +0.006..+0.034; at frac 0.2, dMSE +0.0000 for the SAME dPDS. Best cell k20 lam0.6 frac0.5:
dPDS +0.0024, dMSE +0.0008, J +0.0045 ~ **composite +0.0007**. Below the bar -> not adopted.

**Why it is small, and what it means.** The weak deciles' PDS loss is not hubness, it is rho: at rho 0.008-0.025 our
prediction is essentially orthogonal to the truth there, and no norm-preserving rotation creates signal. So BOTH ends
reduce to rho, and rho is source-bound in both populations:
  - MSE (strong responders, rho 0.32): H1GW same-line strong-target cosine is 0.49 vs our 0.34 -> same-line data hits
    exactly this decile. This is the MSE lever and it is a DATA lever.
  - PDS (weak responders, rho ~0.01): PDS is an equal-weight mean over targets, so 30 % of it is decided on targets
    whose real signal is ~1 % of the panel's noise-corrected total. Weak responders' transferable fraction is 0.12-0.23
    (§12) and ~88 % of their reproducible signal is trans (§20) -- there is signal there we do not capture, but it
    needs more source CONTEXTS, not more geometry.
Geometry operators on the sums (centring, pcr, dis/disw, tamp) are now exhausted: the last untested one is worth +0.0007.

### 25.5 The amplitude axis, re-opened and then closed properly (2026-10-03)

**§25.4's "amplitude is closed" was reached with the WRONG oracle and is retracted.** The greedy per-target MSE oracle
zeroes 125/300 targets, which collapses `B_pred` and with it the rebate (`r` 0.89 -> 0.287); it looked negative for that
reason, not because the allocation was optimal. The correct question holds the budget up. Re-run on H1 (exact offline
evaluator, S16 structure): `s_t ~ rho^0.5` with a global 0.7 gives **dPDS +0.0085, dMSE -0.0678, J +0.0901 ~ composite
+0.015**. So there WAS headroom on this axis.

**It decomposes into two nearly additive levers** (H1, offline):
```
(a) rho-hat MAGNITUDE reshaping, budget matched back   dPDS +0.0023  dMSE -0.0065  J +0.0120   budget-INVARIANT
(b) budget level B/claim 0.89 -> 0.977 (uniform 1.048) dPDS -0.0005  dMSE -0.0144  J +0.0140   budget-CHANGING
(c) both (= token `tz0.25x0.9`)                        dPDS +0.0024  dMSE -0.0211  J +0.0275
```
Mechanism of (b): the MSE rebate is `r = min(1, B_pred/sum claim)` and the shipped pipeline's budget match restores the
BASE recipe's `B_pred`, which leaves `r` at 0.89 -- 11 % of the rebate unclaimed. (a) is new information: `ta` uses
rank(rho-hat)^beta, and the rank transform discards the magnitude spread that (a) exploits.

**REAL SCORER, H1, token `tz0.25x0.9` (= c):**
```
                 offline   REAL      base    arm
pds_cosine        +0.0024  +0.0024   0.9095  0.9119   <- transferred EXACTLY
expr_mse_...      -0.0211  +0.0055   0.7040  0.7095   <- SIGN FLIP
composite         +0.0046  +0.0001   0.3169  0.3170
```
**This is §24.6's rule reproduced on a single panel, and it sharpens it.** Offline PDS is exact; offline MSE is
*qualitatively* wrong whenever the arm moves `B_pred`, because `sch_sumeval.evaluate` models `claim` with the REAL
panel's jackknife as a stand-in for `C_pred_hat`, while the scorer uses the emitted cells' own `C_pred_hat` -- and the
#247 cap binds on every row (§25.1), so the modelled `B/claim` ratio is fictitious. **Rule: price PDS offline, never
price MSE offline for any arm that changes the sum budget.**

=> (b) is dead on the real scorer, on H1 itself, without needing the OR_X confirmation that killed the §24.6 kink
factor. The only survivor is (a), whose PDS gain transfers exactly (+0.0024 raw ~ +0.0004 composite). Ported as the
budget-invariant token `tzb<gamma>` (applied BEFORE the budget match, so `bm` restores `B_pred`); real-scorer runs on
H1 / OR_HCT116X / RN_Jurkat in unit aivc-tzb-20261003. Both tokens default off; `tz`/`tzb` absent leaves every existing
arm bit-identical (parse verified).

### 25.6 `tzb` three-panel result: PDS transfers, the composite does NOT (2026-10-03)

Real scorer, S16 structure, token `tzb0.25` (rho-hat MAGNITUDE reshaping of the per-target sums, applied BEFORE the
budget match so `bm` restores `B_pred` -- budget-INVARIANT):
```
panel            dPDS     dMSE    dNMAE   dreach   dcomposite
H1            +0.0023  -0.0065  -0.0001  +0.0010     +0.0022
OR_HCT116X    -0.0000  -0.0379  +0.0373  -0.0024     -0.0108
RN_Jurkat     +0.0016  +0.0126  -0.0008  +0.0012     -0.0013
```
R21 composite 0.5*(0.75 OR + 0.25 RN) = **-0.0042**: FAIL. 1/3 panels positive on the composite.

What IS consistent: **PDS is positive or flat on 3/3** (+0.0023 / -0.0000 / +0.0016), and the H1 offline prediction
(dPDS +0.0023, dMSE -0.0065) reproduced the real scorer to four decimals on BOTH members -- so for a budget-INVARIANT
sum-side arm the offline evaluator is exact. (Contrast the budget-CHANGING `tz0.25x0.9` in §25.5, where offline said
dMSE -0.0211 and the scorer said +0.0055.)

What is NOT consistent: MSE (2/3 better; RN_Jurkat +0.0126) and NMAE (only OR_X breaks, +0.0373). The OR_X NMAE break
is the §23.11 mechanism -- the DE side keeps the PRE-shrink sums while the emission pins the POST-shrink sums, so the
`cpc` composition pin is unreachable for shrunk targets. The shrink itself is mild (1/263 targets below 0.5x; the
budget match applies a broad ~14 % reduction), so this is a sum/DE MISMATCH, not crushed targets. A DE-following
option `tzb<gamma>x<a>` (DE side multiplied by w**a) was added to test the repair; a=0 is the old behaviour and every
existing arm parses bit-identically.

**Ceiling of this lever, budget held (H1, offline -- trustworthy here):**
```
predictor                gamma    dPDS     dMSE      ~composite
tlrn rho-hat (Sp 0.626)   0.6   +0.0038  -0.0118      +0.0035
ORACLE true rho           1.0   +0.0063  -0.0662      +0.0139   <- ceiling
```
We capture ~25 % of it. Unlike the six failed source-weighting learners, the object to predict here is ONE SCALAR per
target (how well our own delta will do), the training set is 8 harness panels, and the offline evaluator is proven
exact -- so iteration is cheap. Unused uncertainty features we already compute: disagreement among the 7 `min_cells`
ensemble members, the DepMap-neighbour-vs-ours disagreement, source-cell bootstrap stability.

### 25.7 `tzb` CLOSED; the DE repair failed and the mechanism is structural (2026-10-03)

The OR_X NMAE break of §25.6 was diagnosed as a sum/DE mismatch (the DE side keeps the PRE-shrink sums while the
emission pins the POST-shrink ones) and a DE-following option `tzb<gamma>x<a>` was added to repair it. **It made it
worse:**
```
OR_HCT116X                base     tzb0.25   tzb0.25x1.0
  pds_cosine            0.7328    -0.0000      -0.0000
  expr_mse              0.0299    -0.0379      -0.0379
  de_wilcoxon_lfc_nmae  0.8628    +0.0373      +0.0530    <- repair made it WORSE
  composite             0.3942    -0.0108      -0.0155
```
So the diagnosis was wrong. NMAE is a MAGNITUDE metric: any per-target shrink of the predicted LFCs raises it
directly, and making the DE side follow the shrink explicitly just doubles the leak. This is structural, not a
plumbing bug -- the amplitude lever trades MSE against NMAE and there is no version that keeps only the MSE half.

**`tzb` is closed**: R21 composite -0.0042, 1/3 panels positive, repair failed. Both tokens (`tz`, `tzb`) stay in
`sch_emit.py` default-off; every existing arm parses bit-identically (verified).

**What survives and is NOT conditional on `tzb`:** the shipped `ta0.75` consumes the SAME rho-hat, so a better
predictor improves the CURRENT pipeline in its own rank functional form -- that is exactly how S16 was made from S15.
Measured gradient on H1 (budget-invariant, offline exact): ~+0.0028 composite per +0.1 of rho-hat Spearman, with the
oracle at +0.0139. Predictor v2 (already trained, Spearman 0.669 vs v1's 0.626) is worth +0.0048 vs +0.0035 in
magnitude form -- i.e. a better rho-hat already exists and is unused.

**v3 features** (`scripts/ml_amp_learn3.py`): the SPREAD / min / max of the pairwise family agreement (v1-v2 use only
its mean), a family-WEIGHTED agreement, a family SPLIT-HALF consensus cosine, the S4-vs-fusion norm ratio and the
measured-gene count. Single-feature Spearman with the true per-target cosine:
```
            H1     RN_RPE1  RN_K562  OR_X   OR_HCT116
ag_w       0.662    0.693    0.611   0.396    0.478
splithalf  0.659    0.590    0.518   0.163    0.160
ag_max     0.587    0.676    0.431   0.391    0.501
```
`ag_w` alone beats the whole v1 predictor (0.626) on H1. Strong on H1/RN, weak on the shallow OR panels (163 cells).

### 25.8 J is BLIND to NMAE, and that is the amplitude lever's main cost (2026-10-03)

`J = 2.23*dPDS - 1.05*dMSE ~ 6 x composite` was derived for sum-side arms that leave the DE side alone. A per-target
AMPLITUDE arm does not: the emission's `cpc` pin couples the sums to the DE LFCs, so the arm moves NMAE, which J does
not contain. Measured, OR_HCT116X, `tzb0.25`, real scorer:
```
J  = 2.23*(-0.0000) - 1.05*(-0.0379) = +0.0398      <- "strongly positive"
real composite delta                  = -0.0108      <- actually negative
the whole gap is NMAE                 = +0.0373
```
On H1 the same arm has NMAE -0.0001, so there J/6 = +0.0020 against a real +0.0022 -- accurate. **J is accurate
exactly where NMAE does not move.** Every offline amplitude screen in this project (including the ones behind `tamp`
/ S15) used J; S15 survived because the official NMAE behaves like H1's, not because J was right.
**Rule: never screen a per-target amplitude arm on J. Screen PDS and MSE offline (exact when budget-invariant,
§25.6), then pay for a real-scorer run to see NMAE.**

### 25.9 Amplitude predictor v3 features: NULL (2026-10-03)

Leave-one-panel-out over all 8 harness panels, v2 (12 features) vs v3 (19):
```
mean Spearman with the true per-target cosine:  v2 +0.6722   v3 +0.6736   (+0.0014, 5/8 panels)
per panel  H1 .697->.700  KOLF .376->.386  RN_Jurkat .783->.800  RN_HepG2 .802->.801
           RN_RPE1 .783->.787  RN_K562 .736->.722  OR_X .556->.576  OR_HCT116 .644->.616
```
The new features have strong MARGINAL Spearman (ag_w 0.662, splithalf 0.659 on H1 -- `ag_w` alone beats the whole v1
predictor) but are REDUNDANT with the existing mean-agreement column, so the GBM gains nothing. The agreement family
is saturated as a rho-hat signal; a better predictor needs an independent signal (the `min_cells` ensemble members'
disagreement is the obvious untried one -- it was not built here because the members are not cached).

**Amplitude axis, final accounting.** (1) budget half: unreachable (offline MSE sign-flips, real scorer negative,
§25.5); (2) shape half: real and PDS-positive 3/3, but costs NMAE, R21 -0.0042, repair failed (§25.6-25.7);
(3) better rho-hat: v3 null, v2 marginal. The axis is CLOSED -- this time with the right oracle, the right
evaluator and a real-scorer measurement, unlike the §25.4 version that was retracted.

### 25.10 Drop-in rho-hat replacement in the SHIPPED form: ~0 (2026-10-03) -- amplitude axis CLOSED

`ta<beta>` consumes rho-hat as `rank(rho-hat)**beta`, so a better predictor should improve the CURRENT pipeline
without any new token. Leave-one-panel-out over 8 panels, J vs the shipped S16 amplitude (= v1 at rank^0.75):
```
  v2 rank^0.5   J -0.0144 (~composite -0.0024)  3/8        v3 rank^0.5   J -0.0142 (-0.0024)  3/8
  v2 rank^0.75  J +0.0016 (~composite +0.0003)  5/8        v3 rank^0.75  J -0.0003 (-0.0000)  5/8
  v2 rank^1.0   J +0.0061 (~composite +0.0010)  4/8        v3 rank^1.0   J +0.0063 (+0.0010)  5/8
```
Replacing v1 with v2 or v3 at the shipped beta is worth +0.0003 / -0.0000, and J OVERSTATES (it is NMAE-blind,
§25.8). The rank transform compresses the predictors' differences away: +0.02 of Spearman does not reorder enough
targets to matter.

**AMPLITUDE AXIS, CLOSED (final).**
```
  budget half (B/claim 0.89 -> ~1.0)     unreachable: offline MSE sign-flips, real scorer +0.0055 on H1
  magnitude/shape half                   real, PDS +3/3, but NMAE-costly; R21 -0.0042; DE repair made it worse
  better rho-hat in the shipped form     +0.0003 (v2) / -0.0000 (v3), and that is an overstatement
  v3 feature engineering                 null (mean Spearman +0.6722 -> +0.6736)
```
The §25.4 oracle of "+0.0139 composite" was in J units on H1 and is therefore inflated twice over: J is NMAE-blind,
and it is an oracle. Net deliverable gain from this whole line: **zero**. Net deliverable KNOWLEDGE: §25.5 (offline
exact iff budget-invariant), §25.8 (J is NMAE-blind), §25.4 (PDS and MSE lose in disjoint target populations).

**Offline exactness confirmed at a SECOND exponent (H1, `tzb0.6`):** offline dPDS +0.0038 / dMSE -0.0118, real
scorer dPDS +0.0038 / dMSE -0.0118 (NMAE -0.0002, composite +0.0033). Two independent settings now reproduce to four
decimals, so §25.5's rule is established, not a coincidence. H1 alone would have adopted this arm (+0.0033); the
3-panel R21 gate is what stops it, and that is the gate working as designed.

### 25.11 Neighbour-cut ensemble (`nbe`): the axis has no member diversity (2026-10-03)

`sch_emit.py:110` carried the note "NB_TOPK = 10 -- another discontinuous SELECTION boundary, ensemble-able like
min_cells", never built. Built now (token `nbe<k1>.<k2>...`, averages the DepMap-neighbour estimate over several cuts;
`NB_TOPK_ENS` empty => the shipped single cut, bit-identical).

Real scorer, H1, `nbe5.10.20.40`: composite +0.0002 (PDS -0.0000, MSE -0.0007, NMAE -0.0005). Nil.

**The note was wrong, and the measurement says why.** Fraction of the c^2 neighbour weight captured by the top-k cut:
```
   top-3  0.455    top-5  0.645    top-10  0.909 (MEDIAN 1.000)    top-20 / 40 / 80  1.000
```
The neighbour lists are at most ~20 long, so members 20/40/80 are IDENTICAL; for the median target the list is <= 10
and the shipped cut already takes everything; and the c^2 weighting leaves the marginal neighbours almost no mass.
Member diversity is ~0, and by our own measured principle (the ensemble gain tracks member DIVERSITY, not member
quality) the gain must be ~0 -- which is what the scorer says. This CONFIRMS the principle rather than contradicting
it. The only cuts with real diversity (3, 5) are strictly worse estimates (fewer neighbours), not peer alternatives
the way different `min_cells` values are; `nbe3.5.10.20` is being measured for completeness.

Diversity-axis status: `min_cells` (shipped, the whole of S11); `pcr-lambda` / centring averages (provable no-ops);
row bootstrap (INFEASIBLE -- the dominant source replogle_k562_gwps has exactly 1 row per target, 0 % with >1);
neighbour cut (no diversity). The cheap ensemble axes are exhausted.

### 25.12 Round summary: what "improving PDS and MSE" is now blocked on (2026-10-03)

Asked whether PDS/MSE improvement is finished. It is NOT -- S15 moved MSE 0.749 -> 0.696 for +0.0106 official two days
ago, and S10i / S10k / S11 before it were all PDS/MSE changes. What IS finished is the family of operators that
RESHAPE a delta we already have. Everything tested this round came from that family and all of it closed:
```
  per-target amplitude, budget half        unreachable (offline MSE sign-flips; real scorer +0.0055 on H1)
  per-target amplitude, shape half         PDS +3/3 but NMAE-costly; R21 -0.0042; DE repair made it worse
  better rho-hat (v2 / v3 predictors)      +0.0003 / -0.0000 in the shipped form; v3 features null
  strength-gated hubness (`disw`)          +0.0007
  neighbour-cut ensemble (`nbe`)           +0.0002; the axis has no member diversity (lists are <= 20 long)
  continuous DE magnitude calibration (`drp`) -0.0013
```
Ensemble diversity axes are now exhausted: `min_cells` is shipped and was the whole of S11; `pcr-lambda` and centring
averages are provable no-ops; the neighbour cut has no diversity; and a ROW bootstrap is infeasible because the
dominant source `replogle_k562_gwps` has exactly 1 row per target (0 % with >1). A CELL bootstrap would need the
single-cell sources rebuilt and would still leave K562 (fusion weight 8) fixed.

What remains, with honest priors:
  (1) change how the delta is BUILT (cell-level model). 17 architectures already lose 0.04-0.22 via PDS; Codex R20
      judged it does not address PDS. Expensive, low prior.
  (2) more source CONTEXTS (+0.003 each, curve still rising at k=6) and same-line data (+0.02..+0.04 per comparable
      destination). Opens on 10/22 when D/E/F are named.
The binding constraint is unchanged; this round bought no score and three reusable methodological results
(§25.5 offline exactness condition, §25.8 J is NMAE-blind, §25.4 PDS and MSE lose in disjoint target populations).

**`nbe` final (H1):** `nbe5.10.20.40` +0.0002, `nbe3.5.10.20` +0.0005 (PDS -0.0000 / -0.0001, MSE -0.0007 / -0.0018).
Even the setting with real member diversity (top-3/5 differ from top-10 by 0.45 / 0.26 of the c^2 weight) is nil, so
the neighbour-cut axis is closed on measurement as well as on mechanism.

### 25.13 Pseudocount ensemble (`cens`): the one POSITIVE result of the round -- 3/3 panels, R21 +0.0024 (2026-10-03)

User asked whether the TRANSFER itself can be improved. The combination side is saturated (fixed/precision/MODZ/robust
/z-shrink/per-source-amplitude/rank-space all negative; 6 learned source-weighting failures; oracle only +0.017). The
ESTIMATOR side had one untried axis: the LFC pseudocount `c`, fixed at 50 for the fusion part. `_C3` (c in 10/20/100)
was already DEFINED in `ENS["mcwc"]` and had **never been scored** (0 arms).

Built as `delta_s4(..., c_list=)`, token `cens<c1>.<c2>...`, unioned with the `min_cells` members exactly as
`ENS = _MC5 + _W6 + _C3` does. Valid c values are only those with precomputed centring vectors,
`tv2_sources.CS = (1, 5, 10, 20, 50, 100)` -- a parse-time guard now enforces this (c=200 produced `KeyError
center_c200`). `cens` absent => `c_list=None` => the loop and its cache keys are bit-identical to before.

Real scorer, `cens5.10.20.100`:
```
panel          dPDS     dMSE    dNMAE   dreach   dcomposite
H1           +0.0009  +0.0065  -0.0020  +0.0042    +0.0003
OR_HCT116X   +0.0003  +0.0208  -0.0177  +0.0037    +0.0042
RN_Jurkat    +0.0026  -0.0138  -0.0076  +0.0033    +0.0063
R21 = 0.5*(0.75 OR + 0.25 RN) = +0.0024     3/3 panels positive
```
**PDS, NMAE and reach are positive on 3/3** (reach remarkably stable at +0.0033..+0.0042); only MSE is panel-dependent.
R21 +0.0024 is BELOW the +0.005 bar, so this is a RIDE-ALONG candidate, not a submission on its own -- the same class
as `drc 0.8/1.1` (H1 +0.0018, 3/3).

**A premise I used to motivate it was wrong**, and the data says so: I cited "the c sweep spans only 0.0071 pds" as
evidence the members are of EQUAL quality, but that figure is from the OLD parameterisation (c = 0.1..3.0, §4b), not
the current CPM pseudocount 1-100. The MSE degradation on H1/OR_X shows the low-c members are genuinely worse on the
sum side. So this axis is "high diversity + UNEQUAL quality", unlike `min_cells` -- which is why it pays on the DE
members (where low c sharpens the calls on the 84-88 % of genes below CP10K 1) and costs on MSE.

A DE-SIDE-ONLY variant (sums pinned via the `dsrc` pattern, so PDS/MSE are invariant by construction) was estimated
from these members at R21 ~ +0.0033 -- still below the bar, and it would give up the 3/3 PDS gain, so it was NOT built.

**`cens` width swept (H1):** `cens5.10.20.100` +0.0003, `cens1.5.10.20.100` **+0.0000** (MSE +0.0100 vs +0.0065).
Adding c=1 costs more on the sum side than its DE-side sharpening returns, so the adopted ride-along setting is
`cens5.10.20.100`.

### 25.14 Ride-alongs PORTED to the official builder and audited (2026-10-03)

`drc` and `cens` existed only in `sch_emit.py`, so the "ride-along list" could not actually be carried on 10/22.
Both are now in the official path, default off.

- **`cens` needed no builder change.** `parts[].cfg_ens` already takes arbitrary cfg overrides, so the pseudocount
  members join the SAME union as `min_cells` -- at `sum.spec_override.parts[1].cfg_ens`, which is where
  `tv2_final.py --sum-ens` writes. New flag `--sum-cens 5,10,20,100`, with a parse-time guard (only
  `tv2_sources.CS = 1/5/10/20/50/100` have precomputed `center_c<c>`; c=200 raised `KeyError: center_c200`).
- **`drc` is a new builder block** (`de.drc`), placed after the DE-side `pcr` step, plus `--de-drc 0.8,1.1`.
  `lfc_de` is now passed when `drc` is set as well.

**Audit on spec_S16, `audit_only`, contexts A/B/C:**
```
drc formula builder vs harness      max|diff| 0.000e+00, bit-identical
genes modified per target           13,882 / 18,533  (= 5 % x 0.8 + 70 % x 1.1)
override columns                    untouched
pred_dsum_{A,B,C} with/without drc  BIT-IDENTICAL  -> PDS / MSE exactly invariant
cens, 11 members (7 min_cells + 4 c)  builds on A/B/C; sums cosine 0.9938 vs base;
                                      gamma_sum 1.399 -> 1.275..1.283 (budget match absorbing the
                                      ensemble's reduced across-target spread, as documented for S11)
```
Note: the `_dde_` dump is written BEFORE the DE-side transforms (true for `pcr` too), so drc is not visible in it;
that pre-existing dump order was left alone and the port was verified against the harness formula directly.

### 25.15 Final-phase robustness audit (2026-10-04, user "다른 방법은 없어?" -> brainstorming spike)

The question that had never been asked: the shipped hyperparameters were tuned on A/B/C + the harness, but the PRIZE
is the final phase on three UNKNOWN lines. Rehearsals 1-9 all verified plumbing, never hyperparameter transfer. And we
already had one knob with a measured SIGN FLIP by destination (KOLF on the sum side: H1 +0.0079 PDS, OR_HCT116 -0.0126).

**Method (free): re-analyse all 515 already-scored (arm, panel) results.** Arms are token lists, so any two arms
differing in exactly one token give a measured per-panel delta for that knob. The first cut compared panels on
DIFFERENT transitions and produced false alarms; grouping by the EXACT transition fixes it. Validation of the method:
it independently reproduces FINDINGS 24.6 (`gk1.35 -> gk` = +0.1813 on OR_HCT116X, i.e. the -0.181 recorded there).

**Verdict: the shipped knobs are at robust values.**
```
SAFE   K562-denoised DE side  7/7 panels positive (+0.0005..+0.0073)
SAFE   --s4-excl KOLF,FENG    5 panels: 4 favour it, 0 against
SAFE   agc                    inert, |d| <= 0.0006 on 7 panels
SAFE   ct0.5                  the hedge in a real trade (ct1.0 favours stem, absent favours essential panels)
SAFE   k3000                  beats k5000 on 2/3
SENSITIVE  gk 1.17 / mg 0.8   the kink; +0.18 swing on OR_X for a 15 % factor change. Already conservative.
```

**A false alarm I raised and then killed with measurement.** I read "RN_Jurkat +0.709 vs official context A -1.33 for
the SAME cell line" as a 2.04 protocol effect against a 0.755 threshold margin. It is not: the RN panels carry a
reduced ~8k gene axis that keeps only 2 of the 12 stem markers (HepG2 2/12, K562 3/12, RPE1 1/12), so their scores are
averages over a handful of genes. On the six control sets with all 24 markers (H1, KOLF, OR_HCT116, official A/B/C)
the stem-vs-somatic gap is **3.585**, and +1.5 sits 0.755 below the lowest stem and 2.83 above the highest somatic.
The threshold's bias toward "somatic" is correct risk management, because the error costs are asymmetric: calling a
somatic line stem costs -0.0126 PDS, calling a stem line somatic only forgoes +0.0079. A scale-free RANK version of
the score was tested as a fix and FAILS (gap -0.10; RN_Jurkat outranks both stem panels). **No change needed** --
except a guard: `stem_score` now asserts >= 10/12 markers on each panel, so a reduced axis stops the run instead of
silently returning a wrong number (verified: official axis -1.330 passes, truncated axis refuses).

**A real bug found while setting up the next test.** `wmul` (per-family source weights, harness token `wm<FAM>x<m>`)
was wired into `delta_fs` only, NOT into `delta_s4` -- the path every S16-structure arm takes. Any family-weight arm
run on the shipped structure would have been a silent NO-OP, and "family weights do not matter" would have been the
wrong conclusion drawn from it. Now forwarded, with a `_wm<FAM><m>` cache-key suffix so weighted builds cannot collide
with plain ones; `wmul=None` leaves every existing cache key unchanged.

**Rehearsal 10 (2026-10-04) PASS.** The current command including `--sum-tamp 0.5 --sum-cens 5,10,20,100
--de-drc 0.8,1.1` ran end to end on the A/B/C controls without `--build`: `de.drc [0.8, 1.1]` in the spec, `cfg_ens`
11 members (7 min_cells + 4 pseudocounts) in all three contexts, stem scores -1.33 / -5.108 / -2.434 reproducing the
documented values, all three routed somatic with the no-FENG-no-KOLF S4, no broken symlinks.

### 25.16 Source family weights ARE destination-dependent, and the predictor is test-time computable (2026-10-04)

The robustness audit's one real finding. Family weights had NEVER been varied as a standalone cross-panel arm (zero
`wm` arms existed), and the reason turned out to be a bug: `wmul` was wired into `delta_fs` only, not `delta_s4`, so
every S16-structure arm would have ignored it silently. Fixed (cache key suffixed `_wm<FAM><m>`; `wmul=None` leaves
all existing keys unchanged), then measured on the real scorer:
```
panel        lineage                   base     K562 x0.5   K562 x2   prefers
H1           stem                     0.3169     +0.0021    -0.0025   LESS K562
OR_HCT116X   colorectal epithelial    0.3942     -0.0006    -0.0027   LESS K562
RN_Jurkat    lymphoid (haematopoietic) 0.4892    -0.0061    +0.0042   MORE K562
```
K562 is erythroleukemia, i.e. haematopoietic like Jurkat. The sign FLIPS with lineage match and the span on one panel
is 0.0103, so a factor-2 error in the wrong direction costs ~0.006 on a lineage-matched destination -- and the final
phase has three lines. (RN_Jurkat excludes JURKAT/RNJURKAT but NOT K562, so this arm is valid; RN_K562 excludes K562
and cannot be used for it.)

**The predictor exists and needs only controls.** Source control profile = mean over rows of `z["fc"]` (control
composition) as log1p(CPM); destination profile = its own control pseudobulk; Spearman on shared genes:
```
panel         prefers      K562    H1   HCT116  HEK293T  KOLF   CD4   FENG   nearest source
H1 (stem)     LESS K562    0.451  0.999  0.768   0.789  0.868  0.726  0.903   H1
OR_X (colorec) LESS K562   0.669  0.754  0.969   0.850  0.788  0.720  0.784   HCT116  <- correct
RN_Jurkat     MORE K562    0.805  0.453  0.683   0.646  0.559  0.527  0.699   K562    <- correct
```
The K562 column is MONOTONE in the measured preference on 3/3 panels, and the nearest-source identification is right
for every panel. This is the stem score (binary, worth +0.0061 officially) generalised to a CONTINUOUS per-family
lineage weight, computable from D/E/F's controls alone.

**Falsification test running** (unit aivc-wmul2-20261004): CD4's similarity is LOWEST on RN_Jurkat (0.527) and highest
on H1 (0.726), the opposite of the lineage intuition (Jurkat and CD4 are both T-lineage) because primary cells differ
from immortalised lines more than lineage separates them. So the rule predicts **Jurkat prefers CD4 x0.5 while H1 and
OR_X prefer CD4 x2**. If that comes out, the signal is similarity, not a lucky lineage story; if Jurkat prefers more
CD4, the rule is wrong and only the K562 result stands as a one-off.

### 25.17 Preprocessing audit (2026-10-05, user "우리가 정말 전처리를 잘 하고있는게 맞아?")

Preprocessing has produced two real bugs before (the gene axis cut from 17,847 to 5,592; the scrambled
`controls_common` gene order), so it was audited rather than assumed. The decisive check is ON-TARGET KNOCKDOWN: if
rows, gene axis and control matching are all correct, knocking down gene X must drive X's own LFC strongly negative.

```
source          self-KD median   KD < -0.3   genes meas   NaN%   |center|/|row|
H1                  -0.991          96%        18077      2.5        0.38
CD4_Stim8hr         -0.573          77%        18056      2.6        0.22
FENG_GWN            -0.473          71%        18513      0.1        0.31
K562                -0.438          71%         7679     58.6        0.12
HCT116              -0.390          63%        18407      0.7        0.07
KOLF                -0.299          50%        18270      1.4        0.20
HEK293T             -0.284          46%        18407      0.7        0.13
JURKAT_STIM         -0.064           5%        18491      0.2        0.25   <- flagged by the audit itself
JURKAT_UNSTIM       -0.098          14%        18491      0.2        0.10   <- (GSE249595, already killed twice,
                                                                                 and FAMILY weight is 0)
```
**No alignment or control-matching bug.** Every shipped family shows the knockdown; K562's 58.6 % NaN is Replogle
measuring 7,679 of 18,533 genes (known and correct); the audit independently re-found the one bad source we had
already killed, which is the check validating itself. `ChangYe2021` is degenerate (1 row) but is not in `FAMILY`.

**One genuine inconsistency found.** `center_c<c>` is an arithmetic MEAN over the source's own rows
(`tv2_sources.py:79`), so what it removes depends on the screen:
```
  HCT116 0.07, K562 0.12, HEK293T 0.13, DLD1 0.11   (9.8k-row genome-wide screens: a technical baseline)
  KOLF 0.20, CD4 0.22, FENG 0.27-0.31
  H1 0.38                                           (its screen IS the 300 official targets, selected to respond)
  TianKampmann 0.50, WesselsSatija 0.58, PapalexiSatija 0.70   (25-180 rows: the mean is unstable)
```
A mean over 25 rows can be dragged by a couple of huge-effect targets, and for H1 the "genome-wide mean" is really
the mean RESPONSE of the panel we are predicting. A ROBUST (trimmed-mean / median) center has never been tested.
Implemented as `Cfg.center_robust` (trim fraction per tail; >= 0.5 = median), computed on the fly from the source's
own `fp`/`fc` and cached; harness token `crob<q>`; `center_robust=0` is bit-identical to the shipped center
(verified, max|diff| 0.000e+00). At q=0.5 it changes 99.9 % of entries with max|diff| 0.089 and energy ratio 1.009.
Real-scorer test queued on H1 / OR_HCT116X / RN_Jurkat (unit aivc-crob-20261005).

Also added this round and parity-checked: `Cfg.fwins` (family-level winsorize across the families that measured each
(t, g); H1 +0.0015, OR_X -0.0008) and `wmul` finally wired into `delta_s4` (it had been `delta_fs`-only, so every
family-weight arm on the shipped structure would have been a silent no-op).

### 25.18 Preprocessing: audited clean, and the whole per-gene family CLOSED by two oracles (2026-10-05)

User: "전처리도 다양하게 시도해보자 ... 어떤 사람은 bin을 이용해서 처리하고". Preprocessing was decomposed and each
branch priced rather than tried one at a time.

**(1) Health audit -- clean.** On-target knockdown (the check that catches row misalignment, a wrong gene axis and
bad control matching at once) is present in every shipped family: self-KD median -0.284 (HEK293T) to -0.991 (H1),
46-96 % of targets below -0.3. The audit independently re-flagged JURKAT_STIM/UNSTIM (-0.064 / -0.098, 5 % / 14 %),
which is GSE249595 -- already killed twice and carrying FAMILY weight 0. K562's 58.6 % NaN is Replogle measuring
7,679 of 18,533 genes. `ChangYe2021` is degenerate but is not in `FAMILY`.

**(2) BINNING by expression -- closed at the ORACLE.** 8 bins by the destination's own control CPM, each given the
scalar that best matches the truth (budget held, so the evaluator is exact):
```
  oracle per-bin scalars [-0.382 -0.021 0.623 0.641 0.510 0.448 0.415 0.422]
  -> dPDS -0.0008  dMSE -0.0027  J +0.0011  ~composite +0.0002
```
The scalars for the bins that carry weight are nearly uniform (0.41-0.64), so the budget match absorbs them; the only
real information is that the two lowest-expression bins should be zeroed, and they contribute nothing.

**(3) Per-GENE correction -- big oracle, ZERO transfer.** Giving every gene its own oracle scalar (18,533 dof) is
worth **+0.0151** composite on H1 (sign alone +0.0043; 19 % of genes have a NEGATIVE scalar, i.e. we predict their
direction backwards). But learned on two panels and applied to the held-out third it is **negative on 3/3**:
```
  held out H1          scalar ~comp -0.0216   sign-only -0.0129
  held out OR_HCT116X                -0.0758             -0.0075
  held out RN_Jurkat                 -0.0478             -0.0137
  (cross-panel Spearman of the per-gene scalar: 0.05-0.16)
```
So that +0.0151 is destination-specific overfitting, not a systematic source bias.

**(2)+(3) together CLOSE the whole family.** Any per-gene or per-expression-level preprocessing change to the sums --
a variance-stabilising transform instead of log, a binned pseudocount, expression weighting (`ew`), per-gene
calibration -- is a special case of one of these two, and is therefore bounded by +0.0002 (if it is a function of
expression) or dead on transfer (if it is gene-specific). **They do not need to be tried individually.**

**(4) Source-level knobs, which those bounds do NOT cover:**
```
  robust (trimmed / median) per-source center `crob`   H1 +0.0003  OR_X -0.0017  RN_Jurkat +0.0029   ~nil
  family-level winsorize `fwins0.1`                    H1 +0.0015  OR_X -0.0008                      marginal
```
Motivation for `crob` was real -- `center_c<c>` is an arithmetic MEAN and |center|/|row| runs 0.07 (HCT116, 9.8k rows)
to 0.70 (PapalexiSatija, 25 rows), and for H1 the "genome-wide mean" is the mean RESPONSE of the panel we predict --
but the measurement says it does not matter.

**(5) Cell binning / metacells: closed by implication.** The pseudobulk is a SUM, and the sum of metacells is the
same sum; only a robust aggregate over them differs, and median / winsor / trim already lose on 5/5 lines.
Token-style expression binning (scGPT / Geneformer) only serves a transformer, and every joint/conditional
architecture loses 0.04-0.22 composite via PDS.

**Conclusion: preprocessing is not where the problem is.** It is audited clean and its remaining degrees of freedom
are bounded at +0.0002 (expression-level) or do not transfer (gene-level). New code is default-off and parity-checked:
`Cfg.center_robust` (+ token `crob<q>`), `Cfg.fwins` (+ token `fw<q>`), and `wmul` finally reaching `delta_s4`.

### 25.19 THE HARNESS AND THE SHIPPED SPEC USE DIFFERENT SOURCE SETS (2026-10-06)

Found while answering "가장 개선 시도 할곳은". Five built families with real official-panel coverage carry NO weight
in `tv2_final.BASE_NEW`:
```
family      rows   of the official 300   median cells/target   self-KD median
RPE1         836          95                   101                -0.659
RNJURKAT     846          94                    91                -0.635
DLD1        4714          85                   137                -0.407
HEPG2        783          77                    50                -0.616
H1GW        2747         131                     -                -0.255   (H1 same-line; stem destinations only)
```
Their on-target knockdown is **better than every shipped family except H1** (HEK293T -0.284, KOLF -0.299,
HCT116 -0.390, K562 -0.438). The handoff's reason for excluding them -- "RN lines: essential-gene panels
(**0/300 official targets**)" (gap analysis 2026-09-23) -- is **wrong**: they cover 77-95 of the 300, verified by
name (ARID3A, ARPC2, ATF4, ATP6V0B, CAST, CHMP3, COASY, CYCS, CYFIP1, DAXX, ...).

**And the harness already includes three of them.** Part-1 weights:
```
  harness:  CD4 1, FENG 0.25, H1 2, HCT116 1, HEK293T 0.5, K562 8, KOLF 0.25, HEPG2 1, RNJURKAT 1, RPE1 1
  S16 spec: CD4 1, FENG 0.25, H1 2, HCT116 1, HEK293T 0.5, K562 8, KOLF 0.25
```
So every harness measurement has been made on a source set the submissions do not have, and **the harness is
structurally blind to this difference** (it cannot measure adding what it already contains) -- which is why two weeks
of harness work never surfaced it. It is also a candidate contributor to the long-standing pattern that harness gains
over-predict official gains.

**Measurement in flight** (unit aivc-rnout-20261006): remove HEPG2/RNJURKAT/RPE1 from the harness on H1 and
OR_HCT116X. The loss from removing them is what the official build is giving up by leaving `--rn-sources` off; a
weight x2 arm asks whether 1.0 is the right level. `tv2_final.py --rn-sources [--rn-weight W]` already exists
(added 2026-09-27, deliberately left off on the strength of the wrong 0/300 claim).

**25.19 MEASURED.** Removing HEPG2/RNJURKAT/RPE1 from the harness (i.e. scoring the SHIPPED source set):
```
  H1          -0.0023   (PDS -0.0003, MSE +0.0090)
  OR_HCT116X  -0.0013   (PDS -0.0005, MSE +0.0228)
  H1, RN weight x2  -0.0008  -> weight 1.0 is the right level, not 2.0
```
So the shipped build gives up **+0.0013..+0.0023, 2/2 panels**, by leaving `--rn-sources` off. Less than the
+0.003/context pricing because coverage is 26-32 %, not full. Note `TM.SRC` is REDIRECTED by `fh_harness` to the
harness source dir, which is why the RN files appear present under the harness and absent under the official path;
`--rn-sources` builds them for the official run via `tv2_sources.py RN_HepG2 RN_Jurkat RN_RPE1`. The harness copies
were built by `fh_sources.py`, so the two builders' RN outputs must be compared before trusting the transfer.
