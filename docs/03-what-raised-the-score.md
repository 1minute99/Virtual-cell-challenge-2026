# 03 — What raised the score

This file lists every change that raised the **official** validation score, in order of impact, with what the change
was, the official before → after, whether the pair isolates it, and the local evidence behind it. Changes that did not
pay are in `04-what-did-not-work.md`; the components are described in `01-challenge-and-final-system.md` §3; the full
upload list is in `02-submission-history.md`.

## 1. How to read this file

- **Official Δ** is the change of the composite against the nearest officially scored parent, as assigned in
  `output/data-analysis/results/e4_official.csv`. Component changes are given as credits (the normalised scale of
  `01-challenge-and-final-system.md` §1.4); one sixth of a credit change is its contribution to the composite.
- **Clean** means the pair differs by one lever (checked by a diff of the two build specifications from S4mix on;
  earlier pairs rely on notes). A **bundle** changed several levers at once and cannot be attributed to any of them.
- Every official Δ is a single upload with no seed replicate. The scorer is deterministic (b18, a byte-identical copy
  of B16, scored Δ −1e-16), but model-seed noise on A/B/C was never measured, so gains below a few thousandths are
  fragile.
- **Local evidence** is marked R (real scorer on local harness panels) or X (offline or proxy), as in
  `04-what-did-not-work.md`. J is the offline sum-side proxy J = 2.23 ΔPDS − 1.05 ΔMSE (about 6 × Δcomposite). The
  "panel-weighted forecast" is the composite rule 0.5 × (0.75 × OR + 0.25 × RN mean), pre-registered on 2026-09-30
  (for M2 the OR term was the mean of OR_HCT116X and H1).
- **Scope of "clean".** Of the three method components, only **M2** has a clean single-change official measurement
  (S14fk2 → S15, +0.0106, mostly via MSE). Among the pairs checked by a build-spec diff, five other single changes
  also raised the official score — the cpc pin, the neighbour blend, a DE-source bug fix, the learned amplitude and
  the DE recipe switch — none of which is M1, M2 or M3. Four earlier pairs are single changes according to the notes
  only (B1 → B2, B6 → B8, B8 → B16, b18 → b19); B1 → B2, three sources instead of one, is a precursor of M1.

Sources: `output/data-analysis/results/e4_official.csv`; `output/data-analysis/results/e4_official.json`;
`knowledge/synthesis/e4_official-assumptions.md`.

## 2. All official gains, ranked

| # | Change | Pair | Official Δ | Clean? | Main component moves (credit) |
|---|---|---|---|---|---|
| 1 | Six-source mean + stochastic rounding | B2 → B6 | +0.0750 | bundle (2) | fid +0.297, PDS +0.063, NMAE +0.063 |
| 2 | Delta ×3 | B6 → B8 | +0.0447 | single (notes) | fid +0.192, PDS +0.102, NMAE −0.094 |
| 3 | Dual-moment emission (+ blend 0.3) | S8 → S10b | +0.0348 | bundle (2) | MSE +0.221, PDS +0.104, NMAE −0.152 |
| 4 | cpc composition pin | S10b → S10c | +0.0247 | **single** | NMAE +0.163, reach −0.016 |
| 5 | **M2** agreement amplitude | S14fk2 → S15 | +0.0106 | **single** | MSE +0.056, PDS +0.005 |
| 6 | 3-source mean instead of CD4 alone | B1 → B2 | +0.0084 | single (notes) | PDS +0.054, NMAE +0.030 |
| 7 | Fold-change fusion (M1) added at weight 0.5 | S4gp → S4mix | +0.0077 | bundle (part + weight) | NMAE +0.036, reach +0.008 |
| 8 | Sum-side transforms + kink + no-KOLF + residual | S10f → S10k | +0.0076 (+0.0061 vs S10c) | bundle | PDS +0.045 |
| 9 | Covariance denoiser | B8 → B16 | +0.0075 | single (parent inferred) | PDS +0.051 |
| 10 | Context-forcing residual | b18 → b19 | +0.0073 | single (notes), confounded | NMAE +0.045, PDS −0.012 |
| 11 | 4 families + cens + drc + **M3** | S16 → S18 | +0.0070 | bundle (4) | PDS +0.034, NMAE +0.019, MSE −0.005 |
| 12 | Gated emission + promoter prior | S4 → S4gp | +0.0065 | bundle (3) | PDS +0.022 |
| 13 | S4 without FENG + K562-denoised DE side | S12 → S14fk | +0.0060 | bundle (2) | NMAE +0.022 |
| 14 | DepMap neighbour blend 0.2 | S4mix → S8 | +0.0041 | **single** | spread over five components (MSE 0) |
| 15 | min_cells ensemble + ci5 + lsb | S10k → S11 | +0.0030 | bundle (3) | MSE +0.012 |
| 16 | DE-source link fix | S14fk → S14fk2 | +0.0027 | **single** | reach +0.007, fid +0.005 |
| 17 | Learned per-target amplitude | S15 → S16 | +0.0014 | **single** | MSE +0.004, PDS +0.004 |
| 18 | Extended refit + FiLM residual | S1 → S4 | +0.0013 | bundle (2) | PDS +0.006 |
| 19 | Refit baseline instead of b19 | b19 → S1 | +0.0012 | new pipeline | reach +0.005, PDS +0.002 |
| 20 | DE side from the sum-side recipe | S11 → S12 | +0.0008 | **single** | reach +0.005, NMAE −0.004 |

One removal, one share change and one probe pair also measured the value of components that were kept: removing the
learned S4 half cost −0.0046 (S10c → S10e; PDS −0.033), raising the fusion share to 0.75 cost −0.0063 (S4mix → S7a),
and moving the sum budget off the kink cost −0.0462 upwards and −0.0073 downwards (S16 → S17u, S17d).

Sources: `output/data-analysis/results/e4_official.csv`; `vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`; `AGENTS.md` (2026-09-28 S10k).

## 3. The four large steps

**3.1 Six-source mean and stochastic rounding — B2 → B6, +0.0750 (bundle).** The step changed two things: three
sources became six, and deterministic rounding of counts became stochastic. The project history credits the gain to
the rounding; the pair cannot separate the two. The rounding effect is measured offline (X): round-to-nearest carried
only 22 % of the intended delta (0.0034 of 0.0152) because most genes sit at 0–3 counts, while dithered rounding
recovered 0.0153 of 0.0152.

**3.2 Delta ×3 — B6 → B8, +0.0447 (single, from notes).** The same six-source delta scaled ×3 took the score from
0.0771 to 0.1217. The leaderboard was the evidence: ×1 (B6) and ×4 (B10, a bundle) were both much worse, and the local
KOLF harness had peaked at ×2 instead (0.1972 vs 0.1944 at ×3). A later reconstruction sweep (2026-09-22) read the
scale as amplitude restoration for the direction metrics, not calibration (least-squares slope 0.5; composite optimum
2.5, only +0.004 over ×3). The ×3 survived in the learned S4 part to the end. PDS is not scale-invariant here: on the
KOLF harness raw PDS rose 47 % from ×1 to ×4, because the delta passes through non-linear count reconstruction.

**3.3 Dual-moment emission — S8 → S10b, +0.0348 (bundle).** The emission was rebuilt so that group sums are pinned to
the prediction and the DE signal is created by re-allocating counts within cells
(`01-challenge-and-final-system.md` §3.6); the neighbour blend also moved from 0.2 to 0.3. It brought the first
non-zero MSE credit (raw ratio 6.73 → 0.778; credit 0 → 0.221) and PDS raw 0.764 → 0.810, at the cost of NMAE
(credit +0.127 → −0.025). Local (R): +0.064 over the blend-0.3 arm on 6/6 panels, so about 55 % transferred; the
blend change alone was about +0.003 on H1/KOLF.

**3.4 cpc composition pin — S10b → S10c, +0.0247 (single).** With sums pinned, the fold changes of weak genes had
become sign-random; cpc tilts each gene's column by cell depth so that the mean per-cell composition matches the old
cell-shift expectation while the sums stay fixed. Official NMAE raw 1.017 → 0.916 (credit −0.025 → +0.138); PDS and
MSE identical; reach −0.016. Local (R): +0.020 to +0.038 on 6/6 panels (H1 +0.023 / +0.026 on two seeds), so the
official gain matched the harness. Rank moved from 107 to 49.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§3, §6, §14); `AGENTS.md` (§2; 2026-09-22 "NEW BEST OFFICIAL —
S4"; 2026-09-25 S8; 2026-09-26 S10b and S10c); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§2);
`output/data-analysis/results/e4_official.csv`.

## 4. The method components M1–M3 and the budget

**4.1 M2 — per-target amplitude from cross-source agreement (S14fk2 → S15, +0.0106, single).** Each target's sum-side
change is scaled by rank(agreement)^0.5 before the budget match, so the fixed amplitude budget moves from targets
whose sources disagree to targets whose sources agree. PDS is invariant to a target's own amplitude, so the lever acts
on MSE: official MSE raw 0.749 → 0.696, credit +0.056, i.e. +0.0093 of the +0.0106 composite; PDS +0.0008, reach
+0.0006 and NMAE −0.0001 in composite units. Local (R, S14fk structure): H1 +0.0065, OR_HCT116X +0.0215 / +0.0206
(two seeds), RN_Jurkat +0.0045, RN_RPE1 +0.0025, RN_K562 +0.0047, KOLF +0.0013, RN_HepG2 −0.0045. The pre-registered
gate failed on the per-panel floor (HepG2) although the panel-weighted forecast passed (+0.0054); the official gain was
about twice that forecast and 1.6 times the H1 harness. The follow-up S16 replaced the agreement rank with a learned
predictor of each target's accuracy (β 0.75): +0.0014 officially, exactly the H1 forecast, after 8/8 positive real-
scorer runs. The final pipeline kept the agreement rank.

**4.2 M1 — multi-source fold-change fusion transfer (no clean official measurement).** M1 grew in steps, each inside a
bundle or from notes only:
- B1 → B2, +0.0084: three sources instead of CD4 alone (from notes). The equal source weighting was deliberate: a
  response that reproduces across K562 and primary CD4 T cells is better transfer evidence than one from the largest
  source.
- B2 → B6, +0.0750: six sources, bundled with stochastic rounding (§3.1).
- S4gp → S4mix, +0.0077: the fold-change fusion (true fold changes, source centring, K562-weighted families) entered
  at weight 0.5. The H1 harness (R) had predicted a tie (0.2003 vs 0.2000): the fusion alone raised H1 PDS from 0.777
  to 0.865 but cost DE metrics. Officially PDS stayed flat and the gain came from NMAE (raw 0.949 → 0.926) and reach
  (0.259 → 0.265). A higher share (0.75, S7a) lost 0.0063.
- S10k → S11, +0.0030: the min_cells ensemble, bundled with ci5 and lsb. Official raw MSE 0.7647 → 0.7529 carried the
  whole gain; ci5 and lsb were nil or marginal. Offline (X), the ensemble (seven thresholds plus the default) was the
  first change to improve PDS and MSE together on 3/3 panels (H1 J +0.0159, OR_HCT116 +0.0762, KOLF +0.0416).
- S16 → S18: four more families and the pseudocount ensemble, bundled with drc and M3 (§4.3).

**4.3 M3 — control-PC deflation of the weaker half (not isolable).** M3 entered with S18. Its clean official pair would
be S17 → S18, but S17 was built and never uploaded. The scored pair S16 → S18 (+0.0070) bundles four changes:
(1) the HEPG2, RNJURKAT, RPE1 and DLD1 families at weight 1; (2) the pseudocount ensemble (cens); (3) the DE
calibration drc 0.8 / 1.1; (4) M3 (k 8, η 0.75). The official components moved PDS +0.034, NMAE +0.019, MSE −0.005 and
reach −0.005 (credits); which change produced which part cannot be told. Local evidence: the S17 bundle alone was
positive on 5/5 panels (panel-weighted forecast +0.0055, R); removing the three RN families cost H1 −0.0023 and
OR_HCT116X −0.0013 (R); DLD1 at weight 1 gave OR_HCT116X +0.0040. A later offline re-screen on the S16 base (X,
2026-10-28) put the F1 deflation at +0.011 PDS on OR_HCT116X and +0.0008 on H1.

**4.4 The MSE budget (kink).** The kink rule (sums at 0.8 × the estimated kink, floor 1) and the budget match entered
with S10k, inside a bundle (S10k vs S10c +0.0061; raw PDS 0.810 → 0.8205, MSE 0.778 → 0.7647; pre-upload estimate
+0.006). On A/B/C the estimated kink was about 1.0–1.1, so the rule itself was nearly inert there (sum scale 1.02 in
the S10i build); its role was to hold the amplitude budget while the other sum-side transforms reshaped the sums.
Locally the kink had given the RN panels +0.022 (12/12) and H1 0 (R). Its official value was bracketed by two probes from S16:
margin 1.2 (S17u) cut the MSE credit from 0.311 to 0.038 (−0.0462 composite), and margin 0.6 without the floor (S17d)
lost 0.0073. The same S10k bundle carried the stem-versus-somatic source decision (KOLF off the sum side): keeping
KOLF was worth +0.0079 PDS locally on H1 (stem) and −0.0126 on OR_HCT116 (somatic).

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §16.4, §23.8, §23.10, §23.12–23.13, §25.19);
`AGENTS.md` (2026-09-24 S6; 2026-09-28 S10k and S11; 2026-09-29 "FINAL-PHASE RISK REPRICED"; 2026-10-02; 2026-10-03);
`vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md` (§3); `vcc2026_data/reports/PROJECT_HISTORY.md` (§2);
`output/data-analysis/results/e4_official.json` (`M2`, `M3_bundle`, `checks.S17_note`);
`vcc2026_data/reports/final_2026-10-27/README.md`.

## 5. Other measured gains

- **Covariance denoiser (B8 → B16, +0.0075).** A frozen rank-128 correction learned across sources lifted PDS from
  0.504 to 0.555 credit. Locally (R, KOLF) the learned correction gained 4/4 seeds, but a label-shuffled copy scored
  higher (0.2449 vs 0.2359), so the mechanism is covariance denoising, not learned response transfer.
- **Context-forcing residual (b18 → b19, +0.0073, confounded).** The whole gain was NMAE; PDS fell. b19 also changed
  the baseline rank and training population, and the controlled local test said the residual hurt (five-line PDS
  0.8005 at α 0 vs 0.7947 at α 1).
- **Gated emission and promoter prior (S4 → S4gp, +0.0065).** Genes below 5 CPM stopped receiving invented counts, and
  promoter neighbours were silenced multiplicatively. The new H1 Flex harness predicted +0.006, its first exact
  forecast; raw MSE fell from 14.7 to 8.3.
- **S4 without FENG plus K562-denoised DE side (S12 → S14fk, +0.0060).** Forecast +0.0057; the DE side over-delivered
  (NMAE +0.0037 composite) even though the upload's DE side was mis-linked (S14fk2 then fixed it for +0.0027, against
  a forecast of about +0.003).
- **DepMap neighbour blend (S4mix → S8, +0.0041, clean).** Locally +0.0065 on 18/18 paired runs (R); 63 % transferred.
- **Learned FiLM residual (S1 → S4, +0.0013)** and **the refit baseline (b19 → S1, +0.0012)**: see §7.
- **DE side from the sum-side recipe (S11 → S12, +0.0008, clean).** 6/6 panels positive locally (+0.0024); about a
  third transferred.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§7, §9); `vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md`;
`AGENTS.md` (2026-09-25 S8; 2026-09-29 S12; 2026-09-30 S14fk; 2026-10-01 S14fk2); `output/data-analysis/results/e4_official.csv`.

## 6. Same-line data: the largest local lever, never usable officially

On the H1 panel (R), replacing the transferred responses by those of H1GW (GSE295214, a genome-scale CRISPRi screen of
H1 from another lab covering 131 of the 300 targets) raised the composite from 0.3029 to 0.3772 at weight 1
(+0.074) and to 0.3798 at weight 0.5 (+0.077) during the challenge (2026-09-29); about 70 % of the value was on the DE
side. The pre-registered post-challenge re-measurement, which is newer and preferred, gives +0.0690 at weight 1, about
10.6 times M2's +0.0065 on the same panel (`05-why-it-does-not-work.md`, A4). Related lines helped less and could hurt
(KOLF ← H1GW +0.021 at weight 0.5; OR_HCT116 ← DLD-1 +0.005 at 0.5, −0.035 at 1), so the policy capped non-exact
sources at weight 0.5.

The lever never reached an official entry. H1GW was reserved for stem-like final lines; no public low-MOI same-line
screen existed for the validation lines (Jurkat, HeLa, CAL-33), and all three final lines were somatic with no
same-line screen either.

Sources: `AGENTS.md` (2026-09-29 "LOOP" and "MULTI-ANGLE BOTTLENECK REVIEW"); `knowledge/synthesis/results-summary.md`
(H4); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§25.19); `vcc2026_data/reports/final_2026-10-27/README.md`.

## 7. Local → official shrinkage

| Entry | Change | Local or forecast | Official Δ | Ratio |
|---|---|---|---|---|
| S4 | refit_v2 + FiLM residual (learned context model) | +0.0247 (five-line composite, 15/15 seeds, 5/5 lines) | +0.0013 | **0.053** |
| S3 | scratch FiLM residual | +0.0180 (five-line composite) | −0.0005 | invalid |
| S1 | refit baseline vs b19 | +0.0047 five-line PDS | +0.0020 PDS (+0.0012 composite) | — |
| B15 | effective-cell-count weighting | KOLF PDS 0.1570 → 0.1957 | −0.0150 | sign flip |
| S4gp | gated emission + promoter prior | +0.006 (H1) | +0.0065 | ≈ 1.1 |
| S8 | neighbour blend | +0.0065 (18/18) | +0.0041 | 0.63 |
| S10b | dual-moment emission | +0.064 (6/6) | +0.0348 | ≈ 0.55 |
| S10c | cpc | +0.020 to +0.038 (6/6) | +0.0247 | ≈ 1 |
| S10f | three small knobs (blend, down-gene share, K) | +0.006 (H1, 3/3) | −0.0015 | sign flip |
| S12 | DE side from sum recipe | +0.0024 | +0.0008 | ≈ 0.33 |
| S14fk | S4 no-FENG + denoised DE | +0.0057 (forecast) | +0.0060 | ≈ 1.05 |
| S15 | M2 | +0.0054 (forecast); H1 +0.0065 | +0.0106 | ≈ 2.0 / 1.6 |
| S16 | learned amplitude | +0.0014 (H1) | +0.0014 | 1.0 |

The one valid case for a learned context model is **S4: +0.0247 locally became +0.0013 officially, a ratio of 0.053
(about 19 times smaller); the sign held**. Part of the local figure was the extended refit alone (+0.0060). The S3
case (+0.0180 → −0.0005, ratio −0.028) is **invalid**: the local number came from the drug-context pretraining
experiment, and both that experiment and the S3 build used the scrambled gene-context vectors found on 2026-09-21.

The pattern: the learned context model shrank about 19-fold, and the KOLF-gated weighting (B15) and a set of small
knobs (S10f) reversed sign. Structural emission and sum-side changes, forecast on the H1 Flex harness (from
2026-09-24) and with the panel-weighted forecast (from 2026-09-30), landed at 0.3 to 2 times the forecast.

Sources: `output/data-analysis/results/e4_official.json` (`shrinkage`); `knowledge/synthesis/e4_official-assumptions.md`;
`vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`; `vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`;
`vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md`; `vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md`;
`AGENTS.md` (2026-09-14 reassessment; 2026-09-25 to 2026-10-03 sections; "BUG FOUND 2026-09-21").

## 8. Ride-alongs

From 2026-09-27 no slot was spent on an expected official gain below about +0.005. Verified changes that leave PDS
and MSE untouched or were too small to test alone rode along in the next substantial build:

| Ride-along | What it does | Measured size | Shipped in |
|---|---|---|---|
| cens | pseudocount ensemble c 5/10/20/100 joined to the min_cells members | R: H1 +0.0003, OR_HCT116X +0.0042, RN_Jurkat +0.0063 (3/3; PDS, NMAE and reach positive on all three); forecast composite +0.0024 | S18, F1, F2 |
| drc 0.8 / 1.1 | DE calibration: top 5 % of genes by \|DE LFC\| × 0.8, bottom 70 % × 1.1; sums untouched | R: H1 +0.0018, OR_HCT116X +0.0004, RN_Jurkat +0.0005 (H1 seed noise +0.0002); about +0.001 expected officially | S18, F1, F2 |
| ci5 | composition-pin fixed point against final cell totals | R: H1 +0.0022, KOLF +0.0004; official NMAE moved −0.0000165 (nil) | S11 onwards |
| lsb | depth-balanced emitted block | R: H1 forecast for S11 +0.0044 without it, +0.0053 to +0.0064 with it; official fid −0.0001, jac +0.0001, reach +0.0021 (nil or marginal) | S11 onwards |

Sources: `AGENTS.md` ("SUBMISSION POLICY"; 2026-09-27 "Codex debate loop"; 2026-09-28 S11; 2026-10-03 "12-HOUR LOOP"
and "PDS/MSE ROUND"); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§25.13–25.14).

## Sources cited

`output/data-analysis/results/e4_official.csv`; `output/data-analysis/results/e4_official.json`;
`knowledge/synthesis/e4_official-assumptions.md`; `knowledge/synthesis/results-summary.md`; `AGENTS.md`;
`vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`;
`vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md`; `vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md`;
`vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`; `vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`;
`vcc2026_data/reports/final_2026-10-27/README.md`.
