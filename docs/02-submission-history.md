# 02 — Submission history

The team made **38 official submissions**: 36 in the validation phase (contexts A/B/C) and 2 in the final phase
(contexts D/E/F), first F1 and then F2. **F2 is the last final entry and the one that counts.** The best validation
score was **S18, 0.2545**. At the end of the validation phase the team stood **129th of 1,635 (top 7.9 %)**; that is
the final validation standing, not the best one, which was **34th of 1,244 (top 2.7 %) on 2026-09-30**. The
final-phase scores have not been published (§4).

Sources: `output/data-analysis/results/e4_official.json`; `knowledge/synthesis/e4_official-assumptions.md`.

## 1. How to read the tables

- **Labels** follow `output/data-analysis/results/official_record.csv`. For S3 onwards the `label` column holds the
  opaque model name and the S-label is in `s_label`; the S-label is used here. "prior A" and "prior B" are the CSV's
  first two rows, whose labels mean "earlier A/B". S4gp and S4mix are called S5 and S6 in the working notes.
- **Score** is the official composite (mean of the six component credits, `01-challenge-and-final-system.md` §1.4).
  **Rank** is the rank shown when the score was first retrieved, not a leaderboard-snapshot standing (§3).
- **Δ** is the change against the nearest officially scored parent, as assigned in
  `output/data-analysis/results/e4_official.csv`. **Attribution**: *single* = one lever changed (checked by a diff of
  the two build specifications where both exist, otherwise taken from the notes, which is lower confidence);
  **bundle** = several levers changed together, so the gain cannot be attributed to any one of them; *root* = a new
  pipeline with no comparable parent.
- **Dates** are UTC. Upload dates for prior A/B and B1–B14 were not recorded on disk.

Sources: `output/data-analysis/results/official_record.csv`; `output/data-analysis/results/e4_official.csv`;
`knowledge/synthesis/e4_official-assumptions.md`.

## 2. All 38 submissions

| # | Date | Label | What changed (short) | Score | Rank | Δ (parent) | Attribution |
|---|---|---|---|---|---|---|---|
| 1 | Sept 5–7¹ | prior A | Training-free transfer of Replogle K562 responses, deterministic rounding | −0.0228 | 468 | — | root |
| 2 | Sept 5–7¹ | prior B | STATE ST single-cell model trained from scratch (4 lines, ESM-2 features) | −0.1910 | 609 | — | root |
| 3 | ≤ 09-14² | B1 | CD4 atlas alone, mean transfer, deterministic rounding | −0.0064 | 484 | — | root |
| 4 | ≤ 09-14 | B2 | 3-source mean, deterministic rounding | 0.0020 | 467 | +0.0084 (B1) | single (notes) |
| 5 | ≤ 09-14 | B6 | 6-source mean + stochastic rounding | 0.0771 | 342 | +0.0750 (B2) | **bundle** (2) |
| 6 | ≤ 09-14 | B7 | Sign-agreement weighting (stochastic rounding) | −0.0224 | 550 | −0.0994 (B6) | single (notes) |
| 7 | ≤ 09-14 | B8 | B6 with the delta scaled ×3 | 0.1217 | 214 | +0.0447 (B6) | single (notes) |
| 8 | ≤ 09-14 | B10 | Plain mean transfer at ×4 | 0.0868 | 328 | −0.0349 (B8) | **bundle** (2), parent inferred |
| 9 | ≤ 09-14 | B11 | Expression-gated per-context delta, ×2.25 | 0.1156 | 234 | −0.0061 (B8) | **bundle** (2), parent inferred |
| 10 | ≤ 09-14 | B14 | Precision-weighted sources, ×3.5 | 0.0888 | 345 | −0.0329 (B8) | **bundle** (2), parent inferred |
| 11 | 09-14 | B15 | Effective-cell-count source weighting, ×3 | 0.1067 | 287³ | −0.0150 (B8) | single, parent inferred |
| 12 | 09-15 | B16 | Frozen rank-128 covariance denoiser on the 6-source mean, ×3 | 0.1292 | 198³ | +0.0075 (B8) | single, parent inferred |
| 13 | 09-16 | b18 | Package byte-identical to B16 (re-validation) | 0.1292 | 231 | 0.0000 (B16) | determinism check |
| 14 | 09-17 | b19 | Context-forcing residual (rank 32, α 1), ×3 | 0.1365 | 193 | +0.0073 (b18) | single (notes)⁴ |
| 15 | 09-20 | S1 | Refit baseline, the best local predictor, first submitted | 0.1377 | 231 | — (+0.0012 vs b19) | root |
| 16 | 09-20 | S2 | Refit + PCA context interaction | 0.1357 | 241 | −0.0020 (S1) | invalid⁵ |
| 17 | 09-21 | S3 | Refit + 0.25 × scratch FiLM context residual | 0.1372 | 241 | −0.0005 (S1) | invalid⁵ |
| 18 | 09-22 | S4 | Extended-table refit + 0.25 × FiLM residual (38 gene contexts) | 0.1390 | 265 | +0.0013 (S1) | **bundle** (2) |
| 19 | 09-24 | S4gp | Gated cell-shift emission + multiplicative promoter-neighbour prior (S4 deltas ×3) | 0.1455 | 286 | +0.0065 (S4) | **bundle** (3) |
| 20 | 09-24 | S4mix | 0.5 × S4 + 0.5 × fold-change fusion (M1), ×1.6 | 0.1532 | 249 | +0.0077 (S4gp) | **bundle** (part + weight) |
| 21 | 09-25 | S7a | Fusion share 0.75 | 0.1469 | 288 | −0.0063 (S4mix) | single (spec) |
| 22 | 09-25 | S8 | + DepMap co-essentiality neighbour blend, weight 0.2 | 0.1573 | 225 | +0.0041 (S4mix) | single (spec) |
| 23 | 09-26 | S10b | Dual-moment emission; blend 0.3 (via unsubmitted S9) | 0.1922 | 107 | +0.0348 (S8) | **bundle** (2) |
| 24 | 09-26 | S10c | + cpc composition pin | 0.2169 | 49 | +0.0247 (S10b) | single (spec) |
| 25 | 09-27 | S10e | Ablation: learned S4 half replaced by equal-weight fold changes | 0.2123 | 57 | −0.0046 (S10c) | single (spec) |
| 26 | 09-27 | S10f | Blend 0.5 + down-gene share 0.2 + re-allocation K 5000 | 0.2154 | 57 | −0.0015 (S10c) | **bundle** (3) |
| 27 | 09-28 | S10k | Centring 0.5, guarded library-PC removal, budget-matched kink (from unsubmitted S10i) + sum side without KOLF + no-embedding residual 0.25 | 0.2230 | 50 | +0.0076 (S10f)⁶ | **bundle** |
| 28 | 09-28 | S11 | min_cells ensemble + ci5 + lsb | 0.2260 | 41 | +0.0030 (S10k) | **bundle** (3) |
| 29 | 09-29 | S12 | DE side from the sum-side recipe | 0.2268 | 38 | +0.0008 (S11) | single (spec) |
| 30 | 09-30 | S14fk | S4 without FENG for somatic contexts + K562-denoised DE side⁷ | 0.2327 | 34 | +0.0060 (S12) | **bundle** (2) |
| 31 | 10-01 | S14fk2 | DE-source link fix (the intended DE side) | 0.2354 | 35 | +0.0027 (S14fk) | single (spec) |
| 32 | 10-02 | S15 | M2: agreement amplitude, β 0.5 | 0.2460 | 32 | +0.0106 (S14fk2) | single (spec) |
| 33 | 10-03 | S16 | Learned per-target amplitude, β 0.75 | 0.2474 | 33 | +0.0014 (S15) | single (spec) |
| 34 | 10-10 | S17u | Probe: sum budget up (kink margin 1.2) | 0.2012 | 260 | −0.0462 (S16) | single (spec) |
| 35 | 10-11 | S17d | Probe: sum budget down (margin 0.6, floor 0) | 0.2402 | 119 | −0.0073 (S16) | single (spec) |
| 36 | 10-11 | **S18** | S16 + 4 families + cens + drc (= unsubmitted S17) + M3 (k 8, η 0.75) | **0.2545** | 74 | +0.0070 (S16) | **bundle** (4) |
| 37 | 10-27 | F1 | Final D/E/F: final-phase command + S17 families + M3 + cens + drc | hidden | — | — | final, superseded |
| 38 | 10-28 | **F2** | F1 + sum-side expression down-weighting kept | hidden | — | — | final, **counting** |

¹ Era 1 of `vcc2026_data/reports/PROJECT_HISTORY.md` (§5); exact dates not recorded.
² Exact dates not recorded; all of B1–B14 had been scored when the handoff of 2026-09-14 was written.
³ B15's rank is 277 in two records and 287 in two others; B16 is 198 in the record and 199 in one table.
⁴ b19 also changed the baseline rank and training population, so the pair cannot isolate the residual.
⁵ S2 and S3 were built with scrambled gene-context control profiles (bug found 2026-09-21); their scores say nothing
about the methods.
⁶ +0.0061 against S10c, the build's actual base; S10i reverted S10f's knob changes.
⁷ The uploaded S14fk had broken links in its DE-source folder, so its DE side saw only the denoised K562 screen;
S14fk2 shipped the intended DE side.

Sources: `output/data-analysis/results/official_record.csv`; `output/data-analysis/results/e4_official.csv`;
`vcc2026_data/reports/PROJECT_HISTORY.md` (§5–§9, §13, §14); `vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`;
`AGENTS.md` (§2; 2026-09-14 to 2026-10-03 sections; "BUG FOUND 2026-09-21"); `vcc2026_data/reports/final_2026-10-27/README.md`.

### Component credits

Normalised credits (0 = predicting no change, 1 = split-half replicate). MSE stayed at exactly 0 until S10b.

| Label | PDS | MSE | NMAE | fid | reach | jac |
|---|---|---|---|---|---|---|
| prior A | 0.294 | 0.000 | 0.043 | −0.512 | 0.064 | −0.025 |
| prior B | 0.014 | 0.000 | −0.016 | −1.077 | −0.008 | −0.059 |
| B1 | 0.285 | 0.000 | 0.037 | −0.404 | 0.059 | −0.015 |
| B2 | 0.339 | 0.000 | 0.067 | −0.473 | 0.083 | −0.004 |
| B6 | 0.402 | 0.000 | 0.130 | −0.175 | 0.137 | −0.031 |
| B7 | 0.277 | 0.000 | 0.065 | −0.496 | 0.060 | −0.040 |
| B8 | 0.504 | 0.000 | 0.036 | 0.016 | 0.180 | −0.006 |
| B11⁸ | 0.508 | 0.000 | 0.005 | 0.017 | 0.162 | 0.002 |
| B14⁸ | 0.462 | 0.000 | −0.158 | 0.035 | 0.196 | n/r |
| B15 | 0.458 | 0.000 | −0.042 | 0.035 | 0.191 | −0.003 |
| B16 | 0.555 | 0.000 | 0.030 | 0.011 | 0.187 | −0.008 |
| b18 | 0.555 | 0.000 | 0.030 | 0.011 | 0.187 | −0.008 |
| b19 | 0.543 | 0.000 | 0.075 | 0.019 | 0.189 | −0.007 |
| S1⁹ | 0.545 | 0.000 | 0.076 | 0.018 | 0.194 | −0.007 |
| S2⁹ | 0.534 | 0.000 | 0.074 | 0.018 | 0.195 | −0.007 |
| S3 | 0.544 | 0.000 | 0.075 | 0.018 | 0.194 | −0.007 |
| S4 | 0.551 | 0.000 | 0.078 | 0.017 | 0.196 | −0.008 |
| S4gp | 0.573 | 0.000 | 0.087 | 0.018 | 0.203 | −0.008 |
| S4mix | 0.576 | 0.000 | 0.123 | 0.020 | 0.210 | −0.009 |
| S7a | 0.566 | 0.000 | 0.112 | 0.014 | 0.200 | −0.010 |
| S8 | 0.582 | 0.000 | 0.127 | 0.025 | 0.217 | −0.008 |
| S10b | 0.686 | 0.221 | −0.025 | 0.035 | 0.240 | −0.004 |
| S10c | 0.686 | 0.221 | 0.138 | 0.036 | 0.224 | −0.004 |
| S10e | 0.653 | 0.229 | 0.141 | 0.039 | 0.216 | −0.004 |
| S10f | 0.666 | 0.237 | 0.133 | 0.034 | 0.221 | 0.002 |
| S10k | 0.711 | 0.235 | 0.138 | 0.037 | 0.221 | −0.004 |
| S11 | 0.713 | 0.247 | 0.139 | 0.036 | 0.223 | −0.004 |
| S12 | 0.713 | 0.247 | 0.135 | 0.039 | 0.229 | −0.002 |
| S14fk | 0.718 | 0.251 | 0.157 | 0.041 | 0.232 | −0.003 |
| S14fk2 | 0.718 | 0.251 | 0.159 | 0.045 | 0.239 | 0.000 |
| S15 | 0.723 | 0.307 | 0.158 | 0.046 | 0.242 | 0.000 |
| S16 | 0.727 | 0.311 | 0.157 | 0.046 | 0.243 | 0.000 |
| S17u | 0.724 | 0.038 | 0.157 | 0.046 | 0.242 | 0.000 |
| S17d | 0.731 | 0.262 | 0.158 | 0.046 | 0.244 | 0.000 |
| S18 | 0.761 | 0.306 | 0.176 | 0.046 | 0.238 | 0.000 |

No component data were recorded for B10. ⁸ B11 from the handoff's decomposition of that entry (PDS 0.5075, NMAE
0.0050, fid 0.0165, reach 0.1621, jac 0.0023); B14 from its per-metric comparison with B8 (n/r = not recorded).
⁹ From the S1/S2 findings note (the CSV has no components for them).

Sources: `output/data-analysis/results/e4_official.csv`; `output/data-analysis/results/official_record.csv`;
`AGENTS.md` (§1, §2); `vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`.

## 3. Standing on the leaderboard

The validation leaderboard showed **each team's latest entry, not its best**. Ten snapshots were saved; only
aggregate statistics and our own row are reported.

| Snapshot | Teams | Our displayed entry | Our score | Our rank | Top % | Top-10 threshold | 90th pct | Median |
|---|---|---|---|---|---|---|---|---|
| 2026-09-18 | 1,036 | b19 | 0.1365 | 205 | 19.8 % | 0.2286 | 0.1665 | 0.0479 |
| 2026-09-22 | 1,104 | S4 | 0.1390 | 266 | 24.1 % | 0.2419 | 0.1759 | 0.0564 |
| 2026-09-27 | 1,171 | S10f | 0.2154 | 52 | 4.4 % | 0.2533 | 0.1942 | 0.0679 |
| 2026-09-28 | 1,188 | S10k | 0.2230 | 48 | 4.0 % | 0.2555 | 0.1984 | 0.0698 |
| 2026-09-29 | 1,214 | S11 | 0.2260 | 40 | 3.3 % | 0.2599 | 0.2042 | 0.0709 |
| 2026-09-30 | 1,231 | S12 | 0.2268 | 38 | 3.1 % | 0.2634 | 0.2072 | 0.0747 |
| **2026-09-30b** | 1,244 | S14fk | 0.2327 | **34** | **2.7 %** | 0.2646 | 0.2062 | 0.0769 |
| 2026-10-01 | 1,250 | S14fk2 | 0.2354 | 35 | 2.8 % | 0.2729 | 0.2069 | 0.0781 |
| 2026-10-10 | 1,434 | S17u | 0.2012 | 252 | 17.6 % | 0.2942 | 0.2315 | 0.0899 |
| **2026-10-27** | 1,635 | S18 | 0.2545 | **129** | **7.9 %** | 0.3193 | 0.2392 | 0.0807 |

- **Best standing: 34 of 1,244 (top 2.7 %)**, snapshot of 2026-09-30 (second save that day), with S14fk displayed.
  S15 and S16 showed ranks 32 and 33 when their scores were retrieved, but none of the ten saved snapshots (which
  carry the team counts) falls on those days, so they are not used as standings.
- **The latest-entry rule.** Twice the displayed entry was not our best: on 2026-09-27 (S10f, 0.2154, instead of
  S10c, 0.2169; a small effect) and on 2026-10-10, when the board displayed the S17u budget probe (0.2012) and the
  rank fell to 252 of 1,434 (17.6 %). S18 replaced S17u as the displayed entry the next day; S17d and S18 were
  uploaded in that order so that S18 would be the entry shown.
- **End of validation: 129 of 1,635 (top 7.9 %; above 92.1 % of teams)**, with S18 — both our latest and our best
  entry — displayed. The gap to the top-10 threshold (0.3193) was 0.0649. The fall from 34th is the field moving: from
  2026-09-30 to 2026-10-27 the number of teams grew from 1,244 to 1,635 and the top-10 threshold rose from 0.2646 to
  0.3193, while our best score rose from 0.2327 to 0.2545.
- Our submission count on the final validation board was 36, matching the record.

Sources: `output/data-analysis/results/e4_official.json` (`leaderboard_snapshots_aggregate`,
`final_validation_standing`, `checks`); `knowledge/synthesis/e4_official-assumptions.md`;
`vcc2026_data/reports/submission_2026-10-11/upload_s17d_s18.sh`.

## 4. The trajectory by period

**Era 1, Sept 5–7 — foundation models (2 entries).** Both first entries scored below zero, i.e. worse than predicting
no change: a training-free transfer of Replogle K562 responses (prior A, −0.0228) and a STATE ST model trained from
scratch (prior B, −0.1910). The arithmetic method beat the 45.6M-parameter trained model on five of the six
components; both scored 0 on MSE.

**Era 2, to Sept 14 — mean transfer and count rounding (8 entries).** Averaging each target's measured change over
public screens and adding it to the context's control cells gave the first positive scores: CD4 alone −0.0064, three
sources 0.0020, six sources with stochastic rounding 0.0771, and the same delta ×3 0.1217 (B8, rank 214). The
alternatives at ×4, with expression gating at ×2.25 and with precision weights at ×3.5 all lost to B8; each changed
two things at once. Sign-agreement weighting (B7) lost 0.0994.

**Era 3, Sept 14–15 — weighting and denoising (2 entries).** B15 tested effective-cell-count weighting as one change
and lost 0.0150, although the KOLF harness had predicted a gain; KOLF was then demoted as a promotion gate. B16 added a
frozen rank-128 covariance denoiser and reached 0.1292 (+0.0075), the first leaderboard gain from learned
post-processing. A training-label-shuffle control scored at least as well locally, so the mechanism was recorded as
covariance denoising, not a learned response to new contexts.

**Eras 4–5, Sept 15–17 — context models (2 entries).** The context campaign ran locally (`04-what-did-not-work.md`,
group (b)). b18 re-submitted B16's package byte for byte and scored identically (Δ −1e-16): a slot spent, which at
least showed that the scorer is deterministic. b19 added a context-forcing residual and scored 0.1365 (rank 193); its
whole gain was in NMAE, and its PDS was below B16's.

**Eras 6–8, Sept 18–20 — the refit baseline (2 entries).** After sweeps, a data audit and the oracle analysis, the
best local predictor was submitted for the first time: the refit baseline S1, 0.1377. S2 added a PCA context
interaction and lost 0.0020; it was later found to have used scrambled context profiles.

**Sept 21–22 — learned context residuals (2 entries).** S3 (a FiLM residual trained from scratch, also on the
scrambled profiles) lost 0.0005. S4 (extended-table refit plus a FiLM residual with correct profiles) gained 0.0013,
the first context-conditioned model above the refit, at about 1/19 of its local gain.

**Sept 23–27 — transfer and emission rebuilt (8 entries).** The leaderboard's per-component data, read on 2026-09-22,
put our raw MSE ratio at 14.7 against 0.58–0.79 for the top 30, and the gap analysis of 2026-09-23 found defects in
count reconstruction, prediction energy, source weighting and promoter neighbours. Rebuilding the transfer (S4gp
0.1455, gated emission and promoter prior; S4mix 0.1532, the fold-change fusion M1; S8 0.1573, the neighbour blend)
and then the emission (S10b 0.1922, the first non-zero MSE credit; S10c 0.2169, rank 49) added 0.078 between 2026-09-22 and
2026-09-26. Two follow-ups lost (S10e, S10f). After them the user set the rule of 2026-09-27: no slot for an expected
gain below about +0.005; small, verified, PDS/MSE-neutral fixes ride along with the next substantial build.

**Sept 28 – Oct 3 — sum-side refinement and M2 (7 entries).** S10k (centring, library-PC removal, budget-matched kink,
no-KOLF sum side, learned residual) 0.2230; S11 (min_cells ensemble) 0.2260; S12 0.2268; S14fk 0.2327, which gave the
best standing of the phase; S14fk2 (a bug fix) 0.2354; S15 (M2) 0.2460, the largest single-change gain since S10c; S16
(learned amplitude) 0.2474.

**Oct 4–11 — source set, budget probes and M3 (3 entries).** On 2026-10-06 the team found that four built source
families with official-panel coverage (HEPG2, RNJURKAT, RPE1, DLD1) carried no weight in the shipped specification,
while the local harness had been using three of them. Two entries were spent on purpose to measure the official MSE
budget: S17u (margin 1.2) cut the MSE credit from 0.311 to 0.038 (−0.0462 composite); S17d (margin 0.6,
floor 0) lost 0.0073. S17 (S16 + the four families + cens + drc) was built but not uploaded; S18 = S17 + M3 scored
0.2545, the best validation score. Validation closed on 2026-10-22 with S18 displayed.

**Oct 22 – Nov 5 — final phase (2 entries).** F1 was uploaded on 2026-10-27 and F2 on 2026-10-28; F2 is the counting
entry (`01-challenge-and-final-system.md`, §4). As of 2026-11-22 the final-phase scores have not been published; the
rules release the final leaderboard, based on each team's last final entry, when the winners are announced in late
November 2026.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§5–§13); `AGENTS.md` (2026-09-14 to 2026-10-27 sections,
"SUBMISSION POLICY"); `vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/submission_2026-09-21/FINDINGS.md`; `vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`;
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§25.19); `output/data-analysis/results/e4_official.json`;
`vcc2026_data/reports/final_2026-10-27/README.md`.

## 5. Record notes

- **Count.** 38 = 36 validation + 2 final; our count on the final validation board was 36.
- **Missing ids.** B10, B11 and B14 have no entry id on disk; their scores and ranks come from
  `vcc2026_data/reports/PROJECT_HISTORY.md`.
- **Lineage confidence.** The B-series parents come from one-line notes; B15's and B16's parent (B8) is inferred.
  From S4mix onwards the levers were checked by diffing build specifications.
- **No seed replicates.** Every official Δ is one upload. b18 ≡ B16 shows the scorer is deterministic; model-seed
  noise on the official contexts was never measured.
- **Built but never uploaded** (not counted): b9, b12, b13, b17, S9, S10, S10d, S10i, S10j and S17.

Sources: `output/data-analysis/results/e4_official.json` (`checks`); `knowledge/synthesis/e4_official-assumptions.md`
("What stays open"); `vcc2026_data/reports/PROJECT_HISTORY.md` (§13); `AGENTS.md` (2026-09-25, 2026-09-26,
2026-09-27 sections); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§3).

## Sources cited

`output/data-analysis/results/official_record.csv`; `output/data-analysis/results/e4_official.csv`;
`output/data-analysis/results/e4_official.json`; `knowledge/synthesis/e4_official-assumptions.md`; `AGENTS.md`;
`vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/submission_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/submission_2026-09-21/FINDINGS.md`; `vcc2026_data/reports/submission_2026-09-22/FINDINGS.md`;
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`; `vcc2026_data/reports/submission_2026-10-11/upload_s17d_s18.sh`;
`vcc2026_data/reports/final_2026-10-27/README.md`.
