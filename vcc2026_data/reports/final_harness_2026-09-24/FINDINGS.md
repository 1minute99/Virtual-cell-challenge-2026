# Final-phase harness: ML arms on top of the final-phase recipe (2026-09-24)

User request: implement and test the ML ideas that learn along the target/gene axes (per-target amplitude, low-rank
denoising, expression-based transfer, unseen-target prediction) instead of the context axis.

## Harness
- Six held-out panels with real cells, scored with the official cell-eval2 `vcc2026` preset and the leaderboard
  constants (`tv2_score.CONST`):
  - H1 (VCC2025 Flex, 300 targets x 300 cells);
  - KOLF (279 official targets, 200 cells);
  - RN K562 / HepG2 / Jurkat / RPE1 (the same 300 essential targets, 53-189 cells, 800 controls).
- For a held-out panel, its own source family is excluded everywhere: base prediction, program learning, and the
  features and labels of every learned arm. Learned arms train on the other five panels' real pseudobulk only.
  RN panels still see same-lab sibling lines as sources (optimistic).
- Base = the final-phase recipe without S4 (`tv2_final.py` default): 0.5 x equal-weight uncentred fold change
  (c=100) x3 + 0.5 x centred fusion (c=50, K562 8 ...) x1.6, promoter prior + knockdown override,
  shared_cellshift_gated. Source library rebuilt for 2,581 targets including RN HepG2/Jurkat/RPE1 (`fh_sources.py`).
- One emission seed, paired across arms.
- Code: `scripts/fh_prep.py`, `fh_sources.py`, `fh_harness.py`. Results: `scores/<arm>/<panel>/result.json`,
  `summary.json`.

## Results (official composite; columns H1, KOLF, RN K562, HepG2, Jurkat, RPE1)

| arm | H1 | KOLF | K562 | HepG2 | Jurkat | RPE1 | mean | vs base (mean / worst / wins) |
|---|---|---|---|---|---|---|---|---|
| **base** | 0.207 | -0.019 | 0.370 | 0.280 | 0.359 | 0.158 | **0.226** | - |
| amp (per-target amplitude, energy kept) | 0.056 | -0.136 | 0.281 | 0.260 | 0.249 | 0.099 | 0.135 | -0.091 / -0.151 / 0 |
| amp_shrink (only weakens) | 0.139 | -0.091 | 0.325 | 0.280 | 0.308 | 0.113 | 0.179 | -0.047 / -0.072 / 1 |
| expr (per-gene multiplier) | 0.144 | -0.012 | 0.356 | 0.277 | 0.360 | 0.173 | 0.216 | -0.010 / -0.064 / 3 |
| lr64 (64 response programs) | 0.047 | -0.217 | 0.249 | 0.176 | 0.190 | 0.007 | 0.075 | -0.150 / -0.198 / 0 |
| lr128_half | 0.173 | -0.114 | 0.357 | 0.268 | 0.313 | 0.102 | 0.183 | -0.043 / -0.095 / 0 |

PDS: base 0.821 / 0.555 / 0.854 / 0.733 / 0.809 / 0.874. expr is within +-0.01 everywhere; amp, amp_shrink and
lr arms lower it on every panel (lr64 by 0.04-0.11).

- **Per-target amplitude loses on all six panels.** It fixes over-calling on weak targets (Jaccard on RN 0.06 ->
  0.15), but fidelity (k / max(n_pred, n_conf), direction only) drops 0.1-0.2: under-calling costs more than
  over-calling, and the predicted target strength (train slope regression) is too inaccurate. PDS also drops.
- **Low-rank denoising loses on all six.** Projecting onto shared programs strips the target-specific residual that
  PDS discriminates on. It lowers energy, so MSE improves (RPE1 ratio 0.97-1.00, the first < 1), but the DE
  members pay more.
- **The expression multiplier is a wash** (3/6 wins, mean -0.010, H1 -0.064 via NMAE). It lowers MSE on every panel
  (ratio -0.09 to -0.40) but never below 1.

## Unseen targets (panel targets removed from every source; override kept)

| arm | H1 | KOLF | K562 | HepG2 | Jurkat | RPE1 | mean composite |
|---|---|---|---|---|---|---|---|
| PDS, base (targets seen) | 0.821 | 0.555 | 0.854 | 0.733 | 0.809 | 0.874 | 0.226 |
| **PDS, unseen_knn** (DepMap co-essentiality neighbours) | **0.698** | **0.558** | **0.679** | **0.585** | **0.645** | **0.695** | **-0.020** |
| PDS, unseen_mean (mean response) | 0.535 | 0.501 | 0.511 | 0.502 | 0.510 | 0.511 | -0.219 |
| PDS, unseen_zero (override only) | 0.542 | 0.503 | 0.517 | 0.503 | 0.515 | 0.508 | -0.239 |

- **The only positive ML result.** For a target no screen measured, the weighted mean of its DepMap co-essentiality
  neighbours' responses (median 8-10 neighbours, weight corr^2) gives PDS 0.59-0.70, against ~0.51 for the mean or
  zero fallback. That is +0.20 composite on average, 6/6 panels.
- On KOLF it equals the measured-source base (0.558 vs 0.555): where the context is far from every source, direct
  measurement is no better than the neighbour estimate.

## Conclusion and actions for the final phase
- None of the learned transforms beats the recipe on the six-panel harness; keep the recipe as is.
- Use `unseen_knn` as the fallback for final-panel targets that no public screen covers. Replace "mean / zero"
  in `tv2_final.py`. Coverage of the new panel decides how much this matters.
- The extra RN HepG2/Jurkat/RPE1 sources should enter the final source library (base H1 0.2073 with them, vs
  0.2003 for the same recipe on the old library; not a controlled comparison).
- Caveats: single emission seed; RN panels see same-lab sibling lines; H1 is the only Flex panel.

## Follow-up: measured-source + DepMap-neighbour blend (seen targets), 3 emission seeds

Arms (`fh_harness.py`):
- blend20: 0.8 D + 0.2 x neighbour estimate. Neighbours = top-10 DepMap co-essentiality genes with source data, any
  library target except t, weight corr^2.
- blend20_e: blend20 rescaled to the base's total energy.
- blend_adapt: neighbour weight 0.5 / 0.35 / 0.1 by the number of covering source families.

| arm | seed 0 wins | mean vs base | note |
|---|---|---|---|
| blend20 | 1/6 | -0.017 | averaging shrinks the total energy; fidelity drops |
| blend_adapt | 2/6 | -0.005 | almost every target has >= 3 families, so w = 0.1 |
| **blend20_e** | **6/6** | **+0.0072** | |

blend20_e vs base, paired over seeds 0/1/2: **18/18 wins, mean +0.0065 (sd 0.0023, min +0.0040)**; per panel
H1 +0.0053, KOLF +0.0041, K562 +0.0060, HepG2 +0.0076, Jurkat +0.0075, RPE1 +0.0085.
Per member:
- NMAE better on 18/18 (-0.004);
- fidelity +0.005 (18/18);
- reach +0.011 (17/18);
- Jaccard +0.003 (17/18);
- PDS -0.002 (7/18 positive);
- MSE ratio -0.027.

Replacing 20 % of each target's measured delta with its co-essentiality neighbours' delta, at the same total energy,
improves the direction and DE members on every panel and seed. This is the first learned-free transform to beat
the recipe on the six-panel harness.
Seed noise is large in absolute terms (RPE1 base PDS 0.874 / 0.752 / 0.868) but the paired difference is stable.
The size (+0.0065 local) is small; past local gains shrank 4-19x officially.

Action:
- Add blend20_e to the final recipe (`tv2_final.py`).
- It applies to the public panel too, since all official targets are seen, so one leaderboard slot (S6 + blend20_e)
  would test transfer. That needs the user's go-ahead.

## 10-hour sweep of low-overfit change predictors (2026-09-25; user request)

All arms are paired against base on the same six panels; seed 0 unless noted. Composite deltas:

| arm | idea | wins | mean | H1 | KOLF |
|---|---|---|---|---|---|
| gbm_half / rf_half / gbm | tree stacking of the real LFC from family LFCs + features (0.5 blend / RF / alone) | 3 / 3 / 1 of 6 | -0.001 / +0.003 / -0.044 | -0.020 / -0.018 / -0.112 | +0.016 / +0.017 / +0.022 |
| gain | GBM per-(t,g) multiplier on the base (sign fixed) | 0/6 | -0.033 | -0.056 | -0.003 |
| nnls_fam | 10 learned non-negative family weights | 2/6 | -0.007 | -0.014 | -0.002 |
| gene_smooth | 20 % smoothing over co-responding genes | 2/6 | -0.001 | -0.000 | +0.003 |
| resp_knn | 20 % blend with response-similar targets | 5/6 | +0.010 | -0.001 | +0.003 |
| gbm_resid_tc (lambda 0.5 / 0.2 / 0.1) | target-centred GBM residual | 2 / 3 / 4 of 6 | -0.005 / +0.004 / +0.004 | -0.013 / -0.002 / -0.001 | -0.003 / -0.002 / -0.001 |
| mlp_resid_tc | same residual with an MLP | 3/6 | -0.011 | -0.023 | +0.002 |
| b30r10 | blend30_e + 0.1 residual | 6/6 | +0.014 | +0.008 | +0.006 |
| sim_soft / sim_strong | control-similarity family weights | 5 / 3 of 6 | +0.009 / +0.012 | +0.003 / -0.007 | +0.001 / +0.009 |
| blend30_e / blend40_e | neighbour blend 30 / 40 % | 6 / 6 of 6 | +0.010 / +0.010 | +0.008 / +0.009 | +0.006 / +0.009 |
| combo | sim_soft weights + blend20_e | 6/6 | +0.017 | +0.010 | +0.005 |

Three-seed paired summaries (18 comparisons each):

| arm vs ref | wins | mean | H1 mean | KOLF mean |
|---|---|---|---|---|
| combo vs base | 18/18 | +0.0164 | +0.0083 | +0.0040 |
| blend30_e vs base | 18/18 | +0.0087 | +0.0082 | +0.0070 |
| blend20_e vs base (= S8) | 18/18 | +0.0065 | +0.0053 | +0.0041 |
| sim_soft vs base | 15/18 | +0.0097 | +0.0026 | -0.0006 |
| combo vs blend30_e | 12/18 | +0.0076 | +0.0002 | -0.0030 |
| blend30_e vs blend20_e | 16/18 | +0.0022 | +0.0029 | +0.0029 |

Readings:
- Everything that learns a mapping from features to the change (trees, MLP, gain, learned weights) loses on H1, the
  official-like panel. Where it helps, the gain sits on the RN panels, which share lab, protocol and targets.
  - Stacking and gain lower PDS on every panel: they learn cross-target regularities and blur what discriminates
    targets.
  - The target-centred residual is the exception on PDS (+0.012 to +0.093 at lambda 0.5), but it costs NMAE. At small
    lambda the H1/KOLF effect is nil.
- Neighbour information helps only on the TARGET axis and only when the neighbours are chosen from independent
  information (DepMap co-essentiality). Response-similar neighbours (resp_knn) and gene-axis smoothing do not
  transfer to H1.
- Control-similarity weighting is real but mostly an RN same-lab effect (H1 +0.003, KOLF 0). combo's edge over
  blend30_e is entirely RN; on H1 they tie and on KOLF blend30_e is better.
- Recommendation for the official panel and the final phase: **blend30_e** (S6 recipe + 0.3 neighbour blend,
  energy-matched). Expected over S8 about +0.003 on the official-like panels (H1/KOLF, 6/6 seeds). Keep sim_soft
  optional for a final line with a close sibling source.

## Dual-moment emission (2026-09-25; user: "assume the organisers allow it and try")

`scripts/fh_dual.py`; builder mode `dual_mult` (`tv2_emit.emit_dual_iter`, verified cell-for-cell identical to the
harness function).
- **Group sums:** each target's n-cell total per gene is pinned to round(full control mean x exp(gamma x d) x n), with
  the knockdown / promoter override applied multiplicatively. d = blend30_e delta. These sums drive PDS and the
  pseudobulk MSE.
- **Cell distribution:** the top-K genes by |d| among genes the context expresses (>= 5 CPM) are re-allocated within
  their column, sum unchanged, to create Wilcoxon calls:
  - up-genes are flattened towards uniform;
  - down-genes are concentrated into the top 10 % of cells;
  - strength decays from 1 to 0.3 with rank, so the p-value order follows |d| (reach).

H1 grid (seed 0):

| gamma | K | composite | notes |
|---|---|---|---|
| 0.3 | 300 / 1500 | 0.100 / 0.183 | |
| 0.6 | 300 / 1500 / 3000 | 0.131 / 0.203 / 0.230 | |
| 0.8 | 3000 | 0.238 | |
| **1.0** | **3000** | **0.2450** | MSE ratio 0.761 (scaled 0.237), PDS 0.858, fid 0.538, reach 0.241, NMAE 0.962 |
| 1.0 | 5000 | 0.2452 | same as K 3000 |
| 1.2 / 1.5 | 3000 | 0.228 / 0.203 | 1.5 gives MSE 1.05 |
| 0.6 / 1.0 with called-gene boost x2 + deepest-cell concentration | 3000 | 0.158 / 0.111 | NMAE 1.22-1.29: do not boost |

Six panels, gamma 1.0, K 3000 vs blend30_e (existing emission): **6/6, mean +0.061**:
H1 +0.029, KOLF +0.034, K562 +0.009, HepG2 +0.092, Jurkat +0.041, RPE1 +0.163. PDS rises everywhere (pinned sums
carry no sampling noise). MSE ratio falls below 1 on H1 and on all RN panels (not on KOLF, 1.04).

S10 = S9 recipe (S6 + blend30_e) with dual_mult, gamma 1.0, K 3000. Built, validated by `vcc prep`
(sha 208834ad...); pinned sums verified exact in the package. **Not uploaded; rule status unconfirmed (issue
#348 residual; see ORGANIZER_INQUIRY_DRAFT.md).**

Seed check (dual gamma 1.0 K 3000 vs blend30_e, paired): H1 +0.029 / +0.023 / +0.028, KOLF +0.034 / +0.036 / +0.030 (6/6 positive, mean +0.030). PDS and MSE are identical across seeds (pinned sums are deterministic); only the DE members vary.

### Dual refinements (H1 seed 0 unless noted; base dual gamma 1.0 K 3000 = 0.2450)

| variant | H1 | notes |
|---|---|---|
| dd (down-genes into the deepest cells, sums fixed) | 0.1743 | NMAE 1.20: per-cell LFC overshoots; the earlier boost failures were this flag |
| ag (strength ordered by source sign agreement x \|d\|) | 0.2481 | |
| du0.5 / du0.25 (up-genes flattened towards x ~ L^-a) | 0.2477 / 0.2476 | |
| **du0.5 + ag, three seeds** | **0.2511 / 0.2492 / 0.2525** | vs base 0.2450 / 0.2433 / 0.2481 -> +0.006 / +0.006 / +0.004 |
| du0.5 + ag on KOLF, three seeds | 0.0219 / 0.0188 / 0.0233 | vs 0.0215 / 0.0165 / 0.0230 -> +0.000 / +0.002 / +0.000 |

du0.5 + ag is adopted: 6/6 positive, mean +0.003 (H1 +0.005). Builder options `dual.depth_up` and `dual.agree`
(`tv2_agree_prep.py`: family sign agreement with the family mean). **S10b** = S10 + du0.5 + ag.

Final pipeline (`tv2_final.py`):
- `--blend 0.3` (default) runs `tv2_final_knn.py`: DepMap neighbours, a library including the neighbours, knn.npz.
  Targets no screen measured get the neighbour estimate only.
- `--dual 1.0,3000,0.5,ag` switches to dual emission.

Rehearsal on the simulated D/E/F panel (`scripts/tv2_final_sim2.sh`) built through to the prediction h5ad:
297/300 targets had neighbours, a 972-target library was built, and all targets were covered. `vcc prep` rejects the
fake labels as before.
