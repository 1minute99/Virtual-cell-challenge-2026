# Debate loop rounds

## Round 1 (2026-09-27 18:40-19:30 PDT)
- Codex critique accepted: (1) evidence by laboratory (RN / Orion / H1 / KOLF), Orion-X sensitivity only; (2) harness parity
  with the submitted form (un-matched no-KOLF S4, seed-MEAN residual); (3) internal validation leaks (PCA before split,
  lofo_diag inputs) - external panels unaffected.
- S10k official 0.22300 (+0.0061 over S10c) arrived during the round; the +0.006 estimate transferred.
- Codex H1 "train the residual in the scored geometry" (diagonal w = a/(1+a), a = 5e4 q), `sch_geo_eval.py`, parity form:
  J current / geo: H1 +0.018 / +0.013, RN K562 +0.031 / +0.027, HepG2 +0.017 / +0.008, Jurkat +0.041 / +0.037, RPE1 +0.025 /
  +0.026, Orion +0.075 / +0.057 -> RN lab mean 0.0285 -> 0.0245. REJECTED. Geometry check: cos(exact scored delta, LFC)
  0.86-0.96, cos(exact, diagonal) 0.94-0.996 -> little geometry headroom (Codex's own falsification branch).
- Weak-quartile PDS: the residual moves it by -0.000..+0.008 -> the weak-responder gap is not addressed.

## Round 2-3 (cis prior)
- Real CRISPRi repression beyond the 5 kb prior confirmed (`ml_cis_distance.py`), strongest on weak targets, but the oracle
  upper bound (`ml_cis_oracle.py`) is J +0.002..+0.011 (H1/RN) and +0.008..+0.042 (Orion) at 2-17 pairs per target.
  Codex's caveat recorded: that script replaces the c=50 LFC on the ORIGINAL S4, so it is an optimistic diagnostic, not a
  strict bound. Deprioritised.

## Round 4 (transfer ceiling) - the decisive round
- `ml_weak_reach.py`: weak responders' within-line split-half cosine 0.108 vs cross-line 0.013-0.025 (transferable fraction
  0.12-0.23; strong targets 0.49-0.74).
- `ml_transfer_ceiling.py`: our S10k-parity prediction (PDS 0.9078 / MSE 0.3187, weak 0.818) BEATS perfect same-lab transfer
  from one line (0.819-0.886) and from the mean of three (0.8974 / 0.4075, weak 0.788), also after per-target energy
  matching (0.9004). A 50/50 blend with ours is 0.9134 / 0.2992 (J +0.033) -> the only transfer gain left is ensemble
  diversity against lines we already use.
- Codex's objections to be honest about: the within/cross cosines share the control estimate and differ in depth
  attenuation; the oracle rows pass through OUR transforms and use c=50 source LFCs; the baseline is the original S4, not
  the S10k construction; RN lines are essential-gene panels. Conclusion stated as "unconditioned transfer is saturated",
  not as a proof about context-conditioned prediction on D/E/F.
- **Codex's counter-proposal (its highest-value non-transfer experiment): the composition pin's denominator.**
  `tv2_emit.emit_dual_iter` solves mean_c(W_cg / Lc) = CPC with Lc the ORIGINAL control row totals, while the scorer's DE
  functional uses the FINAL emitted row totals. Audited (`ml_cpc_audit.py`, H1, 20 targets, 10,783 genes): realized /
  intended composition median 1.0011 (no bias), within-target sd of log 0.0589, of which 0.047 comes from the denominator
  substitution (the bisection's own residual is 0.0356); final row totals are 0.970 of the control totals.
  The existing `ci` fixed point (already in `sch_emit.emit`, never scored - the `_ci2` dir has no result.json) reduces it:
  sd 0.0569 (ci0) -> 0.0542 (ci2) -> 0.0488 (ci5) (`ml_cpc_ci_check.py`). Implied NMAE gain ~+0.001 composite, below the
  adoption bar; one real-scorer H1 run queued because the offline evaluator cannot see the DE members.

## Round 5 (the composition pin) - Codex's hypothesis CONFIRMED, then bounded
- `ci5` (fixed point against the emitted cells' own row totals, already in sch_emit but never scored) on the REAL scorer,
  H1 seed 0: composite 0.2741 -> 0.2763 (+0.0022); PDS and MSE EXACTLY unchanged (sums pinned), NMAE 0.8743 -> 0.8697,
  reach 0.2630 -> 0.2672, fid/jaccard +0.0001. KOLF and a second H1 seed queued.
- Residual error after ci5 is integer rounding: |log err| by mean counts per cell 0.0194 / 0.0116 / 0.0055 / 0.0022 / 0.0017
  vs the 1/sqrt prediction 0.0194 / 0.0093 / 0.0053 / 0.0029 / 0.0017 (`ml_cpc_residual_src.py`) - irreducible at 400 cells.
- Codex's two further fixes, both real in principle, both measured at ZERO (`ml_cpc_fixes.py`): override-consistent target
  changes scored genes by sd log 6e-5 (H1 has 1 override gene/target at 0.24 % of the control total; and cs_cpm_expected
  already normalises by the SHIFTED total, so its denominator was never Lc); the logL vs log(Lcur) tilt is non-monotone in
  0 of 15,000 solves and the continuous solve is exact either way. Direction closed at its floor.

## Round 6 (where is the remaining PDS?) - the per-target SOURCE WEIGHTING oracle
Codex's independent answer to "how do the top-10 get PDS 0.861 vs our 0.8205": "identifying which source response is
relevant to each target in the destination context". I measured exactly that.
- `ml_fusion_oracle.py`: replacing the recipe's fixed family weights by the per-target NNLS fit against the destination's
  REAL delta (energy-matched, everything else identical) gives H1 PDS 0.9067 -> 0.9550 (+0.048), Orion 0.7409 -> 0.8709
  (+0.130). Controls: weights fitted on a random HALF of the genes keep +0.035 (H1) / +0.129 (Orion); weights shuffled
  across targets collapse to -0.374 / -0.224. So the per-target weighting carries real, target-specific signal.
- But it is NOT reachable at test time:
  * MODZ (LINCS-style mutual-agreement weighting, listed as untried in the original handoff) - `ml_modz_weights.py`:
    worse than fixed weights in the SAME construction on both panels at powers 0.5 / 1 / 2 (H1 0.8589 vs 0.8640,
    Orion 0.7007 vs 0.7117).
  * Leave-one-panel-out GBM on test-time features (family response norm and rank, mutual agreement, coverage, source-vs-
    destination control-expression similarity, DepMap, recipe weight) - `ml_weight_lpo.py`, 6 panels: oracle recovered
    +3 %, -8 %, -20 %, -10 %, -84 %, +14 % (mean ~ -17 %), indistinguishable from the shuffled control.
=> Third independent confirmation that destination-specific information is not present in the sources. The remaining PDS
gap to the top-10 is most consistent with MORE source contexts (rank 1 uses a proprietary 200+ cell-line atlas), not with a
better combination of the contexts we have.

## Codex's final-phase recommendation (round 6, accepted as the next build)
Stop emission tuning after the ci5 confirmations; build ONE blind final-phase rehearsal: several pseudo-finals that give the
pipeline only controls and a new 300-target list, with outcomes inaccessible until predictions are frozen, covering three
coverage regimes (many sources / one or poorly matched source / unmeasured -> neighbour fallback), run through the ACTUAL
final builder and scorer, selecting on laboratory-level performance and regressions rather than the pooled mean.
