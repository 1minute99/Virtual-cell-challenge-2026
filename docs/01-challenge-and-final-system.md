# 01 — The challenge and the final system

This file describes the Virtual Cell Challenge 2026 as the team met it, the public data the predictions were built
from, and the system behind the final-phase entry F2. The upload sequence is in `02-submission-history.md`, the
measured value of each change in `03-what-raised-the-score.md`, and the negative results and their causes in
`04-what-did-not-work.md` and `05-why-it-does-not-work.md`.

## 1. The challenge

### 1.1 Task

The task was to predict the single-cell transcriptome (18,533 genes) after CRISPRi knockdown of each of 300 target
genes, in cell lines for which no perturbation data were released. For each context the organisers released only
non-targeting control cells (18,400 per context) and the target list. There was no challenge training set: every
perturbation response used here came from public screens in other cell lines. The validation ground truth was never
released; it existed only as the leaderboard.

Sources: `RESULTS.md` (section "The blocker you need to know about first"); `AGENTS.md` (§1 "The task, in one
paragraph"); `vcc2026_data/reports/final_2026-10-27/README.md`.

### 1.2 Submission format and quota

- 400 predicted cells per target per context, as raw non-negative integer counts: 3 contexts × 300 targets × 400
  cells = 360,000 cells × 18,533 genes.
- The official CLI caps a file at 400,000 cells and each cell's total count at 1,000,000.
- The 2026 data were generated with 10x Flex chemistry.
- Validation quota: two submissions per day, one in flight at a time.

Sources: `VCC2026_gene_perturbation_report.md` (§1); `vcc2026_data/reports/PROJECT_HISTORY.md` (§1).

### 1.3 Phases

| | Validation phase | Final phase |
|---|---|---|
| Contexts | A, B, C (publicly identified as Jurkat, HeLa, CAL-33) | D, E, F (three new lines) |
| Panel | `vcc2026-val-1`, 300 targets | `vcc2026-test-1`, 300 targets, 129 of them shared with A/B/C |
| Dates | until 2026-10-22 | controls released 2026-10-22; deadline 2026-11-05 23:59 UTC |
| Leaderboard | live; each team's latest entry is shown | scores hidden until the winners are announced; only the last final entry counts; 2 entries per UTC day |

The team identified the final lines from their control cells by Spearman correlation with Human Protein Atlas
cell-line profiles: **D = MOLM-13** (AML, myeloid; ρ 0.849), **E = U-87MG** (ρ 0.878), **F = an hTERT-immortalised
fibroblast** (HBF/TERT88, BJ-hTERT-like; ρ 0.827). A stem score computed from control-cell marker panels classified
all three as somatic (−1.8, −5.9, −6.0). No public low-MOI screen of any of the three lines was found.

Sources: `AGENTS.md` (2026-09-23 "Leaderboard gap analysis"; 2026-10-27 "FINAL PHASE BUILD F1");
`vcc2026_data/reports/final_2026-10-27/README.md`.

### 1.4 The official scorer

Entries were scored by `cell-eval2` (version 0.16.0 as recorded in September) with `--preset vcc2026`. Each of
six components is normalised so that **0 = the context-mean baseline** (predicting no change) and **1 = a split-half
replicate** of the real data. The overall score is their unweighted mean; on every one of our entries with all six
components recorded, the mean reproduces the overall score to 8e-17.

| Short | Scorer key | What it measures | Raw value at 0 | Raw value at 1 |
|---|---|---|---|---|
| PDS | `pds_cosine` | Perturbation discrimination. For each target, whether the prediction is closer (cosine distance) to that target's true change than to the other targets' true changes, scored 1 − rank/(n − 1). All 300 panel target genes are removed first, so a target's own knockdown earns nothing; a fully tied row scores 0.5. | 0.500 | 0.9515 |
| MSE | `expr_mse_unbiased_capped_norm` | Error of the predicted pseudobulk (log-normalised per-target group sums) against the truth, after subtracting a capped allowance for sampling noise, relative to predicting no change. Lower raw is better. The credit was exactly 0 for every entry before S10b, whose raw error was far above the no-change level. | ≈ 0.984 | 0.043 |
| NMAE | `de_wilcoxon_lfc_nmae` | Error of the differential-expression log2 fold changes: mean absolute difference between predicted and true values over the genes that are truly significant, divided by the mean absolute true value and referenced to a split-half replicate. Targets with ≤ 10 truly significant genes are dropped. In practice it is limited by sign accuracy. | 1.0008 | 0.3921 |
| fid | `de_wilcoxon_direction_fidelity_yield_raw` | Direction fidelity: k / max(n_pred, n_conf), with k the genes called in the correct direction and n_pred, n_conf the genes the prediction and the truth call significant. Over-calling costs as much as a wrong sign. | 0.5126 | 0.8065 |
| reach | `de_wilcoxon_direction_reach_raw` | Direction reach: how deep a per-target list of genes, ordered by the prediction's strength, can go while keeping 90 % direction purity. | 0.0790 | 0.9653 |
| jac | `de_wilcoxon_sig_jaccard` | Jaccard overlap between the genes the prediction calls significant and those the truth calls significant; set membership only, no direction. | 0.0303 | 0.4010 |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§1, definitions and anchors); `AGENTS.md` (§1 "Scorer facts that
decide strategy"; §5 H4); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§16.2, §16.3, §23.5);
`knowledge/synthesis/e4_official-assumptions.md` (check E2).

Three properties of the scorer decided the shape of the system (details in `05-why-it-does-not-work.md`, Part B):

1. **Sums versus cells.** PDS and MSE are functions of the 300 per-target group sums only (the comparator is
   log1p(TS · S_g / Σ_g S_g)). The four DE components come from Wilcoxon tests on the emitted cells, with fold changes
   taken from arithmetic means of CPM-normalised cells. The system therefore has a **sum side**, which fixes the
   group sums, and a **DE side**, which decides how those counts are spread over the 400 cells (§3).
2. **Only PDS discriminates targets.** On the early KOLF harness a target-shuffled copy of a prediction kept 79 % of
   the composite but only 3 % of PDS.
3. **MSE has a kink.** The noise allowance grows with the predicted amplitude up to a cap; below that point (the
   "kink") MSE falls as the prediction's overlap with the truth grows, above it MSE rises steeply. On the H1 harness
   our plain squared error was larger than that of predicting the control mean (25,410 vs 24,503), so the MSE credit
   came from the allowance and the sums had to sit just below the kink, which the test data's noise sets (about
   1.0–1.1 on A/B/C).

Sources: `AGENTS.md` (2026-10-03 "SCORER FACT"; §3); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §5).

## 2. Data used

Every response came from public perturbation screens. The final system fuses the families below; the weights are
the family weights of the fusion in F1/F2.

| Family | Public screen | Weight in F2 | Validation targets covered (of 300) |
|---|---|---|---|
| K562 | Replogle 2022, K562 genome-wide Perturb-seq (GWPS) | 8 | 272 |
| H1 | VCC 2025 H1 hESC CRISPRi (10x Flex) | 2 | 25 |
| HCT116 | X-Atlas/Orion | 1 | 299 (Orion store, both lines) |
| HEK293T | X-Atlas/Orion | 0.5 | (same store) |
| CD4 | primary CD4⁺ T-cell CRISPRi atlas (resting, 8 h and 48 h stimulation; 4 donors) | 1 | 294 |
| HEPG2, RNJURKAT | Nadig 2025, HepG2 and Jurkat | 1 each | 77, 94 |
| RPE1 | Replogle 2022, RPE1 | 1 | 95 |
| DLD1 | GSE337988, DLD-1 CRISPRi | 1 | 85 |
| FENG | Feng, HipSci iPSC lines | 0.25 | 151 |
| KOLF | KOLF2.1J iPSC genome-scale screen | 0.25, but dropped for somatic destinations (so unused for D/E/F) | 282 |

K562 carries the largest weight because, when the fusion was rebuilt on 2026-09-23, it was the best single source on
every proxy and its knockdown was stronger than that of the Orion lines (85 %, against 75 % for HCT116 and 52 % for
HEK293T). The fusion settings were screened on a fast pseudobulk evaluator (H1 and KOLF PDS and MSE, CD4 and Orion
PDS). The HEPG2, RNJURKAT, RPE1 and DLD1 families were built during validation but carried no weight in the shipped
specification until the team found, on 2026-10-06, that the local harness had been using three of them all along; they
entered with S17/S18.

Other public resources used by the final system:

- **DepMap 23Q4** CRISPR knockout fitness (1,100 models × 18,442 genes): co-essentiality neighbours for the neighbour
  blend; together with STRING, the fallback for targets that no screen measured.
- **GENCODE v47** transcription start sites: the CRISPRi promoter-neighbour prior.
- **Human Protein Atlas** cell-line expression: identity of the destination line.
- **ESM-2** gene embeddings (fixed projection to 64 dimensions): the learned S4 part's ridge for targets never measured.

Used only in local research, not in the final recipe: the Jiang screens (6 lines × 4 stimuli × 2 replicates),
Tahoe-100M (drugs), LINCS, scLong, STATE SE-600M and UCE. Tested as sources and rejected: GSE249595 (Jurkat, high-MOI),
VIPerturb-seq K562 (10x Flex) and the LINCS L1000 CRISPR-KO consensus (`04-what-did-not-work.md`, group (e)). The H1
genome-scale screen GSE295214 ("H1GW") is usable only as a same-line source for stem-like destinations, which the
final phase did not have (`03-what-raised-the-score.md`, §6).

Rules note. The handoff log of 2026-09-23 records the rules as revised on 2026-09-16 (published results may train
models but must not be submitted or used to revise predictions) and lists, as an open item, asking the organisers
whether direct mean transfer qualifies. The sources contain no reply.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§2); `AGENTS.md` (§6 "Data on disk"; 2026-09-15 "Implemented
follow-up"; 2026-09-23 "Leaderboard gap analysis"; 2026-09-24; 2026-09-29);
`vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md` (§1); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`
(§4, §18.2, §25.19); `vcc2026_data/reports/final_2026-10-27/README.md`.

## 3. The final system (F2)

### 3.1 Overview

F2 is F1 with one sum-side change (expression down-weighting kept, §3.7). F1 carried the structure of the best
validation entry S18 — the validated final-phase command, the four families added in S17, M3, and the ride-alongs
cens and drc. The main difference from S18: M2 used the agreement rank of S15 rather than the learned amplitude of
S16–S18, because the learned variant needed a two-pass build for a measured +0.0014.

For each final context the pipeline (`vcc2026_data/scripts/tv2_final.py`) identified the line from its controls and routed it as stem or
somatic (all three: somatic), rebuilt the source library for the new panel and the learned S4 part for the context,
built the sum-side change (§3.2–3.5) and the DE-side change (§3.6), emitted 400 cells per target from the context's
own control cells (§3.6), and packed and validated the file with `vcc prep`.

Sources: `AGENTS.md` (2026-09-24 "Final-phase pipeline"; 2026-10-02; 2026-10-03; 2026-10-27);
`vcc2026_data/reports/final_2026-10-27/README.md`; `vcc2026_data/reports/final_2026-10-27/run_F1_spec.sh`;
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§23.13).

### 3.2 M1 — multi-source fold-change fusion transfer

For each target, M1 takes its measured response in every public screen that perturbed it and computes a true fold
change from the raw summed counts, with a pseudocount of 50 CPM. Each source's genome-wide mean response is removed
("source centring"). Responses are averaged within a source family and then across families with the fixed weights of
§2; the fused log fold change, scaled ×1.6, is applied multiplicatively to the destination's own control cells. This
is the "fusion" part of the sum side, weighted 0.5.

**Ensembles over min_cells and pseudocounts.** Within a family every (source, target) measurement enters with equal
weight behind a hard minimum-cell cut, so a 20-cell measurement counts as much as a 500-cell one. The fusion is
therefore averaged over its default setting (min_cells 20, pseudocount 50) and 11 variants: min_cells
5/10/15/30/50/100/200 (added in S11) and pseudocounts 5/10/20/100 (`cens`, added in S17/S18). Averaging over cell
thresholds approximates weighting each measurement by an increasing function of its cell count.

Sources: `AGENTS.md` (2026-09-24 "NEW BEST OFFICIAL — S6"; 2026-09-28 S11; 2026-10-03 §25.13);
`vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md` (§2); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`
(§16.4, §25.13–25.14).

### 3.3 The learned S4 part

S4 is the learned model in the recipe and the other half (weight 0.5, ×3) of the sum side. It refits the project's
pseudobulk baseline — per-biology target means pooled with equal biology weight, a ridge on the fixed gene embeddings
for unmeasured targets, and a rank-32 covariance correction learned across biologies — on an extended table of
200,163 rows that adds CD4 donor × state contexts and the Feng targeted-screen lines, nested inside their biologies. To
this it adds 0.25 × a FiLM residual (rank 32) trained on 38 gene contexts with table-derived control profiles. S4 is
rebuilt for each destination from that context's control profile; for somatic destinations the KOLF biology (S10k) and
the FENG biology (S14fk) are excluded.

A second, small learned component sits on the sum side: a residual without gene embeddings (a target-agnostic
reshaping in a 64-PC space, trained with a discriminative loss), added at weight 0.25 (S10k).

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md` (§3 "The refit baseline"); `AGENTS.md` (2026-09-22 "NEW BEST
OFFICIAL — S4"; 2026-09-27 "Multi-angle bias review"; 2026-09-28 S10k; 2026-09-30 "NEW SUM-SIDE LEVER").

### 3.4 Neighbour blend and overrides

- **Neighbour blend.** Each target's change is mixed 0.7 / 0.3 with the energy-matched mean change of its DepMap
  co-essentiality neighbours. For targets that no screen measured, the neighbour estimate (DepMap with STRING at 0.25)
  is used in full.
- **Overrides.** The target's own knockdown and the CRISPRi promoter-neighbour prior (genes whose transcription start
  site lies within 5 kb of the target's) are written as multiplicative overrides.

Sources: `AGENTS.md` (2026-09-24 "Final-phase harness + ML arms"; 2026-09-25 "NEW BEST OFFICIAL — S8"; 2026-09-24
"NEW BEST OFFICIAL — S6"); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4).

### 3.5 Sum-side transforms, M2, M3 and the MSE budget

The combined sum-side change goes through these steps, in this order:

1. **Expression down-weighting (ew)**, kept in F2: each gene is scaled by (1 + cpm/50)^−0.25, where cpm is the gene's
   expression in the destination's controls, then renormalised.
2. **Centring 0.5**: half of the across-target mean change is subtracted (override entries excluded).
3. **Guarded library-PC removal**: 75 % of the projection on the top 3 principal axes of the source library is
   removed, scaled down towards 0 when the prediction's own energy on those axes exceeds 15 % (full stop at 25 %).
4. **M2 — per-target amplitude from cross-source agreement.** Each target's change is multiplied by
   rank(agreement)^β with β = 0.5, where agreement is the mean pairwise cosine among the source families that
   measured the target. Because the budget match (step 6) restores the total amplitude, M2 changes how the amplitude
   is shared among targets, not its total; PDS is invariant to a target's own amplitude, so M2 acts on MSE.
5. **M3 — control-cell-PC deflation of the weaker half of targets.** The top k = 8 principal axes of the destination's
   own control cells (6,000 cells, log1p CP10K, gene-centred) are computed. For the half of the targets with the
   lowest amplitude score (the score M2 uses), η = 0.75 of the change that lies along those axes is removed.
6. **MSE budget (the kink).** The kink is estimated per context from its control cells (a delete-one jackknife,
   ×1.17), and the sums are scaled to 0.8 × that estimate, but never below the recipe's own amplitude (floor 1),
   because overshooting falls off a cliff. After the transforms above, the sums are rescaled to that budget
   ("budget match"). For F2 the resulting sum scale was 1.07–1.09.

Two official probes in validation bracketed the budget: S17u (margin 1.2) and S17d (margin 0.6, floor 0) both scored
below S16, which used the shipped margin 0.8 (`03-what-raised-the-score.md`, §4).

Sources: `vcc2026_data/scripts/tv2_submit_build.py` (sum-side block; `ctrl_defl`, `tamp`, `ew`, `pcr`);
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §5, §23.8); `AGENTS.md` (2026-09-27 "S10i";
2026-10-01 "SUM-SIDE AGREEMENT AMPLITUDE"); `output/data-analysis/results/official_record.csv`;
`vcc2026_data/reports/final_2026-10-27/README.md`.

### 3.6 The DE side and the emission

**DE-side change.** The DE side uses the same recipe as the sums (`--de-follow-sum`, from S12), but built from a
source library in which the K562 screen is denoised: its responses are projected onto its 100 leading response
programmes and only 25 % of the residual is kept (`--de-denoise intended`, from S14fk/S14fk2). The sum-side
transforms (ew, centring, library-PC removal, M2, M3, budget) are not applied to it.

**DE calibration (drc 0.8 / 1.1).** Within each target, the top 5 % of genes by |DE log fold change| are scaled by
0.8 and the bottom 70 % by 1.1. The group sums are untouched, so PDS and MSE cannot change.

**Dual-moment emission.** Each target's 400 cells are emitted from the context's own control cells so that each gene's
group sum equals the sum-side prediction (γ = 1); within that fixed total, the counts of the top 3,000 genes on the
DE side are re-allocated across cells (depth_up 0.5, with an agreement ranking), which creates the per-cell shifts
the Wilcoxon tests detect. Three further steps:

- **cpc (composition pin):** every expressed gene's column is tilted by cell depth so that the group sum stays fixed
  while the mean per-cell composition — the quantity the DE fold changes use — matches the expectation of the earlier
  cell-shift emission. Without it, pinned sums made the fold changes of weak genes sign-random.
- **ci 5:** a fixed-point refinement of that pin, solved against the final rather than the original cell totals.
- **lsb:** a depth-balanced emitted block.

Emitted cells are kept below 8 × 10^5 counts, under the scorer's per-cell limit. The dual-moment emission was first
submitted after the user ruled, on 2026-09-25, that the rules do not constrain how cells are emitted.

Sources: `AGENTS.md` (2026-09-26 "NEW BEST OFFICIAL — S10b" and "— S10c"; 2026-09-27 "Codex debate loop"; 2026-09-29
S12; 2026-10-01 S14fk2; 2026-10-03 "PDS/MSE ROUND"); `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§2, §4,
§21.2–21.3).

### 3.7 F2 versus F1

F1 used the adaptive rule inherited from validation: expression down-weighting was dropped whenever the kink computed
with it fell below 1, which happened on every official context. F2 sets `ew_adaptive` to false and keeps it. Against F1
only the sums differ (cosine 0.988); the DE side, the main delta and the overrides are bit-identical.

The evidence was offline only (exact screens on the S16 base, 2026-10-28). Kept with a budget match, ew moved
OR_HCT116X PDS by +0.0139, KOLF by +0.0066 and the RN panels by +0.01 to +0.10 (J units), but was nil on H1 (+0.0002).
Stacked on the F1 deflation it added +0.012 raw PDS on OR_HCT116X and +0.0007 on H1. The README priced the remaining
stack of levers (ew, learned gene weights, a gated kNN source) at about +0.003 to +0.006 composite if official PDS
followed OR_HCT116X, as it had for S18, and about 0 if it followed H1; F2 took only ew.

Source: `vcc2026_data/reports/final_2026-10-27/README.md` (2026-10-28 sections).

## 4. Build and submission

**F1** (2026-10-27). The validated final-phase command with `--rn-sources --fam-weights DLD1=1 --sum-tamp 0.5
--sum-cens 5,10,20,100 --de-drc 0.8,1.1 --sum-defl 8,0.75,0.5`, plus the S4, neighbour-blend, dual-emission, kink,
centring, ew, library-PC, ensemble, residual, no-KOLF/FENG and denoised-DE settings listed in §3. `--same-line 0.5`
was also set, but no same-line source existed for D/E/F. Three inputs that a disk cleanup on 2026-10-11 had deleted
were restored first: a NaN stub for the old base table (the S4 rebuild reads only its features and gene list), the
learned residual's family files regenerated from the pre-rebuild source library, and the denoised K562 source,
rebuilt with denoising bit-identical to the validation version. The source library was rebuilt with the final targets
added. `vcc prep` passed (360,000 cells, targets verified). Package sha256
`4605456755d4626f71a98c9a6e00d28afadd598244b45ad5c840fb16047c90fd`. Uploaded 2026-10-27, 15:45 UTC, with the user's
go-ahead: entry `HPPRo7pAvkuTyzfik2Xl`, model name `b24-460545`, status published. Superseded by F2.

**F2** (2026-10-28). Specification = F1 with `ew_adaptive` false. `vcc prep` passed. Package sha256
`0c2a01803ee94da365cd8a39f74120b6cd64fe3bde66d212fe8c069e8a8046c8`. Uploaded 2026-10-28, 21:05 UTC, with the user's
go-ahead: entry `oYO6WAS8x4FgQt9ypgmq`, model name `b24-0c2a01`, published, partition test. **F2 is the last final
entry and the one that counts.** Nothing was uploaded after it; two further offline searches on 2026-10-29 found no
change worth an F3.

Only the opaque model name `b24-<first six hex digits of the sha256>` was sent; the leaderboard records carry no
description. Each upload left a receipt (model name, entry id, sha256, UTC time, status) next to its package in the
working repository (`vcc2026_data/reports/submission_2026-09-23/`); the same facts are in the README and `AGENTS.md`.

The pipeline had been validated before the final phase: on 2026-10-03 the final-phase command was run on the H1 panel
as if it were a new context (H1 and H1GW removed from every source) and scored with the real scorer. Its rebuilt
inputs and budget choice were sound; the one gap found (−0.030) came from feeding the rehearsal a different draw of
control cells, which led to the rule "emit from the controls the organisers ship for that context".

Sources: `vcc2026_data/reports/final_2026-10-27/README.md`; `vcc2026_data/reports/final_2026-10-27/run_F1_spec.sh`;
`AGENTS.md` (2026-10-27; 2026-10-03 "12-HOUR LOOP"; 2026-09-30 S14fk; 2026-10-31);
`output/data-analysis/results/official_record.csv`.

## 5. Final-phase score

As of 2026-11-22 the final-phase scores have not been published, and the winners have not been announced. Under the
rules the final leaderboard, based only on each team's last final entry, is released when the winners are announced
(late November 2026). Our counting entry is F2.

Sources: `output/data-analysis/results/e4_official.json` (`final_phase`); `knowledge/synthesis/e4_official-assumptions.md`.

## Sources cited

`AGENTS.md`; `RESULTS.md`; `VCC2026_gene_perturbation_report.md`; `vcc2026_data/reports/PROJECT_HISTORY.md`;
`vcc2026_data/reports/final_2026-10-27/README.md`; `vcc2026_data/reports/final_2026-10-27/run_F1_spec.sh`;
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`; `vcc2026_data/reports/transfer_v2_2026-09-23/FINDINGS.md`;
`vcc2026_data/scripts/tv2_submit_build.py`; `output/data-analysis/results/official_record.csv`;
`output/data-analysis/results/e4_official.json`; `knowledge/synthesis/e4_official-assumptions.md`.
