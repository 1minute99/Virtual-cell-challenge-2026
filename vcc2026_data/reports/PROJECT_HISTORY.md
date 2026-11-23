# AIVC / Virtual Cell Challenge 2026 — complete project history

Everything attempted in this repository, what it measured, and what it scored. Compiled
2026-09-20 from 39 report directories, 284 scripts, ~200 harness runs and 14 official
submissions. Every number is quoted from a committed artifact. Where two artifacts disagree,
both are given and the disagreement is named.

**Reading guide.** Sections 1-4 are context: the challenge, the data, the machinery, and the
measurement rules. Sections 5-12 are chronological, one per campaign. Section 13 is the
scoreboard. Section 14 collects every retraction and bug, which in this project were more often
the deliverable than the models were. Section 15 is what to do next.

---

## 1. The challenge

Given a CRISPR gene knockdown and a cell context, predict the perturbed single-cell expression
profile. The hard part is that evaluation contexts are **unseen**: the model must transfer to
cell lines it has never observed perturbed.

Scoring is `cell-eval2` 0.16.0 with the `vcc2026` preset. Six metrics, each normalised so
**0 = the context-mean baseline** (predict no change) and **1 = a split-half replicate** of the
real data. The overall score is their unweighted mean.

| metric | direction | 0 anchor | 1 anchor |
|---|---|---|---|
| `pds_cosine` (perturbation discrimination) | higher | 0.500 | 0.9515 |
| `expr_mse_unbiased_capped_norm` | lower | ~0.984 | 0.043 |
| `de_wilcoxon_lfc_nmae` | lower | 1.0008 | 0.3921 |
| `de_wilcoxon_direction_fidelity_yield_raw` | higher | 0.5126 | 0.8065 |
| `de_wilcoxon_direction_reach_raw` | higher | 0.0790 | 0.9653 |
| `de_wilcoxon_sig_jaccard` | higher | 0.0303 | 0.4010 |

Three validation contexts A/B/C plus three further final test contexts. Control cells for A/B/C
are given; their perturbation responses are withheld. Submission is 360,000 cells
(3 contexts x 300 targets x 400 cells) x 18,533 genes, as raw integer counts.

Quota: two submissions per day, one in flight at a time.

### The property of the scorer that shaped everything

`pds_cosine` is a **discrimination rank**, not a similarity. For each target it asks whether the
prediction is closer to that target's truth than to every other target's truth, scored as
`1 - rank/(n-1)` by cosine distance. Critically, `exclusion_scope="panel"` removes **every panel
target gene** from the feature space before any distance is computed, and `tie_policy="midrank"`
makes a fully tied row score exactly 0.5.

Consequences, learned repeatedly and painfully:

- A **generic all-target response** — the average thing a knockdown does — inflates ordinary row
  cosine about 6.6x while collapsing PDS.
- A **target-shuffled delta still scores 79%** of the harness composite. So the floor is not zero.
- Of the six metrics, **only PDS meaningfully discriminates targets**.
- **Never select a model on row cosine.** It and PDS routinely move in opposite directions. This
  failure recurred at least four separate times in this project.

---

## 2. The data

### The frozen table

`reports/simple_residual_2026-09-16/data/table/`, **116,217 rows x 5,592 genes**
(103,295 genetic + 12,922 drug), 205 source keys, 18,711 genetic targets, 15 genetic biologies.

| field | shape | meaning |
|---|---|---|
| `y` | (116217, 5592) float32 | log1p(CP10K) treated minus matched untreated |
| `controls` | (116217, 5592) float32 | that row's pooled untreated profile |
| `targets` | (116217,) | perturbation target |
| `bios` | (116217,) | canonical biology |
| `modality` | (116217,) int8 | 0 genetic, 1 drug |
| `features` | (116217, 64) float32 | perturbation embedding |
| `weights` | (116217,) float64 | effective count `1/(1/n_treated + 1/n_control)` |
| `sources` | (116217,) | study/replicate key |

Genetic rows per biology: CD4 35,832 · HEK293T 17,855 · HCT116 17,584 · KOLF 11,266 · K562
10,150 · JURKAT 2,201 · RPE1 2,122 · HEPG2 2,010 · FENG 1,997 · A549 424 · HAP1 424 · HT29 424 ·
MCF7 423 · BXPC3 365 · H1 218.

**Normalisation is log1p(CP10K) everywhere**, computed over measured features only. Unmeasured
genes are NaN, never zero-imputed. Minimum 25 treated cells per row.

**Perturbation features** are ESM-2 5,120-d per gene, L2-normalised, multiplied by a *fixed*
Gaussian projection (seed 20260916), L2-normalised again. No fitted statistics. Drugs use a
512-bit hashed SMILES fingerprint through the same path.

### The six effect stores (the mean-transfer baseline)

| store | biology | pairs | targets |
|---|---|---|---|
| `orion_target_cellline_all` | HEK293T 18,174 + HCT116 18,057 | 36,231 | 18,211 |
| `cd4_target_condition_all` | CD4 (Rest/Stim8hr/Stim48hr) | 36,689 | 12,346 |
| `kolf_target_batch` | KOLF2.1J iPSC (ALPHA/BETA/GAMMA) | 11,568 | 11,568 |
| `replogle_k562_gwps` | K562 genome-wide Perturb-seq | 9,866 | 9,866 |
| `feng_genomewidescreen_nonfitnessgenes_target` | FENG iPSC | 4,521 | 4,521 |
| `vcc2025_target_batch` | H1 hESC | 11,764 | 282 |

`empirical()` averages within a store by effective count, NaN-aware, then takes an **unweighted
mean across stores that measured the gene**. Equal per-source weighting is deliberate: a response
reproducing across K562 and primary CD4 T cells is better transfer evidence than one seen only in
the source with the most cells.

Additional table-only sources: `rn/*` (Replogle K562/RPE1, Nadig HepG2/Jurkat),
`jiang/{TGFB,IFNB,TNFA,INS}` (6 lines x 2 replicates), `tahoe/plate{3,6,14}` (drugs).

### External datasets acquired

Jiang (1,383,236 cells, 4 stimuli), Replogle+Nadig scPerturb (13.7 GB), LINCS GSE106127
(8 lines x 200 shRNA), CCLE/DepMap 23Q4 expression (449,782,440 bytes) and dependency
(1,100 models x 18,442 genes), Tahoe-100M (157 shards, 16.3 GB, 16 lines x 150 drugs, 379
distinct drugs across 1,135 drug-dose conditions), MSigDB Hallmark 2025.1.Hs, DoRothEA 0.99.0,
STATE SE-600M, scLong, UCE.

sci-Plex3 planned, never acquired. Moonen CD4 embargoed to 2027-03-05, past the challenge.

---

## 3. The machinery

### The refit baseline (`SourceBaseline`) — the best local predictor

Five steps, in `simple_residual_data.py:112`:

1. **Pool duplicates** within (biology, target) by effective count.
2. **Equal-biology target means.** Each biology contributes equally regardless of cell count.
3. **Ridge from ESM features** for unseen targets: `lambda = 10x` the mean eigenvalue scale,
   intercept essentially unpenalised.
4. **Rank-32 covariance denoiser.** Per-biology centre, then for each target compute the
   *leave-one-biology-out* mean of the others, regress the residual `specific - others` on a
   rank-32 basis of the centred target means with equal biology weight.
5. **Predict:** ridge for all, empirical mean overriding where the target was seen; centre across
   the query set; add `(out @ b) @ coef`; L1 amplitude-match back to the raw norm; unseen targets
   keep the pure ridge anchor.

`training_anchors()` refits the entire baseline once per biology, excluding that biology's own
rows, so every training row's anchor is leakage-free at the biology level.

**Critical, and a trap that wasted a whole experiment:** `predict` overwrites the feature ridge
with each target's empirical mean for every target seen in training, and all 300 official and 279
KOLF targets are seen in every fold. So permuting or zeroing the gene features changes the refit
prediction by **exactly 0.000**.

The deployed composite is: the 5,592 modelled columns from `SourceBaseline`, the remaining 12,941
official columns from the empirical multi-source mean.

### Count reconstruction

A log value is never added to a count. Per cell:

```
ctrl_log[g] = log1p(n[g] / L * 1e4)
pred_log[g] = ctrl_log[g] + delta[g]
counts[g]   = round(expm1(pred_log[g]) * L / 1e4)
```

implemented fused as `want = ((n/L*1e4 + 1) * exp(delta) - 1) * L/1e4`. Where `delta == 0` this
returns exactly `n[g]`, so the no-change fallback is bit-exact.

**Stochastic rounding is essential.** Round-to-nearest transmits only 22% of the intended delta
(0.0034 of 0.0152, correlation 0.78) because most genes sit at 0-3 counts and a small delta never
crosses a 0.5 boundary. Dithered rounding sends 0.63 to 1 with probability 0.63, so the 400-cell
group mean carries the delta exactly: recovery 0.0153 of 0.0152, correlation 0.995.

**One shared cell sample per context**, reused by every target. An independent per-target draw
leaks sampling noise into the between-target differences PDS scores: agreement between realised
and intended target-to-target distance falls 0.9987 to 0.9868.

### The harness

Five held-out lines: K562, RPE1, HepG2, Jurkat, KOLF. Truth from `rn/<line>` for the four RN
lines and `xl/kolf:*` for KOLF. **Nested** leave-one-biology-out: for outer line O and inner line
I, training excludes *both*, the model predicts I, selection ranks by mean inner PDS, and only
then is O scored. The outer line never participates in selection.

Frozen promotion gate, used throughout:
**gain >= 0.001, >= 3/5 line wins, >= 9/15 paired-seed wins, all valid.**

Seeds 0/1/2, 400 cells, scale 3, 32-fold random-derangement null. Seed noise floor is **±0.010**
composite; treat anything smaller as undecided.

### The delta scale

Submissions use **x3**. The justification in the code is that the transferred delta is too small
and every attenuation lost, with the claim that "pds does not move, because the discrimination
metric is cosine and cosine is scale-invariant."

**That claim is false, and the KOLF leave-one-out bundle contradicts the choice:**

| scale | composite | PDS | fidelity | jaccard |
|---|---|---|---|---|
| x1 | 0.1385 | 0.1198 | 0.4517 | 0.1085 |
| **x2** | **0.1972** | 0.1570 | 0.7653 | 0.0878 |
| x3 | 0.1944 | 0.1656 | 0.7876 | 0.0579 |
| x4 | 0.1867 | 0.1759 | 0.7907 | 0.0437 |

PDS rises 47% from x1 to x4 because the delta passes through `expm1`, a clip at zero and
stochastic rounding before pseudobulking, all nonlinear. Scale invariance holds for the metric on
a delta, not for the delta's route through integer counts. The composite **peaks at x2**; Jaccard
collapses faster than PDS gains. This is a live open question, tempered by the fact that KOLF was
demoted from a promotion gate after it mispredicted b15's direction (section 6).

---

## 4. Measurement rules the project adopted

From `AGENTS.md`, after they were learned the hard way:

1. **Judge on `pds_cosine` only. Never gate on `avg_score`.**
2. Always score a target-shuffled copy of the same delta; report `(pds - shuffled)/(baseline - shuffled)`.
3. Four seeds of an identical delta span composite 0.1972-0.2085, sd 0.0076. Compare only paired
   on seed; treat anything under 0.01 composite as undecided.
4. Harness PDS (0.157) and leaderboard PDS (0.507) live on different scales. Never port thresholds.
5. Historical scientific sources, checkpoints and submissions are **read-only**; verification is
   by frozen hash because the repository has no usable git history.
6. Official uploads require explicit user authorisation. Credentials live in the system keyring,
   never in source, chat or logs.
7. Long jobs run detached via `scripts/run_detached.sh` under a systemd slice capped at 40 GiB
   soft / 48 GiB hard.

---

## 5. Era 1 — foundation models and single-cell training (Sept 5-7)

### scLong, 1.01B parameters, frozen

Two uses, both negative.

**Zero-shot in-silico knockout: no signal.** Masking a target gene's input token on real control
cells moves essentially only that gene's own output slot. All other genes move <= 0.03 log units;
the same 3-4 genes shift for every target; `mean|delta|` over 27,875 genes is **0.0017**;
reconstruction-vs-input Pearson ~0.50. The cause is structural: scLong is a masked-expression
autoencoder with no perturbation head. The paper's perturbation results come from using it as a
frozen encoder under a *trained* GEARS head.

**As a context encoder: worst of three.** Correlating inter-context distance against transfer
quality over 6 context pairs:

| representation | Pearson | Spearman | p |
|---|---|---|---|
| PCA, 50 components | -0.618 | **-0.943** | 0.005 |
| STATE SE, 2,058-d | -0.668 | -0.886 | 0.019 |
| scLong, 200-d | **+0.589** | **+0.029** | 0.957 |

scLong's inter-context cosine distances span 0.006-0.016, barely separating cell lines. Its
nearest-source assignment says HepG2 for all three VCC contexts; PCA and STATE SE agree exactly
with each other (A to Jurkat, B to RPE1, C to HepG2).

**Note the ordering: plain PCA beat both foundation models.** That result recurs on Sept 20.

Throughput measured: 8.08 GB VRAM, 5.21 cells/s on an RTX 5090, compute-bound and flat in batch
size. Caching 5,200 cells took 16.8 minutes and produced a 12.9 GB array.

### scLong + trained GEARS head

34.5M-parameter head on the cached frozen encoder, 165,393 training cells, 585 validation groups.
Baseline (predict control) MSE 0.02929. Epoch 0 gave val unseen-perturbation delta-Pearson
**0.5168** but MSE **0.06336**, i.e. more than twice the control baseline's error. The log stops
mid-epoch 1 with no further validation. **Abandoned; no prediction was ever produced.**

### The two first official submissions

| method | what | overall | rank |
|---|---|---|---|
| A | Training-free Replogle transfer, deterministic rounding | **-0.0228** | 468 |
| B | STATE ST trained from scratch, 4 lines + ESM2 features | **-0.1910** | 609 |

Both are worse than predicting no change. **The training-free arithmetic method beat the trained
45.6M-parameter neural model on all six metrics**, by +0.168 overall and +0.279 on PDS.

Offline diagnostics predicted this correctly before submission. The STATE model's delta-over-null
ratio was only 1.30/1.88/1.49 across contexts A/B/C, and a naive "fraction target-specific"
statistic read 96% for the model while pure control-resampling noise scored 99.8%.

Method B also exposed three upstream traps, each documented with file and line: `state_sm` ships
`gene_decoder_bool: False` causing a silent 3x coverage loss; `state tx infer --tsv` is O(n²) and
produced nothing in 11 minutes where a vectorised pre-pad took 3 seconds; and **STATE ST silently
substitutes the control vector for any perturbation not in its mapping**, so without
`perturbation_features_file` it predicts control and reports nothing.

### Cross-cell-line transfer pilot (K562 to RPE1)

100 perturbations, 2,000 genes, macro MSE:

| method | macro MSE | delta Pearson |
|---|---|---|
| no_change | 0.0162 | — |
| mean_effect | 0.0147 | 0.306 |
| matched_effect | **0.0143** | 0.247 |
| **state_st (35.2M params)** | **0.0350** | **-0.076** |

STATE ST lost to no-change on **100 of 100** perturbations, paired bootstrap CI
[+0.0176, +0.0203]. Diagnostics: its predicted delta L2 was 1.24x too large, only **0.48%** of its
prediction variance was perturbation-specific versus 2.74% in real data, and its source validation
loss bottomed at step 1,999 of 20,000 then rose monotonically. The reverse direction was worse:
RPE1 to K562 transfer actively hurt, better on only 5 of 100.

A noise ceiling was established here: split-half reproducibility median 0.300 in K562 and 0.663 in
RPE1, giving `sqrt(rel_source * rel_target)` = 0.384 against an achieved 0.194.

### The 24-experiment transfer grid

Source-copy beat both floors on delta-Pearson, direction match, DE overlap and discrimination in
**6 of 6** source-target pairs, while **losing to return-control on MAE in all 6**. Pooled over
7,651 scored perturbations only **0.21%** of source-copy predictions reach delta-Pearson > 0.7;
both floors reach it 0.00%. Best pair was RPE1 to HepG2 (0.2886), worst Jurkat to HepG2 (0.2263).

The lesson that set the project's direction: **transfer is real but weak, and arithmetic beats
learning.**

---

## 6. Era 2 — mean transfer, the scale finding, and the first real score (Sept 14)

Submissions b1 through b8 explored source aggregation and count rounding:

| sub | what | overall | rank |
|---|---|---|---|
| B1 | CD4 alone, deterministic rounding | -0.0064 | 484 |
| B2 | 3-source mean, deterministic | +0.0020 | 467 |
| B6 | 6-source mean + stochastic rounding | +0.0771 | 342 |
| B7 | sign-agreement weighting + stochastic | -0.0224 | 550 |
| **B8** | **6-source mean + stochastic + delta x3** | **+0.1217** | **214** |
| B10 | plain mean transfer at x4 | +0.0868 | 328 |
| B11 | expression-gated per-context, x2.25 | +0.1156 | 234 |
| B14 | precision-weighted sources, x3.5 | +0.0888 | 345 |

Two findings dominate. **Stochastic rounding is worth +0.075** (b2 to b6). **Delta x3 is worth
+0.045** (b6 to b8). Together they took the project from roughly zero to 0.1217 without any model
at all.

b7 established that sign-agreement weighting hurts. b11 is interesting: an expression-gated
per-context delta with only **eight fitted numbers** (gene quantile bins with multipliers rising
0.55 to 2.21) improved the local harness composite by +0.0085 with 4/4 paired seeds, and was the
first working context-specific method — but it scored **0.1156, a loss** against b8's 0.1217 on
the leaderboard.

---

## 7. Era 3 — the reassessment that demoted KOLF (Sept 14-15)

### The Jiang external benchmark

Holding out an entire cell line on a 6-line, 4-stimulus panel:

| feature set | mean transfer | centred mean | nested learned | shuffle floor |
|---|---|---|---|---|
| Full measured features | 0.74078 | 0.76885 | **0.78356** | 0.49600 |
| Intersected with official genes | 0.73674 | 0.76485 | **0.77825** | 0.49592 |

Learned wins **6/6 cell lines** and 22/23 condition-line folds over mean transfer, gain +0.0415.
The only negative fold is INS/MCF7 at -0.0017. This looked like proof that a learned operator
beats mean transfer across biologies.

### The KOLF result, and the control that reinterpreted it

On the KOLF reconstructed-cell harness the learned correction gained decisively:

| arm | normalised PDS | composite |
|---|---|---|
| mean transfer | 0.1570 | 0.1972 |
| crsres128 (learned on mean) | **0.2359** | 0.2092 |
| crsres128prec (learned on precision) | **0.2856** | **0.2170** |

4/4 paired seeds for both. **Then the training-label shuffle:** `crsres128_trainshuffle` scored
**0.2449**, *above* the matched-label `crsres128` at 0.2359. The rank-based variant showed the
same inversion (0.1699 shuffled vs 0.1635 matched).

The algebra explains it. For centred sources `S_i = M + E_i` with `sum E_i = 0`, the pooled
residual regression learns `-c/(c-1)^2 sum_i E_i^T E_i` — a **source-disagreement covariance**,
i.e. covariance denoising. So the gain is real but it does **not require correctly matched
training labels**, and therefore is not evidence of learned response prediction for a new context.

### b15 fails and KOLF is demoted

b15 applied effective-cell-count source weighting, which KOLF had predicted would help (raw
precision improved KOLF PDS 0.1570 to 0.1957). On the leaderboard it scored **0.1067 against
b8's 0.1217, a 12.3% loss**, with PDS falling 0.0455.

Recorded decision: **"KOLF must not be used as the sole gate for leaderboard promotion."** It
remains useful for numerical checks, target ordering, permutation controls and ablations.

A naming correction was also recorded: the weighting called "precision" is
`1/(1/n_pert + 1/n_control)`, an **effective cell-count proxy**, not measured gene-wise precision.

### b16 succeeds — the frozen covariance denoiser

Rank-128 signed correction over 17,742 coordinates, trained on 1,661 non-official targets from
HCT116, HEK293T, CD4 and Feng iPSC, applied to the equal-source mean at x3.

**0.12919 against b8's 0.12171, +6.15%.** PDS 0.5037 to **0.5547** (+10.1% relative), which
contributes +0.0085 of the +0.0075 overall. Rank 198.

The conclusion was deliberately qualified: "the local training-label-shuffled residual model
performed similarly, so the supported mechanism remains **covariance denoising**; useful
matched-label supervision and a learned response function for new contexts are not established."

---

## 8. Era 4 — the context campaign (Sept 15)

The largest single campaign. Question: does the untreated state of a cell predict how its
perturbation response deviates?

### Eight architectures on Jiang, all failed

184 outer model-folds, 3,300 cached fits. Against b16's 0.779015:

| model | params | correct-context PDS | vs b16 | vs matched no-context |
|---|---|---|---|---|
| Pairwise DeltaNet | 440-832 | 0.779282 | +0.000267 | **-0.000787** |
| Learned metric | 44-76 | 0.779259 | +0.000244 | +0.000244 |
| Shared-kernel GP | 3 | 0.779236 | +0.000221 | +0.000221 |
| Low-rank hypernetwork | 124-188 | 0.778957 | -0.000059 | -0.001317 |
| Conditional neural process | 1152-1680 | 0.778654 | -0.000361 | -0.000391 |
| Mixture of experts | 201-777 | 0.778266 | -0.000749 | -0.001521 |
| Low-rank tensor | 72-136 | 0.778207 | -0.000808 | -0.002108 |
| Conditional DeepSets | 1016-1408 | 0.777513 | -0.001502 | -0.001844 |

The best gain, +0.000267, is **smaller than its own 0.000470 seed standard deviation**. Every
Holm-adjusted p-value against b16 is 1.0. The GP's *label-shuffled* training scored 0.779511,
above the correctly trained GP. The metric model's *zero-query* control beat its correct query.

The best context-free number in that table, the MoE's selected-no-context 0.782063, then lost on
both broader panels when tested independently (-0.0031 on four biologies, -0.0003 on RN).

### The context-interaction model

Regressing the residual on `[z, z x h]` scored **0.77835 against the denoiser's 0.77902**,
winning 3/6 lines. Correct context beat average wrong context in 5/6 lines but **lost to the
zero-context ablation in 5/6**. Shuffling training context assignments changed average PDS by
~0.00006.

### External data: LINCS, Tahoe, DepMap

**LINCS** (8 lines x 200 shRNA targets, CCLE baseline expression as context):

| model | PDS | correct - no-context | line wins |
|---|---|---|---|
| Source-only denoiser | **0.953222** | — | — |
| Additive ridge | 0.937073 | +0.000967 | 5/8 |
| Bilinear | 0.935176 | -0.000930 | 4/8 |
| Mean transfer | 0.920107 | — | — |

Every context model loses to the context-free denoiser on **every one of the 8 lines**. The
decisive debunk: the ridge's correct-minus-no-context prediction is *exactly a target-independent
vector in every line*. Removing the across-target mean gives **identical** centred PDS (0.951743)
and identical deviation R² (0.014769) for correct and no-context. The apparent gain encodes no
perturbation-by-destination interaction at all.

**Tahoe** (16 lines x 150 drugs, 2.1M cells, 8 model families, 128 outer folds): the largest
increment across all eight families was **0.0000364 PDS**. **All eight predefined gates failed**;
every paired-line interval included zero.

**DepMap** was the one success, and for a different endpoint. Predicting CRISPR knockout
*fitness* from untreated expression gave ridge context R² **0.045245**, a 4.52% error reduction,
bootstrap [0.001024, 0.001224], passing its gates. A synthetic positive control recovered an
implanted interaction at R² > 0.99. But this is fitness, not RNA response.

### The variance decomposition that framed everything after

Reproducible cell-specific variation **exists** and is substantial:

| dataset | shared | cell-specific repeatable | noise/disagreement |
|---|---|---|---|
| Jiang IFNB | 17.7% | **24.0%** | 58.4% |
| Jiang INS | 8.0% | **30.1%** | 61.9% |
| Jiang TGFB | 6.8% | **20.9%** | 72.3% |
| Jiang TNFA | 9.4% | **27.2%** | 63.4% |
| Tahoe (91 repeats) | 43.5% | **22.7%** | 33.8% |
| LINCS technical | 15.0% | **16.4%** | 68.6% |
| LINCS distinct reagent | 2.5% | **8.0%** | 89.5% |

So 8-30% of the response is genuinely cell-specific and repeatable. **Almost none of it was
recovered by any mapping from untreated expression.** The one exception, added in a later
supplement, was a Tahoe full-GP recovery of **1.445%** of repeatable contrast — real, small, and
it did not pass the frozen gate.

A STATE control diagnostic sharpened the point: cell line identification from STATE embeddings is
**100% accurate** (within-line replicate L2 ~0.20-0.27 vs between-line ~0.66-0.69). Cell identity
is trivially distinguishable. That does not make response predictable.

### PerturbedVAE

Seven model variants on Tahoe. The winner was **not** the VAE:

| model | repeat recovery | correct - no-context | improving lines |
|---|---|---|---|
| **Deterministic anchored residual** | **+2.5592%** | **+2.5591%** | **14/16** |
| Anchored residual VAE | +0.0382% | +0.4237% | 11/16 |
| Context-conditioned delta VAE | -6.2100% | +0.4255% | 13/16 |
| Source-faithful full VAE | -12.9501% | -5.3815% | 0/16 |

A conservative deterministic anchored residual beat VAE disentanglement by 2.52 percentage
points, 15/16 lines. Its *unshrunk* recovery is **-1.8740%**, so the shrinkage is the whole
method. Drug-to-gene transfer was not established on either genetic benchmark; frozen transfer
even lost to frozen-random on Jiang (-0.2627 pp).

### Correction strength: the contradiction that stayed unresolved

| endpoint | x1.0 | x1.5 | verdict |
|---|---|---|---|
| KOLF reconstructed cells, 4 seeds | 0.230554 | **0.257194** | 1.5x wins 4/4 |
| Four broader biologies, CP10K | **0.610668** | 0.569491 | 1.5x loses 4/4 |

Resolved against 1.5x on the four-biology evidence. This is the same KOLF-versus-everything-else
split that b15 exposed.

---

## 9. Era 5 — the b19 lineage (Sept 16)

### The refit baseline arrives

`simple_residual_2026-09-16` asked whether adding Tahoe drug data helps a small context residual.
Four-line RN panel, three seeds:

| candidate | PDS | vs refit |
|---|---|---|
| No change | 0.500000 | -0.402890 |
| Historical-source mean | 0.763471 | -0.139419 |
| Archived b16 | 0.854647 | -0.048243 |
| **Source-refitted rank-32 baseline, alpha 0** | **0.902890** | 0 |
| Genetic residual | 0.897102 | -0.005788 |
| Joint genetic + drug | 0.896723 | -0.006167 |
| Drug pretrain then genetic finetune | 0.900244 | -0.002646 |

**No residual beat the alpha-0 baseline.** Drug addition contributed +0.001067 and drug
pretraining -0.000267, neither enough. The refitted baseline itself was the finding.

Then the KOLF check reversed the deployment choice: refit **0.5658** versus legacy b16
**0.5961**, losing 3/3 paired seeds. The frozen pick became archived b16 with alpha 0. That
shipped as **b18**, which scored 0.129193042823 — a difference of **-5.55e-17** from b16, i.e.
the package was **byte-identical**. A wasted submission slot.

### Context forcing: did the residuals learn nothing?

They learned; it did not transfer. All 60 selected checkpoints have nonzero corrections with
Frobenius norms 14.2-45.1% of baseline. Fitted latent residual MSE is 8.9-29.9% below the
zero-residual error.

A decisive algebra finding: **the linear family's context term cancels under query
target-centring**. Replay confirmed context-swap changes of 1-2e-7. So K562's apparent +0.001193
gain in that family could not possibly be context learning. That retracted an earlier reading.

Direct prediction with no empirical mean at all was tested and is far worse: **0.574450 versus
0.786234** on identical gene support, losing **15/15** paired comparisons.

Context controls on the five-line mean:

| treatment | residual PDS |
|---|---|
| Correct query context | 0.795414 |
| Separately trained no-context | 0.793970 |
| **Same model, zero query context** | **0.799576** |
| Same model, swapped query context | 0.785168 |
| Separately trained, shuffled context | 0.742342 |

**Zero context beat correct context.** Correct context won only 2/5 lines.

### The strength sweep, and b19 ships

Multiplying the learned residual:

| alpha | K562 | RPE1 | HepG2 | Jurkat | KOLF | five-line |
|---|---|---|---|---|---|---|
| **0** | 0.882070 | 0.810684 | 0.884560 | 0.859588 | 0.565788 | **0.800538** |
| 1 | 0.883530 | 0.805006 | 0.856429 | 0.859666 | 0.568947 | 0.794715 |
| 1.5 | 0.881319 | 0.799706 | 0.840762 | 0.857202 | 0.571465 | 0.790091 |
| 2 | 0.877521 | 0.794121 | 0.827347 | 0.852962 | 0.574495 | 0.785289 |
| 3 | 0.863831 | 0.783162 | — | 0.838038 | 0.582167 | — |
| 4 | 0.844378 | 0.772843 | — | 0.814589 | 0.584307 | — |

**Monotone degradation on four lines; monotone improvement only on KOLF.** Eight outputs at high
alpha were invalid or unscorable because reconstructed cells exceeded the 1,000,000-count cap or
float32's exact integer range.

Since alpha 0 was deliberately ineligible for a "nonzero residual" submission, alpha 1 shipped as
**b19**:

| metric | b16 | b19 | difference |
|---|---|---|---|
| **Overall** | 0.129193 | **0.136470** | **+0.007277 (+5.63%)** |
| PDS | 0.554671 | **0.543156** | **-0.011514** |
| LFC NMAE | 0.030130 | **0.075024** | **+0.044893** |
| fidelity | 0.010708 | 0.018797 | +0.008089 |
| reach | 0.187440 | 0.188645 | +0.001205 |
| jaccard | -0.007791 | -0.006802 | +0.000988 |

**b19 is the project's best official score and its PDS is worse than b16's.** The entire gain is
LFC-NMAE. The report says plainly that b19 also changed the baseline rank and training population,
so the official comparison cannot isolate the residual's causal effect, and the controlled local
experiment (alpha 0 at 0.800538 versus alpha 1 at 0.794715) says the residual *hurts*.

### b19 optimisation: does any of it matter?

Nested CV, 169 fits, 209 checkpoints.

**Longer training makes the source fit better and transfer worse:**

| schedule | updates | source latent MSE | inner PDS |
|---|---|---|---|
| cosine | 2,500 | 0.416515 | **0.772073** |
| cosine | 10,000 | 0.365971 | 0.770536 |
| cosine | 25,000 | 0.344988 | 0.767051 |

Three of five outer selections chose 2,500 updates; **none chose 25,000**.

**Capacity did not help.** All five selections kept rank 32 / hidden 64 (24,192 parameters).
Rank 64 (30,336 params) passed none of 20 directed gates, so rank 128 was never run. Hidden 128
(72,448 params) lost.

**Gate failed:** gain +0.000385, 2/5 lines, 8/15 paired.

Null-adjusted excess relative to b19: refit **1.006815**, b19 1.000000, selected 0.997363. The
baseline beats the model on this measure too.

---

## 10. Era 6 — sweeps, the data audit, and infrastructure (Sept 18)

### Learning-rate sweep: the grid is flat

30 configurations, 5 learning rates x 2 schedules x 3 budgets.

| by learning rate | inner PDS | by budget | inner PDS |
|---|---|---|---|
| 1e-4 | 0.770606 | 2,500 | **0.770505** |
| 3e-4 | 0.770358 | 10,000 | 0.769519 |
| 1e-3 | 0.769819 | 40,000 | 0.766324 |
| 3e-3 | 0.768754 | | |
| 1e-2 | 0.764377 | | |

Total spread across learning rates 0.0062; across schedules **0.0002**. **More steps is worse.**
Gate failed: -0.001689, 1/5 lines, 5/15 paired.

The context controls again: `zero_query` scored **0.799013**, above the correct-context reference
at 0.795777. Only `context_shuffle` (0.751128) was clearly destructive.

### Sweep 2: the one gate that passed

Arm A swept weight decay x batch size at a fixed 640,000-sample budget. Total spread across the
whole 3x3 grid: **0.0028**. Gate failed.

Arm B tested seed ensembling:

| K | five-line PDS | gate |
|---|---|---|
| 1 | 0.795777 | — |
| **3** | **0.796985** | **PASSED** (+0.0012, 4/5 lines, 10/15 paired) |
| 5 | 0.795818 | failed (+0.00004) |

**This is the only promotion gate that passed in the entire project.** Its gain of +0.0012 is
smaller than the measured per-line seed standard deviation on 4 of 5 lines, and it is negative on
KOLF. It was not deployed.

Reproducibility check: 45 of 45 comparable reference fits bit-exact against sweep 1.

### The data audit: the gene-axis defect

The table assembler intersected **every** source axis, including evaluation panels and the unused
Tahoe drug cache:

| axis | official genes |
|---|---|
| Training sources jointly | **17,847** |
| RN evaluation panels | 8,248 / 8,748 / 8,881 / 9,624 |
| **Resulting training axis** | **5,592** |

The genetic model never uses Tahoe's 12,922 drug rows, yet Tahoe constrained its gene axis.
Official contexts measure all **18,533** features, so **~70% of the scored readout got no learned
component**.

Important scoping correction: **mean transfer was never affected**, because `empirical()` reads
per-source stores on the full official axis with proper masking. Measured identical in both
pipelines, difference 0.0000. Only the refit ridge and the neural residual were restricted.

The narrow-arm control, restricting the same rows back to 5,592 genes, quantified what the extra
biologies are worth:

| line | refit frozen | refit narrow | b19 frozen | b19 narrow |
|---|---|---|---|---|
| K562 | 0.8821 | 0.7977 | 0.8835 | 0.7767 |
| RPE1 | 0.8107 | 0.7465 | 0.8103 | 0.7404 |
| HepG2 | 0.8846 | 0.7888 | 0.8564 | 0.7777 |
| Jurkat | 0.8596 | 0.7890 | 0.8597 | 0.7756 |
| **KOLF** | 0.5658 | **0.5683** | 0.5689 | **0.5712** |
| five-line | **0.8005** | 0.7381 | 0.7958 | 0.7283 |

Dropping the RN and Jiang biologies, 9% of rows but three extra biologies, costs **-0.062 on the
baseline and -0.068 on the model**. **Biologies matter roughly 10x more than genes.** And KOLF is
the exception and the warning: it is the one line where losing them costs nothing, and in fact
gains.

**The fix was built and never deployed.** `table_fullaxis` is 84,752 x 17,282 with bit-identical
values on all 5,592 shared genes. `table_fullaxis_masked` is 103,295 x 17,282 with 9.1% NaN, plus
a mask-aware layer that imputes at target mean for the baseline path and fits each row's rank-32
latent by least squares over measured genes only.

**And it worked.** `masked_baseline/summary.json` records the wide-axis refit at
**0.8015888835448705** against the frozen refit's 0.8005378 — **+0.00105**, with KOLF declining
to 0.5638 from 0.5658. That number is not mentioned in the audit's own findings file, which still
says the masked baseline never produced a result.

### Infrastructure

A VS Code main-process crash (V8 4 GB heap limit, two final GCs each retaining 3,989 MB of a
4,000 MB heap) motivated `scripts/run_detached.sh`: user-systemd supervision independent of the
editor, file logs, duplicate-name protection, and a shared slice capped at 40 GiB soft / 48 GiB
hard. W&B logging was added under a bit-exactness requirement: logged training must be
bit-identical to untraced training, verified by test.

---

## 11. Era 7 — diagnosis (Sept 19)

### Is the data good enough?

Not a preprocessing question but a signal question. All figures are target-specific, each source
centred over its own targets first so the shared generic response cannot inflate them.

| comparison | n targets | r |
|---|---|---|
| CD4 Rest vs Stim8hr (same cells, different state) | 11,873 | **0.270** |
| CD4 Rest vs Stim48hr | 11,817 | 0.187 |
| RN HepG2 vs Jurkat (same source) | 1,901 | 0.174 |
| RN K562 vs RPE1 (same source) | 2,089 | 0.163 |
| Orion HCT116 vs HEK293T (same source) | 17,498 | 0.063 |
| Orion HCT116 vs RN K562 (cross source) | 9,123 | 0.071 |
| CD4 Rest vs RN K562 (cross source) | 8,456 | 0.030 |
| Orion vs Feng / CD4 / KOLF | 1,795-11,434 | 0.013-0.017 |

**Reproducibility ceiling: Jiang Rep1 vs Rep2, same line, same stimulus, median r = 0.209.** About
80% of even that signal is not reproducible.

**Only 2.4% of rows have any replicate.** KOLF's three batches do not qualify: ALPHA/BETA/GAMMA
screen disjoint gene sets, sharing **0, 0 and 1** targets. So for 97.6% of the table, noise cannot
be separated from biology.

Verdict: the data is correctly processed and not short of rows (103,295 genetic rows, 18,711
targets), but close to exhausted as a cross-biology signal source. Transferable component
r ~ 0.02-0.07 against a reproducibility ceiling of 0.21.

### ESM embeddings are a lookup key

16 arms, against b19's 0.795778:

| arm | five-line | vs b19 |
|---|---|---|
| **zero_64 (no gene features at all)** | **0.797836** | **+0.0021** |
| esmc_600m_64 | 0.797351 | +0.0016 |
| esmc_600m_32 | 0.796874 | +0.0011 |
| **random_64 (Gaussian noise)** | 0.796778 | +0.0010 |
| esm2_15b_64 (= frozen b19) | 0.795778 | +0.0000 |
| esm2_15b_128 | 0.794033 | -0.0017 |
| esmc_600m_256 | 0.793390 | -0.0024 |

**Deleting the gene embeddings entirely is the best arm.** Random Gaussian vectors of the same
shape beat a 15-billion-parameter protein language model. There is no dose-response with model
scale: ESM-C 600M is both the second-best and the worst arm depending only on width.

Two mechanisms were measured. First, the gene branch treats its input as a **hash**: the trained
first layer puts comparable weight on embedding and baseline columns (ratios 0.80-1.13), so the
embeddings *are* read — what is used is gene distinctness, not biology. Second, and more
important, **the context branch is a 58-value source lookup and it dominates**: context branch RMS
0.308-0.332 versus gene branch 0.154-0.161, so context supplies 67% of the residual. On the KOLF
fold there are only **58 distinct control profiles across 92,029 rows, one per source**. The table
stores one pooled untreated profile per context and broadcasts it, so the "context input" is a
categorical source label, not a continuous covariate.

That single fact explains every context null in the project up to this point.

A caveat on the record: this experiment's own findings file tabulates only 13 of the 16 arms and
omits `zero_64`, the best one. The full set is in `summary.json`.

### Gene subsets: more genes is simply better

Ranking genes by between-target variance x cross-biology ICC, against the all-genes refit:

| K | ranked | vs all 5,592 | ranked minus size-matched random |
|---|---|---|---|
| 500 | 0.7711 | -0.0294 | **+0.0241** |
| 1,000 | 0.7813 | -0.0192 | **+0.0241** |
| 2,000 | 0.7905 | -0.0101 | +0.0163 |
| 3,000 | 0.7955 | -0.0050 | +0.0095 |
| 5,000 | 0.8000 | -0.0006 | — |
| 5,592 | **0.8005** | 0 | — |

**No subset improves at any size**, and the deficit shrinks monotonically toward zero as K
approaches the full set. But the ranking itself works: given that you must pick 500 genes, the
reliability statistic picks much better ones than chance, +0.024.

The exception is diagnostic: the ranking **does not help on KOLF** (-0.0011 at K=500) while
helping strongly on K562, HepG2 and Jurkat. Cross-biology reliability estimated from source
biologies transfers to lines resembling the sources, and not to the one that does not.

### The context model on refit

`prediction = refit + alpha * residual(context, perturbation)`, gated neural net, two encoders,
LOBO, early stopping on a held-out training biology.

**Best configuration 0.8008 versus 0.8005, +0.0003**, and performance fell monotonically with
alpha, so the optimum is alpha = 0. Substituting another line's context cost 0.042, which is
large — but a **gene-shuffled context scored better than the correct one**, so most of the
ablation spread tracks residual magnitude rather than context accuracy.

The enabling data fix, found here: Replogle-Nadig effects had been computed against per-target
batch-weighted controls but paired by the assembler with the line's **pooled** control (median L2
error 0.646 over 9,624 genes). Restoring the matched control raised distinct contexts from
**205 to 16,263**.

---

## 12. Era 8 — the oracle, and the last word on context (Sept 20)

This session. Five linked experiments that first revived context, then killed it.

### 12.1 Pooled basal state (`basal_context_2026-09-20`)

The previous context analysis had used per-target matched-control fluctuations within a line.
Those are roughly 60% sampling noise. Redone with context defined as the **robust pooled control
profile of a whole biology**: 9 biologies x 11,058 shared perturbations x 5,592 genes.

Decomposition `effect[i,p] = mu[p] + g[i] + s[i,p]`:

**94% of the deviation is perturbation-specific** (mean generic share 0.058). So the half context
might plausibly explain is the small half.

**Only ~28% of the deviation is reproducible at all**, by fitting
`||delta_p||^2 = a + b/effective_n_p` and taking the intercept:

| biology | signal fraction | biology | signal fraction |
|---|---|---|---|
| HepG2 | 0.681 | K562 | 0.180 |
| CD4 | 0.565 | RPE1 | 0.165 |
| Jurkat | 0.344 | HCT116 | 0.163 |
| FENG | 0.199 | **KOLF** | **0.007** |
| HEK293T | 0.192 | **mean** | **0.277** |

**KOLF's deviation is 99.3% noise.** That retrospectively explains the unexplained "+0.019 on
KOLF" anomaly from the earlier modelling run: there is essentially nothing reproducible there to
predict. It also explains why KOLF behaves opposite to every other line throughout this project.

### 12.2 The oracle test — context IS worth something

The upper-bound question: if you already knew a cell line's perturbation effect from an
**independent batch**, how much better could you predict it? Split by batch, never by random
cells, so shared technical noise cannot inflate the answer.

| line | n | oracle | universal | wrong-context | gain | wins |
|---|---|---|---|---|---|---|
| HepG2 | 162 | 0.4037 | 0.2301 | 0.1698 | +0.1736 | 140/162 |
| Jurkat | 253 | 0.2884 | 0.1461 | 0.1275 | +0.1423 | 217/253 |
| RPE1 | 234 | 0.4658 | 0.1983 | 0.1916 | +0.2676 | 212/234 |
| **all** | **649** | **0.3812** | **0.1859** | 0.1618 | **+0.1953** | **569/649** |

**Random-split leakage control:** splitting cells at random instead of by batch gives
0.98x / 0.95x / 1.00x the same gain. So the oracle gain is **not** a technical-noise artifact.

In actual score units, through cell-eval2:

| line | no-context composite | oracle composite | headroom |
|---|---|---|---|
| HepG2 | 0.2175 | 0.5875 | +0.3700 |
| Jurkat | 0.2273 | 0.4956 | +0.2683 |
| RPE1 | 0.1870 | 0.6626 | +0.4755 |

**This reversed the project's standing recommendation to abandon context.** Context-specific
effects are genuine, reproducible and worth a lot. The question became whether they are
*predictable from controls alone*.

### 12.3 The unseen-context benchmark — positive, in cosine

Strict protocol: every fold removes **all** perturbation data for the held-out biology; the model
sees only its pooled control profile; mu, both bases, the context PCA and both ridge maps are
fitted on training biologies only; predictions are scored against the oracle's batch-held-out
replicate.

| arm | correct | 95% CI on gain | wrong | shuffled | gain | oracle headroom recovered |
|---|---|---|---|---|---|---|
| no-context | 0.2377 | — | — | — | — | — |
| none (sanity) | 0.2396 | [+0.0006, +0.0031] | 0.2396 | 0.2396 | +0.0019 | 1.6% |
| raw pooled control | 0.2631 | [+0.0215, +0.0294] | 0.1999 | 0.2260 | +0.0254 | 21.2% |
| **PCA of pooled control** | **0.2645** | [+0.0239, +0.0299] | 0.2246 | 0.2304 | **+0.0268** | **22.3%** |
| oracle ceiling | 0.3578 | — | — | — | +0.1201 | 100% |

PCA wins 3/3 folds, beats wrong context by +0.0399 and shuffled by +0.0342.

Follow-up decomposition confirmed it was not a weak-baseline artifact. The optimal global rescale
fits to **exactly 1.0**, so it is not a magnitude correction. Against a rank-32 denoised baseline,
the same mechanism that makes refit beat mean transfer, context still adds **+0.0227**, positive
on all three lines. The gain splits nearly additively into a generic per-biology offset (+0.0161)
and a genuine target-by-context interaction (+0.0137).

**A prediction stated in advance was refuted here.** The expectation on record was that the
context arms would fail because the constraint is eight training context vectors. They did not
fail.

### 12.4 The same model on the real harness — it does not convert

`prediction = refit + alpha * residual`, scored as PDS on the five-line harness. The alpha=0 arm
reproduces the stored refit to 1e-5 per seed, so the comparison is faithful.

| arm | five-line PDS | vs refit |
|---|---|---|
| refit baseline | **0.8005** | — |
| interaction only, alpha 1 | 0.8011 | +0.0006 |
| **wrong context, alpha 1** | **0.8009** | **+0.0003** |
| both terms, alpha 0.5 | 0.8007 | +0.0001 |
| both terms, alpha 1 | 0.7976 | -0.0029 |
| both terms, alpha 3 | 0.7762 | -0.0243 |
| generic shift only | 0.7961 | -0.0045 |

**Wrong context scores 0.8009, tied with the best correct-context arm and better than the full
correct-context model.** The ablation that passed decisively in cosine fails here.

Per line, the failure is concentrated:

| line | refit | ctx alpha 1 | ctx alpha 3 | deviation/refit norm |
|---|---|---|---|---|
| K562 | 0.8821 | +0.0058 | +0.0087 | 0.140 |
| RPE1 | 0.8107 | +0.0415 | +0.0721 | 0.347 |
| Jurkat | 0.8596 | +0.0065 | +0.0134 | 0.119 |
| KOLF | 0.5658 | +0.0044 | -0.0056 | 0.290 |
| **HepG2** | 0.8846 | **-0.0727** | **-0.2101** | **0.487** |

Four lines gain modestly; HepG2 loses an order of magnitude more. **12 of 15 paired seeds "win"
while the mean is negative** — counting folds won would have called this a success.

**Nested alpha selection makes it worse, not better.** Choosing alpha for each line from the other
four gives 0.7650, worse than any fixed alpha, because HepG2's optimum is the opposite of
everyone else's. A hyperparameter that cannot be chosen out-of-sample is not deployable.

The direct comparison that matters: HepG2 scored **+0.0294 in cosine** (42% of its oracle
headroom) and **-0.0727 on PDS**. Same line, same representation, same folds, opposite sign. This
is the row-cosine trap firing in a case where the ablations looked clean — because they were clean
in the wrong metric.

### 12.5 Cell foundation models as context — no better than PCA

STATE SE-600M and UCE, both frozen, both run through their own published inference code and
weights: `state emb transform` from arc-state 0.11.1 with `se600m_epoch16.ckpt`, and
`eval_single_anndata.py` from snap-stanford/UCE with the 33-layer checkpoint and the paper's
figshare vocabulary. Nothing was gated; no vocabulary was reconstructed.

Interaction-only, since the generic shift caused the HepG2 collapse. Seven biologies, because
CD4 ships only pseudobulk and FENG has no non-targeting control arm. All on a common 6,413-gene
axis, because native panels span 8,248 to 38,584 genes in three assay groups and an FM would
otherwise separate biologies by assay panel rather than biology.

| arm, alpha 1 | five-line PDS | vs refit | correct − shuffled |
|---|---|---|---|
| **UCE, 1280-d** | **0.8039** | **+0.0034** | **+0.0005** |
| PCA of pooled control | 0.8030 | +0.0024 | +0.0059 |
| STATE SE-600M, 2058-d | 0.8012 | +0.0007 | +0.0063 |
| refit baseline | 0.8005 | — | — |

**UCE looks like the winner and is the trap.** Positive on 5/5 lines and 15/15 paired seeds, the
most consistent result in the table. Then shuffling its five context coordinates gives +0.0029
against +0.0034 correct. A shuffle preserves the vector's norm and the multiset of its values
while destroying which biology it names, so UCE's gain tracks **correction magnitude, not context
identity**.

PCA is the only genuinely context-specific arm with a non-trivial gain, and it is one line:
removing RPE1 (+0.0211) turns its five-line mean from +0.0024 to **-0.0023**.

Nothing clears the ±0.010 floor. A 600-million-parameter cell model extracts no more usable
unseen-biology context than a PCA of the same cells. **The binding constraint is six training
context vectors per fold, not encoder richness.**

A bug found and corrected mid-experiment: per-component whitening of n training points projected
onto n-1 PCA components is exactly degenerate — the whitened Gram becomes `n(I - 11^T/n)` for any
input, collapsing all three representations onto the same simplex. The symptom was all three
wrong-context arms scoring identically. Verified numerically at d = 6413, 2058 and 1280, replaced
with a single global scale, and the run redone from scratch.

---

## 13. The scoreboard

### Official submissions, all 14

| # | name | entry id | what it was | overall | rank |
|---|---|---|---|---|---|
| 1 | prior A | `EbsheSl3vTFozhB8iQD5` | Replogle transfer, deterministic rounding | -0.022839 | 468 |
| 2 | prior B | `06DcGcHX5ANL1ozEzIpt` | STATE ST from scratch | -0.191023 | 609 |
| 3 | B1 | `y9ZtKe2wxOH8R18tuvZl` | CD4 alone, deterministic | -0.006427 | 484 |
| 4 | B2 | `sYelQBhif6IU6lrMzJ9x` | 3-source mean, deterministic | +0.002011 | 467 |
| 5 | B6 | `Taz5toxJi7SCwkcTJIJk` | 6-source mean + stochastic rounding | +0.077057 | 342 |
| 6 | B7 | `oVFTLcfjIdXji3vhtNlt` | sign-agreement + stochastic | -0.022366 | 550 |
| 7 | **B8** | `KkMccCXITCezd16oBBr8` | 6-source mean + stochastic + **delta x3** | **+0.121709** | 214 |
| 8 | B10 | — | plain mean transfer at x4 | +0.0868 | 328 |
| 9 | B11 | — | expression-gated per-context, x2.25 | +0.1156 | 234 |
| 10 | B14 | — | precision-weighted sources, x3.5 | +0.0888 | 345 |
| 11 | B15 | `NTQy2s1JzkCMQpJcMSjy` | effective-count source weighting, x3 | +0.106714 | 287 |
| 12 | B16 | `vZ26xQ3Xm0PSGCKSfAJo` | frozen rank-128 covariance denoiser, x3 | +0.129193 | 199 |
| 13 | B18 | `Ax90uDSJhyMEFBacq5WB` | **byte-identical to b16** | +0.129193 | 231 |
| 14 | **B19** | `pamfw9b1glTBvlXvEyQA` | context-forcing residual, rank 32, alpha 1, x3 | **+0.136470** | **193** |

b17 was built and validated but never uploaded. b9, b12, b13 built, never uploaded.

**Current live standing, checked 2026-09-20: rank 225 of 1,071 entries. Our 0.13647 against a
top score of 0.34485.** The rank drifts as entries accumulate; it was 193 at submission and 205 in
the 2026-09-18 snapshot of 1,036 entries.

### The gap to the top

Against the 2026-09-18 snapshot leader at 0.33036, the deficit decomposes roughly as expression
38%, PDS 34%, log-FC 16%. Our expression metric scores **0** because our raw value is 15.5 against
a ~1.0 no-change point — we are far past the point where the metric gives any credit.

### Local five-line harness, everything measured

| predictor | five-line PDS |
|---|---|
| no change | 0.500000 |
| mean transfer | 0.736193 |
| b16 (archived) | 0.754619 |
| b19 reference | 0.795777 |
| seed-ensemble K=3 | 0.796985 |
| **refit baseline** | **0.800538** |
| UCE context on refit | 0.803904 |
| **masked wide-axis refit** | **0.801589** |

Note the last two. The **masked wide-axis refit at 0.801589** is a genuine +0.00105 over the
frozen refit and was never promoted or reported. The **UCE arm at 0.803904** is larger still but
fails its shuffle control, so it is not real.

### What the oracle says is available

| quantity | value |
|---|---|
| Same-context oracle over universal, cosine | +0.1953 |
| Fraction of oracle headroom recovered by best context model | 22.3% in cosine, ~0 in PDS |
| Composite headroom from no-context to oracle | +0.27 to +0.48 |
| Perturbation-specific share of deviation | 94% |
| Reproducible share of deviation | 27.7% mean, 0.7% on KOLF |
| Cross-biology transferable component | r ~ 0.02-0.07 |
| Reproducibility ceiling (Jiang replicates) | r ~ 0.21 |

---

## 14. Every retraction, bug and trap

This is the section most worth reading. The project's reliable output was negative results and
corrections.

### Claims retracted after measurement

1. **"We train in CP10K and are scored in CP50K, and that is why expression accuracy is zero."**
   False. The reconstruction is the exact inverse of log1p(CP10K); a perfect delta reproduces
   cells exactly. What survives is an objective-weighting mismatch: the metric weights
   low-abundance genes ~3x more than the training loss.
2. **"Dropping the x3 scale would recover expression accuracy."** False. At x1 the KOLF expression
   error is still 2.63 against a 1.00 no-change point.
3. **"The KOLF learned gain proves supervised transfer."** False. The training-label shuffle scores
   *higher* (0.2449 vs 0.2359). The mechanism is covariance denoising.
4. **"K562 learned context in the linear family."** Algebraically impossible: the context term
   cancels under query target-centring, confirmed by replay at 1-2e-7.
5. **"Better protein language models will help."** Random Gaussian vectors beat ESM-2 15B, and
   deleting gene features entirely beats both.
6. **"The shared-direction hypothesis explains the embedding null."** Centring dropped inter-gene
   cosine from +0.650 to +0.02 and changed nothing.
7. **"Context is useless; abandon the branch."** Overturned by the oracle test, which showed a
   +0.1953 same-context gain with no leakage.
8. **"The context arms will fail to beat no-context."** Stated in advance, refuted: PCA context
   recovered 22.3% of oracle headroom.
9. **"PDS deletes the generic response, so the generic shift cannot help."** The opposite: the
   generic shift carried nearly the entire K562 gain, and is the term that destroys HepG2.
10. **"Cosine is scale-invariant, so PDS does not move with delta scale."** False. Raw PDS rises
    47% from x1 to x4 because the delta passes through nonlinear reconstruction.
11. **"Replogle is available only as pseudobulk."** False; the `raw_bulk` files are per-cell means.
12. **"Tahoe has ~1,100 distinct drugs."** It has **379 distinct drugs** in 1,135 drug-dose
    conditions.
13. **The TRADE "0.16 to 0.90" figure** was never verified and must not be repeated.
14. **"PerturBench covers 30 independent Jiang cell lines."** Those are 30 biological *states*,
    not 30 lines.

### Bugs that changed results

- **Gene axis cut by evaluation panels.** 5,592 of 17,847 available genes, because the assembler
  intersected every source axis including unused Tahoe drug data. Cost measured at -0.062 baseline
  / -0.068 model. Fixed in `table_fullaxis`; fix never deployed.
- **RN effects paired with the wrong control.** Batch-weighted effects matched to pooled controls,
  median L2 error 0.646. Fixing it raised distinct contexts from **205 to 16,263**.
- **Contrast targets carried group-centring artifacts.** Equal raw effects in two biologies
  produced a latent contrast of 1.407 instead of zero. Caught before protocol freeze.
- **The cyclic-shift null is biased.** Alphabetical cyclic shuffles score 0.54-0.58; 32 random
  derangements score near 0.5, and every cyclic value exceeded all 32 random draws.
- **Per-component whitening of n points on n-1 components is exactly degenerate.** Whitened Gram
  becomes `n(I - 11^T/n)` for any input. Symptom: all wrong-context arms scoring identically.
- **The harness sorts targets; the CSV does not.** 1 of 279 positions matched. Invalidated three
  experiments in one session.
- **HGNC alias trap.** HGNC lists MED1 as an alias of MBD4. A naive alias map silently merges the
  Mediator subunit into MBD4, inflating three intersections from 0 to 3. The rule adopted: never
  remap a token that is itself an Approved symbol.
- **`cell_load` with `num_workers > 0`** leaks ~2 GB per 30 s per worker and OOM-kills training.
- **A local `scripts/coverage.py` shadows the `coverage` package**, breaking numba imports for
  anything run from that directory. Hit twice, most recently this session.
- **The uploader could create a duplicate entry** by misclassifying post-create errors as
  pre-creation refusals.
- **Two OOM kills** from float64 promotion inside the frozen baseline (47 GiB peak), fixed by a
  float32 reimplementation that is bit-identical and peaks at 11 GiB.

### Documents that contradict their own data

- `esm_arms_2026-09-19/FINDINGS.md` tabulates 13 of 16 arms and omits `zero_64`, the best one.
- `data_audit_2026-09-18/FINDINGS.md` says the masked baseline never produced a number; its own
  `masked_baseline/summary.json` contains **0.8015888835448705**.
- The `--scale` help text claims PDS is scale-invariant; the KOLF sweep shows it rising 47%.
- b15's rank is 277 in two records and 287 in two others.

### Things measured and explicitly not to repeat

Global scale tuning; pruning or thresholding the delta (fidelity punishes it structurally);
sign-agreement weighting; empirical-Bayes shrinkage as a ranking; adding other on-disk contexts to
the base table (0-5 of 300 panel coverage); per-source amplitude calibration (-0.003 PDS).

---

## 15. What a successor should do

Ranked by expected value, with the evidence for each.

**1. Score the refit baseline officially. It has never been submitted.** It is the best local
predictor at 0.800538 versus b19's 0.795777, and b19 is what is on the leaderboard. This is a
free experiment that resolves whether the local harness predicts the leaderboard at all — which
is currently unknown and which b15 suggests it may not.

**2. Deploy the wide-axis masked baseline.** It measures **0.801589**, +0.00105 over the frozen
refit, and the artifacts already exist. The audit built the fix and never shipped it. Caveat:
KOLF declines slightly.

**3. Re-examine the delta scale.** We submit at x3; the KOLF bundle peaks at x2 (0.1972 vs
0.1944). The reasoning recorded in the code for choosing x3 is demonstrably wrong on one of its
two claims. Tempered by KOLF's demotion as a gate.

**4. Add independent biologies. This is the real bottleneck.** Every line of evidence converges:
biologies matter ~10x more than genes; the context branch is a 58-value source lookup; six to
eight training context vectors cap every context model; a 600M-parameter cell FM extracts no more
than PCA from the same seven biologies. No architecture change will fix a sample-size limit.

**5. Attack expression accuracy, which is 38% of the gap.** We score exactly 0 on it, with a raw
value of 15.5 against a ~1.0 no-change point. Nothing in this project has targeted it directly.

**6. If context is revisited, use pooled per-biology controls, interaction terms only, and a
per-line trust gate.** Never per-target matched controls (1/cell-count noise). Never the generic
shift (it caused -0.21 on HepG2). The open problem is fold-level robustness: four lines gain and
one catastrophically loses, and nested selection makes it worse because HepG2's optimum is
inverted relative to every other line.

**7. Do not pursue frozen cell foundation models as context on this data.** Measured: neither
STATE SE-600M nor UCE beats a PCA of the same cells, and the better of the two fails its own
shuffle control.

### The methodological rules this project paid for

- Judge on PDS, on all five lines, and look at the **worst fold**, not the mean.
- Always run a **shuffle control**, not just a wrong-context control. Wrong context tests whether
  the model uses *a* context; shuffle tests whether it uses the *right* one. UCE passed the first
  and failed the second.
- Count folds won *and* the mean. Twelve of fifteen paired wins coexisted with a negative mean.
- Select hyperparameters out-of-sample or not at all. Nested selection revealed a method that
  looked positive at every fixed alpha was undeployable.
- A gain smaller than its own seed standard deviation is not a gain.
- Before blaming the pipeline, score the **oracle delta** through it. Three separate
  "the data must be broken" hypotheses died against that one test, which returns PDS 1.000.
