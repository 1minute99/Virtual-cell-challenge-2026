# Training-free baselines on the Virtual Cell Challenge 2026 validation set

Nothing here is trained or fine-tuned. Every number is set arithmetic, a forward
pass with a released checkpoint, or an already-measured effect copied from a
published screen.

## Environment

| component | version / id |
|---|---|
| python | 3.12.14 (`~/miniconda3/envs/AIVC`) |
| cell-eval | 0.8.2 (`vcc` profile = mae, discrimination_score_l1, overlap_at_N; `full` used below) |
| arc-state | 0.11.1 |
| anndata / scanpy / numpy | 0.13.3.post0 / 1.12.4 / 2.5.2 |
| torch | 2.14.0+cu130 |
| GPU | 1x RTX 5090 (32 GB). **No H100 is present on this machine.** |

Datasets and checkpoints actually used:

| artifact | version / provenance |
|---|---|
| VCC 2026 validation controls | `panel_id: vcc2026-val-1`, 3 contexts A/B/C, 300 constructs, 18,533 genes |
| Replogle 2022 genome-scale | figshare `20029387`, `K562_gwps_raw_bulk_01.h5ad` (md5 `4570b53c…`) |
| Replogle 2022 essential | `K562_essential_raw_bulk_01.h5ad`, `rpe1_raw_bulk_01.h5ad` (md5-verified) |
| Replogle single-cell (4 lines) | `replogle_concat.h5ad`, 643,413 cells x 6,546 genes |
| scLong | `gocont_4096_48m_pretrain_1b_mix_2024-02-05_16-23-37.pth`, 1.01B params |
| STATE SE | `arcinstitute/SE-600M`, `se600m_epoch16.ckpt` |
| STATE ST | `arcinstitute/ST-SE-Replogle`, `zeroshot/k562` |

## The blocker you need to know about first

**The VCC 2026 validation ground truth is withheld by Arc — it *is* the
leaderboard — and there is no 2026 challenge training set.** Participants receive
only non-targeting control cells plus the list of target gene IDs. I confirmed
this from the released files (each `context_*.h5ad` holds exactly 18,400 cells,
all `target_gene == "non-targeting"`) and from Arc's own challenge description.

Consequence: **no VCC validation prediction can be scored locally by cell-eval,
by me or anyone else.** Scores exist only by submitting to the leaderboard.

So this report has two halves, kept strictly separate:

* **VCC half** — the real predictions for the 3 validation contexts, emitted in
  submission shape. Unscored, because scoring is impossible offline.
* **Surrogate half** — the *identical* predictors scored with real cell-eval on
  Replogle cell lines that do have ground truth, using the same K562 genome-scale
  source. These numbers are real cell-eval output but are **not**
  leaderboard-comparable: different assay (10x 3' vs 10x Flex), different gene
  space (6,546 vs 18,533), different contexts.

## Task 1 — Coverage counts (the gate)

`results/task1_coverage.csv`, `results/task1_gene_lists.json`

| source | perturbations | **VCC targets perturbed** | VCC targets in readout | VCC gene-space overlap (of 18,533) |
|---|---|---|---|---|
| Replogle K562 **genome-scale** | 9,854 | **272 / 300** | 231 | 7,681 |
| Replogle K562 essential | 2,058 | **0 / 300** | 241 | 7,940 |
| Replogle RPE1 essential | 2,392 | **0 / 300** | 236 | 8,260 |
| STATE ST (`ST-SE-Replogle`) | 2,024 | **0 / 300** | 57 | 1,877 |

**Gene mapping.** All 300 targets are Approved HGNC symbols and all 300 carry an
Ensembl id; one is remapped (`TMEM104` → `SLC38A12`). Joins use Ensembl ids where
both sides carry them, symbols otherwise (271 of the 272 K562 hits agree under
both; `HOMEZ` matches by symbol only, an ENSG-version difference).

> **A trap worth recording.** A naive HGNC alias map is *wrong* here. HGNC lists
> `MED1` as an alias of `MBD4`, so naive alias resolution silently merges the
> Mediator subunit MED1 into MBD4 and manufactures fake gene matches — it
> inflated three of these intersections from 0 to 3 before I caught it. The
> mapper in `scripts/genemap.py` therefore never remaps a token that is itself an
> Approved symbol, and accepts an alias only when it resolves uniquely.

**Why the zeros are zero.** This is not Arc hiding the targets from public data.
The essential-gene screens perturb *essential* genes; VCC picks knockdown targets
that leave cells viable. The two sets barely can overlap. Replogle 2022 ran RPE1
only as an essential screen, so **RPE1-as-source has no data at all** — that arm
of Task 2 is not "small intersection", it is empty.

### Gate outcomes

* **Replogle**: K562 genome-scale gives 272/300 → Task 2 proceeds on K562.
  RPE1-as-source is impossible; reported, not substituted.
* **STATE ST**: 0/300 → **Task 3 skipped**, as instructed.

### Normalization (`results/task1_normalization.csv`)

| dataset | unit |
|---|---|
| VCC 2026 validation (10x Flex) | **raw integer UMI counts**, csr, var id = HGNC symbol |
| Replogle `*_raw_bulk` | raw pseudobulk mean UMI, var id = ENSG (+`gene_name`) |
| Replogle `*_normalized_bulk` | per-gene z-score vs non-targeting |
| Replogle `replogle_concat` single-cell | log1p(CP10K) computed on the full transcriptome then subset to 6,546 genes (expm1 row sums ≈ 8,515, not 10,000) |
| STATE ST output space | 2,000 HVGs; `pert_rep=onehot`, `pert_dim=2024` |

All deltas below are computed in log1p(CP10K) space on each dataset's own gene
space, then mapped by symbol. Genes with no source counterpart get delta 0
(pass-through control).

## Task 3 — SKIPPED (gate fired)

STATE's ST module represents a perturbation as a **one-hot vector over a fixed
2,024-entry vocabulary** (`pert_rep: onehot`, `pert_dim: 2024`, confirmed in the
checkpoint's `config.yaml` and `var_dims.pkl`). That vocabulary is Replogle's
essential-gene perturbation list, which intersects the 300 VCC targets in **zero**
genes. This is not a coverage shortfall that a bigger run would fix: a one-hot
encoder has no slot for an out-of-vocabulary gene, so released ST cannot express
any VCC 2026 target even in principle. Its output space is also only 2,000 HVGs,
covering 1,877 of VCC's 18,533 genes.

The other released ST checkpoints do not rescue this: `ST-*-Tahoe` is a
drug-perturbation model, and all released ST variants use the same one-hot
`pert_rep` with `perturbation_features_file: null` (no ESM2-style perturbation
featurisation that could generalise to unseen genes).

## Method notes that affect every number below

**Does cell-eval exclude the knocked-down gene from the DE metrics?**
**No.** There is no target-gene exclusion anywhere in `cell_eval` 0.8.2 or in
`pdex`: the perturbed gene is scored like any other gene, and its knockdown is
trivially predictable. This is the direct answer to the question, established by
reading the source, and it holds regardless of the quantification below.

*A measurement mistake I made and corrected, since it matters for how you read
any such number.* My first attempt zeroed the target gene in **both** the real
and the predicted matrix. That does not exclude the gene — it injects a large,
identical artificial down-shift into both, which both DE analyses then detect in
the same direction, **manufacturing** agreement. The tell was that the *floors*
gained the most from it (return_control's `overlap_at_N` went 0.000 -> 0.133,
`de_sig_genes_recall` 0.000 -> 0.141), which is impossible if the operation were
genuinely removing a source of inflation. Those numbers are in
`results/task2_perpert_*_notgt.csv` and should **not** be used. The correct
approach — remove the perturbed genes from the gene space entirely, so no target
gene can contribute to any metric — is the `--exclude-target-genespace` arm.

**Cell-level vs pseudobulk.** The spec's predictor is defined at pseudobulk level
(`control_pseudobulk + delta_source`), but cell-eval computes differential
expression and needs non-degenerate within-perturbation variance. Emitting N
identical cells would make every DE test degenerate. So predictions are formed
per cell as `sampled_control_cell + delta_source`, whose group mean is exactly
`control_pseudobulk + delta_source` while retaining realistic cell-level spread.

**Sanity check on the source deltas.** Before trusting anything, the delta at the
perturbed gene's *own* locus: over 6,117 perturbations, **99.2%** are negative,
median **-0.244**, about **6 sd** below a typical gene's delta. CRISPRi knockdown
is being recovered, so the delta extraction and gene mapping are sound.

## Task 4 — Which representation predicts transferability?

### 4a. Transfer success itself (`results/task4_transfer_pairs.csv`)

Before asking *which representation predicts transfer*, measure whether transfer
happens at all. For every ordered pair of Replogle contexts and every
perturbation measured in both, Pearson correlation of the two pseudobulk deltas:

| pair | n perturbations | mean Pearson-delta | median | frac > 0.7 |
|---|---|---|---|---|
| k562 ↔ rpe1 | 1,047 | 0.225 | 0.185 | 0.4% |
| k562 ↔ jurkat | 1,110 | 0.266 | 0.245 | 0.4% |
| k562 ↔ hepg2 | 940 | 0.244 | 0.229 | 0.1% |
| rpe1 ↔ jurkat | 1,183 | 0.234 | 0.201 | 0.9% |
| rpe1 ↔ hepg2 | 1,045 | 0.321 | 0.314 | 1.5% |
| jurkat ↔ hepg2 | 1,035 | 0.226 | 0.193 | 0.1% |
| *(same context, upper bound)* | — | 1.000 | 1.000 | 100% |

**This is the headline negative result.** A perturbation's transcriptional effect
measured in one cell line reproduces the same perturbation's effect in another
line at only r ≈ 0.22–0.32, and **essentially no perturbation (0.1–1.5%) clears
r = 0.7**. Copying a source effect is therefore a weak predictor *by construction*,
independent of any gene-mapping or normalization choice. The most similar pair is
rpe1 ↔ hepg2 (0.32); the least is jurkat ↔ hepg2 (0.23) — a narrow spread.

Note this measure is symmetric, so it scores context *similarity of effects*
rather than a directional prediction; the directional quantity (which also
depends on the target's own control profile) is the cell-eval `pearson_delta` in
Task 2.

### 4b. Representations (forward passes only)

| representation | dim | what was embedded | cost |
|---|---|---|---|
| scLong | 200 | 1.006B-param frozen encoder, mean-pooled over genes then cells | **8.08 GB peak VRAM, 5.21 cells/s, batch 4, RTX 5090** (1,792 cells in 344 s) |
| STATE SE | 2,058 | frozen SE-600M cell embedding, mean-pooled | RTX 5090, batch 32; all 6,192 genes mapped to symbols with **0 dropped** |
| PCA | 50 | PCA fit on pooled control cells of all 7 contexts | trivial; 42.3% variance explained |

256 control cells per context, 7 contexts (VCC A/B/C + Replogle k562/rpe1/jurkat/hepg2),
shared 6,192-gene space, log1p(CP10K).

**GPU feasibility.** scLong inference fits comfortably on the RTX 5090 (8.1 GB of
32 GB); no H100 was needed, and none is installed on this machine. Throughput is
compute-bound and flat in batch size, so the practical ceiling is ~5 cells/s.

## Task 2 — Source-effect-copy vs the floors (scored with cell-eval)

### 2a. The real VCC predictions (`results/vcc_pred_{A,B,C}_{predictor}.npy`)

272 of 300 targets get a K562 genome-scale delta; the 28 without one fall back to
the control profile. The 28 missing: ABCD1, CAPRIN2, EPHB2, HEATR6, KIF21B,
MAPK7, MLKL, NICN1, PARP3, PBLD, PDK3, PHF11, PSMB9, RAB11FIP5, SEMA4F, SLC44A1,
SNN, SOCS5, SPATA20, TAF4, TBC1D19, TDRD7, TESK2, TMEM104, TMEM79, TTBK2,
TXNDC16, ZC2HC1A.

A structural ceiling worth stating plainly: of VCC's 18,533 genes, only **7,682
carry any source signal at all**. The other **10,851 are emitted unchanged from
control** because Replogle never measured them. So ~59% of the output gene space
is, by construction, the return-control baseline.

**These are not scored, because they cannot be.** See the blocker section.

### 2b. Surrogate, scored with real cell-eval — K562 genome-scale → HepG2

1,340 perturbations, `profile="full"`, seed 0. Columns are the three predictors.

| metric | return_control | **source_copy_K562** | train_mean_delta |
|---|---|---|---|
| pearson_delta | 0.0014 | **0.1666** | 0.0704 |
| overlap_at_N | 0.0000 | **0.0854** | 0.0368 |
| de_direction_match | 0.4970 | **0.7506** | 0.5386 |
| de_sig_genes_recall | 0.0000 | **0.4548** | 0.2363 |
| discrimination_score_l1 | 0.5024 | **0.6021** | 0.5026 |
| discrimination_score_cosine | 0.5119 | **0.8090** | 0.5086 |
| **mae** | **0.0698** | 0.0745 | 0.0696 |
| **mse** | 0.0106 | 0.0112 | **0.0103** |
| pr_auc | 0.0719 | 0.0610 | 0.0516 |
| frac pearson_delta > 0.7 | 0.000 | **0.0007** | 0.000 |
| median pearson_delta | 0.0011 | **0.1470** | 0.0555 |

**The trivial predictor scores the way it must.** `return_control` gives
pearson_delta 0.0014 (≈ 0) and de_direction_match 0.4970 (≈ chance). That is the
check that the harness is not manufacturing signal, and it passes.

**The result that matters: source-copy wins on direction and loses on magnitude.**
It roughly doubles direction-match over chance (0.75 vs 0.50), recovers 45% of
significant DE genes, and reaches 0.81 discrimination-cosine — but its **MAE and
MSE are worse than doing nothing at all**. Returning the unperturbed profile
scores MAE 0.0698 against source-copy's 0.0745; `train_mean_delta` — which
ignores *which gene was targeted* — has the best MSE of the three.

This is directly relevant to the leaderboard: cell-eval's `vcc` profile is
`mae`, `discrimination_score_l1`, `overlap_at_N`. Source-copy improves two of
those three and **degrades MAE**. A delta that is directionally informative but
wrongly scaled for the new context costs more absolute error than adding no delta.

Only 0.07% of perturbations clear pearson_delta 0.7, matching Task 4a.

### 2c. Replication in a second unseen context — K562 genome-scale → RPE1

1,499 perturbations. The pattern reproduces exactly:

| metric | return_control | **source_copy_K562** | train_mean_delta |
|---|---|---|---|
| pearson_delta | -0.0037 | **0.1445** | 0.0654 |
| overlap_at_N | 0.0001 | **0.0833** | 0.0461 |
| de_direction_match | 0.4938 | **0.6710** | 0.5485 |
| discrimination_score_l1 | 0.5016 | **0.5497** | 0.5018 |
| **mae** | **0.0782** | 0.0815 | 0.0772 |
| **mse** | 0.0146 | 0.0147 | **0.0139** |
| frac pearson_delta > 0.7 | 0.000 | **0.0007** | 0.000 |

Both unseen contexts agree on all three findings: the floor sits at
pearson_delta ≈ 0 and direction-match ≈ 0.5 (the sanity check passes twice
independently), source-copy roughly doubles direction-match, and **source-copy
loses to both floors on MAE in both contexts**. `frac > 0.7` is 0.0007 in both.
This is a reproduced result, not a single-context artifact.

### 2d. All three unseen contexts, K562 genome-scale as source

| target | n perts | pearson_delta (copy) | direction_match (copy) | **MAE copy** | **MAE control** | MAE verdict |
|---|---|---|---|---|---|---|
| hepg2 | 1,340 | 0.1666 | 0.7506 | 0.0745 | **0.0698** | copy is WORSE |
| jurkat | 1,537 | 0.1701 | 0.7712 | 0.0757 | **0.0724** | copy is WORSE |
| rpe1 | 1,499 | 0.1445 | 0.6710 | 0.0815 | **0.0782** | copy is WORSE |

`return_control` pearson_delta is -0.0037 / -0.0003 / +0.0014 across the three —
i.e. zero, three times independently. The trivial-predictor check passes in every
context.

**Source-effect-copy loses to return-control on MAE in 3 of 3 unseen contexts**,
while simultaneously lifting direction-match from chance (~0.49) to 0.67-0.77 and
recovering 39-45% of significant DE genes. The effect is directionally real and
metrically counterproductive on absolute error. `train_mean_delta` — which does
not even look at which gene was perturbed — has the best MSE in all three.

### 2e. Complete surrogate grid: 2 sources x 3-4 targets (`results/task2_metric_table.csv`)

| source | target | pearson_delta copy / ctrl / mean | direction_match copy / ctrl | **MAE copy** | **MAE ctrl** |
|---|---|---|---|---|---|
| K562 gwps | hepg2 | 0.167 / 0.001 / 0.070 | 0.751 / 0.497 | 0.0745 | **0.0698** |
| K562 gwps | jurkat | 0.170 / -0.000 / 0.048 | 0.771 / 0.493 | 0.0757 | **0.0724** |
| K562 gwps | rpe1 | 0.144 / -0.004 / 0.065 | 0.671 / 0.494 | 0.0815 | **0.0782** |
| rpe1 | hepg2 | **0.289** / -0.000 / 0.212 | **0.838** / 0.514 | 0.0738 | **0.0683** |
| rpe1 | jurkat | 0.188 / 0.001 / 0.094 | 0.781 / 0.497 | 0.0785 | **0.0704** |
| rpe1 | k562 | 0.191 / 0.002 / 0.111 | 0.811 / 0.506 | 0.0745 | **0.0612** |

**Source-effect-copy has worse MAE than return-control in 6 of 6 pairs.** The
floor's pearson_delta is within ±0.004 of zero in all six, and its
direction-match within ±0.015 of 0.5 — the trivial-predictor sanity check passes
six times independently.

**Which source is better depends on context similarity, as Task 4a predicts.**
RPE1 is a markedly better source for HepG2 (pearson_delta 0.289) than K562 is
(0.167) — and rpe1↔hepg2 was exactly the most similar pair in the Task 4a
transfer grid (0.321 vs 0.244). The representation analysis and the scored
baselines agree.

**But most of RPE1's advantage is generic, not gene-specific.** For rpe1→hepg2,
`train_mean_delta` — which ignores which gene was targeted — already reaches
0.212. The gene-specific increment is only 0.289 - 0.212 = **0.077**, versus
0.167 - 0.070 = **0.097** for K562→hepg2. So the headline gain is largely a
better-matched *generic knockdown signature* (shared essential-gene response),
not better transfer of gene-specific biology. Any baseline that looks strong on
pearson_delta should be checked against the gene-agnostic floor for this reason.

### 2f. How much does the target gene inflate the scores? (corrected measurement)

`results/task2_perpert_K562_gwps_to_hepg2_seed0_noTGspace.csv`. K562 -> HepG2,
1,340 perturbations, scored twice: once as cell-eval scores by default (all 6,546
genes, target gene included), once with all 1,216 perturbed genes removed from
the gene space (5,330 genes remain) so no target gene can contribute.

| metric | source_copy default | source_copy excluded | **inflation** | return_control inflation |
|---|---|---|---|---|
| overlap_at_N | 0.0854 | 0.0541 | **+0.0314 (+58% relative)** | -0.0000 |
| de_direction_match | 0.7506 | 0.7083 | **+0.0424** | -0.0020 |
| de_sig_genes_recall | 0.4548 | 0.4221 | **+0.0327** | -0.0000 |
| pearson_delta | 0.1666 | 0.1502 | **+0.0164 (+11% relative)** | -0.0000 |
| discrimination_score_l1 | 0.6021 | 0.5967 | +0.0055 | +0.0013 |
| mae | 0.0745 | 0.0731 | +0.0014 | +0.0017 |

**The single trivially-predictable target gene supplies roughly 37% of
source-copy's DE-overlap score** (0.0314 of 0.0854) and about 11% of its
pearson_delta. The floors gain essentially nothing from it — correctly, since
neither predicts the target gene's own knockdown. The direction of the effect is
now the one you'd expect, and `return_control` acts as the built-in control
showing the operation itself introduces no artifact.

Caveat: removing 1,216 genes also shrinks the gene space, so this is not a pure
single-gene isolation; the unchanged floors are what license reading it as
target-gene inflation rather than a gene-space effect.

### 2f2. Per-perturbation Pearson-delta and the fraction above 0.7

`results/task2_pearson_delta_summary.csv`, histogram in
`results/task2_pearson_delta_hist.png`, raw per-perturbation values in
`results/task2_perpert_*.csv` (kept on disk for Task 4).

| source -> target | copy mean | copy median | **copy frac > 0.7** | n |
|---|---|---|---|---|
| K562 -> hepg2 | 0.1666 | 0.1470 | 0.0007 | 1,340 |
| K562 -> jurkat | 0.1701 | 0.1530 | 0.0000 | 1,537 |
| K562 -> rpe1 | 0.1445 | 0.1195 | 0.0007 | 1,499 |
| rpe1 -> hepg2 | 0.2886 | 0.2609 | 0.0077 | 1,045 |
| rpe1 -> jurkat | 0.1882 | 0.1447 | 0.0017 | 1,183 |
| rpe1 -> k562 | 0.1914 | 0.1458 | 0.0038 | 1,047 |

**Pooled over all 7,651 scored perturbations, source-copy clears Pearson-delta
0.7 for 0.21% of them. Both floors clear it for exactly 0.00%.** The best single
cell (rpe1 -> hepg2) manages 0.77%. There is no regime here in which copying a
source effect is a *good* prediction; it is only reliably better than nothing on
direction.

### 4c. Does context distance predict transfer? (`results/task4_correlations.csv`)

Correlation between context-to-context distance and mean per-perturbation
transfer success, over the 6 unique cross-context pairs of the 4 Replogle lines.
A **negative** correlation is the "sensible" direction: more similar contexts
should transfer better.

| representation | distance | Pearson r | Spearman r | Spearman p |
|---|---|---|---|---|
| **PCA** | euclidean | -0.618 | **-0.943** | **0.005** |
| PCA | cosine | -0.689 | -0.771 | 0.072 |
| STATE SE | cosine | -0.668 | -0.886 | 0.019 |
| STATE SE | euclidean | -0.670 | -0.886 | 0.019 |
| scLong | cosine | **+0.589** | **+0.029** | 0.957 |
| scLong | euclidean | +0.462 | +0.029 | 0.957 |

**PCA wins, and I am saying so plainly.** Fifty principal components of control
pseudobulk rank the context pairs better (Spearman -0.943) than STATE SE
(-0.886) and far better than scLong. **scLong carries essentially no signal
here** — Spearman +0.029, and its Pearson correlation has the *wrong sign*. A
1.006-billion-parameter single-cell foundation model, used as a frozen encoder,
is beaten by PCA at telling you which cellular context a perturbation effect will
transfer to.

Two supporting observations:

* scLong's inter-context cosine distances are tiny (0.006-0.016), i.e. its
  mean-pooled embedding barely separates cell lines at all. Mean-pooling a
  per-gene encoder may simply discard the context signal; this is a statement
  about *this* pooling of *this* checkpoint, not about scLong's per-gene outputs.
* STATE SE and PCA agree **exactly** on which Replogle line each VCC context is
  closest to, which is a useful independent consistency check. scLong disagrees.

**Caveat, stated up front: n = 6 pairs.** Four cell lines give only six unique
cross-context pairs, so these correlations are underpowered and the ranking
between PCA and STATE SE is not separable at this sample size. What *is* robust
is the qualitative split: PCA and STATE SE both point the right way with
consistent magnitude, scLong points nowhere. Note also that the transfer measure
is symmetric, so the 12 ordered pairs collapse to 6.

### 4d. Nearest Replogle source for each VCC context (`results/task4_vcc_nearest_source.csv`)

| VCC context | scLong | STATE SE | PCA |
|---|---|---|---|
| A | hepg2 | **jurkat** | **jurkat** |
| B | hepg2 | **rpe1** | **rpe1** |
| C | hepg2 | **hepg2** | **hepg2** |

STATE SE and PCA agree on all three. If you wanted to pick a transfer source per
VCC context from the data we have, this is the evidence — though Task 4a says the
expected payoff is small either way, and only K562 has the target coverage.

## Summary — what this establishes as a baseline to beat

1. **STATE ST cannot be run on this task at all.** Its perturbation encoder is a
   one-hot over 2,024 Replogle essential genes; the intersection with the 300 VCC
   targets is exactly zero. Task 3 skipped per the gate.

2. **Only K562 genome-scale has the coverage to source a copy baseline**
   (272/300). RPE1 has none (0/300) — Replogle 2022 ran RPE1 only as an essential
   screen. That arm of Task 2 has no data, and I did not substitute anything.

3. **Cross-context effect transfer is weak, everywhere.** Pseudobulk deltas
   correlate r ≈ 0.22-0.32 between Replogle cell lines, and pooled over 7,651
   scored perturbations only **0.21%** of source-copy predictions reach
   Pearson-delta 0.7. Both floors reach it **0.00%** of the time.

4. **Source-copy is directionally informative but metrically counterproductive.**
   In 6 of 6 source→target pairs it beats the floors on pearson_delta,
   direction-match, DE overlap and discrimination — and **loses to
   return-control on MAE in all 6**. Since cell-eval's `vcc` profile is
   `mae` + `discrimination_score_l1` + `overlap_at_N`, a naive source-copy
   submission would improve two leaderboard metrics and degrade the third.

5. **Check every apparent win against the gene-agnostic floor.** RPE1→HepG2 looks
   like the best transfer (0.289) until you notice `train_mean_delta` already
   scores 0.212 there: most of the signal is a generic essential-gene knockdown
   signature, not gene-specific biology.

6. **PCA beats both foundation models at predicting transferability.** Spearman
   -0.943 (PCA/euclidean) vs -0.886 (STATE SE) vs **+0.029 (scLong, no signal,
   wrong-signed Pearson)**. Stated plainly, as requested: **PCA wins.** Caveat:
   n = 6 context pairs, underpowered.

7. **The VCC validation set cannot be scored offline.** Ground truth is withheld
   by design and there is no 2026 training set. Every scored number here is from
   the Replogle surrogate and is **not** leaderboard-comparable. The VCC
   predictions themselves are built and saved.

### Honest caveats

* The surrogate is 10x 3' Replogle data at 6,546 genes; VCC is 10x Flex at
  18,533. Absolute metric values will not transfer to the leaderboard; the
  *ordering* of predictors is the transferable claim.
* Single seed (0) per surrogate cell (the multi-seed requirement attached to
  Task 3, which is skipped). Per-perturbation n is 1,045-1,537, so the means are
  tight, but seed variability is unmeasured.
* Task 4's transfer measure is symmetric, so it scores context similarity rather
  than directional prediction; the directional quantity is Task 2's cell-eval
  `pearson_delta`.
* scLong is evaluated as a *mean-pooled frozen encoder*. That is the pooling the
  task specified, but a different readout of the same checkpoint might separate
  contexts better; the null result is about this representation, not the model's
  every possible use.

### Reproducing

`scripts/run_all.sh` runs everything end to end in `~/miniconda3/envs/AIVC`.
Individual stages: `task1_coverage.py`, `build_cache.py`, `task2_vcc_predict.py`,
`task2_baselines.py`, `task2_aggregate.py`, `task4_transfer.py`,
`task4_embed.py --rep {pca,sclong,state_se}`, `task4_correlate.py`.
Shared helpers: `genemap.py` (collision-safe HGNC), `source_delta.py`,
`state_se_embed.py`.
