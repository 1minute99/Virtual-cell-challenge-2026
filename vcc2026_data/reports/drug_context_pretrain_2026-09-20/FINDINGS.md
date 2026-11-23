# Is context diversity the bottleneck? Drug-response pretraining for unseen-context gene perturbation

**Date:** 2026-09-20, closed 2026-09-21  **Status:** COMPLETE. 192 official-preset runs, 3 seeds x 5
lines, all PDS-parity checked. Scoring was interrupted three times by systemd-oomd (see
`CRASH_RECOVERY_2026-09-20.md`) and finished serially at 06:11 PDT on 2026-09-21.

## The question

Every gene-perturbation fold in this project trains on ~6 independent biological contexts. Every
context-conditioned model built on top of the refit baseline has failed in the same way: it either
ties refit or, when it moves, moves the wrong way, and a *wrong* context scores as well as the
correct one. The open question has been whether that is a property of the problem or a property of
having six contexts.

Drug perturbation data answers it directly, because Tahoe supplies 48 contexts of the same kind of
measurement. If context-dependent modulation is learnable at all, it should become learnable
somewhere between 6 and 48 contexts, and the curve should show it.

The prediction is never replaced, only corrected:

    prediction = refit_gene_effect + alpha * context_specific_residual

`alpha = 0` reproduces the stored refit baseline **bit-identically** (harness PDS
0.8819509476031214 on K562 seed 0, equal to the value stored by the tiny-xattn run).

---

## 1. Data

### Tahoe-100M — already on this machine, nothing downloaded

A processed pseudobulk response table already existed at
`reports/simple_residual_2026-09-16/data/tahoe/table.npz` (1.03 GiB), built from 157 locally
cached parquet shards (16 GiB) covering 4,431,233 cells. Full inventory:
`tahoe_inventory.json`.

| property | value |
|---|---|
| rows | 13,011 pseudobulk conditions |
| contexts | 48 cancer cell lines (CVCL ids present) |
| perturbations | 185 small molecules |
| genes | 18,150 |
| plates | plate3, plate6 (primary) + plate14 (full biological repeat of all 48 lines) |
| controls | both DMSO_TF wells pooled in raw-count space within exact plate and line; 150 control profiles |
| density | 8,508 / 8,880 line x drug cells filled on the primary plates = **95.8 %** |
| dose | 184 drugs at 5 uM only, 1 drug at 0.05 uM only — **dose is not a usable covariate** |
| raw counts | `treated_sum_counts` and `control_sum_counts` are stored as uint32 |

Because raw counts are stored for both arms, the response space was a free choice rather than
inherited. The full single-cell dataset was **not** downloaded; it was not needed.

`HEPG2_FAMILY` (HepG2/C3A, 280 rows) is a Tahoe context and is the same biology as one of the five
scored gene folds. It is removed from drug pretraining on that fold.

### GSE306429 — downloaded from GEO

`gse_manifest.json`. Downloaded the processed pseudobulk plus two QC files; the 45 GiB RAW tar and
the 24 GiB single-cell h5ad were deliberately **not** fetched.

| file | bytes |
|---|---|
| `GSE306429_combined_pseudobulk.h5ad` | 119,871,669 |
| `GSE306429_combined_vscores.h5ad` | 77,748,812 |
| `GSE306429_combined_limma.csv.gz` | 1,757,519,980 |

Schema, verified rather than assumed: 1,737 rows (one per `bio_sample_id`, i.e. one 24-well well)
x 13,290 genes. `layers/counts` holds raw counts (row sums 1.6e7–2.4e7); `X` is already
normalised (row sums ~35,000, max 10.1) and was **not** used. 10 contexts — six immortalised lines
(A375, A549, HA1E, PC3, HEK293T, HepG2) and four primary/donor populations (CD34+, CD8+, HBEC,
PAD). 103 compounds over 11 doses (0.0001–10 uM) at a single 24 h timepoint. 58 plates, each a
single cell line, each carrying DMSO wells. Controls are `compound_name == "DMSO"` (dose is NaN on
exactly those 54 rows).

GSE contains **HepG2 and HEK293T**, both gene-perturbation biologies; both are dropped from drug
pretraining on their own fold. Only 7 compounds are shared with Tahoe, so the two drug datasets are
nearly disjoint in perturbation space.

### Response space (stated explicitly, as required)

    response = log1p(1e4 * treated_counts / treated_library)
             - log1p(1e4 * control_counts / control_library)

recomputed **from raw counts** for both datasets in `dcp_prepare.py`. Nothing is derived from an
already-normalised matrix, so there is no second compression on top of a first. This is the space
the frozen gene table uses and the space the VCC count reconstruction inverts, which is what makes
`refit + alpha * residual` a legal object.

Two independent checks that the preprocessing computes the intended quantity:

* Tahoe: the recomputed response matches the stored `effects` array to **max abs diff 4.8e-7**.
* GSE: for one well (Decitabine 1 uM, A375, replicate 1, plate ILD1-001) our delta against the
  official limma `logFC` over 13,290 genes gives **Spearman 0.756**, Pearson 0.618 on the 11,554
  well-expressed genes. The gap is limma's empirical-Bayes moderation, which we deliberately do not
  reproduce.

Controls are matched within dataset **and** within plate, never across.

### Gene axis

| axis | genes |
|---|---|
| gene-perturbation common axis | 6,413 |
| Tahoe | 18,150 |
| GSE306429 | 13,290 |
| **shared context axis (all three)** | **5,955** |

Responses stay on each dataset's own axis with its own low-rank basis (dataset-specific decoders);
only the *context* representation is forced onto the shared 5,955-gene axis so that the three
datasets' contexts live in one space.

### The datasets separate almost perfectly by assay — and centring fixes it

`preprocessing_diagnostics.json`. Untreated control profiles, before any correction:

| | raw | after per-dataset centring |
|---|---|---|
| nearest-neighbour same-dataset fraction | **0.985** | 0.631 |
| PC1 between/within variance ratio | **133.7** | 6.5e-16 |

Concatenating Tahoe and GSE without this would have modelled the assay, not the biology. All
downstream work uses the control-only per-dataset centred profiles.

### One bug worth recording

Context rows were first keyed by bare name. `A549` exists in both Tahoe and GSE, and `HEK293T` in
both GSE and the gene side, so a name→row dict silently handed those contexts another dataset's
control profile. Measured effect: Tahoe results were unchanged, but GSE's apparent context signal
was **entirely an artefact** — GSE A549 had been borrowing Tahoe's A549 profile and was by far its
best context (correct−wrong +0.150). After switching to dataset-qualified keys
(`tahoe:A549`, `gse:A549`, `gene:HEK293T`) GSE shows no usable context signal at all. Everything
below is post-fix; the pre-fix run is kept in `prekeyfix/`.

---

## 2. Architecture selection — on drug contexts only

Three families, all tiny, swept over rank {16, 32} and weight decay {1e-4, 1e-3}, plus experts
{4, 8} and operator rank {4, 8}: 20 configurations x 6 folds x 2 datasets = 240 runs
(`drug_lodo/results.json`, `sweep_table.json`). Folds hold out complete cell lines. The five-line
gene harness was not consulted.

**Every one of the 20 configurations separates correct from wrong context on 6/6 Tahoe folds.**
That has never happened in this project before.

| family | mean correct−wrong (Tahoe) | mean correct cosine |
|---|---|---|
| **FiLM** | **0.0955** | 0.0954 |
| operator | 0.0893 | 0.0968 |
| mixture | 0.0866 | 0.0912 |

Selected: **`film_r32_wd0.001`** — 7,019 trainable parameters, correct−wrong 0.0938, 6/6 folds.
Weight decay was immaterial (1e-4 and 1e-3 agree to four decimals at most settings).

---

## 3. The context-diversity scaling experiment (the central result)

Best architecture, Tahoe only, identical folds, identical budget; only the number of training
contexts changes. `scaling_table.json`.

| Tahoe contexts requested | actually trained on | correct resid cos | wrong | correct−wrong | folds c>w | **total-response cos gain** | resid MSE gain | folds with a gain |
|---|---|---|---|---|---|---|---|---|
| 6 | 6 | 0.0265 | 0.0016 | +0.0250 | 6/6 | **−0.0247** | −0.00046 | 0/6 |
| 12 | 12 | 0.0671 | 0.0004 | +0.0667 | 6/6 | −0.0277 | −0.00039 | 0/6 |
| 24 | 24 | 0.0940 | 0.0040 | +0.0900 | 6/6 | −0.0080 | −0.00010 | 0/6 |
| 47 | 40 | 0.1054 | 0.0074 | +0.0980 | 6/6 | **+0.0019** | +0.00002 | 3/6 |

(47 requested becomes 40 trained because each fold holds out 8 of the 48 contexts.)

Two things happen, and they are not the same thing:

1. **Context discrimination is present from 6 contexts and grows ~4x to 40.** Correct beats wrong
   on 6/6 folds at every size. The model can always tell contexts apart.
2. **Usefulness lags far behind.** At 6, 12 and 24 contexts the correction makes the *total
   predicted response worse* on every fold. It only crosses zero somewhere near 40 contexts, and
   even there it is break-even: over 48 held-out Tahoe contexts the mean total-response cosine is
   0.5214 with context and 0.5213 without.

Per-context detail at 40 training contexts (`drug_per_context_selected.json`):

| | Tahoe (48 held-out) | GSE (10 held-out) |
|---|---|---|
| correct > wrong | **41/48** | 6/10 |
| correct > zero-context | 36/48 | 1/10 |
| correct > no-context (net useful) | **28/48** | **0/10** |
| correct beats no-context on residual MSE | 24/48 | 1/10 |

This is the cleanest available statement of the problem: **discrimination saturates long before
usefulness does**, and six contexts is deep inside the harmful regime. Every prior null in this
project is consistent with sitting at the left edge of this curve.

---

*(Sections 4–7 follow: gene transfer, official VCC2026 preset results, context-specificity check
and the decision list.)*

## 3b. The gene task sits exactly where the drug curve says it should

The same measurement, run on gene perturbations only (leave-one-scored-biology-out over the five
harness lines, 2 / 4 / 6 training biologies, identical FiLM r32 model). `gene_curve.json`.

| training contexts | source | correct resid cos | wrong | correct-wrong | folds c>w | total-response cos gain |
|---|---|---|---|---|---|---|
| 2 | gene | 0.0270 | 0.0273 | -0.0003 | 3/5 | -0.0066 |
| 4 | gene | 0.0411 | 0.0238 | +0.0173 | 3/5 | -0.0028 |
| 6 | gene | 0.0630 | 0.0449 | +0.0181 | 3/5 | -0.0150 |
| 6 | drug (Tahoe) | 0.0265 | 0.0016 | +0.0250 | 6/6 | -0.0247 |
| 12 | drug (Tahoe) | 0.0671 | 0.0004 | +0.0667 | 6/6 | -0.0277 |
| 24 | drug (Tahoe) | 0.0940 | 0.0040 | +0.0900 | 6/6 | -0.0080 |
| 40 | drug (Tahoe) | 0.1054 | 0.0074 | +0.0980 | 6/6 | +0.0019 |

At six contexts the gene task gives correct-minus-wrong +0.0181 and a total-response cosine change
of -0.0150; the drug task at six contexts gives +0.0250 and -0.0247. Two independent datasets,
different assays, different perturbation modalities, same place on the curve. The gene task is not
a special failure - it is the left-hand end of a general one.

A second detail worth keeping: the gene-side *wrong-context* cosine is 0.0449 at six contexts
against 0.0074 for the drug model at forty. With few contexts most of what the model emits is
generic rather than context-specific, which is exactly the signature every earlier experiment in
this project reported as "wrong context scores as well as correct".


---

## 4. Gene transfer — development metrics (harness PDS, before the official preset)

Same FiLM r32 trunk (7,019 trainable parameters), five leave-one-scored-biology-out folds. The
held-out line's responses enter nothing: not training, not early stopping, not the PCA or
response basis, not alpha selection. Early stopping uses one inner training biology (Jurkat for
K562/RPE1/KOLF, K562 for HepG2/Jurkat); alpha for line L is chosen on the other four lines only
(`dcp_select.py`). Tahoe `HEPG2_FAMILY` and GSE `HepG2`/`HEK293T` are dropped from drug
pretraining on the folds that score those biologies (`gene_protocol.json`).

Arms: `scratch` (gene only), `tahoe_frozen` / `tgse_frozen` (pretrained trunk frozen, gene
adapters only), `tahoe_partial` / `tgse_partial` / `gse_partial` (trunk at 0.1x LR), `joint`
(gene 50 % / Tahoe 25 % / GSE 25 % batches), `joint_ft` (joint then 600 low-LR gene steps),
`ENS_plain` / `ENS_gated` (three-seed `joint` ensemble, ungated and uncertainty-gated).

### 4a. Five-line harness PDS versus refit 0.8005 (`gene_pds_table.json`, `gate_scores.json`)

| arm | alpha | mean PDS | delta vs refit | lines won |
|---|---|---|---|---|
| refit (alpha 0) | 0 | 0.8005 | 0 | – |
| scratch | 0.25 | 0.8005 | −0.0001 | 3/5 |
| tahoe_frozen | 0.25 / 0.5 / 1 | 0.8045 / 0.8064 / 0.8061 | +0.0040 / **+0.0058** / +0.0056 | 3/5 / 3/5 / 2/5 |
| tahoe_partial | 0.25 | 0.8012 | +0.0006 | 2/5 |
| tgse_frozen | 0.25 | 0.7977 | −0.0028 | 2/5 |
| tgse_partial | 0.25 / 0.5 | 0.8036 / 0.8038 | **+0.0030 / +0.0033** | **5/5** / 4/5 |
| gse_partial | 0.25 | 0.7953 | −0.0053 | 2/5 |
| joint | 0.25 | 0.8002 | −0.0003 | 3/5 |
| joint_ft | 0.25 | 0.8000 | −0.0006 | 2/5 |
| ENS_plain | 0.25 | – | −0.0017 | 1/5 |
| ENS_gated | 0.25 | – | −0.0001 | 0/5 |

Every gain is small and every one is carried by RPE1 (tahoe_frozen a1: RPE1 +0.049, K562
−0.017). tgse_partial at alpha 0.25 is the only arm that beats refit on all five lines, by
+0.0012 to +0.0058 per line. Alpha 1 hurts every family except tahoe_frozen. These are the
development numbers that chose the arm map; the ranking below uses the official composite.

### 4b. Did pretraining reduce the early overfitting? No — it made it immediate (`training.json`)

| arm | inner-val best step (5 folds) | early-stopped | folds where residual beats zero |
|---|---|---|---|
| drug pretraining, Tahoe | mean 4,660 of 6,000 (2,000–5,600) | 0/5 | 5/5 |
| drug pretraining, Tahoe+GSE | mean 4,860 of 6,000 | 0/5 | 5/5 |
| scratch (gene) | mean 760 (100–2,800) | 5/5 | 1/5 |
| tahoe_frozen / tgse_frozen / tgse_partial / gse_partial | **100 in 5/5 folds** (the first evaluation) | 5/5 | **0/5** |
| joint | mean 200 (100–300) | 5/5 | 2/5 |

On drug contexts the trunk trains for thousands of steps and beats the zero residual on every
fold. On gene contexts the same pretrained trunk's best inner-validation loss is at the very
first checkpoint and it never beats predicting zero residual (validation 0.156–0.235 against
zero's 0.138). Pretraining did not delay gene-side overfitting; the only improvement in
`beat_zero` comes from `joint`, which keeps drug data in the batch (2/5 folds).

Inner selection picked `joint` for K562/RPE1/Jurkat, `scratch` for HepG2, `tahoe_frozen` for
KOLF (`gene_diagnostics.json`), so the Part 19 controls are built on those.

### 4c. Uncertainty gate (`gate_diagnostics.json`)

The gate multiplies the residual by exp(−distance ratio) using the held-out context's distance to
the nearest training context in the shared PCA space, times ensemble agreement. Every gene
context sits **2.3x farther** from the drug training manifold than training contexts are from
each other (nearest-neighbour distance 16.8–17.2 versus training median 7.3–7.6; KOLF 21.4), so
the gate shrinks the correction to 7 % (KOLF 2 %). `ENS_gated` therefore equals refit to
1e-4 on every line. The gate does what was asked — it prevents harmful corrections on
out-of-distribution contexts — but only by turning the context model off, which is itself the
finding: the gene contexts are out of distribution for a model trained on drug contexts.

---

## 5. Official VCC2026 preset — the ranking that counts

`cell-eval2 run --preset vcc2026` (pdex backend, 400 cells per target, scale 3, the same frozen
real bundles and normalisation anchors as `preset_rescore_2026-09-20`) followed by
`cell-eval2 score` against each line's real bundle. Six metrics in the installed preset:
`pds_cosine`, `expr_mse_unbiased_capped_norm`, `de_wilcoxon_lfc_nmae`,
`de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`,
`de_wilcoxon_sig_jaccard`; the composite is their plain mean after per-line replicate
normalisation. Every run's raw `pds_cosine` equals the harness PDS stored at training time to
1e-9. Alphas are the nested choices from `dcp_select.py`. Full tables: `preset_report.log`,
`final_table.log`, `preset_table.json`, `final_table.json`; per-run receipts under `preset/runs/`.

### 5a. Five-line official composite (refit = 0.2755)

| model | pretraining | transfer | composite | delta vs refit | lines won | paired seeds won |
|---|---|---|---|---|---|---|
| refit (alpha 0) | none | – | 0.2755 | 0 | – | – |
| **scratch** | none | gene only | **0.2935** | **+0.0180** | 3/5 | 10/15 |
| **tgse_frozen** | Tahoe+GSE | frozen trunk | 0.2848 | **+0.0093** | **5/5** | **14/15** |
| tgse_partial | Tahoe+GSE | trunk at 0.1x LR | 0.2788 | +0.0033 | 2/5 | 7/15 |
| joint_ft | Tahoe+GSE | joint, short gene FT | 0.2783 | +0.0028 | 4/5 | 12/15 |
| ENS_plain | Tahoe+GSE | joint x3 seeds | 0.2765 | +0.0010 | 3/5 | 9/15 |
| joint | Tahoe+GSE | joint 50/25/25 | 0.2761 | +0.0006 | 3/5 | 6/15 |
| ENS_gated | Tahoe+GSE | joint x3, gated | 0.2755 | +0.0000 | 3/5 | 8/15 |
| gse_partial | GSE only | trunk at 0.1x LR | 0.2697 | −0.0058 | 4/5 | 12/15 |
| tahoe_frozen | Tahoe | frozen trunk | 0.2663 | −0.0092 | 3/5 | 9/15 |

### 5b. Per line (delta vs refit, seeds won)

| model | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|
| scratch | −0.0004 (1/3) | **+0.0729 (3/3)** | +0.0218 (3/3) | −0.0049 (0/3) | +0.0006 (3/3) |
| tgse_frozen | +0.0080 (3/3) | **+0.0365 (3/3)** | +0.0014 (3/3) | +0.0005 (3/3) | +0.0002 (2/3) |
| tgse_partial | −0.0002 | +0.0148 | −0.0073 | +0.0104 | −0.0011 |
| joint_ft | +0.0067 | +0.0048 | +0.0035 | +0.0035 | −0.0044 |
| joint | +0.0005 | −0.0007 | +0.0046 | +0.0003 | −0.0015 |
| tahoe_frozen | −0.0302 (a1) | +0.0392 | −0.0774 (a1) | +0.0090 | +0.0134 (a1) |
| refit composite | 0.3351 | **0.0553** | 0.4184 | 0.3638 | 0.2048 |

**One biology carries everything.** RPE1 contributes +0.0146 of scratch's +0.0180 and +0.0073 of
tgse_frozen's +0.0093. Without RPE1 the five-line gains are **+0.0043** (scratch) and **+0.0025**
(tgse_frozen). RPE1 is also where refit itself is broken: its normalised NMAE is −1.14 and reach
−0.24, so its composite is 0.055 against 0.20–0.42 elsewhere. Any residual that perturbs the
magnitudes helps there; scratch on RPE1 gains +0.30 NMAE, +0.06 reach, +0.05 PDS.

**Which metrics move.** Five-line normalised deltas (the quantities the composite averages):

| model | PDS | exprMSE | NMAE | fidelity | reach | jaccard |
|---|---|---|---|---|---|---|
| scratch | −0.0031 | 0 | **+0.0649** | +0.0079 | **+0.0341** | +0.0041 |
| tgse_frozen | **−0.0405** | 0 | **+0.0605** | +0.0069 | +0.0263 | +0.0027 |
| tgse_partial | +0.0061 | 0 | +0.0085 | −0.0006 | +0.0056 | +0.0003 |
| joint_ft | −0.0020 | 0 | +0.0062 | +0.0020 | +0.0101 | +0.0008 |

Every composite gain comes from NMAE and reach, with a little fidelity and Jaccard. PDS — the
only metric that discriminates targets (target-shuffled deltas keep 79 % of the composite, §3 of
the handoff) — is flat for scratch and **worse by 0.04** for tgse_frozen. The residuals are
improving magnitude calibration and direction-set metrics, not target identity.

**PDS mis-ranked again.** By harness PDS the best arms were tahoe_frozen (+0.0058) and
tgse_partial (+0.0030, 5/5 lines) and scratch was a tie. On the official composite scratch is
first, tgse_partial is +0.0033 and tahoe_frozen is **−0.0092** — its nested alpha of 1.0 on K562
and HepG2, chosen on PDS, costs −0.030 and −0.077 there. This is the third time in this project
that PDS-only development ranking reversed on the full preset.

## 6. Context-specificity check (Part 19), official composite at alpha 1

For each line's inner-selected family (`joint` for K562/RPE1/Jurkat, `scratch` for HepG2,
`tahoe_frozen` for KOLF): correct context, the same model fed a real different biology's control
profile, and zero context, all at alpha 1 so the three are matched.

| line | family | correct | wrong | zero | refit | scratch a0.25 | correct − wrong | correct − zero | seeds c > w |
|---|---|---|---|---|---|---|---|---|---|
| K562 | joint | 0.3311 | 0.2848 | 0.3286 | 0.3351 | 0.3347 | +0.0464 | +0.0025 | 3/3 |
| RPE1 | joint | 0.0409 | 0.0385 | 0.0044 | 0.0553 | 0.1282 | +0.0023 | +0.0365 | 3/3 |
| HepG2 | scratch | 0.3956 | 0.3762 | 0.3192 | 0.4184 | 0.4402 | +0.0194 | +0.0764 | 3/3 |
| Jurkat | joint | 0.3536 | 0.3022 | 0.3503 | 0.3638 | 0.3588 | +0.0514 | +0.0033 | 3/3 |
| KOLF | tahoe_frozen | 0.2182 | 0.2085 | 0.2123 | 0.2048 | 0.2054 | +0.0097 | +0.0059 | 3/3 |
| mean | | | | | | | **+0.0259** | **+0.0249** | 15/15 |

* `official(correct) > official(wrong)`: **yes, 5/5 lines, 15/15 seeds.** For the first time in
  this project a context model tells a real wrong context from the right one on the official
  aggregate. Wrong context is worse than zero context on 4/5 lines, so the model is using the
  context input and not merely emitting a generic shift.
* `official(correct) > official(scratch)`: **no, 1/5.** Only KOLF's correct-context arm beats
  both refit and scratch; on the other four, alpha 1 with the correct context is *below refit*
  (−0.004 to −0.023). The context information is real and small, and at full strength it is
  outweighed by the harm of the residual itself.
* On KOLF, the one line where correct beats refit (+0.0134), zero context also beats refit
  (+0.0075) and so does the wrong context (+0.0037): 56 % of the "context" gain is present with
  no context at all.

## 7. Decision list

1. **Does unseen drug-context performance improve from 6 → 12 → 24 → 47 Tahoe contexts?**
   Discrimination does, monotonically and ~4x (correct−wrong 0.025 → 0.067 → 0.090 → 0.098, 6/6
   folds at every size). Usefulness barely does: total-response cosine gain is −0.025, −0.028,
   −0.008, **+0.002**. The correction is harmful below ~24 contexts and break-even at 40.
2. **Can the architecture distinguish correct from wrong context on unseen drug contexts?** Yes:
   41/48 Tahoe contexts, every configuration, every fold. On GSE 6/10 and net-useful on 0/10.
3. **Does Tahoe pretraining improve unseen gene-context transfer over the identical scratch
   model?** No. Official composite: scratch +0.0180, tahoe_frozen −0.0092, tahoe_partial not
   scored (lost on nested PDS). Tahoe+GSE frozen is +0.0093 with the best consistency (5/5 lines,
   14/15 seeds) but still below scratch and only +0.0025 without RPE1.
4. **Does GSE306429 add value beyond Tahoe?** As a drug-side context source, no measurable
   context signal after the key fix (net useful on 0/10 GSE contexts). As pretraining data,
   tgse_frozen (+0.0093) beats tahoe_frozen (−0.0092), but the difference is dominated by the
   alpha choice, not by GSE content; gse_partial alone is −0.0058.
5. **Does joint training outperform sequential pretrain → fine-tune?** No. joint +0.0006 and
   joint_ft +0.0028 versus tgse_frozen +0.0093 and tgse_partial +0.0033.
6. **Is catastrophic forgetting visible?** There is nothing to forget: every sequential gene head
   is at its best after the first 100 steps, and the drug trunk is frozen or at 0.1x LR. Joint
   training delays the gene-side turnover from step 100 to 200–300, that is all.
7. **Does drug pretraining reduce the early overfitting?** No, it makes it immediate (best inner
   validation at step 100 in 20/20 pretrained-head runs versus a mean of 760 for scratch).
8. **Does the improvement survive the six-metric aggregate?** The PDS ordering does not (§5b).
   Two arms show positive composite deltas, both dominated by RPE1.
9. **Does correct context beat a real wrong context under the official aggregate?** Yes, 5/5
   lines, 15/15 seeds, +0.026. It does not beat scratch (1/5) or refit (1/5) at alpha 1.
10. **Which metrics improve, which degrade?** NMAE, reach, fidelity and Jaccard improve; PDS is
    flat or degrades (tgse_frozen −0.04); expression MSE is 0 everywhere.
11. **Is any gain driven by a single biology?** Yes, RPE1: it supplies 81 % of scratch's gain and
    78 % of tgse_frozen's, on a line where refit's own calibration is the worst in the panel.
12. **Does uncertainty gating prevent harmful OOD corrections?** Yes, trivially: gene contexts are
    2.3x farther from the drug manifold than any training context, the gate shrinks every
    correction to ≤ 7 %, and ENS_gated equals refit to 1e-4 on all five lines.
13. **Was insufficient context diversity the major bottleneck?** For *discrimination*, yes: the
    curve shows it directly and the gene task sits at its harmful left end. For *usefulness*, more
    contexts of the same kind only reach break-even at 40, on the drug data itself, and drug
    contexts do not transfer to gene contexts because the latter are out of distribution for the
    context encoder. The bottleneck is not the count of drug contexts we can add; it is the count
    of *gene-perturbation* biologies, which is fixed at six.
14. **Next investment?** Not more drug contexts (LINCS), not full single-cell Tahoe, not
    embeddings and not a richer interaction architecture: every one of those raises the right-hand
    end of a curve whose gene-side operating point cannot move. The measured lever on the official
    aggregate is magnitude calibration (NMAE/reach on RPE1-like lines), which is target-blind and
    belongs to the refit baseline, not to a context model.

## Final summary table

| model / strategy | drug data | architecture | drug ctx | gene transfer | PDS | exprMSE | NMAE | fid | reach | jac | **official aggregate** | Δ vs refit | correct − wrong |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| refit baseline | none | – | 0 | – | 0.8005 | 3.780 | 0.8430 | 0.6080 | 0.4401 | 0.0422 | **0.2755** | 0 | – |
| scratch interaction | none | FiLM r32 | 0 | gene only | 0.8005 | 3.982 | 0.8422 | 0.6119 | 0.4598 | 0.0434 | **0.2935** | +0.0180 | +0.0194 (HepG2) |
| Tahoe → frozen | Tahoe | FiLM r32 | 47 | frozen trunk | 0.7989 | 3.890 | 0.8489 | 0.6034 | 0.4267 | 0.0421 | 0.2663 | −0.0092 | +0.0097 (KOLF) |
| Tahoe+GSE → frozen | Tahoe+GSE | FiLM r32 | 55 | frozen trunk | 0.7882 | 3.781 | 0.8365 | 0.6116 | 0.4556 | 0.0429 | 0.2848 | +0.0093 | – |
| Tahoe+GSE → partial FT | Tahoe+GSE | FiLM r32 | 55 | trunk 0.1x LR | 0.8027 | 3.838 | 0.8434 | 0.6078 | 0.4434 | 0.0423 | 0.2788 | +0.0033 | – |
| GSE → partial FT | GSE | FiLM r32 | 8 | trunk 0.1x LR | 0.7953 | 3.787 | 0.8419 | 0.6083 | 0.4458 | 0.0421 | 0.2697 | −0.0058 | – |
| joint | Tahoe+GSE | FiLM r32 | 55 | 50/25/25 batches | 0.8002 | 3.814 | 0.8438 | 0.6087 | 0.4418 | 0.0423 | 0.2761 | +0.0006 | +0.0464/+0.0023/+0.0514 (K562/RPE1/Jurkat) |
| joint → short gene FT | Tahoe+GSE | FiLM r32 | 55 | 600 low-LR steps | 0.8000 | 3.795 | 0.8430 | 0.6091 | 0.4458 | 0.0424 | 0.2783 | +0.0028 | – |
| joint ensemble, ungated | Tahoe+GSE | FiLM r32 x3 | 55 | joint | 0.7988 | 3.816 | 0.8425 | 0.6086 | 0.4453 | 0.0424 | 0.2765 | +0.0010 | – |
| joint ensemble, gated | Tahoe+GSE | FiLM r32 x3 | 55 | joint + OOD gate | 0.8005 | 3.781 | 0.8430 | 0.6080 | 0.4403 | 0.0422 | 0.2755 | +0.0000 | – |

Raw six are five-line means of the un-normalised scorer outputs; the aggregate is the official
replicate-normalised composite. Correct − wrong is the alpha-1 comparison from §6 (mean over the
five lines +0.0259).

## Recommendation

Stop building context models on the gene data. The drug curve says six contexts is the harmful
regime and that drug contexts cannot stand in for gene contexts; the official scores say every
composite gain here is NMAE/reach calibration concentrated on RPE1, with PDS flat or worse, and
the only context-specific effect (+0.026 correct-over-wrong) does not beat refit at any strength.

**Next experiment (one, cheap, decisive):** test whether scratch's +0.018 is pure magnitude
calibration. Apply a *target-blind* per-line shrink/scale of the refit delta on RPE1 and HepG2
(the same family of operation the residual performs, with no context input), score it on the
full preset with the same bundles and three seeds, and compare to scratch a0.25 (RPE1 +0.073,
HepG2 +0.022). If it reproduces the gain, the context model has zero value beyond calibration and
that calibration belongs in the refit baseline, where it must then be checked against the actual
leaderboard before any upload, because the +0.0214 local composite gain of the PCA-interaction
model became −0.0020 officially. If it does not reproduce the gain, the remaining +0.004 outside
RPE1 is the true size of what six gene contexts can teach, and it is below the seed noise.

**Protocol notes.** No official upload was made or authorised. These training runs were logged
to `curves.json`, not to W&B, contrary to the 2026-09-18 rule; any follow-up must use the W&B
trainers. The five-line preset scorer needs ~16 GiB per run and must be run one line at a time on
this host (`CRASH_RECOVERY_2026-09-20.md`).

## Postscript — official test of the scratch residual, 2026-09-21

The user authorised one upload of the top local arm. `reports/submission_2026-09-21/FINDINGS.md`:
refit + 0.25 x scratch FiLM residual scored **0.137164** against the refit baseline's 0.137665,
**−0.000501**, PDS −0.0014, every metric except Jaccard slightly worse. The +0.0180 local composite
gain did not transfer. Refit remains the best official entry; the recommendation above stands.

## CAVEAT added 2026-09-21 — gene-side context vectors were scrambled

The gene contexts' control profiles used here (`contexts.npz` rows `gene:*`, built by
`dcp_prepare.py` from `cellfm_context_2026-09-20/controls_common`) carry sorted gene names over
source-order counts for HEPG2, JURKAT, K562, RPE1 and KOLF; only HCT116 and HEK293T are correct.
See `reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`. Consequences for this report:
§2–3 (drug side, Tahoe/GSE profiles built from their own count matrices) stand. §3b and §4–7 used
gene context vectors that were effectively random keys: mutually distinct (so correct-vs-wrong
context separation survives), but unrelated to the drug-context manifold, which by itself
produces "gene contexts 2.3x outside the training manifold", the gate collapsing to refit, and
the absence of drug→gene transfer. **The gene-transfer conclusions (decision items 3, 6, 7, 12,
13 in part) are not established by this experiment.** The official S3 result (postscript above)
used a model trained on these vectors. A corrected gene-context test with table-derived control
profiles is the Feng+CD4 experiment in `reports/ext_contexts_2026-09-21/`.
