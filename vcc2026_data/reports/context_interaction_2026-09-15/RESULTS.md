# Context interaction experiment — 2026-09-15

The context model is implemented and tested. It does **not** improve the general denoiser on this benchmark: whole-cell-line-average PDS 0.77834981 versus 0.77901524. It wins only 3/6 cell-line averages against that denoiser. The original official b16 result remains the best recorded challenge score, 0.12919304. None of these new experiments has been submitted.

## Architecture and fitting

For each allowed training set, the model starts with a b16-style signed covariance correction fitted from those sources only. Each training destination receives a separate denoiser prediction fitted without that destination. The new model regresses the remaining target-specific residual on `[z, z x h]`: `z` contains 8 or 16 components of the source mean perturbation response, and `h` contains up to two components of untreated expression. An output basis learned from training responses decodes 8 or 16 components back into the gene space. The fitted regression has 192 or 768 coefficients when all components are available. All bases, feature scales and coefficients exclude the outer cell line. A positive ridge penalty constrains the correction. No context-only intercept can create the same correction for every target.

The no-context ablation regresses on `z` alone. The grid uses ranks 8/16 and ridge strengths 0.1/1/10, with separate inner whole-line selections for the context and no-context families. The denoiser is also tested at correction strengths 0/0.5/1/1.5. The outer line's responses cannot enter any fitting or hyperparameter selection. Training targets are shared across Jiang lines; this experiment does not claim unseen-target validation.

## Data and evaluation

All four usable Jiang stimuli are included: IFNB, INS and TNFA have six lines; TGFB has five. There are 23 condition-by-line outer folds, summarized with equal weight for each of six cell lines. Both biological replicates are retained for evaluation and averaged for fitting. Controls are replicate matched; integer raw counts are normalized over the complete source library before restricting features. All screened target-gene coordinates are excluded. This run uses the official feature intersection and inherits the audited balanced target panels from the existing benchmark.

These are **CP10K pseudobulk PDS**, not reconstructed-cell competition scores. Every outer prediction receives 32 target derangements. A held-out-line result is a biological fold; cells, targets, stimuli and reconstruction seeds do not create additional independent cell-line folds. This is an exploratory expansion of an already studied benchmark. The new grid and source hashes were frozen in [spec.json](spec.json) before this run's outer results.

## Results

| Method | PDS | Replicate 1 PDS | Replicate 2 PDS |
|---|---:|---:|---:|
| Mean transfer | 0.736738 | 0.710197 | 0.722201 |
| Previously selected Jiang model | 0.778250 | 0.749370 | 0.762632 |
| Refitted b16-style denoiser | 0.779015 | 0.749229 | 0.763686 |
| Context interaction, inner-selected | 0.778350 | 0.750550 | 0.762754 |
| No-context residual, inner-selected | 0.780934 | 0.752549 | 0.766259 |
| No-context residual, matched hyperparameters | 0.780821 | 0.750729 | 0.765344 |
| Context model, zero query context | 0.780485 | 0.750074 | 0.765151 |
| Context model, wrong query context average | 0.775869 | 0.747218 | 0.760500 |
| Training context assignments shuffled | 0.778289 | 0.748391 | 0.762880 |
| New residual training labels shuffled | 0.777183 | 0.747630 | 0.763516 |
| Inner-selected correction strength | 0.781047 | 0.752205 | 0.767099 |

Correct query context beats the average wrong query context in 5/6 cell-line averages, but it loses to the zero-context ablation in 5/6. Shuffling training context assignments changes the average PDS by only about 0.00006. Thus response to a wrong input does not establish a beneficial learned context function. The context model is not promoted. The no-context residual's small gain also requires broader validation and is not proof of context-dependent supervision.

The residual-label shuffle randomizes only labels for the **new incremental correction**; the general denoiser remains fixed. Its score is expected to retain the denoiser's target information, so it is not a chance-level null for the full prediction. Prediction-target derangements provide the separate target-information floor. Context-assignment and residual-label controls each use one fixed training permutation; they are diagnostics, not exhaustive permutation significance tests.

## Correction strength and broader validation

Jiang favors a stronger correction by a small amount: PDS 0.779015 at current strength and 0.781222 at 1.5x. However, [the broader source benchmark](../broad_strength_2026-09-15/summary.json) reverses this: current 0.610668 versus 1.5x 0.569491, a loss in all four held-out biologies. It uses 1,661 shared targets across HCT116, HEK293T, CD4 and Feng iPSC, with every training-panel target coordinate additionally removed. The operator is refitted without the scored biology. Heterogeneous study protocols remain a limitation.

The frozen official operator also shows why "stronger denoising" needs care: increasing its correction to 1.5x reverses eight eigen-directions of its learned input subspace, while the current strength reverses none. These are response-pattern directions, not individual genes; the measurement alone does not prove the cause of the transfer loss. [Operator diagnostic](operator_diagnostics.json).

The exact reconstructed KOLF strength comparison is complete: four-seed mean PDS is 0.165490 for mean transfer, 0.186620 at 0.5x correction, 0.230554 at current strength, and 0.257194 at 1.5x. The stronger correction wins all four paired reconstruction seeds, with average gain 0.026640 and average ratio above its own shuffle of 1.1312 versus current strength. The four-biology loss prevents promotion. [Complete paired results](../correction_strength_2026-09-15/paired_seeds.csv). It uses paired seeds 0–3 and a separate target-shuffled prediction for each candidate, with the same existing scale 2. The original reconstruction arithmetic/RNG was extracted into a shared generator. Bulk normalization and PDS use the installed scorer. A retained complete prediction gives **zero bulk error and identical per-target PDS** before new settings are evaluated. KOLF is diagnostic; its result does not override the four-biology reversal or predict the official overall score.

## Verification and artifacts

Four scientific tests pass: unseen-response poisoning, recovery of a known synthetic context interaction, no-context invariance, and serialization parity. In real IFNB/A549 data, replacing the held-out responses with random values x1e6 leaves all 17 candidate predictions bit-identical. A separate regression test compares the refactored full builder with the saved pre-refactor builder cell for cell, under both stochastic and nearest rounding, including unchanged genes. [Measured checks](invariants.json).

[Per-fold metrics](fold_metrics.csv), [per-target ranks](target_metrics.csv), [inner selections](selection.csv), [cell-line comparisons](comparison_by_cell_line.csv), [summary](summary.json), trained fold models/predictions, and [source snapshots](sources/manifest.json) preserve the experiment. No new submission or upload was launched.

## Research follow-up

A recent response-decomposition preprint reports that untreated expression is more informative about shared responses than about all perturbation-specific components in its studied screens, and motivates gene-coessentiality features. This is a hypothesis for our next input change, not a proof that our context task is impossible. [Molina and Zhang preprint](https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1.full).

The authors' cross-cell-line code maps DepMap gene embeddings to expression responses with ridge/MLP models; its true target response centroid is used for evaluation and must not be copied into prediction inputs. The inspected repository tree lacks the embedding data advertised in its README. [Author code](https://github.com/xinyizhanglab/perturbation-decomposition/blob/main/benchmarks/cross_cl.py), [data availability audit](research_data_audit.json).

TRADE separately motivates effect uncertainty modelling rather than treating cell count as gene-wise precision. The next uncertainty experiment should use real replicate/donor variability and shared controls. The previous exact 0.16-to-0.90 claim remains unverified. [TRADE](https://www.nature.com/articles/s41588-025-02169-3).


## Acquired input for the next model

A versioned public DepMap 23Q4 gene-fitness matrix has now been downloaded, checked against the publisher's file sizes and MD5 hashes, and converted into a feature table with missing values preserved. It contains 18,442 canonical genes across 1,100 screened cell models. Canonical gene mapping covers **300/300 official targets**, each measured in at least 324 cell models. Of b16's 1,661 training targets, 1656 are present and 1644 have at least 100 measured models. Poorly covered genes require an explicit fallback; missing fitness measurements must not be interpreted as measured zero effects.

[Local data manifest](../../processed/depmap_23q4/manifest.json), [versioned publisher record](https://doi.org/10.25452/figshare.plus.24667905.v2). These are CRISPR knockout fitness measurements in cancer models, not CRISPRi transcriptomic responses. This is a deliberate versioned prior, not a claim to the latest release. No response model has been trained on this new input yet.

The next supervised test can learn conserved target-specific response programs from these external gene features, with evaluated targets excluded from response-label fitting and whole destination biologies held out. Compare direct prediction and mixtures with the denoiser; include shuffled gene-feature/response correspondence. This provides an input change after the failed control-expression interaction experiment, while keeping the successful b16 reference.

[Paired biological comparison figure](comparison.png) is also available as [PDF](comparison.pdf).
