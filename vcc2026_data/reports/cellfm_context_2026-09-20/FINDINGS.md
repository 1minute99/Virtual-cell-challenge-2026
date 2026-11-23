# PCA vs STATE SE-600M vs UCE as unseen-biology context, on the real harness metric

Date: 2026-09-20. Scripts: `cellfm_controls.py`, `cellfm_common_axis.py`, `cellfm_encode.py`,
`cellfm_harness_eval.py`, `cellfm_report.py`. Baseline: frozen five-line refit mean **0.8005**.

## Answer in one line

**No.** Neither Cell FM beats PCA, neither beats the refit baseline by more than noise, and the
one arm that looks most consistent (UCE) scores the same with its context coordinates **shuffled**.

## Setup

Both models run through their OWN published inference code and published weights, frozen:

| model | entry point | weights | dim |
|---|---|---|---|
| STATE SE-600M | `state emb transform` (arc-state 0.11.1) | `se600m_epoch16.ckpt` + `protein_embeddings.pt` from `arcinstitute/SE-600M` | 2058 |
| UCE (33-layer) | `eval_single_anndata.py` (snap-stanford/UCE, unmodified) | `33l_8ep_1024t_1280.torch`, `all_tokens.torch`, `species_chrom.csv`, `species_offsets.pkl`, `protein_embeddings/` from the paper's figshare record | 1280 |

Nothing was gated; no vocabulary or preprocessing was reconstructed. Two documented deviations:
UCE pins `torch==2.1.1`, which predates the RTX 5090's compute capability, so it runs on torch
2.11 with its own code and files untouched; and its hardcoded `model_files/` path is satisfied by
symlinks to the published files rather than by editing its source.

`prediction = refit_baseline + alpha * interaction(context, perturbation)`, **interaction only**.
The generic per-biology shift is excluded by design; it caused the -0.21 HepG2 collapse in
`reports/ctx_harness_2026-09-20`, and excluding it worked -- HepG2's worst arm here is -0.0079.

**Seven biologies, not nine.** CD4 ships only pseudobulk and FENG has no non-targeting control
arm (its only non-perturbed label is `unassigned`, a QC category), so neither can supply control
cells. Both are dropped from EVERY arm so the comparison stays matched. Six training contexts
per fold, context dim 5.

**One common gene axis.** Native panels span 8,248 to 38,584 genes in three assay groups; encoding
those directly would let a FM separate biologies by assay panel rather than biology. All seven
matrices are cut to the 6,413-gene intersection and the PCA arm is rebuilt from the same 2,000
cells on the same genes, so the arms differ only in how they encode identical input.

**Validation.** alpha=0 reproduces the stored refit per seed to ~1e-5 (K562 seed 0: 0.881962 vs
0.881951; five-line 0.8005 vs 0.8005).

## Results, five-line mean PDS

| arm | K562 | RPE1 | HepG2 | Jurkat | KOLF | five-line | vs refit |
|---|---|---|---|---|---|---|---|
| refit (alpha=0) | 0.8821 | 0.8107 | 0.8846 | 0.8596 | 0.5658 | **0.8005** | +0.0000 |
| **uce correct a=1** | 0.8862 | 0.8118 | 0.8857 | 0.8626 | 0.5733 | **0.8039** | **+0.0034** |
| **pca correct a=1** | 0.8823 | 0.8318 | 0.8766 | 0.8638 | 0.5602 | **0.8030** | **+0.0024** |
| uce shuffled a=1 | 0.8668 | 0.8266 | 0.8895 | 0.8600 | 0.5741 | 0.8034 | +0.0029 |
| **se correct a=1** | 0.8818 | 0.8086 | 0.8878 | 0.8602 | 0.5678 | **0.8012** | **+0.0007** |
| pca wrong a=1 | 0.8748 | 0.8178 | 0.8886 | 0.8513 | 0.5626 | 0.7990 | -0.0015 |
| se wrong a=1 | 0.8708 | 0.8105 | 0.8882 | 0.8509 | 0.5607 | 0.7962 | -0.0043 |
| se shuffled a=1 | 0.8755 | 0.7830 | 0.8882 | 0.8603 | 0.5678 | 0.7949 | -0.0056 |
| uce wrong a=1 | 0.8669 | 0.8046 | 0.8862 | 0.8493 | 0.5584 | 0.7931 | -0.0075 |
| pca shuffled a=1 | 0.8801 | 0.8036 | 0.8746 | 0.8596 | 0.5674 | 0.7971 | -0.0035 |

Paired per (line, seed), 15 pairs, against the ±0.010 seed floor:

| arm | mean | sd | min | max | wins | clears floor |
|---|---|---|---|---|---|---|
| uce correct a=1 | +0.0034 | 0.0025 | +0.0009 | +0.0094 | 15/15 | no |
| pca correct a=1 | +0.0024 | 0.0107 | -0.0082 | +0.0212 | 9/15 | no |
| se correct a=1 | +0.0007 | 0.0019 | -0.0023 | +0.0036 | 9/15 | no |

Nested alpha selection (alpha for line L chosen on the other four only):

| rep | nested five-line | gain | lines positive |
|---|---|---|---|
| uce | 0.8039 | +0.0034 | 5/5 |
| se | 0.8008 | +0.0002 | 3/5 |
| pca | 0.8000 | **-0.0006** | 3/5 |

## The decisive ablation

| rep | correct − wrong | correct − shuffled | verdict |
|---|---|---|---|
| pca | +0.0039 (9/15) | **+0.0059 (12/15)** | context-specific |
| se | +0.0050 (9/15) | **+0.0063 (9/15)** | context-specific but gain ≈ 0 |
| uce | +0.0108 (12/15) | **+0.0005 (6/15)** | **NOT context-specific** |

UCE is the only arm that is positive on all five lines and all fifteen paired seeds, which looks
like the strongest result in the table. It is not. Shuffling the five context coordinates leaves
the gain essentially unchanged (+0.0029 vs +0.0034, correct beating shuffled on only 6 of 15
pairs). A shuffle preserves the vector's norm and the multiset of its values but destroys which
biology it identifies, so UCE's gain tracks the **magnitude** of the correction, not the identity
of the context. This is the same failure mode recorded earlier where ablations order by residual
size rather than accuracy.

PCA is the only representation whose gain is genuinely context-specific *and* non-trivial, and it
is driven by one line: removing RPE1 (+0.0211) turns its five-line mean from +0.0024 to **-0.0023**.
UCE without its best line still gives +0.0023, but that survives only because the effect is a
magnitude artifact that no line resists.

## The five questions

1. **Does either Cell FM beat PCA?** No. UCE is +0.0009 over PCA winning 9 of 15 pairs, a tie.
   STATE SE is -0.0017 below PCA. A 600M-parameter cell model and a 1,280-dim UCE embedding do
   not beat a PCA of the pooled control expression profile.
2. **Does either beat refit by more than noise?** No. Best is UCE at +0.0034 against a ±0.010
   floor. Its paired spread is tight (sd 0.0025, 15/15 wins) so it is not random, but it is not a
   context effect either, per the shuffle test.
3. **Does correct context clearly beat wrong/shuffled on PDS?** For PCA yes, on both. For STATE SE
   yes on both but from a gain of +0.0007, so the specificity has nothing to act on. For UCE it
   beats wrong but **not shuffled**, which is the test that matters and which it fails.
4. **Consistent or driven by one or two lines?** PCA's entire gain is RPE1; without it PCA is
   negative. STATE SE is ~0 everywhere. UCE is uniform across all five lines, but uniformly for
   the wrong reason.
5. **Is Cell-FM context worth pursuing on the real metric?** **No.** Neither model beats a PCA of
   the same cells, neither clears the noise floor, and the better of the two fails its own shuffle
   control. Frozen whole-cell embeddings of control cells do not carry usable unseen-biology
   perturbation-response signal at this sample size.

## A bug I introduced and corrected

The first run whitened each context representation per component after PCA. With n training
biologies projected onto n-1 components this is exactly degenerate: the whitened Gram matrix
becomes `n*(I - 11^T/n)` for ANY input, so all three representations' training configurations
collapse onto the same simplex and the arms cannot differ. The symptom was `pca_wrong`, `se_wrong`
and `uce_wrong` all scoring exactly 0.8207 on RPE1. Verified numerically at d = 6413, 2058, 1280,
then replaced with a single global scale, which equalises magnitude across representations for the
shared ridge penalty while preserving each one's geometry. The degenerate run is kept in
`degenerate_whitening_run/`. Every number above is from the corrected run.

## Caveats that limit this

- Six training context vectors per fold. That is the binding constraint on every arm, and it is
  why the honest conclusion is "not at this sample size" rather than "these models are useless".
- 2,000 control cells per biology, mean-pooled. A different pooling (attention, per-cell-type) was
  not tried.
- Only the frozen models were used, as specified. No fine-tuning, no adapter.
- Restricting to a 6,413-gene common axis costs depth for the wide-panel biologies (KOLF's median
  UMI drops from 5,130 to 843), which may handicap the FMs relative to their native input. The
  alternative confounds biology with assay panel, which is worse.

## Recommendation

Stop pursuing frozen Cell-FM context on this data. The refit baseline at **0.8005** remains the
best local predictor and still has never been scored officially. The constraint is the number of
independent biologies, not the richness of the context encoder: a 600M-parameter model given the
same seven biologies extracts no more usable context than a PCA does. Adding independent biologies
remains the highest-value investment.
