# Independent denoising validation

The earlier context-free MoE is a linear denoising correction. Its expert weights are constant within each fitted model: `sum_k a_k (p W_k) = p (sum_k a_k W_k)`. It has no destination-state input. This track evaluates that incremental correction separately from context learning, preserving all historical reports and checkpoints.

## Jiang reproduction and seed repeats

| model | b16_pds | selected_pds | gain_vs_b16 | seed_pds_sd | wins | ties | label_shuffle_pds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| moe | 0.779015 | 0.782063 | 0.003047 | 0.000000 | 5 | 0 | 0.776077 |
| linear | 0.779015 | 0.782063 | 0.003047 | 0.000000 | 5 | 0 | 0.776077 |
| ridge | 0.779015 | 0.782063 | 0.003047 | nan | 5 | 0 | 0.776077 |

All 23 stimulus-by-line folds use the unchanged historical preprocessing, source-only denoiser, ranks 8/16 and inner-selected strengths 0/.25/.5/1. The MoE repeats its independently selected no-context settings, rather than settings selected for its context model. Fresh seed0 training under the original four-thread numerical setting reproduces every historical selected score and rank/strength, including overall PDS 0.7820626596083255. All three implementations obtain the same PDS, with 5/6 pooled-truth line wins and 5/6 lines improved on both replicate scores. Seeds are repeat optimizations, not new biological replicates. Equal discrete PDS does not imply bit-identical continuous predictions.

The direct linear map uses the same effective initial weights as the identically seeded MoE and the same 160-step AdamW objective. Ridge solves the common MSE plus 0.5 residual-energy objective analytically with a fixed 0.001 trace-scaled coefficient penalty. All three use the same inner rank/strength grid. This establishes that the Jiang gain does not require expert routing or a nonlinear network.

The local PDS gain does not improve every agreement metric. MoE row cosine is about 0.048624 versus b16's 0.049274, and absolute log-magnitude error is 0.718059 versus 0.705776 (higher is worse). RMSE improves slightly, 0.031326 versus 0.031418. All of these measurements are in the historical Jiang pseudobulk space, not reconstructed competition counts.

Compute elapsed time: 91.3 seconds; unique cached fits: 560. Detailed fold, per-target, magnitude/error/cosine metrics, inner selections, predictions and checkpoints are retained in this dataset's directory.

## Broader biological transfer

| model | b16_pds | selected_pds | gain_vs_b16 | seed_pds_sd | wins | ties | label_shuffle_pds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| moe | 0.609403 | 0.606299 | -0.003104 | 0.000001 | 0 | 2 | 0.609300 |
| linear | 0.609403 | 0.605879 | -0.003523 | 0.000000 | 0 | 2 | 0.609304 |
| ridge | 0.609403 | 0.606305 | -0.003098 | nan | 0 | 2 | 0.609300 |

Each destination biology is entirely excluded from feature fitting, baseline fitting, residual fitting and inner selection. The target panel is independently partitioned into three fixed folds; the evaluated target fold is also excluded from all fitting and inner selection. Its response observations in permitted source cells are available only as query inputs. This tests new model-fitting targets in an unseen cell, with source observations available. It does not test prediction for a perturbation with no observed source response. The three target folds do not increase the number of independent biological contexts.

The four biological groups are HCT116, HEK293T, primary CD4 (all states pooled), and Feng iPSC. The historical same-target reference is independently replayed exactly at mean PDS 0.610667836910556 before the stricter target-disjoint analysis. The 1,661-target panel excludes the official and KOLF evaluation targets; all benchmark target-gene coordinates are excluded. There are 16,121 measured evaluation genes. No separate replicate labels are available in these pooled profiles.

Compute elapsed time: 397.6 seconds; unique cached fits: 228. Detailed fold, per-target, magnitude/error/cosine metrics, inner selections, predictions and checkpoints are retained in this dataset's directory.

## Independent Replogle/Nadig genetic panel

| model | b16_pds | selected_pds | gain_vs_b16 | seed_pds_sd | wins | ties | label_shuffle_pds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| moe | 0.912828 | 0.912536 | -0.000291 | 0.000000 | 1 | 1 | 0.912737 |
| linear | 0.912828 | 0.912537 | -0.000291 | 0.000000 | 1 | 1 | 0.912737 |
| ridge | 0.912828 | 0.912695 | -0.000132 | nan | 2 | 0 | 0.912750 |

Each destination biology is entirely excluded from feature fitting, baseline fitting, residual fitting and inner selection. The target panel is independently partitioned into three fixed folds; the evaluated target fold is also excluded from all fitting and inner selection. Its response observations in permitted source cells are available only as query inputs. This tests new model-fitting targets in an unseen cell, with source observations available. It does not test prediction for a perturbation with no observed source response. The three target folds do not increase the number of independent biological contexts.

The pilot contains 1,744 shared genetic targets and 4,693 measured evaluation genes across K562, RPE1, HepG2 and Jurkat. All five local scPerturb files were aggregated consistently from count sums. Controls were matched by experimental batch, weighted according to each target's cells across batches. Normalization used the full measured gene library before evaluating its intersection. Targets require at least 25 cells in each primary line. Official and KOLF target panels were excluded. Batch annotations are not assumed to be biological replicates. The separately retained K562 essential screen is one additional screen of the same biology. [scPerturb describes its count matrices](https://www.nature.com/articles/s41592-023-02144-y); the named files are distributed through the [public release](https://zenodo.org/records/13350497).

The denoiser itself transfers strongly on this panel: mean-transfer PDS is 0.818842049 versus 0.912827613 for the permitted-source b16-style model. The additional MoE correction has negative mean gain, so this finding supports the denoising reference, not the extra correction. K562 overlaps the Jiang cell-line roster; the new studies and other three lines supply additional evidence, but the dataset line counts must not be summed as distinct global biologies.

Compute elapsed time: 169.0 seconds; unique cached fits: 228. Detailed fold, per-target, magnitude/error/cosine metrics, inner selections, predictions and checkpoints are retained in this dataset's directory.

## Verification and decision

A post-fit K562 screen check evaluates the identical predictions on 499, 503 and 477 eligible test targets in the three target folds. MoE and the direct linear model are identical to b16 there because inner selection disables their corrections. Ridge gains about 0.000656 PDS on the GWPS screen and 0.000594 on the essential-gene screen; this small one-biology gain does not overturn its negative four-line mean. The b16-style improvement over mean transfer persists in both screens. These are two screens of one cell biology, not two independent destination contexts.

Fresh real IFNB/A549 poisoning checks preserve all unshrunk predictions for MoE, linear and ridge at both ranks. Low-rank score algebra agrees with full gene-space PDS. The installed cell-eval2 scorer was checked separately. New broad preparation passes whole-biology poisoning and an inference-target mutation check: changing held-out-target source observations cannot alter training labels, bases or scales. All frozen reference hashes remain unchanged. MoE's numerical collapse to a single matrix is verified for saved models.

An initial two-thread strict replay was stopped after a one-rank-step inner-score discrepancy; its files and log are retained in `jiang_strict_replay_attempt/`. Repeating with the original four-thread setting reproduces the historical results exactly. The completed Jiang spec's numeric `threads: 4` is authoritative; two descriptive strings carried over from the diagnostic attempt mention two threads and are clarified in `jiang/metadata_clarification.json`.

The denoising lead does not pass the predefined broader-transfer gate. Retain b16. No KOLF reconstruction or official submission was started. The modest Jiang gain is real locally but is insufficient evidence for a new challenge candidate.

These are compact, exploratory comparisons with only four to six independent biological contexts. Negative transfer does not prove that every linear correction or context-free model must fail; it shows that this particular Jiang-selected lead is not yet robust enough to promote. This track makes no claim about destination-cell context, drug-to-gene transfer, or unseen perturbations lacking source observations.
