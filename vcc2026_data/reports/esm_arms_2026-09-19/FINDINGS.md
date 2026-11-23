# Gene embeddings and projection width do not matter — 2026-09-19

Requested: ESM-2 15B vs ESM-C 300M/600M, each at 64 / 128 / 256 dimensions. Added: centred
variants, random/zero null controls, and narrower 32/16 widths. Local reconstructed-count PDS,
five held-out lines, seeds 0/1/2, frozen scoring. No submission.

## Headline

**Random Gaussian vectors score as well as a 15B-parameter protein language model.**
ESM-2 15B 0.7958, random noise of the same shape **0.7968** — random ranks 2nd of 13 arms.
Every embedding/width/centring arm sits inside a band four times narrower than the seed noise
floor. The gene branch uses its input as an arbitrary per-gene identifier, not as biology.

## Finding 0 (before any arm ran): embeddings cannot touch the refit baseline

`SourceBaseline.predict` overwrites the feature-based ridge with each target's empirical mean for
every target seen in training, and all 300 official / 279 KOLF evaluation targets are seen in
every fold. Permuting or zeroing the features changes the refit prediction by **exactly 0.000**,
re-confirmed in-flight: refit was bit-identical (0.882070 on K562) across all 13 arms. So these
arms necessarily test the b19 residual path.

## Results

| arm | five-line b19 | vs baseline | KOLF |
|---|---|---|---|
| esm2_15b_128 | 0.7940 | -0.0017 | 0.5667 |
| esm2_15b_256 | 0.7946 | -0.0012 | 0.5686 |
| esm2_15b_256_centred | 0.7947 | -0.0010 | 0.5709 |
| esm2_15b_64 | 0.7958 | +0.0000 | 0.5690 |
| esm2_15b_64_centred | 0.7952 | -0.0005 | 0.5646 |
| esmc_300m_128 | 0.7946 | -0.0011 | 0.5677 |
| esmc_300m_256 | 0.7966 | +0.0008 | 0.5722 |
| esmc_300m_64 | 0.7962 | +0.0004 | 0.5711 |
| esmc_300m_64_centred | 0.7960 | +0.0002 | 0.5700 |
| esmc_600m_128 | 0.7960 | +0.0002 | 0.5674 |
| esmc_600m_256 | 0.7934 | -0.0024 | 0.5676 |
| esmc_600m_64 | 0.7974 | +0.0016 | 0.5694 |
| random_64 | 0.7968 | +0.0010 | 0.5700 |

13 arms, spread 0.0040, seed noise ~0.010

ESM-C 600M went +0.0016 / +0.0002 / -0.0024 across 64/128/256 — the best AND the worst arm in the
grid. No dose-response with model scale; the +0.0016 was the top of a noise band.

## Why: two mechanisms, both measured

**1. The gene branch treats its input as a hash.** The trained first layer puts comparable weight
per input unit on the embedding columns and the baseline-coordinate columns (ratios 0.80-1.13), so
the embeddings are read, not ignored. But random vectors work equally well, so what is used is
gene DISTINCTNESS, not gene biology.

**2. The context branch is a 58-value source lookup, and it dominates.** From training receipts
(15 fits per arm, stable to 3 decimals): context branch rms 0.308-0.332 vs gene-only branch
0.154-0.161, i.e. **context supplies 0.667-0.674 of the residual**. But on the KOLF fold there are
only **58 distinct control profiles across 92,029 rows, one per source** (58/58 pure). The table
stores one pooled untreated profile per context and broadcasts it, so the context input is
constant within a source: a categorical source label, not a continuous covariate.

Together these explain every null in this project: two thirds of the correction keys on source
identity, the rest on arbitrary gene identity, and neither transfers to an unseen line.

## The shared-direction hypothesis was wrong

Raw inter-gene cosine is +0.650 (ESM-2) and +0.763 (ESM-C); the adapter preserves it at every
width because a random projection preserves geometry. Centring drops it to +0.02 / +0.001. I
predicted centring would unlock a difference. It did not (-0.0005, -0.0010, +0.0002). The
adapter was not the bottleneck; the features simply carry nothing usable.

## The fix this points to (feasible with existing data)

Per-batch controls were collapsed at assembly and are still recoverable:
RN HepG2 **56 batches**, RN Jurkat **55**, Replogle GWPS **585 control rows**, VCC2025 **144** —
all pooled to one profile each. Restoring them would make the context input vary WITHIN a source,
so source identity would no longer determine it and the branch would have to learn a continuous
mapping. Caveats: per-batch controls are noisier, and batch variation is technical rather than the
cell-type variation we must generalise across.

## Assets

`processed/esmc/esmc_300m.pt` (960D) and `esmc_600m.pt` (1152D), 19,084 genes, mean-pooled over
residues, truncated at 2,046. Sequences: UniProt human reviewed, matched by `GN=`, 96.4% coverage.
ESM-C installed in isolated venv `.venv-esmc` — do NOT install `esm` into ct_transfer, it moves torch.
Caveat: Arc's ESM-2 map pooling is undocumented, so part of any ESM-2/ESM-C gap could be pooling —
moot now that random vectors match both.

Validation: the ESM-2 64D arm reproduces frozen b19 EXACTLY on all five lines (0.7958 vs 0.7958).
