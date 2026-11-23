# S3: refit + 0.25 x scratch FiLM context residual — official result 2026-09-21

Authorised by the user ("submit it"). One upload, one of the two daily slots. Build:
`scripts/dcp_submit_build.py` (verified against the stored K562 fold to 2.9e-5),
`scripts/dcp_submit_chain.sh`, `scripts/dcp_submit_upload.sh`; receipts `build_receipt.json`,
`upload_receipt.json`, package `S3_scratch.vcc` sha256 `cbe57f70…`.

| | model | entry | overall | rank |
|---|---|---|---|---|
| best official | S1 refit baseline (b20-wk7f3a) | `3CFK0NSTG9Fw6SfjD3IC` | **0.137665** | 231 |
| S2 | refit + PCA interaction (b21-wk2e9c) | `hlhJKNFjK8UHE1dq2rsV` | 0.135702 | 241 |
| **S3** | **refit + 0.25 x scratch FiLM residual (b22-cbe57f)** | `5BTNATrQRIsC3XvRQayZ` | **0.137164** | 241 |

| metric | S1 refit | S3 | S3 − S1 |
|---|---|---|---|
| overall | 0.137665 | 0.137164 | **−0.000501** |
| pds | 0.545108 | 0.543674 | −0.001434 |
| lfc nmae | 0.075883 | 0.074979 | −0.000904 |
| fidelity | 0.018175 | 0.017514 | −0.000661 |
| reach | 0.193821 | 0.193777 | −0.000044 |
| jaccard | −0.006997 | −0.006957 | +0.000040 |
| expr mse | 0 | 0 | 0 |

## Reading

The local five-line official-preset composite said +0.0180 over refit (10/15 seeds, 80 % from
RPE1). The leaderboard says **−0.0005**: a small loss, every metric except Jaccard slightly worse.
The residual for A/B/C was 6–9 % of the refit norm, 95–97 % correlated across the three contexts
and 73–84 % correlated with the zero-context output, i.e. essentially one generic shift; A/B/C sit
near HEK293T in the context space (nn distance 8.6–9.3 against 23.2 between training biologies).

This is the fourth local-versus-official pair for context-conditioned residuals and the pattern
holds: cosine benchmarks → harness PDS → local full preset → leaderboard, and the gain shrinks or
reverses at each step. Combined with the drug-context experiment, the residual direction on gene
data is closed with an official number: **the refit baseline (0.137665) remains the best entry.**

Local diagnostics that would have predicted this: the full-data model early-stopped at its first
checkpoint (step 100) and did not beat the zero residual on the inner biology (0.239 vs 0.196),
exactly like every leave-one-biology-out fold, and the local gain vanished without RPE1 (+0.0043).
