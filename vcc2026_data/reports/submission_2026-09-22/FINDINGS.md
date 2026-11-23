# S4: extended refit + 0.25 x FiLM residual (38 gene contexts) — official result 2026-09-22

Authorised by the user ("해봐"). One upload. Build `scripts/ext_submit_build.py` (K562 fold
reproduced to 6e-8), chain `ext_submit_chain.sh`, upload `ext_submit_upload.sh`; receipts
`build_receipt.json`, `upload_receipt.json`; package `S4_filmv2.vcc` sha256 `4b7c0496…`;
out-of-distribution caveat recorded beforehand in `OOD_NOTE.md`.

| | model | entry | overall | rank |
|---|---|---|---|---|
| previous best | S1 refit baseline (b20-wk7f3a) | `3CFK0NSTG9Fw6SfjD3IC` | 0.137665 | 231 |
| S3 | refit + scratch FiLM (scrambled profiles) | `5BTNATrQRIsC3XvRQayZ` | 0.137164 | 241 |
| **S4** | **extended refit + FiLM residual (b23-4b7c04)** | `df2PmtmQLtThacyTKu00` | **0.138975** | 265* |

*rank at retrieval; the board moved between submissions.

| metric | S1 refit | S4 | S4 − S1 |
|---|---|---|---|
| overall | 0.137665 | **0.138975** | **+0.001310** |
| pds | 0.545108 | 0.551176 | +0.006068 |
| lfc nmae | 0.075883 | 0.077782 | +0.001899 |
| fidelity | 0.018175 | 0.016888 | −0.001287 |
| reach | 0.193821 | 0.195572 | +0.001751 |
| jaccard | −0.006997 | −0.007567 | −0.000570 |

## Reading

**The project's best official score, and the first context-conditioned model to beat the refit on
the leaderboard.** The local five-line composite gain of +0.0247 became +0.0013 officially: the
direction transferred and the magnitude shrank ~19x, the same pattern as refit-vs-b19 (+0.0047
local PDS → +0.0012 official). Unlike the two earlier context residuals (PCA interaction −0.0020,
scrambled-profile FiLM −0.0005) the sign held. Notably the largest official gain is in **PDS**
(+0.0061), the target-discriminating metric, whereas locally PDS was the one metric that fell
(−0.014); NMAE and reach gained as predicted, fidelity and Jaccard slipped.

Three things changed together relative to S3 and cannot be separated by this one entry: the refit
base is the extended table with CD4/Feng nested (+0.0060 locally on its own), the FiLM's context
profiles are correct (S3's were scrambled for 5/7 training biologies), and 38 gene contexts
instead of 6 (the local scaling curve says this one did not matter). A/B/C sit well outside the
training context cloud (`OOD_NOTE.md`), so the residual extrapolated and still helped.

Two daily slots: one used on 2026-09-22. `official_scores.csv` updated.
