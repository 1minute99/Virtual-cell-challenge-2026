# Tahoe context-transfer pilot results

| family   |   correct |   no_context |   correct_minus_none |   correct_minus_wrong |   correct_minus_zero |   heldout_line_wins |   biological_secondary_gain | passes_ten_context_gates   |
|:---------|----------:|-------------:|---------------------:|----------------------:|---------------------:|--------------------:|----------------------------:|:---------------------------|
| mlp      |  0.998990 |     0.998953 |             0.000036 |              0.000035 |             0.000019 |                   1 |                   -0.000229 | False                      |
| gp_full  |  0.998990 |     0.998954 |             0.000036 |              0.000078 |             0.000031 |                   5 |                   -0.001107 | False                      |
| rrr      |  0.998968 |     0.998954 |             0.000014 |              0.000015 |             0.000017 |                   3 |                   -0.000214 | False                      |
| ridge    |  0.998957 |     0.998954 |             0.000003 |              0.000008 |             0.000003 |                   3 |                    0.000153 | False                      |
| metric   |  0.998953 |     0.998954 |            -0.000001 |              0.000009 |             0.000003 |                   3 |                   -0.000257 | False                      |
| gp       |  0.998951 |     0.998954 |            -0.000003 |              0.000034 |            -0.000003 |                   4 |                   -0.000427 | False                      |
| bilinear |  0.998949 |     0.998954 |            -0.000006 |              0.000008 |            -0.000006 |                   3 |                   -0.000160 | False                      |
| film     |  0.998912 |     0.998949 |            -0.000036 |             -0.000011 |            -0.000038 |                   2 |                   -0.000229 | False                      |

The primary panel contains 16 whole held-out cell lines and 150 distinct drugs. The independent biological repeat covers 91 identical conditions from plates 6/14. All models use the same source-only preprocessing and small nested sweep. Elapsed model time was 1530.3 seconds for 3600 fits. Complete metrics, per-line scores, label/context permutations, target permutations and seed records are in `evaluation/`. GP uncertainty is a latent posterior approximation and excludes baseline and omitted-gene uncertainty. The GP has two shared kernel hyperparameters plus the separately reported posterior coefficient table; those values must not be hidden when comparing complexity.

Success requires all ten predeclared context criteria. A positive raw PDS change from an additive Ridge model can reflect a generic cell-specific shift; `correct_minus_zero_centered_pds` tests whether destination context changes perturbation-specific predictions after removing the across-drug component. Deterministic seeds are not replication. The paired-line bootstrap is descriptive and does not correct the exploratory comparison across seven model families.

`variance_components.csv` estimates shared reproducible, cell-specific reproducible, and replicate/noise energy using only the documented independent plate pair. These signed finite-panel estimates are retained without clipping. The technical cell-count eligibility filter can favor less cytotoxic conditions; abundance was never supplied to the predictor.

This chemical benchmark does not yet answer whether drug-trained context rules transfer to genetic perturbations. No drug effect has been equated with a knockdown effect. No official reconstruction or upload occurred.
