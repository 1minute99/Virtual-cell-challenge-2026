# B19 sweep 2 completed

Local reconstructed-count PDS on previously studied held-out panels; not official scores. No official upload.

## Gates
- Arm A (weight decay x batch size): passed=False gain=-0.0003215342064464723 line_wins=2 paired=7
- Arm B (K=3 seed ensemble): passed=True gain=0.0012078749845809108 line_wins=4 paired=10
- Arm B (K=5 seed ensemble): passed=False gain=4.082128622264019e-05 line_wins=3 paired=9

## Arm A inner landscape (seed 0, 10 scopes x 2 inner lines)

| weight_decay | batch | steps | n | mean PDS | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|---|---|---|---|
| 0.003 | 256 | 2500 | 20 | 0.7722 | 0.8521 | 0.7846 | 0.8091 | 0.8394 | 0.5756 |
| 0.03 | 256 | 2500 | 20 | 0.7720 | 0.8518 | 0.7850 | 0.8078 | 0.8407 | 0.5748 |
| 0.03 | 64 | 10000 | 20 | 0.7717 | 0.8511 | 0.7845 | 0.8061 | 0.8452 | 0.5716 |
| 0.3 | 1024 | 625 | 20 | 0.7712 | 0.8488 | 0.7853 | 0.7946 | 0.8469 | 0.5804 |
| 0.003 | 64 | 10000 | 20 | 0.7703 | 0.8529 | 0.7830 | 0.8063 | 0.8401 | 0.5691 |
| 0.3 | 256 | 2500 | 20 | 0.7703 | 0.8531 | 0.7826 | 0.7985 | 0.8429 | 0.5742 |
| 0.3 | 64 | 10000 | 20 | 0.7701 | 0.8513 | 0.7809 | 0.7985 | 0.8459 | 0.5736 |
| 0.03 | 1024 | 625 | 20 | 0.7697 | 0.8495 | 0.7808 | 0.7950 | 0.8447 | 0.5784 |
| 0.003 | 1024 | 625 | 20 | 0.7694 | 0.8491 | 0.7803 | 0.7953 | 0.8442 | 0.5781 |

## Outer five-line mean PDS (seeds 0-2 unless noted)

| tag | mean | n | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|---|---|
| b16 | 0.754619 | 15 | 0.7861 | 0.7419 | 0.8136 | 0.8354 | 0.5961 |
| b19_reference | 0.795777 | 15 | 0.8835 | 0.8103 | 0.8564 | 0.8597 | 0.5689 |
| ensemble_k3 | 0.796985 | 15 | 0.8849 | 0.8115 | 0.8598 | 0.8618 | 0.5669 |
| ensemble_k5 | 0.795818 | 15 | 0.8842 | 0.8118 | 0.8547 | 0.8614 | 0.5671 |
| mean | 0.736193 | 15 | 0.7604 | 0.8073 | 0.7472 | 0.8114 | 0.5547 |
| no_change | 0.500000 | 15 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| reference_extra_seed | 0.794541 | 60 | 0.8828 | 0.8103 | 0.8537 | 0.8595 | 0.5663 |
| refit | 0.800538 | 15 | 0.8821 | 0.8107 | 0.8846 | 0.8596 | 0.5658 |
| selected | 0.795455 | 15 | 0.8835 | 0.8147 | 0.8509 | 0.8608 | 0.5673 |

## Arm B: ensemble size vs single seed

| K | mean PDS | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|---|
| 1 | 0.795777 | 0.8835 | 0.8103 | 0.8564 | 0.8597 | 0.5689 |
| 3 | 0.796985 | 0.8849 | 0.8115 | 0.8598 | 0.8618 | 0.5669 |
| 5 | 0.795818 | 0.8842 | 0.8118 | 0.8547 | 0.8614 | 0.5671 |

Single-seed PDS standard deviation over 15 reference seeds per line: K562 0.0029, RPE1 0.0022, HepG2 0.0104, Jurkat 0.0043, KOLF 0.0040

## Arm C: six-metric local diagnostic (seed 0, mean over 5 lines; raw cell-eval2 values, not calibrated overall)

| tag | de_wilcoxon_direction_fidelity_yield_raw | de_wilcoxon_direction_reach_raw | de_wilcoxon_lfc_nmae | de_wilcoxon_sig_jaccard | expr_mse_unbiased_capped_norm | pds_cosine |
|---|---|---|---|---|---|---|
| b19_reference | 0.6078 | 0.4542 | 0.8391 | 0.0438 | 3.6909 | 0.7956 |
| selected | 0.6083 | 0.4642 | 0.8385 | 0.0441 | 3.6709 | 0.7969 |
| refit | 0.6079 | 0.4400 | 0.8433 | 0.0422 | 3.7802 | 0.8005 |
| ensemble_k3 | 0.6069 | 0.4542 | 0.8392 | 0.0441 | 3.6915 | 0.7967 |
| ensemble_k5 | 0.6070 | 0.4558 | 0.8389 | 0.0439 | 3.6783 | 0.7969 |

Sweep-1 weight parity checks: 45 bit-exact of 45 comparable reference fits.

Per-line six-metric rows: `full_metric_comparison.json`. Selections: `selections.json`, `deployment_choice.json`. Controls: `context_controls.json`, `target_shuffled_controls.json`.
