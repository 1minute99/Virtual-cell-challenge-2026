#!/usr/bin/env python
"""Task 2 aggregation - {predictor} x {cell line} x {metric} table, histogram of
per-perturbation Pearson-delta, and the fraction of perturbations above 0.7.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task2_aggregate.py
"""
import glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")


def main():
    files = sorted(glob.glob(os.path.join(OUT, "task2_perpert_*.csv")))
    if not files:
        raise SystemExit("no task2 per-perturbation files yet")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df["predictor"] = df["predictor"].str.replace(r"^source_copy_.*", "source_copy", regex=True)
    print(f"loaded {len(files)} files, {len(df)} rows")

    metrics = [c for c in df.select_dtypes("number").columns
               if c not in ("seed",)]
    tab = (df.groupby(["source", "target_cell_line", "predictor", "target_gene_dropped"])[metrics]
             .mean().reset_index())
    tab.to_csv(os.path.join(OUT, "task2_metric_table.csv"), index=False)

    if "pearson_delta" in df.columns:
        frac = (df.groupby(["source", "target_cell_line", "predictor", "target_gene_dropped"])
                  ["pearson_delta"]
                  .agg(mean="mean", median="median",
                       frac_above_0p7=lambda s: float((s > 0.7).mean()),
                       n="size").reset_index())
        frac.to_csv(os.path.join(OUT, "task2_pearson_delta_summary.csv"), index=False)
        print(frac.to_string(index=False))

        sub = df[~df["target_gene_dropped"].astype(bool)]
        pairs = sub.groupby(["source", "target_cell_line"]).ngroups
        fig, axes = plt.subplots(1, max(pairs, 1), figsize=(4.2 * max(pairs, 1), 3.4),
                                 squeeze=False)
        for ax, ((s, t), g) in zip(axes[0], sub.groupby(["source", "target_cell_line"])):
            for pred, gg in g.groupby("predictor"):
                ax.hist(gg["pearson_delta"].dropna(), bins=40, alpha=.55, label=pred)
            ax.axvline(0.7, color="k", ls="--", lw=1)
            ax.set_title(f"{s} -> {t}", fontsize=9)
            ax.set_xlabel("per-perturbation Pearson-delta")
        axes[0][0].set_ylabel("perturbations"); axes[0][0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "task2_pearson_delta_hist.png"), dpi=150)
        print("wrote task2_pearson_delta_hist.png")

    piv = tab.pivot_table(index=["source", "target_cell_line", "predictor"],
                          values=[m for m in ["pearson_delta", "mae", "mse",
                                              "discrimination_score_l1", "overlap_at_N",
                                              "de_direction_match", "pr_auc"]
                                  if m in tab.columns])
    print(piv.to_string())
    piv.to_csv(os.path.join(OUT, "task2_headline_metrics.csv"))


if __name__ == "__main__":
    main()
