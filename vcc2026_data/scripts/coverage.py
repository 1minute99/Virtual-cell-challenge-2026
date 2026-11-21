#!/usr/bin/env python3
"""Target and feature coverage of the official VCC 2026 panel, from acquired data only.

Coverage is reported in the four stages the instruction distinguishes, because
they are genuinely different claims:

  measured_feature   the gene is a column of that assay's readout
  designed_target    the gene is in that screen's guide library
  observed_target    cells carrying a guide for that gene exist in the data
  usable_response    enough cells, and a matched control on the same context,
                     to form a control-relative effect

Only evidence read out of files actually on this disk is counted. A target that
no acquired dataset reaches is reported as unresolved WITH the reason, and
"not in the data we hold" is kept distinct from "not accessible yet".

Outputs:
  metadata/observed_target_inventory.csv
  reports/target_coverage.csv
  reports/target_coverage_summary.csv
  reports/unresolved_targets.csv
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import vcclib as V

ROOT = V.ROOT
EFF = ROOT / "effects"
PB = ROOT / "pseudobulk"


def load_effect_sources():
    out = []
    for d in sorted(EFF.iterdir()):
        if not (d / "pairs.csv").exists():
            continue
        P = pd.read_csv(d / "pairs.csv")
        prov = json.loads((d / "provenance.json").read_text())
        mask = np.load(d / "measured_mask.npy")
        out.append((d.name, P, prov, mask))
    return out


def main():
    targets = V.official_targets()
    features = V.official_features()
    gm = V.GeneMap()
    tset = set(targets)
    canon2orig = {}
    for t in targets:
        canon2orig.setdefault(gm.canon(t), t)

    sources = load_effect_sources()
    print(f"{len(sources)} effect stores on disk: {[s[0] for s in sources]}")

    inv_rows, cov_rows = [], []
    for key, P, prov, mask in sources:
        ctx_cols = [c for c in ("cell_context", "stimulation") if c in P.columns]
        for _, r in P.iterrows():
            inv_rows.append(dict(
                store_key=key, dataset_id=r.dataset_id, experiment_id=r.experiment_id,
                cell_context=r.cell_context, stimulation=r.stimulation, time=r.time,
                donor=r.donor, perturbation_mode=r.perturbation_mode,
                observed_target=r.target_symbol, target_ensembl=r.target_ensembl,
                n_cells_pert=r.n_cells_pert, n_cells_control=r.n_cells_control,
                n_source_groups=r.get("n_source_groups", ""),
                qc_flags=r.get("qc_flags", ""),
                official_features_measured=int(mask.sum()),
                evidence="observed target with a matched control (usable_response)",
                count_semantics=prov.get("count_semantics", prov.get("quantity", "")),
                source_url=prov.get("source_url", "")))
        # per (store, context) coverage of the official panel
        for grp, sub in P.groupby(ctx_cols if ctx_cols else lambda _: 0):
            ctx = grp if isinstance(grp, tuple) else (grp,)
            present = set(sub.target_symbol)
            exact = tset & present
            via_alias = {canon2orig[c] for c in
                         ({gm.canon(x) for x in present} & set(canon2orig)) } - exact
            for t in targets:
                hit = t in exact
                alias = (not hit) and (t in via_alias)
                if not (hit or alias):
                    continue
                name = t if hit else t
                srcname = t if hit else gm.canon(t)
                row = sub[sub.target_symbol == srcname]
                cov_rows.append(dict(
                    target_gene=t, store_key=key, dataset_id=sub.dataset_id.iloc[0],
                    cell_context=ctx[0],
                    stimulation=ctx[1] if len(ctx) > 1 else sub.stimulation.iloc[0],
                    match="exact_symbol" if hit else "reviewed_alias",
                    evidence_stage="usable_response",
                    n_cells_pert=int(row.n_cells_pert.max()) if len(row) else 0,
                    n_cells_control=int(row.n_cells_control.max()) if len(row) else 0,
                    official_features_measured=int(mask.sum()),
                    target_locus_measured=bool(mask[features.index(t)])
                    if t in features else False,
                    perturbation_mode=sub.perturbation_mode.iloc[0],
                    available_form="control-relative effect from pseudobulk"))
    inv = pd.DataFrame(inv_rows)
    inv.to_csv(ROOT / "metadata/observed_target_inventory.csv", index=False)
    print(f"observed target inventory: {len(inv)} rows, "
          f"{inv.observed_target.nunique()} distinct observed targets "
          f"-> metadata/observed_target_inventory.csv")

    cov = pd.DataFrame(cov_rows)
    cov.to_csv(ROOT / "reports/target_coverage.csv", index=False)

    # per-source summary over the official panel
    srows = []
    for key, P, prov, mask in sources:
        sub = cov[cov.store_key == key]
        srows.append(dict(
            store_key=key, dataset_id=P.dataset_id.iloc[0],
            contexts="|".join(sorted(set(P.cell_context.astype(str)))),
            conditions="|".join(sorted(set(P.stimulation.astype(str)))),
            perturbation_mode=P.perturbation_mode.iloc[0],
            observed_targets_total=int(P.target_symbol.nunique()),
            official_targets_covered=int(sub.target_gene.nunique()),
            official_features_measured=int(mask.sum()),
            official_features_unmeasured=int((~mask).sum()),
            median_cells_per_covered_target=float(
                sub.groupby("target_gene").n_cells_pert.max().median()) if len(sub) else 0.0))
    S = pd.DataFrame(srows).sort_values("official_targets_covered", ascending=False)
    S.to_csv(ROOT / "reports/target_coverage_summary.csv", index=False)
    print()
    print(S.to_string(index=False))

    covered = set(cov.target_gene)
    union_feat = np.zeros(len(features), dtype=bool)
    for _, _, _, m in sources:
        union_feat |= m
    print(f"\nUNION over acquired data: {len(covered)}/300 official targets have a usable "
          f"control-relative effect; {int(union_feat.sum())}/{len(features)} official "
          f"features are measured by at least one acquired dataset")

    # unresolved: distinguish "absent from what we hold" from "not yet accessible"
    pending = json.loads((ROOT / "manifests/pending_sources.json").read_text()) \
        if (ROOT / "manifests/pending_sources.json").exists() else {}
    urows = []
    for t in targets:
        if t in covered:
            continue
        urows.append(dict(
            target_gene=t, canonical_symbol=gm.canon(t), ensembl_id=gm.ens(t) or "",
            status="NOT_IN_ACQUIRED_DATA",
            note=("no usable target/control pair in any dataset currently on disk; "
                  "this is not evidence of absence from public data"),
            searched_stores="|".join(sorted(s[0] for s in sources)),
            further_routes_pending=";".join(pending.get("datasets", []))))
    U = pd.DataFrame(urows)
    U.to_csv(ROOT / "reports/unresolved_targets.csv", index=False)
    print(f"unresolved: {len(U)}/300 -> reports/unresolved_targets.csv")
    if len(U):
        print("  " + ", ".join(U.target_gene.tolist()[:40]) +
              (" ..." if len(U) > 40 else ""))


if __name__ == "__main__":
    main()
