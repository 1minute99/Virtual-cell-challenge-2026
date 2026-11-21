"""RQ3 Task 16 -- H3a verdict + summary from the rq3_energy.py runs.  -> results/rq3_h3a.json

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md (H3a, its falsification criterion, operational
definition 7, Task 16). Clarification (knowledge/_meta/log.md, 2026-11-05, "Pre-data operationalisation for Task 16"):
  (1) H3a is decided with the registered TARGET bootstrap (B = 2,000, mean over split seeds 0-4) on H1, k = 8; every
      panel also reports the two-level bootstrap (targets x control cells), used for all secondary-panel statements.
  (2) Negative-control rule in NUMERATOR form: D_neg = (pseudo-target numerator excess per target, N_neg) divided by
      the OBSERVED real-data denominator per target; H3a is rejected if D_neg >= 0.5 x observed D (either D_other
      or D_perm).
ASSUMPTIONS
  S1 Primary input = pred 'prem3_ctrlAB' (Task 14 C_a u C_b rebuild per split seed), space 'gcap' (G-cap, 6,014).
  S2 Verdict rule (registered): REJECTED if the target-bootstrap 95 % CI of D_other or of D_perm at k = 8 is not
     entirely above 0, OR the negative control gives D_neg >= 0.5 x D_obs for D_other or D_perm (point estimates).
     Otherwise CONFIRMED. (The plan names no 'inconclusive' branch for H3a.)
  S3 The negative-control comparison uses the H1 negctrl run with the same splits, subspaces and k.
  S4 Cross-panel counts use the two-level percentile CI (lo > 0); target-bootstrap counts are reported alongside.
  S5 Two-level CIs are conservative on small-control panels (synthetic calibration, Task 16 assumption note);
     target-only CIs are anti-conservative there. Both are reported; neither changes the H1 decision rule.
No randomness here (all resampling happened in rq3_energy.py; seeds listed there).
"""
from __future__ import annotations
import json, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402

R = common.RESULTS
KS = (8, 24, 50)
QS = ("F_dest", "F_other_mean", "F_perm_mean", "F_haar", "D_other", "D_perm")


def row(df, panel, q, k, boot="target", stratum="all"):
    x = df[(df.panel == panel) & (df.quantity == q) & (df.k == k) & (df.boot == boot) & (df.stratum == stratum)]
    if len(x) != 1:
        raise KeyError((panel, q, k, boot, stratum, len(x)))
    return x.iloc[0]


def ci(df, panel, q, k, stratum="all"):
    t = row(df, panel, q, k, "target", stratum)
    out = dict(point=float(t.point), target_ci=[float(t.lo), float(t.hi)], target_se=float(t.se))
    if (df.boot == "two_level").any():
        w = row(df, panel, q, k, "two_level", stratum)
        assert np.isclose(w.point, t.point, rtol=1e-9, atol=1e-12)
        out.update(two_level_ci=[float(w.lo), float(w.hi)], two_level_se=float(w.se))
    return out


def main():
    E = pd.read_parquet(R / "rq3_energy.parquet")
    N = pd.read_parquet(R / "rq3_negctrl.parquet")
    prim = E[(E.pred == "prem3_ctrlAB") & (E.space == "gcap")]
    full = E[(E.pred == "prem3_ctrlAB") & (E.space == "full")]
    panels = [p for p in common.PANELS7 if p in set(prim.panel)]
    out = dict(plan="input/ideas/zero-shot-limits-context-complexity-plan.md", task=16, hypothesis="H3a",
               primary_input="results/preds/<panel>_prem3_ctrlAB_seed{0..4}.npz; G-cap (6,014 genes)",
               bootstrap="target: B = 2,000 over targets, seed-mean of split seeds 0-4 (registered); two_level: same "
                         "target draws x per-seed resampling of C_a and C_b cells (robustness; conservative on small "
                         "control thirds)",
               panels={}, negative_control={}, strata={}, full_space_secondary={})
    # ---- per panel x k ----
    for p in panels:
        d = {}
        for k in KS:
            dk = {q: ci(prim, p, q, k) for q in QS}
            dk["F_other_by_line"] = {q.split(":", 1)[1]: float(row(prim, p, q, k).point)
                                     for q in sorted(set(prim[prim.panel == p].quantity)) if q.startswith("F_other:")}
            dk["den_per_target"] = float(row(prim, p, "den_per_target", k).point)
            dk["n_targets"] = int(row(prim, p, "D_perm", k).n_targets)
            d[f"k{k}"] = dk
        out["panels"][p] = d
    # ---- negative control (numerator form, scaled by the observed denominator) ----
    for p in sorted(set(N.panel), key=common.PANELS7.index):
        if p not in panels:
            continue
        d = {}
        for k in KS:
            den_obs = float(row(prim, p, "den_per_target", k).point)
            dk = {"den_obs_per_target": den_obs,
                  "den_neg_per_target": float(row(N, p, "den_per_target", k).point)}
            for w in ("other", "perm"):
                n_ = row(N, p, f"N_{w}", k)
                d_obs = float(row(prim, p, f"D_{w}", k).point)
                d_neg = float(n_.point) / den_obs
                dk[w] = dict(D_neg_numerator_form=d_neg,
                             D_neg_ci=[float(n_.lo) / den_obs, float(n_.hi) / den_obs],
                             D_obs=d_obs, half_D_obs=0.5 * d_obs,
                             frac_of_obs=d_neg / d_obs if d_obs != 0 else None,
                             fails_50pct_rule=bool(d_neg >= 0.5 * d_obs),
                             rule_applicable=bool(d_obs > 0),   # with D_obs <= 0 the 50 % rule is vacuous (H3a fails anyway)
                             D_neg_ratio_form_unstable=float(row(N, p, f"D_{w}", k).point))
            d[f"k{k}"] = dk
        out["negative_control"][p] = d
    # ---- H3a verdict (H1, k = 8) ----
    h = out["panels"]["H1"]["k8"]
    ng = out["negative_control"]["H1"]["k8"]
    c_other = h["D_other"]["target_ci"][0] > 0
    c_perm = h["D_perm"]["target_ci"][0] > 0
    neg_ok = not (ng["other"]["fails_50pct_rule"] or ng["perm"]["fails_50pct_rule"])
    verdict = "confirmed" if (c_other and c_perm and neg_ok) else "rejected"
    out["H3a"] = dict(
        verdict=verdict,
        rule="Rejected if the 95 % target-bootstrap CI (B = 2,000; mean over 5 split seeds) of D_other or D_perm at "
             "k = 8 is not entirely above 0, or the control-only negative control (numerator form scaled by the "
             "observed denominator) yields D >= 50 % of the observed D.",
        deciding_numbers=dict(
            F_dest=h["F_dest"]["point"], F_other_mean=h["F_other_mean"]["point"], F_perm_mean=h["F_perm_mean"]["point"],
            F_haar=h["F_haar"]["point"],
            D_other=h["D_other"]["point"], D_other_target_ci=h["D_other"]["target_ci"],
            D_other_two_level_ci=h["D_other"].get("two_level_ci"),
            D_perm=h["D_perm"]["point"], D_perm_target_ci=h["D_perm"]["target_ci"],
            D_perm_two_level_ci=h["D_perm"].get("two_level_ci"),
            negctrl_D_other=ng["other"]["D_neg_numerator_form"], negctrl_D_perm=ng["perm"]["D_neg_numerator_form"],
            half_D_other=ng["other"]["half_D_obs"], half_D_perm=ng["perm"]["half_D_obs"]),
        criteria=dict(D_other_ci_above_0=c_other, D_perm_ci_above_0=c_perm, negctrl_below_half=neg_ok),
        robustness_two_level_ci_above_0=dict(
            D_other=h["D_other"]["two_level_ci"][0] > 0, D_perm=h["D_perm"]["two_level_ci"][0] > 0))
    # ---- cross-panel summary ----
    cp = {}
    for k in KS:
        rec = {}
        for b, key in (("two_level", "two_level_ci"), ("target", "target_ci")):
            pos = {q: [p for p in panels if out["panels"][p][f"k{k}"][q][key][0] > 0] for q in ("D_other", "D_perm")}
            rec[b] = dict(n_panels=len(panels), D_other_ci_above_0=pos["D_other"], D_perm_ci_above_0=pos["D_perm"],
                          both=[p for p in pos["D_other"] if p in pos["D_perm"]],
                          n_both=len([p for p in pos["D_other"] if p in pos["D_perm"]]))
        rec["point_positive_both"] = [p for p in panels if out["panels"][p][f"k{k}"]["D_other"]["point"] > 0
                                      and out["panels"][p][f"k{k}"]["D_perm"]["point"] > 0]
        nc = out["negative_control"]
        rec["negctrl_fails_50pct"] = {f"D_{w}": [p for p in nc if nc[p][f"k{k}"][w]["fails_50pct_rule"]] for w in ("other", "perm")}
        rec["negctrl_fails_50pct_where_D_obs_positive"] = {
            f"D_{w}": [p for p in nc if nc[p][f"k{k}"][w]["fails_50pct_rule"] and nc[p][f"k{k}"][w]["rule_applicable"]]
            for w in ("other", "perm")}
        rec["negctrl_frac_of_obs"] = {p: {w: nc[p][f"k{k}"][w]["frac_of_obs"] for w in ("other", "perm")} for p in nc}
        cp[f"k{k}"] = rec
    out["cross_panel"] = cp
    # ---- strata (predicted-amplitude quartiles, Q1 = smallest) ----
    for p in panels:
        meta = json.loads((R / "rq3" / f"real_{p}_prem3_ctrlAB.json").read_text())
        d = {"quartile_edges": meta.get("amp_quartile_edges"), "quartile_n": meta.get("amp_quartile_n")}
        for k in KS:
            d[f"k{k}"] = {s: {q: ci(prim, p, q, k, s) for q in ("F_dest", "D_other", "D_perm")}
                          for s in ("Q1", "Q2", "Q3", "Q4")}
        out["strata"][p] = d
    # ---- secondary: destination-only, full gene space ----
    for p in [p for p in common.PANELS7 if p in set(full.panel)]:
        meta = json.loads((R / "rq3" / f"real_{p}_prem3_ctrlAB_full.json").read_text())
        out["full_space_secondary"][p] = {"n_genes": meta["n_genes"],
                                          **{f"k{k}": {q: ci(full, p, q, k) for q in ("F_dest", "F_perm_mean", "F_haar", "D_perm")}
                                             for k in KS}}
    # ---- runtimes / provenance ----
    rt = {}
    for f in sorted((R / "rq3").glob("*.json")):
        m = json.loads(f.read_text()); rt[f.stem] = round(m.get("runtime_s", float("nan")), 1)
    out["runtime_s"] = rt
    (R / "rq3_h3a.json").write_text(json.dumps(out, indent=1))
    # ---- stdout summary ----
    print(f"H3a (H1, k=8): {verdict.upper()}")
    print(json.dumps(out["H3a"]["deciding_numbers"], indent=1))
    for k in KS:
        print(f"k={k}: two-level both>0 {cp[f'k{k}']['two_level']['n_both']}/{len(panels)} {cp[f'k{k}']['two_level']['both']}; "
              f"target both>0 {cp[f'k{k}']['target']['n_both']}; negctrl fails (D_obs>0) {cp[f'k{k}']['negctrl_fails_50pct_where_D_obs_positive']}")
    return out


if __name__ == "__main__":
    main()
