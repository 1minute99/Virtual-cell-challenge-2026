"""S4 improvement, H1 panel (official chemistry, non-essential targets): leave-H1-out refit on the original and the
denoised tables for the harness H1 panel's 300 targets (refit only -- no FiLM residual exists for H1). Saves
reports/s4_improve_2026-09-27/H1/film_{orig,z4,z1}.npy (harness keys s4Yorig / s4Yz4 / s4Yz1, no energy matching)
and film_{z4,z1}m.npy (per-target energy matched to orig; keys s4Yz4m / s4Yz1m)."""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, "/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import ext_baselines as XB
P = json.loads((XB.D.ROOT / "reports/final_harness_2026-09-24/panels.json").read_text())
targets = P["H1"]["targets"]; OUT = XB.D.ROOT / "reports/s4_improve_2026-09-27/H1"; OUT.mkdir(parents=True, exist_ok=True)
base = XB.A.Data(); indices = {g: i for i, g in enumerate(XB.V.official_features())}
canon = [XB.D.canonical_target(t) for t in targets]
missing = [t for t in canon if t not in base.feature_map]
print("H1 targets", len(targets), "without features", len(missing), flush=True)
feats = np.stack([base.feature_map[t] if t in base.feature_map else np.zeros_like(next(iter(base.feature_map.values()))) for t in canon])
emp, nv = XB.empirical(targets, ["H1"], "v2")
out = {}
for tag, tdir in [("orig", XB.ROOT / "table_v2"), ("z4", XB.D.ROOT / "reports/s4_improve_2026-09-27/table_v2z4"), ("z1", XB.D.ROOT / "reports/s4_improve_2026-09-27/table_v2z1")]:
    t0 = time.time(); table = XB.D.Table.load(tdir); columns = np.array([indices[g] for g in table.genes])
    common, nrows, nbios, seen = XB.refit(table, "H1", targets, feats, canon)
    pred = emp.copy(); pred[:, columns] = common; assert np.isfinite(pred).all(); out[tag] = pred.astype(np.float64)
    np.save(OUT / f"film_{tag}.npy", pred.astype(np.float32))
    print(json.dumps(dict(tag=tag, rows=nrows, bios=nbios, all_seen=bool(seen), seconds=round(time.time() - t0, 1))), flush=True)
n0 = np.linalg.norm(out["orig"], axis=1)
for tag in ("z4", "z1"):
    n1 = np.linalg.norm(out[tag], axis=1); m = out[tag] * np.where(n1 > 0, n0 / np.maximum(n1, 1e-12), 1.0)[:, None]
    np.save(OUT / f"film_{tag}m.npy", m.astype(np.float32))
    print(tag, "cos vs orig", round(float((out[tag] * out["orig"]).sum() / np.linalg.norm(out[tag]) / np.linalg.norm(out["orig"])), 4), flush=True)
