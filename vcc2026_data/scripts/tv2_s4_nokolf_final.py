"""PANEL-AGNOSTIC version for the final phase (2026-09-29): reads S4_ctx{c}.delta.npy + S4_refit_v2.delta.npy from $TV2_S4_OUT
(the tv2_s4_rebuild.py --full output for the panel in $VCC_OFFICIAL_DIR), contexts $TV2_CTX, and writes S4noK_ctx{c} there.
Verified to reproduce the A/B/C S4noK files. Original description:
S4 without the KOLF biology for the official contexts (S10k sum side; 2026-09-27). Reproduces the saved S4 refit_v2 base
first (same code path as tv2_s4_rebuild.full), then refits on table_v2 without KOLF rows and takes the empirical means
without KOLF; the FiLM context residual is the SAVED one (S4_ctx - S4_refit_v2), so only the biology set changes.
Writes reports/submission_2026-09-27/S4noK_ctx{A,B,C}.delta.npy + S4noK_refit_v2.delta.npy + S4noK_build.json."""
import os, sys, json, time, pathlib
assert "TV2_S4_OUT" in os.environ and "VCC_OFFICIAL_DIR" in os.environ, "set TV2_S4_OUT (S4 rebuild dir) and VCC_OFFICIAL_DIR"
sys.path.insert(0, "/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import numpy as np
import tv2_s4_rebuild as RB
S4 = pathlib.Path(os.environ["TV2_S4_OUT"]); O = S4; CTXS = os.environ.get("TV2_CTX", "A,B,C").split(",")
t0 = time.time(); base = RB.A.Data(); table2 = RB.D.Table.load(RB.ROOT / "table_v2")
indices = {g: i for i, g in enumerate(RB.V.official_features())}; columns = np.array([indices[g] for g in table2.genes])
targets = RB.V.official_targets(); feats = base.features(targets); canon = [RB.D.canonical_target(t) for t in targets]
ids2 = RB.F.permitted_ids(table2, [])
refit2, _ = RB.refit_on(table2, ids2, targets, feats, canon); emp2, nv = RB.XB.empirical(targets, [], "v2")
base2 = emp2.copy(); base2[:, columns] = refit2
saved = np.load(S4 / "S4_refit_v2.delta.npy").astype(np.float64)
rep = dict(max_abs=float(np.abs(base2 - saved).max()), cos=float((base2 * saved).sum() / np.linalg.norm(base2) / np.linalg.norm(saved)))
print(json.dumps(dict(event="reproduce_saved_refit_v2", **rep, seconds=round(time.time() - t0, 1))), flush=True)
assert rep["max_abs"] < 1e-3, rep
bios = np.asarray(table2.bios).astype(str); keep = ids2[bios[ids2] != "KOLF"]
refitk, _ = RB.refit_on(table2, keep, targets, feats, canon); empk, nvk = RB.XB.empirical(targets, ["KOLF"], "v2")
basek = empk.copy(); basek[:, columns] = refitk; np.save(O / "S4noK_refit_v2.delta.npy", basek.astype(np.float32))
out = dict(rows_all=int(len(ids2)), rows_noK=int(len(keep)), voters_outside_all=nv, voters_outside_noK=nvk,
           cos_base_noK_vs_saved=float((basek * saved).sum() / np.linalg.norm(basek) / np.linalg.norm(saved)), reproduce=rep, ctx={})
for c in CTXS:
    S = np.load(S4 / f"S4_ctx{c}.delta.npy").astype(np.float64); res = S - saved
    outside = np.ones(S.shape[1], bool); outside[columns] = False; assert np.abs(res[:, outside]).max() < 1e-5
    Sk = basek + res; assert Sk.shape == (len(targets), 18533) and np.isfinite(Sk).all()
    np.save(O / f"S4noK_ctx{c}.delta.npy", Sk.astype(np.float32))
    out["ctx"][c] = dict(cos_vs_S4=float((Sk * S).sum() / np.linalg.norm(Sk) / np.linalg.norm(S)), norm_ratio=float(np.linalg.norm(Sk) / np.linalg.norm(S)))
json.dump(out, open(O / "S4noK_build.json", "w"), indent=1); print(json.dumps(out), flush=True)
