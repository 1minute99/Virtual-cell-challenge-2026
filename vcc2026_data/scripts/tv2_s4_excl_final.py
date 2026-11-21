"""Parametrized S4 biology exclusion for the official / final contexts (2026-09-30; generalises tv2_s4_nokolf_final.py).
Env: TV2_S4_OUT (tv2_s4_rebuild --full dir with S4_ctx{c} + S4_refit_v2), VCC_OFFICIAL_DIR, TV2_CTX (A,B,C),
TV2_S4_EXCL (comma list of table_v2 biologies to drop from the refit rows AND the empirical means; default KOLF),
TV2_S4_TAG (output prefix, default S4noK), TV2_S4_DEST (output dir, default TV2_S4_OUT).
Same code path as tv2_s4_nokolf_final.py (saved FiLM residual S4_ctx - S4_refit_v2 kept, only the biology set changes);
with TV2_S4_EXCL=KOLF it must reproduce the shipped S4noK files bit-exactly."""
import os, sys, json, time, pathlib
assert "TV2_S4_OUT" in os.environ and "VCC_OFFICIAL_DIR" in os.environ, "set TV2_S4_OUT (S4 rebuild dir) and VCC_OFFICIAL_DIR"
sys.path.insert(0, "/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import numpy as np
import tv2_s4_rebuild as RB
S4 = pathlib.Path(os.environ["TV2_S4_OUT"]); O = pathlib.Path(os.environ.get("TV2_S4_DEST", str(S4))); O.mkdir(parents=True, exist_ok=True)
CTXS = os.environ.get("TV2_CTX", "A,B,C").split(",")
EXCL = [b for b in os.environ.get("TV2_S4_EXCL", "KOLF").split(",") if b]; TAG = os.environ.get("TV2_S4_TAG", "S4noK")
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
bios = np.asarray(table2.bios).astype(str); assert all(b in set(bios) for b in EXCL), EXCL; keep = ids2[~np.isin(bios[ids2], EXCL)]
refitk, _ = RB.refit_on(table2, keep, targets, feats, canon); empk, nvk = RB.XB.empirical(targets, list(EXCL), "v2")
basek = empk.copy(); basek[:, columns] = refitk; np.save(O / f"{TAG}_refit_v2.delta.npy", basek.astype(np.float32))
out = dict(rows_all=int(len(ids2)), rows_noK=int(len(keep)), voters_outside_all=nv, voters_outside_noK=nvk,
           cos_base_noK_vs_saved=float((basek * saved).sum() / np.linalg.norm(basek) / np.linalg.norm(saved)), reproduce=rep, ctx={})
for c in CTXS:
    S = np.load(S4 / f"S4_ctx{c}.delta.npy").astype(np.float64); res = S - saved
    outside = np.ones(S.shape[1], bool); outside[columns] = False; assert np.abs(res[:, outside]).max() < 1e-5
    Sk = basek + res; assert Sk.shape == (len(targets), 18533) and np.isfinite(Sk).all()
    np.save(O / f"{TAG}_ctx{c}.delta.npy", Sk.astype(np.float32))
    out["ctx"][c] = dict(cos_vs_S4=float((Sk * S).sum() / np.linalg.norm(Sk) / np.linalg.norm(S)), norm_ratio=float(np.linalg.norm(Sk) / np.linalg.norm(S)))
out["excl"] = EXCL; json.dump(out, open(O / f"{TAG}_build.json", "w"), indent=1); print(json.dumps(out), flush=True)
