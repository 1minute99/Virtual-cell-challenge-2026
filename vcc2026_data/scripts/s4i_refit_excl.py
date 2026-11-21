"""S4 without extra biologies (2026-09-27, data-quality question): leave-one-biology-out refit_v2 exactly as s4i_refit.py,
but with EXTRA biologies removed from the training rows (via ext_baselines.CONFLICT_GENE, no table copy), then the S4
recipe film = refit + 0.25 x (film_ext_a1 - refit_v1) with the saved FiLM residual. Writes
reports/s4_improve_2026-09-27/<line>/{refit_TAG,film_TAG}.npy (harness keys X<TAG>/Y<TAG>).
Usage: s4i_refit_excl.py TAG BIO[,BIO...] [LINES]"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, "/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import ext_baselines as XB
tag, extra = sys.argv[1], sys.argv[2].split(","); EMP = not tag.endswith("_refitonly")
PANEL_TARGETS = None                                           # set by the Orion wrapper (targets of a harness panel)
lines = sys.argv[3].split(",") if len(sys.argv) > 3 else ["K562", "HepG2", "Jurkat", "RPE1"]
EXT = XB.ROOT; OUT = XB.D.ROOT / "reports/s4_improve_2026-09-27"; tdir = EXT / "table_v2"
base = XB.A.Data(); indices = {g: i for i, g in enumerate(XB.V.official_features())}
table = XB.D.Table.load(tdir); columns = np.array([indices[g] for g in table.genes])
bios = set(np.asarray(table.bios).astype(str)); assert all(b in bios for b in extra), (extra, sorted(bios)[:20])
for line in lines:
    t0 = time.time(); key = XB.D.canonical_bio(line); orig = list(XB.CONFLICT_GENE.get(key, []))
    XB.CONFLICT_GENE[key] = orig + [b for b in extra if b not in orig]
    try:
        targets = PANEL_TARGETS or base.targets[line]; canon = [XB.D.canonical_target(t) for t in targets]; feats = base.features(targets)
        emp, nv = XB.empirical(targets, [line] + (extra if EMP else []), "v2")   # genes outside the table: empirical means
        common, nrows, nbios, seen = XB.refit(table, line, targets, feats, canon)
    finally:
        XB.CONFLICT_GENE[key] = orig
    pred = emp.copy(); pred[:, columns] = common
    (OUT / line).mkdir(parents=True, exist_ok=True); np.save(OUT / line / f"refit_{tag}.npy", pred.astype(np.float32))
    if (EXT / line / "film_ext_a1.npy").exists():
        res = np.load(EXT / line / "film_ext_a1.npy").astype(np.float64) - np.load(EXT / line / "refit_v1.npy").astype(np.float64)
        ref = np.load(EXT / line / "refit_v2.npy").astype(np.float64); cs = float((pred * ref).sum() / np.linalg.norm(pred) / np.linalg.norm(ref))
    else:                                                      # no saved FiLM residual for this line: refit only
        res = 0.0; cs = float("nan")
    film = pred.astype(np.float64) + 0.25 * res; np.save(OUT / line / f"film_{tag}.npy", film.astype(np.float32))
    print(json.dumps(dict(line=line, tag=tag, extra=extra, empirical_excl=EMP, rows=nrows, bios=nbios, seconds=round(time.time() - t0, 1),
                          cos_vs_saved_refit_v2=cs)), flush=True)
