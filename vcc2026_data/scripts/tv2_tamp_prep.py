"""Per-target cross-family AGREEMENT for the sum-side amplitude rule (harness token `ta<beta>`, FINDINGS §23.8):
agree_t = mean pairwise cosine of the donor FAMILY mean LFCs (c 50, centred, min_cells 20; families measuring >= 500 genes of t),
families = the harness fusion families (fh_harness.W_NEW > 0) minus the context's exclusions (KOLF for somatic contexts).
Same arithmetic as sch_emit 'ta' / ml_amp_probe. Usage: tv2_tamp_prep.py OUT.npz CTX[:somatic|stem] ... [--targets FILE]
(targets default: the official pert_counts.csv order)."""
import sys, pathlib, numpy as np, pandas as pd
from cell_eval2.moments import jackknife_correction  # noqa
_a = list(sys.argv); sys.argv = ["x"]; sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_model as TM
FAMS = ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1", "HEPG2", "RNJURKAT", "RPE1"]   # fh_harness.W_NEW > 0

def fam_means(targets, fams, min_cells=20):
    out = {}
    for fam in fams:
        num = np.zeros((len(targets), TM.NG)); den = np.zeros((len(targets), TM.NG))
        for s_, f_ in TM.FAMILY.items():
            if f_ != fam or not TM._src_path(s_).exists():
                continue
            z = TM.load(s_); cen = np.where(np.isfinite(z["center_c50"]), z["center_c50"], 0.0)
            for k, t in enumerate(targets):
                i = z["pos"].get(TM.canon(t))
                if i is None or z["ncell"][i] < min_cells or TM.canon(t) in TM.DROP_TARGETS:
                    continue
                v = np.log((z["fp"][i].astype(np.float64) * 1e6 + 50) / (z["fc"][i].astype(np.float64) * 1e6 + 50)) - cen
                m = np.isfinite(v); num[k, m] += v[m]; den[k, m] += 1.0
        if den.any():
            out[fam] = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    return out

def agreement(Fm, fams, T):
    ag = np.full(T, np.nan)
    for k in range(T):
        cols = [f for f in fams if f in Fm and np.isfinite(Fm[f][k]).sum() > 500]; cs = []
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                a, b = Fm[cols[i]][k], Fm[cols[j]][k]; mm = np.isfinite(a) & np.isfinite(b)
                cs.append(float(a[mm] @ b[mm] / (np.linalg.norm(a[mm]) * np.linalg.norm(b[mm]) + 1e-12)))
        if cs: ag[k] = float(np.mean(cs))
    return ag

if __name__ == "__main__":
    out = pathlib.Path(_a[1]); args = _a[2:]; tfile = None
    if "--targets" in args:
        i = args.index("--targets"); tfile = args[i + 1]; args = args[:i] + args[i + 2:]
    if tfile:
        targets = pathlib.Path(tfile).read_text().split()
    else:
        targets = pd.read_csv("/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls/pert_counts.csv").iloc[:, 0].astype(str).tolist()
    res = {"targets": np.array(targets)}
    for a in args:
        ctx, _, kind = a.partition(":"); kind = kind or "somatic"
        ex = set(kind[5:].split(",")) if kind.startswith("excl=") else ({"KOLF"} if kind == "somatic" else set())
        fams = [f for f in FAMS if f not in ex]
        Fm = fam_means(targets, fams); ag = agreement(Fm, fams, len(targets))
        res[f"agree_{ctx}"] = ag
        print(ctx, kind, "families", fams, "targets with >= 2 families", int(np.isfinite(ag).sum()), "/", len(targets),
              "median agreement", round(float(np.nanmedian(ag)), 4), flush=True)
    np.savez(out, **res); print("TAMP_PREP_DONE", out)
