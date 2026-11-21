"""transfer_v2 source library: per-source expression FRACTIONS for the targets we need, from raw sums.

For each source s and target t we store f_pert[t, g] = S_t,g / T_t (full-library denominator) and the
matched control fraction f_ctrl[t, g], on the official 18,533-gene axis (NaN = gene not measured by s).
LFC at fusion time: log((f_pert*1e6 + c) / (f_ctrl*1e6 + c)) with a CPM pseudocount c, i.e. a TRUE
fold change of the expression rate (not a log1p(CP10K) difference), so it can be applied
multiplicatively to any destination's own control expression.

Also stores, per source and per pseudocount c, the mean LFC over EVERY target the source measured
(the source's generic/shared response, for source centering).

Sources that measure 0/300 official targets (RPE1, HepG2, Jurkat, K562-essential) are left out of the DEFAULT
list on purpose: they would help the H1 proxy but cannot help the A/B/C panel. They can still be requested by name
(RN_HepG2 / RN_Jurkat / RN_RPE1) - the final phase has a NEW target panel that may well overlap their 896 targets,
and uncovered targets cost ~0.007 composite per 10 % (FINDINGS 2026-09-27 section 13).
"""
import json, pathlib, sys
import numpy as np, pandas as pd
from scipy import sparse

R = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data")
OUT = R / "reports/transfer_v2_2026-09-23/sources"
sys.path.append(str(R / "scripts"))
import vcclib as V

OFF = pathlib.Path("/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls")
GENES = pd.read_csv(OFF / "gene_names.csv").iloc[:, 0].astype(str).values
NG = len(GENES)
GM = V.GeneMap(); canon = GM.canon
GIC = {}
for i, g in enumerate(GENES):
    GIC.setdefault(canon(g), i)
CS = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]


def needed_targets():
    """Official val panel + H1 proxy panel (+ any panel listed in $TV2_TARGETS, e.g. the final-phase panel)."""
    import os
    off = pd.read_csv(OFF / "pert_counts.csv")["target_gene"].astype(str).tolist()
    if os.environ.get("TV2_TARGETS"):
        off += pathlib.Path(os.environ["TV2_TARGETS"]).read_text().split()
    import anndata as ad
    h1 = ad.read_h5ad(R / "reports/transfer_v2_2026-09-23/h1/real_thin.h5ad", backed="r")
    h1t = sorted(set(h1.obs["target"].astype(str)) - {"non-targeting"})
    h1.file.close()
    return {canon(t) for t in off} | {canon(t) for t in h1t}


def lfc(fp, fc, c):
    return np.log((fp * 1e6 + c) / (fc * 1e6 + c))


class Acc:
    """Collect per-target fractions for needed targets + streaming mean LFC over all targets."""

    def __init__(self, need):
        self.need = need; self.fp = {}; self.fc = {}; self.n = {}
        self.cs = {c: np.zeros(NG) for c in CS}; self.cc = {c: np.zeros(NG) for c in CS}

    def add(self, key, fp, fc, ncell):
        for c in CS:
            v = lfc(fp, fc, c); m = np.isfinite(v)
            self.cs[c][m] += v[m]; self.cc[c][m] += 1
        if key in self.need:
            if key in self.fp:   # duplicate rows for one target (aliases): cell-weighted merge
                w0, w1 = self.n[key], ncell
                self.fp[key] = (self.fp[key] * w0 + fp * w1) / (w0 + w1)
                self.fc[key] = (self.fc[key] * w0 + fc * w1) / (w0 + w1)
                self.n[key] += ncell
            else:
                self.fp[key] = fp.astype(np.float32); self.fc[key] = fc.astype(np.float32); self.n[key] = ncell

    def save(self, name, meta):
        keys = sorted(self.fp)
        np.savez_compressed(OUT / f"{name}.npz", targets=np.array(keys),
                            fp=np.vstack([self.fp[k] for k in keys]).astype(np.float32),
                            fc=np.vstack([self.fc[k] for k in keys]).astype(np.float32),
                            ncell=np.array([self.n[k] for k in keys]),
                            **{f"center_c{c:g}": np.where(self.cc[c] > 0, self.cs[c] / np.maximum(self.cc[c], 1), np.nan) for c in CS},
                            **{f"ncenter_c{c:g}": self.cc[c] for c in CS})
        meta.update(n_needed=len(keys), n_all=int(max(self.cc[CS[0]].max(), 0)))
        (OUT / f"{name}.json").write_text(json.dumps(meta, indent=1))
        print(name, meta, flush=True)


def load_store(store):
    p = R / "pseudobulk" / store
    o = pd.read_csv(p / "obs.csv", low_memory=False)
    z = np.load(p / "sum_counts.npz", allow_pickle=True)
    M = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    o2s = np.load(p / "official_to_source.npy")
    return o, M, o2s


def official_row(M, r, o2s):
    s = np.asarray(M[r].todense()).ravel().astype(np.float64)
    v = np.full(NG, np.nan); ok = o2s >= 0
    v[ok] = s[o2s[ok]]
    return v


def from_store(name, store, need, filt=None, batch_matched=False):
    o, M, o2s = load_store(store)
    sel = np.ones(len(o), bool)
    if filt:
        col, val = filt
        sel = (o[col].astype(str) == val).values
    is_c = o["is_control"].astype(str).str.lower().isin(["true", "1"]).values
    acc = Acc(need)
    if batch_matched:
        bcol = "technical_batch"
        ctrl = {}
        for r in np.flatnonzero(sel & is_c):
            b = str(o[bcol].iloc[r]); v = official_row(M, r, o2s)
            T = float(o["total_counts"].iloc[r])
            if b in ctrl:
                ctrl[b] = (ctrl[b][0] + v, ctrl[b][1] + T)
            else:
                ctrl[b] = (v, T)
        tgt = o["target_symbol"].astype(str).values
        rows = np.flatnonzero(sel & ~is_c)
        by_t = {}
        for r in rows:
            by_t.setdefault(canon(tgt[r]), []).append(r)
        for t, rs in by_t.items():
            S = np.zeros(NG); T = 0.0; fc = np.zeros(NG); n = 0
            for r in rs:
                b = str(o[bcol].iloc[r])
                if b not in ctrl:
                    continue
                v = official_row(M, r, o2s); Tr = float(o["total_counts"].iloc[r])
                S = S + v; T += Tr; fc = fc + Tr * ctrl[b][0] / ctrl[b][1]; n += int(o["n_cells"].iloc[r])
            if T <= 0:
                continue
            acc.add(t, S / T, fc / T, n)
    else:
        rc = np.flatnonzero(sel & is_c)
        Sc = np.zeros(NG); Tc = 0.0
        for r in rc:
            Sc = Sc + official_row(M, r, o2s); Tc += float(o["total_counts"].iloc[r])
        fc = Sc / Tc
        tgt = o["target_symbol"].astype(str).values
        rows = np.flatnonzero(sel & ~is_c)
        by_t = {}
        for r in rows:
            by_t.setdefault(canon(tgt[r]), []).append(r)
        for t, rs in by_t.items():
            S = np.zeros(NG); T = 0.0; n = 0
            for r in rs:
                S = S + official_row(M, r, o2s); T += float(o["total_counts"].iloc[r]); n += int(o["n_cells"].iloc[r])
            if T > 0:
                acc.add(t, S / T, fc, n)
    acc.save(name, dict(store=store, filt=filt, batch_matched=batch_matched))


def from_rn(name, line, need):
    z = np.load(R / f"reports/simple_residual_2026-09-16/data/replogle_nadig/{line}/pseudobulk.npz", allow_pickle=True)
    gpos = np.array([GIC.get(canon(str(g)), -1) for g in z["genes"]])
    ok = gpos >= 0
    S = z["sum_counts"]; T = z["source_full_library_sums"]
    C = z["matched_control_sums"]; TC = z["matched_control_full_library_sums"]
    acc = Acc(need)
    for i, t in enumerate(z["targets"]):
        fp = np.full(NG, np.nan); fc = np.full(NG, np.nan)
        fp[gpos[ok]] = S[i][ok] / T[i]; fc[gpos[ok]] = C[i][ok] / TC[i]
        acc.add(canon(str(t)), fp, fc, int(z["target_cell_counts"][i]))
    acc.save(name, dict(store=f"replogle_nadig/{line}", matched_controls="target-batch"))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    need = needed_targets()
    print("needed targets", len(need), flush=True)
    which = sys.argv[1:] or ["K562", "HCT116", "HEK293T", "CD4_Rest", "CD4_Stim8hr", "CD4_Stim48hr", "KOLF", "FENG_T", "FENG_GWF", "FENG_GWN", "H1"]
    for w in which:
        if w == "K562": from_rn("K562", "K562", need)
        elif w.startswith("RN_"): from_rn(w, w[3:], need)     # RN_HepG2 / RN_Jurkat / RN_RPE1 (final-phase coverage)
        elif w in ("HCT116", "HEK293T"): from_store(w, "orion_target_cellline_all", need, ("cell_context", w))
        elif w.startswith("CD4_"): from_store(w, "cd4_target_condition_all", need, ("stimulation", w[4:]))
        elif w == "KOLF": from_store(w, "kolf_target", need)
        elif w == "FENG_T": from_store(w, "feng_targetedscreen_target", need)
        elif w == "FENG_GWF": from_store(w, "feng_genomewidescreen_fitnessgenes_target", need)
        elif w == "FENG_GWN": from_store(w, "feng_genomewidescreen_nonfitnessgenes_target", need)
        elif w == "H1": from_store(w, "vcc2025_target_batch", need, batch_matched=True)
