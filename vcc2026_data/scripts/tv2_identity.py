"""Cell-line identity inference for hidden contexts (final phase D/E/F; validated on A/B/C).

Reference: Human Protein Atlas cell line RNA (nTPM, ~1,200 lines). Query: pooled control pseudobulk of a
context (CP10K). Spearman correlation over the most variable reference genes (log1p), top-k lines reported
with HPA lineage metadata. Usage: tv2_identity.py CONTEXT_H5AD [CONTEXT_H5AD ...]
"""
import json, pathlib, sys, zipfile
import h5py, numpy as np, pandas as pd
from scipy.stats import rankdata

HPA = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/raw/hpa")
CACHE = HPA / "hpa_celline_log1p_ntpm.npz"


def reference():
    if CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        return z["M"], [str(g) for g in z["genes"]], [str(c) for c in z["lines"]]
    df = pd.read_csv(HPA / "rna_celline.tsv", sep="\t", usecols=["Gene name", "Cell line", "nTPM"])
    piv = df.pivot_table(index="Gene name", columns="Cell line", values="nTPM", aggfunc="mean").fillna(0.0)
    M = np.log1p(piv.values.astype(np.float32))
    np.savez_compressed(CACHE, M=M, genes=piv.index.values.astype(str), lines=piv.columns.values.astype(str))
    return M, list(piv.index.astype(str)), list(piv.columns.astype(str))


def lineage_table():
    p = HPA / "celline_description.tsv.zip"
    try:
        with zipfile.ZipFile(p) as z:
            name = [n for n in z.namelist() if n.endswith(".tsv")][0]
            d = pd.read_csv(z.open(name), sep="\t")
        cols = {c.lower(): c for c in d.columns}
        key = cols.get("cell line") or d.columns[0]
        keep = [c for c in d.columns if any(s in c.lower() for s in ("disease", "primary", "lineage", "organ", "cancer", "tissue", "cell type"))]
        return d.set_index(key)[keep]
    except Exception:
        return None


def context_profile(path):
    import anndata as ad
    a = ad.read_h5ad(path, backed="r"); genes = [str(g) for g in a.var_names]; a.file.close()
    with h5py.File(path, "r") as h:
        X = h["X"]; idx = X["indices"][:]; dat = X["data"][:].astype(np.float64)
    s = np.bincount(idx, weights=dat, minlength=len(genes))
    return np.log1p(s / s.sum() * 1e6), genes


def infer(path, n_var=3000, top=5):
    M, rg, lines = reference()
    q, qg = context_profile(path)
    qpos = {g: i for i, g in enumerate(qg)}
    common = [i for i, g in enumerate(rg) if g in qpos]
    Mr = M[common]; qv = q[[qpos[rg[i]] for i in common]]
    var_idx = np.argsort(-Mr.var(1))[:n_var]
    A = Mr[var_idx]; b = qv[var_idx]
    ra = np.apply_along_axis(rankdata, 0, A); rb = rankdata(b)
    ra = (ra - ra.mean(0)) / ra.std(0); rb = (rb - rb.mean()) / rb.std()
    rho = (ra * rb[:, None]).mean(0)
    order = np.argsort(-rho)
    lin = lineage_table()
    out = []
    for j in order[:top]:
        row = dict(line=lines[j], spearman=round(float(rho[j]), 4))
        if lin is not None and lines[j] in lin.index:
            row.update({k: str(v) for k, v in lin.loc[lines[j]].iloc[:4].items()} if lin.loc[lines[j]].ndim == 1 else {})
        out.append(row)
    margin = float(rho[order[0]] - rho[order[1]])
    return dict(context=str(path), best=lines[order[0]], rho=round(float(rho[order[0]]), 4), margin=round(margin, 4), top=out)


if __name__ == "__main__":
    res = [infer(p) for p in sys.argv[1:]]
    print(json.dumps(res, indent=1))
