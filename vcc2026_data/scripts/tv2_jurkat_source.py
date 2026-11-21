"""GSE249595 (genome-scale CRISPRi Perturb-seq in Jurkat, unstimulated + stimulated) -> transfer_v2 sources.

Per channel: STARsolo matrices over all barcodes (features x barcodes). Cells = barcodes with transcriptome UMI
>= 10% of the 99th-percentile UMI among the top expected cells (the authors' EmptyDrops-like rule, here with a
floor of 1,000 UMI). Guide assignment: exactly one target gene among the guides a cell carries (binary sgRNA
matrix). Condition: cell-hashing counts, 'untreated' vs sum of 'activated*', dominant label with >= 2x margin.
Writes JURKAT_UNSTIM / JURKAT_STIM in the tv2_sources format (fractions on the official axis for the needed
targets + genome-wide mean LFC for centering).
"""
import gzip, io, json, pathlib, sys
import numpy as np, pandas as pd
from scipy import sparse
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_sources as TS

RAW = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/raw/jurkat_gse249595")
OUT = TS.OUT
NG = TS.NG


def read_mtx(path):
    with gzip.open(path, "rt") as f:
        header = []
        while True:
            line = f.readline()
            if not line.startswith("%"):
                shape = tuple(int(x) for x in line.split()[:2]); break
        df = pd.read_csv(f, sep=" ", header=None, names=["r", "c", "v"], dtype={"r": np.int32, "c": np.int32, "v": np.float32})
    return sparse.csc_matrix((df.v.values, (df.r.values - 1, df.c.values - 1)), shape=shape)   # features x barcodes


def official_map(features_path):
    feat = pd.read_csv(features_path, sep="\t", header=None, names=["id", "sym", "type"])
    off = pd.read_csv("/home/wonmin-kim/AIVC/vcc2026_data/metadata/vcc_official/official_features.csv")
    ens = {e: i for i, e in zip(off.feature_index, off.ensembl_id) if isinstance(e, str)}
    pos = np.full(len(feat), -1)
    for k, (eid, sym) in enumerate(zip(feat.id, feat.sym)):
        i = ens.get(str(eid).split(".")[0])
        if i is None:
            i = TS.GIC.get(TS.canon(str(sym)), -1)
        pos[k] = i
    return pos


def main():
    """High-MOI design (~12 real guides per cell at >=3 UMI; singletons are ambient): marginal estimator.
    For target t: cells carrying any guide of t vs ALL cells of the same condition (other guides are random,
    so their effects average out; the shared 'many knockdowns' component is removed later by centering)."""
    need = TS.needed_targets() if len(sys.argv) < 2 else {TS.canon(t) for t in pathlib.Path(sys.argv[1]).read_text().split()}
    rng = np.random.default_rng(0)
    chans = sorted({p.name.split("_")[0] + "_" + p.name.split("_")[1] for p in RAW.glob("GSM*_transcriptome_matrix.mtx.gz")})
    gmap = None; log = []
    acc = {c: None for c in ("UNSTIM", "STIM")}
    for ch in chans:
        T = read_mtx(RAW / f"{ch}_transcriptome_matrix.mtx.gz")
        if gmap is None:
            gmap = official_map(RAW / f"{ch}_transcriptome_features.tsv.gz"); ok = gmap >= 0
            gfeat = pd.read_csv(RAW / f"{ch}_guides_features.tsv.gz", sep="\t", header=None)[0].astype(str).values
            gene_of = np.array([TS.canon(g.rsplit("_", 1)[0]) for g in gfeat])
            all_genes = sorted(set(gene_of) - {TS.canon("non-targeting")})
            extra = set(rng.choice([g for g in all_genes if g not in need], size=2000, replace=False))
            tlist = sorted((need & set(all_genes)) | extra)
            tpos = {g: i for i, g in enumerate(tlist)}
            gi = np.array([tpos.get(g, -1) for g in gene_of])
            rows_g = np.flatnonzero(gi >= 0)
            Gmap = sparse.csr_matrix((np.ones(rows_g.size), (rows_g, gi[rows_g])), shape=(len(gfeat), len(tlist)))
            for c in acc:
                acc[c] = dict(S=np.zeros((len(tlist), NG)), T=np.zeros(len(tlist)), n=np.zeros(len(tlist)), Sc=np.zeros(NG), Tc=0.0, nc=0)
        umi = np.asarray(T.sum(0)).ravel()
        thr = max(1000.0, 0.1 * np.percentile(np.sort(umi)[::-1][:20000], 99))
        cells = np.flatnonzero(umi >= thr)
        G = read_mtx(RAW / f"{ch}_guides_matrix.mtx.gz")[:, cells]
        Gb = (G >= 3).astype(np.float64).T.tocsr()                  # cells x guides
        A = (Gb @ Gmap); A.data[:] = 1.0                            # cells x targets (carries any guide of t)
        H = read_mtx(RAW / f"{ch}_labels_matrix.mtx.gz")[:, cells].toarray()
        hfeat = pd.read_csv(RAW / f"{ch}_labels_features.tsv.gz", sep="\t", header=None)[0].astype(str).values
        act = H[[i for i, h in enumerate(hfeat) if h.startswith("activated")]].sum(0); un = H[[i for i, h in enumerate(hfeat) if h.startswith("untreated")]].sum(0)
        cond = np.where(un >= 2 * np.maximum(act, 1), "UNSTIM", np.where(act >= 2 * np.maximum(un, 1), "STIM", ""))
        Tc = T[:, cells].T.tocsr()                                   # cells x features
        for c in acc:
            m = np.flatnonzero(cond == c)
            if m.size == 0:
                continue
            X = Tc[m]; Am = A[m]
            Sf = (Am.T @ X).toarray()                                # targets x features
            v = np.zeros((Sf.shape[0], NG)); v[:, gmap[ok]] += Sf[:, ok]
            acc[c]["S"] += v; acc[c]["T"] += np.asarray(Am.T @ umi[cells][m]).ravel(); acc[c]["n"] += np.asarray(Am.sum(0)).ravel()
            s = np.asarray(X.sum(0)).ravel(); w = np.zeros(NG); w[gmap[ok]] += s[ok]
            acc[c]["Sc"] += w; acc[c]["Tc"] += float(umi[cells][m].sum()); acc[c]["nc"] += int(m.size)
        log.append(dict(channel=ch, cells=int(len(cells)), thr=float(thr), guides_per_cell_median=float(np.median(np.diff(Gb.indptr))),
                        unstim=int((cond == "UNSTIM").sum()), stim=int((cond == "STIM").sum())))
        print(json.dumps(log[-1]), flush=True)
    measured = np.zeros(NG, bool); measured[gmap[gmap >= 0]] = True
    for c, a in acc.items():
        out = TS.Acc(need)
        fc = np.where(measured, a["Sc"] / a["Tc"], np.nan)
        for k, t in enumerate(tlist):
            if a["n"][k] >= 5 and a["T"][k] > 0:
                out.add(t, np.where(measured, a["S"][k] / a["T"][k], np.nan), fc, int(a["n"][k]))
        out.save(f"JURKAT_{c}", dict(source="GSE249595", condition=c, design="high-MOI marginal (>=3 UMI guide calls) vs all cells",
                                    cells=int(a["nc"]), targets=len(tlist), channels=log))


if __name__ == "__main__":
    main()
