#!/usr/bin/env python
"""One sequential pass over replogle_concat.h5ad -> per-cell-line .npy caches.

The 30 GB file stores cells interleaved across cell lines, so scattered fancy
indexing thrashes. We stream it once in contiguous blocks and write one
float32 matrix per cell line, capping cells per perturbation.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/build_cache.py
"""
import os, sys, time
import numpy as np, h5py

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONCAT = os.path.join(ROOT, "data/replogle/replogle_concat.h5ad")
CACHE = os.path.join(ROOT, "work/cellline_cache"); os.makedirs(CACHE, exist_ok=True)
CONTROL = "non-targeting"
MAX_PER_PERT, MAX_CONTROL, SEED = 128, 6000, 0


def main():
    h = h5py.File(CONCAT, "r")
    def cat(c):
        g = h["obs"][c]
        return np.array([x.decode() for x in g["categories"][:]])[g["codes"][:]]
    cl, gene = cat("cell_line"), cat("gene")
    vg = h["var"]["gene_name_index"]
    genes = ([x.decode() for x in vg[:]] if not isinstance(vg, h5py.Group)
             else np.array([x.decode() for x in vg["categories"][:]])[vg["codes"][:]])
    np.save(os.path.join(CACHE, "genes.npy"), np.asarray(genes, dtype=object))

    rng = np.random.default_rng(SEED)
    keep = np.zeros(len(cl), bool)
    for c in np.unique(cl):
        m = cl == c
        for p in np.unique(gene[m]):
            idx = np.where(m & (gene == p))[0]
            cap = MAX_CONTROL if p == CONTROL else MAX_PER_PERT
            if len(idx) > cap:
                idx = rng.choice(idx, cap, replace=False)
            keep[idx] = True
    print(f"selected {keep.sum()} / {len(keep)} cells")

    rows = np.where(keep)[0]
    out = {c: np.empty((int((cl[rows] == c).sum()), len(genes)), np.float32)
           for c in np.unique(cl)}
    fill = {c: 0 for c in out}
    X = h["X"]; N = X.shape[0]; BLK = 25000
    t0 = time.time()
    for s in range(0, N, BLK):
        e = min(s + BLK, N)
        sel = rows[(rows >= s) & (rows < e)]
        if not len(sel):
            continue
        blk = X[s:e]
        sub = blk[sel - s]
        for c in np.unique(cl[sel]):
            m = cl[sel] == c
            n = int(m.sum())
            out[c][fill[c]:fill[c] + n] = sub[m]
            fill[c] += n
        print(f"  {e}/{N} ({time.time()-t0:.0f}s)", flush=True)
    h.close()

    for c, M in out.items():
        np.save(os.path.join(CACHE, f"X_{c}.npy"), M)
        lab = gene[rows][cl[rows] == c]
        np.save(os.path.join(CACHE, f"pert_{c}.npy"), np.asarray(lab, dtype=object))
        print(f"{c}: {M.shape} perts={len(np.unique(lab))-1}")


if __name__ == "__main__":
    main()
