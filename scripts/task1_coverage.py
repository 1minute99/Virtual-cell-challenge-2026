#!/usr/bin/env python
"""Task 1 - Coverage counts for training-free VCC 2026 baselines (GATES Tasks 2-4).

Counts VCC 2026 validation target genes available as *perturbations* in Replogle
2022 (K562 genome-scale, K562 essential, RPE1 essential) and in the STATE ST
perturbation vocabulary; plus gene-space overlaps and normalization notes.

Uses collision-safe HGNC harmonization (see scripts/genemap.py) and joins on
Ensembl gene id where both sides carry it.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task1_coverage.py
"""
import json, os, pickle, sys
import numpy as np, pandas as pd, anndata as ad, h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genemap import GeneMap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VCC  = os.path.join(ROOT, "data/2026_AIVC _validation_controls")
GWPS = os.path.join(ROOT, "data/replogle_gwps")
OUT  = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)

GM = GeneMap()

# ------------------------------------------------------------------ VCC side
targets   = [l.strip() for l in open(os.path.join(VCC, "pert_counts.csv"))][1:]
targets   = [t for t in targets if t]
vcc_genes = [l.strip() for l in open(os.path.join(VCC, "gene_names.csv"))][1:]
vcc_genes = [g for g in vcc_genes if g]
manifest  = json.load(open(os.path.join(VCC, "manifest.json")))

tgt_c  = {t: GM.canon(t) for t in targets}
tgt_e  = {t: GM.ens(t)   for t in targets}
vcc_c  = set(GM.canon(g) for g in vcc_genes)
vcc_e  = set(e for e in (GM.ens(g) for g in vcc_genes) if e)

with h5py.File(os.path.join(VCC, "context_A.h5ad"), "r") as h:
    d = h["X"]["data"][:200000]
vcc_isint = bool(np.all(np.abs(d - np.round(d)) < 1e-6))

# ------------------------------------------------------------------ Replogle
SOURCES = {
    "Replogle_K562_genomescale": "K562_gwps_raw_bulk_01.h5ad",
    "Replogle_K562_essential":   "K562_essential_raw_bulk_01.h5ad",
    "Replogle_RPE1_essential":   "rpe1_raw_bulk_01.h5ad",
}
rows, lists = [], {}
for name, fn in SOURCES.items():
    A = ad.read_h5ad(os.path.join(GWPS, fn))
    p_sym = set(GM.canon(x.split("_")[1])  for x in A.obs_names)
    p_ens = set(x.split("_")[-1]           for x in A.obs_names)
    g_sym = set(GM.canon(g) for g in A.var["gene_name"].astype(str))
    g_ens = set(A.var_names.astype(str))
    hit = sorted(t for t in targets
                 if tgt_c[t] in p_sym or (tgt_e[t] and tgt_e[t] in p_ens))
    rows.append(dict(source=name, n_perturbations=len(p_ens), n_readout_genes=len(g_ens),
                     vcc_targets_perturbed=len(hit),
                     vcc_targets_in_readout=len(set(tgt_c.values()) & g_sym),
                     vcc_genespace_overlap_symbol=len(vcc_c & g_sym),
                     vcc_genespace_overlap_ensembl=len(vcc_e & g_ens),
                     vcc_genes_dropped=len(vcc_genes) - len(vcc_c & g_sym)))
    lists[name] = hit
    del A

# ------------------------------------------------------------------ STATE ST
st_dir = os.path.join(ROOT, "checkpoints/state/ST-SE-Replogle/zeroshot/k562")
vd = pickle.load(open(os.path.join(st_dir, "var_dims.pkl"), "rb"))
st_vocab = set(GM.canon(str(x)) for x in vd["pert_names"])
st_out   = set(GM.canon(str(x)) for x in vd["gene_names"])
st_hit   = sorted(t for t in targets if tgt_c[t] in st_vocab)
rows.append(dict(source="STATE_ST (ST-SE-Replogle, pert_rep=onehot)",
                 n_perturbations=len(vd["pert_names"]), n_readout_genes=len(st_out),
                 vcc_targets_perturbed=len(st_hit),
                 vcc_targets_in_readout=len(set(tgt_c.values()) & st_out),
                 vcc_genespace_overlap_symbol=len(vcc_c & st_out),
                 vcc_genespace_overlap_ensembl=-1,
                 vcc_genes_dropped=len(vcc_genes) - len(vcc_c & st_out)))
lists["STATE_ST"] = st_hit

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "task1_coverage.csv"), index=False)
json.dump({"vcc_targets": targets, "intersections": lists,
           "genemap_stats": GM.stats()},
          open(os.path.join(OUT, "task1_gene_lists.json"), "w"), indent=1)

norm = pd.DataFrame([
 dict(dataset="VCC 2026 validation (10x Flex)", n_genes=len(vcc_genes),
      unit="raw UMI counts" if vcc_isint else "non-integer",
      detail=f"csr, integer-valued={vcc_isint}; 18,400 non-targeting cells x 3 contexts; var id = HGNC symbol"),
 dict(dataset="Replogle K562 gwps raw_bulk", n_genes=8248, unit="raw pseudobulk mean UMI",
      detail="11,258 perturbation pseudobulks; var id = ENSG (+gene_name); obs id = idx_SYMBOL_transcript_ENSG"),
 dict(dataset="Replogle *_normalized_bulk", n_genes=8248, unit="z-score vs non-targeting",
      detail="Replogle's own per-gene z-scored pseudobulk relative to control"),
 dict(dataset="Replogle concat single-cell", n_genes=6546, unit="log-normalized float (max ~5.9)",
      detail="643,413 cells x 6,546 genes; 4 cell lines; ESSENTIAL screen only"),
 dict(dataset="STATE ST output space", n_genes=len(st_out), unit="HVG subset",
      detail=f"pert_rep=onehot, pert_dim={vd['pert_dim']}, gene_dim={vd['gene_dim']}, hvg_dim={vd['hvg_dim']}"),
])
norm.to_csv(os.path.join(OUT, "task1_normalization.csv"), index=False)

print(df.to_string(index=False)); print()
print(norm.to_string(index=False)); print()
print(f"VCC: {len(targets)} targets, {len(vcc_genes)} genes, contexts={manifest['contexts']}")
print(f"GeneMap: {GM.stats()}")
print(f"All targets are Approved HGNC symbols: {all(t==c for t,c in tgt_c.items())}")
print(f"Targets with Ensembl id: {sum(1 for v in tgt_e.values() if v)}/{len(targets)}")
for k, v in lists.items():
    print(f"  {k}: {len(v)}/300 targets perturbed")
print(f"\nGATE Task2 (Replogle): K562-genomescale={len(lists['Replogle_K562_genomescale'])} -> PROCEED; "
      f"RPE1={len(lists['Replogle_RPE1_essential'])} -> NO DATA")
print(f"GATE Task3 (STATE ST): {len(lists['STATE_ST'])} -> SKIP Task 3")
