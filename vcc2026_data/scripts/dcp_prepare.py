"""Part 2 preprocessing: harmonise Tahoe, GSE306429 and the gene-perturbation contexts.

RESPONSE SPACE (documented explicitly, as required):
    response = log1p(1e4 * treated_counts / treated_library) - log1p(1e4 * control_counts / control_library)
Both terms are recomputed FROM RAW COUNTS in this script. Nothing is derived from an
already-normalised matrix, so there is no second compression on top of a first one. This is the
same space the frozen gene-perturbation table uses ("source-native log1p(CP10K) treated-minus-
matched untreated") and the same space the VCC count reconstruction inverts, which is why the
residual can be added to the refit prediction at all.

CONTEXT SPACE: pooled UNTREATED counts only (Tahoe DMSO_TF, GSE DMSO wells, gene-side
non-targeting controls), CP10K + log1p on the 5,955-gene axis shared by all three datasets,
using the shared-axis library sum as the denominator so that a dataset's native panel size
cannot leak in as a constant offset.

Controls are always matched within dataset AND within plate. Tahoe: plate|CVCL line. GSE: the
DMSO wells of the same plate (every plate carries DMSO and every plate is a single cell line).
"""
from pathlib import Path
import json
import numpy as np
import h5py

ROOT = Path('/home/wonmin-kim/AIVC/vcc2026_data')
OUT = ROOT / 'data/processed/drug_context_pretrain'
REP = ROOT / 'reports/drug_context_pretrain_2026-09-20'
TAHOE = ROOT / 'reports/simple_residual_2026-09-16/data/tahoe/table.npz'
GSE = ROOT / 'data/external/gse306429/GSE306429_combined_pseudobulk.h5ad'
CTRLS = ROOT / 'reports/cellfm_context_2026-09-20/controls_common'
AXIS = ROOT / 'reports/cellfm_context_2026-09-20/common_axis.json'
BIOS7 = ('HCT116', 'HEK293T', 'HEPG2', 'JURKAT', 'K562', 'KOLF', 'RPE1')
# Contexts that are the SAME biology as a scored gene-perturbation fold. They are dropped from
# drug pretraining on that fold so that the held-out biology contributes untreated controls only.
FOLD_CONFLICT = {'HEPG2': {'tahoe': ['HEPG2_FAMILY'], 'gse': ['HepG2']},
                 'HEK293T': {'tahoe': [], 'gse': ['HEK293T']}}


def cp10k_log1p(counts, denom):
    return np.log1p(1e4 * counts / np.maximum(denom, 1.0)[:, None]).astype(np.float32)


def load_gse():
    f = h5py.File(GSE, 'r')
    o = f['obs']

    def cat(k):
        g = o[k]
        c = np.array([x.decode() if isinstance(x, bytes) else x for x in g['categories'][:]])
        return c[g['codes'][:]]
    return f, dict(cell=cat('cell_id'), cmpd=cat('compound_name'), plate=cat('plate'),
                   rep=cat('replicate'), well=cat('well'), lib=cat('library_id'),
                   dose=o['dose_uM'][:]), [x.decode() for x in f['var']['_index'][:]]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    REP.mkdir(parents=True, exist_ok=True)
    diag = {}

    # ---------------- shared context gene axis ----------------
    g_gene = json.load(open(AXIS))['genes']
    z = np.load(TAHOE, allow_pickle=False)
    g_tahoe = [str(x) for x in z['genes']]
    fg, gm, g_gse = load_gse()
    shared = sorted(set(g_gene) & set(g_tahoe) & set(g_gse))
    json.dump(dict(n=len(shared), genes=shared,
                   note='intersection of the gene-perturbation 6,413 common axis, the Tahoe '
                        '18,150 challenge intersection and the GSE306429 13,290 pseudobulk axis'),
              open(OUT / 'context_axis.json', 'w'))
    sidx = {g: i for i, g in enumerate(shared)}
    print(json.dumps(dict(event='axis', shared=len(shared), gene=len(g_gene),
                          tahoe=len(g_tahoe), gse=len(g_gse))), flush=True)

    # ---------------- TAHOE responses ----------------
    ti = np.array([sidx.get(g, -1) for g in g_tahoe])
    t_take = np.flatnonzero(ti >= 0)
    t_put = ti[t_take]
    cpi = z['control_profile_index']
    tre = cp10k_log1p(z['treated_sum_counts'].astype(np.float64), z['treated_full_library_sums'])
    ctr = cp10k_log1p(z['control_sum_counts'].astype(np.float64), z['control_full_library_sums'])
    t_resp = tre - ctr[cpi]
    del tre
    stored = z['effects']
    diag['tahoe_response_reconstruction'] = dict(
        max_abs_diff_vs_stored=float(np.abs(t_resp - stored).max()),
        rms_diff=float(np.sqrt(np.mean((t_resp - stored) ** 2))),
        note='confirms the recomputed-from-raw-counts response matches the stored effects')
    del stored
    t_bg = np.array([str(x) for x in z['biological_groups']])
    t_drug = np.array([str(x) for x in z['drug_names_by_row']])
    t_ctx = sorted(set(t_bg))
    t_pert = sorted(set(t_drug))
    nt, nc = z['treatment_cell_counts'].astype(np.float64), z['control_cell_counts'].astype(np.float64)
    t_w = 1.0 / (1.0 / np.maximum(nt, 1) + 1.0 / np.maximum(nc[cpi], 1))

    # Tahoe context profiles: pooled DMSO counts over every plate of that biological group
    t_cprof = np.zeros((len(t_ctx), len(shared)), np.float64)
    keys = [str(k) for k in z['control_profile_keys']]
    ccounts = z['control_sum_counts'].astype(np.float64)
    line_of_key = {}
    for r in range(len(t_bg)):
        line_of_key[keys[cpi[r]]] = t_bg[r]
    for ki, k in enumerate(keys):
        b = line_of_key.get(k)
        if b is None:
            continue
        t_cprof[t_ctx.index(b), t_put] += ccounts[ki][t_take]
    np.savez(OUT / 'tahoe.npz',
             resp=t_resp.astype(np.float32), genes=np.array(g_tahoe),
             ctx=np.array([t_ctx.index(b) for b in t_bg], np.int32), ctx_names=np.array(t_ctx),
             pert=np.array([t_pert.index(d) for d in t_drug], np.int32), pert_names=np.array(t_pert),
             plate=np.array([str(x) for x in z['plates']]), is_repeat=z['is_repeat'],
             weight=t_w.astype(np.float32),
             ctx_profile_counts=t_cprof.astype(np.float64))
    print(json.dumps(dict(event='tahoe', rows=len(t_resp), ctx=len(t_ctx), pert=len(t_pert),
                          recon=diag['tahoe_response_reconstruction'])), flush=True)
    del t_resp, ccounts

    # ---------------- GSE responses ----------------
    gi = np.array([sidx.get(g, -1) for g in g_gse])
    g_take = np.flatnonzero(gi >= 0)
    g_put = gi[g_take]
    C = np.asarray(fg['layers/counts'][:], np.float64)
    isd = gm['cmpd'] == 'DMSO'
    plates = gm['plate']
    # matched control = pooled DMSO counts of the SAME plate
    pl_u = sorted(set(plates))
    ctrl_counts = np.zeros((len(pl_u), C.shape[1]))
    for pi, p in enumerate(pl_u):
        ctrl_counts[pi] = C[isd & (plates == p)].sum(0)
    pl_idx = np.array([pl_u.index(p) for p in plates])
    treated = np.flatnonzero(~isd)
    g_resp = (cp10k_log1p(C[treated], C[treated].sum(1))
              - cp10k_log1p(ctrl_counts, ctrl_counts.sum(1))[pl_idx[treated]])
    g_cell = gm['cell'][treated]
    g_pname = np.array([f'{c}@{d:g}uM' for c, d in zip(gm['cmpd'][treated], gm['dose'][treated])])
    g_ctx = sorted(set(g_cell))
    g_pert = sorted(set(g_pname))
    g_cprof = np.zeros((len(g_ctx), len(shared)))
    for ci, c0 in enumerate(g_ctx):
        g_cprof[ci, g_put] = C[isd & (gm['cell'] == c0)].sum(0)[g_take]
    np.savez(OUT / 'gse.npz',
             resp=g_resp.astype(np.float32), genes=np.array(g_gse),
             ctx=np.array([g_ctx.index(c) for c in g_cell], np.int32), ctx_names=np.array(g_ctx),
             pert=np.array([g_pert.index(p) for p in g_pname], np.int32), pert_names=np.array(g_pert),
             plate=plates[treated], is_repeat=(gm['rep'][treated] == '2'),
             weight=np.sqrt(C[treated].sum(1)).astype(np.float32) / 1e3,
             compound=gm['cmpd'][treated], dose=gm['dose'][treated],
             ctx_profile_counts=g_cprof)
    print(json.dumps(dict(event='gse', rows=len(g_resp), ctx=len(g_ctx), pert=len(g_pert))), flush=True)

    # ---------------- gene-side contexts (untreated controls only) ----------------
    ge_cprof = np.zeros((len(BIOS7), len(shared)))
    gene_diag = []
    for bi, b in enumerate(BIOS7):
        with h5py.File(CTRLS / f'{b}.h5ad', 'r') as fh:
            gv = [x.decode() if isinstance(x, bytes) else str(x) for x in fh['var']['_index'][:]]
            X = fh['X']
            Xd = np.asarray(X[:], np.float64) if isinstance(X, h5py.Dataset) else None
            if Xd is None:
                import scipy.sparse as sp
                Xd = sp.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]),
                                   shape=tuple(X.attrs['shape'])).toarray().astype(np.float64)
        col = np.array([sidx.get(g, -1) for g in gv])
        tk = np.flatnonzero(col >= 0)
        ge_cprof[bi, col[tk]] = Xd[:, tk].sum(0)
        gene_diag.append(dict(bio=b, cells=int(Xd.shape[0]), native_genes=len(gv),
                              mapped=int(len(tk)), integral=bool(np.all(Xd == np.floor(Xd)))))
    np.savez(OUT / 'gene_contexts.npz', ctx_names=np.array(BIOS7), ctx_profile_counts=ge_cprof)
    diag['gene_control_sources'] = gene_diag

    # ---------------- dataset-separation diagnostic (Part 2.8) ----------------
    prof = np.concatenate([t_cprof, g_cprof, ge_cprof])
    lab = np.array(['tahoe'] * len(t_cprof) + ['gse'] * len(g_cprof) + ['gene'] * len(ge_cprof))
    raw = cp10k_log1p(prof, prof.sum(1))
    diag['dataset_separation'] = separation(raw, lab)
    # control-only per-dataset standardisation: centre each dataset's controls on its own mean
    adj = raw.copy()
    for d in np.unique(lab):
        m = lab == d
        adj[m] -= adj[m].mean(0)
    diag['dataset_separation_after_centering'] = separation(adj, lab)
    # names are DATASET-QUALIFIED. Without the prefix 'A549' collides between Tahoe and GSE and
    # 'HEK293T' between GSE and the gene-perturbation side, and a plain name->row dict then hands
    # one dataset's drug rows another dataset's control profile.
    names = np.array([f'tahoe:{c}' for c in t_ctx] + [f'gse:{c}' for c in g_ctx]
                     + [f'gene:{b}' for b in BIOS7])
    assert len(set(names.tolist())) == len(names), 'context keys must be unique'
    np.savez(OUT / 'contexts.npz', profile_raw=raw, profile_adj=adj, dataset=lab,
             names=names, genes=np.array(shared))
    json.dump(diag, open(REP / 'preprocessing_diagnostics.json', 'w'), indent=2, default=float)
    print(json.dumps(dict(event='done', separation=diag['dataset_separation'],
                          after=diag['dataset_separation_after_centering']), indent=2), flush=True)


def separation(Xp, lab):
    """How separable are the dataset groups? 1.0 = a single PC splits them perfectly."""
    Xc = Xp - Xp.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc = U[:, :2] * S[:2]
    out = dict(pc1_var_frac=float(S[0] ** 2 / (S ** 2).sum()))
    # nearest-neighbour dataset purity, a direct "do they cluster by assay" answer
    Dm = ((Xc[:, None, :] - Xc[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(Dm, np.inf)
    nn = Dm.argmin(1)
    out['nn_same_dataset_fraction'] = float((lab[nn] == lab).mean())
    for d in np.unique(lab):
        m = lab == d
        out[f'pc1_mean_{d}'] = float(pc[m, 0].mean())
    # between/within variance ratio on PC1
    mus = [pc[lab == d, 0].mean() for d in np.unique(lab)]
    out['pc1_between_over_within'] = float(np.var(mus) / max(np.mean(
        [np.var(pc[lab == d, 0]) for d in np.unique(lab)]), 1e-9))
    return out


if __name__ == '__main__':
    main()
