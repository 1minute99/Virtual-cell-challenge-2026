"""Task 9 -- cell-level split-half extraction for the genetic contexts (input stage of rq2_gene_params.py).

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 9, definitions 3-5, Amendments A2/A3.
For every source (one per Count-1 line, plus VCC2025-H1 for sensitivities) this writes
output/data-analysis/results/rq2_gene_cache/<SRC>.npz with, per GROUP (a target, or the control pool 'CTRL';
Feng: '<line>|<target>' and '<line>|CTRL'):
  n[g]           cells (CD4: cells summed over pseudobulk units)
  nA[s, g]       cells in half A for split seed s (s = 0..4)
  full[g, G]     summed raw counts over the gene axis (float64)
  A[s, g, G]     summed raw counts of half A (float32; half B = full - A)
  lib[g], libA[s, g]  summed native libraries (full / half A)
  ctrl_sample    up to 2,000 random control cells (log1p CP10K on the axis; for the control-PC sensitivity)

ASSUMPTIONS (stated; see knowledge/synthesis/rq2_gene_params-assumptions.md)
  E1 Gene axis = official 18,533-gene axis indices (common.to_official, the harness's canon mapping) measured in
     EVERY source's count matrix (intersection; duplicate source genes mapping to one official gene are summed).
  E2 Library = the cell's native total UMI over its full expression matrix (obs total where the matrix is
     gene-filtered: RN 'UMI_count', Orion 'total_counts', KOLF 'total_counts', DLD-1 'total', CD4 unit 'total_counts';
     otherwise the row sum over all features).
  E3 Target universe U = targets with >= 20 cells in at least one of the three Replogle/Nadig ESSENTIAL screens
     (RPE1, HepG2, Jurkat; K562 is read from GWPS = the registry's K562 family). Only U targets + controls are read.
  E4 Split halves: within each group the cells (in the source's fixed reading order) are permuted with
     default_rng([SPLIT_BASE, seed, src_id, group]) and the first floor(n/2) go to half A (= rq2_params.half_assignment
     with a per-source base). Controls are split the same way, so A and B have SEPARATE control halves.
     CD4 has no cell-level data on disk (only per run x donor x guide x condition pseudobulk units): its halves are
     unit-level, balanced greedily by cell count within each target (DEVIATION, logged).
  E5 Raw counts only (KOLF: layers/counts; others X). Cells: RN all; Orion all called cells with a gene_target
     (as the registry store 'all called cells'); KOLF all QC'd cells; DLD-1 singlets, control label 'NTC';
     H1GW the builder's guide call (>= 5 UMIs, runner-up <= 20 %); VCC2025-H1 the Training split; Feng all cells of
     the targeted screen with a single-target call (pb_robust_feng rules).
Deterministic (seeds fixed). Run: python rq2_gene_extract.py plan | extract SRC
"""
from __future__ import annotations
import csv, glob, gzip, json, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py

ROOT = Path('/home/wonmin-kim/AIVC'); V = ROOT / 'vcc2026_data'; RAW = V / 'raw'
sys.path.insert(0, str(ROOT / 'output/data-analysis'))
import common as CM
CACHE = ROOT / 'output/data-analysis/results/rq2_gene_cache'
SPLIT_BASE = 20261108
N_SEEDS = 5
MIN_CELLS = 20
N_CTRL_SAMPLE = 2000
SOURCES = ['K562', 'RPE1', 'HEPG2', 'JURKAT', 'HCT116', 'HEK293T', 'CD4', 'KOLF', 'FENG', 'DLD1', 'H1GW', 'VCC_H1']
SRC_ID = {s: i for i, s in enumerate(SOURCES)}
RN_FILE = {'K562': 'ReplogleWeissman2022_K562_gwps', 'RPE1': 'ReplogleWeissman2022_rpe1',
           'HEPG2': 'NadigOConner2024_hepg2', 'JURKAT': 'NadigOConner2024_jurkat'}
ESSENTIAL = ['RPE1', 'HEPG2', 'JURKAT']


def log(**k):
    print(json.dumps(k, default=str), flush=True)


def canon_target(t):
    canon, _ = CM._canon()
    return canon(str(t))


def _cat(h, key):
    o = h[key]
    if isinstance(o, h5py.Group):
        c = o['categories'][:].astype(str)
        return c[o['codes'][:]]
    return o[:].astype(str)


# ============================================================ metadata per source
# meta(src) -> dict(label=array[str] per stream cell ('CTRL' = control, '' = drop), genes=array[str] source columns,
#                   lib=array per stream cell or None (row sums computed while streaming))

def meta_rn(src):
    f = RAW / f'scperturb/{RN_FILE[src]}.h5ad'
    with h5py.File(f, 'r') as h:
        g = _cat(h, 'obs/gene'); lib = h['obs/UMI_count'][:].astype(np.float64)
        genes = h['var/gene_name'][:].astype(str) if not isinstance(h['var/gene_name'], h5py.Group) else _cat(h, 'var/gene_name')
    lab = np.where(g == 'non-targeting', 'CTRL', g)
    return dict(label=lab, genes=genes, lib=lib)


def meta_orion(src):
    import pyarrow.parquet as pq
    files = sorted(glob.glob(str(RAW / f'orion/data/{src}_Batch*.parquet')))
    labs, libs = [], []
    for fn in files:
        t = pq.read_table(fn, columns=['gene_target', 'total_counts']).to_pandas()
        g = t.gene_target.astype(str).to_numpy()
        labs.append(np.where(g == 'Non-Targeting', 'CTRL', np.where(t.gene_target.isna(), '', g)))
        libs.append(t.total_counts.to_numpy(np.float64))
    gm = pd.read_parquet(RAW / 'orion/metadata/gene_metadata.parquet').sort_values('gene_token_id')
    assert (gm.gene_token_id.to_numpy() == np.arange(len(gm))).all()
    return dict(label=np.concatenate(labs), genes=gm.gene_name.astype(str).to_numpy(), lib=np.concatenate(libs),
                files=files, file_rows=[len(x) for x in labs])


def meta_kolf(src):
    f = RAW / 'kolf/KOLF_Pan_Genome_QC_Filtered.h5ad'
    with h5py.File(f, 'r') as h:
        g = _cat(h, 'obs/gene_target'); lib = h['obs/total_counts'][:].astype(np.float64)
        genes = _cat(h, 'var/_index')
    return dict(label=np.where(g == 'NTC', 'CTRL', g), genes=genes, lib=lib)


def meta_dld1(src):
    D = RAW / 'dld1_gse337988'
    tg, dem, tot = [], [], []
    with gzip.open(D / 'cells_meta.csv.gz', 'rt') as f:
        r = csv.reader(f); next(r)
        for row in r:
            if len(row) < 4:
                tg.append(''); dem.append(''); tot.append(0); continue
            tg.append(','.join(row[1:-2])); dem.append(row[-2]); tot.append(float(row[-1]))
    tg = np.array(tg); sing = np.array(dem) == 'singlet'
    lab = np.where(~sing, '', np.where(tg == 'NTC', 'CTRL', tg))
    genes = pd.read_csv(D / 'genes.csv')['symbol'].astype(str).to_numpy()
    with h5py.File(D / 'GSE337988_sublib2_processed_objects_Low_assays.h5', 'r') as h:
        assert h['assay001'].shape == (len(lab), len(genes))
    return dict(label=lab, genes=genes, lib=np.array(tot, np.float64))


H1GW_FILES = lambda: sorted(glob.glob(str(RAW / 'h1_gw_gse295214/*_filtered_feature_bc_matrix.h5')))


def _h1gw_lane(fn):
    from scipy import sparse
    with h5py.File(fn, 'r') as f:
        m = f['matrix']; shape = tuple(int(v) for v in m['shape'][:])
        X = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=shape)
        ft = m['features/feature_type'][:].astype(str); nm = m['features/name'][:].astype(str)
    gex = np.flatnonzero(ft == 'Gene Expression'); gd = np.flatnonzero(ft == 'CRISPR Guide Capture')
    lab = np.array(['NT' if x.startswith('non-targeting') else x.split('_')[0] for x in nm[gd]])
    ul, inv = np.unique(lab, return_inverse=True)
    G = X[gd, :].tocsr()
    A = (sparse.csr_matrix((np.ones(len(gd)), (inv, np.arange(len(gd)))), shape=(len(ul), len(gd))) @ G).toarray()
    top = A.argmax(0); u1 = A[top, np.arange(A.shape[1])]
    A[top, np.arange(A.shape[1])] = -1; u2 = A.max(0)
    ok = (u1 >= 5) & (u2 <= 0.2 * u1)                          # = ml_h1gw_build.py MIN_UMI / DOM
    call = np.where(ok, ul[top], '')
    call = np.where(call == 'NT', 'CTRL', call)
    E = X[gex, :].T.tocsr()
    return call, E, nm[gex]


def meta_h1gw(src):
    labs, genes = [], None
    for fn in H1GW_FILES():
        call, E, g = _h1gw_lane(fn)
        labs.append(call)
        if genes is None:
            genes = g
        assert (g == genes).all()
    return dict(label=np.concatenate(labs), genes=genes, lib=None, file_rows=[len(x) for x in labs])


def meta_vcc(src):
    with h5py.File(RAW / 'vcc2025/adata_Training.h5ad', 'r') as h:
        g = _cat(h, 'obs/target_gene'); genes = _cat(h, 'var/_index')
    return dict(label=np.where(g == 'non-targeting', 'CTRL', g), genes=genes, lib=None)


def meta_cd4(src):
    with h5py.File(RAW / 'cd4/GWCD4i.pseudobulk_merged.h5ad', 'r') as h:
        gt = _cat(h, 'obs/guide_type'); tg = _cat(h, 'obs/perturbed_gene_name')
        n = h['obs/n_cells'][:].astype(np.int64); lib = h['obs/total_counts'][:].astype(np.float64)
        guide = _cat(h, 'obs/guide_id'); genes = _cat(h, 'var/gene_name')
    lab = np.where(gt == 'non-targeting', 'CTRL', tg)
    return dict(label=lab, genes=genes, lib=lib, unit_cells=n, guide=guide)


def meta_feng(src):
    sys.path.insert(0, str(V / 'scripts'))
    import ingest_feng_targeted as IFT
    m = pd.read_csv(RAW / 'feng/TargetedScreen_Cell-Metadata.tsv.gz', sep='\t')
    m['bc'] = m.Cell_ID.map(IFT.norm_bc)
    call = m.Guide_Call.astype(str); tgt = call.map(IFT.IF.guide_to_gene)
    ntc = tgt.map(lambda s: bool(IFT.IF.NTC_EXACT.match(str(s)))); un = call.str.lower().isin(['unassigned', 'nan'])
    m['target'] = np.where(ntc, 'CTRL', np.where(un, '', tgt.astype(str)))
    with gzip.open(RAW / 'feng/TargetedScreen_RNA-UMI-Counts.csv.gz', 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
    cols = [IFT.norm_bc(c) for c in header[1:]]
    by = dict(zip(m.bc, m.index)); rows = np.array([by.get(c, -1) for c in cols])
    line = np.where(rows >= 0, m.Cell_Line.astype(str).to_numpy()[np.maximum(rows, 0)], '')
    t = np.where(rows >= 0, m.target.to_numpy()[np.maximum(rows, 0)], '')
    lab = np.where((rows >= 0) & (t != ''), np.char.add(np.char.add(line.astype(str), '|'), t.astype(str)), '')
    return dict(label=lab, genes=None, lib=None, line=line)


META = {'K562': meta_rn, 'RPE1': meta_rn, 'HEPG2': meta_rn, 'JURKAT': meta_rn, 'HCT116': meta_orion,
        'HEK293T': meta_orion, 'KOLF': meta_kolf, 'DLD1': meta_dld1, 'H1GW': meta_h1gw, 'VCC_H1': meta_vcc,
        'CD4': meta_cd4, 'FENG': meta_feng}


def feng_genes():
    """Gene ids of the Feng RNA matrix rows (first column), cached in the cache dir."""
    p = CACHE / 'feng_genes.json'
    if p.exists():
        return np.array(json.loads(p.read_text()))
    out = []
    proc = subprocess.Popen(['pigz', '-dc', str(RAW / 'feng/TargetedScreen_RNA-UMI-Counts.csv.gz')], stdout=subprocess.PIPE)
    proc.stdout.readline()
    for line in proc.stdout:
        out.append(line[:line.index(b',')].decode())
    proc.wait()
    sym = [x.split(':')[1] if x.count(':') >= 2 else x for x in out]
    p.write_text(json.dumps(sym))
    return np.array(sym)


# ============================================================ plan: universe, gene axis, per-source counts

def plan():
    CACHE.mkdir(parents=True, exist_ok=True)
    info, axes = {}, {}
    for src in SOURCES:
        t0 = time.time()
        m = META[src](src)
        genes = feng_genes() if src == 'FENG' else m['genes']
        om = CM.to_official(genes)
        axes[src] = set(om[om >= 0].tolist())
        lab = m['label']
        if src == 'CD4':
            keys, inv = np.unique(lab, return_inverse=True); cnt = np.bincount(inv, weights=m['unit_cells'])
        else:
            keys, cnt = np.unique(lab, return_counts=True)
        cnt = dict(zip(keys.tolist(), cnt.astype(int).tolist())); cnt.pop('', None)
        if src == 'FENG':                                    # pooled over lines
            pool = {}
            for k, v in cnt.items():
                t = k.split('|', 1)[1]; pool[t] = pool.get(t, 0) + v
            info[src] = dict(per_line_group=cnt, counts=pool)
        else:
            info[src] = dict(counts=cnt)
        log(event='meta', src=src, cells=int(len(lab)), groups=len(cnt), genes_on_axis=len(axes[src]),
            secs=round(time.time() - t0, 1))
    # target universe (E3), canonical names
    cn = {s: {} for s in SOURCES}
    for s in SOURCES:
        for k, v in info[s]['counts'].items():
            if k == 'CTRL':
                continue
            c = canon_target(k); cn[s][c] = cn[s].get(c, 0) + v
    U = sorted(set().union(*[{t for t, v in cn[s].items() if v >= MIN_CELLS} for s in ESSENTIAL]))
    ax = sorted(set.intersection(*axes.values()))
    out = dict(universe=U, gene_axis=ax, n_universe=len(U), n_axis=len(ax),
               axis_per_source={s: len(a) for s, a in axes.items()},
               ctrl_cells={s: int(info[s]['counts'].get('CTRL', 0)) for s in SOURCES},
               targets_ge20_in_U={s: int(sum(1 for t in U if cn[s].get(t, 0) >= MIN_CELLS)) for s in SOURCES},
               feng_ctrl_per_line={k.split('|')[0]: v for k, v in info['FENG']['per_line_group'].items() if k.endswith('|CTRL')})
    (CACHE / 'plan.json').write_text(json.dumps(out))
    log(event='plan', n_universe=len(U), n_axis=len(ax), targets_ge20_in_U=out['targets_ge20_in_U'], ctrl=out['ctrl_cells'])


# ============================================================ accumulator

class Acc:
    """Group sums for full depth and half A per split seed. Cells are added in stream order with stream index."""

    def __init__(self, src, group_of_cell, group_names, n_axis, unit_cells=None, ctrl_group=None, ctrl_groups=None):
        self.src = src; self.gnames = group_names; self.G = n_axis
        self.code = group_of_cell                                   # per stream cell: group id or -1
        ng = len(group_names); self.ng = ng
        w = np.ones(len(group_of_cell)) if unit_cells is None else unit_cells.astype(np.float64)
        self.w = w
        keep = group_of_cell >= 0
        self.n = np.bincount(group_of_cell[keep], weights=w[keep], minlength=ng)
        # half assignment
        self.inA = np.zeros((N_SEEDS, len(group_of_cell)), bool)
        idx_by_g = [np.flatnonzero(group_of_cell == g) for g in range(ng)]
        for s in range(N_SEEDS):
            for g, ii in enumerate(idx_by_g):
                rng = np.random.default_rng([SPLIT_BASE, s, SRC_ID[src], g])
                if unit_cells is None:
                    perm = rng.permutation(len(ii))
                    self.inA[s, ii] = perm < len(ii) // 2
                else:                                               # CD4 units: greedy balance by cells
                    order = rng.permutation(len(ii)); ca = cb = 0.0
                    for j in order:
                        if ca <= cb:
                            self.inA[s, ii[j]] = True; ca += w[ii[j]]
                        else:
                            cb += w[ii[j]]
        self.nA = np.stack([np.bincount(group_of_cell[keep], weights=(w * self.inA[s])[keep], minlength=ng)
                            for s in range(N_SEEDS)])
        self.full = np.zeros((ng, n_axis)); self.A = np.zeros((N_SEEDS, ng, n_axis))
        self.lib = np.zeros(ng); self.libA = np.zeros((N_SEEDS, ng))
        self.seen = np.zeros(len(group_of_cell), bool)
        # control sample for the control-PC sensitivity
        cg = [g for g, nm in enumerate(group_names) if nm == 'CTRL' or nm.endswith('|CTRL')] if ctrl_groups is None else ctrl_groups
        ci = np.flatnonzero(np.isin(group_of_cell, cg))
        rng = np.random.default_rng([SPLIT_BASE, 99, SRC_ID[src]])
        self.sample_idx = np.sort(rng.choice(ci, min(N_CTRL_SAMPLE, len(ci)), replace=False)) if len(ci) else ci
        self.sample = np.zeros((len(self.sample_idx), n_axis), np.float32)
        self.sample_lib = np.zeros(len(self.sample_idx))
        self._spos = {c: i for i, c in enumerate(self.sample_idx)}

    def add_dense(self, cell_idx, X, lib):
        """cell_idx: stream indices (n,), X: (n x n_axis) counts on the axis, lib: (n,) libraries."""
        code = self.code[cell_idx]; m = code >= 0
        if not m.any():
            return
        cell_idx, X, lib, code = cell_idx[m], X[m], lib[m], code[m]
        from scipy import sparse
        n = len(cell_idx)
        I = sparse.csr_matrix((np.ones(n), (code, np.arange(n))), shape=(self.ng, n))
        X = np.asarray(X, np.float64)
        self.full += I @ X
        self.lib += np.bincount(code, weights=lib, minlength=self.ng)
        for s in range(N_SEEDS):
            a = self.inA[s, cell_idx]
            Ia = sparse.csr_matrix((np.ones(int(a.sum())), (code[a], np.flatnonzero(a))), shape=(self.ng, n))
            self.A[s] += Ia @ X
            self.libA[s] += np.bincount(code[a], weights=lib[a], minlength=self.ng)
        self.seen[cell_idx] = True
        for k, c in enumerate(cell_idx):
            j = self._spos.get(int(c))
            if j is not None:
                self.sample[j] = X[k]; self.sample_lib[j] = lib[k]

    def add_csr(self, cell_idx, X, lib):
        """X: scipy CSR (n x n_axis)."""
        code = self.code[cell_idx]; m = code >= 0
        if not m.any():
            return
        self.add_dense(cell_idx[m], X[np.flatnonzero(m)].toarray(), lib[m])

    def save(self, extra):
        assert self.seen[self.code >= 0].all(), 'some grouped cells were never streamed'
        CACHE.mkdir(parents=True, exist_ok=True)
        np.savez(CACHE / f'{self.src}.npz', groups=np.array(self.gnames), n=self.n, nA=self.nA, full=self.full,
                 A=self.A.astype(np.float32), lib=self.lib, libA=self.libA, ctrl_sample=self.sample,
                 ctrl_sample_lib=self.sample_lib, **extra)


def build_groups(src, m, U, min_cells=MIN_CELLS):
    """Group ids per stream cell: U targets with >= min_cells (canonical name) + controls."""
    lab = m['label']; Uset = set(U)
    canon, _ = CM._canon()
    w = m.get('unit_cells')
    if src == 'FENG':
        line = np.array([x.split('|', 1)[0] if x else '' for x in lab]); tg = np.array([x.split('|', 1)[1] if x else '' for x in lab])
        ct = np.array([t if t in ('', 'CTRL') else canon(t) for t in tg])
        key = np.where((ct == 'CTRL') | np.isin(ct, list(Uset)), np.char.add(np.char.add(line, '|'), ct), '')
        key[lab == ''] = ''
    else:
        ct = np.array([t if t in ('', 'CTRL') else canon(t) for t in lab])
        key = np.where((ct == 'CTRL') | np.isin(ct, list(Uset)), ct, ''); key[lab == ''] = ''
    keys, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, weights=w) if w is not None else np.bincount(inv)
    ok = [(k != '') and (cnt[i] >= min_cells or k.endswith('CTRL')) for i, k in enumerate(keys)]
    names = [k for k, o in zip(keys, ok) if o]
    gid = {k: i for i, k in enumerate(names)}
    code = np.array([gid.get(k, -1) for k in keys])[inv]
    return code, names


def colmap(genes, axis):
    om = CM.to_official(genes); pos = {g: i for i, g in enumerate(axis)}
    return np.array([pos.get(int(o), -1) for o in om])


def _agg_cols(X, cm, n_axis):
    """Sum source columns into axis columns (duplicates summed). X dense (n x n_src)."""
    keep = np.flatnonzero(cm >= 0)
    if len(np.unique(cm[keep])) == len(keep):
        out = np.zeros((X.shape[0], n_axis), np.float64); out[:, cm[keep]] = X[:, keep]; return out
    from scipy import sparse
    M = sparse.csr_matrix((np.ones(len(keep)), (keep, cm[keep])), shape=(len(cm), n_axis))
    return np.asarray(X @ M)


# ============================================================ streamers

def extract(src):
    t0 = time.time()
    P = json.loads((CACHE / 'plan.json').read_text()); U, axis = P['universe'], P['gene_axis']; GA = len(axis)
    m = META[src](src)
    code, names = build_groups(src, m, U)
    acc = Acc(src, code, names, GA, unit_cells=m.get('unit_cells'))
    log(event='groups', src=src, groups=len(names), cells_used=int((code >= 0).sum()), secs=round(time.time() - t0, 1))
    extra = {}
    if src in RN_FILE:
        f = RAW / f'scperturb/{RN_FILE[src]}.h5ad'; cm = colmap(m['genes'], axis)
        with h5py.File(f, 'r') as h:
            X = h['X']; B = X.chunks[0] if X.chunks else 4000; B = max(B, 4000 // B * B)
            need = code >= 0
            for r0 in range(0, X.shape[0], B):
                r1 = min(X.shape[0], r0 + B)
                if not need[r0:r1].any():
                    continue
                blk = X[r0:r1, :]; sel = np.flatnonzero(need[r0:r1])
                acc.add_dense(r0 + sel, _agg_cols(blk[sel], cm, GA), m['lib'][r0 + sel])
                if (r0 // B) % 50 == 0:
                    log(event='rows', src=src, done=r1, of=X.shape[0], secs=round(time.time() - t0, 1))
    elif src in ('HCT116', 'HEK293T'):
        import pyarrow.parquet as pq
        cm = colmap(m['genes'], axis); off = 0
        for i, fn in enumerate(m['files']):
            nr = m['file_rows'][i]; c = code[off:off + nr]
            if (c >= 0).any():
                t = pq.read_table(fn, columns=['gene_token_id', 'gene_expression'])
                tok = t.column('gene_token_id').combine_chunks(); ex = t.column('gene_expression').combine_chunks()
                offs = tok.offsets.to_numpy(); tv = tok.values.to_numpy(); ev = ex.values.to_numpy()
                sel = np.flatnonzero(c >= 0)
                from scipy import sparse
                lens = offs[sel + 1] - offs[sel]
                ii = np.repeat(offs[sel], lens) + (np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens))
                rows = np.repeat(np.arange(len(sel)), lens); cols = cm[tv[ii]]; ok = cols >= 0
                Xs = sparse.csr_matrix((ev[ii][ok], (rows[ok], cols[ok])), shape=(len(sel), GA))
                acc.add_dense(off + sel, Xs.toarray(), m['lib'][off + sel])
            off += nr
            if i % 20 == 0:
                log(event='file', src=src, i=i, of=len(m['files']), secs=round(time.time() - t0, 1))
    elif src == 'KOLF':
        cm = colmap(m['genes'], axis)
        need = code >= 0
        with h5py.File(RAW / 'kolf/KOLF_Pan_Genome_QC_Filtered.h5ad', 'r') as h:
            L = h['layers/counts']; ip = L['indptr'][:]
            cols = np.flatnonzero(cm >= 0)
            # per-cell axis, but streaming genes: accumulate group x axis via bincount on (group, axis col)
            key_full = code.astype(np.int64)
            for b0 in range(0, len(cols), 200):
                bc = cols[b0:b0 + 200]
                for j in bc:
                    lo, hi = int(ip[j]), int(ip[j + 1])
                    if hi == lo:
                        continue
                    idx = L['indices'][lo:hi]; dat = L['data'][lo:hi].astype(np.float64)
                    cg = key_full[idx]; ok = cg >= 0
                    idx, dat, cg = idx[ok], dat[ok], cg[ok]
                    a = cm[j]
                    acc.full[:, a] += np.bincount(cg, weights=dat, minlength=acc.ng)
                    for s in range(N_SEEDS):
                        ia = acc.inA[s, idx]
                        acc.A[s, :, a] += np.bincount(cg[ia], weights=dat[ia], minlength=acc.ng)
                    sp = np.isin(idx, acc.sample_idx)
                    if sp.any():
                        acc.sample[np.searchsorted(acc.sample_idx, idx[sp]), a] += dat[sp]
                log(event='genes', src=src, done=b0 + len(bc), of=len(cols), secs=round(time.time() - t0, 1))
        lib = m['lib']; keep = need
        acc.lib = np.bincount(code[keep], weights=lib[keep], minlength=acc.ng)
        for s in range(N_SEEDS):
            k2 = keep & acc.inA[s]
            acc.libA[s] = np.bincount(code[k2], weights=lib[k2], minlength=acc.ng)
        acc.sample_lib = lib[acc.sample_idx]; acc.seen[:] = True
    elif src == 'DLD1':
        cm = colmap(m['genes'], axis); need = code >= 0
        with h5py.File(RAW / 'dld1_gse337988/GSE337988_sublib2_processed_objects_Low_assays.h5', 'r') as h:
            X = h['assay001']; B = 2000
            for r0 in range(0, X.shape[0], B):
                r1 = min(X.shape[0], r0 + B)
                if not need[r0:r1].any():
                    continue
                blk = X[r0:r1, :]; sel = np.flatnonzero(need[r0:r1])
                acc.add_dense(r0 + sel, _agg_cols(blk[sel].astype(np.float64), cm, GA), m['lib'][r0 + sel])
                if (r0 // B) % 50 == 0:
                    log(event='rows', src=src, done=r1, of=X.shape[0], secs=round(time.time() - t0, 1))
    elif src == 'H1GW':
        off = 0; cm = None
        for i, fn in enumerate(H1GW_FILES()):
            call, E, g = _h1gw_lane(fn)
            if cm is None:
                cm = colmap(g, axis)
            nr = len(call); c = code[off:off + nr]; sel = np.flatnonzero(c >= 0)
            lib = np.asarray(E.sum(1)).ravel().astype(np.float64)
            Es = E[sel]
            keep = np.flatnonzero(cm >= 0)
            from scipy import sparse
            M = sparse.csr_matrix((np.ones(len(keep)), (keep, cm[keep])), shape=(len(cm), GA))
            acc.add_dense(off + sel, (Es @ M).toarray(), lib[sel])
            off += nr
            if i % 10 == 0:
                log(event='lane', src=src, i=i, secs=round(time.time() - t0, 1))
    elif src == 'VCC_H1':
        cm = colmap(m['genes'], axis)
        from scipy import sparse
        keep = np.flatnonzero(cm >= 0)
        M = sparse.csr_matrix((np.ones(len(keep)), (keep, cm[keep])), shape=(len(cm), GA))
        with h5py.File(RAW / 'vcc2025/adata_Training.h5ad', 'r') as h:
            X = h['X']; ip = X['indptr'][:]; ncol = int(X.attrs['shape'][1]); B = 20000
            for r0 in range(0, len(ip) - 1, B):
                r1 = min(len(ip) - 1, r0 + B)
                lo, hi = int(ip[r0]), int(ip[r1])
                Xs = sparse.csr_matrix((X['data'][lo:hi], X['indices'][lo:hi], ip[r0:r1 + 1] - lo), shape=(r1 - r0, ncol))
                lib = np.asarray(Xs.sum(1)).ravel().astype(np.float64)
                sel = np.flatnonzero(code[r0:r1] >= 0)
                acc.add_dense(r0 + sel, (Xs[sel] @ M).toarray(), lib[sel])
    elif src == 'CD4':
        cm = colmap(m['genes'], axis)
        from scipy import sparse
        keep = np.flatnonzero(cm >= 0)
        M = sparse.csr_matrix((np.ones(len(keep)), (keep, cm[keep])), shape=(len(cm), GA))
        need = code >= 0
        with h5py.File(RAW / 'cd4/GWCD4i.pseudobulk_merged.h5ad', 'r') as h:
            X = h['X']; ip = X['indptr'][:]; ncol = int(X.attrs['shape'][1]); B = 2000
            for r0 in range(0, len(ip) - 1, B):
                r1 = min(len(ip) - 1, r0 + B); sel = np.flatnonzero(need[r0:r1])
                if not len(sel):
                    continue
                lo, hi = int(ip[r0]), int(ip[r1])
                Xs = sparse.csr_matrix((X['data'][lo:hi], X['indices'][lo:hi], ip[r0:r1 + 1] - lo), shape=(r1 - r0, ncol))
                acc.add_dense(r0 + sel, (Xs[sel] @ M).toarray(), m['lib'][r0 + sel])
                if (r0 // B) % 20 == 0:
                    log(event='units', src=src, done=r1, of=len(ip) - 1, secs=round(time.time() - t0, 1))
        extra['note'] = np.array('CD4: units = run x donor x guide x condition pseudobulks; halves are unit-level')
    elif src == 'FENG':
        genes = feng_genes(); cm = colmap(genes, axis)
        used = np.flatnonzero(code >= 0); cu = code[used]
        lib_used = np.zeros(len(used))
        Xbuf = []; Gbuf = []
        proc = subprocess.Popen(['pigz', '-dc', str(RAW / 'feng/TargetedScreen_RNA-UMI-Counts.csv.gz')], stdout=subprocess.PIPE,
                                bufsize=1 << 24)
        proc.stdout.readline()
        sp = np.searchsorted(used, acc.sample_idx)
        for gi, line in enumerate(proc.stdout):
            k = line.index(b',')
            v = np.fromstring(line[k + 1:], dtype=np.float64, sep=',')[used]
            lib_used += v
            a = cm[gi]
            if a >= 0:
                acc.full[:, a] += np.bincount(cu, weights=v, minlength=acc.ng)
                for s in range(N_SEEDS):
                    ia = acc.inA[s, used]
                    acc.A[s, :, a] += np.bincount(cu[ia], weights=v[ia], minlength=acc.ng)
                acc.sample[:, a] += v[sp]
            if gi % 2000 == 0:
                log(event='generow', src=src, i=gi, secs=round(time.time() - t0, 1))
        proc.wait()
        acc.lib = np.bincount(cu, weights=lib_used, minlength=acc.ng)
        for s in range(N_SEEDS):
            ia = acc.inA[s, used]
            acc.libA[s] = np.bincount(cu[ia], weights=lib_used[ia], minlength=acc.ng)
        acc.sample_lib = lib_used[sp]; acc.seen[:] = True
    else:
        raise ValueError(src)
    acc.save(extra)
    log(event='done', src=src, secs=round(time.time() - t0, 1), groups=len(names),
        median_cells=float(np.median(acc.n[[i for i, n in enumerate(names) if not n.endswith('CTRL')]])),
        ctrl_cells=float(acc.n[[i for i, n in enumerate(names) if n.endswith('CTRL')]].sum()))


if __name__ == '__main__':
    if sys.argv[1] == 'plan':
        plan()
    elif sys.argv[1] == 'extract':
        for s in sys.argv[2:]:
            extract(s)
