"""Task 12 (RQ1) -- cell-level split-half extraction for the harness-panel targets NOT covered by the Task 9 cache.

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 12 (definition 6), Amendment A3.
The Task 9 cache (results/rq2_gene_cache/, read-only) holds only targets in its universe U (>= 20 cells in an
essential screen). The H1, OR_HCT116X and KOLF panels' targets are mostly outside U (H1 201/300, OR 263/263,
KOLF 279/279). This script re-uses rq2_gene_extract's readers verbatim (META, build_groups, extract, Acc) with
  * universe X = (canonical targets of the H1, OR_HCT116X, KOLF truth panels) minus U   (474 targets)
  * the SAME gene axis as the Task 9 cache (plan.json gene_axis, 6,109 official-axis genes)
  * output dir results/rq1_gene_cache/ (the Task 9 cache is never written)
  * sources K562 (GWPS), HCT116, HEK293T, CD4, KOLF, VCC_H1 (FENG covers only 9 of the 474 targets with a
    >= 20-cell donor group and its reader takes ~1 h: skipped, logged as a deviation; RPE1/HepG2/Jurkat are
    essential screens, < 20 cells for every target outside U by construction of U; H1GW, DLD1 are not
    in any panel's pack families).
CONTROL-HALF ALIGNMENT (required, not cosmetic): a source's targets from the Task 9 cache and from this cache must
share ONE control split, otherwise half-A and half-B responses of the same source share control cells and the
cross-half products are biased. The control group's split therefore uses the Task 9 group index of 'CTRL'
(default_rng([SPLIT_BASE, seed, SRC_ID[src], g_T9])); `verify` checks that the control sums (full and half A,
all 5 seeds) are identical to the Task 9 cache. Target groups use index 100000 + own index (independent halves).
Run: python rq1_gene_extract.py plan | extract SRC [SRC ...] | extract_lean HCT116 HEK293T | verify
(Orion lines use extract_lean: row-batched parquet reads, added 2026-11-05 after the main session's memory limit)
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC')
sys.path.insert(0, str(ROOT / 'output/data-analysis'))
import common as CM
import rq2_gene_extract as RG

T9 = ROOT / 'output/data-analysis/results/rq2_gene_cache'
OUT = ROOT / 'output/data-analysis/results/rq1_gene_cache'
SOURCES = ['VCC_H1', 'CD4', 'KOLF', 'K562', 'HCT116', 'HEK293T']
TARGET_OFFSET = 100000


def t9_groups(src):
    with np.load(T9 / f'{src}.npz', allow_pickle=False) as z:
        return [str(x) for x in z['groups']]


class AlignedAcc(RG.Acc):
    def __init__(self, src, group_of_cell, group_names, n_axis, unit_cells=None, ctrl_group=None, ctrl_groups=None):
        super().__init__(src, group_of_cell, group_names, n_axis, unit_cells, ctrl_group, ctrl_groups)
        g9 = {n: i for i, n in enumerate(t9_groups(src))}
        w = self.w
        keep = group_of_cell >= 0
        self.inA[:] = False
        for g, name in enumerate(group_names):
            ii = np.flatnonzero(group_of_cell == g)
            gi = g9[name] if (name == 'CTRL' or name.endswith('|CTRL')) else TARGET_OFFSET + g
            for s in range(RG.N_SEEDS):
                rng = np.random.default_rng([RG.SPLIT_BASE, s, RG.SRC_ID[src], gi])
                if unit_cells is None:
                    perm = rng.permutation(len(ii))
                    self.inA[s, ii] = perm < len(ii) // 2
                else:
                    order = rng.permutation(len(ii)); ca = cb = 0.0
                    for j in order:
                        if ca <= cb:
                            self.inA[s, ii[j]] = True; ca += w[ii[j]]
                        else:
                            cb += w[ii[j]]
        self.nA = np.stack([np.bincount(group_of_cell[keep], weights=(w * self.inA[s])[keep], minlength=self.ng)
                            for s in range(RG.N_SEEDS)])


def plan():
    OUT.mkdir(parents=True, exist_ok=True)
    P = json.loads((T9 / 'plan.json').read_text())
    U = set(P['universe'])
    X = set()
    for p in ['H1', 'OR_HCT116X', 'KOLF']:
        X |= {RG.canon_target(t) for t in set(CM.labels(p)) - {'non-targeting'}}
    X = sorted(X - U)
    (OUT / 'plan.json').write_text(json.dumps(dict(universe=X, gene_axis=P['gene_axis'], n_universe=len(X),
                                                   n_axis=len(P['gene_axis']), source='rq1_gene_extract.plan',
                                                   note='universe = H1/OR_HCT116X/KOLF panel targets not in the Task 9 universe')))
    RG.log(event='plan', n_extra=len(X), n_axis=len(P['gene_axis']))


def extract_orion_lean(src, batch_rows=4096):
    """Memory-lean Orion reader (HCT116 / HEK293T): same cells, groups, axis mapping and half split as
    rq2_gene_extract.extract, but the parquet files are streamed in row batches of `batch_rows` cells instead of
    whole files (peak RSS kept well below 8 GB on the shared machine)."""
    import time
    import pyarrow.parquet as pq
    from scipy import sparse
    t0 = time.time()
    P = json.loads((OUT / 'plan.json').read_text()); U, axis = P['universe'], P['gene_axis']; GA = len(axis)
    m = RG.META[src](src)
    code, names = RG.build_groups(src, m, U)
    acc = AlignedAcc(src, code, names, GA)
    RG.log(event='groups', src=src, groups=len(names), cells_used=int((code >= 0).sum()), secs=round(time.time() - t0, 1))
    cm = RG.colmap(m['genes'], axis)
    off = 0
    for i, fn in enumerate(m['files']):
        nr = m['file_rows'][i]
        if (code[off:off + nr] >= 0).any():
            pf = pq.ParquetFile(fn)
            assert pf.metadata.num_rows == nr
            o = 0
            for b in pf.iter_batches(batch_size=batch_rows, columns=['gene_token_id', 'gene_expression']):
                n = b.num_rows
                sel = np.flatnonzero(code[off + o:off + o + n] >= 0)
                if len(sel):
                    bt = b.take(sel)
                    tok = bt.column(0); ex = bt.column(1)
                    offs = tok.offsets.to_numpy(); tv = tok.values.to_numpy(); ev = ex.values.to_numpy()
                    eo = ex.offsets.to_numpy()
                    assert np.array_equal(offs - offs[0], eo - eo[0])
                    tv = tv[offs[0]:offs[-1]]; ev = ev[eo[0]:eo[-1]]
                    lens = np.diff(offs)
                    rows = np.repeat(np.arange(len(sel)), lens); cols = cm[tv]; ok = cols >= 0
                    Xs = sparse.csr_matrix((ev[ok].astype(np.float64), (rows[ok], cols[ok])), shape=(len(sel), GA))
                    acc.add_dense(off + o + sel, Xs.toarray(), m['lib'][off + o + sel])
                o += n
            assert o == nr
        off += nr
        if i % 20 == 0:
            RG.log(event='file', src=src, i=i, of=len(m['files']), secs=round(time.time() - t0, 1))
    acc.save({})
    RG.log(event='done', src=src, secs=round(time.time() - t0, 1), groups=len(names))


def verify():
    rep = {}
    for src in SOURCES:
        f = OUT / f'{src}.npz'
        if not f.exists():
            rep[src] = 'missing'; continue
        a, b = np.load(f), np.load(T9 / f'{src}.npz', mmap_mode='r')
        ga, gb = [str(x) for x in a['groups']], [str(x) for x in b['groups']]
        ia, ib = ga.index('CTRL'), gb.index('CTRL')
        r = dict(n=float(a['n'][ia]) == float(b['n'][ib]),
                 nA=bool(np.array_equal(a['nA'][:, ia], b['nA'][:, ib])),
                 full_maxabs=float(np.abs(a['full'][ia] - b['full'][ib]).max()),
                 A_maxabs=float(np.abs(a['A'][:, ia].astype(np.float64) - b['A'][:, ib].astype(np.float64)).max()),
                 lib=bool(np.isclose(a['lib'][ia], b['lib'][ib], rtol=1e-12)),
                 libA=bool(np.allclose(a['libA'][:, ia], b['libA'][:, ib], rtol=1e-12)),
                 n_target_groups=len(ga) - 1)
        r['pass'] = bool(r['n'] and r['nA'] and r['full_maxabs'] == 0 and r['A_maxabs'] == 0 and r['lib'] and r['libA'])
        rep[src] = r
        RG.log(event='verify', src=src, **r)
    (OUT / 'verify.json').write_text(json.dumps(rep, indent=1))


if __name__ == '__main__':
    RG.CACHE = OUT
    RG.Acc = AlignedAcc
    cmd = sys.argv[1]
    if cmd == 'plan':
        plan()
    elif cmd == 'extract':
        for s in (sys.argv[2:] or SOURCES):
            RG.extract(s)
    elif cmd == 'extract_lean':
        for s in sys.argv[2:]:
            extract_orion_lean(s)
    elif cmd == 'verify':
        verify()
