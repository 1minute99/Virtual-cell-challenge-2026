"""Task 9 step 1 (Amendment A2): C_public Count 2 (donor lines), computed BEFORE any genetic C*_min.

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 9, definition 5, Amendment A2.
Rule (A2): every distinct cell line with a usable per-line response = >= 300 targets with >= 20 cells each,
counted separately; kolf_2/kolf_3 merged with KOLF2.1J (same line; KOLF2.1J is already in Count 1).
Count 2 = Count 1 (families.json C_public = 11) - 1 (the pooled HipSci/FENG family) + qualifying non-kolf HipSci lines.

ASSUMPTIONS
  C1 'usable cells' = cells present in the Feng targeted-screen RNA count matrix with a single-target guide call,
     i.e. the n_cells of pseudobulk/feng_targetedscreen_target_line (built from the RNA matrix columns matched to
     metadata). Cross-checked against a direct count from TargetedScreen_Cell-Metadata.tsv.gz.
  C2 a HipSci line = the Cell_Line identifier (donor_clone, e.g. eipl_1). Two clones of one donor are two lines
     under A2's literal wording ("every distinct cell line"); a donor-collapsed count is reported as sensitivity.
  C3 only the targeted screen (19 lines x 445 targets) can qualify; the genome-wide arms have 1-2 cells per
     line-target (families.json), checked here from the _line stores.
No randomness.
"""
import json, sys, datetime
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path('/home/wonmin-kim/AIVC'); V = ROOT / 'vcc2026_data'
RES = ROOT / 'output/data-analysis/results'
MIN_T, MIN_CELLS = 300, 20
KOLF = {'kolf_2', 'kolf_3'}


def per_line_counts(store):
    o = pd.read_csv(V / 'pseudobulk' / store / 'obs.csv')
    p = o[~o.is_control.astype(bool)]
    out = {}
    for ln, d in p.groupby('donor'):
        out[ln] = dict(n_targets=int(d.target_symbol.nunique()),
                       n_targets_ge20=int((d.groupby('target_symbol').n_cells.sum() >= MIN_CELLS).sum()),
                       median_cells=float(d.n_cells.median()))
    c = o[o.is_control.astype(bool)].groupby('donor').n_cells.sum()
    for ln in out:
        out[ln]['n_control_cells'] = int(c.get(ln, 0))
    return out


def metadata_crosscheck():
    sys.path.insert(0, str(V / 'scripts'))
    import ingest_feng_targeted as IFT
    m = pd.read_csv(V / 'raw/feng/TargetedScreen_Cell-Metadata.tsv.gz', sep='\t')
    call = m.Guide_Call.astype(str)
    tgt = call.map(IFT.IF.guide_to_gene)
    ntc = tgt.map(lambda s: bool(IFT.IF.NTC_EXACT.match(str(s))))
    un = call.str.lower().isin(['unassigned', 'nan'])
    m['target'] = np.where(ntc, 'non-targeting', np.where(un, None, tgt))
    m = m[m.target.notna() & (m.target != 'non-targeting')]
    g = m.groupby(['Cell_Line', 'target']).size()
    return {ln: int((s >= MIN_CELLS).sum()) for ln, s in g.groupby(level=0)}


def main():
    fam = json.loads((RES / 'families.json').read_text())
    c1 = int(fam['C_public'])
    T = per_line_counts('feng_targetedscreen_target_line')
    meta = metadata_crosscheck()
    gw = {s: max(v['n_targets_ge20'] for v in per_line_counts(s).values())
          for s in ('feng_genomewidescreen_fitnessgenes_target_line', 'feng_genomewidescreen_nonfitnessgenes_target_line')}
    qual = sorted(l for l, v in T.items() if v['n_targets_ge20'] >= MIN_T)
    qual_nonkolf = [l for l in qual if l not in KOLF]
    donors = sorted({l.split('_')[0] for l in qual_nonkolf})
    count2 = c1 - 1 + len(qual_nonkolf)
    out = dict(
        task='Task 9 step 1 (A2) -- computed before any genetic C*_min estimate',
        computed=datetime.date.today().isoformat(),
        rule='>= %d targets with >= %d cells each (Feng targeted screen, per line); kolf_2/kolf_3 merged with KOLF2.1J' % (MIN_T, MIN_CELLS),
        C_public_count1=c1,
        C_public_count2=count2,
        formula='Count1 (%d) - 1 (pooled HipSci family) + %d qualifying non-kolf HipSci lines' % (c1, len(qual_nonkolf)),
        qualifying_lines=qual, qualifying_nonkolf_lines=qual_nonkolf,
        failing_lines={l: T[l]['n_targets_ge20'] for l in T if l not in qual},
        kolf_lines_merged_with_KOLF2_1J=sorted(KOLF & set(qual)),
        per_line=T,
        metadata_crosscheck_n_targets_ge20=meta,
        metadata_crosscheck_agrees_on_qualification=bool(all((meta.get(l, 0) >= MIN_T) == (l in qual) for l in T)),
        genomewide_arms_max_targets_ge20_any_line=gw,
        sensitivity_donor_collapsed=dict(C_public=c1 - 1 + len(donors), donors=donors,
                                         note='clones of the same donor counted once (not the A2 wording; reported only)'),
    )
    f = RES / 'rq2_params_gene.json'
    prev = json.loads(f.read_text()) if f.exists() else {}
    prev['count2'] = out
    f.write_text(json.dumps(prev, indent=1))
    print(json.dumps({k: out[k] for k in ('C_public_count1', 'C_public_count2', 'qualifying_nonkolf_lines', 'failing_lines',
                                          'metadata_crosscheck_agrees_on_qualification', 'genomewide_arms_max_targets_ge20_any_line',
                                          'sensitivity_donor_collapsed')}, indent=1))


if __name__ == '__main__':
    main()
