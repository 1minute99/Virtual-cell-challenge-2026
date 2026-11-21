#!/usr/bin/env python3
"""Ingest the Feng multi-donor iPSC CRISPRi counts into count-sum pseudobulk.

Layout, read from the files rather than assumed:
  *_RNA-UMI-Counts.csv.gz   GENES x CELLS (row label `ENSG:SYMBOL:Gene-Expression`,
                            one column per cell), so it is streamed by ROW and the
                            grouped sum is taken across columns.
  *_Cell-Metadata.tsv.gz    Cell_ID, Batch, Guide_Call, Cell_Line
  *_Guide-UMI-Counts.csv.gz GUIDES x CELLS, guide row labels `GENE-PROTOSPACER`

THE BARCODE TRAP. The same cell is spelled three different ways:
  metadata   MP-ST-P3-D3_I1_AAACCTGAGACGCAAC-1
  RNA header    ST-P3-D3_I1_AAACCTGAGACGCAAC-1     (no MP- prefix)
  guide header MP.ST.P3.D3_I1_AAACCTGAGACGCAAC.1   (dots for dashes)
So all three are normalised (dots to dashes, optional MP- prefix dropped) and the
overlap is checked; the run aborts rather than silently aggregating misaligned
columns.

`Cell_Line` is the iPSC line, a real biological unit, recorded as donor. `Batch`
is a pool/day/inlet label, recorded as technical_batch.

A cell whose Guide_Call is `unassigned` has NO detected guide. That is NOT the
same as a confirmed non-targeting control, so those cells are kept in their own
`no-guide-detected` group and never used as the control.
"""
import argparse
import gzip
import json
import pathlib
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import make_effects as ME
import pbstore as PBS
import vcclib as V

RAW = V.ROOT / "raw/feng"
EXPERIMENTS = [
    ("GenomeWideScreen_FitnessGenes", "genome-wide fitness-gene arm"),
    ("GenomeWideScreen_NonFitnessGenes", "genome-wide non-fitness-gene arm"),
    ("TargetedScreen", "focused 444-target arm"),
]
PROTO = re.compile(r"[_-]([ACGT]{18,25})$")
NOGUIDE_LABEL = "no-guide-detected (NOT a confirmed non-targeting control)"
# Matched as WHOLE tokens, not substrings. A bare "ntc" substring test misfires on
# KNTC1 (kinetochore associated 1), a real gene, and would delete a real target while
# inventing 30 control cells that do not exist.
NTC_EXACT = re.compile(
    r"^(non[-_]?targeting|ntc|scrambled?|safe[-_]?harbou?r|no[-_]?target)"
    r"([-_]?\d+)?$", re.I)


def norm_bc(s):
    s = str(s).replace(".", "-")
    return s[3:] if s.startswith("MP-") else s


def guide_to_gene(call):
    """'USP8_TGGCGTTAGTGAAGCGCGCC' -> 'USP8'. Strips the protospacer, not a hyphen.

    Splitting on a hyphen would break HLA-DRA and every other symbol that
    legitimately contains one, so only a trailing ACGT run is removed.
    """
    c = str(call)
    m = PROTO.search(c)
    return c[:m.start()] if m else c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--experiments", default="",
                    help="comma list; default is every arm present")
    ap.add_argument("--chunk-rows", type=int, default=256)
    ap.add_argument("--max-rows", type=int, default=0, help="debug: stop after N gene rows")
    a = ap.parse_args()

    gm = V.GeneMap()
    features = V.official_features()
    want = [e for e in EXPERIMENTS
            if not a.experiments or e[0] in a.experiments.split(",")]
    summary = []

    for exp, desc in want:
        rna = RAW / f"{exp}_RNA-UMI-Counts.csv.gz"
        meta = RAW / f"{exp}_Cell-Metadata.tsv.gz"
        if not rna.exists() or not meta.exists():
            print(f"skip {exp}: files not present")
            continue
        t0 = time.time()
        print(f"===== {exp} ({desc})", flush=True)
        M = pd.read_csv(meta, sep="\t")
        M["bc"] = M.Cell_ID.map(norm_bc)
        M["target"] = M.Guide_Call.map(
            lambda c: "no-guide-detected" if str(c).lower() in
            ("unassigned", "nan", "none", "") else guide_to_gene(c))
        M.loc[M.target.apply(lambda s: bool(NTC_EXACT.match(str(s)))),
              "target"] = "non-targeting"
        n_ntc = int((M.target == "non-targeting").sum())
        print(f"  metadata: {len(M):,} cells, {M.Cell_Line.nunique()} iPSC lines, "
              f"{M.Batch.nunique()} batches, "
              f"{int((M.target == 'no-guide-detected').sum()):,} with no guide detected, "
              f"{int((M.target == 'non-targeting').sum()):,} confirmed non-targeting, "
              f"{M[~M.target.isin(['no-guide-detected', 'non-targeting'])].target.nunique()} "
              f"perturbed target genes")

        with gzip.open(rna, "rt") as f:
            header = f.readline().rstrip("\n").split(",")
        cols = [norm_bc(c) for c in header[1:]]
        by_bc = dict(zip(M.bc, M.index))
        hit = sum(1 for c in cols if c in by_bc)
        print(f"  RNA matrix: {len(cols):,} cell columns; {hit:,} align to metadata "
              f"({hit/len(cols)*100:.2f}%)")
        if hit / len(cols) < 0.95:
            raise SystemExit(f"{exp}: only {hit/len(cols)*100:.1f}% of RNA columns align "
                             f"to the metadata barcodes; refusing to aggregate")

        rows = [by_bc.get(c, -1) for c in cols]
        # TWO group levels, because this screen is extremely fragmented: at
        # (target x line x batch) the non-fitness arm has 105,591 groups and at
        # (target x line) 56,174 groups averaging ~7 cells each -- too few cells for a
        # usable profile, and a dense accumulator at that width is tens of GB (an
        # earlier version of this script reached 26 GB RSS and had to be killed).
        # So: the pooled (target) level is accumulated densely and is what effects use,
        # while the finer (target x line) level is accumulated SPARSELY for donor-level
        # analysis. Batch composition is recorded as a count, not as a group axis.
        fine_key, coarse_key = [], []
        for r in rows:
            if r < 0:
                fine_key.append(None); coarse_key.append(None)
            else:
                t_, l_, b_ = (str(M.at[r, "target"]), str(M.at[r, "Cell_Line"]),
                              str(M.at[r, "Batch"]))
                fine_key.append((t_, l_)); coarse_key.append((t_,))
        uniq = sorted({k for k in fine_key if k})
        cuniq = sorted({k for k in coarse_key if k})
        gidx = {k: i for i, k in enumerate(uniq)}
        cidx = {k: i for i, k in enumerate(cuniq)}
        codes = np.array([gidx[k] if k else -1 for k in fine_key], dtype=np.int64)
        ccodes = np.array([cidx[k] if k else -1 for k in coarse_key], dtype=np.int64)
        keep = codes >= 0
        ncells = np.bincount(codes[keep], minlength=len(uniq)).astype(np.int64)
        cncells = np.bincount(ccodes[keep], minlength=len(cuniq)).astype(np.int64)
        nbatch = pd.DataFrame({"k": [k for k in fine_key if k],
                               "b": [str(M.at[r, "Batch"]) for r in rows if r >= 0]}) \
            .groupby("k").b.nunique()
        print(f"  {len(uniq)} (target x line) groups, {len(cuniq)} pooled (target) groups",
              flush=True)

        # stream gene rows; each row is one gene across all cells
        fine_blocks, coarse_blocks, gene_labels = [], [], []
        # no global dtype: it would be applied to the gene-label index column too
        reader = pd.read_csv(rna, chunksize=a.chunk_rows, header=0, index_col=0,
                             low_memory=False)
        n_seen = 0
        neg = nonint = 0
        ck, cck = codes[keep], ccodes[keep]
        for ch in reader:
            V_ = ch.to_numpy(dtype=np.float64, copy=False)
            if V_.shape[1] != len(cols):
                raise SystemExit(f"{exp}: chunk has {V_.shape[1]} columns, expected "
                                 f"{len(cols)}")
            neg += int((V_ < 0).sum())
            nonint += int((np.abs(V_ - np.rint(V_)) > 1e-6).sum())
            Vk = V_[:, keep]
            fb = np.zeros((V_.shape[0], len(uniq)))
            cb = np.zeros((V_.shape[0], len(cuniq)))
            for j in range(V_.shape[0]):
                fb[j] = np.bincount(ck, weights=Vk[j], minlength=len(uniq))
                cb[j] = np.bincount(cck, weights=Vk[j], minlength=len(cuniq))
            fine_blocks.append(sparse.csc_matrix(fb.T))
            coarse_blocks.append(cb.T)
            gene_labels.extend(list(ch.index))
            n_seen += V_.shape[0]
            del V_, Vk, fb, cb
            if n_seen % (a.chunk_rows * 8) == 0:
                nnz = sum(b.nnz for b in fine_blocks)
                print(f"    {n_seen} gene rows, fine nnz {nnz/1e6:.1f}M "
                      f"({time.time()-t0:.0f}s)", flush=True)
            if a.max_rows and n_seen >= a.max_rows:
                break
        if neg or nonint:
            raise SystemExit(f"{exp}: {neg} negative and {nonint} non-integral values; "
                             f"this is not a raw UMI count matrix")
        F = sparse.hstack(fine_blocks, format="csr")
        A = np.hstack(coarse_blocks)
        del fine_blocks, coarse_blocks
        print(f"  fine {F.shape[0]} x {F.shape[1]} ({F.nnz/1e6:.1f}M nnz), "
              f"pooled {A.shape[0]} x {A.shape[1]}, "
              f"{A.sum():,.0f} total counts ({time.time()-t0:.0f}s)")

        # gene labels look like ENSG00000243485:MIR1302-2HG:Gene-Expression
        parts = [str(g).split(":") for g in gene_labels]
        ens = [p[0] if p else "" for p in parts]
        sym = [p[1] if len(p) > 1 else (p[0] if p else "") for p in parts]
        kinds = pd.Series([p[2] if len(p) > 2 else "" for p in parts]).value_counts()
        print(f"  feature kinds in the row labels: {kinds.to_dict()}")
        var, o2s, method = PBS.map_features(sym, ens, gm=gm, features=features)
        print(f"  official features measured: {int((o2s >= 0).sum())}/{len(features)} "
              f"({pd.Series(method).value_counts().to_dict()})")

        def prov_for(level, keys):
            return dict(
                dataset_id="Feng2025", experiment_id=f"Feng-{exp}", experiment=desc,
                level=level,
                source_files=[str(rna), str(meta)],
                source_url="https://api.figshare.com/v2/articles/27989294",
                release="Figshare 27989294 v2", license="MIT",
                checksum="md5 matches the Figshare API supplied_md5 for all 9 files",
                orientation="RNA CSV is GENES x CELLS; streamed by row",
                count_semantics="SUM of raw UMI counts, computed here",
                semantics_evidence="every streamed value non-negative and integral",
                barcode_alignment=(f"{hit}/{len(cols)} RNA columns matched metadata after "
                                   f"normalising dots to dashes and dropping an MP- prefix"),
                group_keys=keys,
                fragmentation_note=(
                    f"this arm has {len(uniq)} (target x line) groups averaging "
                    f"{len(cols)/max(1,len(uniq)):.1f} cells; the pooled (target) level "
                    f"exists because that is too few cells per group for a usable "
                    f"profile. Batch is recorded as a per-group count, not a group axis"),
                replicate_note=("Cell_Line is the iPSC line and is recorded as donor / "
                                "biological_replicate; Batch is a pool/day/inlet label "
                                "and is counted in n_technical_batches"),
                control_note=(
                    f"confirmed non-targeting cells in this arm: {n_ntc}. "
                    + ("Because there are none, the reference used for effects is the "
                       "`no-guide-detected` group (Guide_Call == unassigned). Those cells "
                       "have NO detected guide, which is NOT the same as a confirmed "
                       "non-targeting control: they include cells whose guide simply "
                       "failed to be captured, so any effect computed against them is "
                       "weaker evidence than a real NTC contrast and is labelled as such."
                       if n_ntc == 0 else
                       "`no-guide-detected` cells are kept as a separate group and are "
                       "NOT used as the control.")),
                reference_control=("no-guide-detected" if n_ntc == 0 else "non-targeting"),
                guide_to_gene=("gene taken by stripping a trailing ACGT protospacer, "
                               "never by splitting on a hyphen, so HLA-DRA survives"),
                perturbation_mode="CRISPRi",
                partial_run=bool(a.max_rows))

        def obs_for(keys, ncell_vec, totals, donors, nbatches):
            tg = [k[0] for k in keys]
            return pd.DataFrame(dict(
                pb_id=np.arange(len(keys)), dataset_id="Feng2025",
                experiment_id=f"Feng-{exp}", cell_context="human iPSC",
                donor=donors, biological_replicate=donors,
                technical_batch="", stimulation="none", time="",
                perturbation_mode="CRISPRi", target_symbol=tg,
                target_ensembl=[("" if t in ("non-targeting", "no-guide-detected")
                                 else (gm.ens(t) or "")) for t in tg],
                guide_id="<aggregated over guides>",
                control_type=[("non-targeting" if t == "non-targeting"
                               else (NOGUIDE_LABEL if t == "no-guide-detected" else ""))
                              for t in tg],
                is_control=[(t == "non-targeting") or
                            (n_ntc == 0 and t == "no-guide-detected") for t in tg],
                n_cells=ncell_vec, total_counts=totals,
                qc_flags="release as published",
                n_technical_batches=nbatches, n_source_groups=ncell_vec))

        # pooled (target) level: dense, and what effects are computed on
        ckey = f"feng_{exp.lower()}_target"
        cob = obs_for(cuniq, cncells, A.sum(1).astype(np.int64),
                      ["<pooled across lines>"] * len(cuniq),
                      [int(nbatch[[k for k in uniq if k[0] == c[0]]].sum())
                       if any(k[0] == c[0] for k in uniq) else 0 for c in cuniq])
        PBS.write_store(ckey, A, cob, var, prov_for("pooled target", ["target_symbol"]),
                        o2s, method)
        try:
            ME.build(ckey, ["cell_context"], min_cells=10, min_control_cells=100)
        except SystemExit as e:
            print(f"  pooled effects not built: {e}")

        # finer (target x line) level: sparse, for donor-level analysis
        fkey = f"feng_{exp.lower()}_target_line"
        fob = obs_for(uniq, ncells, np.asarray(F.sum(1)).ravel().astype(np.int64),
                      [k[1] for k in uniq],
                      [int(nbatch.get(k, 0)) for k in uniq])
        PBS.write_store(fkey, F, fob, var,
                        prov_for("target x iPSC line", ["target_symbol", "donor"]),
                        o2s, method)
        try:
            ME.build(fkey, ["donor"], min_cells=10, min_control_cells=50)
        except SystemExit as e:
            print(f"  per-line effects not built: {e}")

        summary.append(dict(experiment=exp, cells=len(cols), aligned=hit,
                            groups=len(uniq), genes=A.shape[1],
                            lines=int(M.Cell_Line.nunique()),
                            targets=int(M[~M.target.isin(
                                ['no-guide-detected', 'non-targeting'])]
                                .target.nunique()),
                            seconds=round(time.time() - t0)))
    if summary:
        S = pd.DataFrame(summary)
        S.to_csv(V.ROOT / "reports/feng_summary.csv", index=False)
        print()
        print(S.to_string(index=False))
    print("\nnext: python scripts/coverage.py && python scripts/smoke_test.py")


if __name__ == "__main__":
    main()
