"""Task 4 of plan zero-shot-limits-context-complexity: freeze panels, source registry, leave-destination-out audit,
official record. Writes results/{panels.json, families.json, official_record.csv, g_cap.json, task04_checks.json}.

ASSUMPTIONS (stated; each is checked below and the outcome written to task04_checks.json):
  S1  The ERA packs' truth = fh_harness.PANELS paths (final_harness_2026-09-24/panels.json). Checked: pack target lists
      equal the truth files' target sets.
  S2  "Excluded sources" for destination P = union of what each frozen-pipeline component leaves out:
      M1 fusion + neighbour library (fh_harness.EXCL[P] via TM.predict_lfc(exclude=) / H.lib_delta; families with zero
      weight are inert), fusion family list (ml_h1fuse.OUTER[P] -> weight 0; = pack `families`), S4 part (ext_baselines.refit
      leave-own-biology-out + CONFLICT_GENE + the S4KEY extra exclusions). Parsed from source with `ast`, not imported.
  S3  A LEAK = a source measured in the destination's own cell line (incl. parental/isogenic line or same donor) or
      the same underlying dataset as the truth that is NOT excluded for that destination.
  S4  C_public (plan definition 5) counts distinct cell lines/types with a genome-scale or essential-scale (>= ~2,000
      targets) CRISPRi/KO Perturb-seq response among M1-registry families with non-zero weight in the frozen F2 spec
      (+ H1GW, the registry's same-line/stem source). Same line from different labs counts once.
  S5  Official submission count = distinct entry ids over all on-disk records, plus history-table rows without an id.
Deterministic; no random numbers."""
import ast, collections, csv, gzip, json, pathlib, re, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as C

RP, SC = C.RP, C.V / "scripts"; OUT = C.RESULTS; OUT.mkdir(parents=True, exist_ok=True)
checks = {}


def lit(path, name):
    """ast.literal_eval of a module-level assignment `name = <literal>` (sets allowed)."""
    for n in ast.parse(pathlib.Path(path).read_text()).body:
        tg = n.targets if isinstance(n, ast.Assign) else []
        for t in tg:
            if isinstance(t, ast.Name) and t.id == name:
                return ast.literal_eval(n.value)
            if isinstance(t, ast.Tuple):                       # A, B = ...
                for k, e in enumerate(t.elts):
                    if isinstance(e, ast.Name) and e.id == name:
                        return ast.literal_eval(n.value.elts[k])
    raise KeyError(name)


EXCL = lit(SC / "fh_harness.py", "EXCL"); W_NEW_H = lit(SC / "fh_harness.py", "W_NEW"); FAMS_H = lit(SC / "fh_harness.py", "FAMS")
OUTER = lit(SC / "ml_h1fuse.py", "OUTER"); S4KEY = lit(SC / "ml_h1fuse.py", "S4KEY"); S4_LINE = lit(SC / "sch_emit.py", "S4_LINE")
CONFLICT = lit(SC / "ext_baselines.py", "CONFLICT_GENE"); FAMILY = lit(SC / "tv2_model.py", "FAMILY")
BASE_NEW = lit(SC / "tv2_final.py", "BASE_NEW"); SOURCE_LINES = lit(SC / "tv2_final.py", "SOURCE_LINES")
SAME_LINE_SRC = lit(SC / "tv2_final.py", "SAME_LINE_SRC")
F2 = json.loads((RP / "final_2026-10-27/F1/F2_spec.json").read_text())
F2W = F2["parts"][1]["cfg"]["weights"]
# S4 extra exclusions by key (s4_improve tags: noKOLFFENG = KOLF+FENG biologies removed; orX = HCT116+HEK293T; orR = HCT116)
S4_EXTRA = {"Yorig": [], "F": [], "YnoKOLFFENG": ["KOLF", "FENG"], "YorXnoKOLFFENG": ["HEK293T", "KOLF", "FENG"],
            "YorRnoKOLFFENG": ["KOLF", "FENG"]}

# ---------- table_v2 (S4 training table): which biologies carry the kolf-background Feng lines ----------
tb = np.load(C.V / "reports/ext_contexts_2026-09-21/table_v2/bios.npy", allow_pickle=True).astype(str)
ts = np.load(C.V / "reports/ext_contexts_2026-09-21/table_v2/sources.npy", allow_pickle=True).astype(str)
kolf_rows = collections.Counter((s, b) for s, b in zip(ts, tb) if "kolf" in s.lower() or b.upper().startswith("FENGL"))
checks["table_v2_kolf_background_rows"] = {f"{s}|bio={b}": n for (s, b), n in sorted(kolf_rows.items())}
checks["CONFLICT_GENE_labels_present_in_table_v2"] = {k: [x for x in v if x in set(tb)] for k, v in CONFLICT.items()}

# ---------- Feng (HipSci) metadata: share of kolf_2 / kolf_3 cells ----------
feng_share = {}
for f in sorted((C.V / "raw/feng").glob("*Cell-Metadata.tsv.gz")):
    cnt = collections.Counter()
    with gzip.open(f, "rt") as fh:
        r = csv.reader(fh, delimiter="\t"); h = next(r); j = h.index("Cell_Line")
        for row in r:
            cnt[row[j]] += 1
    tot = sum(cnt.values())
    feng_share[f.name.split("_Cell")[0]] = dict(cells=tot, lines=len(cnt), kolf_2=cnt.get("kolf_2", 0), kolf_3=cnt.get("kolf_3", 0),
                                                 kolf_frac=round((cnt.get("kolf_2", 0) + cnt.get("kolf_3", 0)) / tot, 4))
checks["feng_kolf_cell_share"] = feng_share

# ---------- destination identity: truth dataset + what counts as same line / same dataset ----------
DEST = {
 "H1": dict(line="H1 hESC", truth_dataset="VCC 2025 H1 training data (Arc Institute, 10x Flex), thinned",
            same_line_or_dataset={"H1": "same dataset (vcc2025_target_batch)", "H1GW": "same line, other lab (GSE295214, Weissman)"}),
 "OR_HCT116X": dict(line="HCT116", truth_dataset="Orion / X-Atlas HCT116 genome-scale CRISPRi (orion_target_cellline_all)",
            same_line_or_dataset={"HCT116": "same dataset", "HEK293T": "same underlying Orion dataset (other line)"}),
 "KOLF": dict(line="KOLF2.1J iPSC", truth_dataset="KOLF2.1J genome-scale CRISPRi (kolf_target store), disjoint-cell panel",
            same_line_or_dataset={"KOLF": "same dataset",
                                  "FENG": "pooled HipSci iPSC lines incl. kolf_2 (HPSI0114i-kolf_2, the parental line of KOLF2.1J) and kolf_3 (same donor)"}),
 "RN_K562": dict(line="K562", truth_dataset="Replogle 2022 K562 GWPS (raw/scperturb/ReplogleWeissman2022_K562_gwps.h5ad), 800 ctrl / <=200 cells",
            same_line_or_dataset={"K562": "same dataset (replogle_nadig/K562 = GWPS, 9,729 targets)", "VIP": "same line, other lab (VIPerturb-seq K562)"}),
 "RN_HepG2": dict(line="HepG2", truth_dataset="Nadig 2025 HepG2 essential-scale CRISPRi",
            same_line_or_dataset={"HEPG2": "same dataset"}),
 "RN_Jurkat": dict(line="Jurkat", truth_dataset="Nadig 2025 Jurkat essential-scale CRISPRi",
            same_line_or_dataset={"RNJURKAT": "same dataset", "JURKAT": "same line, other lab (GSE249595, high-MOI)"}),
 "RN_RPE1": dict(line="hTERT-RPE1", truth_dataset="Replogle 2022 RPE1 essential-scale CRISPRi",
            same_line_or_dataset={"RPE1": "same dataset"}),
}

# ---------- per-panel stats, pack check, exclusions, leak verdict ----------
gcap = C.g_cap(recompute=not (OUT / 'g_cap.json').exists()); gcap_info = json.loads((OUT / "g_cap.json").read_text())
HP = C.harness_panels(); panels = {}
for p in C.PANELS7:
    lab = C.labels(p); vc = collections.Counter(lab); nctrl = vc.pop("non-targeting", 0)
    a = C.open_panel(p); ng = a.n_vars; a.file.close()
    inp = C.pack(p, "inputs"); fams_pack = [str(x) for x in inp["families"]]; tg_pack = [str(x) for x in inp["targets"]]
    checks[f"S1_pack_targets_equal_truth_{p}"] = bool(set(tg_pack) == set(vc)) and len(tg_pack) == len(HP[p]["targets"])
    ex_m1 = sorted(EXCL[p]); ex_fusion = sorted(set(FAMS_H) - set(fams_pack))
    s4key = S4KEY[p]; s4_line = S4_LINE[p]
    s4_ex = sorted({s4_line.upper(), *[c for c in CONFLICT.get(s4_line, [])], *S4_EXTRA.get(s4key, [])})
    s4_ex_effective = sorted({s4_line.upper(), *[c for c in CONFLICT.get(s4_line, []) if c in set(tb)], *S4_EXTRA.get(s4key, [])})
    excluded_all = sorted(set(ex_m1) | set(OUTER[p]))
    leaks = []
    for src, why in DEST[p]["same_line_or_dataset"].items():
        in_fusion = src in fams_pack
        in_m1 = (src not in EXCL[p]) and W_NEW_H.get(src, 0) > 0 or (src not in EXCL[p] and src in FAMS_H)
        s4_has = src.upper() in set(tb) and src.upper() not in s4_ex_effective
        if in_fusion or s4_has:
            leaks.append(dict(source=src, relation=why, in_fusion_part=in_fusion, in_first_part_and_neighbour_library=src in FAMS_H and src not in EXCL[p],
                              in_S4_training=s4_has))
    # special: FENG carries kolf_2/kolf_3 in S4 table_v2 (bio FENG) and in the pooled FENG sources
    panels[p] = dict(
        truth_path=str(C.truth_path(p).relative_to(C.ROOT)), full_path=(str(C.truth_path(p, "full").relative_to(C.ROOT))),
        plan_path_note=C.PLAN_PATH_NOTES.get(p), destination_line=DEST[p]["line"], truth_dataset=DEST[p]["truth_dataset"],
        n_targets=len(vc), cells_per_target_median=float(np.median(list(vc.values()))), n_control_cells=int(nctrl),
        gene_count=int(ng), g_cap_size=int(gcap_info["n"]), genes_on_official_axis_with_ctrl_umi=gcap_info["per_panel"][p],
        harness_n_cells_param=HP[p]["n_cells"],
        packs={k: f"vcc2026_data/reports/era_2026-10-07/packs/{k}_{p}.npz" for k in C.PACK_KINDS},
        pack_families=fams_pack,
        excluded_sources=dict(all=excluded_all, M1_first_part_and_neighbour_library=ex_m1,
                              M1_fusion_part_zero_weight=ex_fusion, S4_part_key=s4key,
                              S4_excluded_biologies_intended=s4_ex, S4_excluded_biologies_effective=s4_ex_effective,
                              M2_rho_hat="leave-this-line-out per external/era-aivc/aivc/tasks/vcc_sum/task.py docstring (not re-derived)",
                              notes="KOLF family weight is 0 for every somatic destination (ml_h1fuse.pieces); H1GW/VIP/JURKAT have zero harness weight"),
        same_line_or_dataset_sources=DEST[p]["same_line_or_dataset"],
        leak=("LEAK" if leaks else "OK"), leak_detail=leaks)

kolf = panels["KOLF"]
if kolf["leak"] == "LEAK":
    kolf["leak_summary"] = ("FENG (pooled HipSci iPSC; FENG_T/GWF/GWN) contains kolf_2 = HPSI0114i-kolf_2, the parental line of "
                            "KOLF2.1J, and kolf_3 (same donor): {:.1%}/{:.1%}/{:.1%} of FENG targeted/GW-fitness/GW-nonfitness cells. "
                            "FENG is NOT excluded for the KOLF destination in the M1 fusion part (weight 0.25 of 15.75 = 1.6 %), the "
                            "equal-weight first part (1 of 9 families) or the neighbour library. In the S4 refit the per-line "
                            "fengline/kolf_2 and kolf_3 contexts are filed under biology KOLF in table_v2 and ARE excluded, but the pooled "
                            "Feng rows (xl/feng, fengpooled/fitness; biology FENG) still contain kolf_2/kolf_3 cells "
                            "(ext_baselines.CONFLICT_GENE['KOLF'] = v1 labels FENGLKOLF2/3, absent from table_v2). Magnitude: kolf cells are "
                            "~0.2-0.3 % of the fusion weight and ~1.3-2 % of the first part -- small, but a same-line source by definition S3.").format(
        feng_share["TargetedScreen"]["kolf_frac"], feng_share["GenomeWideScreen_FitnessGenes"]["kolf_frac"],
        feng_share["GenomeWideScreen_NonFitnessGenes"]["kolf_frac"])
extra = {"OR_HCT116 (8th pack, not a paper panel)": dict(excluded=sorted(EXCL["OR_HCT116"]),
         note="keeps HEK293T from the SAME Orion dataset as the truth (same-dataset, other line); OR_HCT116X exists to remove it")}
(OUT / "panels.json").write_text(json.dumps(dict(
    frozen="2026-11-05", panels=panels, g_cap=dict(n=int(gcap_info["n"]), definition=gcap_info["definition"], file="results/g_cap.json"),
    other_packs=extra, feng_kolf_cell_share=feng_share), indent=1))

# ---------- families.json: the M1 registry as tv2_model.FAMILY + frozen F2 weights ----------
def src_meta(s):
    for d in ("transfer_v2_2026-09-23/sources", "transfer_v2_2026-09-23/sources_new"):
        f = RP / d / f"{s}.json"
        if f.exists():
            m = json.loads(f.read_text()); npz = f.with_suffix(".npz")
            return m, str(npz.relative_to(C.ROOT)) if npz.exists() else None
    return {}, None

ANN = {  # family -> (cell line/type, lab/dataset, scale, perturbation, counts_for_C_public)
 "K562": ("K562 (CML)", "Replogle et al. 2022, Weissman lab: K562 GWPS", "genome-wide", "CRISPRi", True),
 "HCT116": ("HCT116 (colorectal)", "Xaira X-Atlas/Orion 2025", "genome-wide", "CRISPRi", True),
 "HEK293T": ("HEK293T (embryonic kidney)", "Xaira X-Atlas/Orion 2025", "genome-wide", "CRISPRi", True),
 "CD4": ("primary CD4+ T cells (3 states pooled)", "genome-scale CD4 T-cell CRISPRi (citation resolved in Task 2)", "genome-wide", "CRISPRi", True),
 "KOLF": ("KOLF2.1J iPSC", "KOLF2.1J CRISPRi Perturb-seq (citation resolved in Task 2)", "genome-wide", "CRISPRi", True),
 "FENG": ("HipSci iPSC lines, pooled (24-34 donors incl. kolf_2/kolf_3)", "Feng et al., HipSci CRISPRi Perturb-seq (citation resolved in Task 2)",
          "genome-wide (fitness + non-fitness arms; 1-2 cells per line-target) + targeted 445", "CRISPRi", True),
 "H1": ("H1 hESC", "VCC 2025 training data (Arc Institute, 10x Flex)", "other (300 targets)", "CRISPRi", False),
 "H1GW": ("H1 hESC", "GSE295214 (Weissman lab, 10x 5')", "genome-scale (2,747 targets passing QC)", "CRISPRi", True),
 "HEPG2": ("HepG2 (hepatocellular)", "Nadig et al. 2025, Weissman lab", "essential", "CRISPRi", True),
 "RNJURKAT": ("Jurkat (T-ALL)", "Nadig et al. 2025, Weissman lab", "essential", "CRISPRi", True),
 "RPE1": ("hTERT-RPE1", "Replogle et al. 2022, Weissman lab", "essential", "CRISPRi", True),
 "DLD1": ("DLD-1 (colorectal)", "GSE337988", "large (4,714 targets; above essential-scale)", "CRISPRi", True),
 "JURKAT": ("Jurkat (T-ALL)", "GSE249595 (high-MOI; deliberately weight 0)", "other (2,733 targets, high-MOI)", "CRISPR", False),
 "VIP": ("K562 (CML)", "VIPerturb-seq K562, 10x Flex (Zenodo 18460279); rejected as a source", "genome-wide", "CRISPRi", False),
}
fam_rows = collections.defaultdict(list)
for s, fam in FAMILY.items():
    m, npz = src_meta(s); fam_rows[fam].append((s, m, npz))
families = []
for fam, srcs in fam_rows.items():
    a = ANN.get(fam); m0 = srcs[0][1]
    if a is None:   # scPerturb small contexts (sources_new)
        a = (str(m0.get("celltype")), f"scPerturb: {', '.join(s for s, _, _ in srcs)}",
             f"other ({'/'.join(str(m.get('targets')) for _, m, _ in srcs)} targets)", str(m0.get("ptype")), False)
    w_f2 = F2W.get(fam, 0.0); w_first = BASE_NEW.get(fam, None)
    in_registry_active = w_f2 > 0 or fam == "H1GW"
    families.append(dict(
        family=fam, cell_line_or_type=a[0], lab_dataset=a[1], scale=a[2], perturbation=a[3],
        sources=[dict(source=s, file=npz, store=m.get("store") or m.get("source") or m.get("geo"),
                      n_targets=m.get("n_all", m.get("targets")), n_targets_built_for_final_panel=m.get("n_needed")) for s, m, npz in srcs],
        n_targets=max([x for x in [m.get("n_all", m.get("targets")) for _, m, _ in srcs] if isinstance(x, (int, float))] or [None]) if True else None,
        f2_fusion_weight=w_f2, harness_fusion_weight=W_NEW_H.get(fam, 0.0),
        role=("active in frozen F2" if w_f2 > 0 else ("same-line/stem fallback (tv2_final --same-line; weight 0 for somatic D/E/F)" if fam == "H1GW"
              else "registered, weight 0")),
        counts_for_C_public=bool(a[4] and in_registry_active)))
# C_public: distinct lines among counting families
LINE_KEY = {"K562": "K562", "HCT116": "HCT116", "HEK293T": "HEK293T", "CD4": "primary CD4+ T", "KOLF": "KOLF2.1J iPSC",
            "FENG": "HipSci iPSC pool", "H1GW": "H1 hESC", "H1": "H1 hESC", "HEPG2": "HepG2", "RNJURKAT": "Jurkat", "JURKAT": "Jurkat",
            "RPE1": "hTERT-RPE1", "DLD1": "DLD-1", "VIP": "K562"}
lines = collections.defaultdict(list)
for f in families:
    if f["counts_for_C_public"]:
        lines[LINE_KEY[f["family"]]].append(f["family"])
for f in families:
    if f["family"] in LINE_KEY and LINE_KEY[f["family"]] in lines and f["family"] not in lines[LINE_KEY[f["family"]]]:
        lines[LINE_KEY[f["family"]]].append(f["family"] + " (not counted: weight 0 or sub-scale)")
c_pub = len(lines)
sens = dict(main=c_pub,
            without_DLD1_if_4714_targets_not_counted_as_essential_scale=c_pub - 1,
            merging_HipSci_pool_with_KOLF2_1J_iPSC=c_pub - 1,
            both=c_pub - 2,
            counting_each_HipSci_donor_line_separately="not applicable: per-line Feng broad screens have 1-2 cells per line-target (unusable per line); targeted screen is 445 targets")
(OUT / "families.json").write_text(json.dumps(dict(
    frozen="2026-11-05", registry_source="vcc2026_data/scripts/tv2_model.py FAMILY + F2 spec weights (reports/final_2026-10-27/F1/F2_spec.json parts[1].cfg.weights)",
    definition_C_public="plan op. def. 5: distinct cell lines/types with a public genome-scale or essential-scale CRISPRi/KO Perturb-seq response in the M1 registry; same line from different labs counts once",
    C_public=c_pub, C_public_sensitivity=sens, line_to_families=dict(lines), families=families), indent=1))

# ---------- official record ----------
rows = {}
def add(eid, **kw):
    r = rows.setdefault(eid, dict(entry_id=eid)); [r.setdefault(k, v) for k, v in kw.items() if v is not None and v != ""]
    r.setdefault("sources", set()).add(kw.get("_src"))
import pandas as pd
o = pd.read_csv(RP / "official_scores.csv")
for _, r in o.iterrows():
    add(r.entry_id, label=r.submission, note=r.note, rank=r["rank"], score_avg=r.score_avg, score_pds=r.score_pds, score_mse=r.score_mse,
        score_nmae=r.score_nmae, score_fid=r.score_fid, score_reach=r.score_reach, score_jac=r.score_jac, phase="validation", _src="official_scores.csv")
for f in [RP / "simple_residual_2026-09-16/submission/official_scores_new.csv", RP / "context_forcing_strength_2026-09-16/submission/official_scores_new.csv"]:
    for _, r in pd.read_csv(f).iterrows():
        add(r.entry_id, label=r.model_name, model_name=r.model_name, rank=r["rank"], score_avg=r.score_avg, score_pds=r.score_pds,
            score_mse=r.score_mse, score_nmae=r.score_nmae, score_fid=r.score_fid, score_reach=r.score_reach, score_jac=r.score_jac,
            submission_date=r.recorded_utc, phase="validation", _src=str(f.relative_to(RP)))
for f in sorted(RP.glob("submission_2026-09-2*/*status.json")):
    j = json.loads(f.read_text()); lab = f.name.replace(".status.json", "").replace("status.json", "")
    lab = lab or {"submission_2026-09-21": "S3", "submission_2026-09-22": "S4"}.get(f.parent.name, "")
    add(j["entry_id"], label=lab, s_label=lab, model_name=j.get("model_name"), rank=j.get("rank"), score_avg=j.get("score_avg"), score_pds=j.get("score_pds"),
        score_mse=j.get("score_mse"), score_nmae=j.get("score_nmae"), score_fid=j.get("score_fid"), score_reach=j.get("score_reach"),
        score_jac=j.get("score_jac"), submission_date=j.get("submission_date"),
        phase=("final" if j.get("partition") == "test" else "validation"), _src=str(f.relative_to(RP)))
for f in [RP / "submission_2026-09-20/S1_status.log", RP / "submission_2026-09-20/S2_status.log"]:
    t = f.read_text(); g = lambda k: (re.search(rf"{k}\s*:?\s*([^\n]+)", t) or [None, None])[1]
    add(g("entry").strip(), label=f.name[:2], s_label=f.name[:2], model_name=g("model").strip(), rank=float(g("Rank")), score_avg=float(g("Overall")),
        submission_date=g("submitted").strip(), phase="validation", _src=str(f.relative_to(RP)))
for lab, sc, rk in [("B10", 0.0868, 328), ("B11", 0.1156, 234), ("B14", 0.0888, 345)]:   # PROJECT_HISTORY.md "Official submissions, all 14": no entry id on disk
    add(f"NO_ID_{lab}", label=lab, rank=rk, score_avg=sc, phase="validation", note="entry id not recorded on disk; score/rank from PROJECT_HISTORY.md", _src="PROJECT_HISTORY.md")
df = pd.DataFrame([{**{k: v for k, v in r.items() if k != "sources"}, "record_sources": ";".join(sorted(s for s in r["sources"] if s))} for r in rows.values()])
df["submission_date"] = df.get("submission_date")
df = df.sort_values(["phase", "submission_date", "label"], na_position="first", key=lambda s: s.astype(str)).reset_index(drop=True)
cols = ["entry_id", "label", "s_label", "model_name", "phase", "submission_date", "rank", "score_avg", "score_pds", "score_mse", "score_nmae",
        "score_fid", "score_reach", "score_jac", "note", "record_sources"]
df[[c for c in cols if c in df]].to_csv(OUT / "official_record.csv", index=False)
lb = json.loads((RP / "leaderboard_2026-10-27.json").read_text())["entries"]
ours = [e for e in lb if e.get("id") == "etQ6SGWT1uu0wk24xuvQ"]
checks["official_count"] = dict(rows=len(df), with_entry_id=int((~df.entry_id.str.startswith("NO_ID")).sum()),
                                validation=int((df.phase == "validation").sum()), final=int((df.phase == "final").sum()),
                                leaderboard_2026_10_27_submission_count=(ours[0]["submission_count"] if ours else None),
                                leaderboard_2026_10_27_team=(ours[0]["team_name"] if ours else None))
checks["g_cap"] = dict(n=int(gcap_info["n"]), per_panel=gcap_info["per_panel"])
(OUT / "task04_checks.json").write_text(json.dumps(checks, indent=1, default=str))
print(json.dumps(dict(leaks={p: v["leak"] for p, v in panels.items()}, C_public=c_pub, lines=dict(lines), g_cap=int(gcap_info["n"]),
                      official=checks["official_count"]), indent=1))
