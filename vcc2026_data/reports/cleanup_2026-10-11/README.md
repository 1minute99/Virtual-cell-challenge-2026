# Disk cleanup 2026-10-11 (user: "꼭 필요한거 아니면 지워")

`delete_list.txt` (1,332 paths, ~247 GB) = run_cleanup.sh deletes them.
- 41 superseded `.vcc` packages (receipts/specs/sha256 stay; S16.vcc, the best official, is kept)
- leftover `pred.h5ad` of S17 / S17d / rehearsals / profiling
- harness and DE-variant delta CACHES not read for > 5 days (`final_harness_2026-09-24/cache`, `denoise_src/cache_*`);
  they are recomputed on demand
- `src_ablation/` (flagged for deletion since 2026-10-01), `denoise_src/full`, `dnALL`, `rsel`, `disc_residual`
- >200 MB intermediates of finished studies not read for > 5 days, `~/.cache/pip`
Protected: raw data, harness real panels, S4 tables (table_v1/v2, ext_effects), sources / sources_new,
replogle_nadig + mock (tv2_sources / fh_prep), S1_refit.delta.npy + S4_ctx (tv2_s4_rebuild), s4_improve, h1fuse,
denoise_src official_* / dnK100r0.25 (S16 DE side links), ERA packs, model checkpoints, all FINDINGS / scores / receipts.
