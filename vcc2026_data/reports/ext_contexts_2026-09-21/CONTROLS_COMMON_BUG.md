# Gene-order bug in `reports/cellfm_context_2026-09-20/controls{,_common}/*.h5ad` — found 2026-09-21

`scripts/cellfm_controls.py::write()` does

    uniq, inv = np.unique(genes, return_inverse=True)      # SORTED symbols
    if len(uniq) != len(genes):                            # only when duplicates exist ...
        X = collapse(X, inv)                               # ... are the columns reordered
    var = pd.DataFrame(index=uniq)

so for a source whose symbols are already unique the counts stay in SOURCE order while the gene
names are alphabetically sorted. Measured on the written files (top columns by total count):

| biology | source | top "genes" | status |
|---|---|---|---|
| HCT116 | orion (duplicate symbols) | MT-CO2, RPLP1, MT-CO1, PTMA | correct |
| HEK293T | orion (duplicate symbols) | MT-CO1, MT-CO2, RPL39, RPS8 | correct |
| HEPG2 | RN cache | ZSWIM8, TRAIP, ZSCAN29, ZSWIM7 | **scrambled** |
| JURKAT | RN cache | ZSWIM7, ZSWIM8, ZSCAN30, NMT1 | **scrambled** |
| K562 | RN cache | ZSWIM7, ZSCAN30, ZSWIM6, ZW10 | **scrambled** |
| RPE1 | RN cache | ZSWIM7, PGAP2, ZSWIM6, ZW10 | **scrambled** |
| KOLF | kolf h5ad | ZSWIM2, ENSG00000271787, ZSCAN5C | **scrambled** |

Two scrambled real cell lines correlate at r = 0.035 on the 5,955-gene context axis; a correctly
built CD4 profile correlates at 0.994 with the frozen table's own CD4 control row and at −0.01
with the scrambled JURKAT.

## What consumed these files (all after 2026-09-20 ~03:00)

1. `cellfm_encode.py` / cell-FM context experiment (`reports/cellfm_context_2026-09-20`): STATE,
   UCE and the PCA-of-pooled-controls context all encoded these cells. Five of seven biologies
   carried the wrong gene for every count. The conclusion "cell FMs tie PCA" compared three
   encodings of the same scrambled input and should not be cited.
2. `dcp_prepare.py` → `contexts.npz` gene rows → every gene-side context vector of the
   drug-context pretraining experiment (`reports/drug_context_pretrain_2026-09-20`, FINDINGS §3b,
   §4–7). Drug contexts (Tahoe, GSE) were built from their own count matrices and are fine. The
   gene contexts were effectively random keys: still mutually distinct (so "correct beats wrong
   context" survives) but unrelated to the drug-context manifold, which mechanically produces
   "gene contexts sit 2.3x outside the drug manifold", "the gate shrinks every correction to 7 %"
   and the failure of drug→gene transfer. **Those gene-side transfer conclusions are not
   established by that experiment.** The drug-side context-scaling curve (§2–3) is unaffected.
3. `submit_build.py` (S2, refit + PCA interaction, entry hlhJKNFjK8UHE1dq2rsV, −0.0020): its seven
   training-biology context profiles were scrambled for five biologies while A/B/C were built
   correctly from the official controls. Its official loss says nothing about the method.
4. `dcp_submit_build.py` (S3, refit + scratch FiLM residual, entry 5BTNATrQRIsC3XvRQayZ,
   −0.0005): trained on the scrambled gene-context vectors, applied to correctly built A/B/C
   vectors. Same caveat.

Not affected: everything built from the training table's `controls` (basal_context, tiny_xattn,
xattn_context, ctx_* experiments, the refit baseline, mean transfer) and all drug-side results.

## What this experiment does instead

All context profiles come from the extended training table's own pooled control profiles on the
5,592 table genes (`ext_context_assets.py`), the same object `basal_context.assemble` produced;
the context PCA is fitted on gene contexts only. `cellfm_controls.py` is left untouched as the
historical source; a fix is one line (`X = X[:, np.argsort(genes)]` when there are no
duplicates) but the downstream experiments would have to be rerun to mean anything.
