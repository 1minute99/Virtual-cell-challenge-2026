# Source acquisition — manual downloads needed

These A+B sources could not be auto-downloaded (paywalled, bot-blocked, or no
Open-Access copy). Download each **original** PDF via your university VPN /
library proxy, save it **into `/home/wonmin-kim/research-library/pdf/`** under
the **exact** "Save as" filename below, then **re-run acquire-sources** to
reconcile (resolved rows disappear) — or run `ingest-source` on the ones
already present.

Do NOT substitute a preprint, prior version, or book review for the original —
`ingest-source` hard-stops on a missing original.

**2026-11-06 — user decision, option (b)** (see `knowledge/_meta/log.md`): where an open-access
preprint of the SAME work exists it has been acquired under its OWN bibkey (`<library>/pdf/<preprint-key>.pdf`,
BibTeX note "preprint version of <published key>"). Preprint pages quote the preprint only. The published
originals below are still wanted (rows stay open) — the last column records the preprint status.

Last updated: 2026-11-06 · 19 of 70 A+B originals still needed (A: 10, B: 9); maas-2005 original acquired 2026-11-06 (row kept, ticked). 12 of the 19 have an OA preprint on disk under its own key.

| ☐ | Citekey | Author Year | Short title | Grade | DOI / landing | Best download URL(s) | Save as (exact filename) | Why manual | Preprint status (option b, 2026-11-06) |
|---|---------|-------------|-------------|-------|---------------|----------------------|--------------------------|------------|---|
| [ ] | adduri-2026-state | Adduri et al. 2026 | state | A | [10.1016/j.cell.2026.07.052](https://doi.org/10.1016/j.cell.2026.07.052) | — | `adduri-2026-state.pdf` | html_login_page: bot-blocked/paywalled host: cell-landing | preprint acquired: `adduri-2025-state` (bioRxiv 10.1101/2025.06.26.661135 v2, 2025-07-10) |
| [ ] | adduri-2026-virtual-cell-challenge | Adduri et al. 2026 | virtual-cell-challenge | A | [10.1016/j.cell.2026.08.004](https://doi.org/10.1016/j.cell.2026.08.004) | — | `adduri-2026-virtual-cell-challenge.pdf` | html_login_page: bot-blocked/paywalled host: cell-landing | no preprint found (bioRxiv published-DOI map, Crossref posted-content, Europe PMC preprints, arXiv, Zenodo, web). The Cell original is marked "Free article" on PubMed (PMID 42648290; © The Authors, CC BY-NC-ND per Europe PMC): full text at cell.com/cell/article/S0092-8674(26)00931-1/fulltext; bot gets 403 → browser download |
| [ ] | breiman-1983-how-many-variables | Breiman & Freedman 1983 | how-many-variables | A | [10.1080/01621459.1983.10477941](https://doi.org/10.1080/01621459.1983.10477941) | — | `breiman-1983-how-many-variables.pdf` | no_oa_url: Taylor & Francis (JASA) paywall; no OA copy | no preprint found online (UCB Statistics Tech. Rep. No. 1, Nov 1981 / rev. Feb 1982, listed in the Berkeley index but not digitised) |
| [ ] | hastie-2022-ridgeless-interpolation | Hastie et al. 2022 | ridgeless-interpolation | A | [10.1214/21-AOS2133](https://doi.org/10.1214/21-AOS2133) | https://pmc.ncbi.nlm.nih.gov/articles/PMC9481183/ (accepted manuscript) | `hastie-2022-ridgeless-interpolation.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `hastie-2019-ridgeless-interpolation` (arXiv 1903.08560 v5, 2020-12-07) |
| [ ] | li-2026-systematic-comparison | Li et al. 2026 | systematic-comparison | A | [10.1126/sciadv.aed3414](https://doi.org/10.1126/sciadv.aed3414) | https://pmc.ncbi.nlm.nih.gov/articles/PMC13557097/ (OA, published version)<br>https://www.science.org/doi/pdf/10.1126/sciadv.aed3414 | `li-2026-systematic-comparison.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `li-2024-systematic-comparison` (bioRxiv 10.1101/2024.12.23.630036 v3, 2026-08-03) |
| [ ] | nadig-2025-transcriptome-wide | Nadig et al. 2025 | transcriptome-wide | A | [10.1038/s41588-025-02169-3](https://doi.org/10.1038/s41588-025-02169-3) | https://pmc.ncbi.nlm.nih.gov/articles/PMC13063516/ (author manuscript) | `nadig-2025-transcriptome-wide.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `nadig-2024-transcriptome-wide` (bioRxiv 10.1101/2024.07.03.601903 v1, 2024-07-03; different title) |
| [ ] | replogle-2022-genome-scale | Replogle et al. 2022 | genome-scale | A | [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013) | https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/ (accepted manuscript) | `replogle-2022-genome-scale.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `replogle-2021-genome-scale` (bioRxiv 10.1101/2021.12.16.473013 v3, 2022-05-28) |
| [ ] | rosset-2020-random-x | Rosset & Tibshirani 2020 | random-x | A | [10.1080/01621459.2018.1424632](https://doi.org/10.1080/01621459.2018.1424632) | preprint only: https://arxiv.org/abs/1704.08160 | `rosset-2020-random-x.pdf` | no_oa_url: Taylor & Francis (JASA) paywall; only arXiv preprint (not the original) | preprint acquired: `rosset-2017-random-x` (arXiv 1704.08160 v2, 2017-06-10) |
| [ ] | wei-2026-benchmarking-algorithms | Wei et al. 2026 | benchmarking-algorithms | A | [10.1038/s41592-025-02980-0](https://doi.org/10.1038/s41592-025-02980-0) | — | `wei-2026-benchmarking-algorithms.pdf` | no_oa_url: Nature Methods, closed | no preprint found (bioRxiv map, Crossref posted-content, Europe PMC, web search; only code/data on GitHub/Zenodo/Figshare) |
| [ ] | zhang-2026-tahoe | Zhang et al. 2026 | tahoe | A | [10.1016/j.cell.2026.08.035](https://doi.org/10.1016/j.cell.2026.08.035) | — | `zhang-2026-tahoe.pdf` | html_login_page: bot-blocked/paywalled host: cell-landing | preprint acquired: `zhang-2025-tahoe` (bioRxiv 10.1101/2025.02.20.639398 v3, 2025-05-10; different title) |
| [ ] | christiansen-2022-causal-framework | Christiansen et al. 2022 | causal-framework | B | [10.1109/TPAMI.2021.3094760](https://doi.org/10.1109/TPAMI.2021.3094760) | preprint only: https://arxiv.org/abs/2006.07433 | `christiansen-2022-causal-framework.pdf` | no_oa_url: IEEE TPAMI paywall; only arXiv preprint | preprint acquired: `christiansen-2020-causal-framework` (arXiv 2006.07433 v3, 2021-08-17) |
| [ ] | copas-1983-regression-shrinkage | Copas 1983 | regression-shrinkage | B | [10.1111/j.2517-6161.1983.tb01258.x](https://doi.org/10.1111/j.2517-6161.1983.tb01258.x) | — | `copas-1983-regression-shrinkage.pdf` | no_oa_url: JRSS-B 1983, closed | no preprint found |
| [ ] | cui-2024-scgpt | Cui et al. 2024 | scgpt | B | [10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0) | preprint only: https://doi.org/10.1101/2023.04.30.538439 | `cui-2024-scgpt.pdf` | no_oa_url: Nature Methods, closed; only bioRxiv preprint | preprint acquired: `cui-2023-scgpt` (bioRxiv 10.1101/2023.04.30.538439 v2, 2023-07-02) |
| [ ] | dicker-2013-equivariant-prediction | Dicker 2013 | equivariant-prediction | B | [10.1214/13-EJS826](https://doi.org/10.1214/13-EJS826) | https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-7/issue-none/Optimal-equivariant-prediction-for-high-dimensional-linear-models-with-arbitrary/10.1214/13-EJS826.full | `dicker-2013-equivariant-prediction.pdf` | bot_block: Diamond OA (EJS) but Project Euclid serves an Incapsula bot challenge — a normal browser download works | no preprint found (arXiv 1203.4572 / 1102.2952 by the same author have different titles; identity not established, not used). The EJS original is Diamond OA: browser download |
| [ ] | feng-2026-crispri-map | Feng et al. 2026 | crispri-map | B | [10.1016/j.xgen.2025.101076](https://doi.org/10.1016/j.xgen.2025.101076) | https://pmc.ncbi.nlm.nih.gov/articles/PMC12903452/ (OA, published version) | `feng-2026-crispri-map.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `feng-2024-crispri-map` (bioRxiv 10.1101/2024.11.28.625833 v1, 2024-11-28) |
| [x] | maas-2005-sample-sizes | Maas & Hox 2005 | sample-sizes | B | [10.1027/1614-2241.1.3.86](https://doi.org/10.1027/1614-2241.1.3.86) | — | `maas-2005-sample-sizes.pdf` | no_oa_url: Hogrefe (Methodology), closed | ORIGINAL acquired (not a preprint): publisher version from the Utrecht University Repository (hdl 1874/23635) saved as `maas-2005-sample-sizes.pdf` |
| [ ] | mcneish-2016-small-sample | McNeish & Stapleton 2016 | small-sample | B | [10.1007/s10648-014-9287-x](https://doi.org/10.1007/s10648-014-9287-x) | — | `mcneish-2016-small-sample.pdf` | no_oa_url: Springer, closed | no preprint found |
| [ ] | nourreddine-2026-crispri-atlas | Nourreddine et al. 2026 | crispri-atlas | B | [10.1038/s41587-026-03199-w](https://doi.org/10.1038/s41587-026-03199-w) | preprint only: https://doi.org/10.1101/2024.11.03.621734 | `nourreddine-2026-crispri-atlas.pdf` | no_oa_url: Nature Biotechnology, closed; only preprint | preprint acquired: `nourreddine-2024-perturbation-atlas` (bioRxiv 10.1101/2024.11.03.621734 v2, 2024-11-04; different title) |
| [ ] | roohani-2025-virtual-cell-challenge | Roohani et al. 2025 | virtual-cell-challenge | B | [10.1016/j.cell.2025.06.008](https://doi.org/10.1016/j.cell.2025.06.008) | https://www.cell.com/action/showPdf?pii=S0092867425006750 (OA on Cell) | `roohani-2025-virtual-cell-challenge.pdf` | http_403: bot-blocked/paywalled host: cell | no preprint found (the Cell commentary is OA on cell.com: browser download) |
| [ ] | zhu-2026-cd4-perturb-seq | Zhu et al. 2026 | cd4-perturb-seq | B | [10.1016/j.cell.2026.08.002](https://doi.org/10.1016/j.cell.2026.08.002) | https://pmc.ncbi.nlm.nih.gov/articles/PMC13614707/ (author manuscript) | `zhu-2026-cd4-perturb-seq.pdf` | html_login_page|http_403: bot-blocked/paywalled host: europepmc, pmc | preprint acquired: `zhu-2025-cd4-perturb-seq` (bioRxiv 10.64898/2025.12.23.696273 v1, 2025-12-24) |

## Downloaded, but not the publisher version (optional upgrade)

These are on disk and ingestable; replace the file (same filename) if you need the final pagination.

| Citekey | On disk | Preferred original |
|---------|---------|--------------------|
| `rosenfeld-2021-risks-irm` | arXiv 2010.05761 | ICLR 2021 camera-ready: https://openreview.net/forum?id=BbNIbVPJ-42 (OpenReview returned 403 to the bot) |
| `du-2021-few-shot-representation` | arXiv 2002.09434 | ICLR 2021: https://openreview.net/forum?id=pW2Q2xLwIMD |
| `lucas-2021-meta-learning-bounds` | arXiv 2010.07140 | ICLR 2021: https://openreview.net/forum?id=SZ3wtsXfzQR |
| `piao-2026-adaptive-contexts` | arXiv 2602.18885 | ICML 2026: https://openreview.net/forum?id=aJVkQAd7tv (or PMLR when listed) |
| `fu-2026-strand` | arXiv 2602.10156 | NeurIPS 2026 proceedings (not yet published) |
| `rothenhausler-2021-anchor-regression` | ETH Research Collection accepted manuscript | JRSS-B publisher version (paywalled): https://doi.org/10.1111/rssb.12398 |
| `ji-2024-universal-prediction` | bioRxiv v3 (2025-09-23, latest) | — note: the STOP-FLAG quotes were read from v2; check fig. S4 numbering in v3 |
