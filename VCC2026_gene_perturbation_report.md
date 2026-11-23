# VCC 2026 유전자 perturbation 데이터 확보 전략

가장 먼저 확보할 자료는 **Replogle genome-wide, Nadig HepG2/Jurkat, X-Atlas/Orion, KOLF2.1J, 공개 primary CD4 T-cell atlas**다. 이들을 중심으로 표적 유전자 범위를 넓히고, **Feng의 다중 donor iPSC, Jiang/Mixscale의 다중 세포·자극 조건, 공개 VCC 2025**로 문맥 다양성을 보강하는 구성이 합리적이다. 동일 표적의 선행 실험은 중요한 정보지만, 효과의 크기와 방향은 세포 유형, CRISPRi 효율, 자극, 시간에 따라 달라진다. 따라서 “같은 gene이 있다”와 “챌린지에 옮길 수 있는 신뢰할 만한 반응이 있다”를 함께 평가해야 한다.

조사 기준은 **2026년 9월 10일 공개 상태**다. 비교표는 91개 연구·실험 분할·배포본·탐색 자원을, 다운로드 manifest는 946개 파일 항목을 수록한다. **91개 독립 연구 또는 946개 독립 데이터셋이라는 뜻은 아니다.** Parquet shard, 같은 세포의 재배포, protein readout, 분석 결과도 파일 목록에 포함한다. 확인되지 않은 다운로드 크기, 아직 공개되지 않은 자료, 실제 표적 목록이 추출되지 않은 자료는 명시적으로 구분했다.

## 1. 첨부 유전자 목록과 챌린지 규칙

첨부된 텍스트는 `gene_name` 헤더 아래 **18,533개 행, 18,533개 고유 유전자**를 포함한다. 공식 VCC 2026의 출력 feature 수가 바로 18,533개이며, perturbation 표적은 별도 **300개**다. 따라서 첨부 목록을 validation perturbation 목록으로 취급하면 거의 모든 전장 전사체 데이터가 “표적을 포함한다”는 잘못된 결론을 얻는다. 공식 controls 번들의 **`pert_counts.csv`**가 표적 대조에 필요한 파일이고, **`gene_names.csv`**는 출력 순서를 정하는 파일이다. 현재 첨부로는 300개 표적의 정확한 교집합을 확정할 수 없다.[^vcc]

| 구분 | 공식 2026 기준 | 데이터 확보에 주는 의미 |
|---|---|---|
| 과제 | 새 세포 문맥에서 CRISPRi 후 전사체 예측 | CRISPRi를 중심으로 확보하되 문맥 전이를 검증 |
| 데이터 생성 | 10x Flex, Ultima UG100 | 일반 droplet scRNA뿐 아니라 Flex 자료의 가치가 큼 |
| validation | 익명 문맥 A/B/C | 대조군으로 문맥 특성을 파악하되 숨겨진 정답이 있다고 가정하지 않음 |
| final | 새로운 문맥 D/E/F | 특정 validation 문맥에만 맞춘 자료 선별은 위험 |
| 출력 | 18,533개 유전자, 표적당 400개 세포, 3개 문맥 | 300 × 400 × 3 = 360,000개 예측 세포 |
| count 형식 | 유한한 비음수 정수 | Z-score, LFC, log-normalized 값을 그대로 제출하면 안 됨 |
| 2026 학습 정답 | 별도 challenge training set 없음 | 외부 실험자료가 학습의 핵심 |
| 일정 | final 자료 10월 22일, 제출 마감 11월 5일 | final 배포 시 그 번들의 표적·feature·규칙을 다시 고정 |

공식 CLI의 추가 한도에는 최대 400,000개 세포, sparse nonzero 최대 4.75 billion, 세포별 total count 최대 1,000,000이 있다. 데이터 수집 단계에서 feature를 과도하게 잘라버리거나 정규화값만 저장하면 최종 count 생성까지 연결하기 어렵다. 외부 데이터 활용은 허용되지만, 각 데이터의 배포 조건은 별도로 확인한다. final과 validation의 표적 구성이 같다고 미리 가정하지 않는다.[^vcc]

```bash
uv tool install vcc-cli
vcc login
vcc datasets list
vcc datasets download controls -d official_controls
```

로그인은 참가자 계정에서 수행한다. controls ZIP을 풀고 `pert_counts.csv`의 `target_gene` 열을 확인한다. 현재 조사에서는 인증이 필요한 official dataset API를 통해 이 300개 목록을 확보하지 못했다. 2025 표적 목록이나 18,533개 feature 목록으로 대체하지 않았다.

## 2. 우선 확보할 데이터와 비용 대비 가치

아래 GB는 **1 GB = 10⁹ bytes**다. GiB는 2³⁰ bytes이며 비교표에서 별도 계산한다. 웹사이트의 “GB” 표시는 실제로 GiB인 경우가 있어 API의 bytes를 우선했다. 압축 다운로드 크기는 실행 시 RAM이나 변환 중 임시 저장 공간과 다르다.

| 순서 | 데이터 / 추천 형태 | 검증된 다운로드 크기 | 선택 이유 | pseudobulk |
|---|---|---:|---|---|
| 1 | Replogle 2022, scPerturb 3개 raw 기반 H5AD | 11.589 GB | K562 genome-wide와 두 essential 실험 | 저자 raw_bulk 3개 합계 0.550 GB |
| 1 | Nadig HepG2 + Jurkat, scPerturb v1.4 | 2.144 GB | 추가 두 세포 문맥과 수천 개 표적 | 저자 DE/LFC는 있음; count PB는 별도 생성 |
| 1 | 공개 CD4, 우선 pseudobulk | 44.567 GB | primary 세포, 4 donors, Flex, 자극 조건 | **raw UMI sum 제공** |
| 2 | Orion, 공식 HF parquet | 126.260 GB + 0.882 MB feature map | 18,903개 설계 표적, 두 세포주 | ready-made count PB 미발견; streaming 생성 |
| 2 | KOLF pan-genome QC filtered | 189.393 GB | iPSC에서 넓은 CRISPRi 표적 범위 | 전체 count PB 미발견; browser signature는 별도 |
| 2 | Feng count/guide/metadata 9개 gzip | 7.926 GB | iPSC donor·유전적 배경 전이 | LFC by gene/by line 제공 |
| 2 | Jiang 5개 stimulus Seurat object | 20.142 GB | 6개 세포주와 5개 자극 조건 | DE 있음; stimulus bulk를 gene PB로 혼동 금지 |
| 2 | VCC 2025 공개 train/validation/test | 계정별 GCS listing 필요 | H1 CRISPRi 외부 학습자료 | raw로 생성 |
| 표적 일치 시 | PerturbME A375 singlets | 1.431 GB | 새 암 세포 문맥, RNA+protein | ElasticNet 계수는 PB가 아님 |
| 모달리티 보강 | Southard RPE1 / Hs27 singlets | 11.687 / 8.048 GB | TF activation과 fibroblast 문맥 | regression mean 1.898 / 1.720 GB |

이 순서는 자원 확보의 기본안이며 실제 300개 표적 교집합을 반영한 최종 순위는 아니다. 소규모 데이터에만 존재하는 표적이 확인되면 그 데이터는 크기와 상관없이 우선 포함 후보로 올린다. 다만 같은 표적이라는 이유만으로 CRISPRa·KO·RNAi를 CRISPRi와 동일한 라벨로 합치지는 않는다.[^reploglefiles][^scperturb][^orion][^kolf][^cd4][^feng][^jiang][^southard][^perturbme]

주요 배포 레코드에서 확인한 라이선스는 다음과 같다. 원본 데이터의 조건과 코드 repository의 라이선스를 혼동하지 않도록, 다운로드한 레코드의 표시와 버전을 함께 보관한다.

| 배포 레코드 | 표시된 라이선스 |
|---|---|
| Replogle 원본, Nadig DE, KOLF, PerturbME Figshare | CC BY 4.0 |
| scPerturb v1.4, Jiang v2.1, Southard processed Zenodo | CC BY 4.0 |
| Feng counts / effects Figshare | MIT |
| X-Atlas/Orion Figshare와 HF | CC BY-NC-SA 4.0 |
| CD4 / Arc GCS 등 | 데이터별 원문·배포 안내와 코드 라이선스를 각각 확인 |

위 표는 배포 API의 라이선스 필드를 전사한 것이다. 원자료의 기여자·원논문 인용도 함께 유지한다.[^reploglefiles][^targettable][^kolf][^perturbme][^scperturb][^jiang][^southard][^feng][^orioncard]

## 3. Replogle와 Arc Hugging Face 배포본 선택

### 3.1 가장 중요한 구분: genome-wide와 essential

Replogle 2022 원자료에는 **K562 genome-wide**, **K562 essential**, **RPE1** 실험이 있다. 후속 저자 분석 부록의 `TargetingGuides`와 실제 세포 수 양수 조건으로 확인한 표적은 각각 **9,866 / 2,057 / 2,393개**다. “Replogle”라는 이름만 보고 K562 essential만 받으면 genome-wide의 수천 개 표적을 놓친다. 시점과 QC가 다른 숫자를 섞지 않고 실험별로 기록해야 한다.[^replogle][^targettable]

| 원본 파일 | bytes | GB | 용도 |
|---|---:|---:|---|
| `K562_gwps_raw_singlecell_01.h5ad` | 65,830,941,948 | 65.831 | genome-wide raw single-cell |
| `K562_essential_raw_singlecell_01.h5ad` | 10,661,879,995 | 10.662 | essential raw single-cell |
| `rpe1_raw_singlecell_01.h5ad` | 8,700,873,216 | 8.701 | RPE1 raw single-cell |
| `K562_gwps_raw_bulk_01.h5ad` | 374,587,922 | 0.375 | 저자 raw pseudobulk |
| `K562_essential_raw_bulk_01.h5ad` | 79,766,954 | 0.080 | 저자 raw pseudobulk |
| `rpe1_raw_bulk_01.h5ad` | 95,350,546 | 0.095 | 저자 raw pseudobulk |

원본은 [저자 데이터 포털](https://gwps.wi.mit.edu/)과 [Figshare v1](https://doi.org/10.25452/figshare.plus.20029387.v1)에서 받는다. **raw와 normalized를 모두 받을 필요는 없다.** normalized 파일은 실험의 GEM-group 기반 표준화가 반영된 표현으로, raw UMI나 일반적인 log1p counts와 동치가 아니다. normalized bulk를 count sum으로 사용할 수 없다. 원본 전체 “159.71GB” 표시는 API 합계로 환산하면 약 **171.49 decimal GB / 159.71 GiB**다.[^reploglefiles]

다운로드 부담을 줄이려면 **scPerturb Zenodo 13350497, v1.4**의 세 파일을 사용한다. 처리 코드를 확인한 결과 Replogle의 `raw_singlecell`을 읽어 메타데이터를 정리하는 경로이며, Nadig도 raw 파일을 읽어 gzip H5AD로 저장한다. 이것이 큰 용량 차이를 단순히 “세포를 대폭 버린 작은 버전”으로 해석해서는 안 되는 이유다. 최종 객체의 세포 수, feature 수, `.X`, `layers`는 내려받은 뒤 다시 검사한다.[^scprocessing]

### 3.2 Arc와 기타 파생본의 정확한 역할

| 배포 위치 / 버전 | 확인한 내용 | 권장 판단 |
|---|---|---|
| `arcinstitute/Replogle-Nadig` | 조사 시 HF tree에는 `.gitattributes`만 존재 | 이 경로를 현재 데이터 다운로드 경로로 안내하면 실패 |
| `arcinstitute/Replogle-Nadig-Preprint`의 개별 raw 파일 | HepG2 5.615 GB, Jurkat 9.367 GB | 원자료 경로의 대안; 파일 단위 선택 |
| 같은 repo의 K562 essential / RPE1 | 파일명이 `normalized_singlecell`인 버전 | 다른 raw 두 파일과 동일한 count scale로 합치지 않음 |
| 같은 repo의 `replogle.h5ad` | 22.298 GB; 원격 H5에서 **643,413 × 6,642**, raw count layer 없음 | 가공 통합본; 전장 출력 학습의 기본 원자료로 부적합 |
| `arcinstitute/State-Replogle-Filtered/replogle_concat.h5ad` | 30.011 GB; **643,413 × 6,546**; `X_state`, `X_vci` 등 포함 | STATE 재현·embedding 분석용 보조본 |
| `arcinstitute/Replogle-Nadig-Models` | 모델·평가·예측 관련 배포 | 추가 실험 세포 데이터로 세지 않음 |
| `tahoebio/replogle-nadig-de-rhaister` | `all_delta.parquet` 6.503 GB, `all_pdex.parquet` 29.762 GB | 파생 효과/집계 표현; raw single-cell 또는 작은 PB로 오해 금지 |
| GEARS / scPertEval 등 benchmark 변형 | 보통 선택한 표적·HVG·분할·정규화가 다름 | 원문과 실험 계보를 연결하고 중복 제거 |

`State-Replogle-Filtered`의 `.X`는 확인한 각 세포주 샘플에서 비정수의 가공된 발현값이었다. `layers`는 비어 있으며 전체 feature는 6,546개다. 첨부 18,533개 symbol과 **이름을 그대로 비교한 교집합은 6,197개**, 즉 **12,336개는 그대로 대응하지 않는다**. 이름 변경을 정리하면 일부는 회복되지만 6,546개 feature만 가진 객체로 18,533개 전체를 관측 학습할 수는 없다. 더 큰 파일이 반드시 더 풍부한 raw RNA라는 보장은 없으며, embedding과 dense 저장 방식 때문에 커질 수도 있다.[^arcpre][^arcfiltered]

Preprint 통합본과 State filtered는 동일한 연구 계열·세포 수를 가진 가공본이다. 개별 원본, scPerturb, 통합본을 서로 독립된 생물학적 샘플로 취급하지 않는다. K562의 별도 sequencing-platform 표현도 같은 생물학적 라이브러리의 재측정일 수 있으므로 플랫폼 전이에는 쓰되 독립 실험으로 중복 집계하지 않는다.

### 3.3 Replogle pseudobulk 사용법

표적 평균 반응을 빠르게 시험하려면 raw_bulk 3개 합계 **549,705,422 bytes**가 좋은 출발점이다. 다만 이 조사는 원본 bulk H5의 내부 집계식을 직접 감사하지 못했다. 파일명의 “raw_bulk”만으로 **sum counts**, **cell mean**, **guide별 평균 후 gene 평균** 중 어느 것인지 단정하지 않는다. DESeq2용 biological-replicate count sum이 필요하면 single-cell raw counts와 batch 정보를 이용해 새로 만든다. 평균 반응 학습용으로 저자 PB를 활용하더라도 controls, gene order, 표준화식을 먼저 확인한다.

## 4. 넓은 유전자 범위를 담당할 대형 자료

### 4.1 Xaira X-Atlas/Orion

공식 배포는 [Hugging Face](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion)와 [Figshare 29190726 v3](https://api.figshare.com/v2/articles/29190726)다. 설계는 **18,903 protein-coding targets**, HCT116·HEK293T 두 세포 문맥, 약 8 million 실험 세포다. 공개 QC 파일의 세포 수와 전체 실험 세포 수는 다르므로 “8M를 전부 받았다”고 가정하지 않는다.[^orion]

| 형식 | 파일 / 범위 | 정확한 크기 |
|---|---|---:|
| HF parquet | HCT116 109 shards | 46,576,484,789 bytes |
| HF parquet | HEK293T 223 shards | 79,683,336,248 bytes |
| feature map | `metadata/gene_metadata.parquet` | 881,915 bytes |
| 원본 H5AD | HCT116 | 209,354,246,272 bytes |
| 원본 H5AD | HEK293T | 350,164,035,901 bytes |
| guide library | `guide_library.csv` | 3,648,784 bytes |
| guide calls | HCT116 / HEK293T | 69,606,079 / 90,886,568 bytes |

**권장: HF parquet와 gene metadata.** 원본 두 H5AD 합계 559.52 GB보다 전송 부담이 훨씬 낮다. HF commit은 manifest에 고정했다. `gene_expression`과 `gene_token_id`는 길이가 같은 sparse 배열이며, token을 gene metadata의 Ensembl/symbol로 해석한 뒤 CSR을 구성한다. token 숫자를 챌린지 feature 순서의 column index로 바로 사용하면 안 된다.

`gene_target`는 유전자 표적, `guide_target`은 guide identity다. `pass_guide_filter`의 같은 표적을 겨냥한 두 sgRNA는 두 유전자 복합 perturbation과 다르다. 정확한 target coverage는 guide library의 설계 목록과 실제 `pass_guide_filter`를 통과한 세포를 구분해서 계산한다. 원자료 문서의 control 표기를 확인한 뒤 적절한 QC를 적용하고, control까지 모두 탈락시키는 필터를 피한다. `sample`/GEM과 세포주를 보존한다.

바로 쓰는 count pseudobulk는 확인하지 못했다. parquet batch를 순차 처리하면서 `(cell_line, sample, gene_target, guide_target)`별 sum counts와 n_cells를 저장하면 single-cell 전체를 RAM에 올리지 않고 집계할 수 있다. 단순 streaming filter는 서버가 필요한 표적의 bytes만 보내준다는 뜻이 아니다. 모든 shard의 대상 column을 읽어야 할 수 있다. 배포 라이선스는 **CC BY-NC-SA 4.0**으로 표기되어 있어 상업적 재사용과 공유 조건을 모델 배포 계획과 함께 확인해야 한다.[^orioncard]

### 4.2 KOLF2.1J iPSC atlas

[저자 portal](https://y-doctor.github.io/KOLF2.1J_Perturbation_Cell_Atlas/)과 [Figshare 27261219](https://api.figshare.com/v2/articles/27261219)가 기준이다. 2.5M 이상 세포의 pan-genome CRISPRi 자료이며 2026년 정식 논문과 기존 preprint/브라우저 사이에 표적 수 표기가 다르다. 11,739·11,692 같은 설계/버전 수치를 관측 표적의 정확한 수로 바꾸어 사용하지 않는다. 현재 browser에서 추출한 perturbation 목록은 **10,167개 symbol**이며, 이것도 전체 QC H5AD의 표적을 모두 직접 센 결과는 아니다.[^kolf]

| 파일 | GB | 선택 |
|---|---:|---|
| `KOLF_Pan_Genome_QC_Filtered.h5ad` | 189.393 | 전체 표적 coverage 우선의 기본 선택 |
| `KOLF_Strong_Perturbations.h5ad` | 46.719 | 강한 효과 subset; 초기 실험/보조학습 |
| chromatin QC / strong | 4.636 / 1.529 | 해당 실험을 구분해 선택 |
| metabolic QC / strong | 6.984 / 1.326 | 해당 실험을 구분해 선택 |

Strong perturbations만 사용하면 약한 효과나 특정 표적을 누락할 수 있다. coverage 목표라면 full QC를 먼저 확보하고 효과 강도를 학습 가중치나 분석 필터로 활용한다. `Boom5426/KOLF` HF mirror의 주요 파일 크기는 저자 배포와 일치하며 MD5 대조가 가능하다. 저자 배포를 기준으로 삼고 mirror는 전송상의 대안으로 사용한다. 양쪽을 별도 데이터로 합치지 않는다.

전체 raw-count PB 파일은 확인하지 못했다. 브라우저의 NTC 기반 signature는 raw count sum과 다르므로 효과 방향·표적 탐색에 사용한다. 원시 FASTQ에는 KOLF 라인의 Y-chromosome 관련 접근 제한이 명시되어 있으며, 공개 count matrix를 사용하는 경로로 준비하는 편이 실용적이다. 제한된 raw read를 받는 것을 학습의 필수 선행 단계로 만들 필요는 없다.[^kolfrepo]

### 4.3 공개 primary CD4 atlas

**Zhu/Dann/Marson/Pritchard의 공개 CD4 atlas**는 약 22M 세포, 4 donors, Rest·8h stimulated·48h stimulated를 포함한다. 공개 S3와 저자 데이터 사전이 있어 다운로드·집계 의미가 특히 명확하다. CRISPRi, Flex, Ultima라는 점에서 2026 챌린지와 assay 측면의 관련성이 높다. 단, 챌린지 익명 문맥이 CD4라는 뜻은 아니다.[^cd4]

```bash
aws s3 ls s3://genome-scale-tcell-perturb-seq/marson2025_data/ --no-sign-request --recursive --human-readable
aws s3 cp s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.pseudobulk_merged.h5ad data/cd4/ --no-sign-request
```

| 파일/묶음 | 크기 | 내용 |
|---|---:|---|
| `D{1..4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad` 12개 | **1,735,835,115,866 bytes = 1.736 TB** | donor·조건별 single-cell |
| `GWCD4i.pseudobulk_merged.h5ad` | 44.567 GB | raw UMI sum |
| `GWCD4i.DE_stats.h5ad` | 16.786 GB | log_fc / p / adjusted p / SE / baseMean 등 |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29.424 GB | guide별 효과 통계 |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16.866 GB | donor별 효과 통계 |
| `suppl_tables/DE_stats.suppl_table.csv` | 7.277 MB | 표적·조건별 DE 요약과 세포 수 |
| `suppl_tables/sgrna_library_metadata.suppl_table.csv` | 9.935 MB | 큐레이션된 guide → gene 매핑 |
| `suppl_tables/guide_kd_efficiency.suppl_table.csv` | 13.320 MB | guide knockdown 효율 |

원격 H5 schema를 확인한 PB는 **278,684개 행 × 18,129개 RNA features**다. 한 행은 단순히 “유전자 하나”가 아니라 **guide × donor × culture condition × 10x run** 구분을 유지한다. `n_cells`, `keep_for_DE`, `keep_effective_guides`, `keep_min_cells`, `keep_total_counts` 등 필터 정보를 함께 보존한다. count sum을 비교할 때 n_cells 차이와 replicate 구조를 고려해야 한다.[^cd4readme]

원자료의 초기 guide label과 후속 큐레이션의 유전자 매핑이 달라질 수 있으므로 **큐레이션된 `perturbed_gene_name`**을 기준으로 한다. gene 이름을 guide 문자열에서 임의 분해하지 않는다. Rest와 stimulation을 서로의 대조군으로 사용해서는 안 되며 같은 donor·condition·run에 대응하는 NTC를 사용한다. PB만으로도 donor별 평균 효과 학습을 시작할 수 있고, 필요한 경우 해당 문맥의 single-cell만 추가한다.

다른 **Moonen/Claringbould CD4 atlas**는 4.1M 세포, 7,664 promoters·1,032 CRE의 독립 연구다. 그러나 [Zenodo 18876252](https://zenodo.org/records/18876252)는 **2027년 3월 5일까지 embargo**다. 2026 대회 학습용으로 지금 받을 수 있는 데이터 목록에 포함해 저장 공간이나 coverage를 계산하면 안 된다. 공개 예정 파일명에 `hvg`가 있는 8,001-feature 버전도 전장 RNA를 대체하지 못한다.[^moonen]

## 5. 문맥 전이를 보강할 자료

### 5.1 Feng: 다중 donor iPSC

Cell Genomics 논문의 genome-scale 실험은 **26 donors에서 유래한 34 iPSC lines, 7,226개 설계 표적**이다. 초기 분석에서 219,206개의 guide-assigned perturbation 세포와 499,998개의 대조/무guide 세포를 구분했으며, 표적당 10개 이상 기준을 통과한 6,673개가 LFC 분석에 사용된다. 후속 focused 실험은 444개 표적에 집중하며 19개 회수된 세포주에서 635,022개의 표적 할당 세포를 얻었다. 두 실험은 서로 다른 세포를 가지지만 일부 표적이 겹친다.[^fengpaper]

추천은 **Figshare 27989294 v2의 9개 파일, 총 7,925,882,368 bytes**다. Fitness genes / NonFitness genes / Targeted screen 각각에 RNA UMI, guide UMI, cell metadata가 있다. 압축 CSV는 열 방향과 barcode 정렬을 먼저 확인하고 chunk로 sparse matrix를 구성한다. RNA count table만 받고 donor/line 및 guide metadata를 생략하면 이 데이터의 가장 큰 장점을 잃는다.

효과량 배포 **Figshare 26819743 v1**에는 genome-wide by-gene LFC 약 0.795 GB, targeted by-gene 약 0.078 GB, targeted per-line 약 1.920 GB 등이 있다. 이는 유전자 평균 효과를 빠르게 학습하는 보조 자료이지 raw pseudobulk counts가 아니다. “No guide detected”는 “실제로 non-targeting guide를 받은 세포”와 다를 수 있어 controls 신뢰도를 구분한다.[^feng]

### 5.2 Jiang / Mixscale: 세포 × 자극 조건

**Zenodo 14518762 v2.1**을 기준으로 한다. 이전 v1/v2와 중복 다운로드하지 않는다. 6개 세포주 A549·MCF7·HT29·HAP1·BxPC3·K562, 5개 자극 IFNB·IFNG·TGFB·TNFA·INS의 반응을 다루며 약 2.6M 세포다. 논문의 약 1,500 perturbation-context 실험을 1,500개 고유 표적 유전자라고 읽으면 안 된다.[^jiangpaper]

| 추천 파일 | GB |
|---|---:|
| `Seurat_object_IFNB_Perturb_seq.rds` | 4.327 |
| `Seurat_object_IFNG_Perturb_seq.rds` | 2.916 |
| `Seurat_object_INS_Perturb_seq.rds` | 5.601 |
| `Seurat_object_TGFB_Perturb_seq.rds` | 2.642 |
| `Seurat_object_TNFA_Perturb_seq.rds` | 4.656 |

R에서 Seurat 객체를 열고 **RNA assay의 counts**를 가져온다. `data`, `scale.data`, integrated assay, SCT residual을 raw count로 사용하지 않는다. Seurat v5에서는 layer가 나뉘었는지 확인하고, 실제 사용하는 버전에 맞는 `LayerData`/`GetAssayData`로 counts를 얻는다. sparse Matrix Market와 obs/var metadata로 내보내거나 검증된 AnnData 변환을 사용한다. cell line, stimulus, donor/replicate, target, guide를 유지한다. 함께 있는 `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds`는 **자극 bulk RNA**이며 유전자 perturbation pseudobulk가 아니다.[^jiang]

### 5.3 VCC 2025 및 Arc의 이동한 경로

2025의 공개 train·validation·test는 2026 관점에서는 과거 외부 실험 자료다. 현행 [Arc Virtual Cell Atlas README](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md)의 다음 경로를 이용한다.

```text
gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/
gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/validation/
gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/test/
```

Arc는 2026년 3월 31일 기존 `arc-ctc-tahoe100` 및 `arc-scbasecount` 버킷에서 이동했다고 공지했다. 오래된 논문·블로그의 명령을 복사하면 잘못된 경로를 요청할 수 있다. 새 버킷은 **Requester Pays**이며 자신의 billing project를 지정하는 경로를 안내한다. 무료 제공량 문구를 모든 계정·지역에 무조건 적용되는 무제한 무료 전송으로 해석하지 않는다.[^arcatlas]

```bash
gcloud storage ls -l --billing-project YOUR_PROJECT 'gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/**'
gcloud storage cp --billing-project YOUR_PROJECT 'gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad' data/vcc2025/
```

이 조사는 사용자 billing project로 조회하거나 비용을 발생시키는 다운로드를 수행하지 않았다. 파일 크기는 해당 계정에서 목록을 확인한 후 고정한다. 데이터 크기가 미검증이라는 빈칸을 0 GB로 계산하면 안 된다.

## 6. 반드시 놓치지 말아야 할 표적별 보강 자료

### 6.1 scPerturb 전체 인벤토리 활용

**RNA/protein Zenodo v1.4에는 54개 H5AD, 총 43,042,767,207 bytes**가 있다. 처음부터 전부 gene perturbation 데이터라고 생각해서는 안 된다. 같은 세포의 RNA/protein 파일이 분리되거나 chemical/cytokine 데이터가 포함되어 있다. 이번 비교표에는 54개 파일을 모두 별도 검색 가능한 항목으로 넣었으며, 사용 우선순위·원문·크기·perturbation 종류를 연결했다.[^scperturb]

| 데이터군 | 문맥 / perturbation | 활용과 주의 |
|---|---|---|
| Adamson 2016, 3개 실험 | K562 CRISPRi UPR | small/combination arm 분리; 포함 control 수와 unique target 수 구분 |
| Dixit 2016, 7d/13d/high-MOI | K562 CRISPR KO | 시간과 복합 perturbation 유지; 단일 표적으로 임의 평탄화 금지 |
| Datlinger 2017/2021 | Jurkat CROP-seq/scifi, TCR stimulation | 추가 T-cell 문맥; 자극별 control 매칭 |
| Norman 2019 | K562 CRISPRa single/double | activation·유전자 상호작용 보조; CRISPRi로 재라벨링 금지 |
| Frangieh 2021 RNA | melanoma/TIL, CRISPR KO, RNA+CITE | RNA 약 1.459 GB; protein 파일은 같은 세포의 보조 readout |
| Papalexi 2021 | THP1 ECCITE-seq | pooled RNA 0.147 GB, arrayed RNA 0.052 GB; protein과 cell 중복 |
| Shifrut 2018 | primary T-cell CRISPR KO | 0.871 GB; v1.4의 수정된 배포 선택 |
| Tian 2019 | iPSC / day7 neuron CRISPRi | 0.351 / 0.269 GB; 발달 시점과 문맥 차이를 학습 |
| Tian 2021 | iPSC-neuron CRISPRi / CRISPRa | 0.289 / 0.155 GB; genome-wide fitness arm과 scRNA subset 구분 |
| Wessels 2023 | THP1 Cas13 RNA knockdown, 조합 | 0.219 GB; RNA targeting modality 구분 |
| Joung 2023 atlas/combination | hPSC TF/ORF overexpression | 5.804 / 0.834 GB; 원본 GEO에는 추가 raw/subsample 버전 |
| Sunshine 2023 | Calu-3 CRISPRi + viral exposure | 0.744 GB; 원문은 CRISPRi, 일부 가공 코드의 generic Cas9 표기 교정 |
| Xu 2023 | HEK293 CRISPRi, nascent RNA | 0.337 GB; total/nascent RNA readout을 구분 |
| Lara-Astiaso 2023 | mouse hematopoiesis ex vivo/in vivo/leukemia | 세 파일 총 약 3.017 GB; 종·분화상태 조건부 보조 |
| Liang 2023 OSCAR | mouse liver organoid KO | 0.364 GB; organoid와 species 전이용 |
| Santinha 2023 | mouse in-vivo brain KO | 0.939 GB; 직접 human loss로 혼합하지 않음 |

Gasperini 2019, Schraivogel TAP-seq, Xie Mosaic-seq는 **cis-regulatory element / enhancer perturbation**의 비중이 크다. enhancer와 연결된 유전자가 있어도 그 유전자를 직접 CRISPRi한 관측과 같지 않다. 특히 TAP-seq는 output panel이 제한될 수 있어 18,533개 feature 전체의 관측 데이터로 사용하면 안 된다. Liscovitch-Brauer, Pierce, Mimitou의 ATAC 중심 배포도 추가적인 유전적 반응 정보지만 RNA count supervision과 구분한다.[^scperturbpaper]

Aissa, Chang, Gehring, Srivatsan sci-Plex, Zhao, Lotfollahi의 chemical arm, Cui cytokine atlas, Schiebinger/Weinreb의 시간·분화·cytokine 자료는 직접 gene perturbation 학습의 우선 후보가 아니다. McFarland는 혼합 resource이므로 실제 genetic subset이 포함된 실험과 drug-only subset을 분리한다. **“perturbation”이라는 단어만으로 gene perturbation으로 분류하지 않는다.**

### 6.2 별도 원자료에서 추가할 수 있는 human 데이터

| 연구/자료 | 접근 경로 | 전처리·선택 |
|---|---|---|
| Replogle 2020 direct-capture | 논문 DOI `10.1038/s41587-020-0470-y`의 데이터 연결 | single/combinatorial arm과 guide identity 보존 |
| Jost 2020 gene-expression titration | GEO `GSE132080` | guide별 knockdown 효율을 연속 변수로 보존 |
| Hill 2018 CROP-seq | GEO `GSE108699` | MCF10A KO; count·guide alignment 확인 |
| Leng 2022 astrocytes | GEO `GSE182308`; RAW.tar 표시 484.8 Mb | iPSC astrocyte CRISPRi; vehicle와 ITC를 분리; `GSE182307`은 별도 bulk |
| Dräger 2022 microglia | GEO `GSE178317`, Kampmann resources | CROP-seq subset 선택; CRISPRi/a 및 microglia 상태 보존 |
| Schmidt/Steinhart 2022 T cells | Zenodo `5784651`, GEO `GSE190604` | `CRISPRa-Perturb-seq.zip` **3,055,676,205 bytes**; resting/restimulated Seurat counts |
| Freimer 2022 T-cell networks | GEO `GSE171737` | bulk genetic RNA와 대응 controls; single-cell variance로 해석 금지 |
| Southard 2025 TF activation | Zenodo `15213619`, `15200179` | singlets의 `layers['counts']`; regression mean과 분리 |
| Southard RPE1-E150 | Zenodo `15215414` | 별도 CRISPRi benchmark; TF CRISPRa와 혼동 금지 |
| PerturbME 2026 | Figshare `33264225 v3` | 아래 혼합 feature·선별 편향을 반드시 처리 |

GEO의 `RAW.tar`라는 이름은 종종 이미 처리된 count H5/CSV 묶음이지 sequencing FASTQ가 아니다. 가급적 count matrix + guide calls + cell metadata + control annotation을 먼저 받는다. FASTQ부터 통일 재처리하는 것은 표준화 효과와 계산비용을 비교한 뒤 별도 결정한다. 정확한 bytes가 확인되지 않은 위 보강 자료는 비교표에 미검증으로 남겨 두었고, 원본 accession 페이지를 다운로드의 시작점으로 연결했다.[^geoastro][^schmidt]

### 6.3 2026년 새 자료와 공개 상태

**PerturbME**는 2026년 8월 preprint와 Genentech 공식 코드·Figshare를 확인했다. A375 HLA-high/low FACS-enriched CRISPR-Cas9 KO 자료다. 추천 `PerturbME_A375_HLA_MOI1.h5ad`는 277,152 singlets, 1,430,639,556 bytes다. 전체 QC 353,091 cells 파일은 1,843,931,953 bytes이며 singlet 파일과 겹친다.[^perturbme]

가장 큰 함정은 feature가 **18,670 RNA raw counts + 13 normalized CITE protein**이라는 점이다. RNA/protein을 분리한 후 RNA만 count 정규화·feature mapping에 넣는다. HLA 기준 선별로 얻은 분포는 무선별 세포군의 전사체 분포와 같지 않다. `elasticnet_beta_matrix.csv.gz`의 231 targets는 신뢰할 만한 모델 계수 대상이지 genome-wide library 전체 coverage가 아니다. 해당 coefficient matrix는 발현 pseudobulk가 아니다.

**Maffa 2026**는 100 genes × 16 cancer lines의 새로운 CRISPR KO Perturb-seq라는 사실과 논문은 확인했다. 그러나 원자료 파일 묶음의 직접 다운로드 URL·정확한 bytes·PB 여부는 이번 조사에서 확정하지 못했다. 문맥 다양성이 중요해 추적 목록에 남겼지만 현재 확보 완료된 학습자료로 계상하지 않는다.[^maffa]

**PerturbAI wholebrain**는 공식 HF에 parquet 약 **81.285 GB**, metadata 약 **470.306 MB**가 있다. HF 전체 약 204.95 GB를 받으면 H5AD 중복 표현과 약 46.10 GB 분석 결과까지 함께 받는다. 약 7.7M mouse brain nuclei와 약 2,000 disease-associated targets는 신경계 문맥 보조에 유용하나 human CRISPRi의 직접 대체재는 아니다. 원본 species identity와 1:1 ortholog 여부를 남긴다.[^wholebrain]

Jin 2020, Alda-Catalinas 2020, 새로운 spatial Perturb-seq/DBiT 자료도 비교표에 포함했다. RNA 측정 방식·공간 aggregation·species 차이 때문에 낮은 우선순위에서, 실제 표적 공백이나 문맥 전이 개선이 확인되는 경우 도입한다. 기업의 100M·120M 규모 발표, AI-generated response atlas, chemical-only mega-atlas를 실측 공개 gene perturbation 데이터 수에 더하지 않았다.

## 7. 남는 표적을 메울 탐색 자원과 bulk 데이터

| 자원 | 조사에서 확인한 범위 | 어떻게 활용할지 |
|---|---|---|
| PerturBase | 122 datasets / 46 studies, RNA·multi-omics | human + genetic + RNA로 검색 후 원논문·counts 다운로드; Zenodo 17972932의 작은 파일은 전체 5M-cell 자료가 아님 |
| PerturbSeq.db | 189 datasets / 77 studies, RNA와 ATAC | target 검색으로 놓친 실험을 찾고 원본 accession 확인; 전체 다운로드 가용성 미검증 |
| EMBL-EBI Perturbation Catalogue | Perturb-seq와 fitness/MAVE 등 다양한 readout | readout 필터 필수; 공개 GitHub의 저자 보충표·guide 목록으로 metadata-first matching |
| PerturbAtlas | 13 species, 122,801 RNA libraries | 인간 KO/KD의 bulk expression/DE와 matched control 확보; libraries를 unique gene 수로 읽지 않음 |
| GPSAdb 2.0 | 7,665 groups, 42,235 RNA samples, 2,810 genes, 1,063 cell lines | 현재 metadata와 group별 expression/DE를 이용해 target 공백 채우기 |
| KnockTF 2.0 | 1,468 datasets, 612 TF + 172 TcoF, multiple species | TF 표적의 human KD/KO를 찾아 원본 count 또는 microarray 구분 |
| ENCODE | TF/RBP knockdown·CRISPR RNA-seq | shRNA RNA-seq / CRISPR RNA-seq 실험, matched controls, gene quantification 다운로드 |
| LINCS L1000 | shRNA / CRISPR KO / ORF expression signature | full RNA가 부족한 유전자에 보조 effect prior; 측정 978 landmark와 inferred feature 구분 |
| ARCHS4/GEO | 방대한 bulk RNA | accession 단위 perturbation-control 쌍을 검토; 질병 비교를 genetic intervention으로 오분류하지 않음 |

GPSAdb는 홈페이지 최신 통계 7,665 groups와 전체 DGE ZIP 안내의 **6,245 groups / 7.74GB 표시**가 공존한다. 따라서 ZIP 하나로 최신 7,665개를 전부 받았다고 선언하면 안 된다. package 안의 group IDs를 최신 metadata와 비교하고 빠진 실험은 개별 다운로드한다. 공개 zip은 DE 결과이며 raw RNA count 총량과 다르다.[^gpsa]

LINCS에서는 `trt_sh`, `trt_xpr`, `trt_oe`를 유전자 perturbation으로 식별하고 chemical `trt_cp`와 분리한다. `GSE92742`, `GSE70138`의 metadata에서 필요한 세포주·표적·시간을 먼저 고른 뒤 GCTX column subset을 읽는다. Level 3 normalized expression과 Level 5 standardized signature는 다르며, consensus gene signature는 개별 shRNA 결과와 다르다. 978개 측정 landmark 및 추론 feature를 18,533개 raw RNA 관측값으로 취급하지 않는다.[^lincs]

Bulk RNA는 이미 세포 집단의 측정이므로 pseudobulk를 다시 만드는 개념과 다르다. counts·TPM·microarray intensity·LFC를 구분하고 donor/replicate별 matched control을 통해 효과를 추출한다. 단일세포의 dispersion, detection rate, 다봉 분포는 bulk만으로 회복할 수 없다. mean-effect auxiliary loss나 gene embedding supervision으로 쓰고, single-cell distribution은 실제 단일세포 데이터에서 배운다.

ENCODE의 eCLIP/ChIP peaks, DepMap의 CRISPR gene effect fitness, CCLE baseline expression은 각각 결합위치·성장효과·무처리 발현이다. 유전자 perturbation 후 RNA 반응과 직접 동일하지 않다. 카탈로그의 분량을 키우기 위해 모두 supervised perturbation examples로 섞지 않는다.[^knocktf][^catalogs]

## 8. 표적 coverage를 정확히 계산하는 절차

### 8.1 존재 여부의 네 단계

1. **Measured feature:** `var`에 그 유전자가 측정돼 있다. perturbation 경험은 증명하지 못한다.
2. **Designed target:** guide library에 그 유전자가 있다. QC 후 충분한 세포가 있는지는 알 수 없다.
3. **Observed target:** `obs`에 유효한 perturbation assignment를 가진 세포가 있다.
4. **Usable response:** 대응 controls, 충분한 세포·guide, 효율·batch 정보가 있어 학습에 사용할 수 있다.

이번 자료집의 `target_registry.csv`는 **63,252개 표적-문맥 행, 14,553개 고유 원문 symbol**을 포함한다. 별칭 통합 후의 unique HGNC gene 수가 아니다. Replogle/Nadig 저자 부록, CD4 DE-tested metadata, KOLF browser 목록을 수집했으며 증거 종류를 행마다 기록했다. **Orion/Feng/Jiang 및 모든 legacy dataset의 전수 observed-target registry까지 완료한 것은 아니다.** 이들의 guide/obs 목록을 확보하는 다음 작업을 위한 접근 경로를 함께 제공한다.

| 검증 출처 / 문맥 | 고유 표적 symbol 수 | 수치의 정의 |
|---|---:|---|
| Replogle K562 genome-wide | 9,866 | 저자 TargetingGuides + 관측 세포 수 > 0 |
| Replogle K562 essential | 2,057 | 동일 |
| Replogle RPE1 | 2,393 | 동일 |
| Nadig HepG2 | 2,393 | 동일 |
| Nadig Jurkat | 2,393 | 동일 |
| CD4 Rest | 11,287 | 저자 DE-tested 표적-조건 summary |
| CD4 Stim8hr | 11,415 | 동일 |
| CD4 Stim48hr | 11,281 | 동일 |
| KOLF browser | 10,167 | 브라우저에 나오는 perturbation 목록; full H5AD 전수 obs 확인 아님 |

이 수치들은 **VCC 300개와의 교집합 수가 아니다.** CD4 PB의 categorical vocabulary 약 12.7k와 조건별 DE-tested 수가 다른 것도 정상이다. 실험 설계, 실제 관측, DE 분석 포함, 강한 효과 subset은 서로 다른 분모다.[^targettable][^cd4readme][^kolf]

### 8.2 실제 표적 파일을 받은 뒤

```bash
python vcc_data.py coverage --targets official_controls/pert_counts.csv --column target_gene --output coverage
```

도구는 기본적으로 300개 고유 표적을 기대하며 잘못된 18,533개 파일을 거절한다. 결과는 target별 hit, dataset/context별 교집합, 현재 부분 registry에서 미매칭인 target을 제공한다. 미매칭을 “공개 데이터에 존재하지 않음”이라고 해석하지 않는다. 해당 유전자를 Orion guide calls → KOLF full obs → Feng/Jiang → scPerturb/별도 human studies → bulk catalog 순서로 확장 조사한다.

Gene mapping은 Ensembl stable ID를 중심으로 version suffix를 분리하고, 원문 symbol과 현재 승인 symbol을 함께 보존한다. 예를 들어 AARS/AARS1 같은 이름 변경을 검토하되, 여러 유전자에 대응하는 alias는 자동 통합하지 않는다. `HLA-DRA`처럼 정상 gene symbol에 들어 있는 하이픈을 복합 perturbation 구분자로 사용하지 않는다. 같은 유전자에 대한 두 guide와 두 유전자 동시 perturbation도 구분한다.

### 8.3 최종 학습 후보 판정

같은 표적이 발견되면 무조건 보존할 **후보**로 등록하되 실제 학습에는 다음 정보를 함께 반영한다: perturbation modality, cell context, cells-per-target, independent guides, on-target KD, matched controls, stimulation/time, raw RNA 여부, QC, 라이선스, 다른 배포본과의 중복 관계. 표적 평균·guide별 평균·single-cell distribution을 각각 평가하면 small-study의 특이적인 artifact를 방지하는 데 도움이 된다.

문맥 전이는 외부 데이터에서 **한 세포주/문맥 전체를 제외한 검증**으로 시험한다. 모든 문맥의 같은 유전자를 학습에 본 경우와 한 유전자 전체를 제외한 경우를 나누어 본다. 이것은 “같은 표적을 미리 보면 얼마나 도움이 되는가”와 “새 문맥에서도 옮겨가는가”를 별도로 측정하는 방법이다. 근거 없이 모든 추가 데이터가 성능을 올린다고 가정하지 않는다.

## 9. 공통 전처리 설계

### 9.1 원자료와 표준 스키마

원자료는 immutable하게 보존하고 manifest에 URL, accession, release/commit, 원본 filename, bytes, 제공 checksum, 취득일을 기록한다. `.h5ad` 확장자가 같다고 `.X`의 의미가 같지는 않다. `X`, `layers/counts`, `raw/X`, `var`, `obs`를 먼저 조사한다. 정수라는 형식 검사는 필요하지만 그것만으로 raw UMI임을 증명하지 못한다.

표준 obs에는 `dataset_id`, `experiment_id`, `cell_context`, `donor`, `biological_replicate`, `technical_batch`, `stimulation`, `time`, `perturbation_mode`, `target_set`, `guide_id`, `is_control`, `assignment_confidence`, `n_counts`, `n_genes`, `qc_flags`를 권장한다. 원문에 없는 donor나 replicate를 만들어 넣지 않는다. 값이 없으면 unknown으로 남긴다.

feature mapping은 full RNA feature union을 먼저 유지하고, 모델 출력에는 공식 18,533개 순서를 고정한다. 어떤 assay에서 측정하지 않은 gene은 **unmeasured mask**로 처리한다. 그 열을 편의상 zero-fill했다면 관측 zero와 구별하는 loss mask도 필요하다. 모든 데이터를 처음부터 HVG 교집합으로 자르면 나중에 원래 출력 범위를 회복할 수 없다.

### 9.2 QC와 모달리티

guide assignment, doublet, cell QC, 최소 read depth는 각 실험의 기준과 분포를 확인한다. 전 연구에 하나의 mitochondrial threshold를 기계적으로 적용하지 않는다. 효율 낮은 guide를 전부 폐기하면 실제 약한 knockdown 분포를 없앨 수 있어 보존본과 학습 선택본을 분리한다. “유의한 DE가 있다”는 조건만으로 표적을 골라 학습하면 효과가 없는 target이 체계적으로 빠진다.

CRISPRi는 transcriptional repression, CRISPR KO는 DNA damage·copy-number·selection effect까지 포함할 수 있고, CRISPRa/ORF는 activation이다. KO 효과에 단순히 부호를 바꾸면 activation을 얻는 것도 아니다. modality별 데이터셋을 보존한 뒤 condition token, 별도 head, 보조 loss 등의 통합 방법을 실제 검증으로 선택한다.

Count 입력은 raw를 유지하고 학습 view에서 필요에 따라 library-size normalization/log1p를 만든다. log1p를 쓰더라도 raw total count와 original features를 보존한다. 여러 연구를 통째로 batch correction해 perturbation signal 자체를 제거하지 않는다. 대조군을 이용해 기술적 차이를 추정하고 cell context·assay·batch를 모델에 명시하는 방식부터 검증한다.

### 9.3 Pseudobulk의 정확한 정의

집단 b, 유전자 g의 count sum과 cell mean은 다음과 같이 다르다.

\[
S_{bg}=\sum_{i\in b}X_{ig},\qquad \bar X_{bg}=S_{bg}/n_b.
\]

`mean(log1p(X))`, `log1p(mean(X))`, library-size-normalized mean, NTC 대비 LFC, Z-score는 별개의 값이다. 한 파일에 섞어서 같은 “pseudobulk”라는 라벨을 주지 않는다. PB 저장에는 집계 정의, n_cells, target/guide, donor/replicate, condition/time, control matching key를 반드시 포함한다.

실용적인 기본 group key는 `(dataset, cell_context, donor, biological_replicate, stimulation, time, target, guide)`이며 기술 run을 추가한다. 자료에 donor가 없으면 GEM group을 biological replicate라고 명명하지 않는다. DE 추론에서는 실제 실험의 독립 단위를 사용한다. 평균 모델 학습에서는 guide별/target별 aggregation 방식과 셀 수 차이를 명시한다.

제공 스크립트 `pseudobulk_h5.py`는 dense/CSR RNA counts를 row chunk로 읽어 sum과 mean_raw를 만든다. 비정수·음수·결측 group label을 거절한다. **메타데이터와 최종 sparse PB는 RAM에 둔다**는 한계가 있다. CD4 전체 1.736TB를 이 도구 하나로 저메모리 처리한다고 기대하지 말고, 공개 PB를 이용하거나 donor/조건별 partition과 분산 aggregation을 사용한다.

### 9.4 중복 제거와 평가 오염 방지

원본 Figshare/GEO ↔ scPerturb ↔ Arc HF ↔ benchmark 가공본의 계보를 연결한다. 파일 checksum이 달라도 같은 세포를 다른 dtype·정규화·압축으로 저장했을 수 있다. `(study, experiment, donor, GEM/run, original cell barcode)`를 중심으로 원본 identity를 만든다. Strong subset, singlet subset, RNA/protein view, pseudobulk/DE view도 원본 single-cell과 연결한다.

서로 다른 날짜의 진짜 반복 실험과 같은 세포의 재배포는 다르다. train/validation split 전에 원본 identity를 기준으로 그룹을 나누고, 같은 세포의 다른 view가 검증쪽에 섞이지 않도록 한다. 문맥을 통째로 뺀 검증과 perturbation을 통째로 뺀 검증을 병행한다. 챌린지의 대조군을 사용하더라도 숨겨진 perturbation response를 아는 것처럼 데이터를 합성하거나 라벨을 만들어 넣지 않는다.

## 10. 단계별 실행안

**단계 A — 메타데이터와 표적 점검.** 실제 `pert_counts.csv`를 확보하고 부분 target registry와 먼저 대조한다. Orion guide library/calls, KOLF full obs, Feng/Jiang metadata를 추가해 target × context × modality별 cells와 controls를 센다. 이 단계는 대형 RNA matrix를 모두 받기 전에 진행할 수 있다.

**단계 B — 작은 raw 기반 학습 코어.** Replogle·Nadig scPerturb 다섯 파일 **13.733 GB**를 받는다. 필요하면 Replogle raw PB 약 0.55 GB로 mean-effect baseline을 먼저 시험한다. Feng의 7.926 GB와 Jiang의 20.142 GB를 추가하면 이 선택 묶음은 약 **41.801 GB**다. 서로 다른 표적·문맥에 대한 부분 loss와 원본 중복 제거를 먼저 구현한다.

**단계 C — 넓은 genome coverage와 Flex.** 위 묶음에 CD4 PB 44.567 GB, Orion parquet 126.260 GB, KOLF full QC 189.393 GB를 더하면 약 **402.020 GB**의 선택 다운로드가 된다. 소형 metadata와 처리 산출물은 별도다. 이것은 계정별 미확인 VCC 2025, 전체 CD4 single-cell, 중복 원본 H5AD, raw FASTQ를 포함하지 않은 계산이다.

**단계 D — 미매칭 target의 집중 확보.** small human scPerturb·neuron/microglia/astrocyte·T-cell screen과 bulk RNA에서 공백을 메운다. CRISPRa, Cas13, ORF, mouse, spatial 자료는 coverage와 검증 개선이 확인될 때 조건부 보조로 추가한다. Maffa처럼 접근 경로가 확정되지 않은 최신 연구는 별도 추적하며 현재 coverage 분자에 넣지 않는다.

저장 공간은 다운로드 bytes에 더해 변환 전후 파일, PB, 인덱스, 학습 cache를 고려해야 한다. “다운로드 크기의 두 배”는 보장된 upper bound가 아니다. dense materialization은 훨씬 커질 수 있고 gzip CSV를 메모리에 통째로 읽으면 실패할 수 있다. 원본 보존과 중간 산출물 관리 계획을 먼저 정한다.

## 11. 조사에서 남겨 둔 불확실성

**가장 중요한 미해결 입력은 실제 300개 validation perturbation 표적 목록이다.** 이를 받기 전에는 어떤 target을 빠짐없이 확보했다고 선언할 수 없다. 더불어 Orion 등 모든 데이터의 전체 obs를 내려받아 QC 후 유전자별 셀 수까지 완전 감사한 상태는 아니다. 데이터 포털·API metadata·저자 부록·가공 코드·일부 원격 H5 schema를 서로 대조한 조사다.

946개 항목의 전체 원자료 전송을 수행한 것은 아니므로 모든 URL의 장시간 전송 성공을 보증하지 않는다. 인증, 지역, 네트워크, rate limit, Requester Pays 설정에 따라 접근이 달라질 수 있다. 원본 API bytes를 확보한 파일은 manifest로 재검증하고, 미검증 용량은 빈 값으로 유지했다. embargo·미공개·논문만 확인된 후보는 현재 사용 가능한 데이터와 분리했다.

따라서 이 자료집은 대형 데이터의 중요한 누락과 잘못된 버전 선택을 줄이고, 표적별 남은 공백을 드러내는 기반이다. “공개된 모든 gene perturbation 자료를 수학적으로 완전 열거했다”는 주장이나 실제 VCC 표적 coverage 100% 주장은 하지 않는다. 실제 표적 목록이 주어지면 부분 registry와 추가 metadata를 결합해 **gene별 확보 데이터, 추천 버전, 최소 다운로드, 남는 공백**으로 최종 계획을 좁힐 수 있다.

## 출처

[^vcc]: Arc Institute, [Virtual Cell Challenge 2026 announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026); [VCC CLI wiki](https://vcc-cli-wiki.virtualcellchallenge.org/); [vcc-cli PyPI release](https://pypi.org/project/vcc-cli/). 공식 안내 및 공개 CLI 0.2.0의 feature/perturbation validation 구현을 대조.
[^replogle]: Replogle et al., [Mapping information-rich genotype–phenotype landscapes with genome-scale Perturb-seq](https://doi.org/10.1016/j.cell.2022.05.013); [저자 데이터 포털](https://gwps.wi.mit.edu/).
[^reploglefiles]: [Replogle Figshare 20029387 v1](https://doi.org/10.25452/figshare.plus.20029387.v1); [파일 bytes·ID·MD5 API](https://api.figshare.com/v2/articles/20029387).
[^scperturb]: [scPerturb Zenodo 13350497 v1.4](https://zenodo.org/records/13350497); [release API](https://zenodo.org/api/records/13350497).
[^scperturbpaper]: Peidli et al., [scPerturb: harmonized single-cell perturbation data](https://doi.org/10.1038/s41592-023-02144-y); [저자 source repository](https://github.com/sanderlab/scPerturb).
[^scprocessing]: [Replogle raw processing script](https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/ReplogleWeissman2022.py); [Nadig raw processing notebook](https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/notebooks/NadigOConner2024.ipynb).
[^targettable]: Nadig et al. 저자 Supplementary Tables ST2–ST6, [EMBL-EBI에 보존된 원문 부록](https://github.com/EMBL-EBI-ABC/PerturbationCatalogue/blob/main/data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx); [GSE264667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667); [DE 결과 Figshare](https://api.figshare.com/v2/articles/29498366).
[^arcpre]: [Arc Replogle-Nadig-Preprint](https://huggingface.co/datasets/arcinstitute/Replogle-Nadig-Preprint/tree/main); [현재 Replogle-Nadig tree](https://huggingface.co/datasets/arcinstitute/Replogle-Nadig/tree/main). 배포 tree와 H5 header/obs/var range-read 결과를 대조.
[^arcfiltered]: [State-Replogle-Filtered](https://huggingface.co/datasets/arcinstitute/State-Replogle-Filtered/tree/main); [rhaister 파생본](https://huggingface.co/datasets/tahoebio/replogle-nadig-de-rhaister/tree/main). `.X` 표본과 feature 교집합은 본 조사 계산.
[^orion]: [Xaira 공식 HF](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion); [Orion Figshare v3 metadata](https://api.figshare.com/v2/articles/29190726).
[^orioncard]: [Orion schema·loading card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/README.md); [배포 라이선스](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/LICENSE.md).
[^kolf]: [KOLF 저자 browser](https://y-doctor.github.io/KOLF2.1J_Perturbation_Cell_Atlas/); [Figshare 27261219 v1](https://api.figshare.com/v2/articles/27261219).
[^kolfrepo]: [KOLF 저자 코드·접근 안내](https://github.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas); [Boom5426 mirror](https://huggingface.co/datasets/Boom5426/KOLF/tree/main).
[^cd4]: [CD4 저자 분석 repository](https://github.com/emdann/GWT_perturbseq_analysis_2025); [CZI dataset card](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq).
[^cd4readme]: [CD4 데이터 사전](https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/data_sharing_readme.md); [DE-tested metadata](https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/suppl_tables/DE_stats.suppl_table.csv); [curated guide metadata](https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/suppl_tables/sgrna_library_metadata.suppl_table.csv).
[^moonen]: [Moonen/Claringbould Zenodo 18876252](https://zenodo.org/records/18876252); embargo와 공개 예정 schema 확인.
[^fengpaper]: Feng et al., [A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00332-5), DOI [10.1016/j.xgen.2025.101076](https://doi.org/10.1016/j.xgen.2025.101076).
[^feng]: [Feng UMI counts v2](https://api.figshare.com/v2/articles/27989294); [Feng LFC/statistical results v1](https://api.figshare.com/v2/articles/26819743).
[^jiangpaper]: Jiang et al., [Nature Cell Biology 2025](https://doi.org/10.1038/s41556-025-01622-z); [Mixscale repository](https://github.com/satijalab/Mixscale).
[^jiang]: [Jiang Zenodo 14518762 v2.1](https://zenodo.org/records/14518762); [GSE281048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048).
[^arcatlas]: [Arc Virtual Cell Atlas repository](https://github.com/ArcInstitute/arc-virtual-cell-atlas); [VCC 2025 datasets](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md).
[^southard]: Southard et al., [Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states](https://doi.org/10.1038/s41588-025-02284-1); [RPE1 processed](https://zenodo.org/records/15213619); [Hs27 processed](https://zenodo.org/records/15200179).
[^geoastro]: [Leng astrocyte CROP-seq GSE182308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182308); [Sunshine Calu-3 CRISPRi GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240); [Joung GSE217460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217460).
[^schmidt]: [Schmidt/Steinhart code and data archive](https://zenodo.org/records/5784651); [GSE190604](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190604).
[^perturbme]: [Genentech PerturbME repository](https://github.com/Genentech/PerturbME); [Figshare 33264225 v3](https://api.figshare.com/v2/articles/33264225); preprint [10.64898/2026.08.17.745330](https://doi.org/10.64898/2026.08.17.745330).
[^maffa]: Maffa et al., [Dissecting context-dependent cancer vulnerabilities using Perturb-seq](https://www.biorxiv.org/content/10.64898/2026.08.24.746802v1), 2026-08-25.
[^wholebrain]: [PerturbAI wholebrain official HF](https://huggingface.co/datasets/perturbai/wholebrain_crispr_atlas); [저자 분석 코드](https://github.com/jinlabneurogenomics/wholebrainperturbseq).
[^gpsa]: [GPSAdb 현재 홈페이지·Download Manager](https://www.gpsadb.com/). 최신 headline과 legacy ZIP 설명의 group 수 차이를 함께 기록.
[^knocktf]: Feng et al., [KnockTF 2.0](https://doi.org/10.1093/nar/gkad1016); [ENCODE](https://www.encodeproject.org/).
[^lincs]: [LINCS 공식 다운로드 경로 안내](https://lincsproject.org/LINCS/tools/workflows/find-the-best-place-to-obtain-the-lincs-l1000-data); [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742); [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138).
[^catalogs]: [PerturBase](http://www.perturbase.cn), [논문 DOI](https://doi.org/10.1093/nar/gkae858); [PerturbSeq.db](http://bioailab.com/PerturbSeq.db); [EMBL-EBI Perturbation Catalogue](https://www.ebi.ac.uk/perturbation-catalogue); [PerturbAtlas](https://perturbatlas.kratoss.site), [논문 DOI](https://doi.org/10.1093/nar/gkae851); [ARCHS4](https://archs4.org/download).
