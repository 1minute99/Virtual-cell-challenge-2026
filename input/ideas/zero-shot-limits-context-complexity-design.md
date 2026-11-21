---
title: "Is zero-shot prediction of perturbation responses in unseen cell lines limited by data (measurement ceiling, non-transferable interaction, number of contexts) rather than by model capacity?"
type: design
created: 2026-11-05
updated: 2026-11-05
tags: [virtual-cell, perturbation-prediction, zero-shot-transfer, sample-complexity, negative-results, perturb-seq]
status: approved
author: mixed
supersedes: "[[zero-shot-transfer-destination-state-design]]"
---

# 처음 보는 세포주로의 교란 반응 예측은 모델이 아니라 데이터에 막혀 있는가: 상한과 맥락 표본 복잡도

(가제: *Data-limited, not model-limited: ceilings and context sample complexity for zero-shot perturbation
prediction in unseen cell lines*)

이 설계는 [[zero-shot-transfer-destination-state-design]]을 **통합·대체**한다. 옛 설계의 E1(목적지 상태 잔차와
M3)은 RQ3, E2(맥락 수 스케일링)는 RQ2의 경험 부분, E3(같은 세포주 데이터)는 RQ4, E4(공식 성적)는 그대로
유지하고, E5(PIE/STATE 플러그인)는 제외한다. (사용자 결정 2026-11-01 "B로 진행")

## Research Question

**주 질문 (반증 가능):** 공개 스크린으로 처음 보는 세포주의 CRISPRi 반응을 예측할 때, 달성 가능한 정확도는 주로
(i) 측정 잡음 천장과 (ii) 다른 세포주에서 옮겨 올 수 없는 세포주×교란 상호작용 성분에 막혀 있는가? 그렇다면 다음
세 가지가 함께 관찰되어야 한다.

- (a) 최선의 전이 파이프라인이 이미 zero-shot 전이 상한의 큰 부분에 도달해 있다.
- (b) 학습형 맥락 모델의 이득은 현재 공개 맥락 수(약 10개 세포주)에서 맥락이 늘어도 커지지 않으며, 이는 맥락 표본
  복잡도 이론이 예측하는 본전 맥락 수 C*보다 C가 작기 때문이다.
- (c) 남은 오차는 이론이 예측하는 곳 — 약한 표적과 목적지 자신의 control 상태 축 — 에 몰린다.

**반증 조건:** 어떤 방법이 추정된 전이 상한을 넘거나(상한 추정이 틀림), 우리 예측과 상한 사이의 격차가 크고 그
격차가 전이 가능한 성분(a_t + b_c)에 몰려 있으면 "데이터 한계" 주장은 기각되고 "모델 한계"로 보고한다. Tahoe에서
C ≤ 24로 적합한 명제 2가 C = 32, 40의 이득을 예측하지 못하면(예측 구간 밖) 정량 이론은 기각되고 명제 1(식별
가능성)과 경험 곡선만 보고한다.

**보조 질문 (기여 순서):**
1. **RQ2 맥락 표본 복잡도 (주 기여):** 상호작용을 control 기술자로부터 배우는 데 필요한 세포주 수 C*는 얼마이며,
   현재 공개 유전자 교란 맥락 수는 그보다 작은가?
2. **RQ1 상한 대비 위치:** 측정 상한과 zero-shot 전이 상한은 얼마이고, 실제 방법들(우리 파이프라인, 학습형 구조
   17종, 단순 평균 전이)은 그 사이 어디에 있는가?
3. **RQ3 남은 오차의 위치와 보정:** 전이 예측 잔차는 잡음 통제 후에도 목적지 control PC에 몰리는가(다른 세포주 PC,
   섞은 control PC 대비)? 그 축을 쓰는 보정(M3)이 대안 보정보다 PDS를 더 올리는가?
4. **RQ4 같은 세포주 데이터의 가치:** 다른 랩의 같은 세포주 스크린 하나(H1GW)는 같은 지표로 잰 모델링 이득, 맥락
   하나 추가의 한계 가치보다 큰가?

## Relevance / Forschungsstand (brief)

"학습형 모델이 단순 기준선을 못 이긴다"는 이미 보고되어 있다(선형 기준선 Nature Methods 2025; VCC 2025 정리;
AdaPert ICML 2026; STATE는 유전자 교란에서 기준선 대비 +5 %). 반응을 공통·교란·세포주·상호작용 성분으로 잡음
보정 분해하고 control PC가 잔차 분산을 얼마나 잡는지는 Molina & Zhang(bioRxiv 2026-07)이, 상태별 잔차와 잡음
바닥은 Science Advances 벤치마크가 다뤘다. 세포 수 기준의 표본 크기는 control-aware cell quota(bioRxiv 2026-10)가
냈고, "맥락 다양성이 병목"이라는 정성적 입장은 Dibaeinia et al.(ICML 2026 position)이 냈다. **비어 있는 부분:**
(1) 맥락(세포주) 수 기준의 표본 복잡도와 본전 맥락 수의 정량값, (2) 분해를 zero-shot 전이 상한으로 바꿔 실제
방법들을 그 위에 놓은 결과, (3) 반응이 아니라 **전이 예측 잔차**를 다른 세포주 축 귀무와 비교한 잡음 통제 분석과
보정 비교, (4) 같은 세포주 데이터와 모델링 이득의 동일 지표 비교. 공식 대회 채점기로 검증된 파이프라인(VCC 2026
검증 단계 1,635팀 중 129위, 상위 8 %)이 실증 근거를 더한다.

## Method

**표기:** δ_{c,t} = μ + a_t + b_c + i_{c,t} + ε_{c,t} (공통 반응, 보존 교란 성분, 세포주 성분, 상호작용, 측정 잡음;
Molina & Zhang 분해를 채택·인용). zero-shot에서 새 세포주 c*에 대해 쓸 수 있는 것은 소스 C개 세포주의 반응과
c*의 control 세포(기술자 x_{c*} ∈ ℝ^d)뿐이다.

**RQ2 — 이론:**
- *명제 1 (식별 가능성):* 상호작용을 기술자의 함수로 학습할 때 관측되는 기술자는 C개 점뿐이다. 선형 모형에서
  식별 가능한 방향은 소스 기술자들이 펼치는 C−1 차원뿐이며, 그 밖의 성분은 어떤 추정기로도 복원되지 않는다.
  비선형 판본은 매끄러움 가정 아래 표본 수 = C인 비모수 수렴률로 부록에서 일반화한다.
- *명제 2 (본전점):* 선형-가우스 모형에서 학습된 맥락 보정의 기대 이득은
  이득(C) ≈ R²·Var(i) − σ²·d/(C−d−1), 본전 맥락 수 C* = d + 1 + d·σ²/(R²·Var(i)).

**RQ2 — 검증:**
- Tahoe 약물 데이터(48개 세포주)에서 맥락 수 곡선을 C ∈ {6, 8, 12, 16, 24, 32, 40}, seed 3개로 다시 그린다(기존
  FiLM 7k 모델, 동일 fold·예산). **사전 예측:** C ≤ 24로 명제 2를 적합하고 C = 32, 40의 이득을 예측해 실측과 비교.
- 독립 교차 확인: 적합된 d와 control 기술자의 유효 차원(참여비)을 비교(2배 이내면 일치).
- 유전자 교란에 적용: 공개 유전자 교란 세포주들에서 d, σ², R²를 따로 추정해 C*를 계산하고, 공개 C ≈ 10과 비교.
  예측(C < C* → 학습형 맥락 모델 이득 ≤ 0)을 기존 경험 결과(학습형 6→35 평평한 곡선, 17개 구조 비교를 seed 3개로
  재확인)와 대조한다. 대비로 학습 없는 전이 파이프라인의 소스 수 곡선(무작위 부분집합 반복 + leave-one-family-out).

**RQ1 — 상한:** 7개 하네스 패널에서 반쪽 분할로 잡음 보정 분산 성분을 추정(공통 표적은 소스 스크린과의 겹침).
측정 상한 = 반쪽 분할 신뢰도, 전이 상한 = 재현 가능 분산 중 (a + b)의 몫. PDS 판본은 "다른 모든 세포주 정답의 평균"
오라클 예측. 실제 방법의 상한 대비 위치를 잡음 보정 상관(전 패널)과 실제 채점기 PDS(H1)로 보고.

**RQ3 — 잔차와 M3:** 반쪽 분할 교차곱 추정기, control 세 조각 분할(PC 추정용 / 반쪽 A 기준선 / 반쪽 B 기준선).
목적지 control PC(k = 8, 24, 50) vs 다른 세포주 PC, 유전자별 섞은 control PC. 음성 대조(control에서 뽑은 가짜
표적). 이어서 M3 vs 다른 세포주 PC 보정, 섞은 PC 보정, 예측 행렬 PC 보정, Systema류 보정(η 격자 중 최선), 보정
없음을 PDS로 비교(7개 패널 평균 + H1 실제 채점기).

**RQ4 — 같은 세포주 데이터:** H1 패널, 실제 채점기, H1GW 비중 w ∈ {0, 0.25, 0.5, 0.75, 1}. 같은 기준·지표에서 M2,
M3 제거 효과, 맥락 하나의 한계 가치(RQ2 경험 곡선), 최선의 구조 이득과 비교. 관련 세포주 사례(KOLF←H1GW,
OR_HCT116←DLD-1)는 여유가 있을 때만.

**공식 근거 (E4):** 공식 제출 기록(상위 8 %), 단일 변경 제출로만 구성요소 이득 귀속(M2 +0.0106, M3 묶음),
로컬 → 공식 격차 사례(FiLM 로컬 +0.0247 → 공식 +0.0013).

**부록 — drug→gene:** 같은 표적의 약물–CRISPR 일치도(Spearman 0.006–0.038, A549·HT-29, 401조건)로 표적
정체성 매개 전이의 상관 상한을 유도. 맥락 벡터 섞임 버그로 무효가 된 기존 실험은 재실행하지 않고 "검증되지 않음"으로
명시.

**제안 방법(파이프라인) 정의:** M1 다중 소스 전이(소스별 fold change 융합, min_cells·pseudocount 앙상블), M2 일치도
기반 크기 보정(β = 0.5, MSE kink 예산), M3 목적지 control PC deflation(하위 50 % 표적, k = 8, η = 0.75). 대회 최종
시스템(S18/F2) 설정을 고정해 사용.

**왜 이 접근인가:** 새 모델을 만들지 않고 이미 있는 데이터와 결과를 하나의 이론(명제 1·2) 아래에 꿰어 "왜 안
되는가"를 정량으로 답한다. 이론에는 사전 예측(Tahoe C = 32, 40)과 독립 교차 확인(d)을 붙여 사후 적합이 아님을
보인다. 추가 계산은 모두 CPU로 가능하다.

## Data Sources / Material Basis

- **유전자 교란 소스 스크린:** Replogle 2022 K562 GWPS·essential·RPE1, Nadig 2025 HepG2·Jurkat, X-Atlas/Orion
  HCT116·HEK293T, CD4 T세포 genome-scale CRISPRi, KOLF2.1J iPSC, HipSci(Feng), VCC2025 H1, H1GW(Pan et al. 2026,
  GSE295214), DLD-1(GSE337988).
- **하네스 패널(세포 수준 정답 있음):** H1 `vcc2026_data/reports/transfer_v2_2026-09-23/h1/real_thin.h5ad`,
  OR_HCT116X `vcc2026_data/reports/sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad`, KOLF
  `checkpoints/state/kolf_real.h5ad`, RN Jurkat `checkpoints/state/jurkat_real.h5ad`, RN K562·RPE1·HepG2
  `vcc2026_data/raw/scperturb/{ReplogleWeissman2022_K562_essential,ReplogleWeissman2022_rpe1,NadigOConner2024_hepg2}.h5ad`.
  정확 오프라인 평가 팩 `vcc2026_data/reports/era_2026-10-07/packs`, 평가기 `vcc2026_data/scripts/sch_sumeval.py`.
- **약물 데이터(RQ2 검증):** Tahoe-100M 로컬 부분집합(443만 세포, 48개 세포주, 13,011 pseudobulk 조건)
  `vcc2026_data/reports/simple_residual_2026-09-16/data/tahoe/table.npz`, `vcc2026_data/data/processed/drug_context_pretrain/tahoe.npz`;
  GSE306429(10개 맥락); 기존 결과 `vcc2026_data/reports/drug_context_pretrain_2026-09-20/`.
- **일치도(부록):** `vcc2026_data/reports/drug_gene_concordance_2026-09-16/` (Tahoe 150약물 × Jiang CRISPRi, LINCS shRNA).
- **기존 경험 결과(재현·정리 대상):** 17개 구조 비교 `vcc2026_data/reports/arch_screen_2026-09-22/` (104 runs, seed 0),
  학습형 맥락 곡선 `vcc2026_data/reports/ext_contexts_2026-09-21/`, 천장·전이 상한 `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`,
  부정 결과 기록 `vcc2026_data/reports/PROJECT_HISTORY.md`, 공식 점수 `vcc2026_data/reports/official_scores.csv`(+ `submission/official_scores_new.csv`).
- **채점기:** cell-eval2 `vcc2026` preset(공개). 실제 채점기는 회당 약 25분·약 15 GB, 한 번에 하나.
- **핵심 문헌(A):** Molina & Zhang 2026(분해), control-aware cell quota 2026, Science Advances 벤치마크, Dibaeinia
  et al. 2026(position), Ahlmann-Eltze 2025, AdaPert 2026, STATE 2025/2026, Tahoe-100M 2025/2026, Systema 2025,
  VCC 2026 commentary(Cell), Replogle 2022, Nadig 2025, X-Atlas/Orion 2025, Pan 2026.

## Output Target

☐ Book chapter   ☐ Journal article   ☐ Grant proposal   ☐ Talk
☒ Other: **arXiv 프리프린트**(cs.LG 주분류, q-bio.GN 교차 게시; 보증 실패 시 bioRxiv), 본문 8~10쪽 + 부록, 코드
공개(GitHub). 목표 제출 2026-11-28~12-01(박사 지원 12월 마감용). v2 후보: ICLR 2027 워크숍(권장 마감 2027-02-01).

## Expected Results

- **의미 있는 긍정 결과:** Tahoe에서 C ≤ 24로 적합한 명제 2가 C = 32, 40의 이득을 예측 구간 안에서 맞히고, 적합된
  d가 control 기술자의 유효 차원과 2배 이내로 일치한다. 유전자 교란에서 추정한 C*가 공개 C ≈ 10보다 크다(→ 학습형
  맥락 모델의 실패를 정량으로 설명). 우리 파이프라인이 전이 상한의 큰 부분에 도달해 있고 남은 격차가 상호작용·잡음
  성분에 있다. 전이 잔차가 목적지 control PC에 몰리고 M3가 대안 보정보다 낫다. H1GW 하나의 가치가 모델링 이득과
  맥락 하나의 가치보다 크다.
- **귀무 결과의 해석:** 사전 예측이 빗나가면 정량 이론을 접고 명제 1(식별 가능성) + 경험 곡선만 보고한다. 파이프라인과
  상한 사이 격차가 크고 전이 가능한 성분에 몰려 있으면 "모델 한계"로 결론을 뒤집어 보고한다(그 자체로 의미 있는
  결과). RQ3에서 잡음 통제 후 집중이 사라지면 M3를 잡음/허브 억제로 재해석한다.

## Ethical / Methodological Risks

- **선행 연구 겹침:** Molina & Zhang 2026과 RQ1 분해, RQ3 control PC 분석이 겹친다. 원문을 1주차에 읽고 차이
  (전이 예측 잔차 대상, 다른 세포주 축 귀무, 보정 비교, 상한 대비 위치)를 확정한다. 겹침이 더 크면 RQ3를 축소한다.
- **이론의 단순성:** 명제 2는 선형-가우스 모형이다. 사전 예측과 독립 교차 확인으로 방어하고, 유전자 쪽 모수는
  유전자 데이터에서 따로 추정한다(약물 쪽 값을 옮기지 않는다). Tahoe는 1억 세포 전체가 아니라 로컬 부분집합임을 명시.
- **기존 결과의 약점:** 17개 구조 비교는 seed 0뿐 → 핵심 구조를 seed 3개로 재실행. 버그로 무효가 된 drug→gene
  실험은 인용하지 않는다. 학습형 6→35 곡선은 한 가지 맥락 순서·약한 기준선 → 그 한계를 명시.
- **측정 잡음:** 얕은 패널(OR 163셀, RN 53~96셀)은 지표가 불안정 → H1 주 패널, 반쪽 분할로 통제, 잡음 공분산이
  control PC 위에 놓이는 구조적 편향을 control 세 조각 분할과 음성 대조로 분리.
- **오프라인 평가기 범위:** 예산 불변 arm에서만 정확 → 그 밖은 실제 채점기.
- **대회 규정·데이터:** VCC 2026 규정상 결과 공개 가능 여부를 제출 전 확인. 다른 팀은 이름 없이 집계로만 언급.
  데이터는 모두 공개·학술 비영리 사용 가능, 인용 의무 준수. ELEM 데이터 미사용.
- **재현성:** 2026-10-11 정리 때 입력이 지워졌다 복원된 이력이 있으므로 입력 목록을 sha256으로 고정하고 외장
  드라이브에 백업, 코드와 함께 공개.

## Feasibility

- **기간:** 2026-11-05~11-28(약 3.5주), 단독 저자.
  - 1주(11-05~11-11): Molina & Zhang 원문 확인, 이론(명제 1·2) 정리, 입력 고정, Tahoe 곡선 재실행·사전 예측.
  - 2주(11-12~11-18): 유전자 쪽 d·σ²·R² 추정과 C*, 분산 성분·상한(RQ1), RQ3, RQ4, 구조 비교 seed 재실행.
  - 3주(11-19~11-25): 초안.
  - 4주(11-26~12-01): 자체 검토·수정·제출.
- **도구:** 기존 하네스, ERA 오프라인 평가기, cell-eval2(CPU). GPU 불필요.
- **외부 의존:** arXiv 첫 제출 보증 — 1주차에 요청, 11-21까지 없으면 bioRxiv.
- **지연 시 축소 순서:** RQ4 관련 세포주 사례 → 구조 비교 seed 재실행(seed 0 단서로 보고) → RQ3를 H1 + 3개
  패널로 축소. RQ2와 RQ1은 줄이지 않는다.

## Scope Boundaries (what is NOT investigated)

- **조사하지 않음:** 새 딥러닝 구조 개발, 딥러닝 모델 플러그인 검증(PIE/STATE, 옛 E5), drug→gene 신규 실험,
  목적지 라벨을 쓰는 few-shot·미세조정, 약물 교란 예측 자체의 개선, 벤치마크 지표 비판(별도 논문; 본 논문은 부록에
  평가 시 주의점만), 비공개 데이터.

## Next Step

writing-research-plan으로 넘겨 `input/ideas/zero-shot-limits-context-complexity-plan.md`를 작성.
