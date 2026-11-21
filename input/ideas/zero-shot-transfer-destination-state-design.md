---
title: "Is the error of zero-shot cross-line perturbation transfer concentrated in the destination cell's own control-state subspace, and can a destination-state-aware correction exploit it better than existing corrections?"
type: design
created: 2026-10-31
updated: 2026-10-31
tags: [virtual-cell, perturbation-prediction, zero-shot-transfer, perturb-seq, single-cell, evaluation]
status: superseded   # by [[zero-shot-limits-context-complexity-design]] (2026-11-05)
author: mixed
---

# 처음 보는 세포주로의 교란 반응 전이: 오차는 목적지 세포의 control 상태 공간에 있는가?

(가제, 영문 논문 제목 후보: *Transfer, calibrate, deflate: destination-state-aware zero-shot prediction of
perturbation responses in unseen cell lines*)

## Research Question

**주 질문 (반증 가능):** 공개 Perturb-seq 스크린을 다른 세포주에서 옮겨 와 처음 보는 세포주의 CRISPRi 반응을
예측할 때, 전이가 놓친 오차(정답 − 예측)는 목적지 세포 **자신의** unperturbed control 세포의 상위 주성분 공간에
집중되는가? 그리고 그 공간을 이용한 보정(목적지 control-PC deflation)은 다른 세포주의 control PC, 분산 스펙트럼을
맞춘 무작위 부분공간, Systema류 보정(교란 전체 평균 제거)보다 판별력(PDS)을 더 높이는가?

- **반증 조건:** 반쪽 분할로 측정 잡음을 통제한 뒤 목적지 control PC 위의 잔차 에너지 비율이 다른 세포주 PC나
  무작위 부분공간과 통계적으로 다르지 않거나, 목적지 PC 보정이 대안 보정보다 PDS를 더 올리지 못하면 주 가설은 기각.

**보조 질문:**
1. 소스 세포 맥락 수(k)에 따른 성능 곡선은 무작위 부분집합으로 통제했을 때 아직 포화되지 않았는가
   (k = K에서 맥락 하나당 한계 가치의 신뢰구간이 0보다 큰가)?
2. 다른 랩의 같은 세포주 스크린 하나의 가치는 같은 지표로 잰 모델링 이득(구조 변경, M2, M3)과 맥락 하나 추가의
   한계 가치보다 큰가?
3. (플러그인 검증) M2/M3를 딥러닝 모델(STATE, 가능하면 PIE)의 zero-shot 출력에 후처리로 적용하면 그 모델들도
   개선되는가?

## Relevance / Forschungsstand (brief)

처음 보는 세포주로의 zero-shot 예측에서 학습형 모델이 교차 세포주 평균 전이를 이기지 못한다는 점은 이미 여러
연구가 보였다(AdaPert, ICML 2026, Appendix Table 13; Pan et al. 2026; STATE 2025; Molina & Zhang 2026). Molina &
Zhang(2026)은 세포주 고유 패턴과 그 잔차가 해당 세포주의 control PC에 몰리지만 세포주×교란 상호작용은 control만으로
복원되지 않는다고 보고했고, CIPHER(2025/2026)는 control 공분산이 세포 유형을 넘어 전이되지 않음을 보였다.
Dibaeinia et al.(2026)은 관찰 분석으로 "모델 크기보다 맥락 다양성"을 주장했다. **비어 있는 부분:** (a) 전이 예측이
놓친 오차가 다른 세포주의 축보다 목적지 자신의 축에 몰리는지에 대한 잡음 통제 비교, (b) 그 구조를 이용한 보정과
기존 보정의 직접 비교, (c) 무작위 부분집합으로 통제한 소스 맥락 수 스케일링 곡선, (d) 같은 세포주 독립 스크린의
가치를 모델링 이득과 같은 지표로 비교한 결과. 공식 대회 채점기(Virtual Cell Challenge 2026, cell-eval2)로 검증된
방법 성과(검증 단계 1,635팀 중 129위, 상위 8%)가 실증 근거를 더한다.

## Method

**제안 방법 (논문의 방법론 기여):**
- **M1 다중 소스 전이:** 공개 스크린별 표적 fold change(원시 count 합 기반, pseudocount 50 CPM, 소스별 전체 평균 반응
  제거)를 소스 가중 융합하고 하이퍼파라미터 앙상블(min_cells 7개 + pseudocount 4개)로 평균. 분석은 해석이 쉬운 융합
  기반 버전으로 하고, 대회 전체 시스템(학습형 S4 부분 포함) 결과는 E4에서만 보고.
- **M2 일치도 기반 크기 보정:** 표적별 소스 간 일치도 순위의 거듭제곱(β = 0.5)으로 반응 크기를 조정하고, 예측의 표적 간
  분산 총량은 채점 지표의 잡음 모형에서 유도되는 최적점(MSE kink)에 맞춤.
- **M3 목적지 상태 보정:** 신뢰도(M2의 일치도 순위) 하위 50 % 표적의 예측에서 목적지 control 세포(log1p CP10K, 6,000셀)의 상위 k = 8개
  주성분 사영의 η = 0.75배를 제거.

**실험:**
- **E1 (핵심) 전이 오차의 위치:** 하네스 패널(H1, OR_HCT116X, KOLF, RN K562·HepG2·Jurkat·RPE1)마다 잔차의 에너지
  비율을 ① 목적지 control PC(k = 8, 24, 50), ② 다른 세포주들의 control PC, ③ 분산 스펙트럼을 맞춘 무작위 부분공간에서
  비교. 측정 잡음은 정답 세포를 반쪽으로 나눈 교차곱으로, PC 추정 잡음은 control 세포 반쪽 분할로 통제. 표적 강도
  (사분위)별로 층화. 이어서 M3를 다른 세포주 PC 보정, 무작위 부분공간 보정, Systema류 보정(교란 평균/공통 반응 제거),
  보정 없음과 PDS·MSE로 비교(예산 불변 실험은 정확 오프라인 평가기, 핵심 결과는 H1에서 실제 채점기로 확인).
- **E2 소스 맥락 수 스케일링:** 패널마다 소스 패밀리의 무작위 부분집합(k = 1…K, 부분집합 반복)으로 융합 파이프라인을
  재구성해 PDS·MSE 곡선을 그리고, 로그 선형 모형과 포화 모형(AIC)을 비교, 맥락 하나당 한계 가치를 부트스트랩
  신뢰구간과 함께 보고. leave-one-family-out 한계 가치도 함께.
- **E3 같은 세포주 데이터의 가치:** H1 패널에 다른 랩 H1 스크린(H1GW, GSE295214)을 섞는 비중 w를 0…1로 바꾼 곡선,
  관련 세포주 사례(KOLF←H1GW, OR_HCT116←DLD-1). 이를 같은 지표로 잰 모델링 이득(과거 17개 구조 비교, M2, M3)과
  E2의 맥락 하나당 한계 가치와 비교.
- **E4 공식 검증:** VCC 2026 검증 단계 성적과, 공식 제출로 측정한 구성요소별 이득(M2 해당 +0.0106, M3 포함 묶음
  +0.0070 등). 결선 결과가 제출 전에 발표되면 추가.
- **E5 플러그인 검증 (딥러닝 모델):** 로컬 STATE 체크포인트(ST-SE-Replogle + SE-600M 임베딩)로 하네스 패널의
  zero-shot 예측을 만들고, 그 pseudobulk 델타에 M2/M3를 후처리로 적용해 개선 여부를 E1과 같은 지표로 평가. PIE
  (Arc, 공개 체크포인트)는 자원이 허락할 때만 추가.

**왜 이 접근인가:** 학습형 모델이 평균 전이를 못 이긴다는 점은 이미 알려져 있으므로, 새 구조를 더 만드는 대신
(i) 전이의 오차가 어디에 남는지 통제 분석으로 밝히고, (ii) 그 진단에서 나온 원리 있는 보정이 기존 보정보다 낫고
딥러닝 모델에도 붙는지 보이는 것이 박사 지원용 ML 방법론 논문으로 가장 방어 가능하다.

## Data Sources / Material Basis

- **소스 스크린:** Replogle 2022 K562 GWPS, Nadig 2025 HepG2·Jurkat·RPE1, X-Atlas/Orion HCT116·HEK293T, CD4 T세포
  genome-scale CRISPRi, KOLF2.1J iPSC, HipSci(Feng), VCC2025 H1, H1GW(Pan et al. 2026, GSE295214), DLD-1, DepMap 23Q4.
- **하네스 패널 (정답이 있는 목적지):** H1(VCC2025, 10x Flex, 300표적), OR_HCT116X(263표적), KOLF(279표적),
  RN K562·HepG2·Jurkat·RPE1(각 300표적). 정확 오프라인 평가 팩(`vcc2026_data/reports/era_2026-10-07/packs`).
- **공식 결과:** 우리 공식 제출 36회의 점수 기록(`vcc2026_data/reports/official_scores.csv`), 공개 리더보드 스냅샷 10회.
- **채점기:** cell-eval2 `vcc2026` preset(공개).
- **딥러닝 모델:** STATE(`checkpoints/state/ST-SE-Replogle`, `SE-600M`), 선택적으로 PIE(arcinstitute/PIE_*).
- **기존 측정 결과(재현·정리 대상):** 맥락 수 곡선(H1 PDS ≈ 0.796 + 0.048 ln k), H1GW 가치(+0.074 합성), 강한 표적의
  control-PC 에너지 비(H1 PC24 0.188 대 0.070), M3 하네스 7/7 패널 양수.

## Output Target

☐ Book chapter   ☐ Journal article   ☐ Grant proposal   ☐ Talk
☒ Other: **arXiv 프리프린트**(cs.LG 주분류, q-bio.GN 교차 게시; 보증 확보 실패 시 bioRxiv), 본문 8~10쪽 + 부록,
코드 공개(GitHub). 목표 제출일 2026-11-28~12-01(박사 지원 12월 마감용). 워크숍은 선택(v2 후보: ICLR 2027 워크숍,
권장 마감 2027-02-01).

## Expected Results

- **의미 있는 긍정 결과:** 잡음 통제 후에도 전이 잔차의 목적지 control PC 에너지 비율이 다른 세포주 PC와 무작위
  부분공간보다 유의하게 크고, 강한 표적에서 그 차이가 더 크다. M3가 대안 보정들보다 PDS를 더 올린다. 맥락 수 곡선이
  K에서 포화되지 않는다. 같은 세포주 데이터의 가치가 맥락 하나 추가와 최선의 모델링 이득보다 한 자릿수 이상 크다.
  M2/M3가 STATE 출력도 개선한다.
- **귀무 결과의 해석:** 잡음 통제 후 집중이 사라지면 "겉보기 집중은 표본 잡음"으로 보고하고 M3의 효과는 잡음/허브
  억제로 재해석한다. 스케일링이 포화되면 "현재 공개 소스의 다양성은 소진됐다"고 보고한다. STATE에서 개선이 없으면
  M2/M3가 소스 융합 예측에 특화된 보정임을 밝힌다. 어떤 경우든 통제 분석 자체는 보고할 가치가 있다.

## Ethical / Methodological Risks

- **선행 연구와의 겹침:** Molina & Zhang 2026, CIPHER, Dibaeinia 2026, AdaPert 2026, Systema. 서론·관련 연구에서 차이
  ("전이 잔차" 대상, 다른 세포주 축 귀무, 보정 비교, 공식 벤치마크)를 명시.
- **측정 잡음:** 얕은 패널(OR 163셀, RN 53~96셀)은 DE 지표가 불안정. H1을 주 패널로, 반쪽 분할로 통제. 잡음의
  공분산(≈ Σ_ctrl/n)이 control PC 위에 놓이는 구조적 편향을 반드시 분리.
- **오프라인 평가기의 범위:** 예산 불변 실험에서만 정확. 예산이 바뀌는 실험은 실제 채점기로만 판정(회당 약 25분, 약
  15 GB, 한 번에 하나).
- **대회 규정·데이터:** VCC 2026 규정상 결과 공개·공식 데이터 사용 가능 여부를 제출 전 확인. 논문의 공식 근거는 공개
  점수에 한정. 다른 팀은 이름 없이 집계로만 언급.
- **데이터 라이선스:** 사용 소스 모두 학술 비영리 사용 가능, 인용 의무 준수. ELEM 데이터는 사용하지 않음.
- **재현성:** 2026-10-11 정리 때 지워진 입력을 복원한 이력이 있으므로 필요한 입력 목록을 고정하고 코드와 함께 공개.

## Feasibility

- **기간:** 4주(2026-11-01~11-28), 단독 저자. 1주 E1, 2주 E2·E3·E4·E5(STATE 추론), 3주 초안, 4주 자체 검토·수정·제출.
- **도구:** 기존 하네스, ERA 정확 오프라인 평가기, cell-eval2(CPU). E5만 GPU 몇 시간(다른 프로젝트와 공유하므로 빈
  시간에 짧게 실행; 불가하면 CPU 추론 또는 E5 축소).
- **외부 의존:** arXiv 첫 제출 보증(endorsement) — 1주차에 확보 시도, 실패 시 bioRxiv.
- **지연 시 축소 순서:** E3 관련 세포주 사례 → E5의 PIE → E2의 반복 횟수.

## Scope Boundaries (what is NOT investigated)

- **조사하지 않음:** 새 딥러닝 구조 개발, 목적지 라벨을 쓰는 few-shot·미세조정, 약물 교란, 복제 데이터 기반의 본격적인
  반응 분해(v2로 미룸), 벤치마크 지표 비판(별도 논문 B; 본 논문은 부록에 평가 시 주의점만), 비공개 데이터.

## Next Step

writing-research-plan으로 넘겨 `input/ideas/zero-shot-transfer-destination-state-plan.md`를 작성.
