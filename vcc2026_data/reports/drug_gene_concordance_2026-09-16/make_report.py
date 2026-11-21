from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent
s=pd.read_csv(O/'summary.csv');r=pd.read_csv(O/'repeat_summary.csv');a=pd.read_csv(O/'pair_metrics.csv')
order=['lincs_full_local_CGS','jiang_IFNB','jiang_INS','jiang_TGFB','jiang_TNFA'];labels=['LINCS shRNA','Jiang IFNB','Jiang INS','Jiang TGFB','Jiang TNFA']
fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
y=np.arange(5)
for mode,offset,col,label in [('raw',-.13,'#2365a1','Raw effects'),('generic_removed',.13,'#d87a28','Common response removed')]:
 d=s[(s.subset=='all_annotated_inhibitors')&(s['mode']==mode)].set_index('dataset').loc[order]
 for ax,key in zip(axes,['spearman','spearman_advantage']):
  v=d[key].to_numpy();lo=d[key+'_lo'].to_numpy();hi=d[key+'_hi'].to_numpy()
  ax.errorbar(v,y+offset,xerr=np.stack([v-lo,hi-v]),fmt='o',capsize=3,color=col,label=label)
for ax,title in zip(axes,['Matched gene-drug rank correlation','Matched minus wrong-target correlation']):
 ax.axvline(0,color='#777',lw=1);ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_title(title,fontsize=11);ax.grid(axis='x',alpha=.15);ax.set_xlabel('Spearman correlation / difference')
axes[0].legend(fontsize=8,loc='upper right');fig.suptitle('Same named cell line, different assays and conditions',fontsize=13)
fig.savefig(O/'concordance.png',dpi=170);plt.close(fig)
raw=s[(s.subset=='all_annotated_inhibitors')&(s['mode']=='raw')].set_index('dataset')
center=s[(s.subset=='all_annotated_inhibitors')&(s['mode']=='generic_removed')].set_index('dataset')
lines=['# 같은 표적의 유전자 억제–약물 반응 일치도: 탐색적 실측 결과','',
'2026-09-16. **관측된 전사체 반응의 일치도는 평균적으로 매우 낮다.** 같은 표적을 연결하면 잘못 연결한 표적보다 조금 더 비슷한 경향은 있지만, 현재 데이터에서 두 반응을 서로 대체 가능한 것으로 볼 근거는 없다. 이 결과는 실험 조건·측정 방식·잡음이 섞인 관측 비교이며, 순수한 처리 방식의 차이나 공통 생물학의 부재를 입증하지 않는다.','',
'## 범위와 해석','',
'A549와 HT-29 두 세포주. Tahoe 150약물 패널 중 억제제/길항제 주석이 있는 후보를 사용했다. Jiang CRISPRi 110조건(자극별), LINCS shRNA 291쌍, 합계 401조건을 계산했다. 같은 표적/약물/세포가 여러 조건에 재등장하므로 401개 독립 실험이 아니다. 기존 후보 선정은 반응값을 보지 않고 수행됐으며, 이번 순위 상관을 보고 표적 주석을 수정하거나 후보를 선별하지 않았다. 출판사 표적 주석은 별도 약리학적 검증을 완료하지 않은 상태다.','',
'주 지표는 공통 측정 유전자에서의 **Spearman 순위 상관**이다. 1은 변화량 순위가 완전히 같음, 0은 단조 순위 관계가 없음, -1은 순위가 완전히 반대임을 뜻한다. 상관 0.018을 1.8% 일치 또는 98.2% 차이라고 해석하면 안 된다. 원래 측정 단위가 다른 LINCS/Tahoe의 MSE나 절대 크기 차이는 보고하지 않았다.','',
'평균은 약물–세포 내 복수 표적을 먼저 평균한 뒤, 약물별 세포주를 평균하고, 마지막으로 약물에 같은 가중치를 부여했다. 각 쌍의 잘못된 표적 대조군은 같은 데이터셋·세포주의 후보 표적 중 그 약물의 **모든 주석 표적을 제외**한 것이다. 공유 경로가 있을 수 있어 생물학적으로 완전한 음성 대조군이라는 뜻은 아니다.','',
'## 원래 반응에서의 결과','',
'| 비교 | 약물 수 | RNA 좌표 수 | 같은 표적 상관 | 잘못된 표적 평균 | 차이 | 같은 표적의 상대 순위 |','|---|---:|---:|---:|---:|---:|---:|']
for name,label in zip(order,labels):
 v=raw.loc[name];lines.append(f'| {label} | {int(v.n_drugs)} | {int(v.n_genes):,} | {v.spearman:.4f} | {v.wrong_spearman:.4f} | {v.spearman_advantage:+.4f} | {100*v.matched_percentile:.1f}% |')
lines += ['', '상대 순위는 가능한 잘못된 표적과 비교해 더 높은 상관을 얻는 비율의 평균(동률 절반 반영)이다. 50%가 교환 가능한 무작위 대응의 기준이다. 공식 PDS 점수가 아니다. Jiang IFNB에는 후보 표적이 4개뿐이라 개별 상대 순위가 매우 거칠다.', '',
'LINCS의 원래 상관은 0.01775이며 약물 부트스트랩 기술 구간은 0.00667–0.03024이다. 잘못된 표적 대비 차이는 +0.01080, 구간은 +0.00310–+0.01936이다. 효과는 작고, 두 세포주만 관측했으며, 동일 표적을 공유하는 약물 간 의존성까지 반영한 구간은 아니다. 확증적 유의성이나 새로운 세포주에 대한 보장으로 해석하지 않는다.', '',
'## 공통 처리 반응을 제거한 민감도 분석','',
'각 자료·세포주에서 모든 사용 가능한 perturbation의 평균 반응을 뺐다. Tahoe는 150개 약물, Jiang은 각 자극 패널 전체 표적, LINCS는 해당 세포주의 96시간 shRNA CGS 전체 표적을 사용했다. 같은 표적끼리만 평균을 낸 것은 아니다. 이는 관측자료의 서술적 변환이며 학습·예측 성능이 아니다.','',
'| 비교 | 같은 표적 상관 | 잘못된 표적 대비 차이 | 차이의 기술적 95% 구간 |','|---|---:|---:|---|']
for name,label in zip(order,labels):
 v=center.loc[name];lines.append(f'| {label} | {v.spearman:.4f} | {v.spearman_advantage:+.4f} | {v.spearman_advantage_lo:+.4f} ~ {v.spearman_advantage_hi:+.4f} |')
lines += ['', 'LINCS의 작은 표적 대응 이점과 Jiang TNFA의 상대적으로 큰 이점은 이 변환 후에도 남는다. IFNB/INS/TGFB는 차이의 기술 구간에 0을 포함한다. 여러 탐색적 비교와 적은 표적 수를 감안해야 한다.', '',
'![같은 표적 상관과 대조군 차이](concordance.png)', '',
'## 반복 실험 확인','',
'Jiang의 동일 유전자 억제 두 반복 간 상관은 후보쌍 가중 기준으로 자극별 약 0.069–0.115에 불과했다. 원래 실험 자체의 전사체 순위 재현성도 낮으므로, 유전자–약물 상관이 낮은 원인을 생물학적 차이 하나로 돌릴 수 없다. 잡음 보정 상관이나 설명 가능한 생물학 비율은 계산하지 않았다.','',
'| 비교 | Tahoe 반복이 있는 약물 수 | 같은 약물의 반복 간 상관 | 유전자–약물 1차 상관 | 유전자–약물 반복 상관 | 1차의 잘못된 표적 대비 차이 | 반복의 잘못된 표적 대비 차이 |','|---|---:|---:|---:|---:|---:|---:|']
for name,label in zip(order,labels):
 v=r[(r.dataset==name)&(r.repeat_kind=='drug')].iloc[0]
 lines.append(f'| {label} | {int(v.n_drugs)} | {v.drug_repeat_rho:.4f} | {v.cross_primary:.4f} | {v.cross_drug_repeat:.4f} | {v.cross_primary_advantage:+.4f} | {v.cross_drug_repeat_advantage:+.4f} |')
lines += ['', '각 행은 1차/반복 모두 존재하는 **동일 약물 부분집합**끼리 비교한다. 위 전체 결과와 직접 비교해 반복 때문에 감소했다고 말하면 안 된다. 특히 LINCS 반복 하위 집합은 26약물/82쌍이다. LINCS 자체의 독립적인 생물학적 반복은 이번 분석에서 확보되지 않았다. Jiang의 Tahoe 반복 하위 집합은 2–3개 약물뿐이라 전체 결론의 검증으로 쓰기에는 작다.', '',
'## 단일 표적 주석만 남긴 결과','',
'여러 표적이 적힌 약물을 제외해도 LINCS 원래 상관은 29약물에서 0.02128, 잘못된 표적 대비 +0.01928이다. 공통 반응 제거 후 상관 0.01689, 차이 +0.01619로 여전히 작다. 단일 표적 주석이 약물의 실제 선택성을 보증하지는 않는다. Jiang 민감도 분석은 자극별 3–6개 약물에 불과하며 IFNB에는 해당 약물이 없다.', '',
'## 제한과 다음 판단','',
'- 같은 **이름의 세포주**이며 동일 클론·배양·실험실 조건이 아니다.',
'- Tahoe 약물 처리 24시간과 LINCS shRNA 96시간은 다르다. Jiang은 별도 신호 자극 조건이며 자극 24시간을 유전자 억제 기간과 혼동하지 않았다.',
'- CRISPRi, shRNA 및 약물 억제를 동일한 처리로 취급하지 않는다. 약물의 다중 표적·변이 선택성·주석 오류/누락을 이 첫 계산으로 해결하지 못했다.',
'- 유전자 기호 중복은 자료 내 산술평균으로 합쳤고, 후보 억제제의 전체 주석 표적을 RNA 평가 좌표에서 제외했다. 공통 945개 LINCS 좌표에서 922개가 남았고, Jiang은 16,554–16,621개다. 기존 Jiang 패널의 표적 좌표 제외도 유지했다.',
'- 같은 표적 연결의 평균 이점은 일부 존재하나, 약물 반응을 유전자 억제 반응의 정답처럼 직접 대체할 근거는 부족하다. 이 결과로 약물 사전학습 자체의 효용을 판정할 수는 없다.',
'- 다음 확증 단계는 표적 주석을 독립 검토한 소수의 명확한 후보에서 시간·자극·배양 조건을 더 가깝게 맞추는 것이다. 현재 수치는 후속 가설을 만드는 탐색 결과다.', '',
'## 재현 및 검증','',
'`ANALYSIS_PLAN.md`를 상관 계산 전에 저장했다. `analyze.py`는 입력 SHA256 보존, 잘못된 표적의 모든 약물 표적 제외, 유한값, 401개 후보 수, 실제 벡터의 scipy 상관/직접 cosine 교차검사를 수행했다. 별도 `verify_outputs.py`가 10개 데이터셋–세포 조합의 실제 입력 ID와 유전자 축을 다시 읽어 상관을 재계산하고 모든 요약표의 평균을 확인한다. `verification.json`, `independent_verification.json`에 검증 결과를 저장한다.', '',
'`pair_metrics.csv`: 401조건 × 두 변환의 지표. `summary.csv`: 약물 균등 가중 결과와 기술적 부트스트랩 구간. `by_line.csv`, `by_drug.csv`: 세포주/약물별 결과. `repeat_metrics.csv`, `repeat_summary.csv`: 반복 분석. `features.json`, `excluded_target_genes.json`: 정확한 평가 좌표. 반응값/모델/공식 제출은 수정하지 않았다.']
(O/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print('Report and figure written.')
