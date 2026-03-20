# Research Log: DGM Alignment Tuning

## 1) 목표
- 데이터: `DGM_Valid_CAO.csv`의 관측 `DGM`
- 모델: `Unit_WorM3_CAO_All_Calib_Direct.py`
- 이슈: salinity-DGM 관계에서 모델 slope가 관측 slope와 불일치
- 조정 대상: `ko_UVA` 계수, `RHg_UVA` 분자식 계수

## 2) 코드 변경 진척
- `Unit_WorM3_CAO_All_Calib_Direct.py`에 튜닝 포인트 주석 추가
  - `RHg_UVA` 분자식 계수 (`rhg_uva_a`, `rhg_uva_b`, `rhg_uva_c`)
  - `ko_UVA` 계수 (`ko_uva_num`)
- 자동 최적화 함수 `optimize_uva_coefficients(...)` 추가
  - 목적함수: `RMSE + slope_weight * |slope_sim(sal) - slope_obs(sal)|`
- CLI 추가: `--tune-uva`, `--n-starts`, `--n-iters`, `--slope-weight`

## 3) uv run 기반 최적화 실행
실행 방식: multi-seed 탐색 (seed: 42, 7, 123), `n_starts=20`, `n_iters=120`, `slope_weight=2.0`

| seed | ko_uva_num | rhg_uva_a | rhg_uva_b | rhg_uva_c | RMSE | objective |
|---|---:|---:|---:|---:|---:|---:|
| 42 | 0.35040701 | -0.02243776 | 1.49119117 | -23.61252119 | 0.02214649 | 0.02214649 |
| 7 | 0.75518161 | -0.02620034 | 2.09406068 | -36.73345321 | 0.02254645 | 0.02254645 |
| 123 | 0.74158509 | -0.02879055 | 2.22536599 | -38.43784032 | 0.02242088 | 0.02242088 |

선정된 최적 계수 (best seed=42):
- `ko_uva_num = 0.350407009078024`
- `rhg_uva_a = -0.02243776483632441`
- `rhg_uva_b = 1.4911911655729884`
- `rhg_uva_c = -23.612521190306616`

## 4) 성능 비교 (Baseline vs Tuned)
- Baseline: 기존 계수
- Tuned: 위 최적 계수 적용

| Metric | Baseline | Tuned | Delta (Tuned - Baseline) |
|---|---:|---:|---:|
| R2_reg | 0.73023 | 0.78166 | +0.05143 |
| RMSE | 0.02667 | 0.02215 | -0.00453 |
| Slope (obs->sim regression) | 0.66718 | 0.84918 | +0.18199 |
| Intercept | 0.02664 | 0.01310 | -0.01354 |
| PBIAS (%) | 9.98269 | 3.62438 | -6.35831 |
| MeanBias | -0.01142 | -0.00414 | +0.00727 |
| Salinity slope gap (model-obs) | -0.00839037 | ~0.00000000 | +0.00839037 |

해석:
- salinity-DGM slope mismatch가 사실상 0으로 보정됨.
- 전체 오차(RMSE)와 편향(PBIAS)도 동반 개선.

## 5) 노트북 산출물
생성 파일: `DGM_tuned_analysis.ipynb`

노트북 내용:
- 최적 계수 적용 모델 재실행
- Baseline vs Tuned 통계 계산: `R2`, `RMSE`, `r`, `slope`, `intercept`, `PBIAS`, `MeanBias`, salinity slope gap
- 변수별 비교 그래프 (관측 DGM vs tuned DGM):
  - Salinity, Temperature, SeaIce, WindSpeed
- 관측값-모델값 관계 그래프:
  - Baseline/Tuned 산점도 + 1:1선 + 회귀선

노트북/분석 코드에서 생성된 그림:
- `fig_tuned_variable_comparison.png`
- `fig_obs_vs_model_baseline_vs_tuned.png`

## 6) 재실행 명령
```bash
uv run python Unit_WorM3_CAO_All_Calib_Direct.py \
  --tune-uva \
  --input-csv DGM_Valid_CAO.csv \
  --n-starts 20 \
  --n-iters 120 \
  --slope-weight 2.0 \
  --output-summary-csv Hg_budget_summary_tuned_best.csv \
  --output-metrics-csv Hg_model_metrics_tuned_best.csv \
  --output-tuning-csv Hg_uva_tuning_result.csv
```
