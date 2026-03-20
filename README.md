# Mercury-arctic

수은(mercury) 해양-대기 상호작용 모델링을 위한 MATLAB/Python 작업공간입니다.  
파일이 많아져서 `입력/출력/예제/분석노트` 기준으로 폴더를 재정리했습니다.

## 폴더 구조
- `inputs/`: 모델 입력 데이터 CSV
- `outputs/`: 모델 실행 결과 CSV/JSON/PNG
- `examples/`: 실험/연습용 MATLAB 스크립트
- `notebooks/`: 분석 노트북
- `notes/`: 연구 로그/메모
- 루트(`./`): 핵심 실행 스크립트(`Unit_WorM3_*`, `Unit_WorM3_ODE_*`, Python 엔트리 등)

## 실행

### Python (기본)
```bash
uv run python Unit_WorM3_CAO_All_Calib_Direct.py
```
- 기본 입력: `inputs/DGM_Valid_CAO.csv`
- 기본 출력: `outputs/Hg_budget_summary.csv`, `outputs/Hg_model_metrics.csv`

### Python (튜닝)
```bash
uv run python Unit_WorM3_CAO_All_Calib_Direct.py \
  --tune-uva \
  --n-starts 20 \
  --n-iters 120 \
  --slope-weight 2.0
```
- 튜닝 결과 출력: `outputs/Hg_uva_tuning_result.csv`

### MATLAB
핵심 MATLAB 스크립트는 루트에 유지됩니다.  
주요 드라이버는 입력을 `inputs/`에서 읽고 결과를 `outputs/`에 저장하도록 경로를 맞춰 두었습니다.
