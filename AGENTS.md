# Repository Guidelines

## Project Structure & Module Organization
This repository is a flat, script-driven modeling workspace (no `src/` tree yet). Most logic lives in root-level MATLAB files:
- `Unit_WorM3_*.m`: core Arctic Hg model variants.
- `Unit_WorM3_ODE_*.m`: ODE definitions used by scenario scripts.
- `*_Simulation.m` and `*_Sensitivity.m`: batch experiments and perturbation analyses.
- `Example*.m`, `Exercise*.m`, `HW*.m`, `TEST.m`, `Twoelectron_test.m`: small experiments and fitting demos.
- `*.csv`: both inputs (for example `DGM_Valid_CAO.csv`, `DGM_Eval_CAO.csv`) and generated outputs (for example `Hg_budget_summary.csv`).
- `Unit_WorM3_CAO_All_Calib_Direct.py`: Python implementation of the calibration workflow.

## Build, Test, and Development Commands
No build system is configured; run scripts directly from the repository root.
- `matlab -batch "Unit_WorM3_CAO_All_Calib_Direct"`: baseline calibration run; writes `Hg_budget_summary.csv` and `Hg_model_metrics.csv`.
- `matlab -batch "Unit_WorM3_CAO_All_Calib_Direct_Simulation"`: 5x5 phase-space experiment; writes `Hg_phase_space_*`.
- `matlab -batch "Unit_WorM3_CAO_All_Calib_Sensitivity"`: 0.1% perturbation sensitivity run.
- `python Unit_WorM3_CAO_All_Calib_Direct.py`: Python equivalent of baseline run (`numpy`, `pandas` required).

## Coding Style & Naming Conventions
- MATLAB: use 4-space indentation, keep one statement per line, and preserve clear section blocks (`%%`).
- Naming pattern already used: scenario drivers as `Unit_WorM3_<Domain>_<Mode>.m`, ODE helpers as `*_ODE_*.m`.
- Prefer descriptive table column names (`Hg0_water`, `PBIAS_percent`) and explicit units in comments when adding new variables.
- Keep generated data filenames deterministic and scenario-specific (for example `Hg_phase_space_summary_5x5_delta.csv`).

## Testing Guidelines
There is no formal test harness yet. Validate changes by re-running at least one baseline and one perturbation script:
- Baseline: `Unit_WorM3_CAO_All_Calib_Direct.m` (or Python equivalent).
- Sensitivity/fit sanity: `Unit_WorM3_CAO_All_Calib_Sensitivity.m`, `Twoelectron_test.m`.
For model changes, compare key outputs (`RMSE`, `R2`, `PBIAS`) against prior CSVs and document deltas in your PR.

## Commit & Pull Request Guidelines
Git history is currently minimal and does not enforce a strict convention. Use:
- Commit format: `<area>: <imperative summary>` (example: `model: adjust ice-threshold sensitivity`).
- Keep commits focused (model logic, data update, or output refresh separately).
- PRs should include: purpose, files changed, commands run, and a short before/after metrics table. Link issues when available.
