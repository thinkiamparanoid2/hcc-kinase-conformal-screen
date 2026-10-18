# COMMANDBOOK — Pipeline Execution Reference

**Project:** Domain-Calibrated Machine Learning & Structure-Based Virtual Screening of Natural Products Against Hepatocellular Carcinoma Kinases (`EGFR`, `FGFR4`, `VEGFR2`)

## Environment
```powershell
Repo:   hcc-kinase-conformal-screen
Python: F:\miniconda3\envs\hcc-screen\python.exe (or conda activate hcc-screen)
```

## Run Full Pipeline (With Checkpoint/Resume)
```powershell
& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/run_pipeline.py
```
- **Resume from specific step:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/run_pipeline.py --from-step 3`
- **Reset all checkpoints:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/run_pipeline.py --reset`

## Run Steps Individually (`scripts/01` – `07`)
1. **Fetch & curate ChEMBL IC50 data:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/01_fetch_chembl_data.py`
2. **Generate training features (21 descriptors + 2048-bit ECFP4):** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/02_feature_engineering.py`
3. **10-seed Murcko-scaffold benchmark:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/03_scaffold_benchmark.py`
4. **Final XGBoost models, Mondrian conformal AD & NP domain shift:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/04_conformal_ad_models.py`
5. **Full 738k COCONUT natural-product screen:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/05_screen_coconut_library.py`
6. **Orthogonal structure-based docking validation (AutoDock Vina 1.2.5):** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/06_docking_validation.py`
7. **Auto-generated report, PAINS/QED triage & figures:** `& "F:\miniconda3\envs\hcc-screen\python.exe" scripts/07_generate_report.py`

## Key Outputs
- `results/metrics/benchmark_long.csv`, `benchmark_summary.csv`, `benchmark_paired.json`
- `results/metrics/final_model_<TARGET>.json`, `<TARGET>_test_predictions.csv`
- `models/<TARGET>_xgb_final.json`, `<TARGET>_conformal.json`
- `results/predictions/<TARGET>_screen_full.csv.gz` (all 738,823 valid compounds with `pred`, 90% conformal `[lo, hi]`, `sim_train`, `in_AD`, `rediscovery`)
- `results/predictions/consensus_hits.csv`, `consensus_hits_annotated.csv`, `consensus_hits_docked.csv`
- `results/metrics/docking_validation.json` & `results/docking_poses/<TARGET>/*_docked.pdbqt`
- `results/reports/results_report.md`, `results/figures/fig1_benchmark.png`, `fig2_domain_shift.png`, `fig3_docking_validation.png`
