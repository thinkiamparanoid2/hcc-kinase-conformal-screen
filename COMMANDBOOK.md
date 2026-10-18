# COMMANDBOOK — pipeline execution reference

**Project:** Domain-calibrated machine learning and structure-based virtual
screening of natural products against hepatocellular carcinoma kinases
(`EGFR`, `FGFR4`, `VEGFR2`).

All commands are platform-independent. Use `python` from the activated
environment; no interpreter path is hard-coded anywhere in the pipeline.

## Environment

```bash
conda env create -f environment.yml && conda activate hcc-screen   # includes vina, obabel
# or
pip install -r requirements.txt                                     # steps 01-05, 07 only
```

External binaries are found on `PATH`. Override with environment variables when
they live elsewhere:

```bash
export VINA_BIN=/path/to/vina
export OBABEL_BIN=/path/to/obabel
```

## Full pipeline

```bash
python scripts/run_pipeline.py              # run, resuming from the last checkpoint
python scripts/run_pipeline.py --list       # show steps and their prerequisites
python scripts/run_pipeline.py --from-step 3
python scripts/run_pipeline.py --only 4 7
python scripts/run_pipeline.py --reset      # clear checkpoints
```

A step whose prerequisites are missing is skipped with the reason printed, rather
than crashing partway through.

## Individual steps

| Step | Command | Needs |
|---|---|---|
| 1 | `python scripts/01_fetch_chembl_data.py` | ChEMBL network access |
| 1b | `python scripts/01_fetch_chembl_data.py --from-existing` | re-curate local CSVs offline |
| 2 | `python scripts/02_feature_engineering.py` | step 1 output |
| 3 | `python scripts/03_scaffold_benchmark.py` | step 2 output |
| 4 | `python scripts/04_conformal_ad_models.py` | steps 2-3 output |
| 5 | `python scripts/05_screen_coconut_library.py` | step 4 + COCONUT library |
| 6 | `python scripts/06_docking_validation.py` | step 5 + `vina`, `obabel` |
| 7 | `python scripts/07_generate_report.py` | whatever exists; missing sections are marked |

Smoke test a single target and seed:

```bash
python scripts/03_scaffold_benchmark.py FGFR4 0
python scripts/04_conformal_ad_models.py FGFR4
```

## Tests

```bash
pytest -q
```

## Key outputs

- `results/metrics/curation_<TARGET>.json` — curation counts, censoring, label noise
- `results/metrics/benchmark_canonical.json` — canonical split with bootstrap CIs
- `results/metrics/benchmark_temporal.json` — temporal holdout
- `results/metrics/model_selection.json` — the family steps 04-05 use, and why
- `results/metrics/final_model_<TARGET>.json` — conformal bins, `tau` + CI, NP transfer
- `models/<TARGET>_final.{json,pkl}`, `models/<TARGET>_conformal.json`
- `results/predictions/<TARGET>_screen_full.csv.gz` — full screen with novelty flags
- `results/predictions/consensus_hits{,_annotated,_docked}.csv`
- `results/metrics/docking_validation.json` — redocking, enrichment, candidates
- `results/reports/results_report.md`, `results/figures/fig{1,2,3}_*.png`
