# Reproducing this work

Everything except two external inputs is reproducible from a clean checkout. The
two exceptions are named below with exactly what they need.

## 1. Environment

```bash
conda env create -f environment.yml
conda activate hcc-screen
pytest -q
```

The conda environment pins AutoDock Vina 1.2.5 and Open Babel 3.1.1 alongside the
Python stack, so step 06 runs without a separate install. With pip
(`pip install -r requirements.txt`) steps 01–05 and 07 run; step 06 needs the two
binaries on `PATH`, or `VINA_BIN` / `OBABEL_BIN` pointing at them.

Docker, if you would rather not touch the host environment:

```bash
docker build -t hcc-screen .
docker run --rm -v "$PWD:/work" hcc-screen python scripts/run_pipeline.py
```

## 2. External inputs

| Input | Needed by | How to get it |
|---|---|---|
| ChEMBL IC50 records | step 01 | Network access to `www.ebi.ac.uk`. Curated CSVs for the three targets are committed, so step 01 is only needed to refresh them — or to add the censored (`>`) records that populate the inactive class. |
| COCONUT library (~550 MB) | step 05 | Download "CSV Lite" from <https://coconut.naturalproducts.net/download> and place it at `data/ligands_phytochemical/phytochemical_screening_set.csv`. |

Receptor structures (`1M17`, `5JKG`, `4ASD`) are committed under `data/targets/`.

## 3. Run

```bash
python scripts/run_pipeline.py --list    # steps and prerequisites
python scripts/run_pipeline.py           # checkpointed; re-run to resume
```

Approximate cost on 8 cores: step 02 ~10 min, step 03 ~1–2 h (10 seeds x 3 targets
x 6 models), step 04 ~15 min, step 05 ~20–40 min, step 06 hours (dominated by
docking x `N_SEEDS` replicates x exhaustiveness 32), step 07 seconds.

## 4. Determinism

| Component | Reproducible? |
|---|---|
| Structure standardisation, featurisation | Yes, deterministic |
| Canonical scaffold split | Yes, seed-independent by construction |
| Repeated scaffold splits | Yes, seeded (`common.SEEDS`) |
| RF / XGBoost fits | Yes, seeded |
| Conformal quantiles, `tau` | Yes, deterministic given the split |
| Bootstrap CIs | Yes, seeded |
| AutoDock Vina | Yes per seed (`--seed 42 + i`); results are averaged over `N_SEEDS` |

Thread counts can change floating-point reduction order; differences should sit
far inside the reported confidence intervals. If they do not, that is a finding
worth reporting.

## 5. What a re-run changes

Steps 03 and 04 overwrite `results/metrics/` and `models/`. Because step 03 now
selects the model family that steps 04 and 05 use, changing the data can change
the selected model — that is intended, and `model_selection.json` records the
decision and the margin.

Step 05 reads `models/<TARGET>_conformal.json`, so a screen is only valid against
the models in the same run. Re-run 04 and 05 together after any data change.

## 6. Known gaps in a fresh reproduction

- If step 01 cannot reach ChEMBL, the committed CSVs are used as-is. They were
  fetched before censored relations were requested, so they contain no inactive
  class; `results/metrics/curation_*.json` records this
  (`censored_records_available: false`) and the report repeats it as a limitation.
- Re-curating pre-aggregated CSVs (`--from-existing`) cannot recover replicate
  spread, so the label-noise ceiling on R² is reported as unknown rather than as a
  spurious 1.0.
- Step 06's enrichment control needs an inactive class; without one it is skipped
  with that reason recorded, rather than falling back to the potency extremes that
  made the earlier control look stronger than it was.
