# HCC kinase natural-product screen: results report

_Auto-generated from `results/metrics` and `results/predictions`. Every number, table and caveat below is read from those files; nothing is hand-typed._

## 0. Data curation

Structures are salt-stripped, neutralised and tautomer-canonicalised, then deduplicated on the parent InChIKey. Censored measurements are labelled, not silently treated as point values: only uncensored rows are used as regression targets.

| target | final compounds | exact labels | censored (potent bound) | censored (inactive bound) | active prevalence (exact) | R2 ceiling from label noise |
|---|---|---|---|---|---|---|
| FGFR4 | 3608 | 3608 | 0 | 0 | 0.763 | - |
| EGFR | 9782 | 9782 | 0 | 0 | 0.749 | - |
| VEGFR2 | 9181 | 9181 | 0 | 0 | 0.767 | - |

**Curation caveats**

- FGFR4: Replicate spread is not recoverable from a pre-aggregated CSV, so the label-noise ceiling on R2 is unknown here. Re-run step 01 with network access to measure it from the raw ChEMBL records.
- FGFR4: the curated file contains no `>` records, so there is NO inactive class. Every classification metric below therefore describes ranking within a set that is 76% active, not active-versus-inactive discrimination. Re-run step 01 with network access to populate the inactive class.
- EGFR: Replicate spread is not recoverable from a pre-aggregated CSV, so the label-noise ceiling on R2 is unknown here. Re-run step 01 with network access to measure it from the raw ChEMBL records.
- EGFR: the curated file contains no `>` records, so there is NO inactive class. Every classification metric below therefore describes ranking within a set that is 75% active, not active-versus-inactive discrimination. Re-run step 01 with network access to populate the inactive class.
- VEGFR2: Replicate spread is not recoverable from a pre-aggregated CSV, so the label-noise ceiling on R2 is unknown here. Re-run step 01 with network access to measure it from the raw ChEMBL records.
- VEGFR2: the curated file contains no `>` records, so there is NO inactive class. Every classification metric below therefore describes ranking within a set that is 77% active, not active-versus-inactive discrimination. Re-run step 01 with network access to populate the inactive class.

## 1. Benchmark: scaffold-split evaluation

Scaffold splits are generic Bemis-Murcko and strictly disjoint. The canonical split is deterministic (largest scaffold families in train) and carries percentile bootstrap CIs over test molecules. The repeated splits measure how much the result moves with the scaffold assignment: that spread is SPLIT variance, not sampling variance, because the same molecules are reused each time.

_Hyperparameters are fixed a priori at common library defaults and were not tuned on any split reported here; no model selection was performed against test-set performance. Tuning, if added, must run as an inner cross-validation loop on the training scaffolds only._

### 1.1 Canonical deterministic scaffold split (with 95% CIs)

**FGFR4** (train 2887, test 721, test active prevalence 0.77, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.228 | [0.159, 0.293] | 1.087 | 0.485 | 0.758 | [0.717, 0.796] |
| RF-ecfp | 0.449 | [0.383, 0.510] | 0.918 | 0.660 | 0.839 | [0.801, 0.874] |
| RF-all | 0.411 | [0.347, 0.469] | 0.949 | 0.644 | 0.826 | [0.787, 0.860] |
| XGB-all | 0.435 | [0.370, 0.496] | 0.930 | 0.645 | 0.830 | [0.793, 0.867] |
| kNN-Tanimoto | 0.330 | [0.231, 0.422] | 1.013 | 0.599 | 0.801 | [0.758, 0.842] |

**EGFR** (train 7826, test 1956, test active prevalence 0.69, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.188 | [0.143, 0.228] | 1.183 | 0.484 | 0.732 | [0.707, 0.756] |
| RF-ecfp | 0.427 | [0.384, 0.468] | 0.994 | 0.668 | 0.813 | [0.792, 0.834] |
| RF-all | 0.388 | [0.344, 0.428] | 1.027 | 0.646 | 0.799 | [0.778, 0.820] |
| XGB-all | 0.383 | [0.338, 0.427] | 1.031 | 0.639 | 0.796 | [0.774, 0.817] |
| kNN-Tanimoto | 0.220 | [0.142, 0.294] | 1.159 | 0.619 | 0.782 | [0.758, 0.806] |

**VEGFR2** (train 7346, test 1835, test active prevalence 0.69, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.119 | [0.076, 0.156] | 1.097 | 0.393 | 0.676 | [0.650, 0.701] |
| RF-ecfp | 0.404 | [0.365, 0.438] | 0.902 | 0.635 | 0.800 | [0.778, 0.821] |
| RF-all | 0.398 | [0.362, 0.431] | 0.907 | 0.640 | 0.798 | [0.777, 0.820] |
| XGB-all | 0.384 | [0.345, 0.421] | 0.917 | 0.622 | 0.788 | [0.767, 0.809] |
| kNN-Tanimoto | 0.280 | [0.221, 0.336] | 0.992 | 0.596 | 0.770 | [0.747, 0.792] |

### 1.2 Repeated random scaffold splits (mean +/- SD over seeds)

**FGFR4** (10 splits, mean test active prevalence 0.75)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.026 +/- 0.039 | 1.258 +/- 0.052 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.499 +/- 0.067 | 0.877 +/- 0.061 | 0.703 +/- 0.052 | 0.872 +/- 0.033 |
| RF-desc | 0.369 +/- 0.056 | 0.985 +/- 0.046 | 0.599 +/- 0.054 | 0.829 +/- 0.028 |
| RF-ecfp | 0.566 +/- 0.053 | 0.816 +/- 0.035 | 0.740 +/- 0.041 | 0.897 +/- 0.025 |
| RF-all | 0.547 +/- 0.056 | 0.834 +/- 0.038 | 0.731 +/- 0.048 | 0.891 +/- 0.025 |
| XGB-all | 0.556 +/- 0.058 | 0.825 +/- 0.038 | 0.735 +/- 0.046 | 0.888 +/- 0.024 |
| XGB-scramble | -0.130 +/- 0.068 | 1.320 +/- 0.067 | 0.063 +/- 0.098 | 0.535 +/- 0.064 |

**EGFR** (10 splits, mean test active prevalence 0.76)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.009 +/- 0.008 | 1.348 +/- 0.022 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.531 +/- 0.034 | 0.918 +/- 0.031 | 0.733 +/- 0.029 | 0.872 +/- 0.012 |
| RF-desc | 0.247 +/- 0.124 | 1.161 +/- 0.094 | 0.502 +/- 0.105 | 0.781 +/- 0.029 |
| RF-ecfp | 0.574 +/- 0.032 | 0.876 +/- 0.037 | 0.761 +/- 0.026 | 0.893 +/- 0.009 |
| RF-all | 0.547 +/- 0.040 | 0.902 +/- 0.043 | 0.746 +/- 0.030 | 0.889 +/- 0.011 |
| XGB-all | 0.549 +/- 0.040 | 0.901 +/- 0.041 | 0.746 +/- 0.033 | 0.886 +/- 0.011 |
| XGB-scramble | -0.066 +/- 0.038 | 1.385 +/- 0.032 | 0.029 +/- 0.056 | 0.513 +/- 0.031 |

**VEGFR2** (10 splits, mean test active prevalence 0.77)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.016 +/- 0.018 | 1.164 +/- 0.040 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.447 +/- 0.058 | 0.858 +/- 0.051 | 0.677 +/- 0.033 | 0.809 +/- 0.025 |
| RF-desc | 0.226 +/- 0.026 | 1.016 +/- 0.046 | 0.482 +/- 0.036 | 0.712 +/- 0.017 |
| RF-ecfp | 0.509 +/- 0.037 | 0.808 +/- 0.040 | 0.711 +/- 0.027 | 0.832 +/- 0.024 |
| RF-all | 0.496 +/- 0.032 | 0.819 +/- 0.040 | 0.704 +/- 0.029 | 0.825 +/- 0.023 |
| XGB-all | 0.496 +/- 0.033 | 0.819 +/- 0.039 | 0.703 +/- 0.026 | 0.828 +/- 0.021 |
| XGB-scramble | -0.074 +/- 0.043 | 1.197 +/- 0.049 | 0.033 +/- 0.091 | 0.511 +/- 0.047 |

### 1.3 Paired comparisons across splits

| target | comparison | mean_dR2 | wins | n_seeds | wilcoxon_p |
|---|---|---|---|---|---|
| EGFR | XGB-all vs kNN-Tanimoto | 0.0174 | 6 | 10 | 0.1934 |
| EGFR | RF-ecfp vs kNN-Tanimoto | 0.0426 | 9 | 10 | 0.0098 |
| EGFR | RF-all vs RF-ecfp | -0.0264 | 0 | 10 | 0.0020 |
| EGFR | XGB-all vs RF-ecfp | -0.0252 | 1 | 10 | 0.0039 |
| FGFR4 | XGB-all vs kNN-Tanimoto | 0.0572 | 9 | 10 | 0.0039 |
| FGFR4 | RF-ecfp vs kNN-Tanimoto | 0.0674 | 10 | 10 | 0.0020 |
| FGFR4 | RF-all vs RF-ecfp | -0.0190 | 1 | 10 | 0.0039 |
| FGFR4 | XGB-all vs RF-ecfp | -0.0102 | 2 | 10 | 0.1055 |
| VEGFR2 | XGB-all vs kNN-Tanimoto | 0.0496 | 9 | 10 | 0.0059 |
| VEGFR2 | RF-ecfp vs kNN-Tanimoto | 0.0626 | 10 | 10 | 0.0020 |
| VEGFR2 | RF-all vs RF-ecfp | -0.0128 | 1 | 10 | 0.0059 |
| VEGFR2 | XGB-all vs RF-ecfp | -0.0130 | 2 | 10 | 0.0098 |

_paired over scaffold splits of one dataset; measures split noise, not generalisation to new chemistry_

### 1.4 Temporal holdout (train on older literature)

_Not available._ FGFR4: no publication year in curated data; EGFR: no publication year in curated data; VEGFR2: no publication year in curated data

### 1.5 Model selection

Selected family: **RF-ecfp** (highest mean R2 across targets on repeated random scaffold splits). Mean R2 by family: `{"RF-desc": 0.281, "RF-ecfp": 0.55, "RF-all": 0.53, "XGB-all": 0.534}`. Step 04 calibrates and step 05 screens with this family, so the model used downstream is the one this benchmark actually favours.

Controls: `mean` and `XGB-scramble` (shuffled labels) must sit at R2 <= 0 and AUC ~ 0.5. If they do not, the pipeline leaks.

## 2. Applicability domain, conformal calibration and NP transfer

Mondrian split-conformal intervals (90%), with the taxonomy defined by max-Tanimoto to the training split. `tau` is read off an isotonic (non-increasing) fit of half-width against similarity, so a single noisy bin cannot discard usable lower bins, and is reported with a bootstrap CI and a sensitivity analysis.

> **Validity statement.** The 90% guarantee holds for held-out ChEMBL scaffolds, which are exchangeable with the calibration set. It does **not** hold for natural-product libraries screened in step 05: the covariate shift between synthetic calibration chemistry and natural products breaks exchangeability, so intervals attached to library compounds are uncalibrated extrapolations and must be described as such wherever they are reported

### FGFR4 (RF-ecfp; train 2459 / cal 608 / test 541)

- empirical coverage **0.874** [0.847, 0.900] against a nominal 0.90; mean half-width 1.37 pIC50 units
- tau = **0.0** (95% CI 0.0-0.7); held-out molecules inside the domain: 100.0%

- tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give.

**tau at other width budgets**

| max_half_width | tau | binding |
|---|---|---|
| 1.000 | - | - |
| 1.250 | 0.700 | True |
| 1.500 | 0.000 | False |
| 1.650 | 0.000 | False |
| 2.000 | 0.000 | False |

| bin | n_cal | n_test | half_width_90 | half_width_90_isotonic | test_MAE | test_coverage | coverage 95% CI |
|---|---|---|---|---|---|---|---|
| [0.00, 0.70) | 328 | 309 | 1.484 | 1.484 | 0.852 | 0.822 | [0.777, 0.864] |
| [0.70, 1.00) | 280 | 232 | 1.214 | 1.214 | 0.521 | 0.944 | [0.914, 0.970] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.000 | -0.100 | 1.000 | 541.000 | 0.429 | 0.874 |
| 0.000 | -0.050 | 1.000 | 541.000 | 0.429 | 0.874 |
| 0.000 | 0.000 | 1.000 | 541.000 | 0.429 | 0.874 |
| 0.050 | 0.050 | 1.000 | 541.000 | 0.429 | 0.874 |
| 0.100 | 0.100 | 1.000 | 541.000 | 0.429 | 0.874 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 541 | 0.429 | [0.357, 0.494] | 0.933 |  |
| inside domain | 541 | 0.429 | [0.357, 0.494] | 0.933 |  |
| outside domain | 0 | - | - | - | fewer than 20 molecules: no metric is reported because the confidence interval would span most of the plausible range |
| NP-like (Ertl score >= 0) | 20 | 0.029 | [-0.805, 0.451] | 0.871 |  |
| not NP-like | 521 | 0.431 | [0.359, 0.498] | 0.935 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (14 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.950 [0.850, 1.000] (n = 20) versus 0.871 [0.843, 0.898] on the rest


### EGFR (RF-ecfp; train 6818 / cal 1498 / test 1466)

- empirical coverage **0.892** [0.875, 0.907] against a nominal 0.90; mean half-width 1.51 pIC50 units
- tau = **0.7** (95% CI 0.6-0.7); held-out molecules inside the domain: 46.2%

- tau excludes part of the similarity range; compounds below it fall outside the applicability domain.

**tau at other width budgets**

| max_half_width | tau | binding |
|---|---|---|
| 1.000 | - | - |
| 1.250 | - | - |
| 1.500 | 0.700 | True |
| 1.650 | 0.700 | True |
| 2.000 | 0.000 | False |

| bin | n_cal | n_test | half_width_90 | half_width_90_isotonic | test_MAE | test_coverage | coverage 95% CI |
|---|---|---|---|---|---|---|---|
| [0.00, 0.50) | 153 | 413 | 1.708 | 1.761 | 1.018 | 0.852 | [0.818, 0.886] |
| [0.50, 0.60) | 308 | 166 | 1.815 | 1.761 | 0.827 | 0.916 | [0.867, 0.958] |
| [0.60, 0.70) | 460 | 210 | 1.656 | 1.656 | 0.739 | 0.919 | [0.881, 0.957] |
| [0.70, 1.00) | 577 | 677 | 1.273 | 1.273 | 0.552 | 0.901 | [0.876, 0.923] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.600 | -0.100 | 0.605 | 887.000 | 0.604 | 0.905 |
| 0.650 | -0.050 | 0.545 | 799.000 | 0.625 | 0.905 |
| 0.700 | 0.000 | 0.462 | 677.000 | 0.639 | 0.901 |
| 0.750 | 0.050 | 0.357 | 523.000 | 0.683 | 0.914 |
| 0.800 | 0.100 | 0.212 | 311.000 | 0.724 | 0.945 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 1466 | 0.448 | [0.395, 0.493] | 0.973 |  |
| inside domain | 677 | 0.639 | [0.577, 0.689] | 0.738 |  |
| outside domain | 789 | 0.207 | [0.130, 0.280] | 1.136 |  |
| NP-like (Ertl score >= 0) | 54 | 0.072 | [-0.129, 0.306] | 1.056 |  |
| not NP-like | 1412 | 0.443 | [0.392, 0.491] | 0.969 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (29 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.944 [0.870, 1.000] (n = 54) versus 0.890 [0.873, 0.905] on the rest


### VEGFR2 (RF-ecfp; train 6387 / cal 1418 / test 1376)

- empirical coverage **0.881** [0.863, 0.898] against a nominal 0.90; mean half-width 1.46 pIC50 units
- tau = **0.5** (95% CI 0.5-0.6); held-out molecules inside the domain: 71.2%

- tau excludes part of the similarity range; compounds below it fall outside the applicability domain.

**tau at other width budgets**

| max_half_width | tau | binding |
|---|---|---|
| 1.000 | - | - |
| 1.250 | - | - |
| 1.500 | 0.600 | True |
| 1.650 | 0.500 | True |
| 2.000 | 0.000 | False |

| bin | n_cal | n_test | half_width_90 | half_width_90_isotonic | test_MAE | test_coverage | coverage 95% CI |
|---|---|---|---|---|---|---|---|
| [0.00, 0.50) | 274 | 396 | 1.714 | 1.714 | 0.963 | 0.821 | [0.780, 0.856] |
| [0.50, 0.60) | 250 | 167 | 1.574 | 1.574 | 0.793 | 0.868 | [0.814, 0.916] |
| [0.60, 0.70) | 332 | 271 | 1.411 | 1.411 | 0.658 | 0.911 | [0.878, 0.945] |
| [0.70, 1.00) | 562 | 542 | 1.270 | 1.270 | 0.530 | 0.913 | [0.891, 0.935] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.400 | -0.100 | 0.818 | 1125.000 | 0.443 | 0.892 |
| 0.450 | -0.050 | 0.761 | 1047.000 | 0.479 | 0.901 |
| 0.500 | 0.000 | 0.712 | 980.000 | 0.512 | 0.905 |
| 0.550 | 0.050 | 0.646 | 889.000 | 0.546 | 0.911 |
| 0.600 | 0.100 | 0.591 | 813.000 | 0.561 | 0.913 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 1376 | 0.376 | [0.332, 0.418] | 0.934 |  |
| inside domain | 980 | 0.512 | [0.464, 0.556] | 0.805 |  |
| outside domain | 396 | 0.041 | [-0.026, 0.098] | 1.196 |  |
| NP-like (Ertl score >= 0) | 40 | -0.392 | [-1.330, 0.070] | 1.049 |  |
| not NP-like | 1336 | 0.377 | [0.335, 0.420] | 0.931 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (24 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.825 [0.700, 0.925] (n = 40) versus 0.882 [0.866, 0.900] on the rest


## 3. Full-library screen

_Not run. Produce it with:_ `python scripts/05_screen_coconut_library.py`

## 4. Orthogonal structure-based validation (AutoDock Vina)

_The docking outputs in this repository predate the current validation protocol: they report best-of-N rather than top-1 RMSD without symmetry correction, a 12-vs-12 enrichment control drawn from potency extremes, and cross-receptor score comparisons. They are withheld here. Regenerate with:_ `python scripts/06_docking_validation.py` _(needs `vina` and `obabel`; see `docs/REPRODUCE.md`)._

## 5. Limitations

- No inactive class for FGFR4: the curated set contains no `>` records, so reported AUC/MCC describe ranking inside a 76%-active set rather than active-versus-inactive discrimination.

- Conformal coverage is guaranteed only for held-out ChEMBL scaffolds, which are exchangeable with the calibration set. It does not hold for natural-product libraries screened in step 05: the covariate shift between synthetic calibration chemistry and natural products breaks exchangeability, so intervals attached to library compounds are uncalibrated extrapolations and must be described as such wherever they are reported

- insufficient NP-like calibration molecules (14 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- The NP-like held-out subset is small (n = 20 for FGFR4); its metrics carry wide confidence intervals and support a direction, not a point claim.

- The applicability-domain threshold is non-binding for FGFR4: tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give. Any in-domain hit count for this target is therefore a count of the whole library, not of a calibrated-reliable subset.

- No experimental validation. Every candidate is a computational prediction and the paper claims nothing beyond that.
