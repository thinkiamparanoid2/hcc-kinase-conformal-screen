# HCC kinase natural-product screen: results report

_Auto-generated from `results/metrics` and `results/predictions`. Every number, table and caveat below is read from those files; nothing is hand-typed._

## 0. Data curation

Structures are salt-stripped, neutralised and tautomer-canonicalised, then deduplicated on the parent InChIKey. Censored measurements are labelled, not silently treated as point values: only uncensored rows are used as regression targets.

| target | final compounds | exact labels | censored (potent bound) | censored (inactive bound) | active prevalence (exact) | R2 ceiling from label noise |
|---|---|---|---|---|---|---|
| FGFR4 | 4819 | 3609 | 132 | 1078 | 0.763 | 0.993 |
| EGFR | 12484 | 9821 | 1004 | 1659 | 0.743 | 0.949 |
| VEGFR2 | 11269 | 9228 | 375 | 1666 | 0.769 | 0.996 |

## 1. Benchmark: scaffold-split evaluation

Scaffold splits are generic Bemis-Murcko and strictly disjoint. The canonical split is deterministic (largest scaffold families in train) and carries percentile bootstrap CIs over test molecules. The repeated splits measure how much the result moves with the scaffold assignment: that spread is SPLIT variance, not sampling variance, because the same molecules are reused each time.

_Hyperparameters are fixed a priori at common library defaults and were not tuned on any split reported here; no model selection was performed against test-set performance. Tuning, if added, must run as an inner cross-validation loop on the training scaffolds only._

### 1.1 Canonical deterministic scaffold split (with 95% CIs)

**FGFR4** (train 2888, test 721, test active prevalence 0.76, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.236 | [0.172, 0.296] | 1.101 | 0.482 | 0.773 | [0.737, 0.810] |
| RF-ecfp | 0.477 | [0.417, 0.532] | 0.910 | 0.688 | 0.857 | [0.824, 0.890] |
| RF-all | 0.433 | [0.373, 0.488] | 0.948 | 0.659 | 0.843 | [0.811, 0.876] |
| XGB-all | 0.449 | [0.385, 0.508] | 0.934 | 0.662 | 0.844 | [0.811, 0.876] |
| kNN-Tanimoto | 0.380 | [0.296, 0.458] | 0.991 | 0.640 | 0.832 | [0.795, 0.867] |

**EGFR** (train 7859, test 1962, test active prevalence 0.70, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.191 | [0.147, 0.235] | 1.189 | 0.490 | 0.740 | [0.716, 0.766] |
| RF-ecfp | 0.428 | [0.379, 0.472] | 1.000 | 0.679 | 0.817 | [0.796, 0.837] |
| RF-all | 0.404 | [0.353, 0.451] | 1.020 | 0.663 | 0.807 | [0.785, 0.827] |
| XGB-all | 0.364 | [0.311, 0.413] | 1.054 | 0.637 | 0.792 | [0.768, 0.813] |
| kNN-Tanimoto | 0.234 | [0.158, 0.303] | 1.157 | 0.623 | 0.778 | [0.753, 0.802] |

**VEGFR2** (train 7384, test 1844, test active prevalence 0.69, scaffold overlap 0)

| model | R2 | R2 95% CI | RMSE | Spearman | AUC | AUC 95% CI |
|---|---|---|---|---|---|---|
| RF-desc | 0.127 | [0.085, 0.166] | 1.096 | 0.394 | 0.683 | [0.657, 0.709] |
| RF-ecfp | 0.408 | [0.372, 0.442] | 0.902 | 0.634 | 0.799 | [0.777, 0.819] |
| RF-all | 0.400 | [0.365, 0.433] | 0.908 | 0.636 | 0.799 | [0.778, 0.819] |
| XGB-all | 0.391 | [0.353, 0.427] | 0.915 | 0.621 | 0.789 | [0.768, 0.811] |
| kNN-Tanimoto | 0.289 | [0.231, 0.345] | 0.989 | 0.599 | 0.772 | [0.748, 0.795] |

### 1.2 Repeated random scaffold splits (mean +/- SD over seeds)

**FGFR4** (10 splits, mean test active prevalence 0.78)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.008 +/- 0.009 | 1.196 +/- 0.064 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.411 +/- 0.080 | 0.911 +/- 0.055 | 0.652 +/- 0.053 | 0.856 +/- 0.033 |
| RF-desc | 0.302 +/- 0.050 | 0.993 +/- 0.041 | 0.537 +/- 0.039 | 0.795 +/- 0.031 |
| RF-ecfp | 0.501 +/- 0.044 | 0.840 +/- 0.047 | 0.689 +/- 0.043 | 0.872 +/- 0.029 |
| RF-all | 0.483 +/- 0.054 | 0.854 +/- 0.039 | 0.675 +/- 0.053 | 0.863 +/- 0.026 |
| XGB-all | 0.484 +/- 0.047 | 0.854 +/- 0.039 | 0.682 +/- 0.043 | 0.859 +/- 0.026 |
| XGB-scramble | -0.182 +/- 0.049 | 1.294 +/- 0.066 | -0.016 +/- 0.056 | 0.486 +/- 0.043 |

**EGFR** (10 splits, mean test active prevalence 0.72)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.012 +/- 0.015 | 1.369 +/- 0.059 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.528 +/- 0.048 | 0.934 +/- 0.056 | 0.746 +/- 0.027 | 0.871 +/- 0.017 |
| RF-desc | 0.335 +/- 0.058 | 1.109 +/- 0.058 | 0.581 +/- 0.045 | 0.800 +/- 0.027 |
| RF-ecfp | 0.572 +/- 0.047 | 0.889 +/- 0.056 | 0.766 +/- 0.027 | 0.888 +/- 0.018 |
| RF-all | 0.556 +/- 0.057 | 0.905 +/- 0.060 | 0.754 +/- 0.028 | 0.883 +/- 0.018 |
| XGB-all | 0.545 +/- 0.064 | 0.916 +/- 0.068 | 0.746 +/- 0.037 | 0.878 +/- 0.018 |
| XGB-scramble | -0.099 +/- 0.039 | 1.427 +/- 0.049 | -0.024 +/- 0.075 | 0.486 +/- 0.041 |

**VEGFR2** (10 splits, mean test active prevalence 0.78)

| model | R2 | RMSE | Spearman | AUC |
|---|---|---|---|---|
| mean | -0.011 +/- 0.011 | 1.140 +/- 0.027 | 0.000 +/- 0.000 | 0.500 +/- 0.000 |
| kNN-Tanimoto | 0.468 +/- 0.038 | 0.826 +/- 0.032 | 0.680 +/- 0.023 | 0.821 +/- 0.018 |
| RF-desc | 0.210 +/- 0.043 | 1.007 +/- 0.040 | 0.458 +/- 0.041 | 0.709 +/- 0.035 |
| RF-ecfp | 0.514 +/- 0.027 | 0.791 +/- 0.037 | 0.707 +/- 0.015 | 0.842 +/- 0.013 |
| RF-all | 0.502 +/- 0.029 | 0.800 +/- 0.036 | 0.699 +/- 0.024 | 0.834 +/- 0.021 |
| XGB-all | 0.496 +/- 0.031 | 0.805 +/- 0.035 | 0.692 +/- 0.019 | 0.833 +/- 0.018 |
| XGB-scramble | -0.088 +/- 0.029 | 1.182 +/- 0.022 | 0.025 +/- 0.047 | 0.519 +/- 0.038 |

### 1.3 Paired comparisons across splits

| target | comparison | mean_dR2 | wins | n_seeds | wilcoxon_p |
|---|---|---|---|---|---|
| EGFR | XGB-all vs kNN-Tanimoto | 0.0170 | 8 | 10 | 0.1934 |
| EGFR | RF-ecfp vs kNN-Tanimoto | 0.0445 | 7 | 10 | 0.0273 |
| EGFR | RF-all vs RF-ecfp | -0.0159 | 3 | 10 | 0.0273 |
| EGFR | XGB-all vs RF-ecfp | -0.0274 | 2 | 10 | 0.0137 |
| FGFR4 | XGB-all vs kNN-Tanimoto | 0.0730 | 10 | 10 | 0.0020 |
| FGFR4 | RF-ecfp vs kNN-Tanimoto | 0.0902 | 10 | 10 | 0.0020 |
| FGFR4 | RF-all vs RF-ecfp | -0.0182 | 1 | 10 | 0.0488 |
| FGFR4 | XGB-all vs RF-ecfp | -0.0171 | 2 | 10 | 0.1309 |
| VEGFR2 | XGB-all vs kNN-Tanimoto | 0.0274 | 9 | 10 | 0.0039 |
| VEGFR2 | RF-ecfp vs kNN-Tanimoto | 0.0453 | 9 | 10 | 0.0039 |
| VEGFR2 | RF-all vs RF-ecfp | -0.0118 | 2 | 10 | 0.0488 |
| VEGFR2 | XGB-all vs RF-ecfp | -0.0179 | 1 | 10 | 0.0098 |

_paired over scaffold splits of one dataset; measures split noise, not generalisation to new chemistry_

### 1.4 Temporal holdout (train on older literature)

| target | model | cutoff year | n test | R2 | R2 95% CI | AUC |
|---|---|---|---|---|---|---|
| FGFR4 | RF-desc | 2022.000 | 794 | 0.079 | [0.009, 0.142] | 0.620 |
| FGFR4 | RF-ecfp | 2022.000 | 794 | 0.074 | [0.030, 0.113] | 0.611 |
| FGFR4 | RF-all | 2022.000 | 794 | 0.085 | [0.030, 0.135] | 0.645 |
| FGFR4 | XGB-all | 2022.000 | 794 | 0.097 | [0.031, 0.154] | 0.651 |
| FGFR4 | kNN-Tanimoto | 2022.000 | 794 | -0.152 | [-0.240, -0.069] | 0.600 |
| EGFR | RF-desc | 2021.000 | 2779 | -0.135 | [-0.184, -0.087] | 0.669 |
| EGFR | RF-ecfp | 2021.000 | 2779 | -0.852 | [-0.955, -0.755] | 0.490 |
| EGFR | RF-all | 2021.000 | 2779 | -0.537 | [-0.613, -0.463] | 0.525 |
| EGFR | XGB-all | 2021.000 | 2779 | -0.541 | [-0.623, -0.464] | 0.576 |
| EGFR | kNN-Tanimoto | 2021.000 | 2779 | -0.533 | [-0.630, -0.442] | 0.625 |
| VEGFR2 | RF-desc | 2018.000 | 2017 | 0.032 | [0.003, 0.057] | 0.593 |
| VEGFR2 | RF-ecfp | 2018.000 | 2017 | 0.208 | [0.174, 0.241] | 0.757 |
| VEGFR2 | RF-all | 2018.000 | 2017 | 0.160 | [0.130, 0.189] | 0.727 |
| VEGFR2 | XGB-all | 2018.000 | 2017 | 0.270 | [0.238, 0.301] | 0.767 |
| VEGFR2 | kNN-Tanimoto | 2018.000 | 2017 | 0.112 | [0.064, 0.163] | 0.724 |

### 1.5 Model selection

Selected family: **RF-ecfp** (highest mean R2 across targets on repeated random scaffold splits). Mean R2 by family: `{"RF-desc": 0.282, "RF-ecfp": 0.529, "RF-all": 0.514, "XGB-all": 0.508}`. Step 04 calibrates and step 05 screens with this family, so the model used downstream is the one this benchmark actually favours.

Controls: `mean` and `XGB-scramble` (shuffled labels) must sit at R2 <= 0 and AUC ~ 0.5. If they do not, the pipeline leaks.

## 2. Applicability domain, conformal calibration and NP transfer

Mondrian split-conformal intervals (90%), with the taxonomy defined by max-Tanimoto to the training split. `tau` is read off an isotonic (non-increasing) fit of half-width against similarity, so a single noisy bin cannot discard usable lower bins, and is reported with a bootstrap CI and a sensitivity analysis.

> **Validity statement.** The 90% guarantee holds for held-out ChEMBL scaffolds, which are exchangeable with the calibration set. It does **not** hold for natural-product libraries screened in step 05: the covariate shift between synthetic calibration chemistry and natural products breaks exchangeability, so intervals attached to library compounds are uncalibrated extrapolations and must be described as such wherever they are reported

### FGFR4 (RF-ecfp; train 2515 / cal 553 / test 541)

- empirical coverage **0.871** [0.841, 0.900] against a nominal 0.90; mean half-width 1.36 pIC50 units
- tau = **0.0** (95% CI 0.0-0.6); held-out molecules inside the domain: 100.0%

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
| [0.00, 0.60) | 159 | 205 | 1.495 | 1.495 | 0.890 | 0.800 | [0.741, 0.854] |
| [0.60, 0.70) | 188 | 100 | 1.425 | 1.425 | 0.739 | 0.890 | [0.830, 0.950] |
| [0.70, 1.00) | 206 | 236 | 1.219 | 1.219 | 0.522 | 0.924 | [0.886, 0.953] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.000 | -0.100 | 1.000 | 541.000 | 0.410 | 0.871 |
| 0.000 | -0.050 | 1.000 | 541.000 | 0.410 | 0.871 |
| 0.000 | 0.000 | 1.000 | 541.000 | 0.410 | 0.871 |
| 0.050 | 0.050 | 1.000 | 541.000 | 0.410 | 0.871 |
| 0.100 | 0.100 | 1.000 | 541.000 | 0.410 | 0.871 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 541 | 0.410 | [0.338, 0.475] | 0.947 |  |
| inside domain | 541 | 0.410 | [0.338, 0.475] | 0.947 |  |
| outside domain | 0 | - | - | - | fewer than 20 molecules: no metric is reported because the confidence interval would span most of the plausible range |
| NP-like (Ertl score >= 0) | 20 | 0.035 | [-0.906, 0.401] | 0.953 |  |
| not NP-like | 521 | 0.408 | [0.338, 0.475] | 0.947 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (22 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.850 [0.650, 1.000] (n = 20) versus 0.871 [0.843, 0.898] on the rest


### EGFR (RF-ecfp; train 6815 / cal 1533 / test 1473)

- empirical coverage **0.858** [0.840, 0.876] against a nominal 0.90; mean half-width 1.28 pIC50 units
- tau = **0.0** (95% CI 0.0-0.6); held-out molecules inside the domain: 100.0%

- tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give.

**tau at other width budgets**

| max_half_width | tau | binding |
|---|---|---|
| 1.000 | 0.700 | True |
| 1.250 | 0.700 | True |
| 1.500 | 0.600 | True |
| 1.650 | 0.000 | False |
| 2.000 | 0.000 | False |

| bin | n_cal | n_test | half_width_90 | half_width_90_isotonic | test_MAE | test_coverage | coverage 95% CI |
|---|---|---|---|---|---|---|---|
| [0.00, 0.50) | 195 | 401 | 1.642 | 1.642 | 1.034 | 0.845 | [0.810, 0.880] |
| [0.50, 0.60) | 188 | 177 | 1.581 | 1.581 | 0.789 | 0.881 | [0.836, 0.927] |
| [0.60, 0.70) | 445 | 220 | 1.346 | 1.346 | 0.704 | 0.850 | [0.800, 0.895] |
| [0.70, 1.00) | 705 | 675 | 0.954 | 0.954 | 0.496 | 0.862 | [0.836, 0.887] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.000 | -0.100 | 1.000 | 1473.000 | 0.461 | 0.858 |
| 0.000 | -0.050 | 1.000 | 1473.000 | 0.461 | 0.858 |
| 0.000 | 0.000 | 1.000 | 1473.000 | 0.461 | 0.858 |
| 0.050 | 0.050 | 1.000 | 1473.000 | 0.461 | 0.858 |
| 0.100 | 0.100 | 1.000 | 1473.000 | 0.461 | 0.858 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 1473 | 0.461 | [0.400, 0.516] | 0.976 |  |
| inside domain | 1473 | 0.461 | [0.400, 0.516] | 0.976 |  |
| outside domain | 0 | - | - | - | fewer than 20 molecules: no metric is reported because the confidence interval would span most of the plausible range |
| NP-like (Ertl score >= 0) | 47 | 0.043 | [-0.175, 0.229] | 1.304 |  |
| not NP-like | 1426 | 0.458 | [0.397, 0.513] | 0.963 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (9 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.894 [0.809, 0.979] (n = 47) versus 0.857 [0.839, 0.874] on the rest


### VEGFR2 (RF-ecfp; train 6429 / cal 1415 / test 1384)

- empirical coverage **0.892** [0.876, 0.909] against a nominal 0.90; mean half-width 1.45 pIC50 units
- tau = **0.5** (95% CI 0.5-0.5); held-out molecules inside the domain: 68.6%

- tau excludes part of the similarity range; compounds below it fall outside the applicability domain.

**tau at other width budgets**

| max_half_width | tau | binding |
|---|---|---|
| 1.000 | - | - |
| 1.250 | 0.700 | True |
| 1.500 | 0.500 | True |
| 1.650 | 0.500 | True |
| 2.000 | 0.000 | False |

| bin | n_cal | n_test | half_width_90 | half_width_90_isotonic | test_MAE | test_coverage | coverage 95% CI |
|---|---|---|---|---|---|---|---|
| [0.00, 0.50) | 182 | 434 | 1.959 | 1.959 | 0.905 | 0.901 | [0.871, 0.929] |
| [0.50, 0.60) | 201 | 171 | 1.343 | 1.374 | 0.837 | 0.807 | [0.749, 0.865] |
| [0.60, 0.70) | 388 | 251 | 1.405 | 1.374 | 0.629 | 0.912 | [0.876, 0.944] |
| [0.70, 1.00) | 644 | 528 | 1.100 | 1.100 | 0.505 | 0.903 | [0.879, 0.928] |

**tau sensitivity**

| tau | delta | frac_test_in_AD | n_test_in_AD | R2_in_AD | coverage_in_AD |
|---|---|---|---|---|---|
| 0.400 | -0.100 | 0.789 | 1092.000 | 0.477 | 0.890 |
| 0.450 | -0.050 | 0.736 | 1018.000 | 0.513 | 0.889 |
| 0.500 | 0.000 | 0.686 | 950.000 | 0.540 | 0.888 |
| 0.550 | 0.050 | 0.618 | 856.000 | 0.585 | 0.901 |
| 0.600 | 0.100 | 0.563 | 779.000 | 0.605 | 0.906 |

| subset | n | R2 | R2 95% CI | RMSE | note |
|---|---|---|---|---|---|
| all held-out | 1384 | 0.400 | [0.357, 0.441] | 0.917 |  |
| inside domain | 950 | 0.540 | [0.493, 0.583] | 0.792 |  |
| outside domain | 434 | 0.066 | [-0.000, 0.125] | 1.143 |  |
| NP-like (Ertl score >= 0) | 45 | -0.573 | [-1.447, -0.082] | 1.147 |  |
| not NP-like | 1339 | 0.405 | [0.362, 0.445] | 0.908 |  |

**Natural-product transfer**

- insufficient NP-like calibration molecules (28 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- coverage on NP-like held-out molecules: 0.889 [0.800, 0.978] (n = 45) versus 0.892 [0.875, 0.908] on the rest


## 3. Full-library natural-product screen

Library: 738,827 records collapsed to 690,177 unique 2D parents (690,177 featurised). Collapsed on standardised parent InChIKey. 2D descriptors and ECFP4 are achiral, so stereoisomer records produce identical predictions and must not be counted as distinct screened compounds.

> **Interval caveat.** `interval_lo`/`interval_hi` are the step-04 Mondrian half-widths applied by similarity bin. Their 90% coverage was established on held-out ChEMBL scaffolds, which are NOT exchangeable with this natural-product library. Treat them as a calibrated-on-synthetic-chemistry error model extrapolated to natural products, never as guaranteed 90% intervals.

| target | model | tau | tau_is_binding | in_AD_caveat | in_AD | in_AD_frac | median_sim_train | in_AD_active | in_AD_active_novel_scaffold | in_AD_active_pan_kinase_alert | interval_lower_bound_active | interval_lower_bound_active_novel_scaffold | rediscovery_total | known_scaffold_total | spearman_model_vs_nn_baseline_in_AD |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | RF-ecfp | 0.000 | False | tau is non-binding for this target (step 04 found every calibration bin inside the width budget), so `in_AD` admits the whole library and the in-domain counts below are not a calibrated-reliable subset. | 690177 | 1.000 | 0.182 | 643633 | 604014 | 90889 | 3220 | 2977 | 276 | 48475 | 0.239 |
| EGFR | RF-ecfp | 0.000 | False | tau is non-binding for this target (step 04 found every calibration bin inside the width budget), so `in_AD` admits the whole library and the in-domain counts below are not a calibrated-reliable subset. | 690177 | 1.000 | 0.217 | 123659 | 116267 | 21087 | 34 | 1 | 226 | 89541 | 0.317 |
| VEGFR2 | RF-ecfp | 0.500 | True | - | 2078 | 0.003 | 0.206 | 977 | 553 | 309 | 111 | 7 | 151 | 61097 | 0.721 |

### 3.1 Consensus across targets

| any_target_inAD_active | two_or_more_targets_inAD_active | all_three_inAD_active | any_target_interval_active | novel_scaffold_and_not_rediscovery | pan_kinase_alert_among_hits |
|---|---|---|---|---|---|
| 645850 | 122205 | 214 | 3271 | 3063 | 293 |

Cross-target prediction rank correlation: `{"FGFR4_vs_EGFR": 0.363, "FGFR4_vs_VEGFR2": 0.164, "EGFR_vs_VEGFR2": 0.18}`

> The three targets are ATP-site kinases trained on overlapping chemotypes, so their predictions are correlated (see cross_target_prediction_spearman). Agreement across targets indicates pan-kinase chemistry, which is a selectivity liability, and is not independent confirmation of activity.

### 3.2 Top candidates on scaffolds absent from the training sets

| identifier | name | n_targets_interval_active | n_targets_inAD_active | mean_pred | np_likeness | QED | PAINS | pan_kinase_alert |
|---|---|---|---|---|---|---|---|---|
| CNP0573320.0 | todralazine | 2 | 3 | 7.457 | -1.600 | 0.788 | False | False |
| CNP0216199.2 | - | 2 | 3 | 7.373 | 0.900 | 0.380 | False | False |
| CNP0552352.1 | KT-5822 | 2 | 3 | 7.373 | 0.750 | 0.380 | False | False |
| CNP0216199.1 | K252a-Me | 2 | 3 | 7.373 | 0.900 | 0.380 | False | False |
| CNP0599845.0 | Gilteritinib | 2 | 2 | 7.782 | -1.100 | 0.428 | True | False |
| CNP0506664.2 | KT 5720 | 2 | 2 | 7.210 | 0.880 | 0.217 | False | False |
| CNP0506664.3 | KT 5720 | 2 | 2 | 7.210 | 0.880 | 0.217 | False | False |
| CNP0340362.1 | KT-5720 | 2 | 2 | 7.210 | 0.720 | 0.217 | False | False |
| CNP0506664.1 | hexyl (15S,16R,18S)-16-hydroxy-15-methyl-3-oxo-28-oxa-4,14,19-triazaoctacyclo[12.11.2.1¹?,¹?.0²,?.0?,²?.0?,¹³.0¹?,²?.0²?,²?]octacosa-1,6,8,10,12,20,22,24,26-nonaene-16-carboxylate | 2 | 2 | 7.210 | 0.880 | 0.217 | False | False |
| CNP0279633.1 | (2S,3R,6R)-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaene-4,16-dione | 2 | 2 | 7.123 | 1.100 | 0.412 | False | False |
| CNP0433234.2 | 3'-epi-N-acetyl-holyrine A | 2 | 2 | 6.892 | 0.940 | 0.307 | False | False |
| CNP0433234.1 | 3'-N-acetyl-holyrine A | 2 | 2 | 6.892 | 0.940 | 0.307 | False | False |
| CNP0497076.1 | ~{N}-[(2~{S},3~{R},4~{R},6~{R})-6-(12,19-dihydroxy-3,13,23-triazahexacyclo[14.7.0.0^{2,10}.0^{4,9}.0^{11,15}.0^{17,22}]tricosa-1,4,6,8,10,12,15,17,19,21-decaen-3-yl)-3-hydroxy-2-methyl-tetrahydropyran-4-yl]acetamide | 2 | 2 | 6.887 | 0.860 | 0.254 | False | False |
| CNP0482271.1 | ~{N}-[(2~{S},3~{R},4~{R},6~{R})-6-(7,12-dihydroxy-3,13,23-triazahexacyclo[14.7.0.0^{2,10}.0^{4,9}.0^{11,15}.0^{17,22}]tricosa-1,4,6,8,10,12,15,17,19,21-decaen-3-yl)-3-hydroxy-2-methyl-tetrahydropyran-4-yl]acetamide | 2 | 2 | 6.814 | 0.880 | 0.254 | False | False |
| CNP0588398.0 | hydralazine hydrochloride | 1 | 3 | 7.529 | -1.220 | 0.481 | False | False |
| CNP0428917.1 | - | 1 | 3 | 7.231 | 0.690 | 0.370 | False | False |
| CNP0430874.1 | - | 1 | 3 | 7.113 | 0.750 | 0.295 | False | False |
| CNP0112734.7 | Hamigeromycin C | 1 | 3 | 6.989 | 1.810 | 0.577 | True | True |
| CNP0427479.1 | - | 1 | 3 | 6.956 | 1.220 | 0.322 | False | False |
| CNP0436341.1 | - | 1 | 3 | 6.943 | 1.230 | 0.322 | False | False |
| CNP0536972.0 | 4-N-methyl-4-N-(3-methyl-1H-indazol-6-yl)-2-N-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine | 1 | 3 | 6.919 | -1.000 | 0.459 | False | True |
| CNP0243075.0 | 4-(3-amino-4-{[(4-chloro-3-methylphenyl)carbamothioyl]amino}phenoxy)-N-methylpyridine-2-carboxamide | 1 | 3 | 6.894 | -1.700 | 0.339 | False | False |
| CNP0246325.0 | 4-{[(6,7-dimethoxy-2-methyl-4-oxoquinazolin-3(4H)-yl)acetyl]amino}benzamide | 1 | 3 | 6.838 | -1.520 | 0.650 | False | True |
| CNP0505548.0 | ethyl 4-(2-{3-oxo-2-phenyl-1,4,8-triazaspiro[4.5]dec-1-en-8-yl}acetamido)benzoate | 1 | 3 | 6.791 | -1.190 | 0.679 | False | False |
| CNP0465346.1 | - | 1 | 3 | 6.742 | -1.430 | 0.524 | False | False |

### 3.3 Hits on scaffolds already present in training (analogues, positive controls)

| identifier | name | n_targets_inAD_active | mean_pred |
|---|---|---|---|
| CNP0250924.4 | - | 3 | 7.500 |
| CNP0250924.6 | - | 3 | 7.500 |
| CNP0250924.1 | staurosporine | 3 | 7.500 |
| CNP0586714.2 | 62996-74-1 | 3 | 7.500 |
| CNP0250924.5 | staurosporine | 3 | 7.500 |
| CNP0250924.3 | staurosporine | 3 | 7.500 |
| CNP0250924.2 | - | 3 | 7.500 |
| CNP0184469.2 | GNF-Pf-70 | 3 | 7.465 |
| CNP0184469.1 | 3'-Demethylstaurosporine | 3 | 7.465 |
| CNP0184469.3 | - | 3 | 7.465 |
| CNP0109266.1 | - | 3 | 7.382 |
| CNP0554764.1 | 11-Hydroxystaurosporine | 3 | 7.374 |
| CNP0158097.2 | - | 3 | 7.374 |
| CNP0158097.3 | - | 3 | 7.374 |
| CNP0158097.1 | - | 3 | 7.374 |

## 4. Orthogonal structure-based validation (AutoDock Vina)

_Produced by Open Babel 3.1.0 -- Feb 17 2026 -- 06:52:02, AutoDock Vina v1.2.5, 16 threads per docking._

Boxes are centred on the co-crystal ligand heavy-atom centroid in chain A and sized to the ligand (minimum 22.0 A, 8.0 A padding), pH 7.4 protonation, Gasteiger charges. Redocking and candidate docking: exhaustiveness 32, 3 seeds per docking. Enrichment control: exhaustiveness 32, 3 seed(s).

> **Score comparability.** AutoDock Vina scores are empirical and receptor-specific. They are not transferable binding free energies, so a weaker score in one pocket than another is not evidence of selectivity. Compare within a target (see the <T>_vina_z columns), and treat any score weaker than -4.0 kcal/mol as a placement failure to investigate, not as a measured steric clash.

### 4.1 Control 1: co-crystal redocking (pose prediction)

Top-1 RMSD is the control: it is what the method would produce prospectively. Best-of-N is a sampling diagnostic and is reported only beside it.

| target | ligand | box (A) | Vina top-1 (kcal/mol) | top-1 RMSD (A) | pass at 2 A | best-of-N RMSD (A) | N poses | native pose rank | score gap (kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | BLU9931-analog (6LF) | 22.000 | -9.127 | 3.334 | False | 1.811 | 9 | 7 | -0.823 |
| EGFR | Erlotinib (AQ4) | 23.800 | -7.809 | 8.004 | False | 1.283 | 9 | 4 | -0.166 |
| VEGFR2 | Sorafenib (BAX) | 27.000 | -12.201 | 0.577 | True | 0.577 | 9 | 1 | 0.000 |

Pose prediction succeeds at the conventional 2 A criterion for 1 of 3 pockets.

Where top-1 fails, the search still produced a near-native pose: FGFR4 ranks it 7, EGFR ranks it 4, separated from the top-scoring pose by at most 0.82 kcal/mol. Vina's scoring function is not quoted to that precision, so these are ranking failures rather than sampling failures: the pose is found and cannot be told apart from its rivals. This is the reason docking enters this work as an orthogonal check and not as a source of truth about binding mode.

### 4.2 Control 2: active-versus-inactive enrichment

| target | n act / n inact | active mean | inactive mean | ROC-AUC | AUC 95% CI | EF1% | EF5% | BEDROC | beats chance |
|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | 50 / 50 | -9.180 | -8.700 | 0.696 | [0.592, 0.791] | 2.000 | 2.000 | 0.890 | ROC_AUC, EF5pct, BEDROC_alpha20 |
| EGFR | 50 / 49 | -9.390 | -9.490 | 0.444 | [0.327, 0.560] | 1.980 | 1.190 | 0.540 | none |
| VEGFR2 | 50 / 50 | -10.110 | -9.390 | 0.693 | [0.585, 0.792] | 2.000 | 1.600 | 0.899 | ROC_AUC, BEDROC_alpha20 |

> **Verdict.** Docking enrichment is target-dependent and is not established here as a general filter. Global separation (AUC outside its null): FGFR4, VEGFR2. Early recognition (EF5% or BEDROC outside its null): FGFR4, VEGFR2. EGFR scores below 0.5, so the ranking is there if anything inverted. The two are different claims and no pocket has both, so the result is not one effect seen three times. With 12 null comparisons in this table and EF5% and BEDROC strongly correlated, a single metric clearing its null in a single pocket is weak evidence on its own; treat it as a direction to confirm on a larger set, not as a validated filter. These are not averaged. A pocket where docking does not enrich cannot corroborate a prediction for a compound in that pocket.

_Sampling: random within class, MW-matched; NOT potency extremes._

> **How to read these metrics.** Top 1% of 100 compounds is 1 compound, and the active fraction is 0.50, not the ~1% that EF and BEDROC were designed for. Read each metric against its null_baselines entry: exceeds_null is whether the value falls outside the middle 90% of what a random ranking gives on this same design.

Each metric against a label-shuffled null on the same design:

| target | metric | observed | random mean | random 5-95% | p | outside null |
|---|---|---|---|---|---|---|
| FGFR4 | ROC_AUC | 0.696 | 0.500 | 0.4022 - 0.5975 | 0.001 | True |
| FGFR4 | EF1pct | 2.000 | 0.982 | 0.0 - 2.0 | 0.491 | False |
| FGFR4 | EF5pct | 2.000 | 0.999 | 0.4 - 1.6 | 0.028 | True |
| FGFR4 | BEDROC_alpha20 | 0.890 | 0.500 | 0.2585 - 0.7504 | 0.002 | True |
| EGFR | ROC_AUC | 0.445 | 0.499 | 0.4029 - 0.5922 | 0.824 | False |
| EGFR | EF1pct | 1.980 | 1.033 | 0.0 - 1.98 | 0.521 | False |
| EGFR | EF5pct | 1.188 | 0.995 | 0.396 - 1.584 | 0.505 | False |
| EGFR | BEDROC_alpha20 | 0.540 | 0.503 | 0.2572 - 0.7426 | 0.411 | False |
| VEGFR2 | ROC_AUC | 0.693 | 0.499 | 0.4032 - 0.594 | 0.000 | True |
| VEGFR2 | EF1pct | 2.000 | 1.002 | 0.0 - 2.0 | 0.501 | False |
| VEGFR2 | EF5pct | 1.600 | 0.997 | 0.4 - 1.6 | 0.178 | False |
| VEGFR2 | BEDROC_alpha20 | 0.899 | 0.499 | 0.2399 - 0.7438 | 0.003 | True |

> **What the negative class is.** censored '>' records (true inactives). This matters for how the AUC reads: a true inactive is a `>` record, where the assay ran and found nothing, whereas a measured weak binder still occupies the pocket. Discriminating actives from weak binders is the harder task, so an AUC measured this way understates what the same control would give against true non-binders, and is not comparable to a published enrichment that used decoys.

### 4.3 Candidate docking by chemotype

| chemotype | n | FGFR4 mean | FGFR4 z | FGFR4 unplaced | EGFR mean | EGFR z | EGFR unplaced | VEGFR2 mean | VEGFR2 z | VEGFR2 unplaced |
|---|---|---|---|---|---|---|---|---|---|---|
| resorcylic_acid_lactone | 2 | -8.044 | 1.146 | 0 | -8.460 | 1.121 | 0 | -7.471 | 1.102 | 0 |
| macrolactone | 2 | -8.044 | 1.146 | 0 | -8.460 | 1.121 | 0 | -7.471 | 1.102 | 0 |
| catechol | 4 | -7.938 | 1.227 | 0 | -8.613 | 1.009 | 0 | -8.236 | 0.166 | 0 |

Chemotype membership is decided by the SMARTS patterns in `common.CHEMOTYPE_SMARTS`. `unplaced` counts dockings whose best score was weaker than the plausibility floor, i.e. the ligand could not be placed in the box -- those are docking failures to investigate, not evidence of a steric clash.

### 4.4 Candidate shortlist

| identifier | name | mean_pred | np_likeness | QED | PAINS | pan_kinase_alert | FGFR4_vina_kcal | FGFR4_vina_sd | EGFR_vina_kcal | EGFR_vina_sd | VEGFR2_vina_kcal | VEGFR2_vina_sd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CNP0573320.0 | todralazine | 7.457 | -1.600 | 0.788 | False | False | -7.466 | 0.028 | -7.648 | 0.012 | -7.978 | 0.016 |
| CNP0216199.2 | - | 7.373 | 0.900 | 0.380 | False | False | -10.508 | 0.003 | -10.686 | 0.042 | -7.932 | 0.006 |
| CNP0552352.1 | KT-5822 | 7.373 | 0.750 | 0.380 | False | False | -10.244 | 0.010 | -10.814 | 0.011 | -8.090 | 0.010 |
| CNP0216199.1 | K252a-Me | 7.373 | 0.900 | 0.380 | False | False | -9.629 | 0.072 | -10.835 | 0.105 | -7.607 | 0.008 |
| CNP0599845.0 | Gilteritinib | 7.782 | -1.100 | 0.428 | True | False | -9.001 | 0.005 | -9.623 | 0.021 | -7.924 | 0.262 |
| CNP0506664.2 | KT 5720 | 7.210 | 0.880 | 0.217 | False | False | - | - | - | - | - | - |
| CNP0506664.3 | KT 5720 | 7.210 | 0.880 | 0.217 | False | False | -10.189 | 0.006 | -9.856 | 0.007 | -7.917 | 0.017 |
| CNP0340362.1 | KT-5720 | 7.210 | 0.720 | 0.217 | False | False | -9.666 | 0.028 | -10.318 | 0.030 | -8.469 | 0.006 |
| CNP0506664.1 | hexyl (15S,16R,18S)-16-hydroxy-15-methyl-3-oxo-28-oxa-4,14,19-triazaoctacyclo[12.11.2.1¹?,¹?.0²,?.0?,²?.0?,¹³.0¹?,²?.0²?,²?]octacosa-1,6,8,10,12,20,22,24,26-nonaene-16-carboxylate | 7.210 | 0.880 | 0.217 | False | False | - | - | - | - | - | - |
| CNP0279633.1 | (2S,3R,6R)-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaene-4,16-dione | 7.123 | 1.100 | 0.412 | False | False | -10.282 | 0.007 | -11.527 | 0.009 | -7.824 | 0.003 |
| CNP0433234.2 | 3'-epi-N-acetyl-holyrine A | 6.892 | 0.940 | 0.307 | False | False | -12.221 | 0.008 | -11.670 | 0.023 | -8.756 | 0.012 |
| CNP0433234.1 | 3'-N-acetyl-holyrine A | 6.892 | 0.940 | 0.307 | False | False | -12.021 | 0.003 | -11.357 | 0.020 | -8.545 | 0.434 |
| CNP0497076.1 | ~{N}-[(2~{S},3~{R},4~{R},6~{R})-6-(12,19-dihydroxy-3,13,23-triazahexacyclo[14.7.0.0^{2,10}.0^{4,9}.0^{11,15}.0^{17,22}]tricosa-1,4,6,8,10,12,15,17,19,21-decaen-3-yl)-3-hydroxy-2-methyl-tetrahydropyran-4-yl]acetamide | 6.887 | 0.860 | 0.254 | False | False | -11.167 | 0.007 | -10.524 | 0.028 | -9.109 | 0.015 |
| CNP0482271.1 | ~{N}-[(2~{S},3~{R},4~{R},6~{R})-6-(7,12-dihydroxy-3,13,23-triazahexacyclo[14.7.0.0^{2,10}.0^{4,9}.0^{11,15}.0^{17,22}]tricosa-1,4,6,8,10,12,15,17,19,21-decaen-3-yl)-3-hydroxy-2-methyl-tetrahydropyran-4-yl]acetamide | 6.814 | 0.880 | 0.254 | False | False | -9.880 | 0.156 | -10.395 | 0.015 | -9.437 | 0.017 |
| CNP0588398.0 | hydralazine hydrochloride | 7.529 | -1.220 | 0.481 | False | False | -6.699 | 0.010 | -7.016 | 0.008 | -6.791 | 0.003 |
| CNP0428917.1 | - | 7.231 | 0.690 | 0.370 | False | False | -10.770 | 0.002 | -10.659 | 0.022 | -8.067 | 0.010 |
| CNP0430874.1 | - | 7.113 | 0.750 | 0.295 | False | False | -10.518 | 0.016 | -12.127 | 0.005 | -7.928 | 0.186 |
| CNP0112734.7 | Hamigeromycin C | 6.989 | 1.810 | 0.577 | True | True | -7.637 | 0.010 | -8.710 | 0.006 | -7.226 | 0.015 |
| CNP0427479.1 | - | 6.956 | 1.220 | 0.322 | False | False | -10.495 | 0.003 | -11.068 | 0.010 | -8.226 | 0.005 |
| CNP0436341.1 | - | 6.943 | 1.230 | 0.322 | False | False | -10.447 | 0.024 | -11.287 | 0.019 | -7.842 | 0.173 |
| CNP0536972.0 | 4-N-methyl-4-N-(3-methyl-1H-indazol-6-yl)-2-N-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine | 6.919 | -1.000 | 0.459 | False | True | -8.484 | 0.007 | -9.082 | 0.012 | -8.555 | 0.011 |
| CNP0243075.0 | 4-(3-amino-4-{[(4-chloro-3-methylphenyl)carbamothioyl]amino}phenoxy)-N-methylpyridine-2-carboxamide | 6.894 | -1.700 | 0.339 | False | False | -8.508 | 0.166 | -9.031 | 0.013 | -9.506 | 0.019 |
| CNP0246325.0 | 4-{[(6,7-dimethoxy-2-methyl-4-oxoquinazolin-3(4H)-yl)acetyl]amino}benzamide | 6.838 | -1.520 | 0.650 | False | True | -7.922 | 0.024 | -8.159 | 0.165 | -8.570 | 0.003 |
| CNP0505548.0 | ethyl 4-(2-{3-oxo-2-phenyl-1,4,8-triazaspiro[4.5]dec-1-en-8-yl}acetamido)benzoate | 6.791 | -1.190 | 0.679 | False | False | -8.918 | 0.068 | -8.746 | 0.125 | -10.752 | 0.082 |
| CNP0465346.1 | - | 6.742 | -1.430 | 0.524 | False | False | -8.772 | 0.022 | -8.443 | 0.022 | -9.674 | 0.025 |

## 5. Limitations

- Conformal coverage is guaranteed only for held-out ChEMBL scaffolds, which are exchangeable with the calibration set. It does not hold for natural-product libraries screened in step 05: the covariate shift between synthetic calibration chemistry and natural products breaks exchangeability, so intervals attached to library compounds are uncalibrated extrapolations and must be described as such wherever they are reported

- insufficient NP-like calibration molecules (22 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- The NP-like held-out subset is small (n = 20 for FGFR4); its metrics carry wide confidence intervals and support a direction, not a point claim.

- The applicability-domain threshold is non-binding for FGFR4: tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give. Any in-domain hit count for this target is therefore a count of the whole library, not of a calibrated-reliable subset.

- The applicability-domain threshold is non-binding for EGFR: tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give. Any in-domain hit count for this target is therefore a count of the whole library, not of a calibrated-reliable subset.

- The three targets are ATP-site kinases trained on overlapping chemotypes, so their predictions are correlated (see cross_target_prediction_spearman). Agreement across targets indicates pan-kinase chemistry, which is a selectivity liability, and is not independent confirmation of activity.

- Top-1 redocking fails the 2 A criterion for EGFR, FGFR4, so pose-level conclusions in those pockets are not supported.

- Docking does not enrich actives over the negative class for EGFR: no metric clears its label-shuffled null. A docking score therefore carries no validated discriminative weight in those pockets and cannot be used to support or filter a prediction there.

- Where docking does enrich, it does so inconsistently: no pocket shows both global separation and early recognition, so the enrichment control licenses a direction to confirm on a larger set rather than a general-purpose structure-based filter.

- AutoDock Vina scores are empirical and receptor-specific. They are not transferable binding free energies, so a weaker score in one pocket than another is not evidence of selectivity. Compare within a target (see the <T>_vina_z columns), and treat any score weaker than -4.0 kcal/mol as a placement failure to investigate, not as a measured steric clash.

- No experimental validation. Every candidate is a computational prediction and the paper claims nothing beyond that.
