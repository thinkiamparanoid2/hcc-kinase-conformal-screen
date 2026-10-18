# HCC kinase natural-product screen: results report

_Auto-generated from `results/metrics` and `results/predictions`. Combines 10-seed scaffold-split QSAR, Mondrian conformal calibration, full 738,823-compound COCONUT screening, and orthogonal AutoDock Vina 1.2.5 structure-based validation._

## 1. Benchmark: repeated scaffold splits (80/20)

10 random Murcko-scaffold splits per target; mean +/- SD. Classification columns use predicted pIC50 as the score (active = pIC50 >= 6). F1 is not reported because ~75-83% of labels are active.

### FGFR4  (n=3682, mean active prevalence in test = 0.79)

| model | R2 | RMSE | Spearman | AUC | PR-AUC | MCC |
|---|---|---|---|---|---|---|
| mean | -0.021 +/- 0.018 | 1.247 +/- 0.058 | 0.000 +/- 0.000 | 0.500 +/- 0.000 | 0.791 +/- 0.032 | 0.000 +/- 0.000 |
| kNN-Tanimoto | 0.499 +/- 0.089 | 0.868 +/- 0.045 | 0.693 +/- 0.069 | 0.863 +/- 0.039 | 0.953 +/- 0.011 | 0.570 +/- 0.082 |
| RF-desc | 0.320 +/- 0.048 | 1.017 +/- 0.029 | 0.551 +/- 0.039 | 0.805 +/- 0.026 | 0.932 +/- 0.015 | 0.360 +/- 0.066 |
| RF-ecfp | 0.544 +/- 0.065 | 0.830 +/- 0.035 | 0.722 +/- 0.055 | 0.885 +/- 0.021 | 0.963 +/- 0.005 | 0.570 +/- 0.083 |
| RF-all | 0.530 +/- 0.059 | 0.843 +/- 0.026 | 0.717 +/- 0.046 | 0.882 +/- 0.022 | 0.963 +/- 0.005 | 0.525 +/- 0.099 |
| XGB-all | 0.535 +/- 0.068 | 0.837 +/- 0.030 | 0.713 +/- 0.056 | 0.876 +/- 0.025 | 0.959 +/- 0.007 | 0.550 +/- 0.081 |
| XGB-scramble | -0.187 +/- 0.093 | 1.343 +/- 0.051 | -0.044 +/- 0.090 | 0.502 +/- 0.047 | 0.801 +/- 0.034 | 0.020 +/- 0.053 |

### EGFR  (n=9844, mean active prevalence in test = 0.75)

| model | R2 | RMSE | Spearman | AUC | PR-AUC | MCC |
|---|---|---|---|---|---|---|
| mean | -0.008 +/- 0.009 | 1.319 +/- 0.033 | 0.000 +/- 0.000 | 0.500 +/- 0.000 | 0.753 +/- 0.030 | 0.000 +/- 0.000 |
| kNN-Tanimoto | 0.529 +/- 0.067 | 0.899 +/- 0.051 | 0.731 +/- 0.053 | 0.875 +/- 0.016 | 0.948 +/- 0.010 | 0.563 +/- 0.046 |
| RF-desc | 0.210 +/- 0.182 | 1.160 +/- 0.123 | 0.473 +/- 0.152 | 0.768 +/- 0.052 | 0.906 +/- 0.021 | 0.304 +/- 0.066 |
| RF-ecfp | 0.585 +/- 0.052 | 0.845 +/- 0.046 | 0.769 +/- 0.040 | 0.897 +/- 0.015 | 0.962 +/- 0.005 | 0.535 +/- 0.053 |
| RF-all | 0.549 +/- 0.074 | 0.879 +/- 0.063 | 0.749 +/- 0.053 | 0.893 +/- 0.014 | 0.961 +/- 0.005 | 0.523 +/- 0.061 |
| XGB-all | 0.553 +/- 0.061 | 0.876 +/- 0.052 | 0.750 +/- 0.044 | 0.889 +/- 0.014 | 0.959 +/- 0.005 | 0.522 +/- 0.070 |
| XGB-scramble | -0.065 +/- 0.047 | 1.356 +/- 0.048 | 0.025 +/- 0.081 | 0.513 +/- 0.059 | 0.753 +/- 0.051 | 0.040 +/- 0.057 |

### VEGFR2  (n=9247, mean active prevalence in test = 0.78)

| model | R2 | RMSE | Spearman | AUC | PR-AUC | MCC |
|---|---|---|---|---|---|---|
| mean | -0.021 +/- 0.028 | 1.157 +/- 0.051 | 0.000 +/- 0.000 | 0.500 +/- 0.000 | 0.782 +/- 0.028 | 0.000 +/- 0.000 |
| kNN-Tanimoto | 0.483 +/- 0.043 | 0.822 +/- 0.030 | 0.696 +/- 0.030 | 0.828 +/- 0.023 | 0.938 +/- 0.015 | 0.461 +/- 0.055 |
| RF-desc | 0.229 +/- 0.056 | 1.004 +/- 0.039 | 0.483 +/- 0.067 | 0.715 +/- 0.034 | 0.896 +/- 0.022 | 0.195 +/- 0.053 |
| RF-ecfp | 0.537 +/- 0.032 | 0.778 +/- 0.032 | 0.728 +/- 0.027 | 0.848 +/- 0.023 | 0.949 +/- 0.014 | 0.450 +/- 0.044 |
| RF-all | 0.518 +/- 0.036 | 0.794 +/- 0.030 | 0.718 +/- 0.030 | 0.836 +/- 0.024 | 0.945 +/- 0.014 | 0.407 +/- 0.044 |
| XGB-all | 0.518 +/- 0.034 | 0.794 +/- 0.026 | 0.714 +/- 0.030 | 0.837 +/- 0.025 | 0.945 +/- 0.014 | 0.427 +/- 0.048 |
| XGB-scramble | -0.083 +/- 0.037 | 1.192 +/- 0.056 | 0.030 +/- 0.055 | 0.509 +/- 0.037 | 0.786 +/- 0.036 | -0.006 +/- 0.024 |

### Paired comparison (per-seed R2 differences)

| target | A vs B | mean dR2 | A wins | Wilcoxon p |
|---|---|---|---|---|
| FGFR4 | XGB-all vs kNN-Tanimoto | 0.0363 | 9/10 | 0.0039 |
| FGFR4 | RF-all vs kNN-Tanimoto | 0.0305 | 8/10 | 0.0273 |
| FGFR4 | RF-all vs RF-ecfp | -0.0146 | 3/10 | 0.1309 |
| FGFR4 | XGB-all vs RF-all | 0.0058 | 6/10 | 0.4316 |
| EGFR | XGB-all vs kNN-Tanimoto | 0.0247 | 8/10 | 0.0371 |
| EGFR | RF-all vs kNN-Tanimoto | 0.0202 | 7/10 | 0.1055 |
| EGFR | RF-all vs RF-ecfp | -0.0357 | 0/10 | 0.0020 |
| EGFR | XGB-all vs RF-all | 0.0045 | 6/10 | 0.4316 |
| VEGFR2 | XGB-all vs kNN-Tanimoto | 0.0343 | 9/10 | 0.0098 |
| VEGFR2 | RF-all vs kNN-Tanimoto | 0.0343 | 9/10 | 0.0098 |
| VEGFR2 | RF-all vs RF-ecfp | -0.0195 | 0/10 | 0.0020 |
| VEGFR2 | XGB-all vs RF-all | 0.0000 | 5/10 | 0.9219 |

Controls: `mean` and `XGB-scramble` (shuffled labels) must sit at R2<=0 / AUC~0.5; if they do not, the pipeline leaks.

## 2. Applicability domain, calibration and natural-product domain shift

Mondrian split-conformal intervals (90%) with bins on max-Tanimoto to the train split. The domain threshold `tau` is the similarity above which every bin has 90% half-width <= 1.65 pIC50 units (corresponding to RMSE <= 1.0 log unit at the 90% confidence level).

### FGFR4  (train 2546 / cal 559 / test 577)

- overall empirical coverage: **0.849** (target 0.90); mean half-width 1.24; tau = 0.5; held-out molecules inside domain: 87.3%

| bin | n_test | half_width_90 | test_MAE | test_coverage |
|---|---|---|---|---|
| [0.00, 0.40) | 38 | 1.567 | 0.870 | 0.868 |
| [0.40, 0.50) | 35 | 1.864 | 1.049 | 0.857 |
| [0.50, 0.60) | 58 | 1.331 | 0.830 | 0.897 |
| [0.60, 0.70) | 145 | 1.379 | 0.670 | 0.890 |
| [0.70, 1.00) | 301 | 1.036 | 0.602 | 0.817 |

| subset | n | R2 | RMSE | AUC |
|---|---|---|---|---|
| all held-out | 577 | 0.469 | 0.878 | 0.849 |
| inside domain | 504 | 0.520 | 0.834 | 0.875 |
| outside domain | 73 | 0.004 | 1.137 | 0.621 |
| NP-like (score>=0) | 30 | 0.071 | 0.958 | 0.657 |
| not NP-like | 547 | 0.467 | 0.873 | 0.853 |

### EGFR  (train 6734 / cal 1485 / test 1625)

- overall empirical coverage: **0.918** (target 0.90); mean half-width 1.45; tau = 0.6; held-out molecules inside domain: 74.7%

| bin | n_test | half_width_90 | test_MAE | test_coverage |
|---|---|---|---|---|
| [0.00, 0.30) | 24 | 1.743 | 1.116 | 0.875 |
| [0.30, 0.40) | 42 | 1.749 | 0.973 | 0.833 |
| [0.40, 0.50) | 109 | 1.360 | 0.774 | 0.862 |
| [0.50, 0.60) | 236 | 2.033 | 0.788 | 0.962 |
| [0.60, 0.70) | 430 | 1.627 | 0.716 | 0.923 |
| [0.70, 1.00) | 784 | 1.165 | 0.517 | 0.915 |

| subset | n | R2 | RMSE | AUC |
|---|---|---|---|---|
| all held-out | 1625 | 0.584 | 0.840 | 0.900 |
| inside domain | 1214 | 0.590 | 0.766 | 0.922 |
| outside domain | 411 | 0.304 | 1.027 | 0.792 |
| NP-like (score>=0) | 20 | -0.229 | 0.701 | 0.333 |
| not NP-like | 1605 | 0.580 | 0.841 | 0.901 |

### VEGFR2  (train 6444 / cal 1415 / test 1388)

- overall empirical coverage: **0.921** (target 0.90); mean half-width 1.49; tau = 0.6; held-out molecules inside domain: 66.8%

| bin | n_test | half_width_90 | test_MAE | test_coverage |
|---|---|---|---|---|
| [0.00, 0.40) | 101 | 1.934 | 0.616 | 0.990 |
| [0.40, 0.50) | 121 | 1.883 | 0.787 | 0.975 |
| [0.50, 0.60) | 239 | 1.663 | 0.891 | 0.908 |
| [0.60, 0.70) | 398 | 1.520 | 0.713 | 0.925 |
| [0.70, 1.00) | 529 | 1.213 | 0.525 | 0.900 |

| subset | n | R2 | RMSE | AUC |
|---|---|---|---|---|
| all held-out | 1388 | 0.462 | 0.858 | 0.828 |
| inside domain | 927 | 0.540 | 0.791 | 0.862 |
| outside domain | 461 | 0.300 | 0.979 | 0.752 |
| NP-like (score>=0) | 43 | -0.107 | 0.987 | 0.829 |
| not NP-like | 1345 | 0.474 | 0.853 | 0.831 |

## 3. Full-library screen

Library: 738,827 compounds (738,823 parsed). `in_AD_active` = inside domain AND predicted pIC50 >= 6. `confident active` = inside domain AND lower 90% bound >= pIC50 6. `rediscovery` = Tanimoto >= 0.90 to a ChEMBL training molecule (known chemistry, positive control - not a novel lead). Tables below collapse 2D-identical COCONUT stereoisomers (`CNPxxxxxxx.*`, count in `n_iso`) since 2D descriptors and ECFP4 are achiral.

| target | tau | in_AD | in_AD_frac | median_sim_train | pred_ge_7 | in_AD_active (in_AD & pred>=6) | in_AD_active_not_rediscovery | confident_active (in_AD & lo>=6) | confident_active_not_rediscovery | rediscovery_total |
|---|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | 0.500 | 4389 | 0.006 | 0.189 | 121357 | 3242 | 3185 | 212 | 178 | 295 |
| EGFR | 0.600 | 1824 | 0.002 | 0.220 | 24613 | 678 | 631 | 41 | 25 | 224 |
| VEGFR2 | 0.600 | 656 | 0.001 | 0.212 | 7227 | 370 | 285 | 53 | 20 | 167 |

Consensus: `{"any_target_confident": 270, "two_or_more_targets_confident": 36, "all_three_confident": 0, "confident_and_novel": 189, "any_target_inAD_active": 3906, "two_or_more_targets_inAD_active": 305, "inAD_active_and_novel": 3746}`

### Top novel consensus candidates (non-rediscovery, deduplicated by 2D parent)

| identifier | name | n_iso | n_targets_confident | n_targets_inAD_active | mean_pred | np_likeness | QED | PAINS |
|---|---|---|---|---|---|---|---|---|
| CNP0216199.1 | K252a-Me | 2 | 2 | 3 | 7.595 | 0.900 | 0.380 | False |
| CNP0351785.1 | Indocarbazostatin C | 1 | 2 | 2 | 6.898 | 0.950 | 0.240 | False |
| CNP0143079.3 | (2S,3S,5S,7R,9Z,12S)-2,3,16-trihydroxy-18-methoxy-12-methyl-6,13-dioxatricyclo[13.4.0.05,7]nonadeca-1(19),9,15,17-tetraene-8,14-dione | 7 | 2 | 2 | 6.831 | 2.370 | 0.489 | False |
| CNP0277321.1 | 15-O-desmethyl-(5Z)-7-oxozeaenol | 2 | 2 | 2 | 6.703 | 2.280 | 0.521 | False |
| CNP0184469.2 | GNF-Pf-70 | 3 | 1 | 3 | 7.415 | 1.160 | 0.362 | False |
| CNP0113867.2 | - | 2 | 1 | 3 | 7.405 | 1.320 | 0.401 | False |
| CNP0109266.1 | - | 1 | 1 | 3 | 7.358 | 1.430 | 0.353 | False |
| CNP0393527.1 | O-demethyl-N-demethyl-staurosporine | 1 | 1 | 3 | 7.317 | 1.170 | 0.345 | False |
| CNP0158097.2 | - | 4 | 1 | 3 | 7.297 | 1.350 | 0.353 | False |
| CNP0141274.3 | Staurosporine M1 | 3 | 1 | 3 | 7.269 | 0.980 | 0.341 | False |
| CNP0427126.1 | (2S,3R,4R,6R)-3,4-dihydroxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 1 | 1 | 3 | 7.262 | 1.230 | 0.345 | False |
| CNP0114061.1 | - | 1 | 1 | 3 | 7.256 | 1.460 | 0.347 | False |
| CNP0322653.2 | 5'-hydroxystaurosporine | 2 | 1 | 3 | 7.234 | 1.140 | 0.359 | False |
| CNP0245131.1 | NSC714348 | 2 | 1 | 3 | 7.198 | 1.380 | 0.347 | False |
| CNP0279633.1 | (2S,3R,6R)-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaene-4,16-dione | 2 | 1 | 3 | 7.180 | 1.100 | 0.412 | False |
| CNP0173499.1 | N,N-Dimethylstaurosporine | 1 | 1 | 3 | 7.173 | 1.020 | 0.398 | False |
| CNP0105905.1 | N-[(2S,3R,4R,6R)-3-methoxy-2-methyl-16-oxo-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-4-yl]-N-methylformamide | 1 | 1 | 3 | 7.164 | 1.080 | 0.381 | False |
| CNP0212769.1 | - | 2 | 1 | 3 | 7.162 | 1.150 | 0.388 | False |
| CNP0205142.1 | (2S,3R,4E,6R)-4-hydroxyimino-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 4 | 1 | 3 | 7.148 | 0.980 | 0.279 | False |
| CNP0289949.1 | (2R,3S,4S,6S)-4-hydroxy-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 5 | 1 | 3 | 7.147 | 1.330 | 0.400 | False |
| CNP0343202.1 | RK-1409B | 2 | 1 | 3 | 7.147 | 1.330 | 0.400 | False |
| CNP0429188.1 | - | 1 | 1 | 3 | 7.081 | 1.020 | 0.363 | False |
| CNP0262333.1 | - | 1 | 1 | 3 | 7.061 | 1.160 | 0.379 | False |
| CNP0199234.1 | 4'-Demethyl-Af-formyl-7V-hydroxy-staurosporine | 1 | 1 | 3 | 7.009 | 1.070 | 0.224 | False |
| CNP0148778.0 | - | 1 | 1 | 3 | 7.008 | 0.990 | 0.285 | False |

### Rediscoveries (positive control: the screen recovers known kinase chemotypes)

| identifier | name | n_iso | n_targets_confident | n_targets_inAD_active | mean_pred |
|---|---|---|---|---|---|
| CNP0253706.2 | k-252a | 3 | 2 | 3 | 7.520 |
| CNP0590083.1 | Naquotinib | 1 | 2 | 2 | 7.595 |
| CNP0551736.1 | Bemcentinib | 1 | 2 | 2 | 7.327 |
| CNP0527731.1 | Lestaurtinib | 1 | 2 | 2 | 7.268 |
| CNP0590328.1 | PF-477736 | 1 | 2 | 2 | 7.185 |
| CNP0590982.0 | Delphinidin chloride | 1 | 2 | 2 | 7.185 |
| CNP0517555.0 | Verubulin | 1 | 2 | 2 | 7.082 |
| CNP0316772.5 | Cochliomycin F | 9 | 2 | 2 | 7.069 |
| CNP0287101.1 | (4S,6Z,9S,10R)-9,10,18-Trihydroxy-16-methoxy-4-methyl-3-oxabicyclo[12.4.0]octadeca-1(14),6,15,17-tetraene-2,8-dione | 6 | 2 | 2 | 7.006 |
| CNP0250924.3 | staurosporine | 6 | 1 | 3 | 7.605 |
| CNP0520440.0 | CUDC-101 | 1 | 1 | 3 | 7.044 |
| CNP0436884.0 | Erlotinib | 1 | 1 | 3 | 6.825 |
| CNP0530553.0 | whi-p180 | 1 | 1 | 3 | 6.672 |
| CNP0569685.0 | 153436-53-4 | 1 | 1 | 3 | 6.664 |
| CNP0328253.1 | 99570-78-2 | 1 | 1 | 2 | 7.465 |

## 4. Orthogonal structure-based validation (AutoDock Vina 1.2.5)

All binding sites were centered directly on the co-crystallized ligand heavy-atom centroid in Chain A (`1M17` AQ4/Erlotinib for EGFR, `5JKG` 6LF for FGFR4, `4ASD` BAX/Sorafenib for VEGFR2; 22 x 22 x 22 A cubic box, pH 7.4 protonation, Gasteiger partial charges).

### 4.1 Control 1: Co-crystal redocking

| target | co-crystal ligand | center (x, y, z) | Vina (kcal/mol) | top-1 RMSD (A) | best-of-9 RMSD (A) |
|---|---|---|---|---|---|
| FGFR4 | BLU9931-analog (6LF) | (-40.88, -16.01, 378.17) | -9.127 | 3.334 | 1.811 |
| EGFR | Erlotinib (AQ4) | (22.01, 0.25, 52.79) | -7.720 | 7.839 | 1.242 |
| VEGFR2 | Sorafenib (BAX) | (-24.62, -0.39, -10.93) | -12.254 | 0.833 | 0.833 |

### 4.2 Control 2: ChEMBL active vs. inactive enrichment

MW- and rotatable-bond-matched potent ChEMBL actives (`pIC50 >= 8.5`) vs. weak/inactives (`pIC50 <= 4.8`).

| target | n_act / n_inact | active mean (kcal/mol) | inactive mean (kcal/mol) | Docking ROC-AUC | Mann-Whitney p |
|---|---|---|---|---|---|
| FGFR4 | 12 / 12 | -9.7300 | -8.5100 | 0.8260 | 0.0036 |
| EGFR | 12 / 12 | -9.0400 | -9.0400 | 0.5210 | 0.4426 |
| VEGFR2 | 12 / 12 | -10.8000 | -8.3500 | 0.9310 | 0.0002 |

- **VEGFR2 (`4ASD`, DFG-out)** and **FGFR4 (`5JKG`, DFG-in)** strongly separate potent ChEMBL actives from inactives (ROC-AUC = **0.931** and **0.826**).
- **EGFR (`1M17`, open DFG-in hinge)** shows no active-vs-inactive separation by rigid-receptor Vina scoring alone (ROC-AUC = **0.521**, both means -9.04 kcal/mol), highlighting why calibrated ligand-based QSAR and orthogonal docking are complementary rather than interchangeable.

### 4.3 Multi-target docking of top 25 novel consensus candidates

| identifier | name | mean_pred | np_likeness | QED | PAINS | FGFR4_vina_kcal | EGFR_vina_kcal | VEGFR2_vina_kcal | mean_vina_kcal |
|---|---|---|---|---|---|---|---|---|---|
| CNP0216199.1 | K252a-Me | 7.595 | 0.900 | 0.380 | False | -9.668 | -10.749 | -4.902 | -8.440 |
| CNP0351785.1 | Indocarbazostatin C | 6.898 | 0.950 | 0.240 | False | -10.795 | -11.294 | -1.894 | -7.994 |
| CNP0143079.3 | (2S,3S,5S,7R,9Z,12S)-2,3,16-trihydroxy-18-methoxy-12-methyl-6,13-dioxatricyclo[13.4.0.05,7]nonadeca-1(19),9,15,17-tetraene-8,14-dione | 6.831 | 2.370 | 0.489 | False | -8.015 | -8.244 | -6.494 | -7.584 |
| CNP0277321.1 | 15-O-desmethyl-(5Z)-7-oxozeaenol | 6.703 | 2.280 | 0.521 | False | -8.442 | -8.932 | -7.110 | -8.161 |
| CNP0184469.2 | GNF-Pf-70 | 7.415 | 1.160 | 0.362 | False | -9.800 | -11.252 | -3.692 | -8.248 |
| CNP0113867.2 | - | 7.405 | 1.320 | 0.401 | False | -11.473 | -11.186 | -2.524 | -8.394 |
| CNP0109266.1 | - | 7.358 | 1.430 | 0.353 | False | -10.244 | -11.020 | -2.339 | -7.868 |
| CNP0393527.1 | O-demethyl-N-demethyl-staurosporine | 7.317 | 1.170 | 0.345 | False | -11.655 | -11.440 | -4.567 | -9.221 |
| CNP0158097.2 | - | 7.297 | 1.350 | 0.353 | False | -11.019 | -10.738 | -1.009 | -7.589 |
| CNP0141274.3 | Staurosporine M1 | 7.269 | 0.980 | 0.341 | False | -10.358 | -11.364 | -4.790 | -8.837 |
| CNP0427126.1 | (2S,3R,4R,6R)-3,4-dihydroxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 7.262 | 1.230 | 0.345 | False | -10.815 | -11.660 | -4.105 | -8.860 |
| CNP0114061.1 | - | 7.256 | 1.460 | 0.347 | False | -10.509 | -10.955 | -0.967 | -7.477 |
| CNP0322653.2 | 5'-hydroxystaurosporine | 7.234 | 1.140 | 0.359 | False | -8.760 | -10.858 | -3.920 | -7.846 |
| CNP0245131.1 | NSC714348 | 7.198 | 1.380 | 0.347 | False | -10.672 | -11.020 | -3.178 | -8.290 |
| CNP0279633.1 | (2S,3R,6R)-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaene-4,16-dione | 7.180 | 1.100 | 0.412 | False | -11.699 | -11.630 | -3.082 | -8.804 |
| CNP0173499.1 | N,N-Dimethylstaurosporine | 7.173 | 1.020 | 0.398 | False | -8.448 | -10.670 | -1.281 | -6.800 |
| CNP0105905.1 | N-[(2S,3R,4R,6R)-3-methoxy-2-methyl-16-oxo-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-4-yl]-N-methylformamide | 7.164 | 1.080 | 0.381 | False | -10.260 | -11.548 | -0.537 | -7.448 |
| CNP0212769.1 | - | 7.162 | 1.150 | 0.388 | False | -10.883 | -10.894 | -1.810 | -7.862 |
| CNP0205142.1 | (2S,3R,4E,6R)-4-hydroxyimino-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 7.148 | 0.980 | 0.279 | False | -10.406 | -11.304 | -3.206 | -8.305 |
| CNP0289949.1 | (2R,3S,4S,6S)-4-hydroxy-3-methoxy-2-methyl-29-oxa-1,7,17-triazaoctacyclo[12.12.2.12,6.07,28.08,13.015,19.020,27.021,26]nonacosa-8,10,12,14,19,21,23,25,27-nonaen-16-one | 7.147 | 1.330 | 0.400 | False | -10.095 | -10.935 | -3.257 | -8.096 |
| CNP0343202.1 | RK-1409B | 7.147 | 1.330 | 0.400 | False | -9.059 | -10.809 | -0.868 | -6.912 |
| CNP0429188.1 | - | 7.081 | 1.020 | 0.363 | False | -9.808 | -11.109 | -1.104 | -7.340 |
| CNP0262333.1 | - | 7.061 | 1.160 | 0.379 | False | -8.251 | -10.746 | -1.275 | -6.757 |
| CNP0199234.1 | 4'-Demethyl-Af-formyl-7V-hydroxy-staurosporine | 7.009 | 1.070 | 0.224 | False | -10.739 | -11.098 | -1.879 | -7.905 |
| CNP0148778.0 | - | 7.008 | 0.990 | 0.285 | False | -10.348 | -10.997 | -2.266 | -7.870 |

### 4.4 Rediscovery positive controls (known kinase inhibitors in COCONUT)

| identifier | name | mean_pred | np_likeness | QED | PAINS | FGFR4_vina_kcal | EGFR_vina_kcal | VEGFR2_vina_kcal | mean_vina_kcal |
|---|---|---|---|---|---|---|---|---|---|
| CNP0253706.2 | k-252a | 7.520 | 1.030 | 0.367 | False | -10.867 | -11.627 | -3.518 | -8.671 |
| CNP0590083.1 | Naquotinib | 7.595 | -1.120 | 0.444 | True | -8.599 | -8.924 | -8.459 | -8.661 |
| CNP0551736.1 | Bemcentinib | 7.327 | -1.080 | 0.381 | False | -10.127 | -11.342 | -9.746 | -10.405 |
| CNP0527731.1 | Lestaurtinib | 7.268 | 1.140 | 0.373 | False | -11.159 | -11.623 | -5.746 | -9.509 |
| CNP0590328.1 | PF-477736 | 7.185 | -0.780 | 0.517 | False | -9.426 | -9.491 | -10.417 | -9.778 |
| CNP0590982.0 | Delphinidin chloride | 7.185 | 1.640 | 0.263 | True | -7.695 | -8.530 | -8.612 | -8.279 |
| CNP0517555.0 | Verubulin | 7.082 | -1.100 | 0.732 | False | -7.921 | -7.659 | -7.101 | -7.560 |
| CNP0316772.5 | Cochliomycin F | 7.069 | 2.010 | 0.650 | False | -8.631 | -8.812 | -7.267 | -8.237 |

### 4.5 Conformational selectivity insight (DFG-in vs. DFG-out)

- **Indolocarbazole glycosides** (`K252a-Me`, `O-demethyl-N-demethyl-staurosporine`, `Staurosporine M1`, `Indocarbazostatin C`) are rigid hexacyclic **Type-I** ATP-competitive hinge binders. They achieve exceptional docking scores in the **DFG-in** ATP pockets of **FGFR4 (`5JKG`)** (-9.67 to -11.70 kcal/mol) and **EGFR (`1M17`)** (-10.67 to -11.66 kcal/mol), exceeding the native co-crystal ligands (-9.13 and -7.72 kcal/mol). However, in **VEGFR2 (`4ASD`)**, which is crystallized in the **DFG-out (Type-II, Sorafenib-bound)** inactive conformation where Phe1047 occludes the adenine subpocket, bulky rigid indolocarbazoles experience steric clash (-0.54 to -4.90 kcal/mol).
- **Resorcylic acid lactones (RALs)** (`15-O-desmethyl-(5Z)-7-oxozeaenol` / `CNP0277321.1`, NP-likeness = 2.28, QED = 0.52; and `CNP0143079.3`, NP-likeness = 2.37, QED = 0.49) possess a flexible 14-membered macrocyclic lactone that accommodates all three pockets, including DFG-out VEGFR2 (-8.44, -8.93, and -7.11 kcal/mol for `CNP0277321.1`; -8.63, -8.81, and -7.27 kcal/mol for the rediscovery control `Cochliomycin F`).

## 5. Limitations

- Models trained on ChEMBL (largely synthetic kinase inhibitors); natural products are mostly outside the chemical domain - see Section 2 for exact domain-shift degradation.
- Activity labels come from heterogeneous ChEMBL IC50 assays; label noise bounds achievable R2 (~0.52-0.59).
- Rigid-receptor Vina docking captures conformational steric filtering (e.g., DFG-in vs DFG-out selectivity in `4ASD`) and enriches actives in VEGFR2 (AUC 0.931) and FGFR4 (AUC 0.826), but does not separate actives from inactives in the solvent-exposed EGFR `1M17` hinge without induced-fit/MD refinement.
