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

> **REDUCED PROTOCOL: enrich_exhaustiveness 8 instead of 32, enrich_n_seeds 1 instead of 3. These numbers are not the protocol's own and should be regenerated before publication: unset DOCK_ENRICH_EXHAUSTIVENESS DOCK_ENRICH_N_SEEDS and re-run step 06 on hardware with more cores. Scores are insensitive to exhaustiveness in this pocket set (8, 16 and 32 agree to 0.01 kcal/mol on a test ligand), and the enrichment control's uncertainty is the bootstrap CI over compounds rather than the per-ligand search spread, so a reduced enrichment budget costs the replicate error bars and not the ranking itself.**

_Produced by Open Babel 3.1.0 -- Feb 17 2026 -- 06:48:53, AutoDock Vina v1.2.5, 4 threads per docking._

Boxes are centred on the co-crystal ligand heavy-atom centroid in chain A and sized to the ligand (minimum 22.0 A, 8.0 A padding), pH 7.4 protonation, Gasteiger charges. Redocking and candidate docking: exhaustiveness 32, 3 seeds per docking. Enrichment control: exhaustiveness 8, 1 seed(s).

> **Score comparability.** AutoDock Vina scores are empirical and receptor-specific. They are not transferable binding free energies, so a weaker score in one pocket than another is not evidence of selectivity. Compare within a target (see the <T>_vina_z columns), and treat any score weaker than -4.0 kcal/mol as a placement failure to investigate, not as a measured steric clash.

### 4.1 Control 1: co-crystal redocking (pose prediction)

Top-1 RMSD is the control: it is what the method would produce prospectively. Best-of-N is a sampling diagnostic and is reported only beside it.

| target | ligand | box (A) | Vina top-1 (kcal/mol) | top-1 RMSD (A) | pass at 2 A | best-of-N RMSD (A) | N poses | native pose rank | score gap (kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | BLU9931-analog (6LF) | 22.000 | -8.550 | 6.499 | False | 2.314 | 9 | 8 | -0.945 |
| EGFR | Erlotinib (AQ4) | 23.800 | -7.839 | 8.004 | False | 1.286 | 9 | 3 | -0.075 |
| VEGFR2 | Sorafenib (BAX) | 27.000 | -12.223 | 0.561 | True | 0.561 | 8 | 1 | 0.000 |

Pose prediction succeeds at the conventional 2 A criterion for 1 of 3 pockets.

Where top-1 fails, the search still produced a near-native pose: FGFR4 ranks it 8, EGFR ranks it 3, separated from the top-scoring pose by at most 0.94 kcal/mol. Vina's scoring function is not quoted to that precision, so these are ranking failures rather than sampling failures: the pose is found and cannot be told apart from its rivals. This is the reason docking enters this work as an orthogonal check and not as a source of truth about binding mode.

### 4.2 Control 2: active-versus-inactive enrichment

| target | n act / n inact | active mean | inactive mean | ROC-AUC | AUC 95% CI | EF1% | EF5% | BEDROC | beats chance |
|---|---|---|---|---|---|---|---|---|---|
| FGFR4 | 50 / 49 | -8.950 | -8.620 | 0.623 | [0.515, 0.730] | 0.000 | 1.580 | 0.721 | ROC_AUC |
| EGFR | 50 / 50 | -9.220 | -9.370 | 0.416 | [0.304, 0.535] | 2.000 | 1.200 | 0.565 | none |
| VEGFR2 | 49 / 50 | -9.650 | -9.460 | 0.536 | [0.417, 0.658] | 2.020 | 2.020 | 0.790 | EF5pct, BEDROC_alpha20 |

> **Verdict.** Docking enrichment is target-dependent and is not established here as a general filter. Global separation (AUC outside its null): FGFR4. Early recognition (EF5% or BEDROC outside its null): VEGFR2. EGFR scores below 0.5, so the ranking is there if anything inverted. The two are different claims and no pocket has both, so the result is not one effect seen three times. With 12 null comparisons in this table and EF5% and BEDROC strongly correlated, a single metric clearing its null in a single pocket is weak evidence on its own; treat it as a direction to confirm on a larger set, not as a validated filter. These are not averaged. A pocket where docking does not enrich cannot corroborate a prediction for a compound in that pocket.

_Sampling: random within class, MW-matched; NOT potency extremes._

> **How to read these metrics.** Top 1% of 99 compounds is 1 compound, and the active fraction is 0.51, not the ~1% that EF and BEDROC were designed for. Read each metric against its null_baselines entry: exceeds_null is whether the value falls outside the middle 90% of what a random ranking gives on this same design.

Each metric against a label-shuffled null on the same design:

| target | metric | observed | random mean | random 5-95% | p | outside null |
|---|---|---|---|---|---|---|
| FGFR4 | ROC_AUC | 0.623 | 0.500 | 0.4067 - 0.5972 | 0.018 | True |
| FGFR4 | EF1pct | 0.000 | 0.978 | 0.0 - 1.98 | 1.000 | False |
| FGFR4 | EF5pct | 1.584 | 0.998 | 0.396 - 1.584 | 0.174 | False |
| FGFR4 | BEDROC_alpha20 | 0.721 | 0.503 | 0.263 - 0.7431 | 0.071 | False |
| EGFR | ROC_AUC | 0.416 | 0.499 | 0.4012 - 0.5972 | 0.923 | False |
| EGFR | EF1pct | 2.000 | 0.984 | 0.0 - 2.0 | 0.492 | False |
| EGFR | EF5pct | 1.200 | 0.984 | 0.4 - 1.6 | 0.482 | False |
| EGFR | BEDROC_alpha20 | 0.565 | 0.494 | 0.245 - 0.7492 | 0.329 | False |
| VEGFR2 | ROC_AUC | 0.536 | 0.500 | 0.4065 - 0.5941 | 0.273 | False |
| VEGFR2 | EF1pct | 2.020 | 0.987 | 0.0 - 2.0204 | 0.488 | False |
| VEGFR2 | EF5pct | 2.020 | 1.004 | 0.4041 - 1.6163 | 0.028 | True |
| VEGFR2 | BEDROC_alpha20 | 0.790 | 0.497 | 0.2559 - 0.7441 | 0.024 | True |

> **What the negative class is.** measured weak binders; no censored '>' records in the committed data, so re-run step 01 with network access to populate a true inactive class. This matters for how the AUC reads: a true inactive is a `>` record, where the assay ran and found nothing, whereas a measured weak binder still occupies the pocket. Discriminating actives from weak binders is the harder task, so an AUC measured this way understates what the same control would give against true non-binders, and is not comparable to a published enrichment that used decoys.

### 4.3 Candidate docking

_Not available. The two controls above stand on their own: they need only the co-crystal structures and the curated training sets. Docking the candidate shortlist needs a current full-library screen, and there is none in this repository._

_Reason: results/predictions/consensus_hits.csv missing: there is no full-library screen to take candidates from. Run `python scripts/05_screen_coconut_library.py` (needs the COCONUT library; see data/README.md), then step 07 to annotate the shortlist._

## 5. Limitations

- No inactive class for FGFR4: the curated set contains no `>` records, so reported AUC/MCC describe ranking inside a 76%-active set rather than active-versus-inactive discrimination.

- Conformal coverage is guaranteed only for held-out ChEMBL scaffolds, which are exchangeable with the calibration set. It does not hold for natural-product libraries screened in step 05: the covariate shift between synthetic calibration chemistry and natural products breaks exchangeability, so intervals attached to library compounds are uncalibrated extrapolations and must be described as such wherever they are reported

- insufficient NP-like calibration molecules (14 < 60): no NP-specific conformal stratum can be fitted from ChEMBL alone. Intervals transferred to a natural-product library therefore carry NO coverage guarantee, only an extrapolation of the synthetic-chemistry error model.

- The NP-like held-out subset is small (n = 20 for FGFR4); its metrics carry wide confidence intervals and support a direction, not a point claim.

- The applicability-domain threshold is non-binding for FGFR4: tau sits at the lowest similarity bin: at a 1.65 pIC50 width budget every calibration bin qualifies, so the applicability domain excludes nothing and the gate is vacuous for this target. With this calibration set, bins wide enough for a stable 90% quantile are too wide to resolve a threshold. See tau_by_budget for what a tighter budget would give. Any in-domain hit count for this target is therefore a count of the whole library, not of a calibrated-reliable subset.

- Top-1 redocking fails the 2 A criterion for EGFR, FGFR4, so pose-level conclusions in those pockets are not supported.

- Docking does not enrich actives over the negative class for EGFR: no metric clears its label-shuffled null. A docking score therefore carries no validated discriminative weight in those pockets and cannot be used to support or filter a prediction there.

- Where docking does enrich, it does so inconsistently: no pocket shows both global separation and early recognition, so the enrichment control licenses a direction to confirm on a larger set rather than a general-purpose structure-based filter.

- AutoDock Vina scores are empirical and receptor-specific. They are not transferable binding free energies, so a weaker score in one pocket than another is not evidence of selectivity. Compare within a target (see the <T>_vina_z columns), and treat any score weaker than -4.0 kcal/mol as a placement failure to investigate, not as a measured steric clash.

- No experimental validation. Every candidate is a computational prediction and the paper claims nothing beyond that.
