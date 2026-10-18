# Limitations register

What this work does not establish, and what was done about each item. Section 5 of
the auto-generated report restates the live ones from the pipeline's own
measurements; this file is the stable register, including items that are now fixed
and are kept here so reviewers can see the reasoning.

Status: **fixed** (addressed in code), **mitigated** (reduced and reported),
**open** (stated as a limitation, not solved).

## Data and labels

| # | Limitation | Status | What was done |
|---|---|---|---|
| D1 | Only `=`, `<`, `<=` relations were fetched, so there was no inactive class and ~77% of labels were active | fixed in code, **open in the committed data** | Step 01 now requests all five relations and labels each row. The committed CSVs predate this, so `curation_*.json` records `censored_records_available: false` and every classification metric is reported with its prevalence. Re-run step 01 with network access to close it. |
| D2 | Censored values entered the regression as point labels | fixed | `relation` / `censored` / `label_kind` survive aggregation; `load_training` returns uncensored rows only by default. Censored rows remain available for classification and docking enrichment, where a bound is informative. |
| D3 | No structure standardisation; ~2% of each set was salt/free-base or tautomer duplicates with different Murcko scaffolds, a scaffold-split leak | fixed | Salt-strip → neutralise → canonical tautomer, deduplicate on parent InChIKey. 62–74 duplicate parents collapsed per target. |
| D4 | Replicate spread was filtered but never reported, so the label-noise ceiling on R² was unknown | mitigated | Step 01 writes replicate counts, SD and the implied R² ceiling. Unavailable when re-curating a pre-aggregated CSV, and reported as unknown in that case rather than as 1.0. |
| D5 | Only IC50; Ki and Kd binding data are excluded | open | A deliberate scope choice. Including them would need an assay-type covariate. |

## Splitting and modelling

| # | Limitation | Status | What was done |
|---|---|---|---|
| S1 | The scaffold split skipped oversized groups, pushing large scaffold families into train and enriching the test set in rare scaffolds | fixed | Two explicit schemes: a deterministic largest-first Bemis-Murcko split for headline numbers, and an unbiased random-group split for variance. A unit test asserts large families can reach the test set. |
| S2 | Repeated-seed SDs were presented as if they were sampling variance | mitigated | Reported as split variance, with the paired Wilcoxon labelled as a statement about split noise. |
| S3 | The benchmark favoured RF-ecfp while the screen ran on XGBoost | fixed | Step 03 writes `model_selection.json`; steps 04 and 05 read it. |
| S4 | Hyperparameters were fixed with no stated provenance | mitigated | Declared as untuned library defaults in `HYPERPARAMETER_NOTE`, carried into the report. Nested CV is the next step if tuning is added. |
| S5 | No external or prospective validation | mitigated | A temporal holdout (train on older literature) runs whenever publication years are available. Still not prospective validation. |
| S6 | No uncertainty on headline numbers | fixed | Percentile bootstrap CIs on R², AUC, coverage and `tau`. |

## Conformal prediction and applicability domain

| # | Limitation | Status | What was done |
|---|---|---|---|
| C1 | Exchangeability between synthetic calibration data and a natural-product deployment set does not hold, so the 90% guarantee did not transfer | mitigated, **open by nature** | Stated explicitly in `conformal_validity`, carried into every downstream file and the report. Coverage is measured separately on NP-like held-out molecules, and an NP-specific conformal stratum is fitted when enough NP-like calibration data exists. Closing it properly needs a natural-product bioactivity set (NPASS, or the LOTUS ∩ ChEMBL overlap) as a true NP calibration set. |
| C2 | The NP domain-shift claim rested on 20–43 molecules | mitigated | Subsets under 20 report no metric at all; the rest carry bootstrap CIs. The claim is stated as a direction, not a point estimate. |
| C3 | `MIN_CAL_PER_BIN = 30` makes a 90% quantile equal to the maximum residual | fixed | Raised to 150; each bin reports its conformal level, its own coverage and a CI. |
| C4 | `tau` was derived by a top-down scan that assumed monotone half-widths | fixed | Isotonic (non-increasing) fit, bootstrap CI, and a `tau ± 0.1` sensitivity table. |
| C5 | The library sits almost entirely in the lowest, least-calibrated similarity bin | open | Reported as a finding: the AD-admissible fraction of the library is small, and that is the measurement the paper makes. |

## Screening

| # | Limitation | Status | What was done |
|---|---|---|---|
| V1 | Stereoisomer records inflated the compound count and duplicated predictions | fixed | Deduplicated on parent InChIKey before screening; both counts reported. |
| V2 | `confident_active` implied calibrated confidence the intervals do not carry on this library | fixed | Renamed `interval_lower_bound_active`, with the caveat attached in the summary and the report. |
| V3 | Novelty was measured only by Tanimoto, so a staurosporine analogue at 0.7 counted as novel | fixed | Scaffold-level novelty (`known_scaffold`) plus pan-kinase frequent-hitter chemotype flags. |
| V4 | Consensus across three ATP-site kinases was presented as confirmation | fixed | Cross-target prediction correlation is computed and reported; consensus is described as a polypharmacology axis and a selectivity liability. |
| V5 | No baseline — plain similarity search might recover the same hits | fixed | A nearest-neighbour baseline ranks the same library; rank agreement with the model is reported per target. |

## Docking

| # | Limitation | Status | What was done |
|---|---|---|---|
| K1 | Best-of-9 RMSD was reported as the redocking control while top-1 failed the 2 Å criterion in two of three pockets | fixed | Top-1 is the control with an explicit pass/fail; best-of-N is labelled a sampling diagnostic. |
| K2 | RMSD ignored atom-order matching and molecular symmetry | fixed | RDKit `GetBestRMS`. Returns `None` on a failed match rather than a misleading number. |
| K3 | Enrichment used 12 vs 12, taking the most potent against the weakest | fixed | Random within-class sampling, MW-matched, n = 50 per class, with bootstrap CIs, EF1%, EF5% and BEDROC. |
| K4 | A fixed 22 Å box could not place large ligands; near-zero scores were read as steric clash | fixed | Boxes sized to the ligand; implausibly weak scores flagged `placement_suspect` instead of interpreted. |
| K5 | Cross-receptor kcal/mol comparison underpinned the DFG-in/DFG-out claim | fixed | Within-target z-scores; the comparability caveat travels with the data. |
| K6 | One seed, no replicates, inconsistent exhaustiveness | fixed | `N_SEEDS = 3`, exhaustiveness 32 throughout, mean ± SD reported. |
| K7 | Minimal receptor preparation: no explicit waters, no His/Asn/Gln flip check, no rescoring | open | Would need PDB2PQR/PROPKA preparation and MM-GBSA or short MD rescoring of the final poses. |
| K8 | Vina treats macrocyclic rings as rigid and uses a single input conformer | open | Affects macrocyclic candidates specifically; stated where those compounds are discussed. |

## Scope

| # | Limitation | Status |
|---|---|---|
| G1 | No experimental validation of any predicted compound | open |
| G2 | No ADMET, developability or natural-product sourcing assessment beyond PAINS and QED | open |
| G3 | Rigid-receptor docking captures no induced fit | open |
