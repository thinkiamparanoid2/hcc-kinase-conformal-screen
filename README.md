# Domain-Calibrated QSAR and Structure-Based Screening of Natural Products Against HCC Kinases

Which natural products are credible inhibitor candidates against the receptor
tyrosine kinases that drive hepatocellular carcinoma, and how much confidence does
the evidence actually support?

This repository screens the COCONUT natural-product library against three HCC
receptor tyrosine kinases - **EGFR** (`1M17`), **FGFR4** (`5JKG`) and **VEGFR2**
(`4ASD`) - with scaffold-split QSAR, Mondrian split-conformal calibration,
scaffold-level novelty filtering and orthogonal AutoDock Vina validation, and
reports the surviving compounds as prioritised candidates.

The conformal applicability domain is what makes that candidate list defensible
rather than decorative. It states, per compound, whether the model is being asked
to predict inside the chemistry it was calibrated on, and that gate is applied
before any compound is called a candidate. Its limits are measured rather than
assumed: the coverage guarantee established on held-out ChEMBL scaffolds does not
transfer to natural products, so the intervals carried into the screen are an
extrapolated error model and are labelled that way wherever they appear.

> **The candidate list is not in this repository yet.** Steps 05 and 06 - the
> full-library screen and the docking validation - have not been run since the
> pipeline was rebuilt, so those sections of the report are withheld rather than
> filled with numbers the current code never produced. Run them (see
> [`docs/REPRODUCE.md`](docs/REPRODUCE.md)) to produce the results this framing
> promises.

---

## Results

All numbers live in **[`results/reports/results_report.md`](results/reports/results_report.md)**,
which is regenerated from `results/metrics/` and `results/predictions/` by step 07.
Nothing is transcribed into this README, so nothing here can drift away from the
data. Figures are written to `results/figures/`.

The report covers:

| Section | Content |
|---|---|
| 0 | Data curation: censored labels, standardisation, label-noise ceiling |
| 1 | Scaffold-split benchmark with bootstrap CIs, temporal holdout, model selection |
| 2 | Conformal calibration, applicability domain, natural-product transfer |
| 3 | Full-library screen with scaffold-level novelty and a nearest-neighbour baseline |
| 4 | Redocking (top-1 RMSD), enrichment (AUC with CI, EF, BEDROC), candidate docking |
| 5 | Limitations, generated from what the pipeline actually measured |

---

## Method, and the assumptions it rests on

**Labels.** Human binding-assay IC50 from ChEMBL. `standard_relation` is preserved:
`<` / `<=` bound potency from below, `>` / `>=` from above, and only uncensored
records are used as regression targets. Structures are salt-stripped, neutralised
and tautomer-canonicalised, then deduplicated on the parent InChIKey.

**Splitting.** Generic Bemis-Murcko scaffold splits, strictly disjoint. A
deterministic canonical split carries the headline numbers with percentile
bootstrap CIs; repeated random scaffold splits measure split variance, which is not
the same thing as generalisation to new chemistry and is labelled as such.

**Model.** The family used downstream is chosen by the benchmark
(`results/metrics/model_selection.json`), not fixed in advance. Hyperparameters are
library defaults and were not tuned against any reported split.

**Calibration.** Mondrian split-conformal intervals at 90%, with the taxonomy given
by max-Tanimoto to the training set. The domain threshold `tau` is read off an
isotonic fit of half-width against similarity and reported with a bootstrap CI and
a sensitivity analysis, because a threshold that gates an entire screen should not
rest on one noisy bin.

> **Validity.** Split-conformal coverage holds under exchangeability. Held-out
> ChEMBL scaffolds satisfy it; a natural-product library does not. Intervals
> attached to library compounds in step 05 are therefore an error model calibrated
> on synthetic chemistry and extrapolated — they are labelled that way in every
> output file, and `interval_lower_bound_active` is never described as a
> guaranteed-coverage hit.

**Screening.** The library is deduplicated to 2D parents before screening, since
the features are achiral. Hits are scored for novelty at the scaffold level, not
only by similarity, and flagged for pan-kinase frequent-hitter chemotypes. A
similarity-search baseline ranks the same library, so the model's added value is
measured rather than assumed.

**Docking.** Top-1 symmetry-corrected RMSD is the redocking control; best-of-N is
reported beside it as a sampling diagnostic. The enrichment control uses randomly
sampled, property-matched actives and inactives with bootstrap CIs, EF1%, EF5% and
BEDROC. Boxes are sized to the ligand, replicate seeds give error bars, and Vina
scores are compared within a target only.

---

## Quick start

```bash
git clone https://github.com/thinkiamparanoid2/hcc-kinase-conformal-screen.git
cd hcc-kinase-conformal-screen

# Conda (includes AutoDock Vina and Open Babel for step 06)
conda env create -f environment.yml && conda activate hcc-screen

# or: pip install -r requirements.txt   (steps 01-05 and 07 only)

pytest -q                          # 24 unit tests, no external data needed
python scripts/run_pipeline.py --list
python scripts/run_pipeline.py     # full pipeline, checkpointed and resumable
```

Steps that need external inputs are skipped with a named reason rather than
failing: step 05 needs the COCONUT library (see [`data/README.md`](data/README.md))
and step 06 needs `vina` and `obabel` on `PATH` (or `VINA_BIN` / `OBABEL_BIN`).

See [`docs/REPRODUCE.md`](docs/REPRODUCE.md) for the full reproduction path and
[`COMMANDBOOK.md`](COMMANDBOOK.md) for per-step commands.

---

## Repository layout

```text
├── data/
│   ├── ligands_chembl/          # Curated IC50 datasets (tracked)
│   ├── ligands_phytochemical/   # COCONUT library (not tracked; see data/README.md)
│   ├── targets/                 # Co-crystal PDB and prepared PDBQT receptors
│   └── features/                # Regenerated by step 02 (not tracked)
├── models/                      # Final models + conformal calibration per target
├── results/
│   ├── figures/  metrics/  predictions/  docking_poses/  reports/
├── scripts/                     # common.py + steps 01-07 + run_pipeline.py
├── tests/                       # pytest suite (no network, no docking binaries)
├── docs/                        # Reproduction guide, limitations register
├── environment.yml  Dockerfile  requirements.txt
└── CITATION.cff  LICENSE
```

---

## Known limitations

Maintained in [`docs/limitations.md`](docs/limitations.md) and regenerated into
section 5 of the report from the pipeline's own measurements. The load-bearing
ones: conformal validity does not transfer to natural products; the NP-like
evaluation subset available from ChEMBL is small; there is no experimental
validation; and rigid-receptor Vina scores are not binding free energies.

## Licence and citation

MIT — see [`LICENSE`](LICENSE). Citation metadata in [`CITATION.cff`](CITATION.cff).
Archive a release and add its DOI before submission.
