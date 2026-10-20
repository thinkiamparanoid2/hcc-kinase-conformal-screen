# Predictions

This directory holds the screen's own outputs:

| File | Written by | What it is |
|---|---|---|
| `consensus_hits.csv.gz` | step 05 | Every library compound that clears the screen, with per-target prediction, interval lower bound, max-Tanimoto similarity, applicability-domain flag, scaffold novelty and pan-kinase alert. Step 05 writes the plain `consensus_hits.csv`; only the gzipped copy is tracked, because the uncompressed file is 213 MB and GitHub's hard limit is 100 MB. Steps 06 and 07 read the plain CSV, so gunzip it after a fresh clone. |
| `consensus_hits_annotated.csv` | step 07 | The shortlist, with QED, PAINS and chemotype annotation. |
| `consensus_hits_docked.csv` | step 06 | The shortlist with Vina scores per pocket, replicate SDs and within-target z-scores. |
| `<TARGET>_screen_full.csv.gz` | step 05 | Predictions for the whole library, before any filter. Not tracked. |

## Why the shortlist is so large

`tau` came out non-binding for FGFR4 and EGFR (step 04 found every calibration
similarity bin inside the width budget), so `in_AD` admits the whole library on
those two targets and the shortlist rule `n_targets_inAD_active >= 1` selects
645,851 of 690,177 parents -- 94% of the library. The in-AD counts for those two
targets are therefore not a calibrated-reliable subset of anything, which
`screen_summary.json` records per target as `in_AD_caveat`.

The criterion that does discriminate is the conformal interval's lower bound
clearing the activity threshold (`n_targets_interval_active`): 3,271 compounds
across all three targets, two to four orders of magnitude fewer. Read that as
the hit count, not the in-AD column.

**The files that were here before** were
produced before the pipeline was rebuilt: they named hits by a rule this work
withdrew (`confident_active`, now `interval_lower_bound_active`), carried no
scaffold-level novelty, and were not deduplicated to 2D parents, so a
stereoisomer counted as several hits. They were removed rather than left in
place, because a candidate list is the one file a reader will take at face
value. They remain in the git history.

Regenerate them with the COCONUT library in place:

```bash
python scripts/05_screen_coconut_library.py   # writes consensus_hits.csv
python scripts/07_generate_report.py          # writes consensus_hits_annotated.csv
python scripts/06_docking_validation.py       # writes consensus_hits_docked.csv
```

See [`../../data/README.md`](../../data/README.md) for how to obtain the library.
Step 06 refuses to dock a shortlist written by the older step 05, so a stale file
here cannot quietly re-enter the results.
