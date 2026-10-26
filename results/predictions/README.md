# Predictions

This directory holds the screen's own outputs:

| File | Written by | What it is |
|---|---|---|
| `consensus_hits.csv.gz` | step 05 | The shortlist: 3,718 compounds, with per-target prediction, interval lower bound, max-Tanimoto similarity, applicability-domain flag, scaffold novelty and pan-kinase alert. Step 05 writes the plain `consensus_hits.csv` (1.3 MB); only the gzipped copy is tracked. Steps 06 and 07 read the plain CSV, so gunzip it after a fresh clone. |
| `consensus_hits_annotated.csv` | step 07 | The shortlist, with QED, PAINS and chemotype annotation. |
| `consensus_hits_docked.csv` | step 06 | The shortlist with Vina scores per pocket, replicate SDs and within-target z-scores. |
| `<TARGET>_screen_full.csv.gz` | step 05 | Predictions for the whole library, before any filter. Not tracked. |

## What the shortlist rule is, and why `in_AD` is not part of it

A compound is on the shortlist if its conformal interval's lower bound clears the
activity threshold on at least one target (`n_targets_interval_active >= 1`), or
if it is a rediscovery of known chemistry and therefore belongs to step 06's
control (`rediscovery_any`). That gives 3,718 of 690,177 parents.

`in_AD` is deliberately not a criterion. `tau` came out non-binding for FGFR4 and
EGFR (step 04 found every calibration similarity bin inside the width budget), so
`in_AD` admits the whole library on those two targets: including it selected
645,851 parents, 94% of the library, in a 213 MB file whose size was a property
of the calibration rather than of the compounds. The in-AD counts for those two
targets are not a calibrated-reliable subset of anything, which
`screen_summary.json` records per target as `in_AD_caveat`. Across all three
targets 3,271 compounds are interval-active, two to four orders of magnitude
fewer than the in-AD column: read that as the hit count.

Only 76 of the 523 rediscovery compounds are interval-active, so the rediscovery
clause is what keeps step 06's control drawing from all 523 rather than from the
subset the model is already confident about, which would quietly bias it.

Every count in `screen_summary.json` is computed on the unfiltered frame, so
narrowing this file changes no reported number.

**The files that were here before** were produced before the pipeline was
rebuilt: they named hits by a rule this work withdrew (`confident_active`, now `interval_lower_bound_active`), carried no
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
