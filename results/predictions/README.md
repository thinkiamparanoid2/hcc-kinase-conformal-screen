# Predictions

This directory holds the screen's own outputs:

| File | Written by | What it is |
|---|---|---|
| `consensus_hits.csv` | step 05 | Every library compound that clears the screen, with per-target prediction, interval lower bound, max-Tanimoto similarity, applicability-domain flag, scaffold novelty and pan-kinase alert. |
| `consensus_hits_annotated.csv` | step 07 | The shortlist, with QED, PAINS and chemotype annotation. |
| `consensus_hits_docked.csv` | step 06 | The shortlist with Vina scores per pocket, replicate SDs and within-target z-scores. |
| `<TARGET>_screen_full.csv.gz` | step 05 | Predictions for the whole library, before any filter. Not tracked. |

**The directory is currently empty of results.** The files that were here were
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
