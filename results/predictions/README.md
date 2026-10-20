# Predictions

This directory holds the screen's own outputs:

| File | Written by | What it is |
|---|---|---|
| `consensus_hits.csv` | step 05 | The shortlist: every compound whose conformal interval's lower bound clears the activity threshold on at least one target, with per-target prediction, interval lower bound, max-Tanimoto similarity, applicability-domain flag, scaffold novelty and pan-kinase alert. 3,271 rows. |
| `consensus_hits_annotated.csv` | step 07 | The shortlist, with QED, PAINS and chemotype annotation. |
| `consensus_hits_docked.csv` | step 06 | The shortlist with Vina scores per pocket, replicate SDs and within-target z-scores. |
| `<TARGET>_screen_full.csv.gz` | step 05 | Predictions for the whole library, before any filter. Not tracked. |

## What counts as a hit, and why not the applicability domain

`tau` came out non-binding for FGFR4 and EGFR (step 04 found every calibration
similarity bin inside the width budget), so `in_AD` admits the whole library on
those two targets. A shortlist rule of `n_targets_inAD_active >= 1` therefore
selected 645,850 of 690,177 parents -- 94% of the library, a 213 MB file naming
almost everything screened. Its size measured the calibration, not the compounds.

The shortlist is now the criterion that discriminates: the conformal interval's
lower bound clearing the activity threshold on at least one target
(`n_targets_interval_active`), which is 3,271 compounds. The in-AD counts are
still reported in `screen_summary.json`, measured over the whole library rather
than over the shortlist, because the fraction they give is what shows the domain
filter is not biting; each target also carries that as `in_AD_caveat`. Read the
interval column as the hit count, never the in-AD one.

Dropping the in-AD arm removes only the tail of the file. The sort puts every
interval-active compound above every merely in-AD one, so the annotated top 200,
the docked candidate set and the report's tables were already drawn from these
3,271 rows and are unchanged.

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
