# Where this work stands, and what to do next

Written 2026-10-20, after re-running the whole pipeline on a Windows machine with
the two external inputs the cloud container could not reach. Read this first when
you come back to the project.

Companion documents: [`limitations.md`](limitations.md) is the formal limitation
register; [`REPRODUCE.md`](REPRODUCE.md) is how to run things. This file is the
state of play and the plan.

---

## 1. What changed in this run

Branch `claude/project-thread-pv3x3o`, PR #4.

**The inactive class now exists.** The committed ChEMBL CSVs were fetched before
step 01 requested all five `standard_relation` values, so they held no `>`
records — no compound that was tested and found inactive. Every classification
and enrichment number was therefore a potent-vs-weak comparison wearing an
active-vs-inactive label. Re-fetched all three targets:

| Target | parents | exact | censored-inactive | R² ceiling from label noise |
|---|---|---|---|---|
| EGFR | 12,484 | 9,821 | **1,659** | 0.949 |
| FGFR4 | 4,819 | 3,609 | **1,078** | 0.993 |
| VEGFR2 | 11,269 | 9,228 | **1,666** | 0.996 |

That closes limitation D1. Everything downstream was re-run on this data at the
full docking protocol — `protocol_note` in `docking_validation.json` is `null`,
which is the machine-checked statement that no budget was reduced.

Environment: AutoDock Vina 1.2.5, Open Babel 3.1.0, 16 cores. `pytest -q` is 36
passed.

---

## 2. The numbers, all in one place

### Model benchmark (step 03) — mean R² over 10 random scaffold splits

| Target | RF-ecfp | RF-all | XGB-all | kNN-Tanimoto | XGB-scramble |
|---|---|---|---|---|---|
| EGFR | **0.572** | 0.556 | 0.545 | 0.528 | −0.099 |
| FGFR4 | **0.501** | 0.483 | 0.484 | 0.411 | −0.182 |
| VEGFR2 | **0.514** | 0.502 | 0.496 | 0.468 | −0.088 |

RF-ecfp selected, same family as before. The scramble control is negative on all
three, so the models learned real signal rather than memorising noise. But note
how close `kNN-Tanimoto` is — see §4.

Temporal holdout (train on older literature, test on newest 20% by publication
year) collapses to roughly R² 0.08–0.27 depending on target and model. That is
the honest estimate of performance on genuinely new chemistry.

### Conformal calibration (step 04)

| Target | empirical coverage (nominal 0.90) | τ | binding? | test in AD |
|---|---|---|---|---|
| EGFR | 0.858 [0.840, 0.876] | 0.0 | **no** | 100% |
| FGFR4 | 0.871 [0.841, 0.900] | 0.0 | **no** | 100% |
| VEGFR2 | 0.892 [0.876, 0.909] | 0.5 | yes | 68.6% |

EGFR's coverage CI excludes 0.90 — the guarantee is not being met even on
exchangeable held-out ChEMBL scaffolds.

**No NP-specific conformal stratum could be fitted for any target:** 9 (EGFR),
22 (FGFR4), 28 (VEGFR2) natural-product-like calibration molecules against a
minimum of 60. Every interval attached to a COCONUT compound is therefore an
extrapolation of a synthetic-chemistry error model, with no coverage guarantee.
This is stated in the output and must stay stated.

### Screen (step 05)

738,827 library records → **690,177 unique 2D parents** (8 unparsable).

| Target | τ | in AD | in-AD active | **interval-LB active** | …on novel scaffolds | median sim to train |
|---|---|---|---|---|---|---|
| FGFR4 | 0.0 non-binding | 100% | 643,633 | **3,220** | 2,977 | 0.182 |
| EGFR | 0.0 non-binding | 100% | 123,659 | **34** | 1 | 0.217 |
| VEGFR2 | 0.5 binding | 0.30% (2,078) | 977 | **111** | 7 | 0.206 |

Across any target: 3,271 interval-active, versus 645,850 in-AD-active (94% of the
library). **The interval-lower-bound column is the defensible hit count.**

Model-vs-nearest-neighbour-baseline rank agreement over in-domain compounds:
FGFR4 0.24, EGFR 0.32, VEGFR2 0.72.

### Docking validation (step 06)

Redocking control — can Vina re-find a pose already known from crystallography?

| Target | top-1 RMSD | 2 Å | best-of-9 | superimposed | best pose rank | ΔG gap to best pose |
|---|---|---|---|---|---|---|
| EGFR (erlotinib) | 8.004 Å | **FAIL** | 1.283 Å | 1.78 Å | 4 | −0.166 kcal/mol |
| FGFR4 (6LF) | 3.334 Å | **FAIL** | 1.811 Å | 1.443 Å | 7 | −0.823 kcal/mol |
| VEGFR2 (sorafenib) | 0.577 Å | PASS | 0.577 Å | 0.265 Å | 1 | — |

EGFR's correctly placed pose exists at rank 4 and is only 0.166 kcal/mol from the
top-ranked one, so the scoring function cannot separate them. Quoting the 1.78 Å
superimposed figure would present a wrong-subpocket pose as a success.

Enrichment control — 50 actives vs 50 MW-matched true `>` inactives per target,
random within class, not potency extremes:

| Target | AUC | 95% CI | was (no true inactives) | MWU p | exceeds permutation null |
|---|---|---|---|---|---|
| FGFR4 | **0.696** | [0.592, 0.791] | 0.623 | 0.0004 | ROC-AUC, EF5%, BEDROC |
| EGFR | **0.444** | [0.327, 0.560] | 0.416 | 0.83 | **nothing** |
| VEGFR2 | **0.693** | [0.585, 0.792] | 0.536 | 0.0005 | ROC-AUC, BEDROC |

A real inactive class helps two pockets materially. EGFR is below chance, its CI
spans 0.5, and it beats no null on any metric. EGFR is also the pocket that fails
redocking worst — placement and ranking fail together there.

Candidates: 23 novel + 8 rediscovery controls docked in all three pockets, no
placement failures. Mean novel affinity FGFR4 −9.63, EGFR −9.98, VEGFR2 −8.38
kcal/mol. **Not comparable between pockets** — Vina scores are receptor-specific.

---

## 3. Bugs this run exposed, all fixed with tests

1. **Step 03's checkpoint could not tell two datasets apart.** `done_keys()` keyed
   on `(target, seed)`, so re-running after the data changed reported all 10 seeds
   per target as "already done" and reused rows computed from the old data;
   `summarize()` would have picked the downstream model family from a mixture of
   two datasets. Nothing crashed. Rows now carry a `data_sha` and `prune_stale()`
   evicts mismatches. 210 stale rows were dropped and every seed recomputed.

2. **Step 06's candidate staleness guard could never pass.** It required columns
   `interval_lower_bound_active` and `known_scaffold` — keys of
   `screen_summary.json` and of a frame internal to step 05, never columns of the
   CSV it reads (`n_targets_interval_active`, `known_scaffold_all_targets`). So
   candidate docking had never once run, on any machine. The same function already
   used the correct name 30 lines below, and the test asserting the old name
   encoded the bug.

3. **Step 06 died instantly on Windows.** `_mol_from_pdbqt_model()` deleted the
   scratch SDF while RDKit's `SDMolSupplier` still held the handle — permitted on
   POSIX, refused on Windows. The supplier is now closed first and a failed delete
   of that scratch file no longer aborts a multi-hour run.

4. **`standardise_series()` was the pipeline's real bottleneck.** Serial, ~101 ms
   per molecule, which is **~21 hours** for the 738k library and made step 05 look
   hung on one core. Now parallel across processes, **~1.5 h measured**. A test
   asserts the parallel path returns exactly what the serial one does, in order,
   including `None` for unparsable SMILES. Steps 01 and 02 use the same path.

---

## 4. Honest assessment

**The validation machinery is stronger than the science it validates, and that is
the most interesting thing here.**

What genuinely works: the models learned real signal (scramble control negative);
R² 0.50–0.57 on unseen scaffolds is respectable for QSAR; VEGFR2 passes every
control end to end; the enrichment control improved exactly as predicted once a
real inactive class existed.

What does not hold up:

- **The headline methodological contribution is inert on two of three targets.**
  The selling point is a domain-calibrated Mondrian conformal applicability
  domain, and τ = 0 for FGFR4 and EGFR, admitting the whole library. See §5 —
  this is a parameter choice, not a broken method.
- **The model barely beats a trivial baseline.** A few points of R² over "find the
  most similar known compound and guess its potency." And where the model's
  ranking differs from that baseline (agreement 0.24 on FGFR4) nothing
  establishes that it differs in the right direction.
- **Prospective generalisation is poor.** Temporal-holdout R² ≈ 0.08–0.27.
- **Coverage undershoots its own nominal level** before any natural-product shift.
- **The natural-product intervals have no guarantee at all**, and the pipeline
  says so itself. This is the most serious weakness.
- **Structure-based validation validates one pocket of three.**
- **Nothing has been tested experimentally.** "Hit" currently means a model and a
  docking program agreed.

**The one coherent positive story is VEGFR2:** τ binds and admits 0.30% of the
library, redocking passes at 0.577 Å, enrichment is 0.693 with a CI excluding 0.5
and exceeding its permutation null. 111 interval-active hits, 7 on novel
scaffolds. Narrow, but real and defensible.

**What should not be presented:** 3,220 FGFR4 "hits". That target's AD filter is
inert and its pocket fails redocking.

---

## 5. The recoverable finding: τ depends on a parameter nobody reports

τ is chosen against a *width budget* — how much prediction error you are willing
to tolerate. The pipeline used **1.65 log units**, roughly 45-fold uncertainty in
IC50, which is permissive enough that nothing is ever excluded. Step 04 prints τ
at every budget, and the table is the finding:

| width budget (log units) | FGFR4 | EGFR | VEGFR2 |
|---|---|---|---|
| 1.00 | — | 0.7 | — |
| **1.25** | **0.7** | **0.7** | **0.7** |
| 1.50 | 0.0 | 0.6 | 0.5 |
| 1.65 (used) | 0.0 | 0.0 | 0.5 |
| 2.00 | 0.0 | 0.0 | 0.0 |

**At a 1.25-log-unit budget the applicability domain binds on all three targets.**
The machinery was never broken; it was asked to enforce a bar low enough that
everything cleared it. Tightening the budget *strengthens* the control — fewer
compounds admitted, stricter claims.

Two things this does **not** fix: EGFR's docking failure (a hard, flexible,
water-mediated pocket — a known-difficult case in the literature, not a mistake
of yours), and the natural-product calibration gap (needs more NP bioactivity
data, not a parameter).

This has **not** been run. It is the first experiment to do, and it should be
reported as a documented sweep rather than quietly swapping the number.

---

## 6. What to do next, in priority order

### Do first (hours)

1. **Re-run steps 04 → 07 at a 1.25-log-unit width budget.** ~2 h. On the evidence
   in §5 this likely gives a binding applicability domain on all three targets,
   which changes the paper's story. Report the sweep, not just the chosen value.
2. **Make the τ-versus-width-budget sweep a headline figure.** This is the core
   methodological contribution, currently buried in a log line.
3. **Plot the similarity distribution between library and training set.** Median
   0.18–0.22 frames every other result: the library is almost entirely outside
   the training chemistry.
4. **Keep the kNN-Tanimoto baseline prominent.** Reviewers respect it being shown.

### Do next (weeks)

5. **Run the pipeline on 10–15 kinases instead of 3.** Mostly compute; the code
   exists. "Applicability domains can be vacuous" is an anecdote at n=3 and a
   finding at n=15. Highest credibility gain per unit effort.
6. **Add two or three competing AD methods** — distance-to-model, bootstrap
   variance, inductive vs Mondrian conformal — and show whether they agree on
   which compounds are out of domain. A comparative AD study stands on its own.
7. **Add property-matched decoys (DUD-E style)** alongside the ChEMBL inactives,
   so the docking control rests on two independent negative classes.
8. **Report the temporal split prominently** — train pre-2015, test post-2020.
   `first_year` is already stored.

### The highest-value step that needs no laboratory

9. **Measure interval coverage on real natural products.** NPASS, and the
   natural-product subsets of ChEMBL and BindingDB, carry NP–kinase activity
   pairs. A few hundred would let you *measure* whether the conformal intervals
   hold on natural products instead of conceding that you cannot. That converts
   the single biggest limitation from a caveat into a result, and it is the
   control the work is currently missing.

### The step that changes the paper's category

10. **Test a handful of VEGFR2 hits experimentally.** 7 novel-scaffold candidates
    from the one target where every control passed. A commercial kinase assay on
    5–10 purchasable compounds, or a collaboration with anyone running kinase
    assays. One confirmed hit turns a methods paper into a discovery paper backed
    by unusually rigorous validation.

### Do not spend time on

A bigger library. A graph neural network or transformer. More targets in the
*screen*. None of these touch the actual weakness, which is that the uncertainty
estimates cannot presently be shown to mean anything on natural products.

---

## 7. Publication view

As a natural-product hit-discovery paper this would not survive competent review:
no experimental validation, an AD filter that admits the whole library on two
targets, docking that fails its own control on two targets. Venues that would
accept that framing are not worth publishing in.

As a methods and negative-results paper it is credible, and the rigour is its
strongest asset. The frame that works: *applicability domains can be vacuous, and
conformal intervals calibrated on synthetic chemistry do not transfer to natural
products — here is a pipeline strict enough to detect both, with VEGFR2 as the
positive control showing it can work end to end.*

Realistic venues for that framing: Journal of Cheminformatics, Molecular
Informatics, RSC Digital Discovery. JCIM becomes plausible with item 5. If item 10
lands, aim higher.

**Expect this objection and prepare for it:** a reviewer can answer the AD-vacuity
finding with "you chose a bad width budget — that is your error, not the field's."
Items 5 and 6 are the defence, together with the observation that published work
rarely states its AD threshold or width budget at all, which is itself the point.

---

## 8. Open decisions and known state

- **Shortlist rule changed in this run.** `consensus_hits.csv` was 213 MB
  (645,851 rows) because the rule was `n_targets_interval_active >= 1 OR
  n_targets_inAD_active >= 1`, and the in-AD clause admits the whole library where
  τ is non-binding. It is now `n_targets_interval_active >= 1 OR rediscovery_any`
  — about 3,718 rows. The rediscovery clause matters: only 76 of 523 rediscovery
  compounds are interval-active, and without it step 06's rediscovery control
  would draw only from compounds the model is already confident about, which would
  bias it. Every count in `screen_summary.json` is computed on the unfiltered
  frame, so narrowing the file changed no reported number.
- **Figures and the report** are regenerated by step 07 into
  `results/reports/results_report.md` and `results/figures/`.
- **Not tracked in git, by design:** `models/*_final.joblib` (too large,
  deterministic — regenerate with step 04), `results/docking_poses/` (Open Babel
  places hydrogens non-deterministically, so a tracked pose is stale the moment
  step 06 runs), `data/features/*.csv`, `data/ligands_phytochemical/*.csv`, and
  the per-target `*_screen_full.csv.gz` dumps.
- **The receptor-preparation gap is still open** (limitation K7): `obabel -xr
  -p 7.4` adds polar hydrogens at an arbitrary rotation, and those HD atoms are
  exactly the hydrogen-bond donors the scoring function reads. Docking scores are
  reproducible within a run and only approximately across runs. Fixing it properly
  means preparing each receptor once with a deterministic tool (PDB2PQR/PROPKA)
  and committing that.
- **Windows note:** steps 01–05 and 07 need only `pip install -r
  requirements.txt`. Step 06 needs `vina` and `obabel`; the Windows Vina binary is
  on the AutoDock-Vina GitHub releases page and `obabel` comes from the
  `openbabel-wheel` pip package. Point `VINA_BIN` / `OBABEL_BIN` at them.
