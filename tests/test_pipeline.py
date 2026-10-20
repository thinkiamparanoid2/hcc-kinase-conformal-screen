"""
Smoke and correctness tests for the pipeline's load-bearing logic.

These run in seconds on synthetic data and need no ChEMBL access, no COCONUT
library and no docking binaries, so they are safe to run in CI. They cover the
parts where a silent bug would change a published number: split disjointness,
conformal coverage, the tau rule, similarity search, and metric edge cases.

    pytest -q
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "scripts"))
import common as C


# ----------------------------------------------------------------------------
# Structure standardisation
# ----------------------------------------------------------------------------

def test_salt_and_free_base_collapse_to_one_parent():
    """A hydrobromide salt and its free base must give the same InChIKey.

    This is the leak the earlier curation had: deduplicating on raw canonical
    SMILES kept both, and they receive different Murcko scaffolds, so one molecule
    could land on both sides of a scaffold split.
    """
    _, key_salt = C.standardise_mol("Br.CC(Nc1ncnc2[nH]ccc12)c1ccccc1")
    _, key_free = C.standardise_mol("CC(Nc1ncnc2[nH]ccc12)c1ccccc1")
    assert key_salt is not None and key_salt == key_free


def test_standardise_rejects_garbage():
    assert C.standardise_mol("not-a-molecule") == (None, None)


# ----------------------------------------------------------------------------
# Scaffold splitting
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("scheme", ["deterministic", "random"])
def test_scaffold_split_is_disjoint(scheme):
    rng = np.random.RandomState(0)
    groups = rng.randint(0, 60, 1000)
    tr, te = C.scaffold_split(groups, 0.2, seed=1, scheme=scheme)
    assert len(set(groups[tr]) & set(groups[te])) == 0
    assert len(tr) + len(te) == len(groups)
    assert len(te) > 0 and len(tr) > 0


def test_deterministic_split_is_seed_independent():
    groups = np.random.RandomState(0).randint(0, 40, 500)
    a = C.scaffold_split(groups, 0.2, seed=1, scheme="deterministic")
    b = C.scaffold_split(groups, 0.2, seed=99, scheme="deterministic")
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def test_random_split_does_not_systematically_exclude_large_groups():
    """Large scaffold families must be able to reach the test set.

    The previous implementation skipped any group that would overshoot the test
    quota, which pushed every large family into train and quietly enriched the test
    set in rare scaffolds.
    """
    sizes = np.array([200] + [5] * 60)
    groups = np.repeat(np.arange(len(sizes)), sizes)
    big_in_test = 0
    for seed in range(30):
        _, te = C.scaffold_split(groups, 0.35, seed=seed, scheme="random")
        if 0 in set(groups[te]):
            big_in_test += 1
    assert big_in_test > 0


# ----------------------------------------------------------------------------
# Similarity search
# ----------------------------------------------------------------------------

def test_knn_tanimoto_self_exclusion_uses_explicit_indices():
    """Leave-one-out must exclude the query's own training entry, not position i."""
    smiles = ["c1ccccc1", "c1ccccc1C", "Clc1ccccc1Br", "OC(=O)c1ccncc1N",
              "CN1CCN(CC1)c1ccc(cc1)C(=O)O"]
    fps = C.morgan_bitvects(smiles)
    y = np.array([1.0, 2.0, 3.0, 9.0, 9.5])
    # Query the last two molecules; they sit at training positions 3 and 4.
    pred_with_self, maxsim_with_self, _ = C.knn_tanimoto(fps[3:], fps, y, k=2)
    pred_loo, maxsim_loo, _ = C.knn_tanimoto(fps[3:], fps, y, k=2,
                                             self_indices=[3, 4])
    assert maxsim_with_self[0] == pytest.approx(1.0)
    assert maxsim_loo[0] < 1.0
    assert not np.allclose(pred_with_self, pred_loo)


def test_knn_tanimoto_handles_k_larger_than_training_set():
    smiles = ["c1ccccc1", "CCO"]
    fps = C.morgan_bitvects(smiles)
    pred, mx, _ = C.knn_tanimoto(fps, fps, np.array([1.0, 2.0]), k=10)
    assert np.isfinite(pred).all() and np.isfinite(mx).all()


# ----------------------------------------------------------------------------
# Conformal prediction
# ----------------------------------------------------------------------------

def _conformal_module():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "scripts", "04_conformal_ad_models.py")
    spec = importlib.util.spec_from_file_location("conformal04", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_conformal_coverage_is_near_nominal_under_exchangeability():
    """With exchangeable calibration and test data, coverage must hit 1 - alpha.

    This is the property the paper leans on, so it is tested directly rather than
    assumed: draw residuals from one distribution, calibrate on half, measure
    coverage on the other half.
    """
    m = _conformal_module()
    rng = np.random.RandomState(0)
    n = 4000
    residuals = np.abs(rng.normal(0, 1.0, n))
    sims = rng.uniform(0.0, 1.0, n)
    cal, test = slice(0, n // 2), slice(n // 2, n)
    edges = m.merge_bins(sims[cal], m.EDGES, m.MIN_CAL_PER_BIN)
    q = m.conformal_quantiles(residuals[cal], sims[cal], edges, alpha=0.10)
    b = m.bin_of(sims[test], edges)
    covered = residuals[test] <= np.array(q)[b]
    assert 0.86 <= covered.mean() <= 0.94


def test_conformal_level_formula_is_finite_sample_correct():
    m = _conformal_module()
    residuals = np.linspace(0, 1, 200)
    sims = np.full(200, 0.8)
    edges = [0.0, 1.0001]
    q = m.conformal_quantiles(residuals, sims, edges, alpha=0.10)
    # ceil((200+1)*0.9)/200 = 181/200 = 0.905
    assert q[0] == pytest.approx(np.quantile(residuals, 0.905), rel=1e-6)


def test_isotonic_tau_ignores_a_single_noisy_bin():
    """A non-monotone blip must not discard the usable bins beneath it.

    Half-widths [1.5, 1.75, 1.0, 0.9] across rising similarity. Scanning from the
    top, the old rule stopped at the 1.75 bin and returned 0.5, discarding the
    lowest bin even though the underlying trend is comfortably inside the budget.
    The isotonic fit pools the non-monotone pair to 1.625, which passes, so the
    whole similarity range is admitted.
    """
    m = _conformal_module()
    edges = [0.0, 0.4, 0.5, 0.6, 1.0001]
    tau, fitted = m.isotonic_tau(edges, [1.5, 1.75, 1.0, 0.9], max_half_width=1.65)
    assert tau is not None
    assert tau == pytest.approx(0.0)
    assert all(fitted[i] >= fitted[i + 1] - 1e-9 for i in range(len(fitted) - 1))


def test_isotonic_tau_returns_none_when_nothing_meets_the_budget():
    m = _conformal_module()
    tau, _ = m.isotonic_tau([0.0, 0.5, 1.0001], [4.0, 3.5], max_half_width=1.65)
    assert tau is None


def test_merge_bins_respects_minimum_count():
    m = _conformal_module()
    sims = np.concatenate([np.full(5, 0.1), np.full(500, 0.8)])
    edges = m.merge_bins(sims, [0.0, 0.4, 0.5, 0.6, 0.7, 1.0001], min_n=150)
    for lo, hi in zip(edges[:-1], edges[1:]):
        n = np.sum((sims >= lo) & (sims < hi))
        assert n >= 150 or len(edges) == 2


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def test_ranking_metrics_report_prevalence():
    """Prevalence must travel with AUC: 80%-active data makes AUC easy to misread."""
    y = np.concatenate([np.full(80, 7.0), np.full(20, 4.0)])
    p = y + np.random.RandomState(0).normal(0, 0.1, len(y))
    out = C.ranking_metrics(y, p)
    assert out["prevalence"] == pytest.approx(0.8)
    assert out["n"] == 100


def test_ranking_metrics_single_class_does_not_crash():
    y = np.full(50, 7.0)
    out = C.ranking_metrics(y, y)
    assert np.isnan(out["AUC"])


def test_enrichment_factor_and_bedroc_reward_early_recognition():
    n = 1000
    y = np.zeros(n, int)
    y[:50] = 1
    perfect = -np.arange(n, dtype=float)          # actives ranked first
    random_score = np.random.RandomState(0).normal(size=n)
    assert C.enrichment_factor(y, perfect, 0.01) > C.enrichment_factor(y, random_score, 0.01)
    assert C.bedroc(y, perfect) > 0.9
    assert 0.0 <= C.bedroc(y, random_score) <= 0.5


def test_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.RandomState(0)
    y = rng.normal(0, 1, 300)
    p = y + rng.normal(0, 0.5, 300)
    ci = C.r2_ci(y, p, n_boot=300)
    from sklearn.metrics import r2_score
    point = r2_score(y, p)
    assert ci["lo"] < point < ci["hi"]


# ----------------------------------------------------------------------------
# Chemotypes and feature contract
# ----------------------------------------------------------------------------

def test_every_chemotype_smarts_parses():
    from rdkit import Chem
    for name, sma in C.CHEMOTYPE_SMARTS.items():
        assert Chem.MolFromSmarts(sma) is not None, f"{name} is not valid SMARTS"


def test_pan_kinase_alert_flags_a_quinone():
    flags = C.chemotype_flags(["O=C1C=CC(=O)C=C1", "CCO"])
    assert bool(flags.loc[0, "pan_kinase_alert"]) is True
    assert bool(flags.loc[1, "pan_kinase_alert"]) is False


def test_feature_column_contract_is_stable():
    """Step 05 asserts on this ordering; pin it so a silent reorder fails here."""
    assert C.FEATURE_COLS[:len(C.DESC_COLS)] == C.DESC_COLS
    assert C.FEATURE_COLS[len(C.DESC_COLS):] == [f"FP_{i}" for i in range(C.N_FP_BITS)]
    assert len(C.FEATURE_COLS) == len(C.DESC_COLS) + C.N_FP_BITS


def test_descriptors_match_declared_columns():
    from rdkit import Chem
    d = C.physchem_descriptors(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    assert set(d) == set(C.DESC_COLS)


def test_resolve_binary_is_not_required_to_exist():
    assert C.resolve_binary("definitely-not-installed-xyz", required=False) is None


# ----------------------------------------------------------------------------
# Report generation against incomplete or stale inputs
# ----------------------------------------------------------------------------

def _report_module():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "scripts", "07_generate_report.py")
    spec = importlib.util.spec_from_file_location("report07", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_report_survives_metrics_written_by_an_older_step_04(tmp_path, monkeypatch):
    """A stale metrics file must produce a "regenerate this" note, not a crash.

    CI caught this: `final_model_*.json` from before the calibration rewrite has no
    `conformal_validity` key, and the report generator raised KeyError. The report
    has to describe whatever state the repository is actually in.
    """
    import json
    mod = _report_module()
    metrics = tmp_path / "metrics"
    metrics.mkdir()
    (metrics / "final_model_FGFR4.json").write_text(json.dumps({
        "target": "FGFR4", "tau": 0.5, "overall_coverage": 0.85,
        "bin_table": [], "alpha": 0.1,            # no conformal_validity
    }))
    monkeypatch.setattr(C, "METRICS_DIR", str(metrics))
    monkeypatch.setattr(mod.C, "METRICS_DIR", str(metrics))
    section = mod.domain_section()
    assert "04_conformal_ad_models.py" in section
    assert "FGFR4" in section


def test_report_sections_report_missing_inputs_rather_than_failing(tmp_path, monkeypatch):
    mod = _report_module()
    empty = tmp_path / "metrics"
    empty.mkdir()
    monkeypatch.setattr(C, "METRICS_DIR", str(empty))
    monkeypatch.setattr(mod.C, "METRICS_DIR", str(empty))
    monkeypatch.chdir(tmp_path)
    for section in (mod.curation_section(), mod.benchmark_section(),
                    mod.domain_section(), mod.screen_section(),
                    mod.docking_section()):
        assert "Not run" in section or "Regenerate" in section or "regenerate" in section
    assert isinstance(mod.limitations_section(), str)


# ----------------------------------------------------------------------------
# Redocking RMSD: the control must measure placement, not conformation
# ----------------------------------------------------------------------------

def _docking_module():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "scripts", "06_docking_validation.py")
    spec = importlib.util.spec_from_file_location("docking06", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _embedded(smiles, seed=0xf00d):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    assert AllChem.EmbedMolecule(mol, params) == 0
    return Chem.RemoveHs(mol)


def test_redock_rmsd_is_measured_in_place_not_superimposed():
    """A pose translated 8 A away is a failed prediction, however well it
    superimposes. `GetBestRMS` would align it first and call it perfect, which is
    how a pose in the wrong subpocket passes a 2 A criterion."""
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms  # noqa: F401  (import guard)
    import numpy as np

    mod = _docking_module()
    ref = _embedded("CC(=O)Oc1ccccc1C(=O)O")

    moved = Chem.Mol(ref)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 8.0, p.y, p.z))

    in_place = mod.symmetry_corrected_rmsd(ref, moved)
    superimposed = mod.symmetry_corrected_rmsd(ref, moved, superimpose=True)

    assert in_place == pytest.approx(8.0, abs=1e-3), in_place
    assert superimposed < 0.1, superimposed


def test_rmsd_matches_molecules_whose_bond_orders_were_perceived_differently():
    """The crystal reference and the Open Babel pose kekulise differently. Matching
    on the heavy-atom graph is what makes the RMSD computable at all; before this
    the control returned None for every target and reported it as a failed pose."""
    from rdkit import Chem

    mod = _docking_module()
    ref = _embedded("c1ccccc1C(=O)[O-]")

    scrambled = Chem.Mol(ref)
    for b in scrambled.GetBonds():
        b.SetIsAromatic(False)
    for a in scrambled.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)

    rms = mod.symmetry_corrected_rmsd(ref, scrambled)
    assert rms is not None
    assert rms == pytest.approx(0.0, abs=1e-6)


def test_rmsd_returns_none_for_a_different_molecule():
    mod = _docking_module()
    assert mod.symmetry_corrected_rmsd(_embedded("CCO"), _embedded("CCCCCCN")) is None


def test_candidate_docking_names_step_05_when_there_is_no_screen(tmp_path, monkeypatch):
    """With no screen at all the blocker is step 05, not step 07's annotation."""
    mod = _docking_module()
    monkeypatch.chdir(tmp_path)
    with pytest.raises(mod.StaleScreenOutput) as excinfo:
        mod.dock_consensus_candidates()
    assert "05_screen_coconut_library" in str(excinfo.value)


def test_candidate_docking_refuses_a_pre_rewrite_screen_output(tmp_path, monkeypatch):
    """Step 06 reuses whatever screen output it finds. An old one names its hits
    by a rule this work withdrew, so docking it would attach current scores to a
    candidate list the current screen would not produce."""
    import pandas as pd

    mod = _docking_module()
    preds = tmp_path / "results" / "predictions"
    preds.mkdir(parents=True)
    stale = pd.DataFrame({"identifier": ["CNP1"], "smiles": ["CCO"],
                          "n_targets_confident": [3], "rediscovery_any": [False]})
    stale.to_csv(preds / "consensus_hits_annotated.csv", index=False)
    stale.to_csv(preds / "consensus_hits.csv", index=False)
    monkeypatch.chdir(tmp_path)

    with pytest.raises(mod.StaleScreenOutput) as excinfo:
        mod.dock_consensus_candidates()
    assert "n_targets_interval_active" in str(excinfo.value)
    assert "05_screen_coconut_library" in str(excinfo.value)


def test_protocol_note_is_silent_at_the_defaults_and_names_every_override():
    mod = _docking_module()
    assert mod.protocol_note() is None

    import importlib
    os.environ["DOCK_ENRICH_N_SEEDS"] = "1"
    try:
        reduced = _docking_module()
        note = reduced.protocol_note()
        assert note is not None
        assert "enrich_n_seeds 1 instead of 3" in note
        assert "DOCK_ENRICH_N_SEEDS" in note
    finally:
        del os.environ["DOCK_ENRICH_N_SEEDS"]


def test_bedroc_on_a_balanced_set_sits_inside_its_own_null():
    """EF and BEDROC assume a large library with few actives. On a 50/50 control
    a random ranking returns BEDROC near 0.5 with a wide spread, so the metric
    must be read against its null rather than its usual thresholds."""
    import numpy as np

    labels = np.array([1] * 50 + [0] * 49)
    rng = np.random.default_rng(7)
    null = C.permutation_null(labels, rng.normal(size=99),
                              lambda y, s: C.bedroc(y, s, 20.0), n_perm=400)
    assert null is not None
    assert 0.4 < null["null_mean"] < 0.6, null["null_mean"]
    assert null["null_95pct"] > 0.7, null["null_95pct"]
    assert not null["exceeds_null"]


def test_permutation_null_flags_a_ranking_that_really_separates():
    import numpy as np
    from sklearn.metrics import roc_auc_score

    labels = np.array([1] * 50 + [0] * 50)
    perfect = np.concatenate([np.ones(50), np.zeros(50)])
    null = C.permutation_null(labels, perfect, roc_auc_score, n_perm=400)
    assert null["observed"] == pytest.approx(1.0)
    assert null["exceeds_null"]
    assert null["p_one_sided"] < 0.01


# ----------------------------------------------------------------------------
# Benchmark checkpoint
# ----------------------------------------------------------------------------

def _benchmark_module():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "scripts", "03_scaffold_benchmark.py")
    spec = importlib.util.spec_from_file_location("benchmark03", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_benchmark_checkpoint_discards_rows_from_changed_data(tmp_path, monkeypatch):
    """A seed is only "done" for the data it was computed from.

    Keying the checkpoint on (target, seed) alone meant that re-running after
    step 01 or 02 changed the data skipped every seed as already done, and the
    summary then averaged rows from two datasets into one model selection.
    """
    import json as _json

    mod = _benchmark_module()
    jsonl = tmp_path / "benchmark_long.jsonl"
    monkeypatch.setattr(mod, "JSONL", str(jsonl))

    rows = [{"target": "EGFR", "seed": 0, "model": "RF-ecfp", "data_sha": "old"},
            {"target": "EGFR", "seed": 1, "model": "RF-ecfp", "data_sha": "new"},
            {"target": "FGFR4", "seed": 0, "model": "RF-ecfp"}]  # predates the field
    jsonl.write_text("".join(_json.dumps(r) + "\n" for r in rows))

    assert mod.prune_stale({"EGFR": "new", "FGFR4": "new"}) == 2
    kept = [_json.loads(l) for l in jsonl.read_text().splitlines()]
    assert kept == [rows[1]]
    assert mod.done_keys() == {("EGFR", 1, "new")}


def test_parallel_standardisation_matches_the_serial_path(monkeypatch):
    """The pool is an optimisation, not a change in what is computed.

    Serially the 738k screening library takes ~21 h in step 05, which is why this
    is parallelised; the parallel path is only admissible if it returns exactly
    what the serial one does, in the same order.
    """
    smiles = ["CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
              "not_a_molecule", "c1ccccc1", "OC(=O)CC(O)(CC(=O)O)C(=O)O",
              "CC(C)Cc1ccc(cc1)C(C)C(=O)O"] * 8
    monkeypatch.setattr(C, "_PARALLEL_STANDARDISE_MIN", 4)

    serial = C.standardise_series(smiles, n_jobs=1)
    parallel = C.standardise_series(smiles, n_jobs=2)
    assert parallel == serial
    assert len(serial[0]) == len(smiles)
    assert serial[0][2] is None and serial[1][2] is None  # unparsable stays None


def test_candidate_staleness_guard_accepts_the_columns_step_05_writes(tmp_path,
                                                                      monkeypatch):
    """The guard must reject a pre-rewrite screen and pass a current one.

    It named `interval_lower_bound_active` and `known_scaffold` -- keys of
    screen_summary.json and of a frame internal to step 05, never columns of
    consensus_hits_annotated.csv. Every current screen was therefore rejected and
    candidate docking never ran.
    """
    import pandas as pd

    mod = _docking_module()
    monkeypatch.chdir(tmp_path)
    (tmp_path / "results" / "predictions").mkdir(parents=True)

    current = ["identifier", "name", "smiles", "pan_kinase_alert",
               "n_targets_interval_active", "n_targets_inAD_active",
               "known_scaffold_all_targets", "rediscovery_any", "mean_pred"]
    pre_rewrite = ["identifier", "name", "smiles", "pan_kinase_alert",
                   "confident_active", "mean_pred"]

    def write(cols):
        for name in ("consensus_hits.csv", "consensus_hits_annotated.csv"):
            pd.DataFrame(columns=cols).to_csv(
                tmp_path / "results" / "predictions" / name, index=False)

    write(pre_rewrite)
    with pytest.raises(mod.StaleScreenOutput, match="predates the current step 05"):
        mod.dock_consensus_candidates()

    # The current header clears the staleness check; it then fails later, on the
    # empty candidate frame, which is proof the guard is no longer the blocker.
    write(current)
    with pytest.raises(Exception) as exc:
        mod.dock_consensus_candidates()
    assert not isinstance(exc.value, mod.StaleScreenOutput)


# ----------------------------------------------------------------------------
# Screen shortlist rule
# ----------------------------------------------------------------------------

def _screen_module():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "scripts", "05_screen_coconut_library.py")
    spec = importlib.util.spec_from_file_location("screen05", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_shortlist_is_the_interval_criterion_not_the_domain_flag():
    """Being in-domain and predicted active must not put a compound on the list.

    tau is non-binding on two of three targets, so `in_AD` admits the whole
    library there; an `n_targets_inAD_active >= 1` arm selected 94% of the parents
    and the shortlist stopped naming anything. Ordering is unaffected by dropping
    it, which is why the annotated and docked subsets did not move.
    """
    import pandas as pd

    mod = _screen_module()
    cons = pd.DataFrame({
        "identifier": ["keep2", "keep1", "drop", "drop_too"],
        "n_targets_interval_active": [2, 1, 0, 0],
        "n_targets_inAD_active": [3, 1, 3, 0],
        "mean_pred": [8.0, 7.0, 9.9, 6.0],
    })

    hits = mod.shortlist(cons)
    assert list(hits["identifier"]) == ["keep2", "keep1"]
    # the dropped rows outrank nothing that was kept, whatever their in-AD count
    assert hits["n_targets_interval_active"].is_monotonic_decreasing
