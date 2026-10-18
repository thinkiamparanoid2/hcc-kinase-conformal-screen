"""
Step 07 - Results report and figures
-------------------------------------
Every number in the report is read from results/metrics/*.json|csv and
results/predictions/*. Nothing is hand-typed -- including the narrative: chemotype
groupings come from the SMARTS definitions in `common.CHEMOTYPE_SMARTS`, not from
a prose paragraph with numbers pasted into it, and not from thresholding a
continuous NP-likeness score.

Sections whose inputs are missing are written as "not run", with the command that
produces them, instead of being silently skipped.

Outputs: results/reports/results_report.md
         results/predictions/consensus_hits_annotated.csv
         results/figures/fig{1,2,3}_*.png
"""

import os
import sys
import json

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C

REPORT = os.path.join("results", "reports", "results_report.md")
FIG_DIR = os.path.join("results", "figures")


def md_table(df, floatfmt="{:.3f}"):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |",
           "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for v in r:
            if isinstance(v, (float, np.floating)):
                cells.append(floatfmt.format(v) if np.isfinite(v) else "-")
            elif v is None:
                cells.append("-")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def ci_str(ci, fmt="{:.3f}"):
    if not ci or ci.get("lo") is None:
        return "-"
    return f"[{fmt.format(ci['lo'])}, {fmt.format(ci['hi'])}]"


def missing(section, command):
    return (f"## {section}\n\n_Not run. Produce it with:_ `{command}`\n")


# ----------------------------------------------------------------------------
# 0. Data curation
# ----------------------------------------------------------------------------

def curation_section():
    rows, notes = [], []
    for t in C.TARGETS:
        cur = C.load_json(os.path.join(C.METRICS_DIR, f"curation_{t}.json"))
        if not cur:
            continue
        rows.append({
            "target": t,
            "final compounds": cur.get("n_final"),
            "exact labels": cur.get("n_exact"),
            "censored (potent bound)": cur.get("n_censored_lower_bound_potency"),
            "censored (inactive bound)": cur.get("n_censored_upper_bound_potency"),
            "active prevalence (exact)": cur.get("exact_active_prevalence_at_pIC50_6"),
            "R2 ceiling from label noise": cur.get("r2_ceiling_from_label_noise"),
        })
        if cur.get("label_noise_note"):
            notes.append(f"- {t}: {cur['label_noise_note']}")
        if cur.get("censored_records_available") is False:
            notes.append(
                f"- {t}: the curated file contains no `>` records, so there is NO "
                f"inactive class. Every classification metric below therefore "
                f"describes ranking within a set that is "
                f"{cur.get('exact_active_prevalence_at_pIC50_6', float('nan')):.0%} "
                f"active, not active-versus-inactive discrimination. Re-run step 01 "
                f"with network access to populate the inactive class.")
    if not rows:
        return missing("0. Data curation", "python scripts/01_fetch_chembl_data.py")
    lines = ["## 0. Data curation\n",
             "Structures are salt-stripped, neutralised and tautomer-canonicalised, "
             "then deduplicated on the parent InChIKey. Censored measurements are "
             "labelled, not silently treated as point values: only uncensored rows "
             "are used as regression targets.\n",
             md_table(pd.DataFrame(rows)) + "\n"]
    if notes:
        lines.append("**Curation caveats**\n")
        lines.extend(notes)
        lines.append("")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 1. Benchmark
# ----------------------------------------------------------------------------

def benchmark_section():
    path = os.path.join(C.METRICS_DIR, "benchmark_long.csv")
    if not os.path.exists(path):
        return missing("1. Benchmark", "python scripts/03_scaffold_benchmark.py")
    df = pd.read_csv(path)
    canonical = C.load_json(os.path.join(C.METRICS_DIR, "benchmark_canonical.json"), {})
    temporal = C.load_json(os.path.join(C.METRICS_DIR, "benchmark_temporal.json"), {})
    selection = C.load_json(os.path.join(C.METRICS_DIR, "model_selection.json"), {})

    lines = ["## 1. Benchmark: scaffold-split evaluation\n"]
    lines.append(
        "Scaffold splits are generic Bemis-Murcko and strictly disjoint. "
        "The canonical split is deterministic (largest scaffold families in train) "
        "and carries percentile bootstrap CIs over test molecules. The repeated "
        "splits measure how much the result moves with the scaffold assignment: "
        "that spread is SPLIT variance, not sampling variance, because the same "
        "molecules are reused each time.\n")
    if selection.get("hyperparameters"):
        lines.append(f"_{selection['hyperparameters']}_\n")

    lines.append("### 1.1 Canonical deterministic scaffold split (with 95% CIs)\n")
    for t in C.TARGETS:
        block = (canonical.get("splits") or {}).get(t)
        if not block:
            continue
        rows = []
        for name, m in block["models"].items():
            rows.append({"model": name,
                         "R2": m["R2"], "R2 95% CI": ci_str(m.get("R2_CI95")),
                         "RMSE": m["RMSE"], "Spearman": m["Spearman"],
                         "AUC": m.get("AUC"), "AUC 95% CI": ci_str(m.get("AUC_CI95"))})
        lines.append(f"**{t}** (train {block['n_train']}, test {block['n_test']}, "
                     f"test active prevalence {block['test_active_prevalence']:.2f}, "
                     f"scaffold overlap {block['scaffold_overlap']})\n")
        lines.append(md_table(pd.DataFrame(rows)) + "\n")

    lines.append("### 1.2 Repeated random scaffold splits (mean +/- SD over seeds)\n")
    order = ["mean", "kNN-Tanimoto", "RF-desc", "RF-ecfp", "RF-all", "XGB-all",
             "XGB-scramble"]
    for t in C.TARGETS:
        d = df[df.target == t]
        if d.empty:
            continue
        rows = []
        for m in order:
            x = d[d.model == m]
            if x.empty:
                continue
            rows.append({"model": m,
                         "R2": f"{x.R2.mean():.3f} +/- {x.R2.std():.3f}",
                         "RMSE": f"{x.RMSE.mean():.3f} +/- {x.RMSE.std():.3f}",
                         "Spearman": f"{x.Spearman.mean():.3f} +/- {x.Spearman.std():.3f}",
                         "AUC": f"{x.AUC.mean():.3f} +/- {x.AUC.std():.3f}"})
        lines.append(f"**{t}** ({d.seed.nunique()} splits, mean test active "
                     f"prevalence {d.prevalence.mean():.2f})\n")
        lines.append(md_table(pd.DataFrame(rows)) + "\n")

    paired = C.load_json(os.path.join(C.METRICS_DIR, "benchmark_paired.json"), [])
    if paired:
        lines.append("### 1.3 Paired comparisons across splits\n")
        pr = pd.DataFrame(paired)
        keep = [c for c in ["target", "comparison", "mean_dR2", "wins", "n_seeds",
                            "wilcoxon_p"] if c in pr.columns]
        lines.append(md_table(pr[keep], "{:.4f}") + "\n")
        lines.append(f"_{paired[0].get('interpretation', '')}_\n")

    tsplits = (temporal.get("splits") or {})
    ok = {t: b for t, b in tsplits.items() if b.get("status") == "ok"}
    lines.append("### 1.4 Temporal holdout (train on older literature)\n")
    if ok:
        rows = []
        for t, b in ok.items():
            for name, m in b["models"].items():
                rows.append({"target": t, "model": name, "cutoff year": b["cutoff_year"],
                             "n test": b["n_test"], "R2": m["R2"],
                             "R2 95% CI": ci_str(m.get("R2_CI95")),
                             "AUC": m.get("AUC")})
        lines.append(md_table(pd.DataFrame(rows)) + "\n")
    else:
        reasons = {t: b.get("reason") for t, b in tsplits.items()}
        lines.append("_Not available._ " + "; ".join(
            f"{t}: {r}" for t, r in reasons.items() if r) + "\n")

    if selection:
        lines.append(
            f"### 1.5 Model selection\n\nSelected family: **{selection['selected_model']}** "
            f"({selection.get('criterion', '')}). Mean R2 by family: "
            f"`{json.dumps({k: round(v, 3) for k, v in selection['mean_R2_by_model'].items()})}`. "
            f"{selection.get('note', '')}\n")

    lines.append(
        "Controls: `mean` and `XGB-scramble` (shuffled labels) must sit at R2 <= 0 "
        "and AUC ~ 0.5. If they do not, the pipeline leaks.\n")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 2. Conformal calibration and applicability domain
# ----------------------------------------------------------------------------

def domain_section():
    blocks = {t: C.load_json(os.path.join(C.METRICS_DIR, f"final_model_{t}.json"))
              for t in C.TARGETS}
    blocks = {t: b for t, b in blocks.items() if b}
    if not blocks:
        return missing("2. Applicability domain and conformal calibration",
                       "python scripts/04_conformal_ad_models.py")

    # A metrics file written by an older step 04 lacks the keys below. Report that
    # and name the command, rather than crashing on a stale artifact: the report
    # must describe whatever state the repository is actually in.
    stale = [t for t, b in blocks.items() if "conformal_validity" not in b]
    blocks = {t: b for t, b in blocks.items() if "conformal_validity" in b}
    if not blocks:
        return ("## 2. Applicability domain, conformal calibration and NP transfer\n\n"
                f"_Metrics for {', '.join(stale)} predate the current calibration "
                "and carry no validity statement. Regenerate with:_ "
                "`python scripts/04_conformal_ad_models.py`\n")

    first = next(iter(blocks.values()))
    validity = first.get("conformal_validity", {})
    lines = [
        "## 2. Applicability domain, conformal calibration and NP transfer\n",
        "Mondrian split-conformal intervals (90%), with the taxonomy defined by "
        "max-Tanimoto to the training split. `tau` is read off an isotonic "
        "(non-increasing) fit of half-width against similarity, so a single noisy "
        "bin cannot discard usable lower bins, and is reported with a bootstrap CI "
        "and a sensitivity analysis.\n",
        "> **Validity statement.** " +
        f"The {validity.get('nominal_coverage', 0.9):.0%} guarantee holds for "
        f"{validity.get('holds_for', 'held-out scaffolds')}. It does **not** hold "
        f"for {validity.get('does_not_hold_for', 'shifted deployment data')}\n",
    ]
    if stale:
        lines.append(f"> Metrics for {', '.join(stale)} predate the current "
                     f"calibration and are omitted here. Regenerate with "
                     f"`python scripts/04_conformal_ad_models.py`.\n")
    for t in C.TARGETS:
        r = blocks.get(t)
        if not r:
            continue
        cov_ci = ci_str(r.get("conformal_validity", {}).get("empirical_coverage_CI95"))
        tau_ci = r.get("tau_CI95", {})
        lines.append(f"### {t} ({r['model']}; train {r['n_train']} / cal {r['n_cal']} "
                     f"/ test {r['n_test']})\n")
        lines.append(
            f"- empirical coverage **{r['overall_coverage']:.3f}** {cov_ci} "
            f"against a nominal {1 - r['alpha']:.2f}; mean half-width "
            f"{r['mean_half_width']:.2f} pIC50 units\n"
            f"- tau = **{r['tau']}** (95% CI "
            f"{tau_ci.get('lo')}-{tau_ci.get('hi')}); held-out molecules inside the "
            f"domain: {r['frac_test_in_AD'] * 100:.1f}%\n")
        if r.get("tau_interpretation"):
            lines.append(f"- {r['tau_interpretation']}\n")
        if r.get("tau_by_budget"):
            lines.append("**tau at other width budgets**\n")
            lines.append(md_table(pd.DataFrame(r["tau_by_budget"])) + "\n")
        bt = pd.DataFrame(r["bin_table"])
        bt["bin"] = bt.apply(
            lambda x: f"[{x.sim_lo:.2f}, {min(x.sim_hi, 1.0):.2f})", axis=1)
        bt["coverage 95% CI"] = [ci_str(c) for c in
                                 bt.get("test_coverage_CI95", [None] * len(bt))]
        show = ["bin", "n_cal", "n_test", "half_width_90",
                "half_width_90_isotonic", "test_MAE", "test_coverage",
                "coverage 95% CI"]
        lines.append(md_table(bt[[c for c in show if c in bt.columns]]) + "\n")

        if r.get("tau_sensitivity"):
            lines.append("**tau sensitivity**\n")
            lines.append(md_table(pd.DataFrame(r["tau_sensitivity"])) + "\n")

        sub = []
        for key, label in [("test_all", "all held-out"),
                           ("test_in_AD", "inside domain"),
                           ("test_out_AD", "outside domain"),
                           ("test_NP_like", "NP-like (Ertl score >= 0)"),
                           ("test_not_NP_like", "not NP-like")]:
            v = r.get(key, {})
            sub.append({"subset": label, "n": v.get("n"),
                        "R2": v.get("R2", np.nan),
                        "R2 95% CI": ci_str(v.get("R2_CI95")),
                        "RMSE": v.get("RMSE", np.nan),
                        "note": v.get("note", "")})
        lines.append(md_table(pd.DataFrame(sub)) + "\n")

        npt = r.get("np_transfer", {})
        lines.append("**Natural-product transfer**\n")
        lines.append(f"- {npt.get('status', '')}\n")
        if npt.get("coverage_on_NP_like_test") is not None:
            lines.append(
                f"- coverage on NP-like held-out molecules: "
                f"{npt['coverage_on_NP_like_test']:.3f} "
                f"{ci_str(npt.get('coverage_on_NP_like_test_CI95'))} "
                f"(n = {npt['n_test_NP_like']}) versus "
                f"{npt['coverage_on_non_NP_like_test']:.3f} "
                f"{ci_str(npt.get('coverage_on_non_NP_like_test_CI95'))} on the rest\n")
        if npt.get("np_specific_coverage_on_NP_like_test") is not None:
            lines.append(
                f"- after NP-specific recalibration: "
                f"{npt['np_specific_coverage_on_NP_like_test']:.3f} "
                f"{ci_str(npt.get('np_specific_coverage_CI95'))} "
                f"(half-width {npt['np_specific_half_width_90']:.2f})\n")
        lines.append("")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 3. Screen
# ----------------------------------------------------------------------------

def annotate(cons, n=200):
    from rdkit import Chem
    from rdkit.Chem import QED
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog(params)
    top = cons.head(n).copy()
    flags, qeds = [], []
    for s in top["smiles"]:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        flags.append(bool(pains.HasMatch(m)) if m else True)
        qeds.append(QED.qed(m) if m else np.nan)
    top["PAINS"] = flags
    top["QED"] = qeds
    chemo = C.chemotype_flags(top["smiles"])
    for col in chemo.columns:
        if col not in top.columns:
            top[col] = chemo[col].to_numpy()
    return top


def screen_section():
    path = os.path.join(C.METRICS_DIR, "screen_summary.json")
    cons_path = os.path.join("results", "predictions", "consensus_hits.csv")
    if not (os.path.exists(path) and os.path.exists(cons_path)):
        return missing("3. Full-library screen",
                       "python scripts/05_screen_coconut_library.py")
    s = C.load_json(path)
    # A summary written by the pre-rewrite step 05 has different keys, different
    # taus and a different model, and its counts were computed without
    # stereoisomer deduplication or scaffold-level novelty. Presenting it beside
    # the current calibration would state numbers the current pipeline never
    # produced, so it is refused rather than rendered.
    if "library_records" not in s or "interval_caveat" not in s:
        return ("## 3. Full-library natural-product screen\n\n"
                "_The screen outputs in this repository predate the current "
                "pipeline: they were produced with a different model family, "
                "different applicability-domain thresholds, and without "
                "stereoisomer deduplication or scaffold-level novelty. They are "
                "withheld here rather than reported beside the current "
                "calibration. Regenerate with:_ "
                "`python scripts/05_screen_coconut_library.py` "
                "_(needs the COCONUT library; see `data/README.md`)._\n")
    lines = ["## 3. Full-library natural-product screen\n",
             f"Library: {s.get('library_records', 0):,} records collapsed to "
             f"{s.get('library_unique_parents', 0):,} unique 2D parents "
             f"({s.get('library_featurised', 0):,} featurised). "
             f"{s.get('dedup_note', '')}\n",
             f"> **Interval caveat.** {s.get('interval_caveat', '')}\n"]
    rows = [{"target": t, **{k: v for k, v in d.items()
                             if not isinstance(v, (dict, list))}}
            for t, d in s.get("targets", {}).items()]
    if rows:
        lines.append(md_table(pd.DataFrame(rows)) + "\n")
    cons_summary = s.get("consensus", {})
    if cons_summary:
        lines.append("### 3.1 Consensus across targets\n")
        flat = {k: v for k, v in cons_summary.items() if not isinstance(v, (dict, str))}
        lines.append(md_table(pd.DataFrame([flat])) + "\n")
        cross = cons_summary.get("cross_target_prediction_spearman", {})
        if cross:
            lines.append("Cross-target prediction rank correlation: "
                         f"`{json.dumps({k: round(v, 3) for k, v in cross.items()})}`\n")
        if cons_summary.get("consensus_caveat"):
            lines.append(f"> {cons_summary['consensus_caveat']}\n")

    cons = pd.read_csv(cons_path)
    novelty_cols = [c for c in ["rediscovery_any", "known_scaffold_all_targets"]
                    if c in cons.columns]
    novel = cons.copy()
    for c in novelty_cols:
        novel = novel[~novel[c].astype(bool)]
    ann = annotate(novel)
    ann.to_csv(os.path.join("results", "predictions",
                            "consensus_hits_annotated.csv"), index=False)

    cols = [c for c in ["identifier", "name", "n_targets_interval_active",
                        "n_targets_inAD_active", "mean_pred", "np_likeness",
                        "QED", "PAINS", "pan_kinase_alert"] if c in ann.columns]
    lines.append("### 3.2 Top candidates on scaffolds absent from the training sets\n")
    if len(ann):
        lines.append(md_table(ann[cols].head(25)) + "\n")
    else:
        lines.append("_No hit survives both novelty filters (not a near-duplicate of "
                     "training chemistry, and not on a scaffold already present in "
                     "the training set). That is itself the result: the "
                     "AD-admissible fraction of this library is known kinase "
                     "chemistry._\n")

    known = cons[cons.get("known_scaffold_all_targets", False).astype(bool)] \
        if "known_scaffold_all_targets" in cons.columns else cons.iloc[:0]
    lines.append("### 3.3 Hits on scaffolds already present in training "
                 "(analogues, positive controls)\n")
    kcols = [c for c in ["identifier", "name", "n_targets_inAD_active", "mean_pred"]
             if c in known.columns]
    lines.append(md_table(known[kcols].head(15)) + "\n" if len(known)
                 else "_none_\n")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 4. Docking
# ----------------------------------------------------------------------------

def docking_section():
    dock_json = os.path.join(C.METRICS_DIR, "docking_validation.json")
    dock_csv = os.path.join("results", "predictions", "consensus_hits_docked.csv")
    if not os.path.exists(dock_json):
        return missing("4. Structure-based docking validation",
                       "python scripts/06_docking_validation.py")
    dv = C.load_json(dock_json)
    # Same reasoning as the screen section: a pre-rewrite docking run reports
    # best-of-N RMSD as its control, 12-vs-12 enrichment, and cross-receptor
    # kcal/mol comparisons. Those are the numbers this work set out to stop
    # quoting, so they are not rendered as current results.
    if "settings" not in dv or "score_comparability_note" not in dv:
        return ("## 4. Orthogonal structure-based validation (AutoDock Vina)\n\n"
                "_The docking outputs in this repository predate the current "
                "validation protocol: they report best-of-N rather than top-1 "
                "RMSD without symmetry correction, a 12-vs-12 enrichment control "
                "drawn from potency extremes, and cross-receptor score "
                "comparisons. They are withheld here. Regenerate with:_ "
                "`python scripts/06_docking_validation.py` "
                "_(needs `vina` and `obabel`; see `docs/REPRODUCE.md`)._\n")
    settings = dv.get("settings", {})

    lines = ["## 4. Orthogonal structure-based validation (AutoDock Vina)\n"]
    if settings.get("protocol_note"):
        lines.append(f"> **{settings['protocol_note']}**\n")
    versions = settings.get("tool_versions") or {}
    if versions:
        lines.append("_Produced by " + ", ".join(
            f"{k.replace('_', ' ')} {v}" for k, v in sorted(versions.items())) +
            f", {settings.get('vina_cpu_threads')} threads per docking._\n")
    lines += [
        f"Boxes are centred on the co-crystal ligand heavy-atom centroid in chain A "
        f"and sized to the ligand (minimum {settings.get('min_box_A')} A, "
        f"{settings.get('box_padding_A')} A padding), pH 7.4 protonation, Gasteiger "
        f"charges, exhaustiveness {settings.get('exhaustiveness')}, "
        f"{settings.get('n_seeds')} seeds per docking.\n",
        f"> **Score comparability.** {dv.get('score_comparability_note', '')}\n",
        "### 4.1 Control 1: co-crystal redocking (pose prediction)\n",
        "Top-1 RMSD is the control: it is what the method would produce "
        "prospectively. Best-of-N is a sampling diagnostic and is reported only "
        "beside it.\n",
    ]
    rows = []
    for t in C.TARGETS:
        rc = (dv.get("redocking_control") or {}).get(t)
        if not rc:
            continue
        rows.append({
            "target": t, "ligand": rc["ligand"],
            "box (A)": rc["box_size"][0],
            "Vina top-1 (kcal/mol)": rc.get("affinity_top1_kcal_mol"),
            "top-1 RMSD (A)": rc.get("rmsd_top1_A"),
            "pass at 2 A": rc.get("pose_prediction_success_top1"),
            f"best-of-N RMSD (A)": rc.get("rmsd_best_of_n_A"),
            "N poses": rc.get("n_poses"),
            "native pose rank": rc.get("best_pose_rank"),
            "score gap (kcal/mol)": rc.get("affinity_gap_to_best_pose_kcal"),
        })
    if rows:
        lines.append(md_table(pd.DataFrame(rows)) + "\n")
        n_pass = sum(1 for r in rows if r["pass at 2 A"])
        lines.append(f"Pose prediction succeeds at the conventional 2 A criterion "
                     f"for {n_pass} of {len(rows)} pockets.\n")
        ranked = [r for r in rows if r.get("native pose rank")]
        missed = [r for r in ranked if r["native pose rank"] > 1]
        if missed:
            worst = max(abs(r["score gap (kcal/mol)"] or 0) for r in missed)
            lines.append(
                "Where top-1 fails, the search still produced a near-native pose: "
                f"it is ranked {', '.join(str(r['native pose rank']) for r in missed)} "
                f"in {', '.join(r['target'] for r in missed)} respectively, separated "
                f"from the top-scoring pose by at most {worst:.2f} kcal/mol. Vina's "
                "scoring function is not quoted to that precision, so these are "
                "ranking failures rather than sampling failures: the pose is found "
                "and cannot be told apart from its rivals. This is the reason "
                "docking enters this work as an orthogonal check and not as a "
                "source of truth about binding mode.\n")

    lines.append("### 4.2 Control 2: active-versus-inactive enrichment\n")
    erows = []
    for t in C.TARGETS:
        ec = (dv.get("enrichment_control") or {}).get(t)
        if not ec:
            continue
        if ec.get("status") != "ok":
            erows.append({"target": t, "status": ec.get("status"),
                          "reason": ec.get("reason")})
            continue
        erows.append({
            "target": t, "n act / n inact": f"{ec['n_actives']} / {ec['n_inactives']}",
            "active mean": ec["active_mean_kcal"],
            "inactive mean": ec["inactive_mean_kcal"],
            "ROC-AUC": ec["docking_ROC_AUC"],
            "AUC 95% CI": ci_str(ec.get("docking_ROC_AUC_CI95")),
            "EF1%": ec.get("EF1pct"), "EF5%": ec.get("EF5pct"),
            "BEDROC": ec.get("BEDROC_alpha20"),
        })
    if erows:
        lines.append(md_table(pd.DataFrame(erows), "{:.3f}") + "\n")
        sampling = next((ec.get("sampling") for ec in
                         (dv.get("enrichment_control") or {}).values()
                         if isinstance(ec, dict) and ec.get("sampling")), None)
        if sampling:
            lines.append(f"_Sampling: {sampling}._\n")
        inactive_class = next((ec.get("inactive_class") for ec in
                               (dv.get("enrichment_control") or {}).values()
                               if isinstance(ec, dict) and ec.get("inactive_class")), None)
        if inactive_class:
            lines.append(
                f"> **What the negative class is.** {inactive_class}. This matters "
                "for how the AUC reads: a true inactive is a `>` record, where the "
                "assay ran and found nothing, whereas a measured weak binder still "
                "occupies the pocket. Discriminating actives from weak binders is "
                "the harder task, so an AUC measured this way understates what the "
                "same control would give against true non-binders, and is not "
                "comparable to a published enrichment that used decoys.\n")

    cand_status = (dv.get("candidates_summary") or {}).get("status", "ok")
    if cand_status != "ok" or not os.path.exists(dock_csv):
        lines.append("### 4.3 Candidate docking\n")
        lines.append(
            "_Not available. The two controls above stand on their own: they need "
            "only the co-crystal structures and the curated training sets. Docking "
            "the candidate shortlist needs a current full-library screen, and there "
            "is none in this repository._\n")
        if cand_status != "ok":
            lines.append(f"_Reason: {cand_status}_\n")
        return "\n".join(lines)

    df = pd.read_csv(dock_csv)
    lines.append("### 4.3 Candidate docking by chemotype\n")
    chemo_cols = [c for c in df.columns if c.startswith("chemotype_")] or \
                 [c for c in C.CHEMOTYPE_SMARTS if c in df.columns]
    if not chemo_cols and "smiles" in df.columns:
        flags = C.chemotype_flags(df["smiles"])
        for col in C.CHEMOTYPE_SMARTS:
            if col in flags.columns:
                df[col] = flags[col].to_numpy()
        chemo_cols = [c for c in C.CHEMOTYPE_SMARTS if c in df.columns]

    grows = []
    for chem in chemo_cols:
        name = chem.replace("chemotype_", "")
        sub = df[df[chem].astype(bool)]
        if sub.empty:
            continue
        row = {"chemotype": name, "n": len(sub)}
        for t in C.TARGETS:
            col = f"{t}_vina_kcal"
            zcol = f"{t}_vina_z"
            if col in sub.columns:
                row[f"{t} mean"] = float(sub[col].mean())
            if zcol in sub.columns:
                row[f"{t} z"] = float(sub[zcol].mean())
            scol = f"{t}_placement_suspect"
            if scol in sub.columns:
                row[f"{t} unplaced"] = int(sub[scol].sum())
        grows.append(row)
    if grows:
        lines.append(md_table(pd.DataFrame(grows)) + "\n")
        lines.append(
            "Chemotype membership is decided by the SMARTS patterns in "
            "`common.CHEMOTYPE_SMARTS`. `unplaced` counts dockings whose best score "
            "was weaker than the plausibility floor, i.e. the ligand could not be "
            "placed in the box -- those are docking failures to investigate, not "
            "evidence of a steric clash.\n")

    cand_cols = [c for c in ["identifier", "name", "mean_pred", "np_likeness", "QED",
                             "PAINS", "pan_kinase_alert"] if c in df.columns]
    cand_cols += [c for t in C.TARGETS for c in (f"{t}_vina_kcal", f"{t}_vina_sd")
                  if c in df.columns]
    nov = df[df["category"] == "novel_consensus"] if "category" in df.columns else df
    lines.append("### 4.4 Candidate shortlist\n")
    lines.append(md_table(nov[cand_cols]) + "\n")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 5. Limitations, generated from what the pipeline actually measured
# ----------------------------------------------------------------------------

def limitations_section():
    lines = ["## 5. Limitations\n"]
    items = []

    for t in C.TARGETS:
        cur = C.load_json(os.path.join(C.METRICS_DIR, f"curation_{t}.json"))
        if cur and cur.get("censored_records_available") is False:
            items.append(
                f"No inactive class for {t}: the curated set contains no `>` "
                f"records, so reported AUC/MCC describe ranking inside a "
                f"{cur.get('exact_active_prevalence_at_pIC50_6', float('nan')):.0%}-active "
                f"set rather than active-versus-inactive discrimination.")
            break

    r = next((b for b in (C.load_json(os.path.join(C.METRICS_DIR,
                                                   f"final_model_{t}.json"))
                          for t in C.TARGETS)
              if b and "conformal_validity" in b), None)
    if r:
        v = r.get("conformal_validity", {})
        items.append(
            f"Conformal coverage is guaranteed only for {v.get('holds_for', '')}. "
            f"It does not hold for {v.get('does_not_hold_for', '')}")
        npt = r.get("np_transfer", {})
        if npt.get("np_specific_half_width_90") is None:
            items.append(npt.get("status", ""))
        for t in C.TARGETS:
            rr = C.load_json(os.path.join(C.METRICS_DIR, f"final_model_{t}.json"))
            if rr and rr.get("test_NP_like", {}).get("n", 0) < 50:
                items.append(
                    f"The NP-like held-out subset is small "
                    f"(n = {rr['test_NP_like']['n']} for {t}); its metrics carry wide "
                    f"confidence intervals and support a direction, not a point claim.")
                break

    for t in C.TARGETS:
        rr = C.load_json(os.path.join(C.METRICS_DIR, f"final_model_{t}.json"))
        if rr and rr.get("tau_is_binding") is False:
            items.append(
                f"The applicability-domain threshold is non-binding for {t}: "
                f"{rr.get('tau_interpretation', '')} Any in-domain hit count for "
                f"this target is therefore a count of the whole library, not of a "
                f"calibrated-reliable subset.")

    s = C.load_json(os.path.join(C.METRICS_DIR, "screen_summary.json"))
    if s and "interval_caveat" not in s:
        items.append(
            "The library screen and docking validation in this repository predate "
            "the current pipeline and are withheld from this report; sections 3 "
            "and 4 name the commands that regenerate them. No screening or "
            "docking claim in this document is current.")
        s = None
    if s:
        cons = s.get("consensus", {})
        if cons.get("consensus_caveat"):
            items.append(cons["consensus_caveat"])
        for t, d in (s.get("targets") or {}).items():
            rho = d.get("spearman_model_vs_nn_baseline_in_AD")
            if rho is not None and rho > 0.9:
                items.append(
                    f"For {t} the model ranks in-domain library compounds almost "
                    f"identically to a nearest-neighbour baseline (Spearman "
                    f"{rho:.2f}), so the learned model adds little over similarity "
                    f"search on this library.")
                break

    dv = C.load_json(os.path.join(C.METRICS_DIR, "docking_validation.json"))
    if dv and "score_comparability_note" not in dv:
        dv = None
    if dv:
        fails = [t for t, rc in (dv.get("redocking_control") or {}).items()
                 if rc.get("pose_prediction_success_top1") is False]
        if fails:
            items.append(
                f"Top-1 redocking fails the 2 A criterion for {', '.join(fails)}, so "
                f"pose-level conclusions in those pockets are not supported.")
        if dv.get("score_comparability_note"):
            items.append(dv["score_comparability_note"])

    items.append(
        "No experimental validation. Every candidate is a computational "
        "prediction and the paper claims nothing beyond that.")

    for i in items:
        if i:
            lines.append(f"- {i}\n")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

def generate_figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG_DIR, exist_ok=True)

    # Fig 1: canonical-split R2 with bootstrap CIs.
    canonical = C.load_json(os.path.join(C.METRICS_DIR, "benchmark_canonical.json"))
    if canonical:
        splits = canonical.get("splits", {})
        models = sorted({m for b in splits.values() for m in b.get("models", {})})
        fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
        x = np.arange(len(models))
        width = 0.8 / max(1, len(splits))
        for j, (t, b) in enumerate(splits.items()):
            vals, los, his = [], [], []
            for m in models:
                e = b["models"].get(m, {})
                v = e.get("R2", np.nan)
                ci = e.get("R2_CI95") or {}
                vals.append(v)
                los.append(v - (ci.get("lo") if ci.get("lo") is not None else v))
                his.append((ci.get("hi") if ci.get("hi") is not None else v) - v)
            ax.bar(x + (j - (len(splits) - 1) / 2) * width, vals, width,
                   yerr=[los, his], capsize=3, label=t)
        ax.set_xticks(x)
        ax.set_xticklabels(models, rotation=25, ha="right")
        ax.set_ylabel("$R^2$ (canonical scaffold split)")
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
        ax.set_title("Scaffold-split performance with 95% bootstrap CIs")
        ax.legend()
        fig.savefig(os.path.join(FIG_DIR, "fig1_benchmark.png"), dpi=200)
        plt.close(fig)

    # Fig 2: conformal coverage and half-width per similarity bin.
    blocks = {t: C.load_json(os.path.join(C.METRICS_DIR, f"final_model_{t}.json"))
              for t in C.TARGETS}
    blocks = {t: b for t, b in blocks.items() if b}
    if blocks:
        fig, axes = plt.subplots(1, len(blocks), figsize=(4.4 * len(blocks), 4.2),
                                 constrained_layout=True, squeeze=False)
        for j, (t, r) in enumerate(blocks.items()):
            bt = pd.DataFrame(r["bin_table"])
            ax = axes[0][j]
            xi = np.arange(len(bt))
            ax.plot(xi, bt["half_width_90"], "o-", label="90% half-width")
            if "half_width_90_isotonic" in bt.columns:
                ax.plot(xi, bt["half_width_90_isotonic"], "^:",
                        label="isotonic fit (sets tau)")
            ax.plot(xi, bt["test_MAE"], "s--", label="held-out MAE")
            ax.set_xticks(xi)
            ax.set_xticklabels([f"[{a:.1f},{min(b_, 1.0):.1f})"
                                for a, b_ in zip(bt.sim_lo, bt.sim_hi)],
                               rotation=25, ha="right")
            ax.set_ylabel("pIC50 units")
            ax.set_xlabel("max Tanimoto to training set")
            cov = r["overall_coverage"]
            ax.set_title(f"{t}: coverage {cov:.2f}, tau={r['tau']}")
            ax.legend(fontsize=8)
            ax2 = ax.twinx()
            ax2.plot(xi, bt["test_coverage"], "d-", color="gray", alpha=0.6)
            ax2.axhline(1 - r["alpha"], color="gray", linestyle=":", linewidth=0.8)
            ax2.set_ylim(0, 1.05)
            ax2.set_ylabel("empirical coverage", color="gray")
        fig.savefig(os.path.join(FIG_DIR, "fig2_domain_shift.png"), dpi=200)
        plt.close(fig)

    # Fig 3: redocking top-1 vs best-of-N, and enrichment with CIs.
    dv = C.load_json(os.path.join(C.METRICS_DIR, "docking_validation.json"))
    if dv:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
        rc = dv.get("redocking_control", {})
        targets = [t for t in C.TARGETS if t in rc]
        if targets:
            ax = axes[0]
            xx = np.arange(len(targets))
            top1 = [rc[t].get("rmsd_top1_A") or np.nan for t in targets]
            best = [rc[t].get("rmsd_best_of_n_A") or np.nan for t in targets]
            ax.bar(xx - 0.19, top1, 0.38, label="top-1 (pose prediction)")
            ax.bar(xx + 0.19, best, 0.38, label="best-of-N (sampling)")
            ax.axhline(2.0, color="crimson", linestyle="--", linewidth=1,
                       label="2 A success criterion")
            ax.set_xticks(xx)
            ax.set_xticklabels(targets)
            ax.set_ylabel("Heavy-atom RMSD (A, symmetry-corrected)")
            ax.set_title("Control 1: co-crystal redocking")
            ax.legend(fontsize=8)

        ec = dv.get("enrichment_control", {})
        usable = [t for t in C.TARGETS
                  if isinstance(ec.get(t), dict) and ec[t].get("status") == "ok"]
        ax = axes[1]
        if usable:
            xx = np.arange(len(usable))
            aucs = [ec[t]["docking_ROC_AUC"] for t in usable]
            los = [a - (ec[t].get("docking_ROC_AUC_CI95") or {}).get("lo", a)
                   for a, t in zip(aucs, usable)]
            his = [(ec[t].get("docking_ROC_AUC_CI95") or {}).get("hi", a) - a
                   for a, t in zip(aucs, usable)]
            ax.bar(xx, aucs, 0.5, yerr=[los, his], capsize=4)
            ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
            ax.set_xticks(xx)
            ax.set_xticklabels([f"{t}\nn={ec[t]['n_actives']}/{ec[t]['n_inactives']}"
                                for t in usable])
            ax.set_ylim(0, 1.0)
            ax.set_ylabel("Docking ROC-AUC")
            ax.set_title("Control 2: active-vs-inactive enrichment (95% CI)")
        else:
            ax.text(0.5, 0.5, "enrichment control not available\n"
                              "(no inactive class in curated data)",
                    ha="center", va="center")
            ax.set_axis_off()
        fig.savefig(os.path.join(FIG_DIR, "fig3_docking_validation.png"), dpi=200)
        plt.close(fig)


def main():
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    parts = [
        "# HCC kinase natural-product screen: results report\n",
        "_Auto-generated from `results/metrics` and `results/predictions`. Every "
        "number, table and caveat below is read from those files; nothing is "
        "hand-typed._\n",
        curation_section(),
        benchmark_section(),
        domain_section(),
        screen_section(),
        docking_section(),
        limitations_section(),
    ]
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))
    generate_figures()
    print(f"[OK] report written: {REPORT}")
    return {"status": "report generated"}


if __name__ == "__main__":
    main()
