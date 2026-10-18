"""
Step 04 - Final model, Mondrian conformal calibration, applicability domain
---------------------------------------------------------------------------
Per target:

  1. Scaffold-disjoint train / calibration / test split (~70 / 15 / 15).
  2. Fit the model family that step 03 selected (results/metrics/model_selection.json).
  3. Mondrian split-conformal intervals (alpha = 0.10) with the Mondrian taxonomy
     defined by max-Tanimoto-to-train bins.
  4. Reliable-domain threshold `tau`, from an ISOTONIC fit of half-width against
     similarity, with a bootstrap CI and an explicit sensitivity analysis.
  5. Natural-product transfer analysis, including NP-specific recalibration when
     the calibration set contains enough NP-like chemistry.

WHAT CHANGED AND WHY (these are the claims that decide the paper):

* EXCHANGEABILITY IS NOW STATED AND TESTED, NOT ASSUMED.
  Split-conformal coverage is guaranteed only when calibration and test data are
  exchangeable. Here calibration is ChEMBL synthetic kinase chemistry and
  deployment (step 05) is a natural-product library, which is a covariate shift.
  The guarantee therefore does NOT transfer, and this script now measures the
  damage instead of implying it away: coverage is reported separately on NP-like
  held-out molecules, an NP-specific conformal calibration is fitted when there is
  enough NP-like calibration data, and `conformal_validity` in the output states
  in plain words which claims the numbers support.

* BINS ARE WIDE ENOUGH FOR A 90% QUANTILE.
  MIN_CAL_PER_BIN rose from 30 to 150. At n = 30 the finite-sample conformal level
  ceil((n+1)(1-alpha))/n is 1.0, i.e. the quantile is the maximum residual -- the
  single noisiest statistic available. Every bin now reports its own coverage with
  a bootstrap CI so under-coverage is visible rather than averaged away.

* TAU NO LONGER ASSUMES MONOTONICITY IT DOES NOT HAVE.
  The old rule walked bins from the top and stopped at the first one exceeding the
  width budget. Empirical half-widths are not monotone in similarity (noise), so a
  single noisy bin could discard usable lower bins. Half-width is now fitted
  against similarity with isotonic regression (non-increasing), tau is read off the
  fit, and both its bootstrap CI and a tau +/- 0.1 sensitivity table are reported.

Outputs: models/<T>_final.json (+ .pkl for sklearn families),
         models/<T>_conformal.json,
         results/metrics/final_model_<T>.json,
         results/metrics/<T>_test_predictions.csv
"""

import os
import sys
import json

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C


ALPHA = 0.10
EDGES = [0.0, 0.4, 0.5, 0.6, 0.7, 1.0001]
MIN_CAL_PER_BIN = 150        # was 30; a 90% quantile from 30 points is the maximum
MIN_NP_CAL = 60              # minimum NP-like calibration points for NP-specific bins
MAX_HALF_WIDTH = 1.65        # 90% half-width consistent with RMSE <= 1.0 log unit
                             # (z_0.95 * 1.0 = 1.645 pIC50 units)
NP_LIKENESS_CUTOFF = 0.0     # Ertl NP-likeness >= 0 counts as NP-like
SPLIT_SEED = 42
TAU_BOOTSTRAP = 400


def _load_benchmark_module():
    """Import step 03 by file path (its module name starts with a digit)."""
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "03_scaffold_benchmark.py")
    spec = importlib.util.spec_from_file_location("benchmark03", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def selected_model_name():
    sel = C.load_json(os.path.join(C.METRICS_DIR, "model_selection.json"))
    if sel and sel.get("selected_model"):
        return sel["selected_model"]
    print("[!] No model_selection.json found; falling back to RF-ecfp. "
          "Run step 03 first so the screening model is the one the benchmark favours.")
    return "RF-ecfp"


# ----------------------------------------------------------------------------
# Mondrian conformal machinery
# ----------------------------------------------------------------------------

def merge_bins(sims_cal, edges, min_n):
    """Merge adjacent similarity bins until each holds >= min_n calibration points."""
    edges = list(edges)
    i = 0
    while i < len(edges) - 1:
        n = np.sum((sims_cal >= edges[i]) & (sims_cal < edges[i + 1]))
        if n < min_n and len(edges) > 2:
            del edges[i + 1 if i + 1 < len(edges) - 1 else i]
            i = max(i - 1, 0)
        else:
            i += 1
    return edges


def conformal_quantiles(res_cal, sims_cal, edges, alpha):
    """Finite-sample conformal quantile of |residual| per similarity bin."""
    qs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        r = res_cal[(sims_cal >= lo) & (sims_cal < hi)]
        n = len(r)
        if n == 0:
            qs.append(float("inf"))
            continue
        level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
        qs.append(float(np.quantile(r, level)))
    return qs


def bin_of(sims, edges):
    return np.clip(np.digitize(sims, edges[1:-1]), 0, len(edges) - 2)


def isotonic_tau(edges, half_widths, max_half_width):
    """Reliable-domain threshold from a non-increasing isotonic fit.

    Half-width should fall as similarity to the training set rises; empirically it
    does so with noise. Fitting the monotone trend and reading tau off the fit stops
    one noisy bin from discarding usable lower bins, which is what the previous
    top-down scan did.

    Returns (tau, fitted_half_widths). tau is None when no bin meets the budget.
    """
    from sklearn.isotonic import IsotonicRegression
    mids = np.array([(edges[i] + min(edges[i + 1], 1.0)) / 2.0
                     for i in range(len(edges) - 1)])
    hw = np.asarray(half_widths, float)
    finite = np.isfinite(hw)
    if finite.sum() < 2:
        return None, hw
    weights = np.ones(finite.sum())
    iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
    fitted_sub = iso.fit_transform(mids[finite], hw[finite], sample_weight=weights)
    fitted = hw.copy()
    fitted[finite] = fitted_sub
    ok = np.where(np.isfinite(fitted) & (fitted <= max_half_width))[0]
    if ok.size == 0:
        return None, fitted
    return float(edges[ok[0]]), fitted


def tau_by_budget(edges, half_widths, budgets):
    """tau at several width budgets, so the threshold's provenance is visible.

    MAX_HALF_WIDTH is a choice, and it decides how much of a library the domain
    admits. Reporting tau across a range of budgets stops one chosen number from
    reading as a property of the data.
    """
    out = []
    for budget in budgets:
        t, _ = isotonic_tau(edges, half_widths, budget)
        out.append({"max_half_width": budget, "tau": t,
                    "binding": None if t is None else bool(t > edges[0])})
    return out


def bootstrap_tau(res_cal, sims_cal, edges, alpha, max_half_width,
                  n_boot=TAU_BOOTSTRAP, seed=0):
    """Percentile CI for tau by resampling the calibration set."""
    rng = np.random.RandomState(seed)
    n = len(res_cal)
    draws = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        q = conformal_quantiles(res_cal[idx], sims_cal[idx], edges, alpha)
        t, _ = isotonic_tau(edges, q, max_half_width)
        if t is not None:
            draws.append(t)
    if not draws:
        return {"lo": None, "hi": None, "n_boot": 0,
                "frac_no_reliable_domain": 1.0}
    draws = np.array(draws)
    return {"lo": float(np.quantile(draws, 0.025)),
            "hi": float(np.quantile(draws, 0.975)),
            "median": float(np.median(draws)),
            "n_boot": int(len(draws)),
            "frac_no_reliable_domain": float(1.0 - len(draws) / n_boot)}


# ----------------------------------------------------------------------------
# Per-target run
# ----------------------------------------------------------------------------

def _subset_block(y, pred, mask, seed=0):
    out = {"n": int(mask.sum())}
    if mask.sum() < 20:
        out["note"] = ("fewer than 20 molecules: no metric is reported because the "
                       "confidence interval would span most of the plausible range")
        return out
    out.update(C.regression_metrics(y[mask], pred[mask]))
    out.update(C.ranking_metrics(y[mask], pred[mask]))
    out["R2_CI95"] = C.r2_ci(y[mask], pred[mask], seed=seed)
    return out


def run_target(target, model_name):
    bench = _load_benchmark_module()
    df, cols, fp_cols, desc_cols = C.load_training(target)
    y = df["pIC50"].to_numpy(float)
    smiles = df["canonical_smiles"].tolist()
    groups = C.murcko_groups(smiles)

    tr_full, te = C.scaffold_split(groups, 0.15, SPLIT_SEED, scheme="deterministic")
    tr_rel, cal_rel = C.scaffold_split(groups[tr_full], 0.18, SPLIT_SEED + 1,
                                       scheme="random")
    tr, cal = tr_full[tr_rel], tr_full[cal_rel]
    assert not (set(groups[tr]) & set(groups[cal])), "train/cal scaffold overlap"
    assert not (set(groups[tr]) & set(groups[te])), "train/test scaffold overlap"
    assert not (set(groups[cal]) & set(groups[te])), "cal/test scaffold overlap"
    print(f"\n=== {target} [{model_name}]: train {len(tr)} | cal {len(cal)} | "
          f"test {len(te)} ===", flush=True)

    X = C.clean_matrix(df, cols)
    view = bench.model_view(model_name, SPLIT_SEED)
    pos = {c: i for i, c in enumerate(cols)}
    idx = {"all": np.arange(len(cols)),
           "fp": np.array([pos[c] for c in fp_cols]),
           "desc": np.array([pos[c] for c in desc_cols])}[view]
    Xtr, Xcal, Xte = X[tr][:, idx], X[cal][:, idx], X[te][:, idx]
    if model_name.startswith("RF"):
        Xtr, Xcal, Xte = C.median_impute(Xtr, Xcal, Xte)

    model = bench.build_model(model_name, SPLIT_SEED)
    model.fit(Xtr, y[tr])
    pred_cal, pred_te = model.predict(Xcal), model.predict(Xte)

    fps = C.morgan_bitvects(smiles)
    fps_tr = [fps[i] for i in tr]
    _, sim_cal, _ = C.knn_tanimoto([fps[i] for i in cal], fps_tr, y[tr], k=5)
    _, sim_te, _ = C.knn_tanimoto([fps[i] for i in te], fps_tr, y[tr], k=5)

    # ---- Mondrian conformal on similarity bins ----
    edges = merge_bins(sim_cal, EDGES, MIN_CAL_PER_BIN)
    res_cal = np.abs(y[cal] - pred_cal)
    q = conformal_quantiles(res_cal, sim_cal, edges, ALPHA)

    b_te = bin_of(sim_te, edges)
    half = np.array(q)[b_te]
    lo, hi = pred_te - half, pred_te + half
    covered = (y[te] >= lo) & (y[te] <= hi)

    bins = []
    for b in range(len(edges) - 1):
        m_cal = (sim_cal >= edges[b]) & (sim_cal < edges[b + 1])
        m = b_te == b
        entry = {"sim_lo": edges[b], "sim_hi": edges[b + 1],
                 "n_cal": int(m_cal.sum()), "n_test": int(m.sum()),
                 "half_width_90": q[b],
                 "conformal_level_used": float(min(
                     1.0, np.ceil((m_cal.sum() + 1) * (1 - ALPHA)) / m_cal.sum()))
                 if m_cal.sum() else None,
                 "test_MAE": float(np.mean(np.abs(y[te][m] - pred_te[m]))) if m.any() else None,
                 "test_coverage": float(covered[m].mean()) if m.any() else None}
        if m.sum() >= 20:
            entry["test_coverage_CI95"] = C.coverage_ci(covered[m])
        bins.append(entry)

    tau, fitted_hw = isotonic_tau(edges, q, MAX_HALF_WIDTH)
    tau_ci = bootstrap_tau(res_cal, sim_cal, edges, ALPHA, MAX_HALF_WIDTH)
    budget_sweep = tau_by_budget(edges, q, [1.0, 1.25, 1.5, MAX_HALF_WIDTH, 2.0])
    # tau only does work when it excludes something. When every bin meets the
    # budget, tau sits at the lowest edge and the domain admits everything -- the
    # gate is vacuous, and saying so is the honest reading. It is not a licence to
    # narrow the bins until tau looks useful.
    tau_is_binding = bool(tau is not None and tau > edges[0])
    for b, fit in zip(bins, fitted_hw):
        b["half_width_90_isotonic"] = float(fit) if np.isfinite(fit) else None
    in_ad = sim_te >= tau if tau is not None else np.zeros(len(te), bool)

    # tau sensitivity: how much does the domain, and its accuracy, move?
    sensitivity = []
    for delta in (-0.1, -0.05, 0.0, 0.05, 0.1):
        if tau is None:
            break
        t2 = float(np.clip(tau + delta, 0.0, 1.0))
        m = sim_te >= t2
        row = {"tau": round(t2, 3), "delta": delta,
               "frac_test_in_AD": float(m.mean()), "n_test_in_AD": int(m.sum())}
        if m.sum() >= 20:
            row["R2_in_AD"] = float(C.regression_metrics(y[te][m], pred_te[m])["R2"])
            row["coverage_in_AD"] = float(covered[m].mean())
        sensitivity.append(row)

    # ---- natural-product transfer ----
    np_te = C.np_likeness([smiles[i] for i in te])
    np_cal = C.np_likeness([smiles[i] for i in cal])
    np_like_te = np_te >= NP_LIKENESS_CUTOFF
    np_like_cal = np_cal >= NP_LIKENESS_CUTOFF

    np_block = {
        "n_cal_NP_like": int(np_like_cal.sum()),
        "n_test_NP_like": int(np_like_te.sum()),
        "cutoff": NP_LIKENESS_CUTOFF,
    }
    if np_like_te.sum() >= 20:
        np_block["coverage_on_NP_like_test"] = float(covered[np_like_te].mean())
        np_block["coverage_on_NP_like_test_CI95"] = C.coverage_ci(covered[np_like_te])
        np_block["coverage_on_non_NP_like_test"] = float(covered[~np_like_te].mean())
        np_block["coverage_on_non_NP_like_test_CI95"] = C.coverage_ci(covered[~np_like_te])
    # NP-specific recalibration: a conformal quantile fitted on NP-like calibration
    # molecules only. This is the honest way to carry a guarantee into NP chemistry:
    # the taxonomy must contain the stratum you intend to deploy on.
    if np_like_cal.sum() >= MIN_NP_CAL:
        r = res_cal[np_like_cal]
        level = min(1.0, np.ceil((len(r) + 1) * (1 - ALPHA)) / len(r))
        q_np = float(np.quantile(r, level))
        np_block["np_specific_half_width_90"] = q_np
        if np_like_te.sum() >= 20:
            cov_np = (np.abs(y[te][np_like_te] - pred_te[np_like_te]) <= q_np)
            np_block["np_specific_coverage_on_NP_like_test"] = float(cov_np.mean())
            np_block["np_specific_coverage_CI95"] = C.coverage_ci(cov_np)
        np_block["status"] = "NP-specific Mondrian stratum fitted"
    else:
        np_block["np_specific_half_width_90"] = None
        np_block["status"] = (
            f"insufficient NP-like calibration molecules "
            f"({int(np_like_cal.sum())} < {MIN_NP_CAL}): no NP-specific conformal "
            f"stratum can be fitted from ChEMBL alone. Intervals transferred to a "
            f"natural-product library therefore carry NO coverage guarantee, only "
            f"an extrapolation of the synthetic-chemistry error model.")

    # NP-likeness as a shift diagnostic rather than a second small-n R2 claim.
    valid_np = np.isfinite(np_te)
    shift = {
        "np_score_vs_abs_error_spearman": (
            C.spearman(np_te[valid_np], np.abs(y[te] - pred_te)[valid_np])
            if valid_np.sum() > 2 else None),
        "median_np_score_train": float(np.nanmedian(C.np_likeness([smiles[i] for i in tr[:2000]]))),
        "median_np_score_test": float(np.nanmedian(np_te)),
        "median_max_sim_test": float(np.median(sim_te)),
    }

    overall_cov = float(covered.mean())
    validity = {
        "nominal_coverage": 1 - ALPHA,
        "empirical_coverage_held_out_chembl": overall_cov,
        "empirical_coverage_CI95": C.coverage_ci(covered),
        "holds_for": "held-out ChEMBL scaffolds, which are exchangeable with the "
                     "calibration set",
        "does_not_hold_for": "natural-product libraries screened in step 05: the "
                             "covariate shift between synthetic calibration "
                             "chemistry and natural products breaks exchangeability, "
                             "so intervals attached to library compounds are "
                             "uncalibrated extrapolations and must be described as "
                             "such wherever they are reported",
    }

    out = {
        "target": target,
        "model": model_name,
        "n_train": len(tr), "n_cal": len(cal), "n_test": len(te),
        "alpha": ALPHA, "edges": edges, "bin_table": bins,
        "tau": tau, "tau_CI95": tau_ci, "tau_sensitivity": sensitivity,
        "tau_is_binding": tau_is_binding,
        "tau_by_budget": budget_sweep,
        "tau_interpretation": (
            "tau excludes part of the similarity range; compounds below it fall "
            "outside the applicability domain."
            if tau_is_binding else
            f"tau sits at the lowest similarity bin: at a {MAX_HALF_WIDTH} pIC50 "
            f"width budget every calibration bin qualifies, so the applicability "
            f"domain excludes nothing and the gate is vacuous for this target. "
            f"With this calibration set, bins wide enough for a stable 90% "
            f"quantile are too wide to resolve a threshold. See tau_by_budget for "
            f"what a tighter budget would give."),
        "max_half_width_budget": MAX_HALF_WIDTH,
        "min_cal_per_bin": MIN_CAL_PER_BIN,
        "overall_coverage": overall_cov,
        "mean_half_width": float(half.mean()),
        "conformal_validity": validity,
        "test_all": _subset_block(y[te], pred_te, np.ones(len(te), bool)),
        "test_in_AD": _subset_block(y[te], pred_te, in_ad),
        "test_out_AD": _subset_block(y[te], pred_te, ~in_ad),
        "test_NP_like": _subset_block(y[te], pred_te, np_like_te),
        "test_not_NP_like": _subset_block(y[te], pred_te, ~np_like_te),
        "np_transfer": np_block,
        "domain_shift_diagnostics": shift,
        "frac_test_in_AD": float(in_ad.mean()),
        "test_active_prevalence": float((y[te] >= C.ACTIVE_THRESHOLD).mean()),
    }
    C.save_json(out, os.path.join(C.METRICS_DIR, f"final_model_{target}.json"))

    pd.DataFrame({"smiles": [smiles[i] for i in te], "y": y[te], "pred": pred_te,
                  "lo": lo, "hi": hi, "max_sim": sim_te, "np_score": np_te,
                  "in_AD": in_ad, "covered": covered}).to_csv(
        os.path.join(C.METRICS_DIR, f"{target}_test_predictions.csv"), index=False)

    # ---- persist model + calibration ----
    os.makedirs(C.MODELS_DIR, exist_ok=True)
    if hasattr(model, "save_model"):            # XGBoost: compact, portable JSON
        model_path = os.path.join(C.MODELS_DIR, f"{target}_final.json")
        model.save_model(model_path)
    else:                                       # scikit-learn
        # A 300-tree forest over 2048 ECFP bits pickles to 60-170 MB raw, which
        # exceeds GitHub's file limit, so these are compressed and git-ignored and
        # regenerated by re-running this step. joblib is also the format sklearn
        # documents for estimator persistence.
        import joblib
        model_path = os.path.join(C.MODELS_DIR, f"{target}_final.joblib")
        joblib.dump(model, model_path, compress=3)
    C.save_json({
        "model_family": model_name, "model_path": model_path, "feature_view": view,
        "alpha": ALPHA, "edges": edges, "half_widths": q, "tau": tau,
        "tau_CI95": tau_ci, "tau_is_binding": tau_is_binding,
        "np_specific_half_width_90": np_block.get("np_specific_half_width_90"),
        "feature_cols": cols,
        "train_smiles": [smiles[i] for i in tr],
        "train_y": y[tr].tolist(),
        "validity": validity,
    }, os.path.join(C.MODELS_DIR, f"{target}_conformal.json"))

    print(f"  coverage {overall_cov:.3f} "
          f"[{validity['empirical_coverage_CI95']['lo']:.3f}, "
          f"{validity['empirical_coverage_CI95']['hi']:.3f}] (nominal {1-ALPHA})")
    print(f"  tau = {tau} (95% CI {tau_ci.get('lo')}-{tau_ci.get('hi')}) | "
          f"{'binding' if tau_is_binding else 'NON-BINDING: excludes nothing'} | "
          f"test in AD {out['frac_test_in_AD']*100:.1f}%")
    print("  tau by width budget: "
          + ", ".join(f"{b['max_half_width']}->{b['tau']}" for b in budget_sweep))
    print(f"  NP transfer: {np_block['status']}")
    for b in bins:
        cov = b["test_coverage"]
        print(f"   sim [{b['sim_lo']:.2f},{min(b['sim_hi'],1.0):.2f}) "
              f"n_cal={b['n_cal']:4d} n_test={b['n_test']:4d} "
              f"half90={b['half_width_90']:.2f} "
              f"cov={'n/a' if cov is None else round(cov, 2)}")
    return out


def main(targets=None):
    os.makedirs(C.METRICS_DIR, exist_ok=True)
    model_name = selected_model_name()
    print(f"=== Step 04: calibrating {model_name} (chosen by step 03) ===")
    for t in (targets or C.TARGETS):
        run_target(t, model_name)
    return {"status": "final models + conformal AD complete", "model": model_name}


if __name__ == "__main__":
    main(sys.argv[1:] or None)
