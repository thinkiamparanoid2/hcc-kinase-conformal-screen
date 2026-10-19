"""
Step 03 - Scaffold-split benchmark and model selection
------------------------------------------------------
Two evaluations per target, because they answer different questions:

A. CANONICAL SPLIT (deterministic Bemis-Murcko, largest scaffold families into
   train). Seed-independent, comparable with published benchmarks. Every headline
   number from this split is reported with a percentile bootstrap CI over the test
   molecules, so the paper never states an R2 or an AUC without its uncertainty.

B. REPEATED RANDOM SCAFFOLD SPLITS (`common.SEEDS`, 10 by default). These measure
   how much the result moves when the scaffold assignment changes. The spread is
   SPLIT variance, not sampling variance -- the same molecules are reused every
   time -- and the paired Wilcoxon test over seeds is therefore a statement about
   split noise, not about generalisation to new chemistry. Reported as such.

C. TEMPORAL HOLDOUT, when step 01 recorded publication years: train on the older
   literature, test on the newest `TEMPORAL_TEST_FRAC`. This is the closest cheap
   stand-in for prospective validation and is the evaluation reviewers increasingly
   ask for. Skipped with a recorded reason when years are unavailable.

Models:
  mean          train-set mean (floor)
  kNN-Tanimoto  similarity-weighted 5-NN on ECFP4 (strong, cheap QSAR baseline)
  RF-desc       Random Forest on the 21 physicochemical descriptors
  RF-ecfp       Random Forest on the 2048-bit ECFP4
  RF-all        Random Forest on descriptors + ECFP4
  XGB-all       XGBoost on descriptors + ECFP4
  XGB-scramble  XGBoost on shuffled labels (leakage control: must give R2 <= 0)

Hyperparameters are fixed library defaults, NOT tuned -- see HYPERPARAMETER_NOTE.
Only uncensored (point-valued) labels are used for regression.

MODEL SELECTION: the family with the best mean R2 on the repeated splits is
written to results/metrics/model_selection.json, and step 04 reads that file to
pick the model it calibrates and screens with. Previously the benchmark favoured
RF-ecfp while the screen ran on XGBoost; the choice is now made by the data.

Outputs: results/metrics/benchmark_long.{jsonl,csv}, benchmark_summary.csv,
         benchmark_paired.json, benchmark_canonical.json, model_selection.json
"""

import os
import sys
import json
import time
import hashlib

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

JSONL = os.path.join(C.METRICS_DIR, "benchmark_long.jsonl")
TEMPORAL_TEST_FRAC = 0.2

HYPERPARAMETER_NOTE = (
    "Hyperparameters are fixed a priori at common library defaults and were not "
    "tuned on any split reported here; no model selection was performed against "
    "test-set performance. Tuning, if added, must run as an inner cross-validation "
    "loop on the training scaffolds only."
)

# Model family used for the final calibrated model if selection cannot be run.
FALLBACK_MODEL = "RF-ecfp"


def make_models(seed):
    return {
        "RF-desc": ("desc", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.5, min_samples_leaf=1,
            n_jobs=-1, random_state=seed)),
        "RF-ecfp": ("fp", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.3, n_jobs=-1, random_state=seed)),
        "RF-all": ("all", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.3, n_jobs=-1, random_state=seed)),
        "XGB-all": ("all", lambda: XGBRegressor(
            n_estimators=500, max_depth=7, learning_rate=0.05, subsample=0.8,
            colsample_bytree=0.5, tree_method="hist", n_jobs=-1,
            random_state=seed, verbosity=0)),
    }


def build_model(name, seed):
    """Public constructor so step 04 can rebuild whichever family won."""
    if name in make_models(seed):
        return make_models(seed)[name][1]()
    raise ValueError(f"unknown model family {name!r}")


def model_view(name, seed):
    return make_models(seed)[name][0]


def data_fingerprint(target):
    """Short digest of the feature table a benchmark row was computed from."""
    path = os.path.join(C.FEATURES_DIR, f"{target}_training_features.csv")
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def done_keys():
    """Rows already computed, keyed by the data they came from as well as the seed.

    Keying on (target, seed) alone made the checkpoint lie after step 01 or 02
    changed the data: every seed was reported "already done" and skipped, and
    summarize() then averaged rows from two different datasets into one
    model-selection decision. A row with no `data_sha` predates this field and
    can never match, so it is treated as stale.
    """
    if not os.path.exists(JSONL):
        return set()
    with open(JSONL) as f:
        return {(r["target"], r["seed"], r.get("data_sha")) for r in map(json.loads, f)}


def prune_stale(fingerprints):
    """Drop rows whose feature table has changed, and return how many went.

    summarize() averages whatever the file holds, so stale rows have to leave the
    file rather than merely be re-run alongside.
    """
    if not os.path.exists(JSONL):
        return 0
    with open(JSONL) as f:
        rows = [json.loads(l) for l in f]
    keep = [r for r in rows
            if r["target"] not in fingerprints
            or r.get("data_sha") == fingerprints[r["target"]]]
    if len(keep) == len(rows):
        return 0
    with open(JSONL, "w") as f:
        for r in keep:
            f.write(json.dumps(r) + "\n")
    return len(rows) - len(keep)


def _fit_predict(name, seed, X_all, col_idx, tr, te, y_tr):
    view, ctor = make_models(seed)[name]
    idx = col_idx[view]
    Xtr, Xte = X_all[tr][:, idx], X_all[te][:, idx]
    if name.startswith("RF"):
        Xtr, Xte = C.median_impute(Xtr, Xte)
    model = ctor()
    model.fit(Xtr, y_tr)
    return model.predict(Xte)


def _col_index(cols, fp_cols, desc_cols):
    pos = {c: i for i, c in enumerate(cols)}
    return {"all": np.arange(len(cols)),
            "fp": np.array([pos[c] for c in fp_cols]),
            "desc": np.array([pos[c] for c in desc_cols])}


def evaluate_split(target, seed, tr, te, y, fps, X_all, col_idx, split_name):
    """Run every model on one split; return a list of metric rows."""
    y_tr, y_te = y[tr], y[te]
    rows = []

    def record(name, pred, secs, extra=None):
        r = {"target": target, "seed": seed, "split": split_name, "model": name,
             "n_train": len(tr), "n_test": len(te), "test_frac": len(te) / len(y),
             "secs": round(secs, 1)}
        r.update(C.regression_metrics(y_te, pred))
        r.update(C.ranking_metrics(y_te, pred))
        if extra:
            r.update(extra)
        rows.append(r)
        print(f"    {name:13s} R2={r['R2']:.3f} RMSE={r['RMSE']:.3f} "
              f"Sp={r['Spearman']:.3f} AUC={r['AUC']:.3f} MCC={r['MCC']:.3f} "
              f"({secs:.0f}s)", flush=True)
        return r

    record("mean", np.full(len(te), y_tr.mean()), 0.0)

    t0 = time.time()
    pred, mx, _ = C.knn_tanimoto([fps[i] for i in te], [fps[i] for i in tr], y_tr, k=5)
    record("kNN-Tanimoto", pred, time.time() - t0,
           {"test_maxsim_median": float(np.median(mx))})

    for name in make_models(seed):
        t0 = time.time()
        record(name, _fit_predict(name, seed, X_all, col_idx, tr, te, y_tr),
               time.time() - t0)

    t0 = time.time()
    y_shuf = np.random.RandomState(seed).permutation(y_tr)
    scrambled = make_models(seed)["XGB-all"][1]()
    scrambled.fit(X_all[tr], y_shuf)
    record("XGB-scramble", scrambled.predict(X_all[te]), time.time() - t0)
    return rows


def canonical_split_report(target, y, fps, X_all, col_idx, groups):
    """Deterministic scaffold split with bootstrap CIs on every headline metric."""
    tr, te = C.scaffold_split(groups, 0.2, seed=0, scheme="deterministic")
    y_te = y[te]
    y_bin = (y_te >= C.ACTIVE_THRESHOLD).astype(int)
    out = {"target": target, "n_train": int(len(tr)), "n_test": int(len(te)),
           "test_active_prevalence": float(y_bin.mean()),
           "n_scaffolds_train": int(len(set(groups[tr]))),
           "n_scaffolds_test": int(len(set(groups[te]))),
           "scaffold_overlap": int(len(set(groups[tr]) & set(groups[te]))),
           "models": {}}
    for name in list(make_models(0)) + ["kNN-Tanimoto"]:
        if name == "kNN-Tanimoto":
            pred, _, _ = C.knn_tanimoto([fps[i] for i in te], [fps[i] for i in tr],
                                        y[tr], k=5)
        else:
            pred = _fit_predict(name, 0, X_all, col_idx, tr, te, y[tr])
        entry = C.regression_metrics(y_te, pred)
        entry.update(C.ranking_metrics(y_te, pred))
        entry["R2_CI95"] = C.r2_ci(y_te, pred)
        if np.unique(y_bin).size > 1:
            entry["AUC_CI95"] = C.auc_ci(y_bin, pred)
        out["models"][name] = entry
        print(f"    [canonical] {name:13s} R2={entry['R2']:.3f} "
              f"[{entry['R2_CI95']['lo']:.3f}, {entry['R2_CI95']['hi']:.3f}]", flush=True)
    return out


def temporal_split_report(target, df, y, fps, X_all, col_idx):
    """Train on older literature, test on the newest TEMPORAL_TEST_FRAC."""
    if "first_year" not in df.columns:
        return {"status": "skipped", "reason": "no publication year in curated data"}
    years = pd.to_numeric(df["first_year"], errors="coerce").to_numpy(float)
    if np.isfinite(years).sum() < 200:
        return {"status": "skipped",
                "reason": f"only {int(np.isfinite(years).sum())} molecules carry a year"}
    valid = np.where(np.isfinite(years))[0]
    cutoff = float(np.quantile(years[valid], 1 - TEMPORAL_TEST_FRAC))
    tr = valid[years[valid] < cutoff]
    te = valid[years[valid] >= cutoff]
    if len(tr) < 100 or len(te) < 50:
        return {"status": "skipped", "reason": "temporal split too small"}
    out = {"status": "ok", "cutoff_year": cutoff,
           "n_train": int(len(tr)), "n_test": int(len(te)),
           "test_active_prevalence": float((y[te] >= C.ACTIVE_THRESHOLD).mean()),
           "models": {}}
    for name in list(make_models(0)) + ["kNN-Tanimoto"]:
        if name == "kNN-Tanimoto":
            pred, _, _ = C.knn_tanimoto([fps[i] for i in te], [fps[i] for i in tr],
                                        y[tr], k=5)
        else:
            pred = _fit_predict(name, 0, X_all, col_idx, tr, te, y[tr])
        entry = C.regression_metrics(y[te], pred)
        entry.update(C.ranking_metrics(y[te], pred))
        entry["R2_CI95"] = C.r2_ci(y[te], pred)
        out["models"][name] = entry
        print(f"    [temporal] {name:13s} R2={entry['R2']:.3f}", flush=True)
    return out


def summarize():
    rows = [json.loads(l) for l in open(JSONL)]
    df = pd.DataFrame(rows)
    df = df[df.get("split", "repeated") == "repeated"] if "split" in df.columns else df
    df.to_csv(os.path.join(C.METRICS_DIR, "benchmark_long.csv"), index=False)
    metrics = ["R2", "RMSE", "MAE", "Spearman", "AUC", "PR_AUC", "MCC"]
    g = df.groupby(["target", "model"])[metrics]
    summ = g.mean().round(3).astype(str) + " +/- " + g.std().round(3).astype(str)
    summ.insert(0, "n_seeds", df.groupby(["target", "model"]).size())
    summ.to_csv(os.path.join(C.METRICS_DIR, "benchmark_summary.csv"))

    pair = []
    try:
        from scipy.stats import wilcoxon
    except Exception:
        wilcoxon = None
    for t, d in df.groupby("target"):
        piv = d.pivot_table(index="seed", columns="model", values="R2")
        for a, b in [("XGB-all", "kNN-Tanimoto"), ("RF-ecfp", "kNN-Tanimoto"),
                     ("RF-all", "RF-ecfp"), ("XGB-all", "RF-ecfp")]:
            if a not in piv.columns or b not in piv.columns:
                continue
            diff = (piv[a] - piv[b]).dropna()
            p = float("nan")
            if wilcoxon is not None and len(diff) >= 6 and (diff != 0).any():
                p = float(wilcoxon(diff).pvalue)
            pair.append({"target": t, "comparison": f"{a} vs {b}",
                         "mean_dR2": float(diff.mean()),
                         "wins": int((diff > 0).sum()), "n_seeds": int(len(diff)),
                         "wilcoxon_p": p,
                         "interpretation": "paired over scaffold splits of one "
                                           "dataset; measures split noise, not "
                                           "generalisation to new chemistry"})
    C.save_json(pair, os.path.join(C.METRICS_DIR, "benchmark_paired.json"))

    # ---- model selection: best mean R2 across the repeated splits ----
    candidates = [m for m in df.model.unique()
                  if m not in ("mean", "XGB-scramble", "kNN-Tanimoto")]
    per_target, scores = {}, {}
    for m in candidates:
        means = []
        for t in df.target.unique():
            sub = df[(df.target == t) & (df.model == m)]["R2"]
            if len(sub):
                means.append(float(sub.mean()))
                per_target.setdefault(t, {})[m] = float(sub.mean())
        scores[m] = float(np.mean(means)) if means else float("-inf")
    winner = max(scores, key=scores.get) if scores else FALLBACK_MODEL
    selection = {
        "selected_model": winner,
        "mean_R2_by_model": scores,
        "mean_R2_by_target": per_target,
        "criterion": "highest mean R2 across targets on repeated random scaffold splits",
        "note": "Step 04 calibrates and step 05 screens with this family, so the "
                "model used downstream is the one this benchmark actually favours.",
        "hyperparameters": HYPERPARAMETER_NOTE,
    }
    C.save_json(selection, os.path.join(C.METRICS_DIR, "model_selection.json"))

    print("\n=== SUMMARY (mean +/- SD over repeated scaffold splits) ===")
    print(summ.to_string())
    print(f"\nSelected model family for steps 04-05: {winner}")
    return selection


def main(targets=None, seeds=None):
    os.makedirs(C.METRICS_DIR, exist_ok=True)
    targets = targets or C.TARGETS
    seeds = C.SEEDS if seeds is None else seeds
    fingerprints = {t: data_fingerprint(t) for t in targets}
    dropped = prune_stale(fingerprints)
    if dropped:
        print(f"[!] dropped {dropped} benchmark rows computed from a feature table "
              f"that has since changed; those seeds will be re-run", flush=True)
    finished = done_keys()
    canonical, temporal = {}, {}

    for t in targets:
        df, cols, fp_cols, desc_cols = C.load_training(t)
        print(f"\n=== {t}: {len(df)} molecules with point-valued labels, "
              f"{len(cols)} features ===", flush=True)
        groups = C.murcko_groups(df["canonical_smiles"])
        print(f"    {len(set(groups))} generic Murcko scaffolds | "
              f"active prevalence {float((df['pIC50'] >= C.ACTIVE_THRESHOLD).mean()):.3f}",
              flush=True)
        fps = C.morgan_bitvects(df["canonical_smiles"])
        X_all = C.clean_matrix(df, cols)
        y = df["pIC50"].to_numpy(float)
        col_idx = _col_index(cols, fp_cols, desc_cols)

        for s in seeds:
            if (t, s, fingerprints[t]) in finished:
                print(f"  seed {s}: already done, skipping")
                continue
            print(f"  seed {s}:", flush=True)
            tr, te = C.scaffold_split(groups, 0.2, s, scheme="random")
            rows = evaluate_split(t, s, tr, te, y, fps, X_all, col_idx, "repeated")
            with open(JSONL, "a") as f:
                for r in rows:
                    r["data_sha"] = fingerprints[t]
                    f.write(json.dumps(r) + "\n")

        print(f"  canonical deterministic scaffold split:", flush=True)
        canonical[t] = canonical_split_report(t, y, fps, X_all, col_idx, groups)
        temporal[t] = temporal_split_report(t, df, y, fps, X_all, col_idx)

    C.save_json({"splits": canonical,
                 "note": "Deterministic Bemis-Murcko split (largest scaffold "
                         "families in train). CIs are percentile bootstrap over "
                         "test molecules.",
                 "hyperparameters": HYPERPARAMETER_NOTE},
                os.path.join(C.METRICS_DIR, "benchmark_canonical.json"))
    C.save_json({"splits": temporal,
                 "note": "Temporal holdout: train on older literature, test on the "
                         "most recent 20% by first publication year."},
                os.path.join(C.METRICS_DIR, "benchmark_temporal.json"))
    return summarize()


if __name__ == "__main__":
    # smoke test:  python scripts/03_scaffold_benchmark.py FGFR4 0
    if len(sys.argv) > 2:
        main([sys.argv[1]], [int(sys.argv[2])])
    else:
        main()
