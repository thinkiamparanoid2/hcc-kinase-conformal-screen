"""
Step 03 - Rigorous benchmark (replaces the old single-split training script)
-----------------------------------------------------------------------------
For every target, 5 repeated Murcko-scaffold splits (80/20). Models:

  mean          train-set mean (floor)
  kNN-Tanimoto  similarity-weighted 5-NN on ECFP4 (simple, strong QSAR baseline)
  RF-desc       Random Forest on the 21 physicochemical descriptors only
  RF-ecfp       Random Forest on the 2048-bit ECFP4 only
  RF-all        Random Forest on descriptors + ECFP4
  XGB-all       XGBoost on descriptors + ECFP4
  XGB-scramble  XGBoost on shuffled labels (must give R2 <= 0, AUC ~ 0.5)

Regression metrics: R2, RMSE, MAE, Spearman. Ranking view (label = pIC50 >= 6):
ROC-AUC, PR-AUC, MCC. F1 is intentionally not reported (~75% of labels are active).

Resumable: each (target, seed) is appended to results/metrics/benchmark_long.jsonl
and skipped on re-run. Output: benchmark_long.csv, benchmark_summary.csv.
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
import common as C
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

JSONL = os.path.join(C.METRICS_DIR, "benchmark_long.jsonl")


def make_models(seed):
    return {
        "RF-desc": ("desc", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.5, min_samples_leaf=1, n_jobs=-1, random_state=seed)),
        "RF-ecfp": ("fp", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.3, n_jobs=-1, random_state=seed)),
        "RF-all": ("all", lambda: RandomForestRegressor(
            n_estimators=300, max_features=0.3, n_jobs=-1, random_state=seed)),
        "XGB-all": ("all", lambda: XGBRegressor(
            n_estimators=500, max_depth=7, learning_rate=0.05, subsample=0.8,
            colsample_bytree=0.5, tree_method="hist", n_jobs=-1, random_state=seed, verbosity=0)),
    }


def done_keys():
    if not os.path.exists(JSONL):
        return set()
    with open(JSONL) as f:
        return {(r["target"], r["seed"]) for r in map(json.loads, f)}


def run_target_seed(target, seed, df, cols, fp_cols, desc_cols, groups, fps, X_all):
    tr, te = C.scaffold_split(groups, 0.2, seed)
    y = df["pIC50"].to_numpy(float)
    y_tr, y_te = y[tr], y[te]
    col_idx = {"all": np.arange(len(cols)),
               "fp": np.array([cols.index(c) for c in fp_cols]),
               "desc": np.array([cols.index(c) for c in desc_cols])}
    rows = []

    def record(name, pred, secs, extra=None):
        r = {"target": target, "seed": seed, "model": name, "n_train": len(tr),
             "n_test": len(te), "test_frac": len(te) / len(y), "secs": round(secs, 1)}
        r.update(C.regression_metrics(y_te, pred)); r.update(C.ranking_metrics(y_te, pred))
        if extra: r.update(extra)
        rows.append(r)
        print(f"    {name:13s} R2={r['R2']:.3f} RMSE={r['RMSE']:.3f} Sp={r['Spearman']:.3f} "
              f"AUC={r['AUC']:.3f} PR={r['PR_AUC']:.3f} MCC={r['MCC']:.3f} ({secs:.0f}s)", flush=True)

    record("mean", np.full(len(te), y_tr.mean()), 0.0)

    t0 = time.time()
    fps_tr = [fps[i] for i in tr]; fps_te = [fps[i] for i in te]
    pred, mx, mk = C.knn_tanimoto(fps_te, fps_tr, y_tr, k=5)
    record("kNN-Tanimoto", pred, time.time() - t0,
           {"test_maxsim_median": float(np.median(mx))})

    for name, (view, ctor) in make_models(seed).items():
        t0 = time.time()
        idx = col_idx[view]
        Xtr, Xte = X_all[tr][:, idx], X_all[te][:, idx]
        if name.startswith("RF"):
            Xtr, Xte = C.median_impute(Xtr, Xte)
        m = ctor(); m.fit(Xtr, y_tr)
        record(name, m.predict(Xte), time.time() - t0)

    t0 = time.time()
    y_shuf = np.random.RandomState(seed).permutation(y_tr)
    m = make_models(seed)["XGB-all"][1](); m.fit(X_all[tr], y_shuf)
    record("XGB-scramble", m.predict(X_all[te]), time.time() - t0)
    return rows


def summarize():
    rows = [json.loads(l) for l in open(JSONL)]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(C.METRICS_DIR, "benchmark_long.csv"), index=False)
    metrics = ["R2", "RMSE", "MAE", "Spearman", "AUC", "PR_AUC", "MCC"]
    g = df.groupby(["target", "model"])[metrics]
    summ = g.mean().round(3).astype(str) + " +/- " + g.std().round(3).astype(str)
    summ.insert(0, "n_seeds", df.groupby(["target", "model"]).size())
    summ.to_csv(os.path.join(C.METRICS_DIR, "benchmark_summary.csv"))
    # paired comparison: XGB-all vs kNN baseline, per target (R2 difference, wins out of n seeds)
    pair = []
    for t, d in df.groupby("target"):
        a = d[d.model == "XGB-all"].set_index("seed")["R2"]
        b = d[d.model == "kNN-Tanimoto"].set_index("seed")["R2"]
        diff = (a - b).dropna()
        pair.append({"target": t, "XGB_minus_kNN_R2_mean": float(diff.mean()),
                     "XGB_wins": int((diff > 0).sum()), "n_seeds": int(len(diff))})
    C.save_json(pair, os.path.join(C.METRICS_DIR, "benchmark_paired.json"))
    print("\n=== SUMMARY (mean +/- SD over seeds) ===")
    print(summ.to_string())
    print("\nPaired XGB-all vs kNN-Tanimoto:", json.dumps(pair, indent=1))


def main(targets=None, seeds=None):
    os.makedirs(C.METRICS_DIR, exist_ok=True)
    targets = targets or C.TARGETS
    seeds = C.SEEDS if seeds is None else seeds
    finished = done_keys()
    for t in targets:
        df, cols, fp_cols, desc_cols = C.load_training(t)
        print(f"\n=== {t}: {len(df)} molecules, {len(cols)} features ===", flush=True)
        groups = C.murcko_groups(df["canonical_smiles"])
        print(f"    {len(set(groups))} generic Murcko scaffolds", flush=True)
        fps = C.morgan_bitvects(df["canonical_smiles"])
        X_all = C.clean_matrix(df, cols)
        for s in seeds:
            if (t, s) in finished:
                print(f"  seed {s}: already done, skipping"); continue
            print(f"  seed {s}:", flush=True)
            rows = run_target_seed(t, s, df, cols, fp_cols, desc_cols, groups, fps, X_all)
            with open(JSONL, "a") as f:
                for r in rows:
                    f.write(json.dumps(r) + "\n")
    summarize()
    return {"status": "benchmark complete"}


if __name__ == "__main__":
    # optional smoke test:  python scripts/03_scaffold_benchmark.py FGFR4 0
    if len(sys.argv) > 2:
        main([sys.argv[1]], [int(sys.argv[2])])
    else:
        main()
