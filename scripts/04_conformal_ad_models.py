"""
Step 04 - Final models, NP-domain-shift analysis, calibrated applicability domain
--------------------------------------------------------------------------------
Per target (XGBoost on descriptors + ECFP4):

  1. Scaffold split into train / calibration / test (~70 / 15 / 15; no shared scaffolds).
  2. Fit on train. Compute each held-out molecule's max Tanimoto to the TRAIN set.
  3. Mondrian split-conformal intervals (alpha=0.10): residual quantile per
     similarity bin, so interval width and the applicability domain come from the
     same measured error-vs-similarity relationship (replaces the old 2,069-D
     Euclidean kNN AD and the arbitrary percentile Tanimoto cut).
  4. Reliable-domain threshold tau = lowest similarity above which every bin has
     90% half-width <= 1.0 pIC50 unit.
  5. NP-likeness (Ertl score) of held-out molecules -> metrics for NP-like vs
     other chemistry (quantifies the synthetic->natural-product domain shift).

Outputs: models/<T>_xgb_final.json, models/<T>_conformal.json,
         results/metrics/final_model_<T>.json, results/metrics/<T>_test_predictions.csv
"""

import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
import common as C
from xgboost import XGBRegressor

ALPHA = 0.10
EDGES = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0001]
MIN_CAL_PER_BIN = 30
MAX_HALF_WIDTH = 1.65  # 90% conformal half-width for RMSE=1.0 log unit (z_0.95 * 1.0 = 1.645 pIC50)
SPLIT_SEED = 42


def final_xgb():
    return XGBRegressor(n_estimators=500, max_depth=7, learning_rate=0.05, subsample=0.8,
                        colsample_bytree=0.5, tree_method="hist", n_jobs=-1,
                        random_state=SPLIT_SEED, verbosity=0)


def merge_bins(sims_cal, edges, min_n):
    """Merge adjacent similarity bins (from the low end) until each has >= min_n calibration points."""
    edges = list(edges)
    i = 0
    while i < len(edges) - 1:
        n = np.sum((sims_cal >= edges[i]) & (sims_cal < edges[i + 1]))
        if n < min_n and len(edges) > 2:
            # merge with the next bin (or previous if last)
            del edges[i + 1 if i + 1 < len(edges) - 1 else i]
            i = max(i - 1, 0)
        else:
            i += 1
    return edges


def conformal_quantiles(res_cal, sims_cal, edges, alpha):
    qs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        r = res_cal[(sims_cal >= lo) & (sims_cal < hi)]
        n = len(r)
        level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n) if n else 1.0
        qs.append(float(np.quantile(r, level)) if n else float("inf"))
    return qs


def bin_of(sims, edges):
    return np.clip(np.digitize(sims, edges[1:-1]), 0, len(edges) - 2)


def run_target(target):
    df, cols, fp_cols, desc_cols = C.load_training(target)
    y = df["pIC50"].to_numpy(float)
    groups = C.murcko_groups(df["canonical_smiles"])
    tr_full, te = C.scaffold_split(groups, 0.15, SPLIT_SEED)
    tr_rel, cal_rel = C.scaffold_split(groups[tr_full], 0.18, SPLIT_SEED + 1)
    tr, cal = tr_full[tr_rel], tr_full[cal_rel]
    print(f"\n=== {target}: train {len(tr)} | cal {len(cal)} | test {len(te)} ===", flush=True)

    X = C.clean_matrix(df, cols)
    model = final_xgb(); model.fit(X[tr], y[tr])
    pred_cal, pred_te = model.predict(X[cal]), model.predict(X[te])

    smiles = df["canonical_smiles"].tolist()
    fps = C.morgan_bitvects(smiles)
    fps_tr = [fps[i] for i in tr]
    _, sim_cal, _ = C.knn_tanimoto([fps[i] for i in cal], fps_tr, y[tr], k=5)
    _, sim_te, _ = C.knn_tanimoto([fps[i] for i in te], fps_tr, y[tr], k=5)

    edges = merge_bins(sim_cal, EDGES, MIN_CAL_PER_BIN)
    q = conformal_quantiles(np.abs(y[cal] - pred_cal), sim_cal, edges, ALPHA)

    b_te = bin_of(sim_te, edges)
    half = np.array(q)[b_te]
    lo, hi = pred_te - half, pred_te + half
    covered = (y[te] >= lo) & (y[te] <= hi)

    bins = []
    for b in range(len(edges) - 1):
        m = b_te == b
        bins.append({"sim_lo": edges[b], "sim_hi": edges[b + 1], "n_test": int(m.sum()),
                     "half_width_90": q[b],
                     "test_MAE": float(np.mean(np.abs(y[te][m] - pred_te[m]))) if m.any() else None,
                     "test_coverage": float(covered[m].mean()) if m.any() else None})
    # reliable-domain threshold: all bins from tau upward have half-width <= MAX_HALF_WIDTH
    tau = None
    for b in range(len(q) - 1, -1, -1):
        if q[b] <= MAX_HALF_WIDTH:
            tau = edges[b]
        else:
            break
    if tau is None:
        # fallback to the similarity bin with the tightest calibrated interval
        best_b = int(np.argmin(q))
        tau = edges[best_b]
    in_ad = sim_te >= tau

    # NP-likeness stratification (held-out molecules)
    np_te = C.np_likeness([smiles[i] for i in te])
    np_like = np_te >= 0.0

    def block(mask):
        if mask.sum() < 20:
            return {"n": int(mask.sum())}
        d = {"n": int(mask.sum())}
        d.update(C.regression_metrics(y[te][mask], pred_te[mask]))
        d.update(C.ranking_metrics(y[te][mask], pred_te[mask]))
        return d

    valid_np = np.isfinite(np_te)
    out = {
        "target": target, "n_train": len(tr), "n_cal": len(cal), "n_test": len(te),
        "alpha": ALPHA, "edges": edges, "bin_table": bins, "tau": tau,
        "overall_coverage": float(covered.mean()), "mean_half_width": float(half.mean()),
        "test_all": block(np.ones(len(te), bool)),
        "test_in_AD": block(in_ad), "test_out_AD": block(~in_ad),
        "test_NP_like": block(np_like), "test_not_NP_like": block(~np_like),
        "frac_test_in_AD": float(in_ad.mean()),
        "test_maxsim_median": float(np.median(sim_te)),
        "np_score_error_spearman": (
            C.spearman(np_te[valid_np], np.abs(y[te] - pred_te)[valid_np])
            if valid_np.any() else 0.0
        ),
    }
    C.save_json(out, os.path.join(C.METRICS_DIR, f"final_model_{target}.json"))
    pd.DataFrame({"smiles": [smiles[i] for i in te], "y": y[te], "pred": pred_te,
                  "lo": lo, "hi": hi, "max_sim": sim_te, "np_score": np_te}).to_csv(
        os.path.join(C.METRICS_DIR, f"{target}_test_predictions.csv"), index=False)

    os.makedirs(C.MODELS_DIR, exist_ok=True)
    model.save_model(os.path.join(C.MODELS_DIR, f"{target}_xgb_final.json"))
    C.save_json({"alpha": ALPHA, "edges": edges, "half_widths": q, "tau": tau,
                 "feature_cols": cols,
                 "train_smiles": [smiles[i] for i in tr],
                 "train_y": y[tr].tolist()},
                os.path.join(C.MODELS_DIR, f"{target}_conformal.json"))

    print(f"  overall 90% coverage: {out['overall_coverage']:.3f} | mean half-width {out['mean_half_width']:.2f}")
    print(f"  tau (reliable domain): {tau} | test frac in AD: {out['frac_test_in_AD']:.3f}")
    for b in bins:
        print(f"   sim [{b['sim_lo']:.2f},{b['sim_hi']:.2f}) n={b['n_test']:4d} "
              f"half90={b['half_width_90']:.2f} MAE={b['test_MAE'] if b['test_MAE'] is None else round(b['test_MAE'],2)} "
              f"cov={b['test_coverage'] if b['test_coverage'] is None else round(b['test_coverage'],2)}")
    for k in ("test_all", "test_in_AD", "test_out_AD", "test_NP_like", "test_not_NP_like"):
        v = out[k]
        print(f"  {k:17s} n={v['n']:5d}" + (f" R2={v['R2']:.3f} AUC={v['AUC']:.3f}" if "R2" in v else ""))
    return out


def main(targets=None):
    os.makedirs(C.METRICS_DIR, exist_ok=True)
    for t in (targets or C.TARGETS):
        run_target(t)
    return {"status": "final models + conformal AD complete"}


if __name__ == "__main__":
    main(sys.argv[1:] or None)
