"""
Shared utilities for the HCC natural-product screening pipeline.
----------------------------------------------------------------
One place for: feature-column detection, Murcko-scaffold splitting, fingerprints,
Tanimoto-kNN baseline, metrics, and NP-likeness scoring. Import with:

    sys.path.insert(0, "scripts"); import common

Windows note: avoid large numpy.dot / matmul here (OpenBLAS aborts seen on this
machine). Everything below uses RDKit bulk ops, trees, or elementwise numpy.
"""

import os
import sys

_BIN_DIR = os.path.join(sys.prefix, "Library", "bin")
if os.path.isdir(_BIN_DIR):
    os.environ["PATH"] = _BIN_DIR + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_BIN_DIR)

import json
import warnings
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, average_precision_score, matthews_corrcoef)

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

TARGETS = ["FGFR4", "EGFR", "VEGFR2"]
ACTIVE_THRESHOLD = 6.0          # pIC50 >= 6.0  <=>  IC50 <= 1 uM
SEEDS = list(range(10))         # repeated scaffold splits
FEATURES_DIR = os.path.join("data", "features")
METRICS_DIR = os.path.join("results", "metrics")
MODELS_DIR = "models"

_EXCLUDE = {"ligand_id", "canonical_smiles", "pIC50", "pchembl_value",
            "molecule_chembl_id", "scaffold", "scaffold_group"}


# ----------------------------------------------------------------------------
# Data / features
# ----------------------------------------------------------------------------

def load_training(target):
    """Training feature table with a valid label. Returns (df, feature_cols, fp_cols, desc_cols)."""
    df = pd.read_csv(os.path.join(FEATURES_DIR, f"{target}_training_features.csv"))
    df = df.dropna(subset=["pIC50", "canonical_smiles"]).reset_index(drop=True)
    cols = [c for c in df.columns if c not in _EXCLUDE]
    fp_cols = [c for c in cols if c.startswith("FP_")]
    desc_cols = [c for c in cols if not c.startswith("FP_")]
    return df, cols, fp_cols, desc_cols


def clean_matrix(df, cols):
    """float32 matrix; inf -> NaN (XGBoost handles NaN natively; RF gets median-imputed by caller)."""
    X = df[cols].to_numpy(dtype=np.float32)
    X[~np.isfinite(X)] = np.nan
    return X


def median_impute(X_train, *others):
    med = np.nanmedian(X_train, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    out = []
    for X in (X_train,) + others:
        X = X.copy()
        idx = np.where(np.isnan(X))
        X[idx] = med[idx[1]]
        out.append(X)
    return out if len(out) > 1 else out[0]


# ----------------------------------------------------------------------------
# Scaffold split
# ----------------------------------------------------------------------------

def murcko_groups(smiles):
    """Integer group id per molecule from the GENERIC Murcko scaffold (strictest grouping)."""
    ids, lookup = [], {}
    for smi in smiles:
        key = "UNKNOWN"
        try:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                key = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(
                    MurckoScaffold.GetScaffoldForMol(m)))
        except Exception:
            pass
        ids.append(lookup.setdefault(key, len(lookup)))
    return np.asarray(ids)


def scaffold_split(groups, test_frac, seed):
    """Random scaffold-group assignment until test_frac is reached. No scaffold in both sets.
    Returns (train_idx, test_idx)."""
    rng = np.random.RandomState(seed)
    uniq, counts = np.unique(groups, return_counts=True)
    order = rng.permutation(len(uniq))
    target = test_frac * len(groups)
    test_groups, n = set(), 0
    for i in order:
        if n >= target:
            break
        # skip groups that would overshoot massively (keeps the split near test_frac)
        if n + counts[i] > target * 1.25 and n > 0:
            continue
        test_groups.add(uniq[i]); n += counts[i]
    mask = np.isin(groups, list(test_groups))
    return np.where(~mask)[0], np.where(mask)[0]


# ----------------------------------------------------------------------------
# Fingerprints, Tanimoto kNN
# ----------------------------------------------------------------------------

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def morgan_bitvect(smiles):
    m = Chem.MolFromSmiles(smiles)
    return None if m is None else _MORGAN.GetFingerprint(m)


def morgan_bitvects(smiles_list):
    return [morgan_bitvect(s) for s in smiles_list]


def knn_tanimoto(query_fps, train_fps, train_y, k=5, exclude_self=False):
    """Similarity-weighted kNN regression on ECFP4 Tanimoto.
    Returns (pred, max_sim, mean_topk_sim) arrays over queries."""
    train_y = np.asarray(train_y, dtype=float)
    pred = np.empty(len(query_fps)); mx = np.empty(len(query_fps)); mk = np.empty(len(query_fps))
    for i, q in enumerate(query_fps):
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(q, train_fps))
        if exclude_self:
            sims[i] = -1.0
        top = np.argpartition(-sims, k)[:k]
        w = np.clip(sims[top], 1e-6, None)
        pred[i] = np.sum(w * train_y[top]) / np.sum(w)
        mx[i] = sims[top].max(); mk[i] = sims[top].mean()
    return pred, mx, mk


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def pearson(a, b):
    """Elementwise Pearson r (no BLAS: scipy.stats pearsonr/spearmanr abort silently on this Windows setup)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    den = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / den) if den > 0 else 0.0


def spearman(a, b):
    from scipy.stats import rankdata
    return pearson(rankdata(a), rankdata(b))


def regression_metrics(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    return {
        "R2": float(r2_score(y, p)),
        "RMSE": float(np.sqrt(mean_squared_error(y, p))),
        "MAE": float(mean_absolute_error(y, p)),
        "Spearman": spearman(y, p),
        "Pearson": pearson(y, p),
    }



def ranking_metrics(y, p, thr=ACTIVE_THRESHOLD):
    """Classification view of a regression model: label = y >= thr, score = predicted pIC50."""
    y_bin = (np.asarray(y) >= thr).astype(int)
    p = np.asarray(p, float)
    out = {"prevalence": float(y_bin.mean())}
    if y_bin.min() == y_bin.max():
        out.update({"AUC": float("nan"), "PR_AUC": float("nan"), "MCC": float("nan")})
        return out
    out["AUC"] = float(roc_auc_score(y_bin, p))
    out["PR_AUC"] = float(average_precision_score(y_bin, p))
    out["MCC"] = float(matthews_corrcoef(y_bin, (p >= thr).astype(int)))
    return out


# ----------------------------------------------------------------------------
# NP-likeness (Ertl & Schuffenhauer, RDKit Contrib)
# ----------------------------------------------------------------------------

_NP_MODEL = None


def np_likeness(smiles_list):
    """Natural-product-likeness score per SMILES (higher = more NP-like; NaN if unparsable)."""
    global _NP_MODEL
    sys.path.append(os.path.join(RDConfig.RDContribDir, "NP_Score"))
    import npscorer
    if _NP_MODEL is None:
        _NP_MODEL = npscorer.readNPModel()
    out = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        out.append(npscorer.scoreMol(m, _NP_MODEL) if m is not None else np.nan)
    return np.asarray(out)


def save_json(obj, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=float)
