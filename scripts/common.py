"""
Shared utilities for the HCC natural-product screening pipeline.
----------------------------------------------------------------
One place for: structure standardisation, descriptor definitions, feature-column
detection, Murcko-scaffold splitting, fingerprints, Tanimoto-kNN, metrics with
bootstrap confidence intervals, external-binary resolution, and NP-likeness.

Import with:

    sys.path.insert(0, "scripts"); import common

Platform note: external binaries (AutoDock Vina, Open Babel) are resolved via
`resolve_binary()`, which honours environment variables first and then PATH, so
the pipeline runs on Linux, macOS and Windows. Nothing here hard-codes a path.
"""

import os
import shutil
import sys

# On a Windows conda install the RDKit/Open Babel DLLs live in
# <prefix>/Library/bin and are not on PATH by default. Harmless elsewhere.
_WIN_BIN_DIR = os.path.join(sys.prefix, "Library", "bin")
if os.path.isdir(_WIN_BIN_DIR):
    os.environ["PATH"] = _WIN_BIN_DIR + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_WIN_BIN_DIR)

import json
import warnings
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, average_precision_score, matthews_corrcoef)

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

TARGETS = ["FGFR4", "EGFR", "VEGFR2"]
TARGET_CHEMBL_IDS = {"EGFR": "CHEMBL203", "FGFR4": "CHEMBL1974", "VEGFR2": "CHEMBL279"}
ACTIVE_THRESHOLD = 6.0          # pIC50 >= 6.0  <=>  IC50 <= 1 uM
SEEDS = list(range(10))         # repeated scaffold splits
N_BOOTSTRAP = 2000              # resamples for every reported confidence interval
DATA_DIR = "data"
LIGAND_DIR = os.path.join(DATA_DIR, "ligands_chembl")
FEATURES_DIR = os.path.join(DATA_DIR, "features")
METRICS_DIR = os.path.join("results", "metrics")
MODELS_DIR = "models"

# Columns that are metadata, never model inputs.
_EXCLUDE = {"ligand_id", "canonical_smiles", "standardised_smiles", "inchikey",
            "pIC50", "pchembl_value", "molecule_chembl_id", "scaffold",
            "scaffold_group", "relation", "censored", "n_replicates",
            "replicate_sd", "standard_value", "standard_units",
            "standard_relation", "label_kind", "first_year", "original_smiles",
            "n_censored_records"}


# ----------------------------------------------------------------------------
# External binaries (cross-platform)
# ----------------------------------------------------------------------------

def resolve_binary(name, env_var=None, required=True):
    """Locate an external executable.

    Resolution order: explicit environment variable, then PATH (which includes
    <sys.prefix>/Library/bin on Windows conda installs). Returns None when the
    binary is absent and `required` is False, so callers can degrade gracefully
    instead of crashing at import time.
    """
    if env_var:
        override = os.environ.get(env_var)
        if override:
            if os.path.isfile(override):
                return override
            found = shutil.which(override)
            if found:
                return found
            raise FileNotFoundError(f"{env_var}={override!r} is not an executable")
    found = shutil.which(name)
    if found:
        return found
    if required:
        raise FileNotFoundError(
            f"Could not find {name!r} on PATH. Install it (see environment.yml) or "
            f"set {env_var or name.upper()} to its full path.")
    return None


# ----------------------------------------------------------------------------
# Structure standardisation
# ----------------------------------------------------------------------------

_LARGEST_FRAGMENT = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()
_TAUTOMER = rdMolStandardize.TautomerEnumerator()


def standardise_mol(smiles):
    """Salt-strip -> neutralise -> canonical tautomer.

    Returns (standardised_smiles, inchikey) or (None, None) when the input cannot
    be parsed. Deduplicating on the returned InChIKey collapses salt/free-base and
    tautomer pairs that raw canonical SMILES keeps apart -- those pairs otherwise
    receive different Murcko scaffolds and can leak across a scaffold split.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = _LARGEST_FRAGMENT.choose(mol)
        mol = _UNCHARGER.uncharge(mol)
        mol = _TAUTOMER.Canonicalize(mol)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol), Chem.MolToInchiKey(mol)
    except Exception:
        return None, None


def standardise_series(smiles_iterable):
    """Vectorised `standardise_mol` -> (list[smiles|None], list[inchikey|None])."""
    out_smi, out_key = [], []
    for s in smiles_iterable:
        smi, key = standardise_mol(s)
        out_smi.append(smi)
        out_key.append(key)
    return out_smi, out_key


# ----------------------------------------------------------------------------
# Descriptors (single source of truth; step 02 and step 05 both use these)
# ----------------------------------------------------------------------------

PHYSCHEM_DESCRIPTOR_NAMES = [
    "MolWt", "MolLogP", "TPSA",
    "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "NumAromaticRings",
    "NumHeteroatoms", "NumAliphaticRings",
    "FractionCSP3", "HeavyAtomCount",
    "RingCount", "BertzCT",
    "NumValenceElectrons",
    "NumRadicalElectrons",
    "LabuteASA",          # Labute accessible surface area
    "BalabanJ",           # Balaban's J index (topological)
    "Chi0n", "Chi1n",     # Kier-Hall connectivity indices
    "HallKierAlpha",      # Hall-Kier alpha value
]
DESC_COLS = PHYSCHEM_DESCRIPTOR_NAMES + ["LipinskiViolations"]
N_FP_BITS = 2048
FP_COLS = [f"FP_{i}" for i in range(N_FP_BITS)]
FEATURE_COLS = DESC_COLS + FP_COLS

_DESC_CALC = None


def descriptor_calculator():
    global _DESC_CALC
    if _DESC_CALC is None:
        from rdkit.ML.Descriptors import MoleculeDescriptors
        _DESC_CALC = MoleculeDescriptors.MolecularDescriptorCalculator(
            PHYSCHEM_DESCRIPTOR_NAMES)
    return _DESC_CALC


def physchem_descriptors(mol):
    """21 physicochemical descriptors (20 RDKit + Lipinski violation count)."""
    values = dict(zip(PHYSCHEM_DESCRIPTOR_NAMES,
                      descriptor_calculator().CalcDescriptors(mol)))
    values["LipinskiViolations"] = int(
        (values.get("MolWt", 0) > 500) + (values.get("MolLogP", 0) > 5) +
        (values.get("NumHDonors", 0) > 5) + (values.get("NumHAcceptors", 0) > 10))
    return values


# ----------------------------------------------------------------------------
# Data / features
# ----------------------------------------------------------------------------

def load_training(target, exact_only=True):
    """Training feature table. Returns (df, feature_cols, fp_cols, desc_cols).

    `exact_only=True` keeps only uncensored (`standard_relation == '='`) records,
    which are the only ones with a point-valued pIC50 and therefore the only ones
    admissible as regression targets. Censored records (`<`, `<=`, `>`, `>=`) are
    retained in the file and reachable with `exact_only=False` for the
    classification and docking-enrichment analyses, where a bound is informative.
    """
    df = pd.read_csv(os.path.join(FEATURES_DIR, f"{target}_training_features.csv"))
    df = df.dropna(subset=["pIC50", "canonical_smiles"]).reset_index(drop=True)
    if exact_only and "censored" in df.columns:
        df = df[~df["censored"].astype(bool)].reset_index(drop=True)
    cols = [c for c in df.columns if c not in _EXCLUDE]
    fp_cols = [c for c in cols if c.startswith("FP_")]
    desc_cols = [c for c in cols if not c.startswith("FP_")]
    return df, cols, fp_cols, desc_cols


def clean_matrix(df, cols):
    """float32 matrix; inf -> NaN (XGBoost handles NaN natively; RF gets imputed)."""
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
    """Integer group id per molecule from the GENERIC Murcko scaffold (strictest)."""
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


def scaffold_split(groups, test_frac, seed, scheme="deterministic"):
    """Scaffold-disjoint train/test split. No scaffold appears in both sets.

    scheme="deterministic" (default) is the conventional Bemis-Murcko split: groups
    are sorted largest-first and filled into train until the train quota is met,
    so the biggest scaffold families anchor the training set and the test set holds
    the tail. It is seed-independent and is what the literature means by "scaffold
    split", which makes results comparable with published benchmarks.

    scheme="random" assigns whole scaffold groups at random until `test_frac` is
    reached. It is the variant used for repeated-split variance estimates. Unlike
    the earlier implementation it does NOT skip oversized groups, because skipping
    them pushed every large scaffold family into train and silently enriched the
    test set in rare scaffolds; it simply stops once the quota is reached.

    Returns (train_idx, test_idx).
    """
    uniq, counts = np.unique(groups, return_counts=True)
    target_n = test_frac * len(groups)

    if scheme == "deterministic":
        # Largest-first into train; the remainder (smallest families) becomes test.
        order = np.argsort(-counts, kind="stable")
        train_groups, n_train = set(), 0
        train_quota = len(groups) - target_n
        for i in order:
            if n_train >= train_quota:
                break
            train_groups.add(uniq[i])
            n_train += counts[i]
        mask_test = ~np.isin(groups, list(train_groups))
    elif scheme == "random":
        rng = np.random.RandomState(seed)
        order = rng.permutation(len(uniq))
        test_groups, n = set(), 0
        for i in order:
            if n >= target_n:
                break
            test_groups.add(uniq[i])
            n += counts[i]
        mask_test = np.isin(groups, list(test_groups))
    else:
        raise ValueError(f"unknown scheme {scheme!r}")

    return np.where(~mask_test)[0], np.where(mask_test)[0]


# ----------------------------------------------------------------------------
# Fingerprints, Tanimoto kNN
# ----------------------------------------------------------------------------

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=N_FP_BITS)


def morgan_bitvect(smiles):
    m = Chem.MolFromSmiles(smiles)
    return None if m is None else _MORGAN.GetFingerprint(m)


def morgan_bitvects(smiles_list):
    return [morgan_bitvect(s) for s in smiles_list]


def knn_tanimoto(query_fps, train_fps, train_y, k=5, self_indices=None):
    """Similarity-weighted kNN regression on ECFP4 Tanimoto.

    `self_indices[i]` is the position of query i inside `train_fps`, or -1 when the
    query is not in the training set. Passing it excludes each query from its own
    neighbourhood (leave-one-out). The previous API took a bare `exclude_self` flag
    and masked `sims[i]`, which was only correct when query and train lists were
    literally the same object -- silently wrong otherwise.

    Returns (pred, max_sim, mean_topk_sim) arrays over queries.
    """
    train_y = np.asarray(train_y, dtype=float)
    k = int(min(k, max(1, len(train_fps) - 1)))
    pred = np.empty(len(query_fps))
    mx = np.empty(len(query_fps))
    mk = np.empty(len(query_fps))
    for i, q in enumerate(query_fps):
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(q, train_fps))
        if self_indices is not None and self_indices[i] >= 0:
            sims[self_indices[i]] = -1.0
        top = np.argpartition(-sims, k)[:k] if len(sims) > k else np.arange(len(sims))
        w = np.clip(sims[top], 1e-6, None)
        pred[i] = np.sum(w * train_y[top]) / np.sum(w)
        mx[i] = sims[top].max()
        mk[i] = sims[top].mean()
    return pred, mx, mk


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def pearson(a, b):
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
    """Classification view of a regression model: label = y >= thr, score = pred.

    Report these only alongside `prevalence`: when one class dominates, AUC and MCC
    describe ranking within a skewed set, not active-versus-inactive discrimination.
    """
    y_bin = (np.asarray(y) >= thr).astype(int)
    p = np.asarray(p, float)
    out = {"prevalence": float(y_bin.mean()), "n": int(len(y_bin))}
    if y_bin.min() == y_bin.max():
        out.update({"AUC": float("nan"), "PR_AUC": float("nan"), "MCC": float("nan")})
        return out
    out["AUC"] = float(roc_auc_score(y_bin, p))
    out["PR_AUC"] = float(average_precision_score(y_bin, p))
    out["MCC"] = float(matthews_corrcoef(y_bin, (p >= thr).astype(int)))
    return out


def bootstrap_ci(values_fn, n, seed=0, n_boot=N_BOOTSTRAP, alpha=0.05):
    """Percentile bootstrap CI for any statistic over n paired observations.

    `values_fn(idx)` receives a resampled index array and returns a float (or NaN).
    Returns {"lo", "hi", "n_boot"}; NaN draws are dropped.
    """
    rng = np.random.RandomState(seed)
    draws = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        try:
            draws[b] = values_fn(idx)
        except Exception:
            draws[b] = np.nan
    draws = draws[np.isfinite(draws)]
    if draws.size == 0:
        return {"lo": None, "hi": None, "n_boot": 0}
    return {"lo": float(np.quantile(draws, alpha / 2)),
            "hi": float(np.quantile(draws, 1 - alpha / 2)),
            "n_boot": int(draws.size)}


def r2_ci(y, p, seed=0, n_boot=N_BOOTSTRAP):
    y, p = np.asarray(y, float), np.asarray(p, float)

    def stat(idx):
        if np.unique(y[idx]).size < 2:
            return np.nan
        return r2_score(y[idx], p[idx])

    return bootstrap_ci(stat, len(y), seed=seed, n_boot=n_boot)


def auc_ci(y_bin, score, seed=0, n_boot=N_BOOTSTRAP):
    y_bin, score = np.asarray(y_bin, int), np.asarray(score, float)

    def stat(idx):
        if np.unique(y_bin[idx]).size < 2:
            return np.nan
        return roc_auc_score(y_bin[idx], score[idx])

    return bootstrap_ci(stat, len(y_bin), seed=seed, n_boot=n_boot)


def coverage_ci(covered, seed=0, n_boot=N_BOOTSTRAP):
    """Binomial (bootstrap) CI on empirical conformal coverage."""
    covered = np.asarray(covered, bool)
    return bootstrap_ci(lambda idx: covered[idx].mean(), len(covered),
                        seed=seed, n_boot=n_boot)


# ----------------------------------------------------------------------------
# Virtual-screening enrichment metrics
# ----------------------------------------------------------------------------

def enrichment_factor(y_bin, score, frac=0.01):
    """EF at the top `frac` of a ranking (higher score = predicted more active)."""
    y_bin = np.asarray(y_bin, int)
    score = np.asarray(score, float)
    n = len(y_bin)
    n_top = max(1, int(round(frac * n)))
    order = np.argsort(-score, kind="stable")
    hit_rate_top = y_bin[order[:n_top]].mean()
    base = y_bin.mean()
    return float(hit_rate_top / base) if base > 0 else float("nan")


def bedroc(y_bin, score, alpha=20.0):
    """BEDROC (Truchon & Bayly 2007): early-recognition-weighted ranking metric."""
    y_bin = np.asarray(y_bin, int)
    score = np.asarray(score, float)
    n = len(y_bin)
    n_act = int(y_bin.sum())
    if n_act == 0 or n_act == n:
        return float("nan")
    order = np.argsort(-score, kind="stable")
    ranks = np.where(y_bin[order] == 1)[0] + 1
    ra = n_act / n
    rie_num = np.sum(np.exp(-alpha * ranks / n))
    rie_den = (ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1))
    rie = rie_num / rie_den
    rie_max = ((1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha))))
    rie_min = ((1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha))))
    return float((rie - rie_min) / (rie_max - rie_min))


def permutation_null(y_bin, score, metric, n_perm=2000, seed=0):
    """What this metric gives on the same design when the ranking carries no signal.

    EF and BEDROC were built for a large library with a small active fraction.
    On a balanced control set they do not mean what their usual thresholds
    suggest: with 50 actives in 99 compounds, a random ranking returns BEDROC
    (alpha=20) around 0.5 with a 5-95% spread of roughly 0.25 to 0.75, so a
    BEDROC of 0.72 is inside the noise while reading like strong early
    recognition. EF at 1% is worse: the top 1% of 99 compounds is one compound,
    so the statistic has two reachable values.

    Rather than drop the metrics or quote them bare, each is reported against its
    own null, obtained by shuffling the labels. `exceeds_null` is the honest
    summary: whether the observed value falls outside the middle 90% of what
    chance produces on this design.
    """
    y_bin = np.asarray(y_bin, int)
    score = np.asarray(score, float)
    observed = float(metric(y_bin, score))
    rng = np.random.default_rng(seed)
    draws = np.empty(n_perm, float)
    for i in range(n_perm):
        draws[i] = metric(rng.permutation(y_bin), score)
    draws = draws[np.isfinite(draws)]
    if draws.size == 0 or not np.isfinite(observed):
        return None
    lo, hi = (float(np.percentile(draws, 5)), float(np.percentile(draws, 95)))
    return {
        "observed": round(observed, 4),
        "null_mean": round(float(draws.mean()), 4),
        "null_5pct": round(lo, 4),
        "null_95pct": round(hi, 4),
        "p_one_sided": round(float((draws >= observed).mean()), 4),
        "exceeds_null": bool(observed > hi),
        "n_perm": int(draws.size),
    }


# ----------------------------------------------------------------------------
# NP-likeness (Ertl & Schuffenhauer, RDKit Contrib)
# ----------------------------------------------------------------------------

_NP_MODEL = None


def np_likeness(smiles_list):
    """Natural-product-likeness score per SMILES (higher = more NP-like; NaN if bad)."""
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


def load_json(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path) as f:
        return json.load(f)


# ----------------------------------------------------------------------------
# Chemotype and promiscuity annotation
# ----------------------------------------------------------------------------
# Substructure definitions used for chemotype assignment and for flagging
# frequent-hitter kinase chemistry. Assigning chemotypes by SMARTS replaces the
# earlier practice of splitting candidates on an NP-likeness threshold, which is a
# continuous drug-likeness score and carries no structural meaning.

CHEMOTYPE_SMARTS = {
    # Staurosporine / K252a family: fused indolo[2,3-a]carbazole core.
    "indolocarbazole": "c1ccc2c(c1)[nX3]c1c2c2c(c3c1cccc3)[nX3]cc2",
    # Resorcylic acid lactones: resorcinol fused to a macrolactone.
    "resorcylic_acid_lactone": "Oc1cc(O)c(C(=O)O[#6])c(c1)",
    "macrolactone": "[O;R][CX3;R](=O)",
    "flavonoid": "O=C1CC(c2ccccc2)Oc2ccccc21",
    "quinone": "O=C1C=CC(=O)C=C1",
    "catechol": "c1cc(O)c(O)cc1",
}

# Chemotypes that hit many kinases regardless of target. A hit matching one of
# these is a promiscuity liability, not a selective lead, and the paper should say
# so rather than presenting it as a discovery.
PAN_KINASE_ALERTS = ["indolocarbazole", "quinone", "catechol"]


def chemotype_flags(smiles_list):
    """Returns a DataFrame of boolean chemotype columns, one row per input SMILES."""
    patterns = {}
    for name, sma in CHEMOTYPE_SMARTS.items():
        patt = Chem.MolFromSmarts(sma)
        if patt is not None:
            patterns[name] = patt
    rows = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            rows.append({name: False for name in patterns})
            continue
        rows.append({name: m.HasSubstructMatch(p) for name, p in patterns.items()})
    out = pd.DataFrame(rows)
    present = [c for c in PAN_KINASE_ALERTS if c in out.columns]
    out["pan_kinase_alert"] = out[present].any(axis=1) if present else False
    return out


def generic_scaffold(smiles):
    """Generic (element- and bond-order-agnostic) Murcko scaffold SMILES, or None."""
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        return Chem.MolToSmiles(
            MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))
    except Exception:
        return None
