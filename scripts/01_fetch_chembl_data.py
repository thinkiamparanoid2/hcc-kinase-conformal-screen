"""
Step 01 - ChEMBL IC50 download and curation
-------------------------------------------
Fetches human binding-assay IC50 records for EGFR, FGFR4 and VEGFR2 and curates
them into one CSV per target.

Three curation decisions matter for every claim downstream, so they are explicit:

1. CENSORED RECORDS ARE KEPT AND LABELLED, NOT SILENTLY TREATED AS EXACT.
   ChEMBL stores `standard_relation`. `<` / `<=` mean "IC50 is below this value",
   i.e. the compound is at least this potent -- pIC50 is a LOWER bound. `>` / `>=`
   mean the opposite -- pIC50 is an UPPER bound, and these are the genuinely
   inactive compounds. The earlier version of this script requested only
   `=`, `<`, `<=`, which (a) discarded every true inactive, leaving ~75-80% of
   labels "active" with no negative class, and (b) fed lower-bounded values into
   the regression as if they were point measurements. Both are now fixed: all five
   relations are fetched, each row carries `relation` / `censored` / `label_kind`,
   and only uncensored rows are admissible as regression targets
   (`common.load_training(exact_only=True)`, the default).

2. STRUCTURES ARE STANDARDISED BEFORE DEDUPLICATION.
   Salt-strip -> neutralise -> canonical tautomer, then deduplicate on InChIKey.
   Deduplicating on raw canonical SMILES (the old behaviour) left ~2% of each set
   as salt/free-base or tautomer pairs of the same parent. Those pairs get
   different generic Murcko scaffolds, so the same molecule could appear on both
   sides of a scaffold split.

3. REPLICATE SPREAD IS REPORTED, NOT JUST FILTERED.
   Compounds whose exact replicates disagree by more than `MAX_REPLICATE_SD`
   log units are dropped, as before, but the per-compound replicate count and SD
   are written to the CSV and the distribution is summarised to
   results/metrics/curation_<TARGET>.json. That summary is what bounds the
   achievable R2, so the paper can state the label-noise ceiling instead of
   hand-waving at it.

Output: data/ligands_chembl/<TARGET>_IC50_training_data.csv
        results/metrics/curation_<TARGET>.json
"""

import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C

if hasattr(sys.stdout, "reconfigure"):      # non-ASCII safety on Windows consoles
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MAX_RETRIES = 3
RETRY_DELAY = 30                 # seconds between retries
MAX_REPLICATE_SD = 1.0           # log units; drop compounds noisier than this
RELATIONS = ["=", "<", "<=", ">", ">="]

# relation -> (censored?, what the bound means)
_LABEL_KIND = {
    "=":  (False, "exact"),
    "<":  (True, "lower_bound_potency"),    # at least this potent
    "<=": (True, "lower_bound_potency"),
    ">":  (True, "upper_bound_potency"),    # at most this potent -> inactive
    ">=": (True, "upper_bound_potency"),
}


def pic50_from_nm(value_nm):
    """pIC50 = 9 - log10(IC50 in nM). Computed consistently for every relation.

    ChEMBL's own `pchembl_value` is populated only for uncensored records, so
    deriving the label here is what lets censored rows carry a usable bound.
    """
    v = pd.to_numeric(value_nm, errors="coerce")
    v = v.where(v > 0)
    return 9.0 - np.log10(v)


def fetch_records(target_chembl_id, activity_type="IC50"):
    from chembl_webresource_client.new_client import new_client
    res = new_client.activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type=activity_type,
        standard_relation__in=RELATIONS,
        standard_units="nM",
        target_organism="Homo sapiens",
        assay_type="B",
    ).only([
        "molecule_chembl_id", "canonical_smiles", "standard_value",
        "standard_units", "pchembl_value", "assay_type", "standard_relation",
        "assay_chembl_id", "document_chembl_id", "document_year",
    ])
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return list(res)
        except Exception as err:
            print(f"[!] Attempt {attempt}/{MAX_RETRIES} failed: {str(err)[:200]}")
            if attempt < MAX_RETRIES:
                print(f"    Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
    return None


def _first_year(g):
    """Earliest publication year across a compound's records, for temporal splits.

    A temporal split (train on older literature, test on newer) is the closest
    cheap stand-in for prospective validation, so the year is carried through
    curation even though ChEMBL does not always populate it.
    """
    if "document_year" not in g.columns:
        return np.nan
    years = pd.to_numeric(g["document_year"], errors="coerce").dropna()
    return float(years.min()) if len(years) else np.nan


def curate(target_name, records):
    """Raw activity records -> curated per-parent table + curation summary."""
    df = pd.DataFrame.from_records(records)
    n_raw = len(df)
    if "document_year" not in df.columns:
        df["document_year"] = np.nan
    df = df.dropna(subset=["canonical_smiles", "standard_value"])
    df = df[df["standard_units"] == "nM"]
    df["relation"] = df["standard_relation"].fillna("=")
    df = df[df["relation"].isin(RELATIONS)]
    df["pIC50"] = pic50_from_nm(df["standard_value"])
    df = df.dropna(subset=["pIC50"])

    # --- standardise structures, then deduplicate on the parent InChIKey ---
    uniq = df["canonical_smiles"].drop_duplicates()
    std_smi, std_key = C.standardise_series(uniq)
    smi_map = dict(zip(uniq, std_smi))
    key_map = dict(zip(uniq, std_key))
    df["standardised_smiles"] = df["canonical_smiles"].map(smi_map)
    df["inchikey"] = df["canonical_smiles"].map(key_map)
    n_unparsed = int(df["inchikey"].isna().sum())
    df = df.dropna(subset=["inchikey", "standardised_smiles"])

    df["censored"] = df["relation"].map(lambda r: _LABEL_KIND[r][0])
    df["label_kind"] = df["relation"].map(lambda r: _LABEL_KIND[r][1])

    rows = []
    for key, g in df.groupby("inchikey", sort=False):
        exact = g[~g["censored"]]
        if len(exact):
            # Uncensored measurements win: they are the only point-valued labels.
            sd = float(exact["pIC50"].std(ddof=1)) if len(exact) > 1 else 0.0
            rows.append({
                "inchikey": key,
                "canonical_smiles": exact["standardised_smiles"].iloc[0],
                "original_smiles": exact["canonical_smiles"].iloc[0],
                "molecule_chembl_id": exact["molecule_chembl_id"].iloc[0],
                "pIC50": float(exact["pIC50"].mean()),
                "relation": "=",
                "censored": False,
                "label_kind": "exact",
                "n_replicates": int(len(exact)),
                "replicate_sd": sd,
                "n_censored_records": int(len(g) - len(exact)),
                "first_year": _first_year(exact),
            })
            continue
        # Censored only. Keep the least extreme bound so the row is conservative:
        # the tightest lower bound on potency (min) for `<`, the tightest upper
        # bound (max) for `>`.
        kind = g["label_kind"].iloc[0]
        agg = g["pIC50"].min() if kind == "lower_bound_potency" else g["pIC50"].max()
        rows.append({
            "inchikey": key,
            "canonical_smiles": g["standardised_smiles"].iloc[0],
            "original_smiles": g["canonical_smiles"].iloc[0],
            "molecule_chembl_id": g["molecule_chembl_id"].iloc[0],
            "pIC50": float(agg),
            "relation": g["relation"].iloc[0],
            "censored": True,
            "label_kind": kind,
            "n_replicates": int(len(g)),
            "replicate_sd": np.nan,
            "n_censored_records": int(len(g)),
            "first_year": _first_year(g),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(
            f"{target_name}: no usable records after curation. Check that the "
            f"input carries standard_value in nM (or pIC50) and parsable SMILES.")
    n_parents = len(out)

    noisy = out["censored"].eq(False) & out["replicate_sd"].gt(MAX_REPLICATE_SD)
    n_noisy = int(noisy.sum())
    out = out[~noisy].reset_index(drop=True)

    exact_rows = out[~out["censored"]]
    multi = exact_rows[exact_rows["n_replicates"] > 1]["replicate_sd"]
    # Label-noise ceiling: with replicate SD s, the best achievable RMSE for a
    # perfect model is s, so R2_max = 1 - s^2 / var(y).
    noise_sd = float(multi.median()) if len(multi) else float("nan")
    var_y = float(exact_rows["pIC50"].var(ddof=1)) if len(exact_rows) > 1 else float("nan")
    r2_ceiling = float(1.0 - (noise_sd ** 2) / var_y) if np.isfinite(noise_sd) and var_y > 0 else float("nan")

    summary = {
        "target": target_name,
        "n_raw_records": int(n_raw),
        "n_unparsable_structures": n_unparsed,
        "n_unique_parents_before_noise_filter": int(n_parents),
        "n_dropped_replicate_sd_above_threshold": n_noisy,
        "max_replicate_sd": MAX_REPLICATE_SD,
        "n_final": int(len(out)),
        "n_exact": int((~out["censored"]).sum()),
        "n_censored_lower_bound_potency": int((out["label_kind"] == "lower_bound_potency").sum()),
        "n_censored_upper_bound_potency": int((out["label_kind"] == "upper_bound_potency").sum()),
        "exact_active_prevalence_at_pIC50_6": float((exact_rows["pIC50"] >= C.ACTIVE_THRESHOLD).mean()),
        "replicate_sd_median_multi_measured": noise_sd,
        "replicate_sd_n_multi_measured": int(len(multi)),
        "label_variance_exact": var_y,
        "r2_ceiling_from_label_noise": r2_ceiling,
    }
    return out, summary


def recurate_existing(target_name):
    """Re-curate an already-downloaded CSV without touching the network.

    Applies decisions 2 and 3 (standardisation, parent-level deduplication,
    replicate reporting) to a file produced by an earlier run. Decision 1 cannot be
    applied retroactively: a CSV fetched before censored relations were requested
    simply has no `>` records in it, so the curation summary records
    `censored_records_available: false` and the paper must state that the inactive
    class is absent until step 01 is re-run with network access.
    """
    path = os.path.join(C.LIGAND_DIR, f"{target_name}_IC50_training_data.csv")
    if not os.path.exists(path):
        return None, None
    raw = pd.read_csv(path)
    smi_col = "original_smiles" if "original_smiles" in raw.columns else "canonical_smiles"
    # A file already written by this script carries pIC50 but not the raw
    # standard_value it came from. Invert the transform so re-curation is
    # idempotent instead of silently dropping every row.
    if "standard_value" not in raw.columns and "pIC50" in raw.columns:
        raw = raw.copy()
        raw["standard_value"] = 10 ** (9.0 - pd.to_numeric(raw["pIC50"],
                                                           errors="coerce"))
        raw["standard_units"] = "nM"
    records = []
    for _, r in raw.iterrows():
        records.append({
            "molecule_chembl_id": r.get("molecule_chembl_id"),
            "canonical_smiles": r.get(smi_col),
            "standard_value": r.get("standard_value"),
            "standard_units": r.get("standard_units", "nM"),
            "standard_relation": r.get("relation", r.get("standard_relation", "=")),
            "pchembl_value": r.get("pchembl_value"),
            "document_year": r.get("document_year"),
        })
    out, summary = curate(target_name, records)
    summary["source"] = "recurated_from_existing_csv"
    # An already-curated CSV holds one averaged row per compound, so the original
    # replicate spread is gone. Reporting a label-noise ceiling from it would give
    # a spurious R2_max of ~1.0; null it out and say why instead.
    for k in ("replicate_sd_median_multi_measured", "replicate_sd_n_multi_measured",
              "r2_ceiling_from_label_noise"):
        summary[k] = None
    summary["label_noise_note"] = (
        "Replicate spread is not recoverable from a pre-aggregated CSV, so the "
        "label-noise ceiling on R2 is unknown here. Re-run step 01 with network "
        "access to measure it from the raw ChEMBL records.")
    summary["censored_records_available"] = bool(
        raw.get("relation", raw.get("standard_relation", pd.Series(["="]))).astype(str)
        .isin([">", ">="]).any())
    return out, summary


def main(targets=None, from_existing=False):
    os.makedirs(C.LIGAND_DIR, exist_ok=True)
    os.makedirs(C.METRICS_DIR, exist_ok=True)
    print("=== Step 01: ChEMBL IC50 download and curation ===")
    print("    Homo sapiens | binding assays | all relations kept and labelled")
    print("    Structures standardised; duplicates collapsed on parent InChIKey")
    if from_existing:
        print("    MODE: re-curating existing CSVs offline (no censored records added)")
        for name in (targets or C.TARGETS):
            out, summary = recurate_existing(name)
            if out is None:
                print(f"[!] {name}: no existing CSV to re-curate")
                continue
            path = os.path.join(C.LIGAND_DIR, f"{name}_IC50_training_data.csv")
            out.to_csv(path, index=False)
            C.save_json(summary, os.path.join(C.METRICS_DIR, f"curation_{name}.json"))
            print(f"[OK] {name}: {summary['n_raw_records']} rows -> {len(out)} parents "
                  f"(collapsed {summary['n_raw_records'] - len(out)}); "
                  f"R2 ceiling {summary['r2_ceiling_from_label_noise']:.3f}")
        return {"status": "recurated existing csvs"}
    for name in (targets or C.TARGETS):
        cid = C.TARGET_CHEMBL_IDS[name]
        print(f"\n[+] {name} ({cid})...", flush=True)
        records = fetch_records(cid)
        if not records:
            print(f"[!] No records retrieved for {name}; leaving existing file untouched.")
            continue
        out, summary = curate(name, records)
        path = os.path.join(C.LIGAND_DIR, f"{name}_IC50_training_data.csv")
        out.to_csv(path, index=False)
        C.save_json(summary, os.path.join(C.METRICS_DIR, f"curation_{name}.json"))
        print(f"[OK] {len(out)} parents -> {path}")
        print(f"     exact {summary['n_exact']} | "
              f"censored-potent {summary['n_censored_lower_bound_potency']} | "
              f"censored-inactive {summary['n_censored_upper_bound_potency']} | "
              f"dropped noisy {summary['n_dropped_replicate_sd_above_threshold']}")
        print(f"     active prevalence among exact labels: "
              f"{summary['exact_active_prevalence_at_pIC50_6']:.3f} | "
              f"R2 ceiling from label noise: {summary['r2_ceiling_from_label_noise']:.3f}")
    return {"status": "chembl fetch complete"}


if __name__ == "__main__":
    argv = [a for a in sys.argv[1:] if a != "--from-existing"]
    main(argv or None, from_existing="--from-existing" in sys.argv[1:])
