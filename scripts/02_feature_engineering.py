"""
Step 02 - Feature engineering
-----------------------------
Builds the ML-ready table for each target: 21 physicochemical/topological
descriptors + a 2048-bit ECFP4 (Morgan radius 2) fingerprint, carried alongside
the curation metadata from step 01 (pIC50, relation, censored, replicate stats).

Descriptor and fingerprint definitions live in `common.py` so that this script and
the streaming screen in step 05 cannot drift apart.

Two paths were deleted from the earlier version, deliberately:

- The optional `docking_csv` argument, which appended a `DockingScore` feature.
  Docking scores are the orthogonal validation in step 06; feeding them into the
  QSAR features would make that validation circular, and the column silently
  changed the feature count that step 05 asserts on.
- The full phytochemical feature table (~3 GB on disk). Step 05 streams the
  library and featurises in-process, so that file was written and never read.

Output: data/features/<TARGET>_training_features.csv

References
  Rogers & Hahn (2010) Extended-Connectivity Fingerprints. J. Chem. Inf. Model.
  Bemis & Murcko (1996) The properties of known drugs. J. Med. Chem.
"""

import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C

# Metadata carried through from step 01 so downstream code can filter on label
# quality without re-reading the ligand CSV.
CARRY_COLS = ["molecule_chembl_id", "inchikey", "pIC50", "relation", "censored",
              "label_kind", "n_replicates", "replicate_sd", "first_year"]


def featurise_smiles(smiles):
    """One molecule -> {descriptor: value, FP_i: bit} or None if unparsable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        features = C.physchem_descriptors(mol)
    except Exception:
        return None
    bits = C._MORGAN.GetFingerprintAsNumPy(mol)
    for i, bit in enumerate(bits):
        features[f"FP_{i}"] = int(bit)
    return features


def build_feature_table(ligand_csv, output_path):
    if not os.path.exists(ligand_csv):
        print(f"[!] Input not found: {ligand_csv}")
        return None

    df_mols = pd.read_csv(ligand_csv)
    smiles_col = "canonical_smiles" if "canonical_smiles" in df_mols.columns else "smiles"

    print(f"\n[+] Featurising {len(df_mols)} compounds from {ligand_csv}")
    print(f"    {len(C.DESC_COLS)} descriptors + {C.N_FP_BITS}-bit ECFP4")

    records, skipped = [], 0
    for _, row in tqdm(df_mols.iterrows(), total=len(df_mols), desc="features"):
        smi = row.get(smiles_col)
        feat = featurise_smiles(smi) if isinstance(smi, str) else None
        if feat is None:
            skipped += 1
            continue
        feat["canonical_smiles"] = smi
        feat["ligand_id"] = str(row.get("molecule_chembl_id", row.get("ligand_id", ""))
                                ).replace("/", "_").replace("\\", "_")
        for col in CARRY_COLS:
            if col in df_mols.columns:
                feat[col] = row[col]
        if "censored" not in df_mols.columns:
            # Legacy CSV from before censoring was tracked: every row is an exact
            # measurement by construction, because censored relations were not
            # requested. Recorded explicitly so no downstream step has to guess.
            feat["censored"] = False
            feat["label_kind"] = "exact"
        records.append(feat)

    out = pd.DataFrame(records)
    # Stable, explicit column order: metadata, then descriptors, then fingerprint.
    meta = [c for c in ["ligand_id", "canonical_smiles", "molecule_chembl_id",
                        "inchikey", "pIC50", "relation", "censored", "label_kind",
                        "n_replicates", "replicate_sd", "first_year"]
            if c in out.columns]
    out = out[meta + C.DESC_COLS + C.FP_COLS]

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    out.to_csv(output_path, index=False)
    n_exact = int((~out["censored"].astype(bool)).sum()) if "censored" in out else len(out)
    print(f"[OK] {len(out)} molecules x {len(C.FEATURE_COLS)} features -> {output_path}")
    print(f"     usable as regression targets (uncensored): {n_exact}")
    if skipped:
        print(f"     skipped {skipped} unparsable structures")
    return out


def main(targets=None):
    print("=== Step 02: Molecular feature engineering ===")
    for t in (targets or C.TARGETS):
        build_feature_table(
            os.path.join(C.LIGAND_DIR, f"{t}_IC50_training_data.csv"),
            os.path.join(C.FEATURES_DIR, f"{t}_training_features.csv"))
    return {"status": "feature engineering complete"}


if __name__ == "__main__":
    main(sys.argv[1:] or None)
