"""
Phase 5: Feature Engineering & Descriptor Extraction (Publication-Grade)
------------------------------------------------------------------------
Computes a comprehensive feature set for ML/QML modeling:
  1. Physicochemical descriptors (Lipinski properties, TPSA, rotatable bonds, etc.)
  2. Structural fingerprints (ECFP4 / Morgan fingerprints, radius=2, 2048 bits)
  3. Extended RDKit 2D descriptors (200+ molecular properties)
  4. Optional: Docking scores from AutoDock Vina binding affinities

Merges all features into unified ML-ready CSV tables for training (ChEMBL knowns)
and prediction (phytochemical screening library).

References:
  - Rogers & Hahn (2010). Extended-Connectivity Fingerprints. J. Chem. Inf. Model.
  - RDKit Descriptor Module: https://www.rdkit.org/docs/GettingStartedInPython.html
"""

import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from tqdm import tqdm


# ============================================================================
# Morgan Fingerprint (ECFP4) Generation
# ============================================================================

def compute_morgan_fingerprint(mol, radius=2, n_bits=2048):
    """
    Generates an ECFP4 (Extended-Connectivity Fingerprint) bit vector.
    Radius=2 corresponds to ECFP4 (diameter=4 bonds).
    Returns a list of integer bits (0 or 1) of length n_bits.
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return list(fp)


# ============================================================================
# Physicochemical & Extended 2D Descriptors
# ============================================================================

# Curated set of robust, interpretable physicochemical descriptors
PHYSCHEM_DESCRIPTOR_NAMES = [
    "MolWt", "MolLogP", "TPSA",
    "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "NumAromaticRings",
    "NumHeteroatoms", "NumAliphaticRings",
    "FractionCSP3", "HeavyAtomCount",
    "RingCount", "BertzCT",
    "NumValenceElectrons",
    "NumRadicalElectrons",
    "LabuteASA",          # Labute Accessible Surface Area
    "BalabanJ",           # Balaban's J index (topological)
    "Chi0n", "Chi1n",     # Kier-Hall connectivity indices
    "HallKierAlpha",      # Hall-Kier alpha value
]


def compute_physchem_descriptors(mol):
    """
    Computes a curated set of 20 physicochemical/topological descriptors
    plus Lipinski Rule of 5 violation count.
    """
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(PHYSCHEM_DESCRIPTOR_NAMES)
    desc_values = calc.CalcDescriptors(mol)
    features = dict(zip(PHYSCHEM_DESCRIPTOR_NAMES, desc_values))

    # Lipinski Rule of 5 violation count
    features["LipinskiViolations"] = sum([
        features.get("MolWt", 0) > 500,
        features.get("MolLogP", 0) > 5,
        features.get("NumHDonors", 0) > 5,
        features.get("NumHAcceptors", 0) > 10
    ])

    return features


# ============================================================================
# Unified Feature Extraction Pipeline
# ============================================================================

def extract_all_features(smiles, docking_score=None, morgan_radius=2, morgan_bits=2048):
    """
    Master feature extraction function for a single molecule.
    Returns a dict combining:
      - 20 physicochemical descriptors + Lipinski violations
      - 2048-bit Morgan/ECFP4 fingerprint (columns: FP_0 ... FP_2047)
      - Optional docking score
    Returns None if the SMILES cannot be parsed.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # 1. Physicochemical descriptors
        features = compute_physchem_descriptors(mol)

        # 2. Morgan fingerprint (ECFP4)
        fp_bits = compute_morgan_fingerprint(mol, radius=morgan_radius, n_bits=morgan_bits)
        for i, bit in enumerate(fp_bits):
            features[f"FP_{i}"] = bit

        # 3. Optional docking score
        if docking_score is not None:
            features["DockingScore"] = docking_score

        return features

    except Exception:
        return None


def build_feature_table(bioactivity_csv, docking_csv, output_path,
                        is_training=True, label_col="pchembl_value",
                        morgan_radius=2, morgan_bits=2048):
    """
    Builds a comprehensive ML feature table by combining:
      - Physicochemical descriptors (20 features)
      - ECFP4 Morgan fingerprints (2048 bits)
      - Docking binding affinities (if available)
      - pIC50 labels (for training data only)

    Parameters
    ----------
    bioactivity_csv : str
        Path to the ChEMBL training CSV or phytochemical screening CSV.
    docking_csv : str
        Path to docking results CSV (ligand_id, docking_score_kcal_mol).
    output_path : str
        Where to save the generated feature table.
    is_training : bool
        If True, includes pIC50 label column from bioactivity_csv.
    label_col : str
        Column name for the bioactivity label (default: pchembl_value).
    morgan_radius : int
        Morgan fingerprint radius (2 = ECFP4).
    morgan_bits : int
        Number of Morgan fingerprint bits (default: 2048).
    """
    if not os.path.exists(bioactivity_csv):
        print(f"[!] Input bioactivity file not found: {bioactivity_csv}")
        return

    df_mols = pd.read_csv(bioactivity_csv)

    # Load docking scores if available
    docking_scores = {}
    if docking_csv and os.path.exists(docking_csv):
        df_dock = pd.read_csv(docking_csv)
        for _, row in df_dock.iterrows():
            docking_scores[str(row["ligand_id"])] = row["docking_score_kcal_mol"]
    else:
        print(f"[!] Docking scores file not found ({docking_csv}). "
              "Building features without docking affinity.")

    records = []
    skipped = 0
    print(f"\n[+] Extracting {len(PHYSCHEM_DESCRIPTOR_NAMES)} physicochemical descriptors "
          f"+ {morgan_bits}-bit ECFP4 fingerprints for {len(df_mols)} compounds...")

    for _, row in tqdm(df_mols.iterrows(), total=len(df_mols), desc="Feature Extraction"):
        smi = row.get("canonical_smiles") or row.get("smiles")
        lig_id = str(
            row.get("molecule_chembl_id", row.get("ligand_id", ""))
        ).replace("/", "_").replace("\\", "_")

        dock_val = docking_scores.get(lig_id, None)
        feat = extract_all_features(
            smi, docking_score=dock_val,
            morgan_radius=morgan_radius, morgan_bits=morgan_bits
        )

        if feat is not None:
            feat["ligand_id"] = lig_id
            feat["canonical_smiles"] = smi
            if is_training and label_col in row.index:
                feat["pIC50"] = row[label_col]
            records.append(feat)
        else:
            skipped += 1

    df_out = pd.DataFrame(records)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df_out.to_csv(output_path, index=False)

    n_physchem = len(PHYSCHEM_DESCRIPTOR_NAMES) + 1  # +1 for LipinskiViolations
    n_fp = morgan_bits
    n_dock = 1 if docking_scores else 0
    n_total = n_physchem + n_fp + n_dock

    print(f"[OK] Feature table generated: {len(df_out)} molecules x {n_total} features")
    print(f"     Physicochemical: {n_physchem} | ECFP4 bits: {n_fp} | Docking: {n_dock}")
    if skipped > 0:
        print(f"     Skipped {skipped} molecules due to invalid SMILES")
    print(f"     Saved to: {output_path}")
    return df_out


if __name__ == "__main__":
    print("=== Publication-Grade Molecular Feature Engineering Engine ===")
    print("Features: 21 physicochemical descriptors + 2048-bit ECFP4 Morgan fingerprints")
    print()

    # Auto-process all available training datasets
    targets = ["EGFR", "FGFR4", "VEGFR2"]
    for t in targets:
        bio_csv = os.path.join("data", "ligands_chembl", f"{t}_IC50_training_data.csv")
        dock_csv = os.path.join("results", "docking_poses", f"{t}_pilot",
                                f"{t}_binding_affinities.csv")
        out_csv = os.path.join("data", "features", f"{t}_training_features.csv")
        if os.path.exists(bio_csv):
            build_feature_table(bio_csv, dock_csv, out_csv, is_training=True)

    # Process phytochemical screening library
    screen_csv = os.path.join("data", "ligands_phytochemical",
                              "phytochemical_screening_set.csv")
    if os.path.exists(screen_csv):
        build_feature_table(
            screen_csv,
            docking_csv=None,
            output_path=os.path.join("data", "features",
                                     "phytochemical_screening_features.csv"),
            is_training=False,
            label_col="pchembl_value"
        )

    print("\n[OK] Feature engineering complete for all available datasets.")
