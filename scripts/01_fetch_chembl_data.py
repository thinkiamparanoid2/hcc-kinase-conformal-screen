import os
import sys
import time
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Fix Windows terminal encoding for non-ASCII characters
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

OUTPUT_DIR = os.path.join("data", "ligands_chembl")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGETS = {
    "EGFR": "CHEMBL203",
    "FGFR4": "CHEMBL1974",
    "VEGFR2": "CHEMBL279",
}

MAX_RETRIES = 3
RETRY_DELAY = 30  # seconds between retries

def fetch_target_bioactivity(target_name, target_chembl_id, activity_type="IC50"):
    print(f"\n[+] Fetching bioactivity data for {target_name} ({target_chembl_id}) - {activity_type}...")
    activity = new_client.activity
    
    # Filter constraints applied:
    # - target_organism="Homo sapiens": exclude non-human assays
    # - assay_type="B": only include binding assays (not functional or ADMET assays)
    # - standard_relation__in=["=", "<", "<="]: include potent compounds with threshold relations
    res = activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type=activity_type,
        standard_relation__in=["=", "<", "<="],
        target_organism="Homo sapiens",
        assay_type="B"
    ).only([
        "molecule_chembl_id", "canonical_smiles", "standard_value", 
        "standard_units", "pchembl_value", "assay_type", "standard_relation"
    ])
    
    # Retry logic for transient server errors (ChEMBL API can return HTTP 500)
    records = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            records = list(res)
            break
        except Exception as net_err:
            err_msg = str(net_err)[:200]  # Truncate to avoid HTML dump
            print(f"[!] Attempt {attempt}/{MAX_RETRIES} failed: {err_msg}")
            if attempt < MAX_RETRIES:
                print(f"    Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
            else:
                print(f"[!] All {MAX_RETRIES} attempts failed for {target_name}. Skipping.")
                return None

    if not records:
        print(f"[-] No quantitative {activity_type} records found for {target_name}.")
        return None
        
    df = pd.DataFrame.from_records(records)
    df = df.dropna(subset=["canonical_smiles", "pchembl_value"])
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df.dropna(subset=["pchembl_value"])
    
    print(f"    Raw records after filtering: {len(df)}")
    
    # Group by SMILES for deduplication
    # Calculates mean and standard deviation (std) for duplicate entries
    grouped = df.groupby("canonical_smiles").agg(
        molecule_chembl_id=("molecule_chembl_id", "first"),
        pchembl_value_mean=("pchembl_value", "mean"),
        pchembl_value_std=("pchembl_value", "std"),
        standard_value=("standard_value", "first"),
        standard_units=("standard_units", "first"),
        standard_relation=("standard_relation", "first")
    ).reset_index()
    
    # Conflict detection: 
    # For compounds with only 1 entry, std is NaN, fill with 0.0
    grouped["pchembl_value_std"] = grouped["pchembl_value_std"].fillna(0.0)
    
    initial_len = len(grouped)
    # Drop conflicting data where standard deviation > 1.0 log unit
    df = grouped[grouped["pchembl_value_std"] <= 1.0].copy()
    dropped_count = initial_len - len(df)
    if dropped_count > 0:
        print(f"[!] Dropped {dropped_count} records due to conflicting duplicate values (SD > 1.0).")
    
    # Clean up column names after aggregation
    df = df.rename(columns={"pchembl_value_mean": "pchembl_value"})
    df = df.drop(columns=["pchembl_value_std"])
    
    output_path = os.path.join(OUTPUT_DIR, f"{target_name}_{activity_type}_training_data.csv")
    df.to_csv(output_path, index=False)
    print(f"[OK] Saved {len(df)} curated training records to {output_path}")
    return df

if __name__ == "__main__":
    print("=== ChEMBL Bioactivity Training Data Downloader ===")
    print("    Filters: Homo sapiens | Binding assays only | Conflict detection (SD>1.0)")
    for name, cid in TARGETS.items():
        try:
            df = fetch_target_bioactivity(name, cid, activity_type="IC50")
        except Exception as e:
            err_str = str(e)[:200]
            print(f"[!] Error fetching {name}: {err_str}")
    
    print("\n[OK] Data fetching script complete. Inspect 'data/ligands_chembl/' for output CSVs.")
