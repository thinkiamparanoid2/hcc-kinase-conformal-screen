"""
Master Pipeline Orchestrator with Checkpoint/Resume Support
============================================================
This script runs the ENTIRE ML pipeline end-to-end with automatic
checkpoint saving after every step. If power goes out (loadshedding)
or the process is interrupted, simply re-run this script and it will
RESUME from the last completed step.

Usage:
    python scripts/run_pipeline.py              # Run full pipeline (resume if interrupted)
    python scripts/run_pipeline.py --reset       # Clear all checkpoints and start fresh
    python scripts/run_pipeline.py --from-step 3 # Force start from step 3

Environment:
    Activate: conda activate hcc-screen
    Repo:     hcc-kinase-conformal-screen
"""

import os
import sys
import json
import time
import datetime
import traceback

_BIN_DIR = os.path.join(sys.prefix, "Library", "bin")
if os.path.isdir(_BIN_DIR):
    os.environ["PATH"] = _BIN_DIR + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_BIN_DIR)

# Fix Windows terminal encoding
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ============================================================================
# Configuration — ALL DECISIONS PRE-MADE (no ambiguity for next agent)
# ============================================================================

CONFIG = {
    # Targets to process (order: smallest first for fast debugging)
    "targets": ["FGFR4", "EGFR", "VEGFR2"],

    # Target ChEMBL IDs
    "target_ids": {
        "EGFR": "CHEMBL203",
        "FGFR4": "CHEMBL1974",
        "VEGFR2": "CHEMBL279",
    },

    # Feature engineering settings
    "morgan_radius": 2,
    "morgan_bits": 2048,

    # ML training settings
    "pIC50_threshold": 6.0,      # Active/inactive cutoff (6.0 = 1 uM IC50)
    "test_fraction": 0.2,        # 80/20 scaffold split
    "random_state": 42,

    # Applicability domain: calibrated empirically in scripts/04_conformal_ad_models.py
    # (Mondrian conformal by max-Tanimoto bin; see constants at top of that file).

    # Paths
    "data_dir": "data",
    "features_dir": os.path.join("data", "features"),
    "models_dir": "models",
    "results_dir": "results",
    "checkpoint_file": "pipeline_checkpoint.json",

    # Python executable
    "python_exe": r"F:\miniconda3\envs\hcc-screen\python.exe",
}


# ============================================================================
# Checkpoint Manager
# ============================================================================

class CheckpointManager:
    """Tracks completed pipeline steps and allows resume after interruption."""

    def __init__(self, checkpoint_path):
        self.path = checkpoint_path
        self.state = self._load()

    def _load(self):
        if os.path.exists(self.path):
            with open(self.path, "r") as f:
                return json.load(f)
        return {"completed_steps": [], "step_results": {}, "started_at": None}

    def save(self):
        with open(self.path, "w") as f:
            json.dump(self.state, f, indent=2, default=str)

    def is_done(self, step_id):
        return step_id in self.state["completed_steps"]

    def mark_done(self, step_id, result=None):
        if step_id not in self.state["completed_steps"]:
            self.state["completed_steps"].append(step_id)
        self.state["step_results"][step_id] = {
            "completed_at": datetime.datetime.now().isoformat(),
            "result": result
        }
        self.save()

    def reset(self):
        self.state = {"completed_steps": [], "step_results": {}, "started_at": None}
        self.save()

    def set_started(self):
        if not self.state["started_at"]:
            self.state["started_at"] = datetime.datetime.now().isoformat()
            self.save()


# ============================================================================
# Directory Setup
# ============================================================================

def setup_directories():
    """Create all output directories upfront."""
    dirs = [
        CONFIG["features_dir"],
        CONFIG["models_dir"],
        os.path.join(CONFIG["results_dir"], "predictions"),
        os.path.join(CONFIG["results_dir"], "reports"),
        os.path.join(CONFIG["results_dir"], "figures"),
        os.path.join(CONFIG["results_dir"], "docking_poses"),
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    print("[OK] All output directories created.")


# ============================================================================
# Pipeline Steps
# ============================================================================

def step_1_verify_training_data():
    """Verify that ChEMBL training CSVs exist and have valid data."""
    import pandas as pd
    results = {}
    for target in CONFIG["targets"]:
        csv_path = os.path.join(CONFIG["data_dir"], "ligands_chembl",
                                f"{target}_IC50_training_data.csv")
        if not os.path.exists(csv_path):
            raise FileNotFoundError(
                f"Training data missing: {csv_path}\n"
                f"Run: python scripts/01_fetch_chembl_data.py"
            )
        df = pd.read_csv(csv_path)
        n = len(df)
        has_smiles = "canonical_smiles" in df.columns
        has_pchembl = "pchembl_value" in df.columns
        if not has_smiles or not has_pchembl:
            raise ValueError(f"{csv_path} missing required columns")
        results[target] = n
        print(f"    {target}: {n} records [OK]")

    # Verify phytochemical screening set
    screen_path = os.path.join(CONFIG["data_dir"], "ligands_phytochemical",
                               "phytochemical_screening_set.csv")
    if not os.path.exists(screen_path):
        raise FileNotFoundError(f"Screening set missing: {screen_path}")
    df_screen = pd.read_csv(screen_path, nrows=10)
    print(f"    Phytochemical screening set: present [OK]")

    # Verify PDB structures
    for pdb in ["egfr_1m17.pdb", "fgfr4_5jkg.pdb", "vegfr2_4asd.pdb"]:
        pdb_path = os.path.join(CONFIG["data_dir"], "targets", pdb)
        if not os.path.exists(pdb_path):
            print(f"    [WARN] PDB file missing: {pdb} (docking step will be skipped)")

    return results


def step_2_feature_engineering_training():
    """Generate feature tables for all training datasets."""
    sys.path.insert(0, "scripts")
    if "scripts.02_feature_engineering" in sys.modules:
        del sys.modules["scripts.02_feature_engineering"]

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "feat_eng", "scripts/02_feature_engineering.py")
    feat_eng = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(feat_eng)

    results = {}
    for target in CONFIG["targets"]:
        bio_csv = os.path.join(CONFIG["data_dir"], "ligands_chembl",
                               f"{target}_IC50_training_data.csv")
        dock_csv = os.path.join(CONFIG["results_dir"], "docking_poses",
                                f"{target}_binding_affinities.csv")
        out_csv = os.path.join(CONFIG["features_dir"],
                               f"{target}_training_features.csv")

        if os.path.exists(out_csv):
            import pandas as pd
            existing = len(pd.read_csv(out_csv))
            print(f"    {target}: Feature file already exists ({existing} rows). Skipping.")
            results[target] = existing
            continue

        df = feat_eng.build_feature_table(
            bio_csv, dock_csv, out_csv,
            is_training=True,
            morgan_radius=CONFIG["morgan_radius"],
            morgan_bits=CONFIG["morgan_bits"]
        )
        results[target] = len(df) if df is not None else 0

    return results


def _load_module(name, path):
    """Load a numerically-prefixed script (e.g. 03_scaffold_benchmark.py) as a module."""
    import importlib.util
    sys.path.insert(0, "scripts")          # so scripts can `import common`
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def step_3_benchmark():
    """Repeated scaffold-split benchmark: baselines + RF + XGBoost + scramble control."""
    return _load_module("bench", "scripts/03_scaffold_benchmark.py").main()


def step_4_final_models_ad():
    """Final XGBoost models, Mondrian conformal intervals, empirical AD, NP-likeness analysis."""
    return _load_module("final", "scripts/04_conformal_ad_models.py").main()


def step_5_screen_library():
    """Screen the FULL natural-product library (all ~738k compounds)."""
    import subprocess
    subprocess.check_call([sys.executable, "-u", "scripts/05_screen_coconut_library.py"])
    with open(os.path.join(CONFIG["results_dir"], "metrics", "screen_summary.json")) as f:
        s = json.load(f)
    return {"status": "full library screened", "n": s.get("library_valid", 0)}


def step_6_docking_validation():
    """Orthogonal 3D structure-based validation with AutoDock Vina 1.2.5."""
    return _load_module("dock_val", "scripts/06_docking_validation.py").main()


def step_7_report():
    """Generate the results report and figures from metrics/prediction files."""
    return _load_module("report", "scripts/07_generate_report.py").main()


# ============================================================================
# Pipeline Runner
# ============================================================================

PIPELINE_STEPS = [
    ("step_1_verify_data",        "Verify all input datasets exist and are valid",
     step_1_verify_training_data),
    ("step_2_features_training",  "Feature engineering for training data (21 descriptors + 2048-bit ECFP4)",
     step_2_feature_engineering_training),
    ("step_3_benchmark",          "Repeated scaffold-split benchmark (baselines, RF, XGBoost, scramble control)",
     step_3_benchmark),
    ("step_4_final_models_ad",    "Final models + conformal intervals + calibrated applicability domain",
     step_4_final_models_ad),
    ("step_5_screen_library",     "Full-library natural-product screen (~738k compounds)",
     step_5_screen_library),
    ("step_6_docking_validation", "Orthogonal structure-based docking validation (AutoDock Vina 1.2.5)",
     step_6_docking_validation),
    ("step_7_report",             "Generate results report and publication figures",
     step_7_report),
]


def run_pipeline(from_step=None, reset=False):
    """Run the full pipeline with checkpoint/resume support."""
    ckpt = CheckpointManager(CONFIG["checkpoint_file"])

    if reset:
        ckpt.reset()
        print("[OK] All checkpoints cleared. Starting fresh.\n")

    ckpt.set_started()
    setup_directories()

    print("\n" + "=" * 70)
    print("  HCC PHYTOMEDICINE ML PIPELINE — MASTER ORCHESTRATOR")
    print("  Checkpoint/Resume enabled. Safe to interrupt and restart.")
    print("=" * 70)

    total = len(PIPELINE_STEPS)
    for i, (step_id, description, func) in enumerate(PIPELINE_STEPS, 1):
        # Skip steps if from_step is specified
        if from_step and i < from_step:
            continue

        # Skip already completed steps
        if ckpt.is_done(step_id) and not from_step:
            print(f"\n[SKIP] Step {i}/{total}: {description}")
            print(f"       (Already completed. Delete pipeline_checkpoint.json to rerun)")
            continue

        print(f"\n{'='*70}")
        print(f"  STEP {i}/{total}: {description}")
        print(f"  Step ID: {step_id}")
        print(f"  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*70}\n")

        try:
            result = func()
            ckpt.mark_done(step_id, result)
            print(f"\n[OK] Step {i}/{total} complete: {step_id}")
        except Exception as e:
            print(f"\n[FAILED] Step {i}/{total} failed: {step_id}")
            print(f"         Error: {str(e)[:500]}")
            traceback.print_exc()
            print(f"\n[!] Pipeline paused. Fix the issue and re-run this script.")
            print(f"    It will automatically resume from step {i}.")
            sys.exit(1)

    print(f"\n{'='*70}")
    print(f"  PIPELINE COMPLETE!")
    print(f"  All {total} steps finished successfully.")
    print(f"  Results saved to: {CONFIG['results_dir']}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    args = sys.argv[1:]

    if "--reset" in args:
        run_pipeline(reset=True)
    elif "--from-step" in args:
        idx = args.index("--from-step")
        step_num = int(args[idx + 1])
        run_pipeline(from_step=step_num)
    else:
        run_pipeline()
