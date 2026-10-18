"""
Pipeline orchestrator with checkpoint/resume
=============================================
Runs steps 01-07 in order, saving a checkpoint after each so an interrupted run
resumes where it stopped.

Cross-platform: steps are imported and called in-process rather than shelled out
to a hard-coded interpreter path, and external binaries are resolved through
`common.resolve_binary` (PATH, or VINA_BIN / OBABEL_BIN). Nothing here assumes a
particular operating system or conda prefix.

Usage:
    python scripts/run_pipeline.py                # run (resume if interrupted)
    python scripts/run_pipeline.py --reset        # clear checkpoints, start fresh
    python scripts/run_pipeline.py --from-step 3  # force start at step 3
    python scripts/run_pipeline.py --only 3 4     # run only these steps
    python scripts/run_pipeline.py --list         # show steps and prerequisites
"""

import argparse
import datetime
import importlib.util
import json
import os
import sys
import time
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)
import common as C

CHECKPOINT_FILE = os.path.join(REPO_ROOT, "pipeline_checkpoint.json")

STEPS = [
    (1, "01_fetch_chembl_data", "Download and curate ChEMBL IC50 data",
     ["network access to ChEMBL"]),
    (2, "02_feature_engineering", "Physicochemical descriptors + ECFP4", []),
    (3, "03_scaffold_benchmark", "Scaffold benchmark and model selection", []),
    (4, "04_conformal_ad_models", "Final model, conformal calibration, AD", []),
    (5, "05_screen_coconut_library", "Full natural-product library screen",
     ["data/ligands_phytochemical/phytochemical_screening_set.csv"]),
    (6, "06_docking_validation", "AutoDock Vina structure-based validation",
     ["vina on PATH or VINA_BIN", "obabel on PATH or OBABEL_BIN"]),
    (7, "07_generate_report", "Report and figures", []),
]


def load_module(name):
    path = os.path.join(SCRIPT_DIR, f"{name}.py")
    spec = importlib.util.spec_from_file_location(f"step_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Checkpoint:
    def __init__(self, path=CHECKPOINT_FILE):
        self.path = path
        self.state = {"completed": [], "history": []}
        if os.path.exists(path):
            try:
                with open(path) as f:
                    self.state = json.load(f)
            except Exception:
                pass

    def done(self, step):
        return step in self.state.get("completed", [])

    def mark(self, step, name, secs, result):
        if step not in self.state["completed"]:
            self.state["completed"].append(step)
        self.state["history"].append({
            "step": step, "name": name, "seconds": round(secs, 1),
            "finished_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "result": result if isinstance(result, (dict, str, type(None))) else str(result),
        })
        self.save()

    def reset(self):
        self.state = {"completed": [], "history": []}
        self.save()

    def save(self):
        with open(self.path, "w") as f:
            json.dump(self.state, f, indent=2, default=str)


def check_prereqs(step):
    """Return a list of unmet prerequisites, so a step fails early and clearly."""
    unmet = []
    for _, _, _, prereqs in [s for s in STEPS if s[0] == step]:
        for p in prereqs:
            if p.startswith("data/") or p.startswith("results/"):
                if not os.path.exists(os.path.join(REPO_ROOT, p)):
                    unmet.append(f"missing file: {p}")
            elif " on PATH" in p:
                binary = p.split()[0]
                env_var = p.split("or ")[-1] if "or " in p else None
                try:
                    C.resolve_binary(binary, env_var=env_var, required=True)
                except Exception:
                    unmet.append(f"missing executable: {p}")
    return unmet


def run(from_step=1, only=None, reset=False):
    os.chdir(REPO_ROOT)
    cp = Checkpoint()
    if reset:
        cp.reset()
        print("[+] Checkpoints cleared.")

    selected = set(only) if only else {s[0] for s in STEPS if s[0] >= from_step}
    failures = []

    for step, module_name, title, _ in STEPS:
        if step not in selected:
            continue
        if cp.done(step) and not only:
            print(f"[=] Step {step} ({title}) already complete; skipping.")
            continue

        unmet = check_prereqs(step)
        if unmet:
            print(f"[!] Step {step} ({title}) cannot run:")
            for u in unmet:
                print(f"      - {u}")
            print("      Skipping. Downstream steps that need it will also be skipped.")
            failures.append((step, title, "; ".join(unmet)))
            continue

        print(f"\n{'=' * 70}\n[>] Step {step}: {title}\n{'=' * 70}", flush=True)
        t0 = time.time()
        try:
            module = load_module(module_name)
            result = module.main()
            cp.mark(step, title, time.time() - t0, result)
            print(f"[OK] Step {step} finished in {(time.time() - t0) / 60:.1f} min")
        except Exception as exc:
            traceback.print_exc()
            failures.append((step, title, str(exc)[:200]))
            print(f"[!] Step {step} failed: {exc}")
            break

    print("\n" + "=" * 70)
    print("Completed steps:", sorted(cp.state.get("completed", [])))
    if failures:
        print("Not completed:")
        for step, title, reason in failures:
            print(f"  step {step} ({title}): {reason}")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reset", action="store_true", help="clear checkpoints")
    parser.add_argument("--from-step", type=int, default=1)
    parser.add_argument("--only", type=int, nargs="+")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    if args.list:
        for step, _, title, prereqs in STEPS:
            extra = f"  (needs: {', '.join(prereqs)})" if prereqs else ""
            print(f"  {step}. {title}{extra}")
        return 0

    failures = run(from_step=args.from_step, only=args.only, reset=args.reset)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
