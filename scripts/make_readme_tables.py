"""Print the README's two summary tables from results/metrics/.

The README states that no number is transcribed into it by hand. These tables are
the exception that keeps that true: they are generated here and checked by
tests/test_readme_tables.py, which fails if the README drifts from the metrics.

Usage: python scripts/make_readme_tables.py            # print both blocks
       python scripts/make_readme_tables.py --check     # compare against README
"""
import argparse
import io
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
METRICS = REPO_ROOT / "results" / "metrics"
README = REPO_ROOT / "README.md"
TARGETS = ("EGFR", "FGFR4", "VEGFR2")
BUDGETS = (1.0, 1.25, 1.5, 1.65, 2.0)
USED_BUDGET = 1.65
NULL_LABELS = {"ROC_AUC": "ROC-AUC", "EF1pct": "EF1%", "EF5pct": "EF5%", "BEDROC_alpha20": "BEDROC"}


def load(name):
    return json.loads((METRICS / name).read_text())


def evidence_table():
    screen = load("screen_summary.json")["targets"]
    dock = load("docking_validation.json")
    rows = ["| Target | AD threshold | Library in domain | Interval-active hits | Redocking top-1 RMSD | Enrichment ROC-AUC (95% CI) | Clears its null |",
            "|---|---|---|---|---|---|---|"]
    for t in TARGETS:
        s = screen[t]
        tau = f"`tau` = {s['tau']:.1f}" + ("" if s["tau_is_binding"] else " (non-binding)")
        redock = dock["redocking_control"][t]
        rmsd = f"{redock['rmsd_top1_A']:.3f} Å " + ("pass" if redock["pose_prediction_success_top1"] else "**fail**")
        enrich = dock["enrichment_control"][t]
        auc = (f"{enrich['docking_ROC_AUC']:.3f} "
               f"[{enrich['docking_ROC_AUC_CI95']['lo']:.3f}, {enrich['docking_ROC_AUC_CI95']['hi']:.3f}]")
        cleared = [NULL_LABELS.get(k, k) for k, v in enrich["null_baselines"].items() if v.get("exceeds_null")]
        rows.append(f"| {t} | {tau} | {100 * s['in_AD_frac']:.2f}% | {s['interval_lower_bound_active']:,} "
                    f"({s['interval_lower_bound_active_novel_scaffold']:,} on novel scaffolds) | {rmsd} | {auc} "
                    f"| {', '.join(cleared) if cleared else '**none**'} |")
    return "\n".join(rows)


def tau_sweep_table():
    sweeps = {t: {row["max_half_width"]: row for row in load(f"final_model_{t}.json")["tau_by_budget"]}
              for t in TARGETS}
    rows = ["| Width budget (pIC50) | " + " | ".join(TARGETS) + " |", "|---|" + "---|" * len(TARGETS)]
    for budget in BUDGETS:
        cells = []
        for t in TARGETS:
            row = sweeps[t].get(budget)
            if row is None or row.get("tau") is None:
                cells.append("no fit")
            elif row.get("binding"):
                cells.append(f"{row['tau']:.1f} binding")
            else:
                cells.append("non-binding")
        label = f"**{budget} (used)**" if budget == USED_BUDGET else f"{budget}"
        rows.append(f"| {label} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


BLOCKS = {"evidence-table": evidence_table, "tau-sweep": tau_sweep_table}


def rendered(name):
    return f"<!-- BEGIN:{name} -->\n{BLOCKS[name]()}\n<!-- END:{name} -->"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="compare the README against the metrics")
    args = ap.parse_args()
    if not args.check:
        # The tables carry Å, which a cp1252 console would refuse.
        out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", newline="\n")
        out.write("\n\n".join(rendered(name) for name in BLOCKS) + "\n")
        out.flush()
        return 0
    text = README.read_text(encoding="utf-8")
    stale = []
    for name in BLOCKS:
        found = re.search(rf"<!-- BEGIN:{name} -->.*?<!-- END:{name} -->", text, re.S)
        if found is None:
            stale.append(f"{name}: block missing from README")
        elif found.group(0).replace("\r\n", "\n").strip() != rendered(name).strip():
            stale.append(f"{name}: README does not match the metrics")
    for line in stale:
        print(line, file=sys.stderr)
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
