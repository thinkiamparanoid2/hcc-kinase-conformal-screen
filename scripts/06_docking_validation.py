"""
Step 06 - Orthogonal Structure-Based Validation with AutoDock Vina 1.2.5
------------------------------------------------------------------------
1. Prepares clean receptor PDBQTs (Chain A, polar H at pH 7.4, Gasteiger charges)
   from data/targets/{egfr_1m17,fgfr4_5jkg,vegfr2_4asd}.pdb and extracts exact
   binding-site centroids from co-crystallized inhibitors (AQ4, 6LF, BAX).
2. Control 1 (Redocking): Redocks each native co-crystal ligand into its own
   pocket and computes heavy-atom RMSD (top-1 pose and best-of-9 poses).
3. Control 2 (Active vs Weak/Inactive Enrichment): Docks MW-matched potent
   ChEMBL actives (pIC50 >= 8.0) vs weak/inactives (pIC50 <= 5.0) per target
   and computes Docking ROC-AUC and Mann-Whitney U p-value.
4. Candidate Docking: Docks the top 25 novel consensus natural products (from
   consensus_hits_annotated.csv) and top 8 rediscovery positive controls across
   all 3 HCC kinase pockets (EGFR, FGFR4, VEGFR2).

Outputs:
  results/docking_poses/<TARGET>/<ID>_docked.pdbqt
  results/metrics/docking_validation.json
  results/predictions/consensus_hits_docked.csv
"""

import os
import sys
import json
import time
import subprocess
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
import common as C
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

OBABEL = os.path.join(sys.prefix, "Library", "bin", "obabel.exe")
VINA = os.path.join(sys.prefix, "Library", "bin", "vina.exe")
POSE_DIR = os.path.join("results", "docking_poses")

POCKETS = {
    "EGFR": {
        "pdb": "data/targets/egfr_1m17.pdb",
        "pdbqt": "data/targets/egfr_1m17.pdbqt",
        "lig_res": "AQ4",
        "lig_name": "Erlotinib (AQ4)",
        "box_size": [22.0, 22.0, 22.0],
    },
    "FGFR4": {
        "pdb": "data/targets/fgfr4_5jkg.pdb",
        "pdbqt": "data/targets/fgfr4_5jkg.pdbqt",
        "lig_res": "6LF",
        "lig_name": "BLU9931-analog (6LF)",
        "box_size": [22.0, 22.0, 22.0],
    },
    "VEGFR2": {
        "pdb": "data/targets/vegfr2_4asd.pdb",
        "pdbqt": "data/targets/vegfr2_4asd.pdbqt",
        "lig_res": "BAX",
        "lig_name": "Sorafenib (BAX)",
        "box_size": [22.0, 22.0, 22.0],
    },
}


def parse_pdbqt_models(pdbqt_path):
    """Returns list of (affinity_kcal, heavy_atom_coords[N,3], elements[N]) per docked model."""
    models = []
    cur_coords, cur_elems, cur_aff = [], [], None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    cur_aff = float(parts[3])
            elif line.startswith(("ATOM", "HETATM")):
                atype = line[77:79].strip() if len(line) >= 78 else line[12:16].strip()
                if not atype.startswith("H") and atype != "HD":
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    cur_coords.append([x, y, z])
                    cur_elems.append(atype[0].upper())
            elif line.startswith("ENDMDL"):
                if cur_coords:
                    models.append((cur_aff, np.array(cur_coords), cur_elems))
                cur_coords, cur_elems, cur_aff = [], [], None
    if cur_coords:
        models.append((cur_aff, np.array(cur_coords), cur_elems))
    return models


def prepare_receptors_and_redock(exhaustiveness=16):
    os.makedirs( os.path.join(POSE_DIR, "controls"), exist_ok=True)
    redock = {}
    for t, cfg in POCKETS.items():
        pdb_path = cfg["pdb"]
        rec_pdbqt = cfg["pdbqt"]
        lig_res = cfg["lig_res"]
        clean_prot = f"data/targets/{t}_clean_receptor.pdb"
        lig_pdb = os.path.join(POSE_DIR, "controls", f"{t}_{lig_res}_crystal.pdb")
        lig_pdbqt = os.path.join(POSE_DIR, "controls", f"{t}_{lig_res}_crystal.pdbqt")
        out_pdbqt = os.path.join(POSE_DIR, "controls", f"{t}_{lig_res}_redocked.pdbqt")

        coords = []
        with open(pdb_path) as fin, open(clean_prot, "w") as fp, open(lig_pdb, "w") as fl:
            for line in fin:
                if line.startswith("ATOM") and line[21].strip() in ("A", ""):
                    fp.write(line)
                elif line.startswith("HETATM") and line[17:20].strip() == lig_res and line[21].strip() in ("A", ""):
                    fl.write(line)
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            fp.write("END\n")
            fl.write("END\n")

        center = np.round(np.mean(coords, axis=0), 2).tolist()
        cfg["center"] = center

        subprocess.run([OBABEL, "-ipdb", clean_prot, "-opdbqt", "-O", rec_pdbqt,
                        "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([OBABEL, "-ipdb", lig_pdb, "-opdbqt", "-O", lig_pdbqt,
                        "-p", "7.4", "--partialcharge", "gasteiger"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        box = cfg["box_size"]
        cmd = [
            VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
            "--exhaustiveness", str(exhaustiveness), "--seed", "42", "--out", out_pdbqt
        ]
        t0 = time.time()
        if not (os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100):
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ref_models = parse_pdbqt_models(lig_pdbqt)
        dock_models = parse_pdbqt_models(out_pdbqt)
        c_ref = ref_models[0][1]
        rmsds = [float(np.sqrt(np.mean(np.sum((c_ref - m[1]) ** 2, axis=1)))) for m in dock_models]
        top1_aff = dock_models[0][0]
        redock[t] = {
            "ligand": cfg["lig_name"],
            "center": center,
            "box_size": box,
            "affinity_kcal_mol": top1_aff,
            "rmsd_top1_A": round(rmsds[0], 3),
            "rmsd_best_of_9_A": round(min(rmsds), 3),
            "n_heavy_atoms": int(len(c_ref)),
            "secs": round(time.time() - t0, 1),
        }
        print(f"  [Redock {t}] {cfg['lig_name']} @ {center}: {top1_aff:.2f} kcal/mol | "
              f"top-1 RMSD = {rmsds[0]:.2f} A | best-of-9 RMSD = {min(rmsds):.2f} A ({redock[t]['secs']}s)",
              flush=True)
    return redock


def smiles_to_pdbqt(smiles, out_pdbqt):
    """RDKit ETKDGv3 3D embed + MMFF94 minimization -> OpenBabel pH 7.4 PDBQT."""
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
        return True
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = False
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    except Exception:
        pass
    tmp_sdf = out_pdbqt.replace(".pdbqt", ".sdf")
    w = Chem.SDWriter(tmp_sdf)
    w.write(mol); w.close()
    res = subprocess.run([OBABEL, "-isdf", tmp_sdf, "-opdbqt", "-O", out_pdbqt,
                          "-p", "7.4", "--partialcharge", "gasteiger"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if os.path.exists(tmp_sdf):
        os.remove(tmp_sdf)
    return res.returncode == 0 and os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100


def dock_pdbqt(target, lig_pdbqt, out_pdbqt, exhaustiveness=8, cpu=4):
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
        models = parse_pdbqt_models(out_pdbqt)
        if models and models[0][0] is not None:
            return models[0][0]
    cfg = POCKETS[target]
    center, box = cfg["center"], cfg["box_size"]
    cmd = [
        VINA, "--receptor", cfg["pdbqt"], "--ligand", lig_pdbqt,
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
        "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu), "--seed", "42",
        "--out", out_pdbqt
    ]
    res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if res.returncode != 0 or not os.path.exists(out_pdbqt):
        return np.nan
    models = parse_pdbqt_models(out_pdbqt)
    return models[0][0] if (models and models[0][0] is not None) else np.nan


def run_enrichment_control(n_per_class=12, exhaustiveness=8):
    """Dock potent ChEMBL actives (pIC50 >= 8.5) vs weak/inactives (pIC50 <= 4.5), MW-matched."""
    lig_dir = os.path.join(POSE_DIR, "enrichment_ligands")
    os.makedirs(lig_dir, exist_ok=True)
    enrichment = {}
    for t in C.TARGETS:
        t_dir = os.path.join(POSE_DIR, f"{t}_enrichment")
        os.makedirs(t_dir, exist_ok=True)
        df, _, _, _ = C.load_training(t)
        df = df[(df["MolWt"] >= 320) & (df["MolWt"] <= 520) & (df["NumRotatableBonds"] <= 8)].copy()
        act = df[df["pIC50"] >= 8.5].sort_values("pIC50", ascending=False).head(n_per_class)
        inact = df[df["pIC50"] <= 4.8].sort_values("pIC50", ascending=True).head(n_per_class)
        sub = pd.concat([act.assign(label=1), inact.assign(label=0)], ignore_index=True)

        scores, labels = [], []
        t0 = time.time()
        for _, r in sub.iterrows():
            lid = str(r["ligand_id"])
            lpdbqt = os.path.join(lig_dir, f"{lid}.pdbqt")
            opdbqt = os.path.join(t_dir, f"{lid}_docked.pdbqt")
            if not smiles_to_pdbqt(r["canonical_smiles"], lpdbqt):
                continue
            sc = dock_pdbqt(t, lpdbqt, opdbqt, exhaustiveness=exhaustiveness, cpu=8)
            if np.isfinite(sc):
                scores.append(sc)
                labels.append(int(r["label"]))
        scores = np.array(scores); labels = np.array(labels)
        act_sc = scores[labels == 1]; inact_sc = scores[labels == 0]
        # More negative affinity (kcal/mol) = stronger binding -> negate for ROC-AUC
        auc = float(roc_auc_score(labels, -scores))
        pval = float(mannwhitneyu(-act_sc, -inact_sc, alternative="greater").pvalue)
        enrichment[t] = {
            "n_actives": int((labels == 1).sum()),
            "n_inactives": int((labels == 0).sum()),
            "active_mean_kcal": round(float(act_sc.mean()), 2),
            "active_median_kcal": round(float(np.median(act_sc)), 2),
            "inactive_mean_kcal": round(float(inact_sc.mean()), 2),
            "inactive_median_kcal": round(float(np.median(inact_sc)), 2),
            "docking_ROC_AUC": round(auc, 3),
            "mannwhitney_p": round(pval, 4),
            "secs": round(time.time() - t0, 1),
        }
        print(f"  [Enrichment {t}] actives={act_sc.mean():.2f} kcal/mol vs inactives={inact_sc.mean():.2f} kcal/mol | "
              f"Docking AUC = {auc:.3f} (p={pval:.4f}, {enrichment[t]['secs']}s)", flush=True)
    return enrichment


def dock_consensus_candidates(n_novel=25, n_redisc=8, exhaustiveness=8):
    ann = pd.read_csv("results/predictions/consensus_hits_annotated.csv").head(n_novel).copy()
    ann["category"] = "novel_consensus"
    cons = pd.read_csv("results/predictions/consensus_hits.csv")
    # deduplicate rediscoveries by parent_id
    red = cons[cons["rediscovery_any"]].copy()
    red["parent_id"] = red["identifier"].astype(str).str.split(".").str[0]
    red["_has_name"] = red["name"].notna() & (red["name"].astype(str).str.strip() != "")
    p_order = red["parent_id"].drop_duplicates().tolist()
    red = (red.sort_values("_has_name", ascending=False)
              .drop_duplicates("parent_id")
              .set_index("parent_id").loc[p_order].reset_index().head(n_redisc))
    red["category"] = "rediscovery_control"

    from rdkit.Chem import QED
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    p = FilterCatalogParams(); p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog(p)
    r_pains, r_qed = [], []
    for s in red["smiles"]:
        m = Chem.MolFromSmiles(s)
        r_pains.append(bool(pains.HasMatch(m)) if m else True)
        r_qed.append(QED.qed(m) if m else np.nan)
    red["PAINS"] = r_pains
    red["QED"] = r_qed

    keep_cols = ["identifier", "name", "smiles", "category", "n_targets_confident",
                 "n_targets_inAD_active", "mean_pred", "np_likeness",
                 "FGFR4_pred", "FGFR4_lo", "FGFR4_inAD",
                 "EGFR_pred", "EGFR_lo", "EGFR_inAD",
                 "VEGFR2_pred", "VEGFR2_lo", "VEGFR2_inAD",
                 "QED", "PAINS"]
    keep_cols = [c for c in keep_cols if c in ann.columns and c in red.columns]

    batch = pd.concat([ann[keep_cols], red[keep_cols]], ignore_index=True)
    lig_dir = os.path.join(POSE_DIR, "candidate_ligands")
    os.makedirs(lig_dir, exist_ok=True)
    for t in C.TARGETS:
        os.makedirs(os.path.join(POSE_DIR, t), exist_ok=True)
        batch[f"{t}_vina_kcal"] = np.nan

    print(f"\n[+] Docking {len(batch)} candidates ({n_novel} novel consensus + {len(red)} rediscovery controls) across 3 kinases...",
          flush=True)
    t0 = time.time()
    for idx, r in batch.iterrows():
        cid = str(r["identifier"]).replace(".", "_")
        lpdbqt = os.path.join(lig_dir, f"{cid}.pdbqt")
        if not smiles_to_pdbqt(r["smiles"], lpdbqt):
            print(f"    [!] 3D prep failed for {r['identifier']}", flush=True)
            continue
        for t in C.TARGETS:
            opdbqt = os.path.join(POSE_DIR, t, f"{cid}_{t}_docked.pdbqt")
            sc = dock_pdbqt(t, lpdbqt, opdbqt, exhaustiveness=exhaustiveness, cpu=8)
            batch.loc[idx, f"{t}_vina_kcal"] = sc
        if (idx + 1) % 5 == 0 or (idx + 1) == len(batch):
            print(f"    Docked {idx+1}/{len(batch)} compounds ({time.time()-t0:.1f}s)...", flush=True)

    batch["mean_vina_kcal"] = batch[[f"{t}_vina_kcal" for t in C.TARGETS]].mean(axis=1)
    out_csv = os.path.join("results", "predictions", "consensus_hits_docked.csv")
    batch.to_csv(out_csv, index=False)
    print(f"[OK] Saved docked candidates to {out_csv}", flush=True)
    return batch


def main():
    print("=== Step 06: Structure-Based Docking Validation (AutoDock Vina 1.2.5) ===", flush=True)
    redock = prepare_receptors_and_redock(exhaustiveness=16)
    enrichment = run_enrichment_control(n_per_class=12, exhaustiveness=8)
    batch = dock_consensus_candidates(n_novel=25, n_redisc=8, exhaustiveness=8)

    nov = batch[batch["category"] == "novel_consensus"]
    red = batch[batch["category"] == "rediscovery_control"]
    summary = {
        "redocking_control": redock,
        "enrichment_control": enrichment,
        "candidates_summary": {
            "n_novel_docked": int(nov["mean_vina_kcal"].notna().sum()),
            "n_rediscovery_docked": int(red["mean_vina_kcal"].notna().sum()),
            "novel_mean_kcal": round(float(nov["mean_vina_kcal"].mean()), 2),
            "rediscovery_mean_kcal": round(float(red["mean_vina_kcal"].mean()), 2),
            "per_target_novel_mean_kcal": {
                t: round(float(nov[f"{t}_vina_kcal"].mean()), 2) for t in C.TARGETS
            },
        },
    }
    C.save_json(summary, os.path.join(C.METRICS_DIR, "docking_validation.json"))
    print("\n" + json.dumps(summary, indent=2), flush=True)
    return {"status": "docking validation complete"}


if __name__ == "__main__":
    main()
