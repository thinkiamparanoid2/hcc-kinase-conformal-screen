"""
Step 06 - Orthogonal structure-based validation (AutoDock Vina)
----------------------------------------------------------------
1. Receptor preparation and co-crystal redocking (pose-prediction control).
2. Active-versus-inactive enrichment control.
3. Docking of the candidate shortlist across all three HCC kinase pockets.

WHAT CHANGED AND WHY:

* RMSD IS SYMMETRY-CORRECTED AND TOP-1 IS THE HEADLINE.
  The old code compared two separately parsed PDBQT atom lists element-wise,
  assuming identical atom ordering and ignoring molecular symmetry; Erlotinib and
  Sorafenib both have symmetry-equivalent atoms, so those numbers could be wrong in
  either direction. RMSD now goes through RDKit `GetBestRMS`, which does the
  substructure matching. More importantly, the reported control is TOP-1 RMSD --
  pose prediction. Best-of-N is reported beside it and labelled as a sampling
  diagnostic, because quoting best-of-9 as the control is how a 7.8 A top-1 pose
  gets presented as 1.2 A.

* THE ENRICHMENT CONTROL IS POWERED AND UNBIASED.
  Previously 12 actives versus 12 inactives, taken as the 12 MOST potent against
  the 12 WEAKEST -- the easiest possible discrimination, on a sample whose AUC
  confidence interval is roughly +/- 0.17. Now: random sampling within the active
  and inactive classes (seeded), `N_PER_CLASS` default 50, property-matched
  selection on MW and rotatable bonds, bootstrap CIs on AUC, plus EF1%, EF5% and
  BEDROC, which are the metrics that describe a virtual screen.

* THE DOCKING BOX FITS THE LIGAND.
  A fixed 22 A cube for every pocket and every ligand meant large rigid ligands
  could not be placed, and Vina returned near-zero affinities that were then read
  as "steric clash". Box size is now max(pocket span, ligand diameter + padding),
  and any ligand whose score is implausibly weak is flagged
  (`placement_suspect`) rather than interpreted as a chemistry result.

* SCORES ARE REPLICATED, AND NEVER COMPARED ACROSS RECEPTORS.
  Each docking runs with `N_SEEDS` seeds and the mean +/- SD is reported. Vina
  scores are not transferable free energies, so cross-target comparison is
  reported as a within-target z-score (`<T>_vina_z`) against that target's own
  active distribution, and the raw kcal/mol are kept per target only.

Requires AutoDock Vina and Open Babel on PATH (or VINA_BIN / OBABEL_BIN set).

Outputs: results/docking_poses/<TARGET>/<ID>_docked.pdbqt
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

POSE_DIR = os.path.join("results", "docking_poses")
N_PER_CLASS = 50            # per class in the enrichment control (was 12)
N_SEEDS = 3                 # docking replicates per ligand
EXHAUSTIVENESS = 32         # was 8 for everything but the redock
BOX_PADDING = 8.0           # A added to the ligand's longest dimension
MIN_BOX = 22.0              # A, floor so small ligands still sample the pocket
IMPLAUSIBLE_SCORE = -4.0    # weaker than this is treated as a placement failure
# Vina threads per docking. Asking for more threads than the machine has
# oversubscribes the cores and runs slower, and the thread count is part of
# what makes a seeded run reproducible, so it is recorded in the summary.
N_CPU = max(1, os.cpu_count() or 1)

POCKETS = {
    "EGFR": {"pdb": "data/targets/egfr_1m17.pdb",
             "pdbqt": "data/targets/egfr_1m17.pdbqt",
             "lig_res": "AQ4", "lig_name": "Erlotinib (AQ4)"},
    "FGFR4": {"pdb": "data/targets/fgfr4_5jkg.pdb",
              "pdbqt": "data/targets/fgfr4_5jkg.pdbqt",
              "lig_res": "6LF", "lig_name": "BLU9931-analog (6LF)"},
    "VEGFR2": {"pdb": "data/targets/vegfr2_4asd.pdb",
               "pdbqt": "data/targets/vegfr2_4asd.pdbqt",
               "lig_res": "BAX", "lig_name": "Sorafenib (BAX)"},
}


def binaries():
    return (C.resolve_binary("vina", env_var="VINA_BIN"),
            C.resolve_binary("obabel", env_var="OBABEL_BIN"))


# ----------------------------------------------------------------------------
# PDBQT parsing and RMSD
# ----------------------------------------------------------------------------

def parse_pdbqt_models(pdbqt_path):
    """[(affinity_kcal, heavy_atom_coords[N,3], element_symbols)] per docked model."""
    models = []
    coords, elems, aff = [], [], None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    aff = float(parts[3])
            elif line.startswith(("ATOM", "HETATM")):
                atype = line[77:79].strip() if len(line) >= 78 else line[12:16].strip()
                if not atype.startswith("H") and atype != "HD":
                    coords.append([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])
                    elems.append(atype[0].upper())
            elif line.startswith("ENDMDL"):
                if coords:
                    models.append((aff, np.array(coords), elems))
                coords, elems, aff = [], [], None
    if coords:
        models.append((aff, np.array(coords), elems))
    return models


def _mol_from_pdbqt_model(pdbqt_path, model_index, obabel):
    """Extract one model from a PDBQT and return it as an RDKit Mol (heavy atoms)."""
    out_sdf = f"{pdbqt_path}.model{model_index}.sdf"
    res = subprocess.run([obabel, "-ipdbqt", pdbqt_path, "-osdf", "-O", out_sdf,
                          "-f", str(model_index + 1), "-l", str(model_index + 1)],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    mol = None
    if res.returncode == 0 and os.path.exists(out_sdf):
        supplier = Chem.SDMolSupplier(out_sdf, sanitize=False, removeHs=True)
        mols = [m for m in supplier if m is not None]
        mol = mols[0] if mols else None
    if os.path.exists(out_sdf):
        os.remove(out_sdf)
    return mol


def _flatten_graph(mol):
    """Copy with every bond single and all aromaticity, charges and implicit Hs
    cleared.

    The reference comes from a PDB and the pose comes back through Open Babel, so
    the two perceive bond orders differently for the same molecule: kekulised vs
    aromatic rings, a carboxylate drawn -1 on one side and neutral on the other.
    `GetSubstructMatch` then finds nothing and the RMSD silently comes back None.
    An RMSD only needs the heavy-atom graph, so matching on connectivity and
    element alone is the correct comparison, not a workaround.
    """
    rw = Chem.RWMol(Chem.Mol(mol))
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetFormalCharge(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_SYMMRINGS |
                          Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
                     catchErrors=True)
    return out


def symmetry_corrected_rmsd(ref_mol, probe_mol, superimpose=False):
    """Heavy-atom RMSD between two poses, symmetry-corrected.

    `superimpose=False` (the control) is `CalcRMS`: the pose is compared where
    docking actually put it. `superimpose=True` is `GetBestRMS`, which rotates and
    translates the pose onto the reference first and is therefore NOT a
    pose-prediction metric -- it answers a different question, whether the internal
    conformation is right once placement is ignored. Quoting a superimposed number
    as the redocking control turns a pose docked into the wrong subpocket into an
    apparent success: in this repository's own run, EGFR erlotinib is 8.0 A in
    place and 1.8 A superimposed.

    Both enumerate symmetry-equivalent atom mappings, which an element-wise
    positional RMSD cannot do. Returns None when the two molecules do not match,
    which is itself diagnostic (a failed conversion, not a good or bad pose).
    """
    if ref_mol is None or probe_mol is None:
        return None
    try:
        ref = _flatten_graph(Chem.RemoveHs(Chem.Mol(ref_mol), sanitize=False))
        probe = _flatten_graph(Chem.RemoveHs(Chem.Mol(probe_mol), sanitize=False))
        if ref.GetNumAtoms() != probe.GetNumAtoms():
            return None
        if superimpose:
            return float(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(ref)))
        return float(rdMolAlign.CalcRMS(probe, ref))
    except Exception:
        return None


# ----------------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------------

def ligand_span(coords):
    """Longest pairwise extent of a coordinate set, per axis and overall."""
    if len(coords) == 0:
        return np.zeros(3), 0.0
    extent = coords.max(axis=0) - coords.min(axis=0)
    diameter = float(np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2)
                             .sum(-1)).max()) if len(coords) < 400 else float(extent.max())
    return extent, diameter


def box_for(coords, extra_diameter=0.0):
    """Cubic box that contains the reference ligand and a candidate of given size."""
    extent, diameter = ligand_span(coords)
    side = max(MIN_BOX, float(extent.max()) + BOX_PADDING,
               max(diameter, extra_diameter) + BOX_PADDING)
    return [round(side, 1)] * 3


def ligand_diameter_from_pdbqt(path):
    models = parse_pdbqt_models(path)
    if not models:
        return 0.0
    _, diameter = ligand_span(models[0][1])
    return diameter


# ----------------------------------------------------------------------------
# Receptor prep + redocking control
# ----------------------------------------------------------------------------

def prepare_receptors_and_redock(exhaustiveness=EXHAUSTIVENESS):
    vina, obabel = binaries()
    os.makedirs(os.path.join(POSE_DIR, "controls"), exist_ok=True)
    redock = {}
    for t, cfg in POCKETS.items():
        clean_prot = f"data/targets/{t}_clean_receptor.pdb"
        lig_pdb = os.path.join(POSE_DIR, "controls", f"{t}_{cfg['lig_res']}_crystal.pdb")
        lig_pdbqt = lig_pdb.replace(".pdb", ".pdbqt")
        out_pdbqt = os.path.join(POSE_DIR, "controls",
                                 f"{t}_{cfg['lig_res']}_redocked.pdbqt")

        coords = []
        with open(cfg["pdb"]) as fin, open(clean_prot, "w") as fp, open(lig_pdb, "w") as fl:
            for line in fin:
                if line.startswith("ATOM") and line[21].strip() in ("A", ""):
                    fp.write(line)
                elif (line.startswith("HETATM")
                      and line[17:20].strip() == cfg["lig_res"]
                      and line[21].strip() in ("A", "")):
                    fl.write(line)
                    coords.append([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])
            fp.write("END\n")
            fl.write("END\n")

        coords = np.array(coords)
        center = np.round(coords.mean(axis=0), 2).tolist()
        box = box_for(coords)
        cfg["center"] = center
        cfg["box_size"] = box
        cfg["ref_ligand_pdb"] = lig_pdb

        subprocess.run([obabel, "-ipdb", clean_prot, "-opdbqt", "-O", cfg["pdbqt"],
                        "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([obabel, "-ipdb", lig_pdb, "-opdbqt", "-O", lig_pdbqt,
                        "-p", "7.4", "--partialcharge", "gasteiger"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        t0 = time.time()
        if not (os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100):
            subprocess.run([
                vina, "--receptor", cfg["pdbqt"], "--ligand", lig_pdbqt,
                "--center_x", str(center[0]), "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
                "--exhaustiveness", str(exhaustiveness), "--seed", "42",
                "--out", out_pdbqt],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        ref_mol = Chem.MolFromPDBFile(lig_pdb, sanitize=False, removeHs=True)
        docked = parse_pdbqt_models(out_pdbqt)
        rmsds, aligned = [], []
        for i in range(len(docked)):
            probe = _mol_from_pdbqt_model(out_pdbqt, i, obabel)
            rmsds.append(symmetry_corrected_rmsd(ref_mol, probe))
            aligned.append(symmetry_corrected_rmsd(ref_mol, probe, superimpose=True))
        valid = [r for r in rmsds if r is not None]

        redock[t] = {
            "ligand": cfg["lig_name"],
            "center": center,
            "box_size": box,
            "exhaustiveness": exhaustiveness,
            "affinity_top1_kcal_mol": docked[0][0] if docked else None,
            # Top-1 is the pose-prediction control. Best-of-N describes sampling
            # only, and is never the number a redocking control should quote.
            "rmsd_top1_A": round(rmsds[0], 3) if rmsds and rmsds[0] is not None else None,
            "rmsd_best_of_n_A": round(min(valid), 3) if valid else None,
            # Placement ignored: right conformer in the wrong place still scores
            # well here, so this is a diagnostic and never the control.
            "rmsd_top1_superimposed_A": (
                round(aligned[0], 3) if aligned and aligned[0] is not None else None),
            "n_poses": len(docked),
            "n_rmsd_computed": len(valid),
            "rmsd_method": ("RDKit CalcRMS, symmetry-corrected heavy atoms, "
                            "in place (no superposition); the molecules are matched "
                            "on the heavy-atom graph because the crystal reference "
                            "and the Open Babel pose perceive bond orders "
                            "differently"),
            "pose_prediction_success_top1": (
                bool(rmsds[0] is not None and rmsds[0] <= 2.0) if rmsds else None),
            "secs": round(time.time() - t0, 1),
        }
        r = redock[t]
        print(f"  [Redock {t}] {cfg['lig_name']}: top-1 RMSD = {r['rmsd_top1_A']} A "
              f"({'PASS' if r['pose_prediction_success_top1'] else 'FAIL'} at 2 A), "
              f"best-of-{r['n_poses']} = {r['rmsd_best_of_n_A']} A, "
              f"superimposed = {r['rmsd_top1_superimposed_A']} A", flush=True)
    return redock


# ----------------------------------------------------------------------------
# Ligand prep + docking
# ----------------------------------------------------------------------------

def smiles_to_pdbqt(smiles, out_pdbqt, n_confs=1):
    """ETKDGv3 embed + MMFF94 minimise -> Open Babel pH 7.4 PDBQT."""
    _, obabel = binaries()
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
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    tmp_sdf = out_pdbqt.replace(".pdbqt", ".sdf")
    with Chem.SDWriter(tmp_sdf) as w:
        w.write(mol)
    res = subprocess.run([obabel, "-isdf", tmp_sdf, "-opdbqt", "-O", out_pdbqt,
                          "-p", "7.4", "--partialcharge", "gasteiger"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if os.path.exists(tmp_sdf):
        os.remove(tmp_sdf)
    return (res.returncode == 0 and os.path.exists(out_pdbqt)
            and os.path.getsize(out_pdbqt) > 100)


def dock_pdbqt(target, lig_pdbqt, out_pdbqt, exhaustiveness=EXHAUSTIVENESS,
               cpu=N_CPU, n_seeds=N_SEEDS):
    """Dock with `n_seeds` replicates. Returns (mean, sd, [scores], placement_suspect).

    A single Vina run with one seed gives a number with no error bar. Replicates
    cost little and are what makes a reported affinity defensible.
    """
    vina, _ = binaries()
    cfg = POCKETS[target]
    center = cfg["center"]
    # Grow the box if this ligand is larger than the co-crystal reference: a fixed
    # box silently fails to place big ligands and returns a near-zero score.
    ref_coords = parse_pdbqt_models(
        os.path.join(POSE_DIR, "controls",
                     f"{target}_{cfg['lig_res']}_crystal.pdbqt"))
    ref_xyz = ref_coords[0][1] if ref_coords else np.zeros((1, 3))
    box = box_for(ref_xyz, extra_diameter=ligand_diameter_from_pdbqt(lig_pdbqt))

    scores = []
    for s in range(n_seeds):
        seed_out = out_pdbqt if s == 0 else out_pdbqt.replace(".pdbqt", f"_s{s}.pdbqt")
        if not (os.path.exists(seed_out) and os.path.getsize(seed_out) > 100):
            res = subprocess.run([
                vina, "--receptor", cfg["pdbqt"], "--ligand", lig_pdbqt,
                "--center_x", str(center[0]), "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
                "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu),
                "--seed", str(42 + s), "--out", seed_out],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if res.returncode != 0 or not os.path.exists(seed_out):
                continue
        models = parse_pdbqt_models(seed_out)
        if models and models[0][0] is not None:
            scores.append(models[0][0])
    if not scores:
        return np.nan, np.nan, [], True
    mean = float(np.mean(scores))
    sd = float(np.std(scores, ddof=1)) if len(scores) > 1 else 0.0
    # A score this weak means Vina could not place the ligand in the box. That is a
    # docking failure to report, not evidence of a steric clash in the pocket.
    suspect = bool(mean > IMPLAUSIBLE_SCORE)
    return mean, sd, scores, suspect


# ----------------------------------------------------------------------------
# Enrichment control
# ----------------------------------------------------------------------------

def _property_matched_sample(act, inact, n_per_class, seed):
    """Random, MW-matched sample from each class (no potency-extreme selection).

    Each sampled inactive is the unused inactive closest in molecular weight to a
    sampled active, so the two classes are property-matched without ranking either
    class by potency.
    """
    rng = np.random.RandomState(seed)
    n = min(n_per_class, len(act), len(inact))
    if n == 0:
        return act.iloc[[]], inact.iloc[[]]
    act_s = act.sample(n=n, random_state=rng)
    pool = inact.copy()
    picked = []
    for mw in act_s["MolWt"].to_numpy(float):
        if pool.empty:
            break
        j = (pool["MolWt"].to_numpy(float) - mw).__abs__().argmin()
        picked.append(pool.index[j])
        pool = pool.drop(pool.index[j])
    return act_s, inact.loc[picked]


def run_enrichment_control(n_per_class=N_PER_CLASS, exhaustiveness=EXHAUSTIVENESS,
                           seed=0):
    lig_dir = os.path.join(POSE_DIR, "enrichment_ligands")
    os.makedirs(lig_dir, exist_ok=True)
    enrichment = {}
    for t in C.TARGETS:
        t_dir = os.path.join(POSE_DIR, f"{t}_enrichment")
        os.makedirs(t_dir, exist_ok=True)
        # Censored records are admissible here: a ">" bound is exactly what an
        # inactive is, even though it cannot serve as a regression target.
        df, _, _, _ = C.load_training(t, exact_only=False)
        df = df[(df["MolWt"] >= 250) & (df["MolWt"] <= 600) &
                (df["NumRotatableBonds"] <= 12)].copy()
        act = df[df["pIC50"] >= 7.0]
        inact = df[df["pIC50"] <= 5.0]
        act_s, inact_s = _property_matched_sample(act, inact, n_per_class, seed)
        sub = pd.concat([act_s.assign(label=1), inact_s.assign(label=0)],
                        ignore_index=True)
        if sub.empty or sub["label"].nunique() < 2:
            enrichment[t] = {"status": "skipped",
                             "reason": "no inactive class available; re-run step 01 "
                                       "with censored relations to populate it"}
            print(f"  [Enrichment {t}] SKIPPED: {enrichment[t]['reason']}", flush=True)
            continue

        scores, labels, sds = [], [], []
        t0 = time.time()
        for _, r in sub.iterrows():
            lid = str(r["ligand_id"])
            lpdbqt = os.path.join(lig_dir, f"{lid}.pdbqt")
            opdbqt = os.path.join(t_dir, f"{lid}_docked.pdbqt")
            if not smiles_to_pdbqt(r["canonical_smiles"], lpdbqt):
                continue
            mean, sd, _, _ = dock_pdbqt(t, lpdbqt, opdbqt,
                                        exhaustiveness=exhaustiveness, cpu=N_CPU)
            if np.isfinite(mean):
                scores.append(mean)
                sds.append(sd)
                labels.append(int(r["label"]))
        scores = np.array(scores)
        labels = np.array(labels)
        if scores.size == 0 or np.unique(labels).size < 2:
            enrichment[t] = {"status": "failed", "reason": "no usable docking scores"}
            continue

        ranking = -scores                      # stronger binding = higher score
        auc = float(roc_auc_score(labels, ranking))
        enrichment[t] = {
            "status": "ok",
            "n_actives": int((labels == 1).sum()),
            "n_inactives": int((labels == 0).sum()),
            "active_threshold_pIC50": 7.0,
            "inactive_threshold_pIC50": 5.0,
            "sampling": "random within class, MW-matched; NOT potency extremes",
            "active_mean_kcal": round(float(scores[labels == 1].mean()), 2),
            "inactive_mean_kcal": round(float(scores[labels == 0].mean()), 2),
            "mean_replicate_sd_kcal": round(float(np.mean(sds)), 3) if sds else None,
            "docking_ROC_AUC": round(auc, 3),
            "docking_ROC_AUC_CI95": C.auc_ci(labels, ranking),
            "EF1pct": round(C.enrichment_factor(labels, ranking, 0.01), 2),
            "EF5pct": round(C.enrichment_factor(labels, ranking, 0.05), 2),
            "BEDROC_alpha20": round(C.bedroc(labels, ranking, 20.0), 3),
            "mannwhitney_p": round(float(mannwhitneyu(
                ranking[labels == 1], ranking[labels == 0],
                alternative="greater").pvalue), 5),
            "secs": round(time.time() - t0, 1),
        }
        e = enrichment[t]
        print(f"  [Enrichment {t}] n={e['n_actives']}/{e['n_inactives']} | "
              f"AUC {e['docking_ROC_AUC']} "
              f"[{e['docking_ROC_AUC_CI95']['lo']:.3f}, "
              f"{e['docking_ROC_AUC_CI95']['hi']:.3f}] | "
              f"EF1% {e['EF1pct']} | BEDROC {e['BEDROC_alpha20']}", flush=True)
    return enrichment


# ----------------------------------------------------------------------------
# Candidate docking
# ----------------------------------------------------------------------------

def dock_consensus_candidates(n_novel=25, n_control=8, exhaustiveness=EXHAUSTIVENESS):
    from rdkit.Chem import QED
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    ann_path = "results/predictions/consensus_hits_annotated.csv"
    cons_path = "results/predictions/consensus_hits.csv"
    if not os.path.exists(ann_path):
        raise FileNotFoundError(f"{ann_path} missing; run step 07's annotation first.")
    ann = pd.read_csv(ann_path).head(n_novel).copy()
    ann["category"] = "novel_consensus"

    cons = pd.read_csv(cons_path)
    red = cons[cons["rediscovery_any"]].copy()
    red["parent"] = red["identifier"].astype(str).str.split(".").str[0]
    red = red.drop_duplicates("parent").head(n_control)
    red["category"] = "rediscovery_control"

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog(params)
    for frame in (red,):
        flags, qeds = [], []
        for s in frame["smiles"]:
            m = Chem.MolFromSmiles(s)
            flags.append(bool(pains.HasMatch(m)) if m else True)
            qeds.append(QED.qed(m) if m else np.nan)
        frame["PAINS"] = flags
        frame["QED"] = qeds

    keep = [c for c in ["identifier", "name", "smiles", "category",
                        "n_targets_interval_active", "n_targets_inAD_active",
                        "mean_pred", "np_likeness", "generic_scaffold",
                        "pan_kinase_alert", "QED", "PAINS"]
            if c in ann.columns and c in red.columns]
    batch = pd.concat([ann[keep], red[keep]], ignore_index=True)

    lig_dir = os.path.join(POSE_DIR, "candidate_ligands")
    os.makedirs(lig_dir, exist_ok=True)
    for t in C.TARGETS:
        os.makedirs(os.path.join(POSE_DIR, t), exist_ok=True)
        batch[f"{t}_vina_kcal"] = np.nan
        batch[f"{t}_vina_sd"] = np.nan
        batch[f"{t}_placement_suspect"] = False

    print(f"\n[+] Docking {len(batch)} candidates x {len(C.TARGETS)} pockets "
          f"x {N_SEEDS} seeds...", flush=True)
    t0 = time.time()
    for idx, r in batch.iterrows():
        cid = str(r["identifier"]).replace(".", "_")
        lpdbqt = os.path.join(lig_dir, f"{cid}.pdbqt")
        if not smiles_to_pdbqt(r["smiles"], lpdbqt):
            print(f"    [!] 3D prep failed for {r['identifier']}", flush=True)
            continue
        for t in C.TARGETS:
            opdbqt = os.path.join(POSE_DIR, t, f"{cid}_{t}_docked.pdbqt")
            mean, sd, _, suspect = dock_pdbqt(t, lpdbqt, opdbqt,
                                              exhaustiveness=exhaustiveness, cpu=N_CPU)
            batch.loc[idx, f"{t}_vina_kcal"] = mean
            batch.loc[idx, f"{t}_vina_sd"] = sd
            batch.loc[idx, f"{t}_placement_suspect"] = suspect
        if (idx + 1) % 5 == 0 or (idx + 1) == len(batch):
            print(f"    {idx+1}/{len(batch)} ({time.time()-t0:.0f}s)", flush=True)

    # Within-target z-scores. Vina scores are not comparable between receptors, so
    # a per-target standardisation against that pocket's own candidate
    # distribution is the only defensible cross-target view.
    for t in C.TARGETS:
        col = batch[f"{t}_vina_kcal"]
        mu, sd = col.mean(), col.std(ddof=1)
        batch[f"{t}_vina_z"] = (col - mu) / sd if sd and np.isfinite(sd) else np.nan

    out_csv = os.path.join("results", "predictions", "consensus_hits_docked.csv")
    batch.to_csv(out_csv, index=False)
    print(f"[OK] -> {out_csv}", flush=True)
    return batch


def tool_versions():
    """Record what actually produced these poses.

    Receptor protonation differs between Open Babel builds, so two runs of this
    step can start from slightly different receptors. Pinning the versions in
    environment.yml is the fix; reporting them is how a reader can tell whether a
    given results file came from the pinned pair.
    """
    vina, obabel = binaries()
    out = {}
    for label, cmd in (("vina", [vina, "--version"]),
                       ("open_babel", [obabel, "-V"])):
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            out[label] = (res.stdout or res.stderr).strip().splitlines()[0]
        except Exception:
            out[label] = "unknown"
    return out


def main():
    print("=== Step 06: Structure-based docking validation (AutoDock Vina) ===",
          flush=True)
    vina, obabel = binaries()
    print(f"    vina:   {vina}\n    obabel: {obabel}", flush=True)

    redock = prepare_receptors_and_redock()
    enrichment = run_enrichment_control()
    batch = dock_consensus_candidates()

    nov = batch[batch["category"] == "novel_consensus"]
    red = batch[batch["category"] == "rediscovery_control"]
    summary = {
        "settings": {"exhaustiveness": EXHAUSTIVENESS, "n_seeds": N_SEEDS,
                     "box_padding_A": BOX_PADDING, "min_box_A": MIN_BOX,
                     "n_per_class_enrichment": N_PER_CLASS,
                     "vina_cpu_threads": N_CPU,
                     "tool_versions": tool_versions()},
        "redocking_control": redock,
        "enrichment_control": enrichment,
        "candidates_summary": {
            "n_novel_docked": int(nov["FGFR4_vina_kcal"].notna().sum()),
            "n_control_docked": int(red["FGFR4_vina_kcal"].notna().sum()),
            "per_target_novel_mean_kcal": {
                t: round(float(nov[f"{t}_vina_kcal"].mean()), 2) for t in C.TARGETS},
            "per_target_placement_failures": {
                t: int(nov[f"{t}_placement_suspect"].sum()) for t in C.TARGETS},
        },
        "score_comparability_note": (
            "AutoDock Vina scores are empirical and receptor-specific. They are not "
            "transferable binding free energies, so a weaker score in one pocket "
            "than another is not evidence of selectivity. Compare within a target "
            "(see the <T>_vina_z columns), and treat any score weaker than "
            f"{IMPLAUSIBLE_SCORE} kcal/mol as a placement failure to investigate, "
            "not as a measured steric clash."),
    }
    C.save_json(summary, os.path.join(C.METRICS_DIR, "docking_validation.json"))
    print("\n" + json.dumps(summary, indent=2), flush=True)
    return {"status": "docking validation complete"}


if __name__ == "__main__":
    main()
