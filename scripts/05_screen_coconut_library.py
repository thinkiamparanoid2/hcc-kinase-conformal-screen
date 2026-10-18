"""
Step 05 - Full-library screen (all ~738k COCONUT natural products)
-------------------------------------------------------------------
Streams the whole library in parallel chunks (no 3 GB feature CSV): featurise ->
predict pIC50 with the final XGB models -> attach Mondrian-conformal 90% interval
and applicability-domain flag. For each target writes:

  results/predictions/<T>_screen_full.csv.gz
      identifier, name, smiles, np_likeness, pred, lo, hi, sim_train, sim_all,
      nn_id, nn_pIC50, in_AD, rediscovery

  rediscovery = max Tanimoto >= 0.90 to ANY ChEMBL training molecule for that
                target (known chemistry present in COCONUT -> positive control,
                NOT a novel lead)
  in_AD       = max Tanimoto to the model's train split >= tau (from step 04)

Then results/predictions/consensus_hits.csv and results/metrics/screen_summary.json.
"""

import os
import sys
import time
import json
import numpy as np
import pandas as pd
import multiprocessing as mp

sys.path.insert(0, "scripts")
import common as C

LIB_CSV = os.path.join("data", "ligands_phytochemical", "phytochemical_screening_set.csv")
OUT_DIR = os.path.join("results", "predictions")
CHUNK = 2000
REDISCOVERY_SIM = 0.90

# worker globals
_W = {}


def _init_worker(train_smiles_by_target, in_train_by_target, desc_cols):
    sys.path.insert(0, "scripts")
    import importlib.util
    from rdkit.ML.Descriptors import MoleculeDescriptors
    spec = importlib.util.spec_from_file_location("feat_eng", "scripts/02_feature_engineering.py")
    fe = importlib.util.module_from_spec(spec); spec.loader.exec_module(fe)
    import common as CC
    _W["fe"] = fe
    _W["C"] = CC
    _W["desc_cols"] = desc_cols
    _W["physchem_names"] = list(fe.PHYSCHEM_DESCRIPTOR_NAMES)
    _W["calc"] = MoleculeDescriptors.MolecularDescriptorCalculator(_W["physchem_names"])
    _W["fps"] = {t: CC.morgan_bitvects(s) for t, s in train_smiles_by_target.items()}
    _W["in_train"] = in_train_by_target


def _featurize(args):
    """smiles chunk -> (start, ok_mask, desc[n,D], fp[n,2048], sim_all, nn_idx, sim_tr)."""
    start, smiles = args
    from rdkit import Chem, DataStructs
    CC, calc, names, desc_cols = _W["C"], _W["calc"], _W["physchem_names"], _W["desc_cols"]
    n = len(smiles)
    ok = np.zeros(n, bool)
    desc = np.full((n, len(desc_cols)), np.nan, np.float32)
    fp = np.zeros((n, 2048), np.uint8)
    targets = list(_W["fps"].keys())
    sim_all_c = {t: np.full(n, np.nan, np.float32) for t in targets}
    sim_tr_c = {t: np.full(n, np.nan, np.float32) for t in targets}
    nn_idx_c = {t: np.full(n, -1, np.int32) for t in targets}
    for i, smi in enumerate(smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            vals = calc.CalcDescriptors(mol)
            d = dict(zip(names, vals))
            d["LipinskiViolations"] = int(
                (d.get("MolWt", 0) > 500) + (d.get("MolLogP", 0) > 5) +
                (d.get("NumHDonors", 0) > 5) + (d.get("NumHAcceptors", 0) > 10)
            )
            desc[i] = [d[c] for c in desc_cols]
            bv = CC._MORGAN.GetFingerprint(mol)
            fp[i] = CC._MORGAN.GetFingerprintAsNumPy(mol)
            for t in targets:
                s = np.asarray(DataStructs.BulkTanimotoSimilarity(bv, _W["fps"][t]), dtype=np.float32)
                idx_max = int(s.argmax())
                sim_all_c[t][i] = s[idx_max]
                nn_idx_c[t][i] = idx_max
                sim_tr_c[t][i] = s[_W["in_train"][t]].max()
            ok[i] = True
        except Exception:
            continue
    return start, ok, desc, fp, sim_all_c, nn_idx_c, sim_tr_c


def main():
    t_all = time.time()
    from xgboost import XGBRegressor
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- load per-target artifacts from step 04 ----
    art = {}
    for t in C.TARGETS:
        with open(os.path.join(C.MODELS_DIR, f"{t}_conformal.json")) as f:
            cj = json.load(f)
        df, cols, fp_cols, desc_cols = C.load_training(t)
        smi = df["canonical_smiles"].tolist()
        in_train = np.isin(smi, cj["train_smiles"])
        m = XGBRegressor(); m.load_model(os.path.join(C.MODELS_DIR, f"{t}_xgb_final.json"))
        art[t] = {"cj": cj, "smiles": smi, "ids": df["ligand_id"].astype(str).tolist(),
                  "y": df["pIC50"].to_numpy(float), "in_train": in_train, "model": m,
                  "cols": cj["feature_cols"]}
    cols = art[C.TARGETS[0]]["cols"]
    assert all(a["cols"] == cols for a in art.values()), "feature columns differ across targets"
    desc_cols = [c for c in cols if not c.startswith("FP_")]
    assert cols == desc_cols + [f"FP_{i}" for i in range(2048)], "unexpected feature order"

    # ---- library ----
    lib = pd.read_csv(LIB_CSV, usecols=["identifier", "name", "canonical_smiles", "np_likeness"])
    lib = lib.dropna(subset=["canonical_smiles"]).reset_index(drop=True)
    print(f"[+] Library: {len(lib):,} compounds", flush=True)
    smiles = lib["canonical_smiles"].tolist()
    chunks = [(s, smiles[s:s + CHUNK]) for s in range(0, len(smiles), CHUNK)]

    n = len(lib)
    pred = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    sim_train = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    sim_all = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    nn_idx = {t: np.full(n, -1, np.int32) for t in C.TARGETS}
    valid = np.zeros(n, bool)

    train_smiles_by_target = {t: art[t]["smiles"] for t in C.TARGETS}
    in_train_by_target = {t: art[t]["in_train"] for t in C.TARGETS}
    nproc = max(1, min(14, mp.cpu_count() - 2))
    print(f"[+] Featurising + similarity with {nproc} workers, {len(chunks)} chunks", flush=True)
    t0 = time.time()
    with mp.Pool(nproc, initializer=_init_worker,
                 initargs=(train_smiles_by_target, in_train_by_target, desc_cols)) as pool:
        for k, (start, ok, desc, fp, sim_all_c, nn_idx_c, sim_tr_c) in enumerate(pool.imap(_featurize, chunks)):
            e = start + len(ok)
            valid[start:e] = ok
            X = np.concatenate([desc, fp.astype(np.float32)], axis=1)
            X[~np.isfinite(X)] = np.nan
            for t in C.TARGETS:
                pred[t][start:e] = art[t]["model"].predict(X)
                sim_all[t][start:e] = sim_all_c[t]
                nn_idx[t][start:e] = nn_idx_c[t]
                sim_train[t][start:e] = sim_tr_c[t]
            if (k + 1) % 20 == 0 or k + 1 == len(chunks):
                el = time.time() - t0
                print(f"    {k+1}/{len(chunks)} chunks | {el/60:.1f} min | ETA {(el/(k+1))*(len(chunks)-k-1)/60:.1f} min",
                      flush=True)

    # ---- assemble per-target tables ----
    summary = {"library_total": int(n), "library_valid": int(valid.sum()), "targets": {}}
    tables = {}
    for t in C.TARGETS:
        a = art[t]; cj = a["cj"]
        edges, hw, tau = cj["edges"], np.array(cj["half_widths"]), cj["tau"]
        b = np.clip(np.digitize(sim_train[t], edges[1:-1]), 0, len(edges) - 2)
        half = hw[b]
        nn = np.where(nn_idx[t] >= 0, nn_idx[t], 0)
        d = pd.DataFrame({
            "identifier": lib["identifier"], "name": lib["name"], "smiles": lib["canonical_smiles"],
            "np_likeness": lib["np_likeness"],
            "pred": pred[t], "lo": pred[t] - half, "hi": pred[t] + half,
            "sim_train": sim_train[t], "sim_all": sim_all[t],
            "nn_id": np.array(a["ids"])[nn], "nn_pIC50": a["y"][nn],
        })
        d["in_AD"] = (d["sim_train"] >= tau) if tau is not None else False
        d["rediscovery"] = d["sim_all"] >= REDISCOVERY_SIM
        d = d[valid].sort_values("pred", ascending=False).reset_index(drop=True)
        d.to_csv(os.path.join(OUT_DIR, f"{t}_screen_full.csv.gz"), index=False)
        tables[t] = d
        conf = d["in_AD"] & (d["lo"] >= C.ACTIVE_THRESHOLD)
        in_ad_act = d["in_AD"] & (d["pred"] >= C.ACTIVE_THRESHOLD)
        summary["targets"][t] = {
            "tau": tau, "in_AD": int(d["in_AD"].sum()),
            "in_AD_frac": float(d["in_AD"].mean()),
            "median_sim_train": float(d["sim_train"].median()),
            "pred_ge_7": int((d["pred"] >= 7).sum()),
            "in_AD_active (in_AD & pred>=6)": int(in_ad_act.sum()),
            "in_AD_active_not_rediscovery": int((in_ad_act & ~d["rediscovery"]).sum()),
            "confident_active (in_AD & lo>=6)": int(conf.sum()),
            "confident_active_not_rediscovery": int((conf & ~d["rediscovery"]).sum()),
            "rediscovery_total": int(d["rediscovery"].sum()),
        }
        print(f"  {t}: in_AD {d['in_AD'].sum():,} ({d['in_AD'].mean()*100:.2f}%) | "
              f"in_AD&pred>=6: {int(in_ad_act.sum())} (novel {int((in_ad_act & ~d['rediscovery']).sum())}) | "
              f"confident (lo>=6): {int(conf.sum())} (novel {int((conf & ~d['rediscovery']).sum())}) | "
              f"rediscovery {int(d['rediscovery'].sum())}", flush=True)

    # ---- consensus across targets ----
    keyed = {}
    for t, d in tables.items():
        dd = d.set_index("identifier")
        keyed[t] = dd
    ids = tables[C.TARGETS[0]]["identifier"]
    cons = pd.DataFrame({"identifier": ids.values})
    base = keyed[C.TARGETS[0]]
    cons["name"] = base.loc[cons["identifier"], "name"].values
    cons["smiles"] = base.loc[cons["identifier"], "smiles"].values
    cons["np_likeness"] = base.loc[cons["identifier"], "np_likeness"].values
    n_conf = np.zeros(len(cons), int)
    n_inad_act = np.zeros(len(cons), int)
    any_red = np.zeros(len(cons), bool)
    for t in C.TARGETS:
        dd = keyed[t].loc[cons["identifier"]]
        cons[f"{t}_pred"] = dd["pred"].values
        cons[f"{t}_lo"] = dd["lo"].values
        cons[f"{t}_sim"] = dd["sim_train"].values
        cons[f"{t}_inAD"] = dd["in_AD"].values
        n_conf += (dd["in_AD"].values & (dd["lo"].values >= C.ACTIVE_THRESHOLD))
        n_inad_act += (dd["in_AD"].values & (dd["pred"].values >= C.ACTIVE_THRESHOLD))
        any_red |= dd["rediscovery"].values
    cons["n_targets_confident"] = n_conf
    cons["n_targets_inAD_active"] = n_inad_act
    cons["rediscovery_any"] = any_red
    cons["mean_pred"] = cons[[f"{t}_pred" for t in C.TARGETS]].mean(axis=1)
    cons = cons[(cons["n_targets_confident"] >= 1) | (cons["n_targets_inAD_active"] >= 1)].sort_values(
        ["n_targets_confident", "n_targets_inAD_active", "mean_pred"], ascending=False)
    cons.to_csv(os.path.join(OUT_DIR, "consensus_hits.csv"), index=False)
    summary["consensus"] = {
        "any_target_confident": int((cons["n_targets_confident"] >= 1).sum()),
        "two_or_more_targets_confident": int((cons["n_targets_confident"] >= 2).sum()),
        "all_three_confident": int((cons["n_targets_confident"] == 3).sum()),
        "confident_and_novel": int(((cons["n_targets_confident"] >= 1) & ~cons["rediscovery_any"]).sum()),
        "any_target_inAD_active": int((cons["n_targets_inAD_active"] >= 1).sum()),
        "two_or_more_targets_inAD_active": int((cons["n_targets_inAD_active"] >= 2).sum()),
        "inAD_active_and_novel": int(((cons["n_targets_inAD_active"] >= 1) & ~cons["rediscovery_any"]).sum()),
    }
    summary["runtime_min"] = round((time.time() - t_all) / 60, 1)
    C.save_json(summary, os.path.join(C.METRICS_DIR, "screen_summary.json"))
    print("\n" + json.dumps(summary, indent=1), flush=True)
    return {"status": "full library screened", "n": int(valid.sum())}


if __name__ == "__main__":
    mp.freeze_support()
    main()
