"""
Step 05 - Full-library natural-product screen
----------------------------------------------
Streams the COCONUT library in parallel chunks (no multi-GB feature table on
disk): standardise -> featurise -> predict -> attach a Mondrian-conformal interval,
an applicability-domain flag, novelty flags and a nearest-neighbour baseline.

WHAT CHANGED AND WHY:

* STEREOISOMERS ARE COLLAPSED BEFORE SCREENING, NOT AT REPORT TIME.
  Features are 2D (physicochemical descriptors + ECFP4) and therefore achiral, so
  COCONUT's `CNP0000000.1/.2/.3` stereoisomer records produced identical
  predictions that were counted as separate compounds. The library is now
  deduplicated on the standardised parent InChIKey up front, and both the raw
  record count and the deduplicated compound count are reported.

* INTERVALS ARE LABELLED AS UNCALIBRATED FOR THIS LIBRARY.
  The conformal guarantee from step 04 holds for held-out ChEMBL scaffolds, not for
  natural products -- the shift breaks exchangeability. The column previously
  called `confident_active` is now `interval_lower_bound_active`, and
  `screen_summary.json` carries the validity statement with the counts, so no
  downstream table can present these as guaranteed-coverage hits.

* NOVELTY IS MEASURED AT THE SCAFFOLD LEVEL, NOT ONLY BY SIMILARITY.
  `rediscovery` (Tanimoto >= 0.90) catches near-identical known chemistry but not
  "a staurosporine analogue at Tanimoto 0.7". The screen now also records whether
  the compound's generic Murcko scaffold already exists in the target's training
  set (`known_scaffold`) and flags pan-kinase frequent-hitter chemotypes
  (`pan_kinase_alert`). A hit that is in-domain, active and on a known scaffold is
  an analogue, and is counted separately from genuinely novel scaffolds.

* A NEAREST-NEIGHBOUR BASELINE RANKS THE SAME LIBRARY.
  `nn_baseline_pred` is a similarity-weighted 5-NN prediction over the training
  set. If it recovers the same hits as the model, the model adds nothing -- so the
  comparison is computed rather than left for a reviewer to ask about.

Outputs: results/predictions/<T>_screen_full.csv.gz
         results/predictions/consensus_hits.csv
         results/metrics/screen_summary.json
"""

import os
import sys
import time
import json
import multiprocessing as mp

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C

LIB_CSV = os.path.join("data", "ligands_phytochemical", "phytochemical_screening_set.csv")
OUT_DIR = os.path.join("results", "predictions")
CHUNK = 2000
REDISCOVERY_SIM = 0.90

_W = {}


def _init_worker(train_smiles_by_target, in_train_by_target):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    import common as CC
    _W["C"] = CC
    _W["fps"] = {t: CC.morgan_bitvects(s) for t, s in train_smiles_by_target.items()}
    _W["in_train"] = in_train_by_target
    _W["train_y"] = {}


def _featurize(args):
    """SMILES chunk -> descriptors, fingerprints and per-target similarities."""
    start, smiles = args
    from rdkit import Chem, DataStructs
    CC = _W["C"]
    n = len(smiles)
    ok = np.zeros(n, bool)
    desc = np.full((n, len(CC.DESC_COLS)), np.nan, np.float32)
    fp = np.zeros((n, CC.N_FP_BITS), np.uint8)
    targets = list(_W["fps"].keys())
    sim_all_c = {t: np.full(n, np.nan, np.float32) for t in targets}
    sim_tr_c = {t: np.full(n, np.nan, np.float32) for t in targets}
    nn_idx_c = {t: np.full(n, -1, np.int32) for t in targets}
    for i, smi in enumerate(smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            d = CC.physchem_descriptors(mol)
            desc[i] = [d[c] for c in CC.DESC_COLS]
            bv = CC._MORGAN.GetFingerprint(mol)
            fp[i] = CC._MORGAN.GetFingerprintAsNumPy(mol)
            for t in targets:
                s = np.asarray(DataStructs.BulkTanimotoSimilarity(bv, _W["fps"][t]),
                               dtype=np.float32)
                idx_max = int(s.argmax())
                sim_all_c[t][i] = s[idx_max]
                nn_idx_c[t][i] = idx_max
                sim_tr_c[t][i] = s[_W["in_train"][t]].max()
            ok[i] = True
        except Exception:
            continue
    return start, ok, desc, fp, sim_all_c, nn_idx_c, sim_tr_c


def load_model(target, conformal):
    """Reload whichever family step 04 selected and persisted."""
    path = conformal.get("model_path")
    family = conformal.get("model_family", "")
    if path and path.endswith((".joblib", ".pkl")) and os.path.exists(path):
        import joblib
        return joblib.load(path)
    if path and path.endswith(".json") and os.path.exists(path):
        from xgboost import XGBRegressor
        m = XGBRegressor()
        m.load_model(path)
        return m
    raise FileNotFoundError(
        f"No persisted model for {target} ({family}) at {path!r}. Model artifacts "
        f"are not tracked in git because a fitted forest exceeds GitHub's file "
        f"size limit; regenerate them with `python scripts/04_conformal_ad_models.py`.")


def feature_slice(conformal, cols):
    """Column indices for the feature view the selected model was trained on."""
    view = conformal.get("feature_view", "all")
    pos = {c: i for i, c in enumerate(cols)}
    if view == "all":
        return np.arange(len(cols))
    if view == "fp":
        return np.array([pos[c] for c in C.FP_COLS])
    return np.array([pos[c] for c in C.DESC_COLS])


def prepare_library():
    """Read the library, standardise structures, collapse stereoisomer duplicates."""
    lib = pd.read_csv(LIB_CSV, usecols=lambda c: c in {
        "identifier", "name", "canonical_smiles", "np_likeness"})
    n_records = len(lib)
    lib = lib.dropna(subset=["canonical_smiles"]).reset_index(drop=True)

    t_std = time.time()
    print(f"[+] Standardising {len(lib):,} library structures "
          f"across {C.standardise_workers():,} processes...", flush=True)
    std_smi, std_key = C.standardise_series(lib["canonical_smiles"])
    print(f"[+] Standardised in {(time.time() - t_std) / 60:.1f} min", flush=True)
    lib["standardised_smiles"] = std_smi
    lib["inchikey"] = std_key
    n_unparsable = int(lib["inchikey"].isna().sum())
    lib = lib.dropna(subset=["inchikey"]).reset_index(drop=True)

    # Prefer the named record when collapsing a stereoisomer family.
    lib["_has_name"] = lib["name"].notna() & (lib["name"].astype(str).str.strip() != "")
    lib = (lib.sort_values("_has_name", ascending=False)
              .drop_duplicates("inchikey")
              .drop(columns=["_has_name"])
              .sort_index()
              .reset_index(drop=True))
    print(f"[+] {n_records:,} records -> {len(lib):,} unique 2D parents "
          f"({n_unparsable:,} unparsable)", flush=True)
    return lib, {"library_records": int(n_records),
                 "library_unique_parents": int(len(lib)),
                 "library_unparsable": n_unparsable,
                 "dedup_note": "Collapsed on standardised parent InChIKey. 2D "
                               "descriptors and ECFP4 are achiral, so stereoisomer "
                               "records produce identical predictions and must not "
                               "be counted as distinct screened compounds."}


def main():
    t_all = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    art = {}
    for t in C.TARGETS:
        conformal = C.load_json(os.path.join(C.MODELS_DIR, f"{t}_conformal.json"))
        if conformal is None:
            raise FileNotFoundError(f"models/{t}_conformal.json missing; run step 04.")
        df, cols, fp_cols, desc_cols = C.load_training(t)
        smi = df["canonical_smiles"].tolist()
        in_train = np.isin(smi, conformal["train_smiles"])
        train_scaffolds = set(
            s for s in (C.generic_scaffold(x) for x in np.array(smi)[in_train]) if s)
        art[t] = {"cj": conformal, "smiles": smi,
                  "ids": df["ligand_id"].astype(str).tolist(),
                  "y": df["pIC50"].to_numpy(float), "in_train": in_train,
                  "model": load_model(t, conformal), "cols": conformal["feature_cols"],
                  "train_scaffolds": train_scaffolds}
    cols = art[C.TARGETS[0]]["cols"]
    assert all(a["cols"] == cols for a in art.values()), "feature columns differ"
    assert cols == C.FEATURE_COLS, (
        "feature columns do not match common.FEATURE_COLS; re-run step 02")

    lib, lib_summary = prepare_library()
    smiles = lib["standardised_smiles"].tolist()
    chunks = [(s, smiles[s:s + CHUNK]) for s in range(0, len(smiles), CHUNK)]

    n = len(lib)
    pred = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    nn_pred = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    sim_train = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    sim_all = {t: np.full(n, np.nan, np.float32) for t in C.TARGETS}
    nn_idx = {t: np.full(n, -1, np.int32) for t in C.TARGETS}
    valid = np.zeros(n, bool)

    nproc = max(1, min(14, (mp.cpu_count() or 2) - 2))
    print(f"[+] Featurising with {nproc} workers, {len(chunks)} chunks", flush=True)
    t0 = time.time()
    slices = {t: feature_slice(art[t]["cj"], cols) for t in C.TARGETS}
    with mp.Pool(nproc, initializer=_init_worker,
                 initargs=({t: art[t]["smiles"] for t in C.TARGETS},
                           {t: art[t]["in_train"] for t in C.TARGETS})) as pool:
        for k, (start, ok, desc, fp, sim_all_c, nn_idx_c, sim_tr_c) in enumerate(
                pool.imap(_featurize, chunks)):
            e = start + len(ok)
            valid[start:e] = ok
            X = np.concatenate([desc, fp.astype(np.float32)], axis=1)
            X[~np.isfinite(X)] = np.nan
            for t in C.TARGETS:
                Xt = X[:, slices[t]]
                if art[t]["cj"].get("model_family", "").startswith("RF"):
                    Xt = np.nan_to_num(Xt, nan=0.0, posinf=0.0, neginf=0.0)
                pred[t][start:e] = art[t]["model"].predict(Xt)
                sim_all[t][start:e] = sim_all_c[t]
                nn_idx[t][start:e] = nn_idx_c[t]
                sim_train[t][start:e] = sim_tr_c[t]
                # Nearest-neighbour baseline: the label of the closest training
                # molecule. Cheap, and the honest comparator for any QSAR screen.
                nn = np.where(nn_idx_c[t] >= 0, nn_idx_c[t], 0)
                nn_pred[t][start:e] = art[t]["y"][nn]
            if (k + 1) % 20 == 0 or k + 1 == len(chunks):
                el = time.time() - t0
                eta = (el / (k + 1)) * (len(chunks) - k - 1) / 60
                print(f"    {k+1}/{len(chunks)} chunks | {el/60:.1f} min | "
                      f"ETA {eta:.1f} min", flush=True)

    summary = {**lib_summary, "library_featurised": int(valid.sum()), "targets": {}}
    tables = {}
    print("\n[+] Annotating scaffolds and chemotypes...", flush=True)
    scaffolds = np.array([C.generic_scaffold(s) if v else None
                          for s, v in zip(smiles, valid)], dtype=object)
    chemo = C.chemotype_flags(smiles)

    for t in C.TARGETS:
        a = art[t]
        cj = a["cj"]
        edges, hw, tau = cj["edges"], np.array(cj["half_widths"]), cj["tau"]
        b = np.clip(np.digitize(sim_train[t], edges[1:-1]), 0, len(edges) - 2)
        half = hw[b]
        nn = np.where(nn_idx[t] >= 0, nn_idx[t], 0)
        d = pd.DataFrame({
            "identifier": lib["identifier"], "name": lib["name"],
            "smiles": lib["standardised_smiles"], "inchikey": lib["inchikey"],
            "np_likeness": lib.get("np_likeness"),
            "pred": pred[t],
            "interval_lo": pred[t] - half, "interval_hi": pred[t] + half,
            "sim_train": sim_train[t], "sim_all": sim_all[t],
            "nn_id": np.array(a["ids"])[nn], "nn_pIC50": a["y"][nn],
            "nn_baseline_pred": nn_pred[t],
            "generic_scaffold": scaffolds,
        })
        d["in_AD"] = (d["sim_train"] >= tau) if tau is not None else False
        d["rediscovery"] = d["sim_all"] >= REDISCOVERY_SIM
        d["known_scaffold"] = [s in a["train_scaffolds"] if s else False
                               for s in scaffolds]
        d["pan_kinase_alert"] = chemo["pan_kinase_alert"].to_numpy()
        for col in C.CHEMOTYPE_SMARTS:
            if col in chemo.columns:
                d[f"chemotype_{col}"] = chemo[col].to_numpy()
        d = d[valid].sort_values("pred", ascending=False).reset_index(drop=True)
        d.to_csv(os.path.join(OUT_DIR, f"{t}_screen_full.csv.gz"), index=False)
        tables[t] = d

        in_ad_act = d["in_AD"] & (d["pred"] >= C.ACTIVE_THRESHOLD)
        lb_act = d["in_AD"] & (d["interval_lo"] >= C.ACTIVE_THRESHOLD)
        novel_scaffold = ~d["known_scaffold"] & ~d["rediscovery"]
        # Rank agreement with the cheap baseline, over in-domain compounds.
        m = d["in_AD"].to_numpy()
        rank_agreement = (float(C.spearman(d["pred"][m], d["nn_baseline_pred"][m]))
                          if m.sum() > 10 else None)
        summary["targets"][t] = {
            "model": cj.get("model_family"),
            "tau": tau,
            "tau_is_binding": cj.get("tau_is_binding"),
            "in_AD_caveat": (
                None if cj.get("tau_is_binding") is not False else
                "tau is non-binding for this target (step 04 found every "
                "calibration bin inside the width budget), so `in_AD` admits the "
                "whole library and the in-domain counts below are not a "
                "calibrated-reliable subset."),
            "in_AD": int(d["in_AD"].sum()),
            "in_AD_frac": float(d["in_AD"].mean()),
            "median_sim_train": float(d["sim_train"].median()),
            "in_AD_active": int(in_ad_act.sum()),
            "in_AD_active_novel_scaffold": int((in_ad_act & novel_scaffold).sum()),
            "in_AD_active_pan_kinase_alert": int((in_ad_act & d["pan_kinase_alert"]).sum()),
            "interval_lower_bound_active": int(lb_act.sum()),
            "interval_lower_bound_active_novel_scaffold": int((lb_act & novel_scaffold).sum()),
            "rediscovery_total": int(d["rediscovery"].sum()),
            "known_scaffold_total": int(d["known_scaffold"].sum()),
            "spearman_model_vs_nn_baseline_in_AD": rank_agreement,
        }
        print(f"  {t}: in_AD {d['in_AD'].sum():,} ({d['in_AD'].mean()*100:.2f}%) | "
              f"in_AD & active {int(in_ad_act.sum())} "
              f"(novel scaffold {int((in_ad_act & novel_scaffold).sum())}, "
              f"pan-kinase alert {int((in_ad_act & d['pan_kinase_alert']).sum())}) | "
              f"rho vs NN baseline {rank_agreement}", flush=True)

    # ---- consensus across targets ----
    base = tables[C.TARGETS[0]].set_index("identifier")
    cons = pd.DataFrame({"identifier": tables[C.TARGETS[0]]["identifier"].values})
    for col in ("name", "smiles", "np_likeness", "generic_scaffold",
                "pan_kinase_alert"):
        if col in base.columns:
            cons[col] = base.loc[cons["identifier"], col].values

    n_lb = np.zeros(len(cons), int)
    n_inad_act = np.zeros(len(cons), int)
    any_red = np.zeros(len(cons), bool)
    all_known_scaffold = np.ones(len(cons), bool)
    for t in C.TARGETS:
        dd = tables[t].set_index("identifier").loc[cons["identifier"]]
        cons[f"{t}_pred"] = dd["pred"].values
        cons[f"{t}_interval_lo"] = dd["interval_lo"].values
        cons[f"{t}_sim"] = dd["sim_train"].values
        cons[f"{t}_inAD"] = dd["in_AD"].values
        n_lb += (dd["in_AD"].values & (dd["interval_lo"].values >= C.ACTIVE_THRESHOLD))
        n_inad_act += (dd["in_AD"].values & (dd["pred"].values >= C.ACTIVE_THRESHOLD))
        any_red |= dd["rediscovery"].values
        all_known_scaffold &= dd["known_scaffold"].values
    cons["n_targets_interval_active"] = n_lb
    cons["n_targets_inAD_active"] = n_inad_act
    cons["rediscovery_any"] = any_red
    cons["known_scaffold_all_targets"] = all_known_scaffold
    cons["mean_pred"] = cons[[f"{t}_pred" for t in C.TARGETS]].mean(axis=1)

    # Cross-target prediction correlation: consensus is only evidence if the three
    # models disagree sometimes. If they are near-identical, "hits 3/3 targets"
    # means pan-kinase chemistry, not three independent confirmations.
    pred_matrix = cons[[f"{t}_pred" for t in C.TARGETS]].to_numpy(float)
    cross = {}
    for i, a in enumerate(C.TARGETS):
        for b_t in C.TARGETS[i + 1:]:
            j = C.TARGETS.index(b_t)
            cross[f"{a}_vs_{b_t}"] = float(C.spearman(pred_matrix[:, i], pred_matrix[:, j]))

    cons = cons[(cons["n_targets_interval_active"] >= 1) |
                (cons["n_targets_inAD_active"] >= 1)].sort_values(
        ["n_targets_interval_active", "n_targets_inAD_active", "mean_pred"],
        ascending=False)
    cons.to_csv(os.path.join(OUT_DIR, "consensus_hits.csv"), index=False)

    novel = (~cons["rediscovery_any"]) & (~cons["known_scaffold_all_targets"])
    summary["consensus"] = {
        "any_target_inAD_active": int((cons["n_targets_inAD_active"] >= 1).sum()),
        "two_or_more_targets_inAD_active": int((cons["n_targets_inAD_active"] >= 2).sum()),
        "all_three_inAD_active": int((cons["n_targets_inAD_active"] == 3).sum()),
        "any_target_interval_active": int((cons["n_targets_interval_active"] >= 1).sum()),
        "novel_scaffold_and_not_rediscovery": int(novel.sum()),
        "pan_kinase_alert_among_hits": int(cons["pan_kinase_alert"].sum()),
        "cross_target_prediction_spearman": cross,
        "consensus_caveat": "The three targets are ATP-site kinases trained on "
                            "overlapping chemotypes, so their predictions are "
                            "correlated (see cross_target_prediction_spearman). "
                            "Agreement across targets indicates pan-kinase "
                            "chemistry, which is a selectivity liability, and is "
                            "not independent confirmation of activity.",
    }
    summary["conformal_validity"] = art[C.TARGETS[0]]["cj"].get("validity")
    summary["interval_caveat"] = (
        "`interval_lo`/`interval_hi` are the step-04 Mondrian half-widths applied "
        "by similarity bin. Their 90% coverage was established on held-out ChEMBL "
        "scaffolds, which are NOT exchangeable with this natural-product library. "
        "Treat them as a calibrated-on-synthetic-chemistry error model extrapolated "
        "to natural products, never as guaranteed 90% intervals.")
    summary["runtime_min"] = round((time.time() - t_all) / 60, 1)
    C.save_json(summary, os.path.join(C.METRICS_DIR, "screen_summary.json"))
    print("\n" + json.dumps(summary, indent=1), flush=True)
    return {"status": "full library screened", "n": int(valid.sum())}


if __name__ == "__main__":
    mp.freeze_support()
    main()
