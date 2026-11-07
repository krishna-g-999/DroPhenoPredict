"""Package the predictable DGRP drug/toxin traits as model drug-v1.0 (extension of v1.1).

Inclusion (fixed before building): treated or drug-effect traits from Findings 21 with
FDR q < 0.05 in the genotype analysis or the multi-modal analysis
(data/processed/chemical_predictability.csv, script 49).

Per trait, the deployed set-up is the passing analysis with the higher accuracy:
  - 'geno'  : GRM on all genotyped lines with the trait (accuracy = CV r, 10 repeats)
  - 'multi' : the screen's deployed kernel (geno / exp_F / exp_M) on lines with both
              data types (accuracy = nested adaptive r)
Choosing the better of two passing analyses adds a small optimistic bias to the
reported r; the p and q reported are those of the chosen analysis.

Intervals: k-fold split-conformal (gblup.conformal_coverage, 90%), coverage measured,
as for v1.1. Artefacts are written in the v1.1 multimodal format so
DroPhenoPredictor("drug-v1.0") loads them unchanged:
  models/trained/drophenopredict_drug-v1.0.json
  models/trained/drophenopredict_drug-v1.0_kernels.npz
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import expression, gblup, genotypes, multimodal  # noqa: E402

warnings.filterwarnings("ignore")
VERSION = "drug-v1.0"
ALPHA = 0.10
Q_MAX = 0.05
PAIRS_SCRIPT = Path(__file__).with_name("49_chemical_predictability.py")
# Measurements that are proportions or percentages (checked against the data ranges in
# Findings 21). Clipping an interval to the valid range cannot change its coverage.
BOUNDS = {
    "parathion_survival_48h_B": [0, 1],
    "lead_viability_B": [0, 1],
    "malathion_12h_mortality_M": [0, 1],
    "malathion_12h_mortality_F": [0, 1],
    "octanoic_acid_survival_M": [0, 1],
    "ethanol_dev_viability_B": [0, 100],
}


def load_pairs() -> dict:
    """Reuse the exact treated/control pairs defined in script 49."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("s49", PAIRS_SCRIPT)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod.PAIRS


def main() -> None:
    screen = pd.read_csv("data/processed/chemical_predictability.csv")
    cat = pd.read_csv("data/processed/chemical_traits_catalog.csv").set_index("trait_id")
    long = pd.read_csv("data/processed/chemical_traits_long.csv")
    wide = long.pivot(index="line", columns="trait_id", values="value")
    pairs = load_pairs()
    for eff, (t, c) in pairs.items():
        wide[eff] = wide[t] - wide[c]

    sel = screen[screen["kind"].isin(["treated", "effect"])
                 & ((screen["q_geno"] < Q_MAX) | (screen["q_adaptive"] < Q_MAX))]

    Kg_full, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna = {"F": expression.load_rnaseq("female"), "M": expression.load_rnaseq("male")}

    manifest = {"model": "DroPhenoPredict", "version": VERSION, "extends": "v1.1",
                "method": "per-trait GBLUP (genotype GRM or RNA-seq expression kernel); "
                          "k-fold split-conformal 90% intervals",
                "inclusion": f"treated/effect traits with FDR q < {Q_MAX} in Findings 21",
                "traits": {}}
    arrays = {}
    for _, s in sel.iterrows():
        tid = s["trait_id"]
        geno_ok = s["q_geno"] < Q_MAX
        multi_ok = s["q_adaptive"] < Q_MAX
        use_multi = multi_ok and (not geno_ok or s["r_adaptive"] > s["r_geno"])
        y_s = wide[tid].dropna()
        lines = [l for l in g_iids if l in y_s.index]
        sex = tid.rsplit("_", 1)[1]
        if use_multi:
            sexes = ["F", "M"] if sex == "B" else [sex]
            lines = [l for l in lines if all(l in rna[x].columns for x in sexes)]
            modality = s["deployed"]
            r, p, q = s["r_adaptive"], s["p_adaptive"], s["q_adaptive"]
        else:
            modality = "geno"
            r, p, q = s["r_geno"], s["p_geno"], s["q_geno"]
        if modality == "geno":
            gi = [g_iids.index(l) for l in lines]
            K = Kg_full[np.ix_(gi, gi)]; K = K / np.mean(np.diag(K))
        else:
            K = multimodal.expression_kernel(rna[modality.split("_")[1]], lines)
        y = y_s.loc[lines].to_numpy(float)
        fit = gblup.reml(K, y)
        conf = gblup.conformal_coverage(K, y, alpha=ALPHA, n_repeats=10)
        src_trait = pairs[tid][0] if tid in pairs else tid
        info = cat.loc[src_trait]
        arrays[f"{tid}__K"] = K
        arrays[f"{tid}__y"] = y
        manifest["traits"][tid] = {
            "trait": tid.rsplit("_", 1)[0], "sex": sex, "kind": s["kind"],
            "compound": info["compound"], "higher_is": info["higher_is"],
            "unit": "treated - control (source units)" if tid in pairs else "source units",
            "source": info["source"],
            "control_trait": pairs[tid][1] if tid in pairs else None,
            "n_lines": len(lines), "lines": lines, "analysis": "multi" if use_multi else "geno",
            "deployed_modality": modality, "delta": round(fit.delta, 5),
            "genomic_h2": round(float(s["h2_geno"]), 3),
            "cv_r": round(float(r), 3), "perm_p": round(float(p), 4), "fdr_q": round(float(q), 4),
            "predictable": True,
            "interval_halfwidth_90": round(conf["mean_interval_width"] / 2, 4),
            "interval_coverage_90": round(conf["empirical_coverage"], 3),
            "y_mean": round(float(y.mean()), 4), "y_sd": round(float(y.std(ddof=1)), 4),
            # valid range of the measurement; the predictor clips intervals to it
            "bounds": BOUNDS.get(tid),
        }
        m = manifest["traits"][tid]
        # How informative are the intervals? Width of each line's (clipped) interval, as a
        # share of the observed or valid range, using the same leave-one-out point the API serves.
        pts = multimodal.loo_predict(K, y, fit.delta)
        hw = m["interval_halfwidth_90"]
        lo, hi = pts - hw, pts + hw
        b = m["bounds"] or [float(y.min()), float(y.max())]
        if m["bounds"]:
            lo, hi = np.maximum(lo, b[0]), np.minimum(hi, b[1])
        rel = float(np.median((hi - lo) / (b[1] - b[0])))
        m["interval_median_share_of_range"] = round(rel, 3)
        if rel >= 0.85:
            m["interval_note"] = ("interval usually spans most of the possible range: use the "
                                  "prediction to rank lines, not as a calibrated value")
        print(f"{tid:38s} {m['analysis']:5s} {modality:5s} n={len(lines):3d} r={r:.3f} q={q:.3f} "
              f"cover90={m['interval_coverage_90']:.3f} halfwidth={m['interval_halfwidth_90']:.3g} "
              f"(y sd {m['y_sd']:.3g})")

    out = Path("models/trained")
    np.savez_compressed(out / f"drophenopredict_{VERSION}_kernels.npz", **arrays)
    (out / f"drophenopredict_{VERSION}.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nSaved {VERSION}: {len(manifest['traits'])} trait models.")


if __name__ == "__main__":
    main()
