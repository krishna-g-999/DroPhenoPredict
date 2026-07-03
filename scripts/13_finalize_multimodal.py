"""Finalize DroPhenoPredict v1.1 — multi-modal, per-trait-adaptive.

(1) Nested-CV modality selection (unbiased adaptive accuracy + selection counts)
(2) Deployed-modality LOO validation + jackknife conformal 90% intervals + perm p
(3) Serialized artifact (manifest json + per-model kernel/y npz) for the predictor

Outputs:
  models/trained/drophenopredict_v1.1.json
  models/trained/drophenopredict_v1.1_kernels.npz
  data/processed/multimodal_selection.csv
Run: ./.venv/Scripts/python.exe scripts/13_finalize_multimodal.py
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, expression, gblup, multimodal  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

VERSION = "v1.1"
ALPHA = 0.10
TRAITS = {
    "climbing":   (1543, "mn_NegGeotaxis_CTRL", "index"),
    "startle":    (2799, "StartleRes", "score"),
    "starvation": (2798, "StarvationRes", "hours"),
    "sleep":      (1483, "mn_SleepLength_Night", "minutes"),
}
OUTDIR = Path("models/trained")
PROC = Path("data/processed")


def main() -> None:
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna = {s: expression.load_rnaseq(s) for s in ("female", "male")}
    sex_map = {"F": "female", "M": "male"}

    manifest = {"model": "DroPhenoPredict", "version": VERSION,
                "method": "multi-modal per-trait-adaptive GBLUP (genotype GRM | RNA-seq expression), "
                          "modality chosen by nested CV; jackknife conformal intervals",
                "traits": {}}
    arrays = {}
    rows = []

    for trait, (pid, phen, unit) in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in ("F", "M"):
            if sex not in set(vals["sex"].unique()):
                continue
            y_series = dgrpool.line_means(vals, sex=sex, harmonize=True)
            common, kernels, yv = multimodal.build_candidates(
                Kg, g_iids, rna[sex_map[sex]], y_series)
            if len(common) < 110:
                continue
            key = f"{trait}_{sex}"

            fixed = multimodal.fixed_modality_r(kernels, yv)
            nested = multimodal.nested_modality_cv(kernels, yv, n_repeats=10)
            deployed = multimodal.deployed_modality(kernels, yv)
            Kdep = kernels[deployed]
            fit = gblup.reml(Kdep, yv)
            loo = multimodal.loo_predict(Kdep, yv, fit.delta)   # point-prediction mechanism
            loo_r = float(np.corrcoef(yv, loo)[0, 1])
            # HONEST intervals: k-fold split-conformal with verified coverage (NOT LOO residuals)
            conf = gblup.conformal_coverage(Kdep, yv, alpha=ALPHA, n_repeats=10)
            q = conf["mean_interval_width"] / 2.0
            # significance against the unbiased adaptive R2
            perm = gblup.permutation_test(Kdep, yv, nested["adaptive_R2"], n_perm=200)

            arrays[f"{key}__K"] = Kdep
            arrays[f"{key}__y"] = yv
            manifest["traits"][key] = {
                "trait": trait, "sex": sex, "unit": unit, "phenotype_id": pid,
                "n_lines": len(common), "lines": common,
                "deployed_modality": deployed, "delta": round(fit.delta, 5),
                "genomic_h2": round(gblup.reml(kernels["geno"], yv).h2, 3),
                "r_geno": round(fixed["geno"], 3), "r_exp": round(fixed["exp"], 3),
                "r_combined": round(fixed["combined"], 3),
                "adaptive_nested_r": round(nested["adaptive_r"], 3),
                "adaptive_nested_r_sd": round(nested["adaptive_r_sd"], 3),
                "selection_counts": nested["selections"],
                "loo_r_inpanel": round(loo_r, 3),
                "perm_p": round(perm["p_value"], 4),
                "predictable": bool(perm["p_value"] < 0.05),
                "interval_halfwidth_90": round(q, 3),
                "interval_coverage_90": round(conf["empirical_coverage"], 3),
                "y_mean": round(float(yv.mean()), 3), "y_sd": round(float(yv.std(ddof=1)), 3),
            }
            m = manifest["traits"][key]
            rows.append({"trait_sex": key, "n": m["n_lines"], "geno": m["r_geno"],
                         "exp": m["r_exp"], "combined": m["r_combined"],
                         "deployed": deployed, "nested_adaptive_r": m["adaptive_nested_r"],
                         "perm_p": m["perm_p"], "predictable": m["predictable"],
                         "cover90": m["interval_coverage_90"]})
            print(f"{key:13s} n={m['n_lines']:3d} geno={m['r_geno']:+.3f} exp={m['r_exp']:+.3f} "
                  f"comb={m['r_combined']:+.3f} -> {deployed:4s} nested_r={m['adaptive_nested_r']:+.3f} "
                  f"p={m['perm_p']:.3f} cover90={m['interval_coverage_90']:.2f} "
                  f"{'PRED' if m['predictable'] else 'n.s.'}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTDIR / f"drophenopredict_{VERSION}_kernels.npz", **arrays)
    (OUTDIR / f"drophenopredict_{VERSION}.json").write_text(json.dumps(manifest, indent=2))
    pd.DataFrame(rows).to_csv(PROC / "multimodal_selection.csv", index=False)
    npred = sum(t["predictable"] for t in manifest["traits"].values())
    print(f"\nSaved {VERSION}: {len(manifest['traits'])} models, {npred} predictable.")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
