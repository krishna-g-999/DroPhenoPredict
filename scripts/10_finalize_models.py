"""Finalize the DroPhenoPredict model: per-trait GBLUP with leave-one-out
validation, jackknife conformal intervals, and predictability classification.

Decision (from Findings 03-05 + covariate test): the validated model is
genotype-only GBLUP on the GRM. Transcriptome, marker-selection, and
inversion/Wolbachia covariates were each tested and did NOT improve out-of-sample
accuracy, so the final model omits them (documented, not assumed).

For every trait x sex with coverage it stores: REML delta/h2, leave-one-out
predictions, LOO predictive r/R2, a jackknife conformal half-width (90% PI),
permutation p, and a predictability flag (perm p < 0.05).

Outputs:
  models/trained/drophenopredict_v1.0_targets.npz   (iids + per-trait y vectors)
  models/trained/drophenopredict_v1.0.json          (per-trait params + metrics)
Run: ./.venv/Scripts/python.exe scripts/10_finalize_models.py
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

MODEL_VERSION = "v1.0"
ALPHA = 0.10                     # 90% prediction intervals
TRAITS = {
    "climbing":   (1543, "mn_NegGeotaxis_CTRL", "index"),
    "startle":    (2799, "StartleRes", "score"),
    "starvation": (2798, "StarvationRes", "hours"),
    "sleep":      (1483, "mn_SleepLength_Night", "minutes"),
    "lifespan":   (1315, "mn_Longevity", "days"),
}
OUTDIR = Path("models/trained")


def loo_predictions(K, y, delta):
    """Leave-one-out GBLUP predictions (fixed variance ratio delta)."""
    n = len(y)
    pred = np.empty(n)
    idx = np.arange(n)
    for i in range(n):
        train = idx[idx != i]
        pred[i] = gblup.gblup_predict(K, y, train, np.array([i]), delta=delta)[0]
    return pred


def main() -> None:
    K_full, iids, gmeta = genotypes.load_grm("data/processed/grm.npz")
    n_master = len(iids)
    pos = {l: i for i, l in enumerate(iids)}

    manifest = {"model": "DroPhenoPredict", "version": MODEL_VERSION,
                "method": "genotype-only GBLUP on DGRP GRM (1.49M autosomal SNPs)",
                "interval": f"{int((1-ALPHA)*100)}% jackknife conformal",
                "grm": "data/processed/grm.npz", "traits": {}}
    target_arrays = {"iids": np.array(iids)}

    for trait, (pid, phen_name, unit) in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in sorted(s for s in vals["sex"].unique() if s in ("F", "M")):
            y_series = dgrpool.line_means(vals, sex=sex, harmonize=True)
            common = [l for l in iids if l in set(y_series.index)]
            if len(common) < 110:
                continue
            sub = [pos[l] for l in common]
            K = K_full[np.ix_(sub, sub)]
            yv = y_series.loc[common].to_numpy(float)

            fit = gblup.reml(K, yv)
            loo = loo_predictions(K, yv, fit.delta)
            resid = np.abs(yv - loo)
            loo_r = float(np.corrcoef(yv, loo)[0, 1])
            loo_R2 = float(1 - np.sum((yv - loo) ** 2) / np.sum((yv - yv.mean()) ** 2))
            q = float(np.quantile(resid, min(1.0, np.ceil((len(resid) + 1) * (1 - ALPHA)) / len(resid))))
            perm = gblup.permutation_test(K, yv, loo_R2, n_perm=200)
            key = f"{trait}_{sex}"

            # store y aligned to master iids (NaN where missing)
            yfull = np.full(n_master, np.nan)
            for l, v in zip(common, yv):
                yfull[pos[l]] = v
            target_arrays[key] = yfull

            manifest["traits"][key] = {
                "trait": trait, "sex": sex, "phenotype_id": pid,
                "phenotype_name": phen_name, "unit": unit, "n_lines": len(common),
                "genomic_h2": round(fit.h2, 3), "h2_at_bound": bool(fit.at_bound),
                "delta": round(fit.delta, 5),
                "loo_r": round(loo_r, 3), "loo_R2": round(loo_R2, 4),
                "perm_p": round(perm["p_value"], 4),
                "predictable": bool(perm["p_value"] < 0.05),
                "interval_halfwidth_90": round(q, 3),
                "y_mean": round(float(yv.mean()), 3), "y_sd": round(float(yv.std(ddof=1)), 3),
            }
            m = manifest["traits"][key]
            print(f"{key:14s} n={m['n_lines']:3d} h2={m['genomic_h2']:.2f} "
                  f"LOO_r={m['loo_r']:+.3f} p={m['perm_p']:.3f} "
                  f"{'PREDICTABLE' if m['predictable'] else 'not-predictive':14s} "
                  f"PI90=±{m['interval_halfwidth_90']:.2f}{unit[:3]}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTDIR / f"drophenopredict_{MODEL_VERSION}_targets.npz", **target_arrays)
    (OUTDIR / f"drophenopredict_{MODEL_VERSION}.json").write_text(json.dumps(manifest, indent=2))
    n_pred = sum(t["predictable"] for t in manifest["traits"].values())
    print(f"\nSaved model {MODEL_VERSION}: {len(manifest['traits'])} trait-sex models, "
          f"{n_pred} predictable (perm p<0.05).")
    print(f"  -> {OUTDIR/f'drophenopredict_{MODEL_VERSION}.json'}")


if __name__ == "__main__":
    main()
