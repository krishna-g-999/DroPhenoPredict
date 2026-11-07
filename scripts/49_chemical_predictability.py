"""Predictability screen for DGRP drug/toxin response traits (Findings 21).

Traits: data/processed/chemical_traits_long.csv (script 48), plus derived
drug-EFFECT traits = treated - control for each paired design (computed per line,
only for lines measured in both arms).

Analysis, using the same verified engine as the v1.1 predictor:
  A. Genotype GBLUP on every trait: 5-fold CV x 10 repeats (mean r, R2),
     200-permutation test (gblup.permutation_test), Benjamini-Hochberg FDR across traits.
  B. Multi-modal: genotype vs RNA-seq expression (GSE117850) by nested CV
     (multimodal.nested_modality_cv), on lines with both data types. Sex-specific
     traits use the same-sex expression; traits with sexes pooled ('B') offer both
     female and male expression as candidates. Significance as in v1.1: permutation
     of the deployed kernel against the nested adaptive R2.
Seed 42. Output: data/processed/chemical_predictability.csv
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import expression, gblup, genotypes, multimodal  # noqa: E402

warnings.filterwarnings("ignore")
MIN_LINES = 80
PAIRS = {  # effect trait -> (treated, control)
    "climbing_MSB_effect_F": ("climbing_MSB_F", "climbing_MSB_ctrl_F"),
    "climbing_MSB_effect_M": ("climbing_MSB_M", "climbing_MSB_ctrl_M"),
    "startle_MSB_effect_F": ("startle_MSB_F", "startle_MSB_ctrl_F"),
    "startle_MSB_effect_M": ("startle_MSB_M", "startle_MSB_ctrl_M"),
    "MeHg10_effect_B": ("MeHg10_survival_B", "MeHg_ctrl_survival_B"),
    "caffeine_on_MeHg10_effect_B": ("MeHg10_caffeine_survival_B", "MeHg10_survival_B"),
    "lead_dev_time_effect_F": ("lead_dev_time_F", "lead_dev_time_ctrl_F"),
    "lead_dev_time_effect_M": ("lead_dev_time_M", "lead_dev_time_ctrl_M"),
    "lead_viability_effect_B": ("lead_viability_B", "lead_viability_ctrl_B"),
    "lead_activity_effect_F": ("lead_activity_F", "lead_activity_ctrl_F"),
    "lead_activity_effect_M": ("lead_activity_M", "lead_activity_ctrl_M"),
    "methylphenidate_food_intake_effect_M": ("methylphenidate_food_intake_M", "methylphenidate_food_intake_ctrl_M"),
    "methylphenidate_activity_effect_M": ("methylphenidate_activity_M", "methylphenidate_activity_ctrl_M"),
    "ethanol_dev_viability_effect_B": ("ethanol_dev_viability_B", "ethanol_dev_viability_ctrl_B"),
    "ethanol_dev_time_effect_B": ("ethanol_dev_time_B", "ethanol_dev_time_ctrl_B"),
    "paraquat_climbing_effect_B": ("paraquat_climbing_B", "paraquat_climbing_ctrl_B"),
    "rapamycin_pupation_effect_B": ("rapamycin_pupation_time_B", "rapamycin_pupation_time_ctrl_B"),
}


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main() -> None:
    long = pd.read_csv("data/processed/chemical_traits_long.csv")
    cat = pd.read_csv("data/processed/chemical_traits_catalog.csv").set_index("trait_id")
    wide = long.pivot(index="line", columns="trait_id", values="value")
    for eff, (t, c) in PAIRS.items():
        wide[eff] = wide[t] - wide[c]

    Kg_full, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna = {"F": expression.load_rnaseq("female"), "M": expression.load_rnaseq("male")}

    rows = []
    t0 = time.time()
    for trait in wide.columns:
        y_s = wide[trait].dropna()
        lines = [l for l in g_iids if l in y_s.index]
        if len(lines) < MIN_LINES:
            continue
        sex = trait.rsplit("_", 1)[1]
        kind = "effect" if trait in PAIRS else cat.loc[trait, "role"]
        compound = cat.loc[PAIRS[trait][0] if trait in PAIRS else trait, "compound"]

        # A. genotype on all genotyped lines
        gi = [g_iids.index(l) for l in lines]
        Kg = Kg_full[np.ix_(gi, gi)]; Kg = Kg / np.mean(np.diag(Kg))
        y = y_s.loc[lines].to_numpy(float)
        cv = gblup.cross_validate(Kg, y, n_repeats=10, seed=42)
        r2 = float(cv.r2_per_repeat.mean())
        perm = gblup.permutation_test(Kg, y, r2, n_perm=200)
        h2 = gblup.reml(Kg, y).h2
        row = {"trait_id": trait, "compound": compound, "kind": kind, "sex": sex,
               "n_geno": len(lines), "h2_geno": round(h2, 3),
               "r_geno": round(float(cv.r_per_repeat.mean()), 3), "R2_geno": round(r2, 3),
               "p_geno": perm["p_value"]}

        # B. multi-modal on lines with expression
        sexes = ["F", "M"] if sex == "B" else [sex]
        common = [l for l in lines if all(l in rna[s].columns for s in sexes)]
        if len(common) >= MIN_LINES:
            gi2 = [g_iids.index(l) for l in common]
            Kg2 = Kg_full[np.ix_(gi2, gi2)]; Kg2 = Kg2 / np.mean(np.diag(Kg2))
            kernels = {"geno": Kg2}
            for s in sexes:
                kernels[f"exp_{s}"] = multimodal.expression_kernel(rna[s], common)
            y2 = y_s.loc[common].to_numpy(float)
            fixed = {k: float(gblup.cross_validate(K, y2, n_repeats=10, seed=42).r_per_repeat.mean())
                     for k, K in kernels.items()}
            nested = multimodal.nested_modality_cv(kernels, y2, n_repeats=5)
            dep = max(fixed, key=fixed.get)
            p2 = gblup.permutation_test(kernels[dep], y2, nested["adaptive_R2"], n_perm=200)["p_value"]
            row.update({"n_multi": len(common),
                        "r_exp_F": round(fixed.get("exp_F", np.nan), 3),
                        "r_exp_M": round(fixed.get("exp_M", np.nan), 3),
                        "deployed": dep, "r_adaptive": round(nested["adaptive_r"], 3),
                        "R2_adaptive": round(nested["adaptive_R2"], 3), "p_adaptive": p2})
        rows.append(row)
        print(f"[{time.time() - t0:6.0f}s] {trait:38s} n={len(lines):3d} h2={h2:.2f} "
              f"r_geno={row['r_geno']:+.3f} p={row['p_geno']:.3f} | "
              f"adaptive={row.get('r_adaptive', np.nan):+.3f} ({row.get('deployed', '-')}) "
              f"p={row.get('p_adaptive', np.nan):.3f}", flush=True)

    out = pd.DataFrame(rows)
    out["q_geno"] = bh(out["p_geno"].to_numpy())
    has = out["p_adaptive"].notna()
    out.loc[has, "q_adaptive"] = bh(out.loc[has, "p_adaptive"].to_numpy())
    out["best_r"] = out[["r_geno", "r_adaptive"]].max(axis=1)
    out = out.sort_values("best_r", ascending=False)
    out.to_csv("data/processed/chemical_predictability.csv", index=False)
    print(f"\n{len(out)} traits screened. Genotype FDR<0.05: {(out.q_geno < 0.05).sum()}; "
          f"adaptive FDR<0.05: {(out.q_adaptive < 0.05).sum()}")


if __name__ == "__main__":
    main()
