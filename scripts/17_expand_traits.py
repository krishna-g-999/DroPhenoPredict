"""Feature expansion — discover additional predictable DGRP traits.

Scans DGRPool for biologically-distinct, high-coverage numeric-continuous traits
beyond the current 4 families, and runs each through the validated pipeline
(genotype vs RNA-seq expression CV r + permutation test) to find new predictable
traits to add to the resource. Verified, no assumptions.

Output: data/processed/expanded_traits.csv
"""
from __future__ import annotations

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

# New biologically-distinct families to probe (keyword -> match in phenotype name)
NEW_FAMILIES = {
    "oxidative_paraquat": ["paraquat"],
    "oxidative_msb": ["msb", "menadione"],
    "chill_coma": ["chill", "coma", "ccrt"],
    "activity": ["activity", "waking"],
    "food_intake": ["cafe", "food intake", "feeding"],
    "body_weight": ["weight", "mass"],
    "triglyceride": ["triglyceride", "tag "],
    "glucose": ["glucose"],
    "ethanol": ["ethanol", "etoh"],
    "lifespan2": ["longevity", "lifespan"],
}
MIN_COVER = 150
PROC = Path("data/processed")


def main() -> None:
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna = {s: expression.load_rnaseq(s) for s in ("female", "male")}
    sex_map = {"F": "female", "M": "male"}
    gset = set(g_iids)

    cat = dgrpool.fetch_phenotype_catalog()
    cont = dgrpool.numeric_continuous_phenotypes(cat)
    rows = []

    for fam, kws in NEW_FAMILIES.items():
        # candidate phenotype ids whose name matches a keyword
        mask = cont["name"].str.lower().apply(lambda s: any(k in s for k in kws))
        cands = cont[mask]
        if cands.empty:
            continue
        # pick best-covered candidate (overlap with genotype + rnaseq)
        best = None
        for pid, row in cands.iterrows():
            try:
                vals = dgrpool.fetch_phenotype_values(int(pid))
            except Exception:
                continue
            for sex in ("F", "M"):
                if sex not in set(vals["sex"].unique()):
                    continue
                y = dgrpool.line_means(vals, sex=sex, harmonize=True)
                common = [l for l in g_iids if l in set(rna[sex_map[sex]].columns)
                          and l in set(y.index)]
                if len(common) >= MIN_COVER and (best is None or len(common) > best[3]):
                    best = (int(pid), row["name"], sex, len(common), y)
        if best is None:
            print(f"{fam}: no candidate with >= {MIN_COVER} covered lines")
            continue
        pid, name, sex, n, y = best
        common, kernels, yv = multimodal.build_candidates(Kg, g_iids, rna[sex_map[sex]], y)
        r_geno = float(gblup.cross_validate(kernels["geno"], yv, n_repeats=10, seed=42).r_per_repeat.mean())
        r_exp = float(gblup.cross_validate(kernels["exp"], yv, n_repeats=10, seed=42).r_per_repeat.mean())
        best_mod = "exp" if r_exp > r_geno else "geno"
        r_best = max(r_geno, r_exp)
        perm = gblup.permutation_test(kernels[best_mod], yv, r_best ** 2 if r_best > 0 else 0,
                                      n_perm=200)
        rows.append({"family": fam, "phenotype": name, "phenotype_id": pid, "sex": sex,
                     "n": n, "r_geno": round(r_geno, 3), "r_exp": round(r_exp, 3),
                     "best_modality": best_mod, "perm_p": round(perm["p_value"], 4),
                     "predictable": bool(perm["p_value"] < 0.05)})
        print(f"{fam:18s} [{name[:24]:24s}] {sex} n={n:3d} geno={r_geno:+.3f} "
              f"exp={r_exp:+.3f} -> {best_mod} p={perm['p_value']:.3f} "
              f"{'PRED' if perm['p_value']<0.05 else 'n.s.'}")

    df = pd.DataFrame(rows)
    df.to_csv(PROC / "expanded_traits.csv", index=False)
    npred = int(df["predictable"].sum()) if len(df) else 0
    print(f"\nProbed {len(df)} new trait families; {npred} newly predictable (perm p<0.05).")
    if len(df):
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
