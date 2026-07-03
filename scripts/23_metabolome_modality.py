"""Tier 1 — test the metabolome as a 3rd modality + VALIDATE vs Morgante/Heredity 2021.

Builds a metabolomic relationship matrix (MRM) from the NMR bins (script 22) and
compares genotype | expression | metabolome | combined prediction for male traits.
CRITICAL VALIDATION: the metabolome should predict LOCOMOTOR ACTIVITY at r≈0.4
(Heredity 2021) vs genotype <0.1 — reproducing that confirms the NMR processing.
NMR is on MALE pools, so metabolome predictions are valid for male traits only.

Output: data/processed/metabolome_modality.csv
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
from drophenopredict import dgrpool, genotypes, expression, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

PROC = Path("data/processed")
# male behavioural traits we already model + activity (paper's headline metabolome trait)
CORE = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}


def _norm(K):
    return K / np.mean(np.diag(K))


def kernel_from_features(X):
    keep = np.isfinite(X).all(0) & (X.std(0) > 1e-12)
    Z = X[:, keep]; Z = (Z - Z.mean(0)) / Z.std(0)
    return _norm((Z @ Z.T) / Z.shape[1])


def cvr(K, y):
    return float(gblup.cross_validate(K, y, n_repeats=10, seed=42).r_per_repeat.mean())


def find_activity_phenotypes():
    cat = dgrpool.numeric_continuous_phenotypes(dgrpool.fetch_phenotype_catalog())
    hits = cat[cat["name"].str.contains("activ|locomot", case=False, na=False)]
    return [(int(i), r["name"]) for i, r in hits.iterrows()]


def main() -> None:
    met = pd.read_parquet(PROC / "metabolome_matrix.parquet")     # lines x bins
    print(f"Metabolome: {met.shape[0]} lines x {met.shape[1]} bins")
    Kg_full, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna_m = expression.load_rnaseq("male")

    targets = [(f"{t}_M", pid, "M") for t, pid in CORE.items()]
    for pid, name in find_activity_phenotypes():
        targets.append((f"activity[{name[:18]}]", pid, "M"))

    rows = []
    for label, pid, sex in targets:
        try:
            y = dgrpool.line_means(dgrpool.fetch_phenotype_values(pid), sex=sex, harmonize=True)
        except Exception:
            continue
        common = [l for l in met.index if l in g_iids and l in rna_m.columns and l in set(y.index)]
        if len(common) < 90:
            continue
        yv = y.loc[common].to_numpy(float)
        gi = [g_iids.index(l) for l in common]
        Kg = _norm(Kg_full[np.ix_(gi, gi)])
        Km = kernel_from_features(met.loc[common].to_numpy(float))
        Ke = kernel_from_features(rna_m[common].to_numpy(float).T)
        r_g, r_m, r_e = cvr(Kg, yv), cvr(Km, yv), cvr(Ke, yv)
        comb = gblup.cross_validate_multi([Kg, Km], yv, n_repeats=5, seed=42)["r_mean"]
        best = max([("geno", r_g), ("exp", r_e), ("metab", r_m)], key=lambda z: z[1])[0]
        rows.append({"trait": label, "n": len(common), "geno": round(r_g, 3),
                     "exp": round(r_e, 3), "metab": round(r_m, 3),
                     "geno+metab": round(float(comb), 3), "best": best})
        print(f"{label:24s} n={len(common):3d} geno={r_g:+.3f} exp={r_e:+.3f} "
              f"metab={r_m:+.3f} g+m={comb:+.3f} -> {best}")

    df = pd.DataFrame(rows)
    df.to_csv(PROC / "metabolome_modality.csv", index=False)
    print("\n=== METABOLOME AS 3rd MODALITY (male) ===")
    print(df.to_string(index=False))
    act = df[df["trait"].str.startswith("activity")]
    if len(act):
        print(f"\nVALIDATION (locomotor activity): metab r up to {act['metab'].max():+.3f} "
              f"vs geno {act['geno'].max():+.3f}  (Heredity 2021: metab>0.4, geno<0.1)")


if __name__ == "__main__":
    main()
