"""DECISIVE TEST — does RNA-seq expression break the r~0.2 genotype barrier?

Reproduces the Morgante et al. 2020 setup with the RNA-seq line-mean expression
(GSE117850) they used. Compares, leakage-free on the same lines/folds:
  GENO_lin   : additive genomic GBLUP (current baseline)
  EXP_lin    : linear expression kernel (RNA-seq)
  EXP_rbf    : Gaussian/RKHS expression kernel (captures non-additivity)
  GENO+EXP   : two-component REML (genomic + best expression kernel)
Morgante reported expression r=0.28(F)/0.38(M) for starvation vs genotype 0.07/0.15.

Output: data/processed/rnaseq_prediction.csv
"""
from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, expression, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAITS = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
N_REPEATS = 10
OUT = Path("data/processed")


def gaussian_kernel(X):
    D2 = squareform(pdist(X, "sqeuclidean"))
    return np.exp(-D2 / np.median(D2[D2 > 0]))


def _norm(K):
    return K / np.mean(np.diag(K))


def cv_r(K, y):
    return float(gblup.cross_validate(K, y, n_repeats=N_REPEATS, seed=42).r_per_repeat.mean())


def main() -> None:
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rna = {s: expression.load_rnaseq(s) for s in ("female", "male")}
    sex_map = {"F": "female", "M": "male"}
    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in ("F", "M"):
            if sex not in set(vals["sex"].unique()):
                continue
            y = dgrpool.line_means(vals, sex=sex, harmonize=True)
            e = rna[sex_map[sex]]
            common = [l for l in g_iids if l in set(e.columns) and l in set(y.index)]
            if len(common) < 110:
                continue
            yv = y.loc[common].to_numpy(float)
            gi = [g_iids.index(l) for l in common]
            G = _norm(Kg[np.ix_(gi, gi)])
            Xe = e[common].to_numpy(float).T
            keep = np.isfinite(Xe).all(0) & (Xe.std(0) > 1e-9)
            Xe = Xe[:, keep]
            Xe = (Xe - Xe.mean(0)) / Xe.std(0)
            E_lin = _norm((Xe @ Xe.T) / Xe.shape[1])
            E_rbf = _norm(gaussian_kernel(Xe))
            r_geno = cv_r(G, yv)
            r_elin = cv_r(E_lin, yv)
            r_erbf = cv_r(E_rbf, yv)
            best_E = E_rbf if r_erbf >= r_elin else E_lin
            ge = gblup.cross_validate_multi([G, best_E], yv, n_repeats=5, seed=42)["r_mean"]
            row = {"trait": trait, "sex": sex, "n": len(common), "n_genes": int(keep.sum()),
                   "GENO_lin": round(r_geno, 3), "EXP_lin": round(r_elin, 3),
                   "EXP_rbf": round(r_erbf, 3), "GENO+EXP": round(float(ge), 3)}
            rows.append(row)
            print(f"{trait:11s} {sex} n={row['n']:3d}  GENO={row['GENO_lin']:+.3f}  "
                  f"EXP_lin={row['EXP_lin']:+.3f}  EXP_rbf={row['EXP_rbf']:+.3f}  "
                  f"GENO+EXP={row['GENO+EXP']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "rnaseq_prediction.csv", index=False)
    print("\n=========== RNA-seq EXPRESSION PREDICTION ===========")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
