"""Can non-linear (RKHS Gaussian) kernels break the r~0.2 barrier?

Drosophila quantitative traits are epistasis-dominated (Huang 2012 PNAS), and
Morgante 2020 reported expression predicts starvation far better than genotype
using a Gaussian-kernel (RKHS) model. My Findings 05 used a LINEAR expression
kernel on MICROARRAY data and found ~0. This script disentangles METHOD
(linear vs Gaussian/RKHS) from everything else, on the same lines/folds:

  GRM_linear  (additive genomic, GBLUP)         <- current baseline
  GRM_rbf     (genomic, Gaussian kernel = epistasis-capturing)
  TRM_linear  (additive expression)
  TRM_rbf     (expression, Gaussian kernel)

All via the SAME validated gblup engine (any PSD kernel = RKHS BLUP).
Output: data/processed/nonlinear_kernels.csv
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


def gaussian_kernel(X: np.ndarray) -> np.ndarray:
    """RKHS Gaussian kernel with the median-distance heuristic bandwidth."""
    D2 = squareform(pdist(X, "sqeuclidean"))
    med = np.median(D2[D2 > 0])
    return np.exp(-D2 / med)          # gamma = 1/median(sq dist)


def _norm(K):
    return K / np.mean(np.diag(K))


def cv_r(K, y):
    return float(gblup.cross_validate(K, y, n_repeats=N_REPEATS, seed=42).r_per_repeat.mean())


def main() -> None:
    Kg_lin, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    snp = np.load("data/processed/snp_matrix.npz", allow_pickle=True)
    X_snp, snp_iids = snp["X"], [str(x) for x in snp["iids"]]
    expr = {s: expression.load_expression(s) for s in ("female", "male")}
    sex_map = {"F": "female", "M": "male"}

    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in ("F", "M"):
            if sex not in set(vals["sex"].unique()):
                continue
            y = dgrpool.line_means(vals, sex=sex, harmonize=True)
            e = expr[sex_map[sex]]
            common = [l for l in g_iids if l in set(snp_iids) and l in set(e.columns)
                      and l in set(y.index)]
            if len(common) < 110:
                continue
            yv = y.loc[common].to_numpy(float)
            gi = [g_iids.index(l) for l in common]
            si = [snp_iids.index(l) for l in common]
            # genomic kernels
            GRM_lin = _norm(Kg_lin[np.ix_(gi, gi)])
            Xs = (X_snp[si] - X_snp[si].mean(0)) / (X_snp[si].std(0) + 1e-9)
            GRM_rbf = _norm(gaussian_kernel(Xs))
            # expression kernels
            Xe = e[common].to_numpy(float).T
            Xe = (Xe - Xe.mean(0)) / (Xe.std(0) + 1e-9)
            TRM_lin = _norm((Xe @ Xe.T) / Xe.shape[1])
            TRM_rbf = _norm(gaussian_kernel(Xe))

            row = {"trait": trait, "sex": sex, "n": len(common),
                   "GRM_linear": round(cv_r(GRM_lin, yv), 3),
                   "GRM_rbf": round(cv_r(GRM_rbf, yv), 3),
                   "TRM_linear": round(cv_r(TRM_lin, yv), 3),
                   "TRM_rbf": round(cv_r(TRM_rbf, yv), 3)}
            rows.append(row)
            print(f"{trait:11s} {sex} n={row['n']:3d}  GRM_lin={row['GRM_linear']:+.3f} "
                  f"GRM_rbf={row['GRM_rbf']:+.3f}  TRM_lin={row['TRM_linear']:+.3f} "
                  f"TRM_rbf={row['TRM_rbf']:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "nonlinear_kernels.csv", index=False)
    print("\n============ LINEAR vs GAUSSIAN(RKHS) KERNELS ============")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
