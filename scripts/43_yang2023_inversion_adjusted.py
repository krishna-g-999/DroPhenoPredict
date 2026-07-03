"""Inversion/Wolbachia-adjusted robustness check for the Yang 2023 AD
degeneration GBLUP results (script 41). The ORIGINAL AUTHORS' OWN PCA on this
162-line subset found In(3R)Mo and In(2L)t as dominant structure (their own
finding, not ours) -- re-using our existing, tested `covariates.py`
infrastructure (same one used for our own DGRP GWAS adjustment) to check
whether the genotype-only GBLUP null (script 41: 0/5 significant) or any
future signal is robust to inversion/Wolbachia confounding, exactly as we
already did for our own traits.

Output: data/processed/yang2023_gblup_inversion_adjusted.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")
from drophenopredict import genotypes, gblup, covariates  # noqa: E402

PROC = Path("data/processed")
PRIMARY_TRAITS = ["nn.mean", "nn.sd", "cc.mean", "cc.sd", "mean.s.area"]


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    blup = pd.read_csv(PROC / "yang2023_blup_lines.csv", index_col=0)
    blup.index = [f"line_{n}" for n in blup.index]
    K_full, iids, _ = genotypes.load_grm("data/processed/grm.npz")
    cov = covariates.load_covariates()

    common = [l for l in iids if l in blup.index]
    idx = [iids.index(l) for l in common]
    K = K_full[np.ix_(idx, idx)]; K = K / np.mean(np.diag(K))
    X = covariates.design_matrix(cov, common)
    print(f"n={len(common)}, covariates: intercept + {X.shape[1]-1} cols "
          f"({list(cov.columns)})")

    rows = []
    for trait in PRIMARY_TRAITS:
        y = blup.loc[common, trait].to_numpy(float)
        cv_plain = gblup.cross_validate(K, y, n_repeats=10, seed=42)
        cv_adj = gblup.cross_validate(K, y, n_repeats=10, seed=42, X=X)
        r2_adj = float(cv_adj.r2_per_repeat.mean())
        perm_adj = gblup.permutation_test(K, y, r2_adj, n_perm=200)
        rows.append({
            "trait": trait, "n": len(common),
            "cv_r_plain": round(float(cv_plain.r_per_repeat.mean()), 3),
            "cv_r_inversion_adjusted": round(float(cv_adj.r_per_repeat.mean()), 3),
            "perm_p_adjusted": round(perm_adj["p_value"], 4),
        })
        print(f"  {trait:14s} plain_r={rows[-1]['cv_r_plain']:+.3f}  "
              f"adjusted_r={rows[-1]['cv_r_inversion_adjusted']:+.3f}  "
              f"adjusted_p={rows[-1]['perm_p_adjusted']:.3f}")

    df = pd.DataFrame(rows)
    df["fdr_adjusted"] = bh(df["perm_p_adjusted"].to_numpy())
    df.to_csv(PROC / "yang2023_gblup_inversion_adjusted.csv", index=False)
    print("\n=== Inversion/Wolbachia-ADJUSTED robustness check ===")
    print(df.to_string(index=False))
    n_sig = int((df.fdr_adjusted < 0.05).sum())
    print(f"\n{n_sig}/{len(df)} primary traits significant after inversion/Wolbachia "
          f"adjustment (FDR<0.05).")


if __name__ == "__main__":
    main()
