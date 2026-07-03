"""GENOTYPE-ONLY GBLUP test: does our validated GBLUP framework predict real
AD-transgene (Abeta42+Tau) eye-degeneration severity across DGRP genetic
backgrounds? Uses the SAME verified engine (gblup.py) as every other trait in
this project — REML, line-aware GroupKFold-style CV, permutation significance,
split-conformal calibrated intervals. No new methodology; new phenotype.

Targets: the 5 "primary" BLUP traits from script 40, pre-committed (not
post-hoc selected) per docs/EDITORIAL_RISK_ANALYSIS.md #4. BH-corrected across
the 5 within this analysis.

Output: data/processed/yang2023_gblup_results.csv
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
from drophenopredict import genotypes, gblup  # noqa: E402

PROC = Path("data/processed")
PRIMARY_TRAITS = ["nn.mean", "nn.sd", "cc.mean", "cc.sd", "mean.s.area"]


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    blup = pd.read_csv(PROC / "yang2023_blup_lines.csv", index_col=0)
    blup.index = [f"line_{n}" for n in blup.index]
    K_full, iids, meta = genotypes.load_grm("data/processed/grm.npz")
    print(f"GRM: {meta['n_snps_used']} SNPs, {len(iids)} lines")

    common = [l for l in iids if l in blup.index]
    print(f"Common lines (GRM ∩ Yang2023 phenotyped): {len(common)} / {len(blup)} phenotyped, "
          f"{len(iids)} genotyped")
    idx = [iids.index(l) for l in common]
    K = K_full[np.ix_(idx, idx)]
    K = K / np.mean(np.diag(K))

    rows = []
    for trait in PRIMARY_TRAITS:
        y = blup.loc[common, trait].to_numpy(float)
        fit = gblup.reml(K, y)
        cv = gblup.cross_validate(K, y, n_repeats=10, seed=42)
        r2 = float(cv.r2_per_repeat.mean())
        perm = gblup.permutation_test(K, y, r2, n_perm=200)
        conf = gblup.conformal_coverage(K, y, alpha=0.10, n_repeats=10)
        rows.append({
            "trait": trait, "n": len(common), "genomic_h2": round(fit.h2, 3),
            "h2_at_bound": fit.at_bound,
            "cv_R2": round(r2, 4), "cv_r": round(float(cv.r_per_repeat.mean()), 3),
            "perm_p": round(perm["p_value"], 4),
            "cover90": round(conf["empirical_coverage"], 3),
        })
        print(f"  {trait:14s} n={len(common):3d} genomic_h2={fit.h2:.3f} "
              f"CV_r={rows[-1]['cv_r']:+.3f} p={rows[-1]['perm_p']:.3f} "
              f"cover90={rows[-1]['cover90']:.2f}")

    df = pd.DataFrame(rows)
    df["fdr"] = bh(df["perm_p"].to_numpy())
    df.to_csv(PROC / "yang2023_gblup_results.csv", index=False)
    print("\n=== GENOTYPE-ONLY GBLUP: Yang 2023 AD eye-degeneration (F1 x DGRP) ===")
    print(df.to_string(index=False))
    n_sig = int((df.fdr < 0.05).sum())
    print(f"\n{n_sig}/{len(df)} primary traits predictable at FDR<0.05 "
          f"(BH-corrected across the 5 pre-committed primary traits).")


if __name__ == "__main__":
    main()
