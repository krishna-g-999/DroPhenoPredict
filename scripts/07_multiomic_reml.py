"""Option 1 (fair test) — does transcriptome ADD predictive value over genotype?

Two-variance-component REML model (genomic G + transcriptomic T), variance
components estimated from data and refit per CV fold (Morgante et al. 2020, G3).
For each trait x sex: genomic-only CV r vs G+T CV r, plus the transcriptomic
variance fraction p_T estimated on full data.

Output: data/processed/multiomic_reml.csv
Run: ./.venv/Scripts/python.exe scripts/07_multiomic_reml.py
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

TRAITS = {"climbing": 1543, "startle": 2799, "starvation": 2798,
          "sleep": 1483, "lifespan": 1315}
MIN_LINES = 110
N_REPEATS = 5
OUT = Path("data/processed")


def _norm(K):
    return K / np.mean(np.diag(K))


def _subset(K, iids, lines):
    idx = [iids.index(l) for l in lines]
    return K[np.ix_(idx, idx)]


def main() -> None:
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    trms = {s: expression.load_trm(OUT / f"trm_{s}.npz") for s in ("female", "male")}
    sex_map = {"F": "female", "M": "male"}
    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in ("F", "M"):
            if sex not in set(vals["sex"].unique()):
                continue
            y = dgrpool.line_means(vals, sex=sex, harmonize=True)
            Kt, t_iids, _ = trms[sex_map[sex]]
            common = [l for l in g_iids if l in set(t_iids) and l in set(y.index)]
            if len(common) < MIN_LINES:
                continue
            yv = y.loc[common].to_numpy(float)
            G = _norm(_subset(Kg, g_iids, common))
            T = _norm(_subset(Kt, t_iids, common))

            # variance partition on full data
            vc = gblup.reml_multi([G, T], yv)
            tot = vc.sum()
            p_g, p_t = vc[0] / tot, vc[1] / tot

            geno = gblup.cross_validate(G, yv, n_repeats=N_REPEATS, seed=42)
            gt = gblup.cross_validate_multi([G, T], yv, n_repeats=N_REPEATS, seed=42)
            row = {
                "trait": trait, "sex": sex, "n_lines": len(common),
                "var_frac_geno": round(float(p_g), 3),
                "var_frac_trans": round(float(p_t), 3),
                "geno_r": round(float(geno.r_per_repeat.mean()), 3),
                "GplusT_r": round(float(gt["r_mean"]), 3),
                "delta_r": round(float(gt["r_mean"] - geno.r_per_repeat.mean()), 3),
            }
            rows.append(row)
            print(f"{trait:11s} {sex}  n={row['n_lines']:3d}  varT={row['var_frac_trans']:.2f}  "
                  f"geno_r={row['geno_r']:+.3f}  G+T_r={row['GplusT_r']:+.3f}  "
                  f"Δr={row['delta_r']:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "multiomic_reml.csv", index=False)
    print("\n============== MULTI-OMIC (G vs G+T, REML) ==============")
    print(df.to_string(index=False))
    print(f"\nmean Δr (G+T − G) = {df['delta_r'].mean():+.3f}")
    print(f"Saved -> {OUT/'multiomic_reml.csv'}")


if __name__ == "__main__":
    main()
