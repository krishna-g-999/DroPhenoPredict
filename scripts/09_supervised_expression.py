"""Option 1 (final check) — SUPERVISED expression model vs GBLUP.

The TRM (unsupervised) and G+T REML tests found no transcriptomic gain. This
closes the loop with the most favourable expression method: within-fold select
the top-k genes most correlated with the trait (training lines only -> no
leakage), fit RidgeCV, evaluate on held-out lines. If even supervised gene
selection fails to beat GBLUP, the no-gain conclusion is robust.

Output: data/processed/supervised_expression.csv
"""
from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, expression, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAITS = {"climbing": 1543, "startle": 2799, "starvation": 2798,
          "sleep": 1483, "lifespan": 1315}
TOP_K = 200
N_REPEATS, N_SPLITS = 10, 5
OUT = Path("data/processed")


def supervised_cv(X, y, seed=42):
    n = y.size
    rs = []
    for rep in range(N_REPEATS):
        rng = np.random.default_rng(seed + rep)
        folds = np.array_split(rng.permutation(n), N_SPLITS)
        yt_all, yp_all = [], []
        for k in range(N_SPLITS):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(N_SPLITS) if j != k])
            Xtr, ytr = X[train], y[train]
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr = (Xtr - mu) / sd
            yc = ytr - ytr.mean()
            corr = (Ztr.T @ yc) / (len(train) * (yc.std() + 1e-9))
            sel = np.argsort(-np.abs(corr))[:TOP_K]
            m = RidgeCV(alphas=np.logspace(-1, 4, 12)).fit(Ztr[:, sel], ytr)
            yp_all.append(m.predict(((X[test] - mu) / sd)[:, sel]))
            yt_all.append(y[test])
        yt = np.concatenate(yt_all); yp = np.concatenate(yp_all)
        rs.append(np.corrcoef(yt, yp)[0, 1])
    return float(np.mean(rs))


def main() -> None:
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
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
            common = [l for l in g_iids if l in set(e.columns) and l in set(y.index)]
            if len(common) < 110:
                continue
            yv = y.loc[common].to_numpy(float)
            X = e[common].to_numpy(float).T          # lines x genes
            sup_r = supervised_cv(X, yv)
            gi = [g_iids.index(l) for l in common]
            gb = gblup.cross_validate(Kg[np.ix_(gi, gi)], yv, n_repeats=N_REPEATS, seed=42)
            row = {"trait": trait, "sex": sex, "n_lines": len(common),
                   "gblup_r": round(float(gb.r_per_repeat.mean()), 3),
                   "sup_expr_r": round(sup_r, 3),
                   "delta_r": round(sup_r - float(gb.r_per_repeat.mean()), 3)}
            rows.append(row)
            print(f"{trait:11s} {sex} n={len(common):3d}  GBLUP_r={row['gblup_r']:+.3f}  "
                  f"sup_expr(top{TOP_K})_r={row['sup_expr_r']:+.3f}  Δr={row['delta_r']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "supervised_expression.csv", index=False)
    print("\n========== SUPERVISED EXPRESSION vs GBLUP ==========")
    print(df.to_string(index=False))
    print(f"\nmean Δr (sup-expr − GBLUP) = {df['delta_r'].mean():+.3f}")


if __name__ == "__main__":
    main()
