"""Option 2 — marker-effect models vs GBLUP (does sparse QTL architecture help?).

GBLUP assumes infinitesimal (all SNPs, equal-variance). This tests the opposite:
a few large-effect SNPs. To avoid selection leakage, SNP selection happens INSIDE
each training fold (marginal GWAS on training lines only), then a regularized
model is fit on the selected SNPs and evaluated on the held-out lines.

Compares, per trait x sex: GBLUP genomic r  vs  top-k-SNP RidgeCV r.

Output: data/processed/marker_effects.csv
Run: ./.venv/Scripts/python.exe scripts/08_marker_effects.py
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
from drophenopredict import dgrpool, genotypes, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAITS = {"climbing": 1543, "startle": 2799, "starvation": 2798,
          "sleep": 1483, "lifespan": 1315}
TOP_K = 500
N_REPEATS = 10
N_SPLITS = 5
OUT = Path("data/processed")
SNP_CACHE = OUT / "snp_matrix.npz"
BED = "data/raw/dgrp2/dgrp2.bed"


def get_snp_matrix():
    if SNP_CACHE.exists():
        d = np.load(SNP_CACHE, allow_pickle=True)
        return d["X"], [str(x) for x in d["iids"]]
    X, iids, idx = genotypes.build_snp_matrix(BED, target_snps=100_000)
    np.savez_compressed(SNP_CACHE, X=X, iids=np.array(iids), snp_index=idx)
    return X, iids


def marker_cv(X, y, *, seed=42):
    """Within-fold GWAS-select top-k SNPs -> RidgeCV. Leakage-free."""
    n = y.size
    r2s, rs = [], []
    for rep in range(N_REPEATS):
        rng = np.random.default_rng(seed + rep)
        idx = rng.permutation(n)
        folds = np.array_split(idx, N_SPLITS)
        yt_all, yp_all = [], []
        for k in range(N_SPLITS):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(N_SPLITS) if j != k])
            Xtr, Xte = X[train], X[test]
            ytr = y[train]
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr = (Xtr - mu) / sd
            yc = ytr - ytr.mean()
            corr = (Ztr.T @ yc) / (len(train) * (yc.std() + 1e-9))   # marginal r per SNP
            sel = np.argsort(-np.abs(corr))[:TOP_K]
            Zte = (Xte - mu) / sd
            model = RidgeCV(alphas=np.logspace(-1, 4, 12))
            model.fit(Ztr[:, sel], ytr)
            yp_all.append(model.predict(Zte[:, sel]))
            yt_all.append(y[test])
        yt = np.concatenate(yt_all); yp = np.concatenate(yp_all)
        r2s.append(1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2))
        rs.append(np.corrcoef(yt, yp)[0, 1])
    return float(np.mean(rs)), float(np.mean(r2s))


def main() -> None:
    X, iids = get_snp_matrix()
    print(f"SNP matrix: {X.shape[0]} lines x {X.shape[1]} SNPs")
    Kg, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")
    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex in ("F", "M"):
            if sex not in set(vals["sex"].unique()):
                continue
            y = dgrpool.line_means(vals, sex=sex, harmonize=True)
            common = [l for l in iids if l in set(y.index)]
            if len(common) < 110:
                continue
            xi = [iids.index(l) for l in common]
            yv = y.loc[common].to_numpy(float)
            mk_r, mk_r2 = marker_cv(X[xi], yv)
            # GBLUP on same common lines
            gi = [g_iids.index(l) for l in common]
            Gs = Kg[np.ix_(gi, gi)]
            gb = gblup.cross_validate(Gs, yv, n_repeats=N_REPEATS, seed=42)
            row = {"trait": trait, "sex": sex, "n_lines": len(common),
                   "gblup_r": round(float(gb.r_per_repeat.mean()), 3),
                   "marker_r": round(mk_r, 3), "delta_r": round(mk_r - float(gb.r_per_repeat.mean()), 3)}
            rows.append(row)
            print(f"{trait:11s} {sex} n={len(common):3d}  GBLUP_r={row['gblup_r']:+.3f}  "
                  f"marker(top{TOP_K})_r={row['marker_r']:+.3f}  Δr={row['delta_r']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "marker_effects.csv", index=False)
    print("\n============== MARKER-EFFECT vs GBLUP ==============")
    print(df.to_string(index=False))
    print(f"\nmean Δr (marker − GBLUP) = {df['delta_r'].mean():+.3f}")
    print(f"Saved -> {OUT/'marker_effects.csv'}")


if __name__ == "__main__":
    main()
