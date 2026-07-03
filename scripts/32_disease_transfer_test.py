"""HONEST transfer test: does the DGRP natural-variation gene->climbing map
predict the DIRECTION of disease-model expression change?

Logic (no fabrication, no locked behaviour labels needed):
  - DGRP: for each gene, Pearson r between its female RNA-seq expression and
    female climbing across ~180 DGRP lines  =  "climbing weight" (natural variation).
  - BCM-DMAS: for each disease model, log2 fold-change of head expression vs the
    Elav-Gal4 control (age-matched)  =  "disease effect".
  - All 5 disease models LOSE climbing with age. If the DGRP gene->climbing map
    transfers, genes positively associated with climbing should be systematically
    DOWN-regulated in disease -> NEGATIVE Spearman(climbing_weight, disease_log2FC).
  Significance by gene-label permutation.

Caveat (stated honestly): DGRP expression is whole-body; BCM-DMAS is head; both
female. Tissue mismatch weakens but does not invalidate a directional test.

Output: data/processed/disease_transfer_test.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, expression  # noqa: E402

warnings.filterwarnings("ignore")
PROC = Path("data/processed")


def dgrp_climbing_weights() -> pd.Series:
    """Per-gene Pearson r(expression, climbing) across DGRP females (FBgn -> r)."""
    y = dgrpool.line_means(dgrpool.fetch_phenotype_values(1543), sex="F", harmonize=True)
    expr = expression.load_rnaseq("female")                 # genes x lines
    common = [l for l in expr.columns if l in set(y.index)]
    X = expr[common].to_numpy(float)
    keep = np.isfinite(X).all(1) & (X.std(1) > 1e-9)
    X = X[keep]; genes = expr.index[keep]
    yv = y.loc[common].to_numpy(float); yc = yv - yv.mean()
    Xc = X - X.mean(1, keepdims=True)
    r = (Xc @ yc) / (np.sqrt((Xc ** 2).sum(1)) * np.sqrt((yc ** 2).sum()) + 1e-12)
    return pd.Series(r, index=genes)


def disease_log2fc(counts: pd.DataFrame, samp: pd.DataFrame) -> dict:
    """log2FC of each disease model vs Elav-Gal4 control, age-matched, per gene."""
    cpm = counts / counts.sum(0) * 1e6
    logcpm = np.log2(cpm + 1)
    ctrl = samp[samp["model"] == "Elav-Gal4"]
    out = {}
    for model in ["Abeta42", "Tau", "aSyn", "HTT128Q", "HTT200Q"]:
        dsamp = samp[samp["model"] == model]
        ages = sorted(set(dsamp["ageDays"]) & set(ctrl["ageDays"]))
        fc = []
        for a in ages:
            d = logcpm[dsamp[dsamp["ageDays"] == a].index].mean(1)
            c = logcpm[ctrl[ctrl["ageDays"] == a].index].mean(1)
            fc.append(d - c)
        out[model] = pd.concat(fc, axis=1).mean(1)          # mean log2FC across ages
    return out


def main() -> None:
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    samp = pd.read_csv(PROC / "bcm_dmas_samples.csv", index_col=0)
    weights = dgrp_climbing_weights()
    print(f"DGRP climbing weights: {len(weights)} genes")
    fcs = disease_log2fc(counts, samp)

    rng = np.random.default_rng(42)
    rows = []
    for model, fc in fcs.items():
        genes = weights.index.intersection(fc.index)
        w = weights.loc[genes].to_numpy()
        f = fc.loc[genes].to_numpy()
        ok = np.isfinite(w) & np.isfinite(f)
        w, f = w[ok], f[ok]
        rho, _ = stats.spearmanr(w, f)
        # permutation null: shuffle gene labels of weights
        null = np.array([stats.spearmanr(rng.permutation(w), f)[0] for _ in range(500)])
        p = (1 + np.sum(np.abs(null) >= abs(rho))) / (1 + len(null))
        disease = samp[samp["model"] == model]["disease"].iloc[0]
        rows.append({"disease": disease, "model": model, "n_genes": int(ok.sum()),
                     "spearman_rho": round(float(rho), 4), "perm_p": round(float(p), 4),
                     "direction": "as predicted (neg)" if rho < 0 else "opposite (pos)"})
        print(f"{disease:8s} {model:9s} n={int(ok.sum()):5d}  rho={rho:+.3f}  p={p:.3f}  "
              f"{'NEG (good-climbing genes down in disease)' if rho<0 else 'POS'}")

    df = pd.DataFrame(rows)
    df.to_csv(PROC / "disease_transfer_test.csv", index=False)
    print("\n=== DGRP climbing-map -> disease-signature transfer ===")
    print(df.to_string(index=False))
    neg_sig = df[(df.spearman_rho < 0) & (df.perm_p < 0.05)]
    print(f"\n{len(neg_sig)}/{len(df)} disease models show the predicted negative transfer "
          f"(good-climbing genes down-regulated in disease, p<0.05).")


if __name__ == "__main__":
    main()
