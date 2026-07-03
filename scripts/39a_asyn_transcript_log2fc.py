"""Compute aSyn-vs-control transcript log2FC at TWO resolutions, both saved as
permanent artifacts, for later comparison against the TMT proteomics (script 38
must be run first to get the protein data; this script needs only what's
already in hand).

Two estimates, each with a real trade-off, both reported for honesty:
  (a) age-adjusted, all 6 timepoints (days 2-21), n=18 per group — more
      statistical power, but NOT timepoint-matched to the proteomics (day 10 only).
  (b) day-10-only, n=3 per group — exactly timepoint-matched to the proteomics,
      but much noisier (a simple two-sample comparison, no age adjustment needed
      since there is only one timepoint).

Output: data/processed/asyn_log2fc_ageadjusted.csv
        data/processed/asyn_log2fc_day10only.csv
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
warnings.filterwarnings("ignore")
PROC = Path("data/processed")


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    samp = pd.read_csv(PROC / "bcm_dmas_samples.csv", index_col=0)
    cpm = counts / counts.sum(0) * 1e6
    expressed = (cpm >= 1).sum(1) >= (0.5 * counts.shape[1])
    logcpm = np.log2(cpm[expressed] + 1)
    print(f"expressed genes: {expressed.sum()}")

    ctrl = samp[samp["model"] == "Elav-Gal4"]
    dis = samp[samp["model"] == "aSyn"]

    # ---- (a) age-adjusted, all shared timepoints ----
    ages = sorted(set(dis["ageDays"]) & set(ctrl["ageDays"]))
    s = pd.concat([dis[dis["ageDays"].isin(ages)], ctrl[ctrl["ageDays"].isin(ages)]])
    dvec = (s["model"] == "aSyn").astype(float).to_numpy()
    agedum = pd.get_dummies(s["ageDays"].astype(int), prefix="age", drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(s)), dvec, agedum])
    Yt = logcpm[s.index.tolist()].to_numpy(float).T
    XtX_inv = np.linalg.inv(X.T @ X); Beta = XtX_inv @ X.T @ Yt
    resid = Yt - X @ Beta; dof = X.shape[0] - X.shape[1]
    se = np.sqrt((resid ** 2).sum(0) / dof * XtX_inv[1, 1])
    t = Beta[1] / se; p = 2 * stats.t.sf(np.abs(t), dof)
    df_a = pd.DataFrame({"log2FC": Beta[1], "p": p, "fdr": bh(p)}, index=logcpm.index)
    df_a.to_csv(PROC / "asyn_log2fc_ageadjusted.csv")
    print(f"(a) age-adjusted (n={len(s)}, days {ages}): "
          f"{int((df_a.fdr < 0.05).sum())} FDR<0.05 genes")

    # ---- (b) day-10-only, simple two-sample t-test ----
    ctrl10 = ctrl[ctrl["ageDays"] == 10]
    dis10 = dis[dis["ageDays"] == 10]
    print(f"day-10 samples: control n={len(ctrl10)}, aSyn n={len(dis10)}")
    Yc = logcpm[ctrl10.index.tolist()].to_numpy(float)
    Yd = logcpm[dis10.index.tolist()].to_numpy(float)
    log2fc = Yd.mean(1) - Yc.mean(1)
    tstat, pval = stats.ttest_ind(Yd, Yc, axis=1, equal_var=False)
    df_b = pd.DataFrame({"log2FC": log2fc, "p": pval, "fdr": bh(pval)}, index=logcpm.index)
    df_b.to_csv(PROC / "asyn_log2fc_day10only.csv")
    print(f"(b) day-10-only (n=3 vs 3): {int((df_b.fdr < 0.05).sum())} FDR<0.05 genes")

    # sanity: the two estimates should correlate reasonably (same underlying biology,
    # different statistical power / timepoint scope)
    common = df_a.index.intersection(df_b.index)
    r, pr = stats.pearsonr(df_a.loc[common, "log2FC"], df_b.loc[common, "log2FC"])
    print(f"\nConcordance between the two transcript estimates: "
          f"Pearson r={r:.3f} (p={pr:.2e}), n={len(common)} genes")
    print("(expected: positive and reasonably strong, since both measure the same underlying "
          "aSyn-vs-control effect; day-10-only is a subset/noisier version of the full model)")


if __name__ == "__main__":
    main()
