"""Cross-disease molecular convergence in the BCM-DMAS 5-model fly panel.

Tests the v2 premise (Mallik & Mukhopadhyay 2021: shared neurodegeneration
biology) directly, in flies, on MATCHED data. No DGRP bridge, no behaviour
labels needed.

Steps (each verified):
  1. Per disease model, differential expression vs the driver-matched control
     (Elav-Gal4/+), age-adjusted, via a linear model on log2-CPM.
     -> VERIFIED against statsmodels OLS on a random gene before use.
  2. Convergence: how many genes are DE in >=k models, is that MORE than chance
     (column-permutation null), and are shared genes DIRECTIONALLY concordant.
  3. Shared pathways: GO-BP enrichment (hypergeometric) of the convergent genes.

Control = Elav-Gal4/+ (group A) isolates the disease transgene from the Gal4 driver.
Outputs: data/processed/shared_biology_*.csv
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
from drophenopredict import pathways  # noqa: E402

warnings.filterwarnings("ignore")
PROC = Path("data/processed")
MODELS = ["Abeta42", "Tau", "aSyn", "HTT128Q", "HTT200Q"]
FDR, LFC = 0.05, 1.0


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def de_one(logcpm, samp, model):
    """Age-adjusted DE (disease vs Elav-Gal4 control). Returns df: log2FC,p,fdr + design for verify."""
    ctrl = samp[samp["model"] == "Elav-Gal4"]
    dis = samp[samp["model"] == model]
    ages = sorted(set(dis["ageDays"]) & set(ctrl["ageDays"]))
    s = pd.concat([dis[dis["ageDays"].isin(ages)], ctrl[ctrl["ageDays"].isin(ages)]])
    ids = s.index.tolist()
    dvec = (s["model"] == model).astype(float).to_numpy()
    agedum = pd.get_dummies(s["ageDays"].astype(int), prefix="age", drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(s)), dvec, agedum])     # col 1 = disease
    Yt = logcpm[ids].to_numpy(float).T                        # samples x genes
    XtX_inv = np.linalg.inv(X.T @ X)
    Beta = XtX_inv @ X.T @ Yt                                 # p x genes
    resid = Yt - X @ Beta
    dof = X.shape[0] - X.shape[1]
    sigma2 = (resid ** 2).sum(0) / dof
    se = np.sqrt(sigma2 * XtX_inv[1, 1])
    t = Beta[1] / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    df = pd.DataFrame({"log2FC": Beta[1], "t": t, "p": p}, index=logcpm.index)
    df["fdr"] = bh(df["p"].to_numpy())
    return df, (X, Yt, ids)


def verify_de(logcpm, samp):
    """Confirm vectorized OLS == statsmodels for one gene/model."""
    import statsmodels.api as sm
    df, (X, Yt, ids) = de_one(logcpm, samp, "aSyn")
    g = 100                                                   # arbitrary gene index
    m = sm.OLS(Yt[:, g], X).fit()
    mine_t = df["t"].iloc[g]; mine_fc = df["log2FC"].iloc[g]
    print(f"[verify] gene {logcpm.index[g]}: mine log2FC={mine_fc:+.4f} t={mine_t:+.4f} | "
          f"statsmodels coef={m.params[1]:+.4f} t={m.tvalues[1]:+.4f}")
    assert np.allclose(mine_fc, m.params[1], atol=1e-6) and np.allclose(mine_t, m.tvalues[1], atol=1e-4)
    print("[verify] vectorized DE matches statsmodels -> OK")


def main() -> None:
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    samp = pd.read_csv(PROC / "bcm_dmas_samples.csv", index_col=0)
    cpm = counts / counts.sum(0) * 1e6
    expressed = (cpm >= 1).sum(1) >= (0.5 * counts.shape[1])
    logcpm = np.log2(cpm[expressed] + 1)
    print(f"genes: {counts.shape[0]} total -> {expressed.sum()} expressed (CPM>=1 in >=50% samples)")

    verify_de(logcpm, samp)

    # per-model DE
    de = {m: de_one(logcpm, samp, m)[0] for m in MODELS}
    lfc = pd.DataFrame({m: de[m]["log2FC"] for m in MODELS})
    sig = pd.DataFrame({m: (de[m]["fdr"] < FDR) & (de[m]["log2FC"].abs() > LFC) for m in MODELS})
    for m in MODELS:
        print(f"  {m:9s}: {int(sig[m].sum()):4d} DE genes (FDR<{FDR}, |log2FC|>{LFC})")

    # convergence: genes DE in >=k models
    k = sig.sum(1)
    print("\ngenes DE in exactly k models:", {int(i): int((k == i).sum()) for i in range(1, 6)})
    obs = int((k >= 3).sum())

    # permutation null (independently shuffle each model's DE labels)
    rng = np.random.default_rng(42)
    S = sig.to_numpy()
    null = []
    for _ in range(2000):
        Sp = np.column_stack([rng.permutation(S[:, j]) for j in range(S.shape[1])])
        null.append(int((Sp.sum(1) >= 3).sum()))
    null = np.array(null)
    p_conv = (1 + np.sum(null >= obs)) / (1 + len(null))
    print(f"genes DE in >=3 models: observed={obs}, null mean={null.mean():.1f}, "
          f"perm p={p_conv:.4f}")

    # directional concordance across models (genome-wide log2FC correlation)
    print("\ncross-disease log2FC Spearman correlation (genome-wide):")
    corr = lfc.corr(method="spearman")
    print(corr.round(2).to_string())

    # shared convergent genes (DE in >=3 models) + their direction
    shared = k[k >= 3].index
    shared_df = lfc.loc[shared].copy()
    shared_df["n_models_DE"] = k.loc[shared]
    shared_df["mean_log2FC"] = shared_df[MODELS].mean(1)
    shared_df = shared_df.sort_values("n_models_DE", ascending=False)
    shared_df.to_csv(PROC / "shared_biology_convergent_genes.csv")

    # GO-BP enrichment of convergent genes (hypergeometric vs expressed background)
    gsets = pathways.load_go_genesets(aspect="P", min_genes=10, max_genes=500)
    background = set(logcpm.index)
    hits = set(shared) & background
    rows = []
    for go, (name, genes) in gsets.items():
        gb = genes & background
        if len(gb) < 10:
            continue
        overlap = len(gb & hits)
        if overlap < 3:
            continue
        # hypergeometric: overlap of hits with pathway
        pval = stats.hypergeom.sf(overlap - 1, len(background), len(gb), len(hits))
        rows.append({"GO_id": go, "pathway": name, "n_path": len(gb),
                     "overlap": overlap, "p": pval})
    enr = pd.DataFrame(rows).sort_values("p")
    enr["fdr"] = bh(enr["p"].to_numpy())
    enr.head(40).to_csv(PROC / "shared_biology_pathways.csv", index=False)

    print(f"\nconvergent genes (DE in >=3 models): {len(shared)}")
    print("top shared pathways (GO-BP, hypergeometric):")
    print(enr.head(12)[["pathway", "n_path", "overlap", "p", "fdr"]].to_string(index=False))
    print(f"\nSaved -> shared_biology_convergent_genes.csv + shared_biology_pathways.csv")


if __name__ == "__main__":
    main()
