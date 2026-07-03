"""Strengthen the shared-biology result: genome-wide preranked GSEA per disease
model (verified engine, scripts/34) + leave-one-out robustness of the gene-level
convergence test (Findings 13).

GSEA is far better powered than the 27-gene hypergeometric test because it uses
every expressed gene's log2FC, not just genes crossing an arbitrary FDR/LFC cut.

Output: data/processed/gsea_per_model.csv, gsea_convergent_pathways.csv,
        data/processed/convergence_leave_one_out.csv
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
from drophenopredict import pathways, gsea  # noqa: E402

warnings.filterwarnings("ignore")
PROC = Path("data/processed")
MODELS = ["Abeta42", "Tau", "aSyn", "HTT128Q", "HTT200Q"]


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def de_log2fc_all(logcpm, samp, model):
    """Age-adjusted log2FC vs Elav-Gal4 control for ALL expressed genes (no cutoff)."""
    ctrl = samp[samp["model"] == "Elav-Gal4"]
    dis = samp[samp["model"] == model]
    ages = sorted(set(dis["ageDays"]) & set(ctrl["ageDays"]))
    s = pd.concat([dis[dis["ageDays"].isin(ages)], ctrl[ctrl["ageDays"].isin(ages)]])
    dvec = (s["model"] == model).astype(float).to_numpy()
    agedum = pd.get_dummies(s["ageDays"].astype(int), prefix="age", drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(s)), dvec, agedum])
    Yt = logcpm[s.index.tolist()].to_numpy(float).T
    Beta = np.linalg.inv(X.T @ X) @ X.T @ Yt
    return pd.Series(Beta[1], index=logcpm.index)


def main() -> None:
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    samp = pd.read_csv(PROC / "bcm_dmas_samples.csv", index_col=0)
    cpm = counts / counts.sum(0) * 1e6
    expressed = (cpm >= 1).sum(1) >= (0.5 * counts.shape[1])
    logcpm = np.log2(cpm[expressed] + 1)
    print(f"expressed genes: {expressed.sum()}")

    lfc = {m: de_log2fc_all(logcpm, samp, m) for m in MODELS}
    gsets = pathways.load_go_genesets(aspect="P", min_genes=10, max_genes=300)
    print(f"GO-BP gene sets tested: {len(gsets)}")

    # GSEA per model per pathway (verified engine)
    rows = []
    for m in MODELS:
        scores = lfc[m].to_dict()
        for go, (name, genes) in gsets.items():
            r = gsea.gsea_preranked(scores, set(genes), n_perm=500, seed=42)
            if r["n_hit"] >= 5:
                rows.append({"model": m, "GO_id": go, "pathway": name, **r})
    gdf = pd.DataFrame(rows)
    gdf["fdr"] = gdf.groupby("model")["p_value"].transform(lambda p: bh(p.to_numpy()))
    gdf.to_csv(PROC / "gsea_per_model.csv", index=False)

    for m in MODELS:
        sub = gdf[gdf.model == m]
        print(f"  {m:9s}: {len(sub)} pathways tested, "
              f"{int((sub.fdr < 0.1).sum())} at FDR<0.1")

    # convergence: pathways significant (FDR<0.1) in >=3 models, SAME direction
    sig = gdf[gdf.fdr < 0.1]
    piv_sig = sig.pivot_table(index="GO_id", columns="model", values="ES", aggfunc="first")
    n_sig_models = piv_sig.notna().sum(1)
    conv_ids = n_sig_models[n_sig_models >= 3].index
    conv = piv_sig.loc[conv_ids].copy()
    conv["pathway"] = [gsets[i][0] for i in conv.index]
    conv["n_models_sig"] = n_sig_models.loc[conv_ids]
    same_dir = conv[MODELS].apply(lambda row: (row.dropna() > 0).all() or (row.dropna() < 0).all(), axis=1)
    conv["concordant_direction"] = same_dir
    conv = conv.sort_values("n_models_sig", ascending=False)
    conv.to_csv(PROC / "gsea_convergent_pathways.csv")
    print(f"\npathways significant (FDR<0.1) in >=3 models: {len(conv)} "
          f"({int(same_dir.sum())} same-direction across all models where significant)")
    print(conv[["pathway", "n_models_sig", "concordant_direction"] + MODELS].head(15)
         .to_string())

    # --- leave-one-out robustness of the gene-level convergence test (Findings 13) ---
    print("\n=== leave-one-out robustness (gene-level convergence, FDR<0.05,|log2FC|>1) ===")
    def de_sig(model):
        fc = lfc[model]
        # recompute p-values with same design for the sig test (reuse earlier logic)
        ctrl = samp[samp["model"] == "Elav-Gal4"]; dis = samp[samp["model"] == model]
        ages = sorted(set(dis["ageDays"]) & set(ctrl["ageDays"]))
        s = pd.concat([dis[dis["ageDays"].isin(ages)], ctrl[ctrl["ageDays"].isin(ages)]])
        dvec = (s["model"] == model).astype(float).to_numpy()
        agedum = pd.get_dummies(s["ageDays"].astype(int), prefix="age", drop_first=True).to_numpy(float)
        X = np.column_stack([np.ones(len(s)), dvec, agedum])
        Yt = logcpm[s.index.tolist()].to_numpy(float).T
        XtX_inv = np.linalg.inv(X.T @ X); Beta = XtX_inv @ X.T @ Yt
        resid = Yt - X @ Beta; dof = X.shape[0] - X.shape[1]
        se = np.sqrt((resid ** 2).sum(0) / dof * XtX_inv[1, 1])
        t = Beta[1] / se; p = 2 * stats.t.sf(np.abs(t), dof)
        q = bh(p)
        return pd.Series((q < 0.05) & (np.abs(Beta[1]) > 1.0), index=logcpm.index)

    sig_full = pd.DataFrame({m: de_sig(m) for m in MODELS})
    rng = np.random.default_rng(42)
    for excluded in MODELS:
        subset = [m for m in MODELS if m != excluded]
        k = sig_full[subset].sum(1)
        obs = int((k >= 3).sum())
        S = sig_full[subset].to_numpy()
        null = np.array([int((np.column_stack([rng.permutation(S[:, j])
                                                for j in range(S.shape[1])]).sum(1) >= 3).sum())
                         for _ in range(500)])
        p_conv = (1 + np.sum(null >= obs)) / (1 + len(null))
        print(f"  excluding {excluded:9s}: genes DE in >=3/{len(subset)} remaining models = "
              f"{obs:3d} (null mean {null.mean():.2f}), perm p={p_conv:.4f}")


if __name__ == "__main__":
    main()
