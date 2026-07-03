"""Complete the GSEA convergence analysis from the already-computed
data/processed/gsea_per_model.csv (script 35's expensive part succeeded; a
downstream pivot bug crashed on Tau, which had ZERO significant pathways at
FDR<0.1 and so was silently absent as a pivot column — reindex fixes it).

Also runs the leave-one-out robustness check (gene-level, Findings 13) that
never got to execute due to the same crash.
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


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    gdf = pd.read_csv(PROC / "gsea_per_model.csv")
    print("per-model pathways tested / significant (FDR<0.1):")
    for m in MODELS:
        sub = gdf[gdf.model == m]
        print(f"  {m:9s}: {len(sub)} tested, {int((sub.fdr < 0.1).sum())} significant")

    sig = gdf[gdf.fdr < 0.1]
    piv_sig = sig.pivot_table(index="GO_id", columns="model", values="ES", aggfunc="first")
    piv_sig = piv_sig.reindex(columns=MODELS)          # <-- the fix: keep all 5 columns even if empty
    n_sig_models = piv_sig.notna().sum(1)
    conv_ids = n_sig_models[n_sig_models >= 3].index
    conv = piv_sig.loc[conv_ids].copy()

    name_map = {}
    gsets = None
    if len(conv_ids):
        gsets = pathways.load_go_genesets(aspect="P", min_genes=10, max_genes=300)
        name_map = {go: gsets[go][0] for go in conv_ids if go in gsets}
    conv["pathway"] = [name_map.get(i, i) for i in conv.index]
    conv["n_models_sig"] = n_sig_models.loc[conv_ids]
    same_dir = conv[MODELS].apply(
        lambda row: (row.dropna() > 0).all() or (row.dropna() < 0).all(), axis=1)
    conv["concordant_direction"] = same_dir
    conv = conv.sort_values("n_models_sig", ascending=False)
    conv.to_csv(PROC / "gsea_convergent_pathways.csv")
    print(f"\npathways significant (FDR<0.1) in >=3 models: {len(conv)} "
          f"({int(same_dir.sum())} same-direction where significant)")
    if len(conv):
        print(conv[["pathway", "n_models_sig", "concordant_direction"] + MODELS]
             .head(20).to_string())
    else:
        print("(none reached the >=3-model threshold at FDR<0.1 genome-wide)")

    # ------------------------------------------------- leave-one-out (gene-level) --
    print("\n=== leave-one-out robustness (gene-level convergence, Findings 13) ===")
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    samp = pd.read_csv(PROC / "bcm_dmas_samples.csv", index_col=0)
    cpm = counts / counts.sum(0) * 1e6
    expressed = (cpm >= 1).sum(1) >= (0.5 * counts.shape[1])
    logcpm = np.log2(cpm[expressed] + 1)

    def de_sig(model):
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
    loo_rows = []
    for excluded in MODELS:
        subset = [m for m in MODELS if m != excluded]
        k = sig_full[subset].sum(1)
        obs = int((k >= 3).sum())
        S = sig_full[subset].to_numpy()
        null = np.array([int((np.column_stack([rng.permutation(S[:, j])
                                                for j in range(S.shape[1])]).sum(1) >= 3).sum())
                         for _ in range(500)])
        p_conv = (1 + np.sum(null >= obs)) / (1 + len(null))
        loo_rows.append({"excluded": excluded, "n_remaining_models": len(subset),
                         "genes_ge3": obs, "null_mean": round(float(null.mean()), 2),
                         "perm_p": round(float(p_conv), 4)})
        print(f"  excluding {excluded:9s}: genes DE in >=3/{len(subset)} remaining models = "
              f"{obs:3d} (null mean {null.mean():.2f}), perm p={p_conv:.4f}")
    pd.DataFrame(loo_rows).to_csv(PROC / "convergence_leave_one_out.csv", index=False)
    print(f"\nSaved -> gsea_convergent_pathways.csv + convergence_leave_one_out.csv")


if __name__ == "__main__":
    main()
