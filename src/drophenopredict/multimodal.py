"""Multi-modal, per-trait-adaptive prediction (DroPhenoPredict v1.1).

Per trait x sex, the predictor chooses the best DATA MODALITY:
  - geno : genomic relationship (GRM, additive genotype)
  - exp  : RNA-seq expression relationship (linear kernel, GSE117850)
  - combined : two-component REML of both (Morgante 2020 style)
The choice is made by NESTED cross-validation (inner CV selects the modality on
training data only; outer CV evaluates) so the reported accuracy of the adaptive
strategy is unbiased — the selection never sees the test fold.

Rationale (Findings 06): RNA-seq expression beats genotype for some traits
(climbing, starvation male); genotype wins for others (startle, sleep). No single
modality dominates, so adaptive selection is the principled design. Deep learning
/ non-linear kernels were tested and did not help (Findings 06) and are excluded.
"""
from __future__ import annotations

import numpy as np

from . import gblup


def _norm(K: np.ndarray) -> np.ndarray:
    return K / np.mean(np.diag(K))


def expression_kernel(expr_df, lines: list[str]) -> np.ndarray:
    """Linear (additive) expression relationship over `lines`, mean-diag 1."""
    X = expr_df[lines].to_numpy(float).T            # lines x genes
    keep = np.isfinite(X).all(0) & (X.std(0) > 1e-9)
    X = X[:, keep]
    X = (X - X.mean(0)) / X.std(0)
    return _norm((X @ X.T) / X.shape[1])


def build_candidates(Kg_full, g_iids, expr_df, y_series):
    """Common lines + genomic and expression kernels + aligned y, for one trait."""
    common = [l for l in g_iids if l in set(expr_df.columns) and l in set(y_series.index)]
    gi = [g_iids.index(l) for l in common]
    Kg = _norm(Kg_full[np.ix_(gi, gi)])
    Ke = expression_kernel(expr_df, common)
    yv = y_series.loc[common].to_numpy(float)
    return common, {"geno": Kg, "exp": Ke}, yv


# ------------------------------------------------------- nested modality CV ----
def _inner_select(kernels: dict, y: np.ndarray, train: np.ndarray, seed: int) -> str:
    """Pick the modality with best inner-CV r, using ONLY training lines."""
    best, best_r = None, -np.inf
    yt = y[train]
    for name, K in kernels.items():
        Ktt = K[np.ix_(train, train)]
        r = gblup.cross_validate(Ktt, yt, n_splits=5, n_repeats=2, seed=seed).r_per_repeat.mean()
        if r > best_r:
            best, best_r = name, r
    return best


def nested_modality_cv(kernels: dict, y: np.ndarray, *, n_splits: int = 5,
                       n_repeats: int = 10, seed: int = 42) -> dict:
    """Unbiased CV of the adaptive (inner-selected) modality + selection counts."""
    n = y.size
    rs, r2s = [], []
    selections: dict[str, int] = {k: 0 for k in kernels}
    yt_pool_last, yp_pool_last = None, None
    for rep in range(n_repeats):
        rng = np.random.default_rng(seed + rep)
        folds = gblup._kfold_indices(n, n_splits, rng)
        yt_all, yp_all = [], []
        for k in range(n_splits):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(n_splits) if j != k])
            choice = _inner_select(kernels, y, train, seed=seed + 100 + rep)
            selections[choice] += 1
            pred = gblup.gblup_predict(kernels[choice], y, train, test)
            yt_all.append(y[test]); yp_all.append(pred)
        yt = np.concatenate(yt_all); yp = np.concatenate(yp_all)
        rs.append(np.corrcoef(yt, yp)[0, 1])
        r2s.append(1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2))
        yt_pool_last, yp_pool_last = yt, yp
    return {"adaptive_r": float(np.mean(rs)), "adaptive_r_sd": float(np.std(rs)),
            "adaptive_R2": float(np.mean(r2s)), "selections": selections,
            "y_true": yt_pool_last, "y_pred": yp_pool_last}


# ------------------------------------------------- per-modality fixed CV r ------
def fixed_modality_r(kernels: dict, y: np.ndarray, *, seed: int = 42,
                     n_repeats: int = 10) -> dict:
    out = {}
    for name, K in kernels.items():
        out[name] = float(gblup.cross_validate(K, y, n_repeats=n_repeats, seed=seed)
                          .r_per_repeat.mean())
    # combined (two-kernel REML)
    out["combined"] = gblup.cross_validate_multi(list(kernels.values()), y,
                                                 n_repeats=5, seed=seed)["r_mean"]
    return out


def deployed_modality(kernels: dict, y: np.ndarray, seed: int = 42) -> str:
    """Modality chosen on ALL data (what the shipped model uses)."""
    r = {name: float(gblup.cross_validate(K, y, n_repeats=10, seed=seed).r_per_repeat.mean())
         for name, K in kernels.items()}
    return max(r, key=r.get)


# ------------------------------------------------------------- LOO + intervals --
def loo_predict(K: np.ndarray, y: np.ndarray, delta: float) -> np.ndarray:
    n = len(y); idx = np.arange(n); pred = np.empty(n)
    for i in range(n):
        pred[i] = gblup.gblup_predict(K, y, idx[idx != i], np.array([i]), delta=delta)[0]
    return pred
