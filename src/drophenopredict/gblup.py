"""GBLUP: genomic best linear unbiased prediction via a mixed model.

Model:  y = 1*mu + g + e,   g ~ N(0, sigma_g^2 * K),   e ~ N(0, sigma_e^2 * I)
where K is the genomic relationship matrix (genotypes.py).

This is the canonical quantitative-genetics method for predicting a quantitative
trait from genome-wide genotypes. We:
  - estimate genomic heritability h2 = sigma_g^2/(sigma_g^2+sigma_e^2) by REML,
  - predict held-out lines by BLUP,
  - validate by repeated k-fold CV (lines are the unit -> no leakage),
  - test significance by permutation,
  - quantify uncertainty by split-conformal prediction intervals.

All numbers are computed; nothing is assumed. delta = sigma_e^2/sigma_g^2, so
h2 = 1/(1+delta).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.linalg import cho_factor, cho_solve


# ---------------------------------------------------------------- alignment ---
def align(K: np.ndarray, iids: list[str], y) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Restrict GRM and target to their common lines, consistent order.

    y: pandas Series indexed by line id. Returns (K_sub, y_vec, lines).
    """
    common = [l for l in iids if l in set(y.index)]
    idx = [iids.index(l) for l in common]
    K_sub = K[np.ix_(idx, idx)]
    y_vec = y.loc[common].to_numpy(dtype=float)
    return K_sub, y_vec, common


# ---------------------------------------------------------------- REML --------
# Mixed model with fixed effects X (defaults to an intercept; pass extra columns
# for covariates such as inversion karyotype / Wolbachia status).
def _reml_neg_ll(log_delta: float, lam: np.ndarray, Uty: np.ndarray, UtX: np.ndarray) -> float:
    """REML negative log-likelihood as a function of log(delta), general X."""
    delta = np.exp(log_delta)
    inv = 1.0 / (lam + delta)
    XtViX = (UtX * inv[:, None]).T @ UtX        # q x q
    XtViy = (UtX * inv[:, None]).T @ Uty        # q
    beta = np.linalg.solve(XtViX, XtViy)
    r = Uty - UtX @ beta
    rVr = float(np.sum(r * r * inv))
    n, q = lam.size, UtX.shape[1]
    sigma_g2 = rVr / (n - q)
    _, logdet_XViX = np.linalg.slogdet(XtViX)
    return 0.5 * ((n - q) * np.log(sigma_g2) + np.sum(np.log(lam + delta)) + logdet_XViX)


@dataclass
class REMLFit:
    h2: float
    delta: float
    sigma_g2: float
    sigma_e2: float
    at_bound: bool = False     # True if delta hit an optimizer bound (h2 not literal)


def _as_X(X, n):
    return np.ones((n, 1)) if X is None else np.asarray(X, dtype=float).reshape(n, -1)


def reml(K: np.ndarray, y: np.ndarray, X: np.ndarray | None = None) -> REMLFit:
    """REML estimate of variance components / genomic heritability (fixed effects X)."""
    lo, hi = np.log(1e-5), np.log(1e5)
    Xm = _as_X(X, y.size)
    lam, U = np.linalg.eigh(K)
    lam = np.clip(lam, 0, None)
    Uty = U.T @ y
    UtX = U.T @ Xm
    res = minimize_scalar(_reml_neg_ll, bounds=(lo, hi), args=(lam, Uty, UtX),
                          method="bounded")
    at_bound = bool(res.x < lo + 0.05 or res.x > hi - 0.05)
    delta = float(np.exp(res.x))
    inv = 1.0 / (lam + delta)
    XtViX = (UtX * inv[:, None]).T @ UtX
    beta = np.linalg.solve(XtViX, (UtX * inv[:, None]).T @ Uty)
    r = Uty - UtX @ beta
    sigma_g2 = float(np.sum(r * r * inv) / (y.size - Xm.shape[1]))
    return REMLFit(h2=1.0 / (1.0 + delta), delta=delta,
                   sigma_g2=sigma_g2, sigma_e2=delta * sigma_g2, at_bound=at_bound)


# ---------------------------------------------------------------- BLUP --------
def gblup_predict(K: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray,
                  delta: float | None = None, X: np.ndarray | None = None) -> np.ndarray:
    """Predict test lines by BLUP with fixed effects X (default intercept)."""
    Ktt = K[np.ix_(train, train)]
    Kst = K[np.ix_(test, train)]
    yt = y[train]
    Xm = _as_X(X, y.size)
    Xtr, Xte = Xm[train], Xm[test]
    if delta is None:
        delta = reml(Ktt, yt, Xtr).delta
    Ainv = np.linalg.inv(Ktt + delta * np.eye(Ktt.shape[0]))
    beta = np.linalg.solve(Xtr.T @ Ainv @ Xtr, Xtr.T @ Ainv @ yt)   # GLS fixed effects
    resid = yt - Xtr @ beta
    return Xte @ beta + Kst @ (Ainv @ resid)


# ----------------------------------------------- multi-kernel mixed model -----
# y = 1*mu + sum_k u_k + e, u_k ~ N(0, sigma_k^2 K_k), e ~ N(0, sigma_e^2 I).
# Variance components by REML; used for multi-omic (genomic + transcriptomic)
# prediction (Morgante et al. 2020 approach). This is the FAIR test of whether
# transcriptome adds predictive value over genotype.

def _neg_reml_multi(log_vc, Ks, y, X):
    vc = np.exp(log_vc)
    n = y.size
    V = sum(v * K for v, K in zip(vc[:-1], Ks)) + vc[-1] * np.eye(n)
    try:
        c, low = cho_factor(V, lower=True)
    except np.linalg.LinAlgError:
        return 1e12
    Vinv_y = cho_solve((c, low), y)
    Vinv_X = cho_solve((c, low), X)
    XtVinvX = X.T @ Vinv_X
    beta = np.linalg.solve(XtVinvX, X.T @ Vinv_y)
    r = y - X @ beta
    Vinv_r = cho_solve((c, low), r)
    logdet_V = 2.0 * np.sum(np.log(np.diag(c)))
    _, logdet_XVX = np.linalg.slogdet(XtVinvX)
    return 0.5 * (logdet_V + logdet_XVX + r @ Vinv_r)


def reml_multi(Ks: list[np.ndarray], y: np.ndarray) -> np.ndarray:
    """REML variance components [sigma_1^2, ..., sigma_K^2, sigma_e^2] (>=0)."""
    var_y = float(np.var(y))
    start = np.log(np.full(len(Ks) + 1, var_y / (len(Ks) + 1) + 1e-6))
    res = minimize(_neg_reml_multi, start, args=(Ks, y, np.ones((y.size, 1))),
                   method="Nelder-Mead",
                   options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-4})
    return np.exp(res.x)


def predict_multi(Ks_full, y, train, test, vc):
    Ks_tt = [K[np.ix_(train, train)] for K in Ks_full]
    Ks_st = [K[np.ix_(test, train)] for K in Ks_full]
    yt = y[train]
    V = sum(v * K for v, K in zip(vc[:-1], Ks_tt)) + vc[-1] * np.eye(len(train))
    c, low = cho_factor(V, lower=True)
    one = np.ones(len(train))
    Vinv_y = cho_solve((c, low), yt)
    Vinv_1 = cho_solve((c, low), one)
    mu = (one @ Vinv_y) / (one @ Vinv_1)
    a = cho_solve((c, low), yt - mu)
    Kst = sum(v * K for v, K in zip(vc[:-1], Ks_st))
    return mu + Kst @ a


def cross_validate_multi(Ks: list[np.ndarray], y: np.ndarray, *, n_splits: int = 5,
                         n_repeats: int = 5, seed: int = 42) -> dict:
    """CV for the multi-kernel model; variance components refit per fold."""
    n = y.size
    r2s, rs = [], []
    for rep in range(n_repeats):
        rng = np.random.default_rng(seed + rep)
        folds = _kfold_indices(n, n_splits, rng)
        yt_all, yp_all = [], []
        for k in range(n_splits):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(n_splits) if j != k])
            vc = reml_multi([K[np.ix_(train, train)] for K in Ks], y[train])
            yp_all.append(predict_multi(Ks, y, train, test, vc))
            yt_all.append(y[test])
        yt = np.concatenate(yt_all); yp = np.concatenate(yp_all)
        r2s.append(1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2))
        rs.append(np.corrcoef(yt, yp)[0, 1])
    return {"R2_mean": float(np.mean(r2s)), "r_mean": float(np.mean(rs))}


# ---------------------------------------------------------------- CV ----------
@dataclass
class CVResult:
    y_true: np.ndarray            # pooled out-of-fold truth (one repeat's worth, concatenated)
    y_pred: np.ndarray
    r2_per_repeat: np.ndarray
    r_per_repeat: np.ndarray
    mae_per_repeat: np.ndarray
    rmse_per_repeat: np.ndarray
    null_mae_per_repeat: np.ndarray


def _kfold_indices(n: int, n_splits: int, rng: np.random.Generator) -> list[np.ndarray]:
    idx = rng.permutation(n)
    return [a for a in np.array_split(idx, n_splits)]


def cross_validate(K: np.ndarray, y: np.ndarray, *, n_splits: int = 5,
                   n_repeats: int = 10, seed: int = 42,
                   refit_delta_per_fold: bool = True,
                   X: np.ndarray | None = None) -> CVResult:
    n = y.size
    r2s, rs, maes, rmses, null_maes = [], [], [], [], []
    last_true, last_pred = None, None
    for rep in range(n_repeats):
        rng = np.random.default_rng(seed + rep)
        folds = _kfold_indices(n, n_splits, rng)
        yt_all, yp_all = [], []
        for k in range(n_splits):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(n_splits) if j != k])
            delta = None if refit_delta_per_fold else reml(K, y, X).delta
            pred = gblup_predict(K, y, train, test, delta=delta, X=X)
            yt_all.append(y[test]); yp_all.append(pred)
        yt = np.concatenate(yt_all); yp = np.concatenate(yp_all)
        ss_res = np.sum((yt - yp) ** 2)
        ss_tot = np.sum((yt - yt.mean()) ** 2)
        r2s.append(1 - ss_res / ss_tot)
        rs.append(np.corrcoef(yt, yp)[0, 1])
        maes.append(np.mean(np.abs(yt - yp)))
        rmses.append(np.sqrt(np.mean((yt - yp) ** 2)))
        null_maes.append(np.mean(np.abs(yt - yt.mean())))
        last_true, last_pred = yt, yp
    return CVResult(
        y_true=last_true, y_pred=last_pred,
        r2_per_repeat=np.array(r2s), r_per_repeat=np.array(rs),
        mae_per_repeat=np.array(maes), rmse_per_repeat=np.array(rmses),
        null_mae_per_repeat=np.array(null_maes),
    )


# ------------------------------------------------------- permutation test -----
def permutation_test(K: np.ndarray, y: np.ndarray, observed_r2: float, *,
                     n_perm: int = 200, n_splits: int = 5, seed: int = 1000) -> dict:
    """Permute y across lines; recompute single-pass CV R^2 to build the null."""
    delta = reml(K, y).delta   # fix delta for speed; conservative
    null = np.empty(n_perm)
    base_rng = np.random.default_rng(seed)
    for i in range(n_perm):
        rng = np.random.default_rng(seed + 1 + i)
        yp_perm = y[base_rng.permutation(y.size)]
        folds = _kfold_indices(y.size, n_splits, rng)
        yt_all, yhat_all = [], []
        for k in range(n_splits):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(n_splits) if j != k])
            yhat_all.append(gblup_predict(K, yp_perm, train, test, delta=delta))
            yt_all.append(yp_perm[test])
        yt = np.concatenate(yt_all); yh = np.concatenate(yhat_all)
        null[i] = 1 - np.sum((yt - yh) ** 2) / np.sum((yt - yt.mean()) ** 2)
    p = (1 + np.sum(null >= observed_r2)) / (1 + n_perm)
    return {"p_value": float(p), "null_mean": float(null.mean()),
            "null_q95": float(np.quantile(null, 0.95)), "null": null}


# --------------------------------------------------- conformal intervals ------
def conformal_coverage(K: np.ndarray, y: np.ndarray, *, alpha: float = 0.10,
                       n_splits: int = 5, n_repeats: int = 10, seed: int = 7) -> dict:
    """Split-conformal: empirical coverage & width of (1-alpha) prediction intervals.

    Per outer fold: split train -> proper-fit + calibration; quantile of
    calibration |residuals| sets the half-width; measure coverage on test.
    """
    covers, widths = [], []
    for rep in range(n_repeats):
        rng = np.random.default_rng(seed + rep)
        folds = _kfold_indices(y.size, n_splits, rng)
        for k in range(n_splits):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(n_splits) if j != k])
            rng.shuffle(train)
            cut = int(0.7 * train.size)
            fit, calib = train[:cut], train[cut:]
            delta = reml(K[np.ix_(fit, fit)], y[fit]).delta
            cal_pred = gblup_predict(K, y, fit, calib, delta=delta)
            resid = np.abs(y[calib] - cal_pred)
            q_level = min(1.0, np.ceil((calib.size + 1) * (1 - alpha)) / calib.size)
            q = np.quantile(resid, q_level)
            test_pred = gblup_predict(K, y, fit, test, delta=delta)
            covers.append(np.mean(np.abs(y[test] - test_pred) <= q))
            widths.append(2 * q)
    return {"target_coverage": 1 - alpha,
            "empirical_coverage": float(np.mean(covers)),
            "mean_interval_width": float(np.mean(widths))}
