"""Regression tests for the local preranked-GSEA engine (src/drophenopredict/gsea.py).

These lock in the two properties verified in scripts/34_verify_gsea.py so a future
change cannot silently reintroduce the one-sided p-value bug that was caught and
fixed there (which capped p-values near 0.5 and inflated the false-positive rate).
Uses smaller n/perm counts than the full verification script for test speed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import gsea  # noqa: E402


def test_planted_top_set_is_detected():
    rng = np.random.default_rng(0)
    n = 2000
    genes = [f"g{i}" for i in range(n)]
    scores = {g: v for g, v in zip(genes, np.linspace(3, -3, n))}
    top_set = set(genes[:50])
    r = gsea.gsea_preranked(scores, top_set, n_perm=500, seed=1)
    assert r["ES"] > 0.3
    assert r["p_value"] < 0.01


def test_null_pvalues_are_calibrated_not_capped_at_half():
    """The bug this guards against: one-sided-by-observed-sign p-values are
    mechanically bounded near 0.5 under the null. A correct two-sided test
    should give p-values spanning up past 0.5 towards 1.0."""
    rng = np.random.default_rng(0)
    n = 1500
    genes = [f"g{i}" for i in range(n)]
    null_scores = {g: v for g, v in zip(genes, rng.standard_normal(n))}
    pvals = []
    for t in range(60):
        random_set = set(rng.choice(genes, size=60, replace=False))
        r = gsea.gsea_preranked(null_scores, random_set, n_perm=200, seed=100 + t)
        pvals.append(r["p_value"])
    pvals = np.array(pvals)
    assert pvals.max() > 0.6, (
        f"p-values capped at {pvals.max():.3f} — the one-sided p-value bug may have "
        "reappeared (see scripts/34_verify_gsea.py)"
    )
    fpr = float((pvals < 0.05).mean())
    assert fpr < 0.20, f"false-positive rate too high ({fpr:.3f}); test may be miscalibrated"


def test_rank_sensitivity():
    genes = [f"g{i}" for i in range(2000)]
    scores = {g: v for g, v in zip(genes, np.linspace(3, -3, 2000))}
    top_set = set(genes[:50])
    mid_set = set(genes[975:1025])
    r_top = gsea.gsea_preranked(scores, top_set, n_perm=200, seed=1)
    r_mid = gsea.gsea_preranked(scores, mid_set, n_perm=200, seed=2)
    assert abs(r_top["ES"]) > abs(r_mid["ES"])
