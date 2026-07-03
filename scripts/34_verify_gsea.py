"""Verify the local GSEA implementation against known ground-truth properties
BEFORE trusting it on real data (gseapy could not be installed for a direct
cross-check, so we validate against the algorithm's defining properties instead).

Tests:
  1. POSITIVE CONTROL: a gene set planted entirely in the top of the ranked list
     must give a strong positive ES with p < 0.01.
  2. NEGATIVE CONTROL / CALIBRATION: gene sets with NO true association should
     give p-values that are approximately UNIFORM(0,1) across many trials (the
     defining property of a valid permutation test) -> checked via KS test
     against Uniform(0,1), and a coarse false-positive-rate check at alpha=0.05.
  3. RANK SENSITIVITY: a set planted in the middle of the list should give a
     smaller |ES| than one planted at the very top (monotonicity sanity check).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import gsea  # noqa: E402


def main() -> None:
    rng = np.random.default_rng(0)
    n_genes = 5000
    genes = [f"g{i}" for i in range(n_genes)]

    # ------------------------------------------------------------- test 1 ----
    scores = {g: v for g, v in zip(genes, rng.standard_normal(n_genes)[::-1] +
                                   np.linspace(3, -3, n_genes))}  # already trending
    top_set = set(genes[:100])                                   # top of ranked list
    r1 = gsea.gsea_preranked(scores, top_set, n_perm=2000, seed=1)
    print(f"[TEST 1: positive control] ES={r1['ES']:+.3f}  p={r1['p_value']:.4f}  "
          f"n_hit={r1['n_hit']}")
    assert r1["ES"] > 0.3 and r1["p_value"] < 0.01, "FAILED: planted top set not detected"
    print("  PASS: planted top-ranked gene set detected as strongly, significantly enriched.\n")

    # ------------------------------------------------------- test 2 (calib) --
    n_trials = 300
    null_scores = {g: v for g, v in zip(genes, rng.standard_normal(n_genes))}
    pvals = []
    for t in range(n_trials):
        random_set = set(rng.choice(genes, size=80, replace=False))
        r = gsea.gsea_preranked(null_scores, random_set, n_perm=300, seed=100 + t)
        pvals.append(r["p_value"])
    pvals = np.array(pvals)
    ks_stat, ks_p = stats.kstest(pvals, "uniform")
    fpr = float((pvals < 0.05).mean())
    print(f"[TEST 2: calibration] {n_trials} random (null) gene sets")
    print(f"  KS test vs Uniform(0,1): stat={ks_stat:.3f}, p={ks_p:.3f} "
          f"(p>0.05 means 'consistent with uniform' -> well-calibrated)")
    print(f"  empirical false-positive rate at alpha=0.05: {fpr:.3f} (expect ~0.05)")
    assert ks_p > 0.01, "FAILED: null p-values are not approximately uniform -> miscalibrated test"
    assert fpr < 0.15, "FAILED: false-positive rate far too high"
    print("  PASS: null p-values are approximately uniform -> test is calibrated.\n")

    # ------------------------------------------------------------- test 3 ----
    mid_set = set(genes[n_genes // 2 - 50: n_genes // 2 + 50])    # planted in the middle
    r3 = gsea.gsea_preranked(scores, mid_set, n_perm=1000, seed=3)
    print(f"[TEST 3: rank sensitivity] top-set ES={r1['ES']:+.3f} vs mid-set ES={r3['ES']:+.3f}")
    assert abs(r1["ES"]) > abs(r3["ES"]), "FAILED: top-ranked set should show larger |ES|"
    print("  PASS: enrichment score correctly decays for a set planted mid-list.\n")

    print("ALL VERIFICATION TESTS PASSED. The local GSEA implementation is trustworthy "
          "for use on the real BCM-DMAS data.")


if __name__ == "__main__":
    main()
