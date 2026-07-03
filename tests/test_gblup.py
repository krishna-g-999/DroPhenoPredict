"""Engine-correctness tests for the GBLUP pipeline + a guarded scientific fact.

Two things are validated here:
  1. The ENGINE is mathematically correct: on a kernel with genuine relatedness
     structure, REML recovers a known heritability and GBLUP cross-validation
     predicts strongly. (test_engine_*)
  2. The DGRP is an UNRELATED association panel: on the real GRM, even a highly
     heritable simulated trait is ~unpredictable out-of-sample (in-sample fit is
     high, CV R2 ~ 0). This is why real-startle CV R2 is at the noise floor — a
     data property, not a bug. (test_dgrp_unrelated_*)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import gblup, genotypes  # noqa: E402

GRM_PATH = Path("data/processed/grm.npz")


def _sim(K, true_h2, seed):
    """y = g + e with g ~ N(0, h2*Kc), e ~ N(0, 1-h2); Kc has unit avg diagonal."""
    rng = np.random.default_rng(seed)
    n = K.shape[0]
    Kc = K / np.mean(np.diag(K))
    L = np.linalg.cholesky(Kc + 1e-6 * np.eye(n))
    g = L @ rng.standard_normal(n)
    g *= np.sqrt(true_h2) / g.std()
    e = rng.standard_normal(n)
    e *= np.sqrt(1 - true_h2) / e.std()
    return g + e


def _structured_kernel(nfam=10, per=25, seed=0):
    """Synthetic GRM with real relatedness: lines cluster into 'families'."""
    rng = np.random.default_rng(seed)
    n = nfam * per
    anc = rng.standard_normal((n, 50))
    for f in range(nfam):
        anc[f * per:(f + 1) * per] += 3 * rng.standard_normal((1, 50))
    Z = (anc - anc.mean(0)) / anc.std(0)
    return Z @ Z.T / Z.shape[1]


# ----------------------------------------------------------- engine is correct
def test_engine_reml_recovers_h2_on_structured_kernel():
    K = _structured_kernel()
    for true_h2 in (0.3, 0.6, 0.8):
        ests = [gblup.reml(K, _sim(K, true_h2, s)).h2 for s in range(5)]
        assert abs(np.mean(ests) - true_h2) < 0.15, \
            f"h2={true_h2}: REML got {np.mean(ests):.3f}"


def test_engine_cv_tracks_heritability_on_structured_kernel():
    K = _structured_kernel()
    r2_noise = gblup.cross_validate(K, _sim(K, 0.001, 11), n_repeats=4).r2_per_repeat.mean()
    r2_signal = gblup.cross_validate(K, _sim(K, 0.8, 11), n_repeats=4).r2_per_repeat.mean()
    assert r2_noise < 0.05, f"noise CV R2 should be ~0, got {r2_noise:.3f}"
    assert r2_signal > 0.40, f"h2=0.8 CV R2 should be high, got {r2_signal:.3f}"


def test_multikernel_downweights_noise_kernel():
    """G+T REML: adding a pure-noise second kernel must not destroy prediction;
    its variance fraction should be small and CV r stay close to signal-only."""
    Ksig = _structured_kernel(seed=0)
    Knoise = _structured_kernel(seed=999)              # unrelated structure
    y = _sim(Ksig, 0.8, seed=4)
    vc = gblup.reml_multi([Ksig, Knoise], y)
    frac_noise = vc[1] / vc.sum()
    r_signal = gblup.cross_validate(Ksig, y, n_repeats=3).r_per_repeat.mean()
    r_multi = gblup.cross_validate_multi([Ksig, Knoise], y, n_repeats=3)["r_mean"]
    assert frac_noise < 0.30, f"noise kernel got variance frac {frac_noise:.2f}"
    assert r_multi > r_signal - 0.15, f"multi r {r_multi:.2f} vs signal {r_signal:.2f}"


def test_engine_conformal_calibrated_on_structured_kernel():
    K = _structured_kernel()
    conf = gblup.conformal_coverage(K, _sim(K, 0.5, 9), alpha=0.10, n_repeats=6)
    assert 0.82 <= conf["empirical_coverage"] <= 0.97


def test_engine_permutation_null_centered():
    K = _structured_kernel()
    perm = gblup.permutation_test(K, _sim(K, 0.5, 3), observed_r2=0.0, n_perm=50)
    # Out-of-sample CV R2 under no signal is ~0 or slightly negative (shuffled
    # labels do worse than the global mean on held-out folds) — never positive.
    assert -0.15 < perm["null_mean"] < 0.03


# ----------------------------------------- DGRP is an unrelated panel (science)
@pytest.mark.skipif(not GRM_PATH.exists(), reason="GRM not built")
def test_dgrp_grm_is_an_unrelated_panel():
    K, _, _ = genotypes.load_grm(GRM_PATH)
    Kc = K / np.mean(np.diag(K))
    off = Kc[~np.eye(len(Kc), dtype=bool)]
    assert abs(off.mean()) < 0.02 and off.std() < 0.06, \
        f"expected near-zero relatedness; mean={off.mean():.3f} sd={off.std():.3f}"


@pytest.mark.skipif(not GRM_PATH.exists(), reason="GRM not built")
def test_dgrp_high_h2_trait_is_unpredictable_out_of_sample():
    """The crux: in-sample fit is high but CV R2 ~ 0 on the real GRM."""
    K, _, _ = genotypes.load_grm(GRM_PATH)
    y = _sim(K, 0.8, 11)
    idx = np.arange(len(K))
    insample = gblup.gblup_predict(K, y, idx, idx)
    r2_in = 1 - np.sum((y - insample) ** 2) / np.sum((y - y.mean()) ** 2)
    r2_cv = gblup.cross_validate(K, y, n_repeats=4).r2_per_repeat.mean()
    assert r2_in > 0.30, f"in-sample fit unexpectedly low ({r2_in:.3f})"
    assert r2_cv < 0.10, f"DGRP CV R2 should be ~0 ({r2_cv:.3f})"
