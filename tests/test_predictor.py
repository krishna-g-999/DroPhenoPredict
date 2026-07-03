"""Tests for the packaged DroPhenoPredict inference API."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

MODEL = Path("models/trained/drophenopredict_v1.0.json")
pytestmark = pytest.mark.skipif(not MODEL.exists(), reason="model not built (run scripts/10)")


def _pred():
    from drophenopredict.predictor import DroPhenoPredictor
    return DroPhenoPredictor("v1.0")


def test_line_id_normalization():
    p = _pred()
    assert p.normalize_line("DGRP_021") == "line_21"
    assert p.normalize_line("DGRP-21") == "line_21"
    assert p.normalize_line("357") == "line_357"
    assert p.normalize_line("line_75") == "line_75"


def test_prediction_has_calibrated_interval_and_flag():
    p = _pred()
    r = p.predict("line_21", "starvation_F")
    assert r["pi90_lower"] < r["predicted"] < r["pi90_upper"]
    assert isinstance(r["predictable"], bool)
    assert r["predictable"] is True            # starvation_F is a predictable model
    assert "unit" in r


def test_novel_genotype_is_refused_not_faked():
    p = _pred()
    r = p.predict("line_99999", "starvation_F")
    assert "error" in r and "predicted" not in r


def test_profile_covers_all_models():
    p = _pred()
    prof = p.predict_profile("line_357")
    assert len(prof["predictions"]) == len(p.available_traits())
