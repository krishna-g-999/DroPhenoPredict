"""Tests for the v1.1 multi-modal DroPhenoPredict API."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

MODEL = Path("models/trained/drophenopredict_v1.1.json")
pytestmark = pytest.mark.skipif(not MODEL.exists(), reason="v1.1 model not built")


def _pred():
    from drophenopredict.predictor import DroPhenoPredictor
    return DroPhenoPredictor("v1.1")


def test_multimodal_mode_loads():
    p = _pred()
    assert p.mode == "multimodal"
    assert len(p.available_traits()) >= 6


def test_prediction_reports_modality_and_interval():
    p = _pred()
    r = p.predict("line_21", "starvation_M")
    assert r["modality"] in ("geno", "exp")
    assert r["pi90_lower"] < r["predicted"] < r["pi90_upper"]
    assert "model_cv_r" in r and isinstance(r["predictable"], bool)


def test_expression_trait_uses_expression_modality():
    p = _pred()
    # climbing_M is selected as 'exp' by nested CV (Findings 06)
    meta = p.models["climbing_M"]
    assert meta["deployed_modality"] == "exp"


def test_line_outside_panel_is_refused():
    p = _pred()
    r = p.predict("line_99999", "starvation_M")
    assert "error" in r and "predicted" not in r


def test_profile_covers_all_models():
    p = _pred()
    prof = p.predict_profile("line_357")
    assert len(prof["predictions"]) == len(p.available_traits())
