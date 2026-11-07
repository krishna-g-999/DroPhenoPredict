"""Tests for the drug-v1.0 drug/toxin extension model and its input tables (Findings 21)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

MODEL = Path("models/trained/drophenopredict_drug-v1.0.json")
SCREEN = Path("data/processed/chemical_predictability.csv")
LONG = Path("data/processed/chemical_traits_long.csv")
pytestmark = pytest.mark.skipif(not MODEL.exists(), reason="drug-v1.0 model not built")


@pytest.fixture(scope="module")
def pred():
    from drophenopredict.predictor import DroPhenoPredictor
    return DroPhenoPredictor("drug-v1.0")


def test_loads_in_multimodal_format(pred):
    assert pred.mode == "multimodal"
    assert len(pred.available_traits()) == 12


def test_included_traits_match_screen_exactly(pred):
    """Model contents must be exactly the treated/effect traits with FDR q < 0.05."""
    s = pd.read_csv(SCREEN)
    passing = s[s["kind"].isin(["treated", "effect"])
                & ((s["q_geno"] < 0.05) | (s["q_adaptive"] < 0.05))]
    assert set(passing["trait_id"]) == set(pred.available_traits())
    # no control trait may be packaged as a drug response
    assert all(m["kind"] in ("treated", "effect") for m in pred.models.values())


def test_reported_accuracy_and_coverage(pred):
    for t, m in pred.models.items():
        assert m["fdr_q"] < 0.05, t
        assert m["cv_r"] > 0.2, t
        assert 0.85 <= m["interval_coverage_90"] <= 0.95, t


def test_intervals_clipped_to_valid_range(pred):
    for t, m in pred.models.items():
        if not m["bounds"]:
            continue
        lo_b, hi_b = m["bounds"]
        for line in m["lines"][:25]:
            r = pred.predict(line, t)
            assert lo_b <= r["pi90_lower"] <= r["pi90_upper"] <= hi_b, (t, line)


def test_interval_note_flags_uninformative_intervals(pred):
    for t, m in pred.models.items():
        assert ("interval_note" in m) == (m["interval_median_share_of_range"] >= 0.85), t
    r = pred.predict(pred.models["parathion_survival_48h_B"]["lines"][0], "parathion_survival_48h_B")
    assert "interval_note" in r


def test_prediction_carries_drug_fields(pred):
    r = pred.predict("DGRP_357", "paraquat_resistance_F")
    assert r["compound"] == "paraquat" and r["kind"] == "treated"
    assert r["modality"] == "exp_F"
    eff = pred.models["methylphenidate_food_intake_effect_M"]
    assert eff["kind"] == "effect" and eff["control_trait"] == "methylphenidate_food_intake_ctrl_M"


def test_out_of_panel_line_is_refused(pred):
    r = pred.predict("DGRP_999", "lead_viability_B")
    assert "error" in r and "predicted" not in r


def test_trait_table_is_clean():
    d = pd.read_csv(LONG)
    assert not d.duplicated(["trait_id", "line"]).any()
    panel = set(np.load("data/processed/grm.npz", allow_pickle=True)["iids"].tolist())
    assert set(d["line"]) <= panel
    assert np.isfinite(d["value"]).all()


def test_v11_unchanged_by_extension():
    from drophenopredict.predictor import DroPhenoPredictor
    v = DroPhenoPredictor("v1.1")
    r = v.predict("DGRP_357", "starvation_M")
    assert "compound" not in r and r["pi90_lower"] < r["predicted"] < r["pi90_upper"]
