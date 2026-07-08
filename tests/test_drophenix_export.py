"""Tests for the DroPhenix export bridge (schema correctness, never simulated values)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import drophenix_export as dpx  # noqa: E402


def test_external_line_id_pads_to_three_digits():
    assert dpx.external_line_id("line_21") == "DGRP_021"
    assert dpx.external_line_id("line_357") == "DGRP_357"
    assert dpx.external_line_id("DGRP_021") == "DGRP_021"


def test_external_line_id_rejects_no_digits():
    with pytest.raises(ValueError):
        dpx.external_line_id("no_number_here")


def test_climbing_template_matches_drophenix_schema_exactly():
    df = dpx.climbing_template("line_357", "F")
    assert list(df.columns) == dpx.CLIMBING_COLUMNS
    assert (df["Genotype"] == "DGRP_357").all()
    assert (df["Sex"] == "Female").all()
    # the measured value is never fabricated -- always blank for the user to fill in
    assert df["Count"].isna().all()
    assert not df["n"].isna().any()


def test_dams_template_matches_drophenix_schema_exactly():
    df = dpx.dams_template("line_21", "M")
    assert list(df.columns) == dpx.DAMS_COLUMNS
    assert (df["Genotype"] == "DGRP_021").all()
    assert (df["Sex"] == "Male").all()
    assert df["Activity"].isna().all()
    assert df["FlyID"].nunique() == 8


def test_template_for_trait_family_dispatch():
    c = dpx.template_for_trait_family("climbing", "line_5", "F")
    assert list(c.columns) == dpx.CLIMBING_COLUMNS
    s = dpx.template_for_trait_family("sleep", "line_5", "F")
    assert list(s.columns) == dpx.DAMS_COLUMNS


def test_unsupported_trait_family_raises_not_silently_guesses():
    # startle/starvation have no verified matching DroPhenix assay -- must fail loudly,
    # never silently fabricate a schema for them
    with pytest.raises(ValueError):
        dpx.template_for_trait_family("startle", "line_5", "F")
    with pytest.raises(ValueError):
        dpx.template_for_trait_family("starvation", "line_5", "F")
