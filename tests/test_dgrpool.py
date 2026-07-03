"""Tests for DGRPool line-ID harmonization and phenotype aggregation.

to_grm_line_id is load-bearing: every join between DGRPool phenotypes and the
GRM/expression matrices depends on it. A silent bug here would misalign labels
without raising an error (a false-negative-shaped bug), so it gets a direct test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import dgrpool  # noqa: E402


def test_to_grm_line_id_strips_leading_zeros():
    assert dgrpool.to_grm_line_id("DGRP_021") == "line_21"
    assert dgrpool.to_grm_line_id("DGRP_357") == "line_357"
    assert dgrpool.to_grm_line_id("DGRP_005") == "line_5"


def test_to_grm_line_id_handles_double_underscore():
    # Wolbachia sheet uses 'DGRP__21' (double underscore) — must still resolve
    assert dgrpool.to_grm_line_id("DGRP__21") == "line_21"


def test_to_grm_line_id_zero_edge_case():
    assert dgrpool.to_grm_line_id("DGRP_000") == "line_0"


def test_line_means_averages_replicates_and_harmonizes():
    df = pd.DataFrame({
        "DGRP_line": ["DGRP_021", "DGRP_021", "DGRP_357"],
        "sex": ["F", "F", "F"],
        "value": [10.0, 20.0, 5.0],
    })
    s = dgrpool.line_means(df, sex="F", harmonize=True)
    assert s.loc["line_21"] == 15.0        # mean of replicates
    assert s.loc["line_357"] == 5.0
    assert set(s.index) == {"line_21", "line_357"}


def test_line_means_sex_filter():
    df = pd.DataFrame({
        "DGRP_line": ["DGRP_021", "DGRP_021"],
        "sex": ["F", "M"],
        "value": [10.0, 999.0],
    })
    s = dgrpool.line_means(df, sex="F", harmonize=True)
    assert s.loc["line_21"] == 10.0
