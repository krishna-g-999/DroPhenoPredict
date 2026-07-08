"""Export DroPhenix-compatible experiment templates, for validating a DroPhenoPredict prediction.

DroPhenoPredict predicts trait values for DGRP lines from genotype/expression (no new experiment
needed). DroPhenix (github.com/krishna-g-999/DroPhenix, GPL-3.0) analyzes REAL experimental assay
data a researcher has already collected on any genotype. This module is the validation bridge
between them: given a DroPhenoPredict prediction, it generates an empty, correctly-schemaed CSV
template a researcher can fill in after running the real assay in DroPhenix — enabling a genuine
predicted-vs-measured comparison for the same DGRP line.

No DroPhenix code is imported or embedded here (this project does not adopt DroPhenix's GPL-3.0
license as a result) — only two column-name schemas below, read from DroPhenix's public source and
reproduced here as plain lists, not copied implementation:
  - modules/data_parser.py  REQUIRED_LONG_COLS = {"Genotype","Sex","Treatment","Replicate","Time","Count","n"}
  - modules/dams.py         canonical columns    Genotype, Treatment, Sex, FlyID, Time_min, Activity

Only traits with a verified, schema-compatible DroPhenix assay module are supported:
  - climbing -> DroPhenix's climbing / negative-geotaxis assay (data_parser.py)
  - sleep    -> DroPhenix's DAMS activity/sleep assay (dams.py)
startle and starvation have no verified matching DroPhenix assay module (DroPhenix's documented
assays are climbing, lifespan/survival, Y-maze, activity/DAMS, larval crawling, grooming) — no
export template is offered for those rather than guessing an unverified schema.
"""
from __future__ import annotations

import pandas as pd

CLIMBING_COLUMNS = ["Genotype", "Sex", "Treatment", "Replicate", "Time", "Count", "n"]
DAMS_COLUMNS = ["Genotype", "Treatment", "Sex", "FlyID", "Time_min", "Activity"]

SUPPORTED_TRAIT_FAMILIES = ("climbing", "sleep")

_SEX_LABEL = {"F": "Female", "M": "Male"}


def external_line_id(internal_line_id: str) -> str:
    """'line_357' -> 'DGRP_357'; 'line_21' -> 'DGRP_021' (Mackay et al. nomenclature)."""
    digits = "".join(ch for ch in str(internal_line_id) if ch.isdigit())
    if not digits:
        raise ValueError(f"Cannot parse a line number from '{internal_line_id}'")
    return f"DGRP_{int(digits):03d}"


def _sex_label(sex_suffix: str) -> str:
    return _SEX_LABEL.get(sex_suffix.upper(), sex_suffix)


def climbing_template(line_id: str, sex: str, *, n_replicates: int = 3,
                       timepoints_sec: tuple[int, ...] = (0, 5, 10, 15, 20),
                       n_flies_per_replicate: int = 20,
                       treatment: str = "Vehicle") -> pd.DataFrame:
    """Empty DroPhenix climbing-assay (negative geotaxis) template for one DGRP line.

    Matches DroPhenix's modules/data_parser.py REQUIRED_LONG_COLS exactly. `Count` is left blank
    for the user to fill in after running the real assay — this only lays out rows/columns, it
    never simulates or invents a measurement. Replicate count / timepoints / n are a common-default
    skeleton to edit to match the user's actual protocol, not a fixed requirement.
    """
    line = external_line_id(line_id)
    rows = [
        {"Genotype": line, "Sex": _sex_label(sex), "Treatment": treatment,
         "Replicate": rep, "Time": t, "Count": None, "n": n_flies_per_replicate}
        for rep in range(1, n_replicates + 1) for t in timepoints_sec
    ]
    return pd.DataFrame(rows, columns=CLIMBING_COLUMNS)


def dams_template(line_id: str, sex: str, *, n_flies: int = 8,
                   duration_min: int = 1440, interval_min: int = 30,
                   treatment: str = "Vehicle") -> pd.DataFrame:
    """Empty DroPhenix DAMS activity/sleep template for one DGRP line.

    Matches DroPhenix's modules/dams.py canonical columns exactly. `Activity` is left blank for
    the user to fill in after running the real DAMS assay.
    """
    line = external_line_id(line_id)
    rows = [
        {"Genotype": line, "Treatment": treatment, "Sex": _sex_label(sex),
         "FlyID": f"fly_{fly}", "Time_min": t, "Activity": None}
        for fly in range(1, n_flies + 1) for t in range(0, duration_min + 1, interval_min)
    ]
    return pd.DataFrame(rows, columns=DAMS_COLUMNS)


def template_for_trait_family(trait_family: str, line_id: str, sex: str) -> pd.DataFrame:
    """Dispatch to the right DroPhenix template by DroPhenoPredict trait-family name."""
    if trait_family == "climbing":
        return climbing_template(line_id, sex)
    if trait_family == "sleep":
        return dams_template(line_id, sex)
    raise ValueError(f"No verified DroPhenix assay module for trait '{trait_family}'. "
                      f"Supported: {SUPPORTED_TRAIT_FAMILIES}")
