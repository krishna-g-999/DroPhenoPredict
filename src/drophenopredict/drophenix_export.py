"""Export DroPhenix-compatible experiment templates, for validating a DroPhenoPredict prediction.

DroPhenoPredict predicts trait values for DGRP lines from genotype/expression (no new experiment
needed). DroPhenix (github.com/krishna-g-999/DroPhenix, GPL-3.0) analyzes REAL experimental assay
data a researcher has already collected on any genotype. This module is the validation bridge
between them: given a DroPhenoPredict prediction, it generates an empty, correctly-schemaed CSV
template a researcher can fill in after running the real assay in DroPhenix — enabling a genuine
predicted-vs-measured comparison for the same DGRP line.

No DroPhenix code is imported or embedded here (this project does not adopt DroPhenix's GPL-3.0
license as a result) — only the column-name schemas below, read from DroPhenix's public source and
reproduced here as plain lists, not copied implementation:
  - modules/data_parser.py  REQUIRED_LONG_COLS = {"Genotype","Sex","Treatment","Replicate","Time","Count","n"}
                            DESIGN_KEY_COLS    includes VialID (optional, but see below)
  - modules/dams.py         canonical columns    Genotype, Treatment, Sex, FlyID, Time_min, Activity

Only traits with a verified, schema-compatible DroPhenix assay module are supported:
  - climbing -> DroPhenix's climbing / negative-geotaxis assay (data_parser.py)
  - sleep    -> DroPhenix's DAMS activity/sleep assay (dams.py)
startle and starvation have no verified matching DroPhenix assay module (DroPhenix's documented
assays are climbing, lifespan/survival, Y-maze, activity/DAMS, larval crawling, grooming) — no
export template is offered for those rather than guessing an unverified schema.

WHY THE TEMPLATES CARRY MORE THAN THE REQUIRED COLUMNS (revised 2026-08-08)
--------------------------------------------------------------------------
The required columns are the minimum DroPhenix will *accept*. They are not enough for the
comparison this bridge exists to support, for three measured reasons:

  VialID   DroPhenix clusters on the physical vial. `Replicate` is not one of its VialID
           aliases, so a table without VialID falls back to treating Replicate as the cluster
           with a warning, or -- worse -- treats every fly as independent. In the authoring
           laboratory's own climbing data the vial-level ICC is 0.289 (design effect 5.55), so
           60 flies in 3 vials carry about 11 independent observations. A predicted-vs-measured
           comparison built on the wrong denominator produces a falsely tight measured value.

  1-MINUTE Drosophila sleep is defined as >= 5 minutes of inactivity, and DroPhenix implements
  BINS     that as ceil(5 / bin_minutes) empty bins. At 30-minute bins that is ONE bin, so a
           single empty bin becomes a 30-minute sleep bout. Measured on a simulated fly whose
           true sleep is 968 min/day: 1-min bins recover 968, 5-min 905, 15-min 780, 30-min 570
           -- a 41% under-estimate. The DGRP sleep values this tool predicts were measured on
           standard 1-minute DAM records, so a 30-minute template would systematically
           "disconfirm" correct predictions. The default is now 1 minute.

  MONITOR/ Activity monitors differ in beam sensitivity and incubator position. If a whole
  CHANNEL  group sits in one monitor, the monitor effect and the group effect are the same
           column of numbers and no analysis separates them. Emitting Monitor and Channel lets
           DroPhenix's own design checks read the layout; without those columns they cannot.

The allocation rule below (alternate monitors, then fill channels) is written here from the
design principle, NOT ported from DroPhenix's planner -- importing that code would create
exactly the GPL-3.0 entanglement this bridge exists to avoid. It is deliberately simple; the
authoritative check remains DroPhenix's, run on the exported file.
"""
from __future__ import annotations

import pandas as pd

CLIMBING_COLUMNS = ["Genotype", "Sex", "Treatment", "Replicate", "VialID", "Time", "Count", "n"]
DAMS_COLUMNS = ["Genotype", "Treatment", "Sex", "FlyID", "Monitor", "Channel",
                "Time_min", "Activity"]

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


def vial_id(line: str, sex_label: str, treatment: str, replicate: int) -> str:
    """A vial label unique across lines, sexes and treatments, so several exported
    templates can be concatenated without two different vials colliding on one id."""
    safe = lambda s: "".join(ch if ch.isalnum() else "-" for ch in str(s))   # noqa: E731
    return f"{safe(line)}_{safe(sex_label)[:1]}_{safe(treatment)}_R{int(replicate)}"


def climbing_template(line_id: str, sex: str, *, n_replicates: int = 3,
                       timepoints_sec: tuple[int, ...] = (0, 5, 10, 15, 20),
                       n_flies_per_replicate: int = 20,
                       treatment: str = "Vehicle") -> pd.DataFrame:
    """Empty DroPhenix climbing-assay (negative geotaxis) template for one DGRP line.

    Carries DroPhenix's modules/data_parser.py REQUIRED_LONG_COLS plus `VialID`, which is the
    unit DroPhenix clusters on -- see the module docstring for why the required columns alone
    are not sufficient here. One Replicate is one physical vial, and VialID says so explicitly
    instead of leaving DroPhenix to infer it.

    `Count` is left blank for the user to fill in after running the real assay — this only lays
    out rows/columns, it never simulates or invents a measurement. Replicate count / timepoints
    / n are a common-default skeleton to edit to match the user's actual protocol.
    """
    line = external_line_id(line_id)
    sx = _sex_label(sex)
    rows = [
        {"Genotype": line, "Sex": sx, "Treatment": treatment, "Replicate": rep,
         "VialID": vial_id(line, sx, treatment, rep),
         "Time": t, "Count": None, "n": n_flies_per_replicate}
        for rep in range(1, n_replicates + 1) for t in timepoints_sec
    ]
    return pd.DataFrame(rows, columns=CLIMBING_COLUMNS)


def dams_channel_map(line_id: str, sex: str, *, n_flies: int = 8, n_monitors: int = 2,
                     n_channels: int = 32, treatment: str = "Vehicle") -> pd.DataFrame:
    """One row per fly: which monitor and channel it occupies.

    Flies are dealt alternately across monitors so that no group is nested inside one
    instrument. With a single line and sex there is nothing yet to confound, but templates for
    several lines are meant to be concatenated, and by then the layout is fixed and unfixable.
    Written from the design principle, not ported from DroPhenix (see module docstring).
    """
    if n_flies > n_monitors * n_channels:
        raise ValueError(f"{n_flies} flies exceed {n_monitors} monitors x {n_channels} "
                         f"channels = {n_monitors * n_channels} tubes.")
    line = external_line_id(line_id)
    rows = []
    next_channel = {m: 1 for m in range(1, n_monitors + 1)}
    for fly in range(1, n_flies + 1):
        mon = (fly - 1) % n_monitors + 1
        ch = next_channel[mon]
        next_channel[mon] += 1
        rows.append({"Genotype": line, "Treatment": treatment, "Sex": _sex_label(sex),
                     "FlyID": f"fly_{fly}", "Monitor": f"M{mon}", "Channel": ch})
    return pd.DataFrame(rows)


def combine_dams_templates(templates: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate per-line DAMS templates and RE-DEAL the tubes.

    Each template numbers its own channels from 1, so concatenating two of them puts two
    different flies in monitor 1 channel 1 -- a physically impossible layout that no
    downstream check would flag, because both rows look individually valid.

    This re-deals every fly across the available monitors in round-robin order of the groups,
    so each genotype/sex/treatment appears in every monitor rather than being nested inside
    one. That is the property that makes a monitor effect separable from a group effect; with
    a single template there is nothing to cross, which is why this only matters once several
    are combined. FlyID is rewritten to include the group, so ids are unique too.
    """
    if not templates:
        raise ValueError("No templates to combine.")
    per_fly, blocks = [], []
    for t in templates:
        missing = {"Genotype", "Treatment", "Sex", "FlyID", "Monitor", "Channel"} - set(t.columns)
        if missing:
            raise ValueError(f"Template is missing {sorted(missing)}; "
                             "combine_dams_templates expects dams_template output.")
        blocks.append(t)
        per_fly.append(t.drop_duplicates("FlyID")[["Genotype", "Treatment", "Sex", "FlyID"]])

    n_monitors = max(int(t["Monitor"].str.lstrip("M").astype(int).max()) for t in templates)
    groups = [g.reset_index(drop=True) for _, g in
              pd.concat(per_fly, ignore_index=True).groupby(
                  ["Genotype", "Treatment", "Sex"], sort=False)]

    # Interleave: take one fly from each group in turn, rather than emptying one group
    # before starting the next. Two things fall out of that. Monitors alternate, so no group
    # is nested in one instrument; and channels of the same group are non-adjacent, so a
    # temperature or light gradient across the rack cannot line up with a group boundary.
    # The monitor is offset by the group index as well as the position. Without that offset,
    # a straight round-robin aliases whenever the number of groups is a multiple of the
    # number of monitors -- with 2 groups and 2 monitors every fly of group 0 lands in M1 and
    # every fly of group 1 in M2, which is the exact confound this is meant to prevent.
    assign, next_ch = {}, {m: 1 for m in range(1, n_monitors + 1)}
    deepest = max(len(g) for g in groups)
    for i in range(deepest):
        for gi, g in enumerate(groups):
            if i >= len(g):
                continue
            row = g.iloc[i]
            mon = (i + gi) % n_monitors + 1
            key = (row["Genotype"], row["Treatment"], row["Sex"], row["FlyID"])
            assign[key] = (f"M{mon}", next_ch[mon])
            next_ch[mon] += 1

    out = []
    for t in blocks:
        t = t.copy()
        keys = list(zip(t["Genotype"], t["Treatment"], t["Sex"], t["FlyID"]))
        t["Monitor"] = [assign[k][0] for k in keys]
        t["Channel"] = [assign[k][1] for k in keys]
        t["FlyID"] = [f"{k[0]}_{str(k[2])[:1]}_{k[3]}" for k in keys]
        out.append(t)
    combined = pd.concat(out, ignore_index=True)[list(templates[0].columns)]

    tubes = combined.drop_duplicates("FlyID")[["Monitor", "Channel"]]
    if tubes.duplicated().any():                      # invariant, not decoration
        raise AssertionError("re-dealing produced two flies in one tube")
    return combined


def dams_template(line_id: str, sex: str, *, n_flies: int = 8,
                   duration_min: int = 1440, interval_min: int = 1,
                   n_monitors: int = 2, n_channels: int = 32,
                   treatment: str = "Vehicle") -> pd.DataFrame:
    """Empty DroPhenix DAMS activity/sleep template for one DGRP line.

    Carries DroPhenix's modules/dams.py canonical columns plus `Monitor` and `Channel`, so the
    layout is visible to a design check rather than lost.

    `interval_min` defaults to 1 because Drosophila sleep is defined as >= 5 minutes of
    inactivity: coarser bins cannot express that rule and under-count sleep badly (30-minute
    bins recovered 570 of a true 968 min/day in a simulated record). Raising it is allowed but
    will not produce sleep values comparable with the DGRP measurements this tool predicts.

    `Activity` is left blank for the user to fill in after running the real DAMS assay.
    """
    if interval_min < 1:
        raise ValueError("interval_min must be at least 1 minute.")
    cmap = dams_channel_map(line_id, sex, n_flies=n_flies, n_monitors=n_monitors,
                            n_channels=n_channels, treatment=treatment)
    times = range(0, duration_min + 1, interval_min)
    rows = [
        {**fly_row, "Time_min": t, "Activity": None}
        for fly_row in cmap.to_dict("records") for t in times
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
