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
    # DroPhenix's REQUIRED_LONG_COLS must all still be present
    assert {"Genotype", "Sex", "Treatment", "Replicate", "Time", "Count", "n"} <= set(df.columns)
    assert (df["Genotype"] == "DGRP_357").all()
    assert (df["Sex"] == "Female").all()
    # the measured value is never fabricated -- always blank for the user to fill in
    assert df["Count"].isna().all()
    assert not df["n"].isna().any()


def test_climbing_template_names_the_clustering_unit():
    """DroPhenix clusters on VialID; 'Replicate' is not one of its VialID aliases.

    Without this column the vial structure is invisible and every fly is treated as an
    independent observation. At the measured climbing ICC of 0.289 that overstates the
    effective sample size roughly five-fold, which would make a predicted-vs-measured
    comparison look far tighter than the data support.
    """
    df = dpx.climbing_template("line_357", "F", n_replicates=3)
    assert "VialID" in df.columns
    assert df["VialID"].notna().all()
    # one VialID per replicate, and it does not vary within a replicate
    assert df["VialID"].nunique() == 3
    assert (df.groupby("Replicate")["VialID"].nunique() == 1).all()


def test_vial_ids_do_not_collide_across_lines_sexes_or_treatments():
    """Templates are meant to be concatenated; two different vials must never share an id."""
    frames = [dpx.climbing_template(line, sex, treatment=tr)
              for line in ("line_21", "line_357")
              for sex in ("F", "M")
              for tr in ("Vehicle", "Paraquat")]
    ids = set()
    for f in frames:
        new = set(f["VialID"])
        assert not (ids & new), "VialID collision between templates"
        ids |= new
    assert len(ids) == 2 * 2 * 2 * 3


def test_dams_template_matches_drophenix_schema_exactly():
    df = dpx.dams_template("line_21", "M", duration_min=60)
    assert list(df.columns) == dpx.DAMS_COLUMNS
    assert {"Genotype", "Treatment", "Sex", "FlyID", "Time_min", "Activity"} <= set(df.columns)
    assert (df["Genotype"] == "DGRP_021").all()
    assert (df["Sex"] == "Male").all()
    assert df["Activity"].isna().all()
    assert df["FlyID"].nunique() == 8


def test_dams_template_bins_at_one_minute_so_the_sleep_rule_is_expressible():
    """Sleep is >= 5 minutes of inactivity, implemented as ceil(5/bin) empty bins.

    At a 30-minute bin that threshold is ONE bin, so a single empty bin becomes a
    30-minute bout: a simulated fly with 968 min/day of true sleep scored 570 at 30-minute
    bins. The DGRP sleep values this tool predicts come from 1-minute records, so anything
    coarser is not comparable.
    """
    df = dpx.dams_template("line_21", "M", duration_min=120)
    step = df.loc[df["FlyID"] == "fly_1", "Time_min"].diff().dropna().unique()
    assert list(step) == [1], f"bin size is {step}, must be 1 minute"


def test_dams_template_records_monitor_and_channel():
    """Without these columns a design check cannot see which monitor a fly sat in, and a
    group nested inside one instrument is indistinguishable from an instrument effect."""
    df = dpx.dams_template("line_21", "M", n_flies=8, duration_min=10)
    assert {"Monitor", "Channel"} <= set(df.columns)
    per_fly = df.drop_duplicates("FlyID")
    # every fly has exactly one monitor/channel, and no two flies share a tube
    assert (df.groupby("FlyID")[["Monitor", "Channel"]].nunique() == 1).all().all()
    assert not per_fly.duplicated(["Monitor", "Channel"]).any()
    # flies are split across monitors, not all loaded into one
    assert per_fly["Monitor"].nunique() == 2
    counts = per_fly["Monitor"].value_counts()
    assert counts.max() - counts.min() <= 1


def test_dams_template_refuses_more_flies_than_tubes():
    with pytest.raises(ValueError):
        dpx.dams_template("line_21", "M", n_flies=100, n_monitors=2, n_channels=32,
                          duration_min=10)


def test_dams_template_rejects_sub_minute_bins():
    with pytest.raises(ValueError):
        dpx.dams_template("line_21", "M", interval_min=0, duration_min=10)


def test_combining_templates_never_puts_two_flies_in_one_tube():
    """Each template numbers channels from 1, so a naive concatenation is physically
    impossible -- two flies in monitor 1, channel 1 -- and each row still looks valid on
    its own, so nothing downstream would catch it."""
    import pandas as pd
    a = dpx.dams_template("line_357", "F", duration_min=5)
    b = dpx.dams_template("line_21", "F", duration_min=5)

    # a fly is identified by its group AND its id -- 'fly_1' exists in both templates
    fly_key = ["Genotype", "Treatment", "Sex", "FlyID"]
    naive = pd.concat([a, b], ignore_index=True).drop_duplicates(fly_key)
    assert naive.duplicated(["Monitor", "Channel"]).any(), \
        "the collision this test guards against no longer occurs by construction"

    combined = dpx.combine_dams_templates([a, b])
    tubes = combined.drop_duplicates("FlyID")[["Monitor", "Channel"]]
    assert not tubes.duplicated().any()
    assert combined["FlyID"].nunique() == 16
    assert len(combined) == len(a) + len(b)


@pytest.mark.parametrize("lines,sexes", [
    (("line_357", "line_21"), ("F",)),                     # 2 groups, 2 monitors -- aliases
    (("line_357", "line_21"), ("F", "M")),                 # 4 groups
    (("line_357", "line_21", "line_45"), ("F", "M")),      # 6 groups, odd count
])
def test_combining_templates_crosses_every_factor_with_monitor(lines, sexes):
    """A group nested inside one monitor is indistinguishable from a monitor effect.

    The 2-group case is the one that matters: a plain round-robin over groups aliases
    whenever the group count is a multiple of the monitor count, sending every fly of
    group 0 to M1 and every fly of group 1 to M2 while still looking interleaved.
    """
    ts = [dpx.dams_template(l, s, duration_min=5) for l in lines for s in sexes]
    combined = dpx.combine_dams_templates(ts)
    per_fly = combined.drop_duplicates("FlyID")
    for factor in ("Genotype", "Sex"):
        if per_fly[factor].nunique() < 2:
            continue
        table = per_fly.groupby([factor, "Monitor"]).size().unstack(fill_value=0)
        assert (table > 0).all().all(), f"{factor} is confined to one monitor:\n{table}"
        assert (table.max(axis=1) - table.min(axis=1) <= 1).all(), \
            f"{factor} is unevenly split across monitors:\n{table}"


def test_combining_templates_does_not_block_groups_on_adjacent_channels():
    """Contiguous channels put a group on one part of the rack, where a temperature or
    light gradient lines up with the group boundary."""
    ts = [dpx.dams_template(l, "F", duration_min=5) for l in ("line_357", "line_21")]
    per_fly = dpx.combine_dams_templates(ts).drop_duplicates("FlyID")
    for (_, _), sub in per_fly.groupby(["Genotype", "Monitor"]):
        ch = sorted(sub["Channel"])
        assert not (len(ch) > 1 and ch[-1] - ch[0] + 1 == len(ch)), \
            f"a genotype occupies the contiguous block {ch}"


def test_combine_rejects_frames_that_are_not_templates():
    import pandas as pd
    with pytest.raises(ValueError):
        dpx.combine_dams_templates([pd.DataFrame({"Genotype": ["x"]})])
    with pytest.raises(ValueError):
        dpx.combine_dams_templates([])


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
