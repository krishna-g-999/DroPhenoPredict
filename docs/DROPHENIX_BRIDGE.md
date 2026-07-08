# DroPhenix validation bridge — what it is, what was verified, what was deliberately left out

**Date added:** 2026-07-08 · **Code:** `src/drophenopredict/drophenix_export.py`,
`app.py` ("Validate with DroPhenix" page) · **Tests:** `tests/test_drophenix_export.py` (6 tests).

## The question that prompted this

The user asked whether DroPhenix (a separate, existing Drosophila phenotyping tool by the same
author, github.com/krishna-g-999/DroPhenix) and DroPhenoPredict should be merged into one module.

## What was checked before building anything (no assumption)

1. **License.** DroPhenix's GitHub repository page states it is licensed **GPL-3.0**. GPL-3.0 is
   copyleft: importing or embedding its code into DroPhenoPredict would require DroPhenoPredict to
   adopt GPL-3.0 too, propagating to anyone downstream. This is a real, largely irreversible
   licensing decision, not a technical detail — flagged to the user before any code was written.
2. **Data model.** DroPhenoPredict predicts trait values for the 205 fixed DGRP reference lines from
   genotype/expression (no new experiment needed). DroPhenix analyzes real experimental assay data a
   researcher has already collected, on any genotype they name. These are different questions
   ("what would this line's value likely be" vs. "what did I just measure"); a literal merge risks
   conflating **predicted** and **measured** values in one interface — the exact overclaiming risk
   this project has spent multiple sessions eliminating (see `AUDIT.md` §5–6).
3. **User decision.** Given the above, the user was asked to choose between a full code-level merge
   (GPL-3.0 adoption required) and a lightweight, license-free validation bridge. **They chose the
   bridge.**

## What was verified about DroPhenix's actual source (not the README) before building the bridge

Read directly from `github.com/krishna-g-999/DroPhenix`:
- `modules/data_parser.py` — canonical climbing/negative-geotaxis schema:
  `REQUIRED_LONG_COLS = {"Genotype", "Sex", "Treatment", "Replicate", "Time", "Count", "n"}`.
  `Genotype` is confirmed free-text with no fixed strain/DGRP vocabulary check.
- `modules/dams.py` — canonical DAMS activity/sleep schema: `Genotype, Treatment, Sex, FlyID,
  Time_min, Activity` (with documented column aliases). Also free-text `Genotype`.
- DroPhenix's documented assay set (climbing, lifespan/survival, Y-maze, activity/DAMS, larval
  crawling, grooming) has **no verified module for startle response or starvation resistance** —
  two of DroPhenoPredict's four shipped trait families. No template is offered for those traits.
  This is a stated, honest gap, not a silent omission — guessing an unverified schema for them
  would have violated the project's no-fabrication rule.

## What the bridge actually does

`drophenix_export.py` reproduces only the two verified **column-name lists** above (not any
DroPhenix implementation code) and generates an empty, correctly-shaped CSV template for a chosen
DGRP line/trait/sex: `Genotype` pre-filled with the external line ID (`DGRP_357`, converted from the
internal `line_357` id — same nomenclature as Mackay et al. 2012), the measured-value column
(`Count` or `Activity`) always blank, and reasonable default skeleton values (replicate count,
timepoints, fly count) the user is told to edit to match their own protocol.

The app's new "Validate with DroPhenix" page shows DroPhenoPredict's prediction and 90% interval
**next to** the download button, for reference — but that value is never written into the exported
file. The intended workflow: predict a line/trait here → export the template → run the real assay →
analyze it in DroPhenix → compare DroPhenix's real, measured output against this tool's prediction.
This is a genuine external-validation loop this project has repeatedly noted as missing (e.g.
`AUDIT.md` §7 gap 11 on the missing behavioural ground-truth for the BCM-DMAS panel is a different,
unrelated gap — this bridge does not close that one; it creates a *new* validation path for the core
DGRP predictor specifically).

## What this deliberately does NOT do

- No DroPhenix code, model, or data is imported, embedded, or redistributed.
- No claim that DroPhenix has been run, or that any validation has actually occurred — this ships a
  blank template generator, not a result.
- No template offered for startle/starvation (unverified DroPhenix compatibility).
- No change to DroPhenoPredict's license (remains unlicensed/all-rights-reserved by default, as
  before — a separate, later decision for the user to make if desired).

## Verification performed before shipping

- 6 new unit tests (`tests/test_drophenix_export.py`): schema-exact column matching against the two
  verified DroPhenix schemas, line-ID nomenclature conversion, blank-measured-column invariant (the
  test suite explicitly asserts `Count`/`Activity` are never populated), and that unsupported trait
  families raise loudly rather than silently guessing.
- Full app boot-tested in a live browser for both supported traits (climbing, sleep): confirmed the
  rendered table's columns and values match the verified schemas exactly, and that the displayed
  DroPhenoPredict prediction matches the already-verified Predict-page value for the same line/trait.
- Full test suite re-run after the addition: 47/47 passing (41 prior + 6 new).
