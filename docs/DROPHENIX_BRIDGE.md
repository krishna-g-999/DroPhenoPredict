# DroPhenix validation bridge — what it is, what was verified, what was deliberately left out

**Date added:** 2026-07-08 · **Revised:** 2026-08-08 (see §Revision) · **Code:**
`src/drophenopredict/drophenix_export.py`, `app.py` ("Validate with DroPhenix" page) ·
**Tests:** `tests/test_drophenix_export.py` (18 tests; suite total 59).

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

---

## Revision, 2026-08-08 — the schema was right and the *design* was wrong

An audit of both projects' calculations re-examined this bridge. The column names were still
correct: DroPhenix's `REQUIRED_LONG_COLS` and DAMS canonical columns are unchanged, and the
templates matched them exactly. But matching the **required** columns turned out not to be enough
to produce a comparison worth making. Three defects, each measured rather than argued:

**1. The climbing template did not name the clustering unit.** It shipped `Replicate`, which is not
one of DroPhenix's `VialID` aliases (`vial, vial_id, tube, tube_id, chamber, container, …`).
DroPhenix's `design_schema` therefore fell back to `Replicate` with the caveat "correct only if one
Replicate == one physical vial". In the authoring laboratory's own climbing data the vial-level ICC
is **0.289** (design effect 5.55), so 60 flies in 3 vials carry about 11 independent observations. A
predicted-vs-measured comparison built on the wrong denominator returns a falsely tight measured
value — it would appear to confirm or refute a prediction far more sharply than the data allow.
**Fixed:** `VialID` is now emitted, one per replicate, and constructed to stay unique when templates
for several lines/sexes/treatments are concatenated. Verified: DroPhenix now reports
`cluster_unit = VialID` on the exported file.

**2. The DAMS template used 30-minute bins, which cannot express the sleep definition.** Drosophila
sleep is ≥ 5 minutes of inactivity; DroPhenix implements it as `ceil(5 / bin_minutes)` empty bins,
which at a 30-minute bin is **one** bin — so a single empty bin became a 30-minute sleep bout. On a
simulated fly whose true sleep is 968 min/day:

| bin size | sleep recovered |
|---|---|
| 1 min | **968** (exact) |
| 5 min | 905 |
| 15 min | 780 |
| 30 min | **570** (−41%) |

The DGRP sleep values this tool predicts were measured on standard 1-minute DAM records, so the old
template would have systematically "disconfirmed" correct predictions. **Fixed:** `interval_min`
defaults to 1 and sub-minute values are refused. Verified: DroPhenix infers `bin_min = 1`.

**3. The DAMS template carried no `Monitor`/`Channel`, so the layout was invisible.** Monitors differ
in beam sensitivity and incubator position; a group nested inside one monitor is arithmetically
indistinguishable from a monitor effect, and that is unfixable after the run. **Fixed:** both columns
are emitted and flies are dealt across monitors.

**4. Found while testing the fix — concatenated templates put two flies in one tube.** Each template
numbered its channels from 1, so combining a template for one line with another produced two
different flies in monitor 1, channel 1. Every row still looked individually valid, so nothing
downstream would have caught it. **Fixed:** `combine_dams_templates()` re-deals the tubes. Its first
implementation then failed its own test in a second, subtler way: a plain round-robin over groups
aliases whenever the group count is a multiple of the monitor count, sending every fly of group 0 to
M1 and every fly of group 1 to M2 while still *looking* interleaved. The monitor index is now offset
by the group index as well as the position; the 2-group case is a parametrised regression test.

### Why the allocation logic is written here rather than imported

The obvious move — call DroPhenix's `dam_planner.plan_run` — would import GPL-3.0 code into this
project and force it to adopt that licence, which is precisely the entanglement §"What was checked
before building anything" exists to avoid. The rule implemented here (deal across monitors,
interleave positions, offset by group) is written from the design principle and is deliberately
simple. DroPhenix's planner remains the authority; it is run **on the exported file**, not imported.

### Cross-repository verification performed

The exported files were generated in this project's environment and read back in DroPhenix's, which
is the actual workflow rather than a mock of it:

- climbing template → `design_schema.design_report` reports `cluster_unit = VialID`, no blocking
  defects.
- four combined DAMS templates (2 lines × 2 sexes, 32 flies, 46,112 rows) →
  `dam_planner.check_plan` returns `usable: True`, `{'Genotype': 'crossed', 'Sex': 'crossed'}`,
  **no warnings and nothing blocking**; `dams.sleep_per_fly` reads `bin_min = 1` over 32 flies with
  8 in each group; genotypes split 8/8 across both monitors.
- `pytest tests/` — **59 passed** (41 prior + 18 bridge).

### Still not offered

No template for startle or starvation: DroPhenix has no verified matching assay module, and guessing
a schema for them would violate the no-fabrication rule. Unchanged from the original design.
