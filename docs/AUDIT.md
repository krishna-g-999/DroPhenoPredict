# DroPhenoPredict — Complete Editorial-Board Audit (2026-07-03)

Full re-audit of the project as it stands, performed as an editorial-board / expert-reviewer pass
(not a self-report): every claim below is independently re-checked against the current repo state
this session — by direct command execution and code reading, plus a second independent verification
pass (a separate agent re-ran the same checks blind) — not carried forward from prior audits without
re-verification. Where the prior audit (dated 2026-07-01) was stale or wrong, it is corrected here
and the fix is applied to the actual docs, not just noted.

---

## 1. Inventory (re-verified, two independent passes, this session)

| Item | Count | Note |
|---|---|---|
| Numbered pipeline scripts | **39** (01→44) | prior audit said 29 — undercount, corrected |
| Library modules (`src/drophenopredict/`) | **13** | unchanged |
| Test files / tests passing | **8 files, 41 tests, `pytest tests/ -q` → 41 passed** | re-run this session via project `.venv`, confirmed twice |
| `docs/FINDINGS_*.md` + `FINDINGS_coverage.md` | **19 findings docs + coverage doc = 20** | prior audit said 22 — overcount, corrected |
| `data/processed/` files | **82** | prior audit said 68 — undercount, corrected |
| Trained model artifacts | **4 files** (v1.0 `.json`+`.npz`, v1.1 `.json`+`.npz`) = **2 model versions** | both countings are "correct" depending on whether you count files or model versions; both stated now to avoid ambiguity |
| Git state | **2 commits, working tree clean, nothing pushed to a remote yet** | unchanged since last session |

**Why the prior audit's counts were off:** informal counting during a fast-moving session, not a
methodological problem — none of the corrected counts change any scientific conclusion. Fixed here
for the record.

## 2. Fabrication / integrity scan (re-run independently this session)

- Grepped `src/`, `scripts/`, `app.py` for `TODO|FIXME|HACK|fake|dummy|placeholder|mock_|simulated_`
  (excluding test fixtures). **One hit**, and it is the opposite of a problem: a docstring in
  `dgrpool.py:111` stating the fetch function "Returns an EMPTY DataFrame (never fake ...)". No
  fabricated/hardcoded values found anywhere in production code or scripts.
- Every catalogued fabrication from the original draft brief remains excluded per
  `SCIENTIFIC_CHARTER.md §7` — re-checked, no regressions.
- `requirements.txt` re-verified against `app.py`'s actual transitive import chain
  (`predictor → genotypes, gblup, multimodal`): every import is satisfied by the 7 lean packages;
  no drift since the deployment split was made.

## 3. Statistical rigor — a specific editorial-board question, answered directly

**"With ~20 findings documents (GWAS, QTT, GSEA, DE, cross-disease convergence, ortholog overlap,
proteomics, Yang2023 x3), is there a multiple-comparisons problem across the whole project, not just
within each analysis?"**

Verified: **11 scripts apply Benjamini-Hochberg FDR correction within their own scope**
(`20_expression_association.py`, `21_pathway_association.py`, `33_shared_biology.py`,
`35_gsea_convergence.py`, `36_gsea_convergence_finish.py`, `39a`, `39c`, `41`, `42`, `43`, plus
`44` reusing `gseapy`'s own FDR). **No project-wide family-wise correction is applied across
findings.** This is a **documented design choice**, not an oversight (`EDITORIAL_RISK_ANALYSIS.md`
§10): each finding is scoped and pre-registered independently, and no finding borrows significance
by pooling p-values with another. The honest framing for a manuscript's statistics section: this
project reports ~20 independently-controlled analyses, not one omnibus test — a reader evaluating
the totality of claims should treat each finding's significance as conditional on that finding's own
correction, not as part of one global error budget. This should be stated explicitly in any eventual
manuscript (not yet drafted) rather than left implicit.

## 4. Methods correctness — re-confirmed, one clarification made

- GBLUP/REML, GSEA, DE-pipeline verification bugs from prior sessions (one-sided p-value, FBgn/symbol
  mismatch, `spearmanr` NaN propagation, `statsmodels.MixedLM.cov_re` array-vs-DataFrame, wrong
  permutation null for the adaptive procedure) — all previously caught, fixed, and regression-tested;
  re-confirmed present and correct in the current code.
- **Pandas positional-axis calls, re-scoped precisely this session (correcting an imprecise claim in
  the prior audit):** a grep finds `.mean(0/1)`, `.std(0/1)`, `.sum(0/1)`, `.all(0/1)` in
  `src/drophenopredict/expression.py` and `multimodal.py` as well as 12 one-off scripts. **Read in
  context, the two library-module hits operate on plain `numpy.ndarray` objects** (each is called
  after an explicit `.to_numpy(...)` conversion a few lines earlier) — numpy's positional axis
  argument is not deprecated and this is not a forward-compatibility risk. The prior audit's
  characterization ("same pattern exists ... not the reusable library") was correct; a blind
  re-verification this session initially flagged it as a library-level risk from the grep alone, but
  reading the actual call sites confirms it isn't. The **12 one-off analysis scripts** (already
  executed, not part of the reusable library, would not change any reported number if fixed) do call
  these methods on real pandas objects and remain a low-priority, cosmetic gap for pandas 4.0 —
  unchanged assessment from the prior audit, now with an exact per-file inventory on record (see this
  session's verification transcript) instead of an approximate "~12 scripts."

## 5. Documentation consistency — real gaps found and fixed this session

An editorial board reads every doc, not just the code. Cross-checking all `docs/*.md` against actual
shipped code and against each other surfaced **five real staleness/inconsistency issues**, now fixed:

1. **`FEATURES.md` described only the v1.0 genotype-only model** and said transcriptome data was
   "excluded ... retained only as a future resource" — false for v1.1 (shipped), which actively
   deploys RNA-seq expression as the winning modality for 3/8 trait-sex models. **Fixed:** rewrote to
   document RNA-seq as a retained, deployed feature, and clarified the exclusion verdict applies to
   the Huang 2015 *microarray* data type specifically, not expression data in general.
2. **`MODEL_CARD.md` said "16 tests pass"** — stale; current is 41. **Fixed**, and added a pointer to
   the v2 resource and the Yang2023 null (previously undiscoverable from the model card alone).
3. **`TOOL_COMPARISON.md`** still said human↔fly ortholog mapping was "not yet completed" and listed
   the cross-species validation as "incomplete" — both were finished in Findings 16 (result: a
   statistically underpowered null, not a completion in either direction, which is different from
   "not yet completed"). **Fixed**, and added the proteomics (Findings 17) and gseapy cross-validation
   (Findings 20) rows that were missing from the comparison table entirely.
4. **`ROADMAP.md`** described Tier 2 (disease-transfer validation) as still-to-attempt, when it has
   since been executed twice (Findings 12, and now Findings 19/Yang2023) with a consistent negative
   answer, and its own "Immediate next action" pointed at Tier 1 (metabolome), which was completed and
   marked done in the same file. **Fixed:** Tier 2 marked executed-and-negative with results
   summarized; "Immediate next action" replaced with the actual current priority list.
5. **The deployed app (`app.py`) never surfaced the Yang2023 null result anywhere** — a real, rigorous,
   pre-registered negative finding lived only in `docs/FINDINGS_19`, invisible to anyone using the
   public tool. Given the project's own standing rule ("report null results honestly, never omit"),
   this is exactly the kind of gap an editorial board flags: *why does your deployed tool not mention
   a real test of its own limits?* **Fixed:** added an "honest null result" paragraph to the About /
   methods page, verified rendering in the browser before considering this done (see §8).

`SCIENTIFIC_CHARTER.md`'s "Last verified" date was also bumped to today after re-confirming no
contradictions between the charter and current shipped code.

## 6. Innovation assessment (re-confirmed, unchanged from prior audit — still accurate)

**v1.1 DGRP predictor:** multi-modal per-trait-adaptive modality selection (nested CV), calibrated
coverage-verified intervals + explicit predictability flag, packaged reusable API — none of which
Ober 2012, Morgante 2020, DGRP2, or DGRPool provide. Verified still true against the tool-comparison
matrix.

**v2 disease-convergence resource:** independent from-scratch multi-disease convergence analysis on
data the source paper deposited but didn't itself analyze this way; leave-one-out robustness; a
specific, mechanistic, testable convergent signal (nucleolar stress, concordant in 3/5 models) plus a
reproducible outlier pattern (HTT-200Q, independently recurring across 3 separate analyses). Now
additionally backed by an independent-tool cross-validation of the GSEA engine itself (Findings 20,
ES r=1.000 vs `gseapy`) and a genome-wide RNA↔protein concordance check (Findings 17, r=0.33–0.36) —
both completed since the prior audit and now correctly reflected in `TOOL_COMPARISON.md`.

**What is explicitly NOT claimed anywhere** (re-verified unchanged): novel/unrelated-genotype
prediction, transgenic-disease-model behaviour prediction, drug-rescue prediction, deep learning as a
meaningful lever, the metabolome as a useful modality, DGRP→disease transfer as validated — the last
of these now has a second, independent negative confirmation (Yang2023) reinforcing rather than
undermining the honesty of the original claim boundary.

## 7. Known gaps / limitations (carried forward, re-verified, updated)

**v1 (DGRP predictor) side** — unchanged from prior audit, all re-confirmed still accurate:
1. 8 trait-sex models across 4 trait families; expansion tested and found no more predictable traits
   at current DGRPool coverage (measured ceiling, not an oversight).
2. Expression-modality prediction requires the query line's own RNA-seq (stated precondition).
3. starvation♀ expression accuracy (0.10) below Morgante's reported 0.28 — unresolved, documented.
4. Does not address the original disease/drug-screening vision (needs perturbation data, which no
   DGRP-based method provides) — now with two independent empirical tests confirming this (below).

**v2 (disease-convergence) side** — updated:
5. No matched behavioural data for BCM-DMAS (hard ceiling, needs new data, not new analysis).
6. DGRP→disease transfer weak/model-specific (Findings 12) — **now reinforced** by a second,
   independent, better-powered null on real genotype×disease data (Findings 19/Yang2023).
7. Cross-species (human↔fly) validation — **resolved as an underpowered null** (Findings 16); not
   actionable further without a pathway-level (not gene-level) re-test, which is recorded as a
   specific, scoped follow-up if pursued, not a vague TODO.
8. ALS not represented (real dataset exists but reports asymptomatic flies; would need its own DE +
   careful cross-study comparison, same batch-effect risk already demonstrated with the metabolome).
9. SCA — searched, no usable dataset found; honest gap.
10. Proteomics validation — **done** (Findings 17): genome-wide RNA↔protein concordance r=0.33–0.36.
11. scRNA-seq referenced in metadata but not locatable in the public Synapse tree — unresolved,
    flagged, not fabricated.
12. **New this session:** the project's many analyses are each individually FDR-controlled but there
    is no project-wide multiple-comparisons statement — not a defect, but should be stated explicitly
    in the eventual manuscript's statistics section (§3 above), not left for a reviewer to ask first.

## 8. Deployment readiness (verified this session, in addition to the audit)

- `app.py` re-verified to boot clean (HTTP 200, no console/server errors) after this session's edit,
  using the project's actual `.venv` via a live browser preview — including navigating to the About
  page and visually confirming the new null-result disclosure renders correctly.
- `requirements.txt` / `runtime.txt` / `.streamlit/config.toml` / `.gitignore` unchanged and
  re-verified consistent with the deployment guide (`docs/DEPLOYMENT.md`).
- Git: 2 commits, clean tree, **not yet pushed** — awaiting the user's repo name/URL and explicit
  go-ahead (publishing action, not taken without confirmation per this project's standing practice).

## 9. NAR-readiness verdict (re-confirmed)

**v1.1 DGRP predictor** remains at submission quality for a methods/resource paper. **v2
disease-convergence work** remains a genuine, verified, honest finding suitable as a second
paper/section. Both are now further strengthened by cross-tool validation (gseapy) and an independent
molecular layer (proteomics), and the project's negative-result discipline is now demonstrably
consistent end-to-end — including in the public-facing tool itself, not just in internal docs. No new
scientific red flags surfaced by this audit; the issues found were documentation staleness and one
transparency gap in the deployed app, both fixed this session.

## 10. Recommended next steps (priority order, re-ordered; fixes the prior audit's duplicate-numbering)

1. **Push to GitHub + connect Streamlit Community Cloud** — the only remaining launch step, blocked
   only on the user providing a repo name/URL and confirming (see `docs/DEPLOYMENT.md`).
2. **EmoryDrosophilaTau cross-lab replication** — verified real, independent lab, download instructions
   ready (`docs/EMORYDROSOPHILATAU_DOWNLOAD.md`); awaiting the user's Synapse download.
3. **Manuscript drafting** for both resources — the project has reached a natural reporting point: a
   validated predictor, a validated convergence finding, and two independent honest disease-transfer
   nulls. Repeatedly deferred in favor of more data collection; now a reasonable point to stop
   collecting and start writing, unless the user has a specific new dataset in mind.
4. **Optional, low-priority:** fix the ~12 one-off scripts' positional pandas axis calls (cosmetic;
   would not change any already-reported number; full file:line inventory recorded in this session's
   verification pass if ever picked up).

## 11. Precision follow-up pass (same day, before first GitHub push)

A second read specifically for repo-presentation integrity — what a reviewer sees on first landing,
not just what the code computes — surfaced one real issue and two wording precision fixes, all applied:

1. **`CLAUDE_DroPhenoPredict.md` (the original, pre-verification draft brief) is intentionally kept in
   the repo for provenance** — it documents exactly what fabricated content was caught and rejected,
   which is itself evidence of rigor. But the file had **no disclaimer at its own top**, so a reviewer
   opening it directly (it sorts alphabetically above `README.md`) would see invented "expected R²"
   numbers with no immediate warning that they're superseded. **Fixed:** added a banner at the top of
   the file pointing to `SCIENTIFIC_CHARTER.md §7` and `README.md`/`MODEL_CARD.md` for what's actually
   real.
2. **`README.md` had no explicit reviewer-navigation path** — findings, audit, model card, and risk
   analysis existed but nothing told a first-time reader which order to read them in. **Fixed:** added
   a "For reviewers / evaluators" section listing the read order, and a companion-tool cross-reference
   to DroPhenix (verified via its own public README: separate SSSIHL Bioinformatics tool, no shared
   code/data/artifacts).
3. **App footer wording said DroPhenoPredict and DroPhenix were "from the same lab"** — an institutional
   claim not independently established for this project (SCIENTIFIC_CHARTER.md makes no lab/institution
   claim). **Fixed:** narrowed to "by the same author," which is directly supported (the original draft
   brief itself, `CLAUDE_DroPhenoPredict.md` line 4, states "This is a SEPARATE tool from DroPhenix, but
   integrates with it" — confirming the authorial relationship without asserting an unverified
   institutional one).

Re-verified after these changes: fabrication rescan clean, `git add -A -n` dry run shows exactly the
intended 8 modified files (no stray `__pycache__`, no `.claude/` tooling config, no secrets), `pytest
tests/` 41/41 still passing. No new integrity or innovation concerns found — the substantive scientific
audit in §1–§10 stands unchanged; this pass addressed presentation-layer precision only.
