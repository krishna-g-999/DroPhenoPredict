# Findings 16 — Human↔fly ortholog mapping + cross-species overlap test

**Date:** 2026-07-02 · **Code:** `scripts/37_panther_ortholog_mapping.py`
**Data:** Ruffini et al. 2020 139-gene human shared-NDD list (obtained Findings 15/AUDIT gap #7, via
EuropePMC) + BCM-DMAS 27-gene fly convergent set (Findings 13).

## What was completed this session
The ortholog-mapping barrier flagged in `AUDIT.md` §7 (5 sources tried, all blocked/dead) is now
resolved: **PANTHER's ortholog REST API** (pantherdb.org, a peer-reviewed, GO-Consortium-affiliated
resource) is reachable and correct. Verified on spot checks before bulk use: MAPT→**tau**, APP→**Appl**,
SNCA→**no fly ortholog** (matches our existing curated table in `disease.py`, which encodes SNCA as
human-transgene-only — independent cross-check, not just a new data point).

**Engineering note:** PANTHER's batch endpoint silently errors above ~10–14 genes per request
("more items than supported"); empirically confirmed the limit sits between 10 and 15. Fixed to
chunks of 10, with per-chunk checkpointing (a transient server error mid-run does not lose completed
work) and hard failure (not silent empty-result) on unrecoverable errors.

## Result
- **139 human genes queried → 174 fly-ortholog relationships resolved** (many human genes have
  multiple fly paralogs, e.g. VAMP2→Syb *and* nSyb) — **173 unique fly genes**.
- **Cross-species overlap with our BCM-DMAS convergent set (27 genes, DE in ≥3/5 disease models):**

| | value |
|---|---|
| Background (expressed genes, BCM-DMAS) | 9,385 |
| Ruffini-ortholog genes in background | 136 |
| Our convergent genes in background | 27 |
| **Observed overlap** | **0** |
| Expected overlap by chance | 0.391 |
| Hypergeometric p-value | 1.0 |

## Honest interpretation — this is a NULL result, not a negative one
**This test has essentially no statistical power.** With an expected overlap of 0.39 genes, observing
0 is the single most likely outcome under the null hypothesis *regardless of whether a real biological
relationship exists* — a hypergeometric test cannot distinguish "no true signal" from "a real but modest
signal invisible at gene-set sizes this small." **We do not claim this disproves cross-species
convergence.** We report: *at current gene-set sizes, no detectable overlap; the test is underpowered
to conclude anything stronger.*

Plausible, non-mutually-exclusive reasons for zero overlap, none of which we can adjudicate with this
test alone:
1. **Different biological contexts.** Ruffini's 139 genes converge across *human* AD/PD/HD/ALS from
   177 heterogeneous clinical studies (multiple tissues, platforms, disease stages). Our 27 genes
   converge across *engineered Drosophila overexpression models* in one tissue (head), one lab, one
   platform. Natural human disease progression and acute transgene overexpression are different
   perturbation classes.
2. **Fly-specific gene composition.** Several of our convergent genes (immune-response *IM23*,
   antimicrobial-peptide-like *Drsl4*, detoxification *Cyp6a14/Cyp6a17*) belong to fly-specific gene
   families with no clean 1:1 mammalian ortholog — inherently unable to appear in a human-gene-derived
   list regardless of true biological relevance.
3. **Small-N inherent limit.** 27 genes is a strict, high-confidence convergent set (DE in ≥3/5 models,
   FDR<0.05, |log2FC|>1); a larger, looser convergent set (e.g. the 60 GSEA-convergent *pathways*,
   Findings 15) would be the fairer comparison unit but pathway-level cross-species comparison needs a
   different (GO-term-based, not gene-based) test — not yet run.

## What this does NOT undermine
The **within-fly** convergence finding (Findings 13, 15 — 27 genes and 60 pathways converging across
the 5 Drosophila models, both significant by permutation, robust to leave-one-out) is unaffected by
this result. That finding rests on its own internal statistics, not on agreement with the human list.

## Honest status update to AUDIT.md §7 item 7
Ortholog mapping is **no longer blocked** — PANTHER is a working, verified, scriptable source
(`src/drophenopredict` could gain a `panther.py` module if this becomes a recurring need). The
cross-species *validation* itself is inconclusive (underpowered null), not completed successfully in
either direction. Recommended follow-up if pursued further: a pathway-level (GO-term) cross-species
comparison instead of gene-level, which would have far more statistical power than comparing two small,
strict gene lists directly.
