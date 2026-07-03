# Findings 12 — DGRP→disease transfer test on real multi-disease fly data

**Date:** 2026-07-01 · **Code:** `scripts/31_build_bcm_dmas_matrix.py`, `scripts/32_disease_transfer_test.py`
**Data (real, matched):** BCM-DMAS, Synapse syn34767207 (Farrell/Shulman, npj PD 2025) — head RNA-seq,
**5 disease models + 2 controls** on one Elav-Gal4 background, female, timecourse days 2–57:
AD (Aβ42, Tau2N4R), PD (αSyn), HD (HTT-128Q, HTT-200Q). 17,478 genes × 147 samples.
DGRP side: per-gene Pearson r(female expression, female climbing) across ~180 lines.

## Question
Does the DGRP natural-variation gene→climbing map transfer to disease models — i.e., are
"good-climbing" genes systematically down-regulated in disease (all models lose climbing)?

## Result
Spearman ρ between DGRP climbing-weight and disease-vs-control log2FC (age-matched), 11,078 shared genes:

| model | ρ | perm p |
|---|---|---|
| PD αSyn | +0.147 | 0.002 |
| AD Aβ42 | +0.082 | 0.002 |
| AD Tau | +0.065 | 0.002 |
| HD HTT-128Q | +0.027 | 0.004 |
| HD HTT-200Q | −0.101 | 0.002 |

## Interpretation (honest)
1. **A real but weak, model-specific relationship exists.** All 5 are significant (perm p<0.005) but
   effect sizes are small (|ρ| 0.03–0.15); significance is driven partly by n≈11k.
2. **The naive directional hypothesis is NOT supported.** 4/5 models show the OPPOSITE sign
   (good-climbing genes tend to be *up*, not down); only the most severe model (HTT-200Q) matches
   the predicted negative transfer. So "disease = loss of the natural good-climbing program" is too
   simple.
3. **Caveats that must be resolved before any claim:** DGRP expression is whole-body while BCM-DMAS is
   head (tissue mismatch); no control for baseline expression level / gene length confounds; ρ this
   small is not a usable predictor. This is a first-pass directional probe, not a validated transfer.

## What this means for the vision
- The pipeline **works end-to-end on real multi-disease data** (a genuine milestone — DGRP natural
  variation connected to 5 disease signatures with real code, no fabrication).
- But a **direct, simple DGRP→disease behavioural predictor is not supported** by this test — consistent
  with everything we've measured (unrelated-panel ceiling; cross-context transfer is hard).
- The **more robust use of this dataset** is the shared-biology analysis it was built for: do the 5
  disease models share DE genes/pathways (testing Ruffini et al. 2020's cross-disease convergence directly in
  flies, on matched data)? That does NOT depend on the DGRP bridge and is the recommended next step.

## Status of the v2 "universal disease model"
The shared-biology premise is real (Ruffini et al. 2020) and the matched multi-disease omics is real
(BCM-DMAS). What is NOT real is a supervised behaviour predictor: the behavioural values are locked in
figures (not data), and the DGRP→disease transfer is weak/complex. So the honest v2 contribution is a
**multi-disease molecular-convergence + pathway resource**, not a behaviour oracle. No fabricated
KNOWN_PHENOTYPES / ODE params / effect sizes were used (all excluded per Charter §7).
