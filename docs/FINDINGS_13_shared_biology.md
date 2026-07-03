# Findings 13 — Cross-disease molecular convergence in the fly (v2 premise, tested)

**Date:** 2026-07-01 · **Code:** `scripts/33_shared_biology.py` (DE math verified vs statsmodels)
**Data:** BCM-DMAS matched panel (syn34767207) — head RNA-seq, 5 disease models + driver control,
9,385 expressed genes. Control = Elav-Gal4/+ (isolates the disease transgene from the Gal4 driver).

## Question
Does the shared-neurodegeneration premise (Ruffini et al. 2020 (Cells 9:2642)) hold in flies, on matched data?
I.e., do AD (Aβ42, Tau), PD (αSyn), HD (HTT-128Q, HTT-200Q) converge on shared DE genes/pathways?

## Method (each step verified)
- **DE per model** vs driver control, age-adjusted linear model on log2-CPM.
  **Verified**: vectorized OLS reproduces statsmodels OLS exactly (log2FC & t identical) before use.
- **Convergence**: genes DE (FDR<0.05, |log2FC|>1) in ≥k models; significance by column-permutation
  null (independently shuffle each model's DE labels).
- **Direction**: genome-wide Spearman correlation of log2FC across models.
- **Pathways**: hypergeometric GO-BP enrichment of convergent genes vs expressed background.

## Results (measured)
- DE genes per model: Aβ42 41, Tau 33, αSyn 70, HTT-128Q 124, HTT-200Q 220.
- **Convergence is significant:** 27 genes DE in ≥3 models vs **0.1 expected** under independence
  (perm p = 0.0005); 4 genes DE in *all 5*.
- **Directional concordance for 4/5 models:** genome-wide log2FC Spearman 0.29–0.54 among Aβ42, Tau,
  αSyn, HTT-128Q. **HTT-200Q is a distinct outlier** (correlations ≈ 0 / slightly negative) — the same
  model that was the outlier in the DGRP transfer test (Findings 12).
- **Convergent genes** include genuine shared-NDD markers: immune/neuroinflammation **IM23, Drsl4**;
  oxidative/detox **Cyp6a14/Cyp6a17**; cytoskeleton/axonal transport **Map205, Dhc98D** — plus some
  cuticular/odorant genes (Cpr, Obp, Amy-p) likely reflecting head-tissue variation.
- **Pathway enrichment:** only *defense response* (immune) emerges, and only suggestively
  (p=0.068, FDR=0.068, not significant) — underpowered at 27 genes.

## Interpretation (honest)
1. **The shared-biology premise IS supported in flies on matched data** — significant, directionally
   concordant molecular convergence across 4 of 5 neurodegeneration models. This is the real,
   defensible core of the "one shared model" idea, at the *molecular* level.
2. **It is not pristine:** HTT-200Q (aggressive full-length polyQ) has a distinct signature; the
   convergent gene set is a mix of bona-fide NDD mechanisms (immune, oxidative, transport) and some
   likely tissue/technical genes; single-pathway enrichment is underpowered.
3. **This is molecular convergence, NOT a behaviour predictor.** It supports a shared *pathway/feature*
   layer; it does not, by itself, predict climbing/lifespan (behaviour values remain locked in figures).

## What this means for v2
The honest, real v2 contribution is a **matched multi-disease molecular-convergence resource** in
Drosophila that empirically supports shared-neurodegeneration biology and identifies a convergent gene
set — reproducing Ruffini et al. 2020's cross-disease premise in a controlled fly panel. All computed from real
data; the DE engine is verified; nothing fabricated (no KNOWN_PHENOTYPES / ODE params / effect sizes).
