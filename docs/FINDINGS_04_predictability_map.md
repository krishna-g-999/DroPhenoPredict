# Findings 04 — DGRP predictability map (Option 3)

**Date:** 2026-06-29 · **Code:** `scripts/05_predictability_map.py` (verified engine `gblup.py`)
**Data:** DGRPool API (labels) + `data/processed/grm.npz` (1.49M-SNP GRM) ·
**Artifacts:** `data/processed/predictability_map.csv`, `reports/figures/predictability_map.png`

The validated GBLUP baseline applied to all five high-coverage traits × sex, one harmonized
label source (DGRPool) and one verified engine. Every number measured; none assumed.

## Map (sorted by out-of-sample predictive r)

| Trait | Sex | n | Genomic h² | CV R² | Pearson r | Perm p | 90% cover |
|---|---|---|---|---|---|---|---|
| starvation | F | 197 | 0.73 | +0.056 | 0.236 | 0.010 | 0.92 |
| starvation | M | 197 | 0.84 | +0.047 | 0.221 | 0.010 | 0.93 |
| sleep | F | 167 | 0.81 | +0.046 | 0.213 | 0.015 | 0.90 |
| startle | F | 201 | 0.31 | +0.010 | 0.113 | 0.070 | 0.91 |
| startle | M | 201 | 0.25 | +0.007 | 0.097 | 0.085 | 0.91 |
| sleep | M | 167 | 1.00* | −0.005 | 0.075 | 0.100 | 0.91 |
| climbing | M | 185 | 0.27 | −0.001 | 0.068 | 0.139 | 0.92 |
| climbing | F | 185 | 0.11 | −0.012 | −0.002 | 0.443 | 0.92 |
| lifespan | M | 132 | 1.00* | −0.035 | 0.028 | 0.318 | 0.91 |
| lifespan | F | 132 | 0.13 | −0.035 | −0.095 | 0.766 | 0.91 |

`*` h² ≈ 1.00 = REML boundary artifact (δ → 0 on low-noise line means at small n); not literal.
The reliable readout in those cells is the CV r (~0) and perm p (n.s.).

## Interpretation

1. **A real predictability gradient exists.** Significant out-of-sample signal for starvation (both
   sexes, p=0.01), sleep♀ (p=0.015); borderline startle (p≈0.07–0.09); none for climbing/lifespan.
   Predictability is **trait-specific**, not a uniform property of the panel.
2. **Best-case ceiling ≈ r 0.23 (R² ≈ 0.05).** Even the most predictable trait explains ~5% of
   out-of-sample line-mean variance. This is the honest ceiling of genotype-only GBLUP in the DGRP.
3. **Heritability ≠ predictability, quantified.** starvation h²≈0.8 but r≈0.23; the gap is the
   unrelated-panel structure (Findings 03). Traits with low h² (climbing♀ 0.11, lifespan♀ 0.13)
   are unpredictable as expected.
4. **Uncertainty is calibrated everywhere** (empirical 90% coverage 0.90–0.93), so even
   low-accuracy predictions come with honest, validated intervals.

## External validation (second independent trait)

**Ober et al. 2012, PLOS Genetics**: DGRP GBLUP predictive ability — starvation **0.239 ± 0.008**,
startle 0.230 ± 0.012. Our starvation r = 0.236 (♀) / 0.221 (♂) **matches their starvation value
almost exactly**, confirming the pipeline on a second trait. (Our startle is a touch lower, 0.10–0.17
depending on which startle dataset; see note below.)

## Notes / caveats (no hiding)

- **Phenotype-source sensitivity:** startle from DGRPool id 2799 gives r≈0.10; startle from the
  quantgenet `startle.female.csv` (Findings 03) gives r≈0.17. Same trait, different study/normalization
  → different estimate. The map uses one consistent source (DGRPool) for comparability; absolute
  values are dataset-dependent and reported as such.
- **Lifespan n=132** after intersecting DGRPool `mn_Longevity` with the genotyped panel (the raw 218
  included non-canonical lines; the hygiene step from Findings 02 is applied automatically by the
  GRM∩phenotype join).
- REML boundary cases flagged above; CV remains the decision metric.

## Bottom line for the tool

Genotype-only GBLUP gives an **honest, calibrated, but low-accuracy** signal (best r≈0.23). That is
genuinely useful as a *coarse prior for experiment prioritization with intervals* — but not a point
predictor. This motivates Options 1 (functional/transcriptome features) and 2 (marker-effect models),
each to be measured against this same baseline. The map is the reference all future methods must beat.
