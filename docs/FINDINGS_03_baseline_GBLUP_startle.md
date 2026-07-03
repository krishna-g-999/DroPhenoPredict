# Findings 03 — Single-trait GBLUP baseline (startle response)

**Date:** 2026-06-29 · **Code:** `src/drophenopredict/{genotypes,phenotypes,gblup}.py`,
`scripts/03_build_grm.py`, `scripts/04_baseline_startle.py` · **Tests:** `tests/test_gblup.py` (6 pass)
**Artifacts:** `data/processed/grm.npz`, `data/processed/baseline_startle_{female,male}.json`,
`reports/figures/baseline_startle_{female,male}.png`

## Method (all from verified real data, no fabrication)

- **Target:** startle-induced locomotion, DGRP line means per sex. Source: quantgenet.msu.edu
  (`startle.female.csv` / `startle.male.csv`, Mackay/Huang DGRP resource).
- **Features → GRM:** DGRP Freeze 2.0 PLINK genotypes (Zenodo 5582846; NCSU host offline).
  Streaming VanRaden GRM over **1,493,351 autosomal SNPs** passing QC (MAF≥0.05, call-rate≥0.8),
  from 4,438,427 total variants × 205 lines.
- **Model:** GBLUP mixed model `y = 1μ + g + e`, `g ~ N(0, σ²_g K)`. REML for heritability;
  BLUP for prediction (δ = σ²_e/σ²_g refit per training fold).
- **Validation:** 5-fold × 10-repeat CV (line is the unit → no leakage); permutation test
  (200 perms); split-conformal 90% prediction intervals.

## Results (measured)

| Trait | n lines | Genomic h² (REML) | CV R² | Pearson r | MAE (null) | Perm p | 90% PI coverage |
|---|---|---|---|---|---|---|---|
| startle ♀ | 204 | 0.33 | +0.028 ± 0.013 | 0.169 | 5.11 (5.16) | 0.020 | 0.90 |
| startle ♂ | 204 | 0.27 | +0.025 ± 0.011 | 0.161 | 5.14 (5.20) | 0.020 | 0.91 |

Permutation null R² mean ≈ −0.029 (q95 = 0.018); observed exceeds it → signal is real.
r² (0.169²=0.029) ≈ R² → predictions unbiased, low-correlation. Conformal coverage matches target.

## Interpretation — the central, non-obvious finding

1. **Real but weak genomic predictability.** Startle is genomically predictable above chance
   (p=0.02) but only weakly: **r ≈ 0.16, R² ≈ 0.03**.
2. **This is NOT a bug — it is the DGRP's structure.** The GRM off-diagonal mean ≈ −0.005,
   SD ≈ 0.031: the DGRP is an **unrelated association panel**. GBLUP predicts held-out lines via
   genetic *relatedness*; with almost none, in-sample fit is high (R²≈0.65 on a simulated h²=0.8
   trait) but **out-of-sample CV R² ≈ 0**. On a synthetic *related* kernel the same engine recovers
   h²=0.72 and CV R²=0.78 — proving the engine is correct and the ceiling is the data
   (Daetwyler/Goddard: accuracy² ≈ h²·n/(n+Mₑ); Mₑ is huge in high-recombination *Drosophila*).
3. **Heritability ≠ predictability.** GREML h² ≈ 0.3 (moderate additive variance) coexists with
   R² ≈ 0.03 (weak out-of-sample prediction). Conflating the two is the error behind inflated
   claims (e.g. the multi-omics startle paper's in-sample "R²=0.95", which is model fit on 156
   lines, *not* cross-validated prediction).

## External validation (literature benchmark)

**Ober et al. 2012, PLOS Genetics 8(5):e1002685** — same method (GBLUP+GRM), same trait (startle),
same panel (DGRP) — reported cross-validated predictive ability **r = 0.230 ± 0.012** (starvation
0.239 ± 0.008). Our r = 0.16 is the same order of magnitude; the gap is consistent with our
smaller per-fold training fraction (5-fold vs their scheme). Ober's conclusions ("marker selection
does not help; little gain beyond 150k SNPs") match our unrelated-panel explanation. **Our pipeline
reproduces the canonical DGRP genomic-prediction result.**

## Consequences for DroPhenoPredict (honest)

- The brief's premise — *enter any genotype → predict phenotype* — is **structurally limited** in the
  DGRP via relatedness-based genomic prediction. Within-panel accuracy ceilings at **r ≈ 0.2**; for a
  *novel, unrelated* genotype (the actual use case) relatedness-based transfer is ≈ 0.
- Replaces the fabricated "R²=0.55–0.75 / 74% accuracy" with a **measured, literature-consistent
  ceiling**.
- Paths that could raise predictability (each its own future Finding, none assumed to work):
  (a) functional/pathway features (DGRP transcriptome — currently access-blocked; needs sourcing);
  (b) marker-effect models that capture large-effect QTL (BayesB/LASSO) — Ober found little gain;
  (c) larger panels / multi-trait borrowing; (d) reframing output as a calibrated coarse prior with
  honest intervals rather than a point predictor.

## Reproduce
```
./.venv/Scripts/python.exe scripts/03_build_grm.py        # builds data/processed/grm.npz (once)
./.venv/Scripts/python.exe scripts/04_baseline_startle.py # metrics + figures
./.venv/Scripts/python.exe -m pytest tests/test_gblup.py  # 6 engine/structure tests
```
