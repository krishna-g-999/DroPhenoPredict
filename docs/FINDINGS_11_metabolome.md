# Findings 11 — Metabolome modality: obtained, heritable, but NOT shipped (honest negative)

**Date:** 2026-07-01 · **Code:** `scripts/22_build_metabolome.py`, `scripts/23_metabolome_modality.py`
**Data:** MetaboLights **MTBLS2060** — 673 Bruker ¹H-NMR spectra (800 MHz), 170 DGRP lines, male pools.

## What we did (rigorously)
- Downloaded all 673 sample zips; extracted each authors'-processed spectrum (`pdata/1/1r` + `procs`)
  — using their FT/phasing/baseline, not a re-derivation.
- Binned (0.02 ppm, 0.5–10 ppm, water 4.5–5.1 excluded), **PQN-normalized + log-transformed**,
  averaged replicates → **170 lines × 475 bins**. Line ids harmonized (`DGRP-021`→`line_21`).
- Built a metabolomic relationship matrix (MRM); tested genotype | expression | **metabolome** |
  combined for male traits + 9 locomotor-activity phenotypes; validated vs Heredity 2021.

## Verification results
1. **Heritable / reproducible:** within-line replicate correlation **+0.20** vs between-line −0.02 —
   the metabolome distinguishes lines reproducibly (consistent with the paper).
2. **Strong line structure:** MRM off-diagonal SD ≈ **0.42** (far more structure than the GRM's 0.03).
3. **But orthogonal to the phenotypes:** corr(MRM, GRM) ≈ 0.01; metabolome REML h² for the traits ≈ 0;
   **metabolome CV r is NEGATIVE for every trait** (−0.13 to −0.24), across three MRM variants
   (standardized / intensity-weighted / top-variance bins).
4. **Failed the paper's validation:** locomotor activity metab r ≈ −0.16 vs the reported >0.4.

## Interpretation (honest)
The metabolome has strong, heritable line structure that is **orthogonal to the cross-study DGRPool
behavioural phenotypes** — so GBLUP makes confident, wrong predictions (hence negative r). The likely
reasons the paper's PA>0.4 did not reproduce:
- **Matched-sample requirement:** the paper predicted phenotypes measured on the SAME flies/experiment
  as the metabolome; DGRPool phenotypes are from DIFFERENT studies/conditions.
- **Condition lability:** the metabolome is far more environment/state-sensitive than the
  transcriptome. RNA-seq expression (also single-study) DID transfer to cross-study traits
  (starvation r=0.38); the metabolome did not.
- Possibly the authors' exact pipeline (referencing, spectral alignment, metabolite quantification)
  extracts more transferable signal than standard binning.

## Decision
**The metabolome is NOT added as a predictive modality** — we could not demonstrate it improves the
validated benchmark, and shipping an unhelpful modality would violate the project's integrity rule.
The data, processing code, and this negative result are retained.

## Scientific value of this negative
It sharpens the project's central lesson: **molecular-state predictors are condition/study-specific;
cross-study transfer is not guaranteed and must be tested.** The metabolome — heritable yet
non-transferable here — is a concrete cautionary case, and directly informs the disease-transfer
design (Roadmap Tier 2): external prediction must be validated on matched samples, never assumed.
