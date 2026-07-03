# Model Card — DroPhenoPredict v1.1 (multi-modal)

## Overview
DroPhenoPredict predicts organismal phenotypes of **Drosophila Genetic Reference Panel (DGRP)**
lines and returns **calibrated 90% prediction intervals**, the **data modality** used, and an honest
**per-trait predictability flag**. v1.1 is **multi-modal**: per trait it adaptively uses the best of
genotype (genomic GBLUP) or RNA-seq gene expression, chosen by **nested cross-validation**.

- **Modalities:** `geno` = VanRaden GRM (1,493,351 autosomal SNPs); `exp` = RNA-seq expression
  relationship (GSE117850, 15.7–20.4k genes). Selection by nested CV (inner CV picks modality on
  training data only → unbiased).
- **Targets:** 8 trait-sex models — climbing (negative geotaxis), startle locomotion, starvation
  resistance, night sleep.
- **Intervals:** k-fold split-conformal; empirical coverage verified at 0.90–0.92 for the 90% target.
- **Artifacts:** `models/trained/drophenopredict_v1.1.json`, `..._kernels.npz`; API
  `drophenopredict.predictor.DroPhenoPredictor`.

## Performance (nested-CV adaptive r — unbiased headline metric)
| Model | n | deployed modality | geno r | exp r | **adaptive r** | perm p | 90% coverage |
|---|---|---|---|---|---|---|---|
| starvation ♂ | 192 | exp | 0.268 | 0.376 | **0.343** | 0.005 | 0.92 |
| climbing ♂ | 180 | exp | 0.103 | 0.293 | **0.273** | 0.005 | 0.91 |
| sleep ♀ | 162 | geno | 0.252 | neg | **0.249** | 0.015 | 0.90 |
| starvation ♀ | 192 | geno | 0.227 | 0.099 | **0.159** | 0.025 | 0.90 |
| startle ♀ | 197 | geno | 0.159 | neg | **0.149** | 0.025 | 0.92 |
| climbing ♀ | 180 | exp | 0.009 | 0.166 | **0.134** | 0.050 | 0.91 |
| startle ♂ | 197 | geno | 0.125 | neg | 0.102 | 0.119 | 0.92 |
| sleep ♂ | 162 | geno | 0.105 | neg | 0.093 | 0.144 | 0.91 |

6/8 models significantly predictable (perm p<0.05). The adaptive r is reported, NOT the optimistic
leave-one-out reconstruction. Multi-modal gain is largest where genotype fails: **climbing♂
0.10→0.27, climbing♀ 0.01→0.13** — driven by RNA-seq expression.

## Validation / provenance (what makes it trustworthy)
- **Unbiased accuracy:** nested CV — modality selection never sees the test fold.
- **Calibrated uncertainty:** split-conformal with measured coverage (≈ nominal).
- **Significance:** 200-permutation test per model vs the adaptive R².
- **Engine verified** on known-truth simulations (`tests/`, 41 tests pass, verified 2026-07-03).
- **Reproduces 3 literature benchmarks:** Ober 2012 (genotype startle/starvation r≈0.23/0.24),
  Morgante 2020 (expression starvation♂ r≈0.38). All inputs public (`docs/DATA_SOURCES.md`).

## Intended use
- **Primary:** prioritize *which DGRP lines to phenotype* — rank lines by predicted value with honest
  intervals before committing fly work. One RNA-seq profile (already public for the panel) informs
  multiple behavioral predictions.
- **NOT for:** novel/unrelated genotypes, transgenic disease models, or drug responses (Limitations).
  The API refuses out-of-panel lines instead of fabricating a value.

## Limitations (honest)
1. **Modest accuracy (adaptive r ≈ 0.09–0.34).** A *calibrated prior*, not a precise point predictor —
   always read the interval.
2. **In-panel only.** Lines unrelated to the panel are ≈ trait mean; the API does not predict them.
   Expression-modality predictions also require the line's RNA-seq profile.
3. **Natural variation only.** Says nothing about disease-transgene or drug-perturbation phenotypes
   (the original aim; needs perturbation data — a CMap-for-flies resource — not implemented).
4. **Deep learning does not help** at this scale (tested; Findings 06); the gain is the RNA-seq
   *modality*, integrated linearly.

## Ethics / reproducibility
*Drosophila* only; designed to reduce animal/resource use. Seeded (42). Rebuild:
`scripts/03_build_grm.py` → `scripts/13_finalize_multimodal.py`; validate `pytest tests/` (41 tests).

## Related resource (separate model, same tool)
A second, separately-scoped resource ships in the same app: a multi-disease (AD/PD/HD)
**molecular-convergence** analysis on the BCM-DMAS fly panel — not a phenotype predictor, no
shared model artifact with the above. See "Disease convergence (v2)" in `app.py` and
`docs/TOOL_COMPARISON.md` Part 2. A related, honest **null result** — this predictor's own
GBLUP/multi-modal engine applied to a real DGRP×AD-transgene degeneration phenotype (Yang et al.
2023) — found 0/5, 0/10, 0/5 significant across three robustness checks (`docs/FINDINGS_19`).
Reported for completeness: this predictor's accuracy claims above are for the 8 trait-sex models
it actually ships, not a claim that the same approach works for every DGRP-adjacent phenotype.
