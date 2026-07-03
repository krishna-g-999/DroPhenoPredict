# Tool & Method Landscape — where DroPhenoPredict sits

Verified comparison of DroPhenoPredict against existing DGRP resources and the two
genomic-prediction papers it builds on. Claims about other tools reflect their
documented capabilities (checked 2026-06-29).

## Capability matrix

| Capability | DGRP2 web (NCSU) | DGRPool (EPFL) | Ober 2012 | Morgante 2020 | **DroPhenoPredict** |
|---|---|---|---|---|---|
| Type | GWAS web tool | phenotype repository + GWAS/PheWAS | research paper | research paper | **packaged predictor + API** |
| Harmonized DGRP phenotypes | submit-your-own | **903 phenotypes, 142 studies** | — | — | consumes DGRPool |
| Single-marker GWAS | ✓ | ✓ | — | — | ✓ (architecture analysis) |
| Phenotype **prediction** for a line | ✗ | ✗ | ✓ (genotype) | ✓ (geno+expr) | ✓ |
| Genotype modality | ✓ (assoc only) | ✓ (assoc only) | ✓ GBLUP | ✓ | ✓ GBLUP |
| Transcriptome modality | ✗ | ✗ | ✗ | ✓ (RKHS) | ✓ (RNA-seq) |
| **Per-trait adaptive modality** | ✗ | ✗ | ✗ | ✗ | **✓ (nested CV)** |
| **Calibrated prediction intervals** | ✗ | ✗ | ✗ | ✗ | **✓ (conformal, coverage-verified)** |
| **Per-trait predictability flag** | ✗ | ✗ | ✗ | ✗ | **✓ (permutation test)** |
| Programmatic API / reusable model | partial | API (data) | ✗ | ✗ | **✓ (`DroPhenoPredictor`)** |
| Status | offline (2026) | live | — | — | local, reproducible |

## What is genuinely new here
1. **Multi-modal, per-trait-adaptive integration** of genotype and transcriptome with the modality
   chosen by nested CV — neither Ober (genotype) nor Morgante (fixed combination) does adaptive
   per-trait selection.
2. **Calibrated, coverage-verified uncertainty + an explicit predictability flag** — no prior DGRP
   prediction work ships honest intervals or tells the user when a trait is *not* predictable.
3. **A packaged, reusable resource** (model artifact + API), not a one-off analysis.

## What is NOT claimed (editor's guard)
- We do not beat the information-theoretic ceiling; accuracy is modest (adaptive r ≈ 0.09–0.34) and
  reported as such.
- We do not predict novel/unrelated genotypes, transgenic disease models, or drug responses — those
  need perturbation data, which no DGRP-based method (ours included) provides.
- Deep learning was tested and does **not** help at this scale (Findings 06); we do not claim it.

## Reproduction of prior work (validation, not novelty)
- Ober 2012 (genotype): starvation r≈0.24 — we measure 0.23/0.27. Reproduced.
- Morgante 2020 (expression): starvation♂ r=0.38 — we measure 0.376. Reproduced.
- Honest discrepancies: startle genotype (we 0.13–0.16 vs 0.23) and starvation♀ expression
  (we 0.10 vs 0.28) are lower; likely CV-scheme/normalization differences. Reported, not hidden.

---

## Part 2 — Disease-convergence side (v2): where it sits

This is a **separate resource** from the DGRP predictor above (different data, different question:
molecular convergence across neurodegeneration models, not trait prediction for a fly line).

| Capability | Ruffini et al. 2020 (human, *Cells*) | Farrell/Shulman *npj PD* 2025 (fly, source data) | **DroPhenoPredict v2** |
|---|---|---|---|
| Organism | human (177 studies, 1M+ patients) | *Drosophila* (5 models, 1 lab) | *Drosophila* (same data as npj PD 2025) |
| Diseases | AD, PD, HD, ALS | PD only (αSyn) | **AD, PD, HD** (5 models) from the same matched panel |
| Convergence method | DEG intersection across studies | not the paper's focus (single-disease) | **DEG intersection (verified) + genome-wide GSEA (locally implemented, verified)** |
| Cross-study batch handling | meta-analysis across heterogeneous studies | single lab/platform (no batch issue) | single lab/platform (same as npj PD 2025) — avoids the batch-effect risk of pooling heterogeneous studies |
| Robustness testing | not reported | not applicable | **leave-one-out reported** (convergence survives excluding any model) |
| Independent DE validation | — | — | **cross-checked against the paper's own DESeq2 output** (r=0.77, 7,018 genes) |
| Human↔fly cross-species check | N/A (human-only) | N/A | 139-gene human list obtained (Ruffini S1, via EuropePMC); ortholog mapping not yet completed (network access barrier this session, documented in AUDIT.md) |

## What is genuinely new on the disease side
1. **A verified, from-scratch multi-disease convergence analysis in flies** — the npj PD 2025 paper
   deposited the data but analyzed only its own PD arm; we independently computed DE for all 5 models
   and tested convergence, which the original paper does not report.
2. **Leave-one-out robustness** — not done in either comparator paper for this question.
3. **A specific, mechanistic, testable convergent signal** (nucleolar stress/rRNA biogenesis
   down-regulation, concordant in 3/5 models) plus a **recurring outlier finding** (HTT-200Q diverges
   directionally in 3 independent analyses) — genuinely novel observations from this matched panel.

## What is NOT claimed (disease side)
- This is **not** a behaviour predictor — matched climbing/turning/stumbling data for this panel exists
  only in the source paper's Figure 1B, not as downloadable data.
- The DGRP→disease transfer test (Findings 12) found only a weak, model-specific relationship — we do
  **not** claim the DGRP predictor extends to these disease models.
- The human-gene cross-species validation is incomplete (data in hand, ortholog mapping pending).
