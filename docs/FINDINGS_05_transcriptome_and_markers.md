# Findings 05 — Transcriptome (Option 1) & marker-effect models (Option 2)

**Date:** 2026-06-29 · **Code:** `scripts/06_transcriptome_prediction.py`,
`07_multiomic_reml.py`, `08_marker_effects.py`, `09_supervised_expression.py`;
engine `gblup.py` (multi-kernel REML) · **Tests:** `tests/test_gblup.py` (7 pass)
**Data added:** DGRP per-line microarray expression (18,140 genes × ~184 lines × 2 sexes),
recovered from the Internet Archive (see `DATA_SOURCES.md` #3); thinned 100k-SNP matrix
(`data/processed/snp_matrix.npz`).

All compared against the genotype-only GBLUP baseline (Findings 03–04) on the SAME common lines
and SAME CV folds. The question both options answer: **can we beat genotype-only prediction?**

## Option 1 — does the transcriptome add predictive value? (three independent methods)

| Method | What it tests | Result |
|---|---|---|
| TRM-GBLUP (unsupervised) | whole-transcriptome similarity kernel | r ≈ 0 everywhere; never significant |
| G+T two-component REML (Morgante 2020) | optimal genomic+transcriptomic weighting | Δr (G+T − G) ≤ 0 in ~all cells; mean Δr < 0 |
| Supervised top-200-gene Ridge (within-fold) | most favourable: select trait-correlated genes | Δr < 0 in ~all cells (e.g. startle♂ 0.204→0.073) |

**Conclusion: baseline transcriptome adds no out-of-sample predictive value over genotype** for these
five DGRP traits, by any of the three standard methods. The transcriptome absorbs 0–21% of *in-sample*
variance under REML, but this does not generalize (it overfits). This nuances the original brief's
premise that multi-omics expression is THE predictive layer, and is method-fair w.r.t. Morgante et al.
2020 (whose reported gains were small and trait-specific).

*Important scope:* this concerns **baseline** (unperturbed, healthy-fly) expression predicting
**natural-variation** traits. It does NOT bear on disease-model or drug-perturbation contexts, where
expression changes are the signal — that is the (separate, currently data-limited) extension layer.

## Option 2 — do marker-effect (sparse QTL) models beat infinitesimal GBLUP?

Within-fold GWAS selection of the top-500 SNPs (training lines only → leakage-free) + RidgeCV:

| trait·sex | GBLUP r | marker r | Δr |
|---|---|---|---|
| starvation F | 0.236 | **0.312** | **+0.075** |
| starvation M | 0.221 | 0.202 | −0.019 |
| startle F | 0.113 | 0.095 | −0.019 |
| startle M | 0.097 | 0.064 | −0.034 |
| sleep F | 0.213 | 0.116 | −0.097 |
| (others) | ~0 | ~0 | mixed/noise |

**Conclusion: sparse marker models do not systematically beat GBLUP** (mean Δr < 0), confirming
Ober et al. 2012 ("neither implicit nor explicit marker selection substantially improves predictive
ability" in the DGRP). **One honest exception:** starvation♀ (marker r=0.31 vs GBLUP 0.24), hinting
that starvation resistance carries some larger-effect loci — consistent with known starvation QTL —
though starvation♂ did not replicate the gain. Reported as a lead, not a conclusion.

## Synthesis across Options 1–3

1. **Genotype-only GBLUP is the best general predictor in the DGRP** for these traits; ceiling
   **r ≈ 0.2** (starvation, startle, sleep♀ significant), reproducing Ober 2012.
2. **Neither transcriptome nor sparse markers raise the ceiling** in general. The limit is structural
   — an unrelated panel of n≈200 (Findings 03), not a missing feature layer or a wrong model.
3. **Uncertainty is calibrated throughout** (90% conformal coverage ≈ 0.90), so predictions are honest
   even when weak.
4. **Trait-specific leads worth pursuing:** starvation resistance (highest predictability + a
   marker-model gain) is the best candidate for a useful single-trait predictor.

## What this means for the tool (honest design consequence)

DroPhenoPredict's defensible product is **not** a high-accuracy point predictor. It is a
**calibrated, validated coarse prior with intervals** — most useful for ranking/prioritization on the
few predictable traits (starvation foremost), with explicit uncertainty. Raising accuracy needs
*more lines* (larger panels / related individuals to lift the relatedness ceiling), not more omics
layers on these 200 unrelated lines. The disease/drug-perturbation use case remains a separate,
data-limited extension (Charter §6), not addressed by these natural-variation baselines.
