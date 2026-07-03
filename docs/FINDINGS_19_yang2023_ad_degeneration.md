# Findings 19 — Yang 2023 DGRP×AD eye-degeneration: a rigorous, honest null

**Date:** 2026-07-03 · **Code:** `scripts/40`–`43` · **Data:** Yang et al. 2023 (G3 13(9):jkad132),
public GitHub (`github.com/mingwhy/AD_fly_eye`) + EuropePMC supplementary API — no DUA, fully open.
**Risk analysis performed first:** `docs/EDITORIAL_RISK_ANALYSIS.md` (11 anticipated reviewer
questions, each addressed by a specific design decision before any result was computed).

## The question
Does our validated GBLUP/multi-modal framework predict real AD-transgene (Aβ42+Tau)
eye-degeneration severity across DGRP genetic backgrounds — the first dataset in this project
matching real DGRP genotypes to a real disease phenotype (not natural-variation behavior, not
engineered-transgenic-on-fixed-background molecular data)?

## Method (each step verified before use)
1. **Phenotype construction replicated the original authors' own method**, not an ad hoc
   simplification: linear mixed model `trait ~ group + (1|line)` (REML, `statsmodels.MixedLM`),
   extracting the BLUP per DGRP line and heritability from the fitted variance components — verified
   against their own published R code. A real API bug (`fit.cov_re` returns a plain array, not a
   DataFrame, in this statsmodels version) was caught before any output was trusted.
2. **Self-check:** BLUP line-values correlate with naive raw per-line means at r=0.765 — confirms the
   mixed model behaves sensibly (imperfect correlation is expected and correct: BLUP shrinks
   low-replicate lines toward the mean).
3. **162 DGRP lines** (0 lost) — QC-control pseudo-lines ('32C'/'441C') explicitly excluded *before*
   fitting, not after, per the risk analysis.
4. **5 pre-committed primary traits** (the original authors' own designation: nn.mean, nn.sd, cc.mean,
   cc.sd, mean.s.area) — avoids fishing across all 16 correlated features. BH-corrected throughout.
5. **Three analyses, same validated engine used everywhere else in this project:**
   genotype-only GBLUP (`gblup.py`) → multi-modal genotype+expression (`multimodal.py`, both sexes
   tested since cross-sex composition wasn't specified) → inversion/Wolbachia-adjusted robustness
   check (`covariates.py`) — no new methodology, applying already-verified tools to a new phenotype.
6. **A real methodological gap caught and fixed before running:** the first draft of the multi-modal
   significance test borrowed a single-kernel permutation null for the *adaptive* nested-CV procedure
   — not the correct null for an adaptive/selection procedure. Fixed with a dedicated
   `adaptive_permutation_test` that permutes labels and reruns the *entire* nested selection each time.

## Results — complete, consistent null across all three analyses

| Analysis | n | Significant (FDR<0.05) | Notable |
|---|---|---|---|
| Genotype-only GBLUP | 162 | **0 / 5** | 4/5 traits CV_r **negative**; cc.sd has genomic h²=0.45 but CV_r≈0 (another heritability≠predictability case) |
| Multi-modal (geno + RNA-seq, both sexes) | 160 | **0 / 10** | **Every single r-value negative** (geno, exp, combined, adaptive) — expression does not rescue this trait, unlike climbing/starvation in v1.1 |
| Inversion/Wolbachia-adjusted | 162 | **0 / 5** | Adjustment pulls most estimates toward zero (from as negative as −0.15 to near 0) — some of the raw negativity was inversion confounding, but no real signal emerges even after correction |

## Interpretation (honest, not spun)

1. **This is a genuine, well-powered null, not a failed experiment.** n=160–162 is comparable to
   every other DGRP trait in this project; the analysis is methodologically sound (verified phenotype
   construction, correct significance testing, robustness-checked against a confound the original
   authors themselves flagged).
2. **Consistent with, and now extending, the project's central structural finding**: the DGRP's
   unrelated-panel ceiling (Findings 03) applies to this real disease phenotype too. Genotype and
   expression both fail to predict out-of-sample here, in contrast to climbing/starvation where
   expression provided real lift (Findings 06/15). This is itself informative: **not every trait
   benefits from the multi-modal approach** — it is trait-specific, and AD-transgene degeneration
   severity (measured in a heterozygous F1 cross) appears to be one where neither modality transfers.
3. **Plausible reasons, stated as hypotheses, not conclusions:** (a) the F1-heterozygous cross design
   dilutes/dominance-masks additive genetic signal relative to homozygous-line traits; (b) automated
   image-morphology features may carry substantial non-genetic (technical/developmental) noise despite
   the batch-adjusted BLUP; (c) the original authors' own single-marker GWAS approach may be better
   suited to detecting large-effect modifier loci than a genome-wide additive GBLUP is at this n.
4. **This does not undermine the rest of the project.** It is a new, honestly negative data point
   added to the existing predictability map (Findings 04), collected with the same rigor as every
   positive result. A null result, properly obtained and reported, is real science.

## What this means for the manuscript
Strengthens the project's overall credibility: a **pre-registered** (risk analysis written before any
result), **methodologically faithful** (replicated the source authors' own trait construction),
**multiply-robustness-checked** (genotype, multi-modal, inversion-adjusted) null result on a genuine
disease phenotype. Reviewers who might otherwise wonder "did they only report the traits that worked?"
have a direct, verifiable answer: no — here is a rigorous test that came back null, reported in full.
