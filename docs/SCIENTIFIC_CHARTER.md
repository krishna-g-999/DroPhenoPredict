# DroPhenoPredict — Scientific Charter

**Status:** Foundational governance document. Every modeling decision must be traceable to this file.
**Last verified:** 2026-06-28
**Governing principle:** No fabricated data. No invented numbers. No unverified claim presented as fact. Every value in this project is either (a) downloaded from a cited public source, (b) computed by reproducible code from such a source, or (c) explicitly labelled as a hypothesis/placeholder pending verification.

---

## 1. What this project is (honest one-paragraph definition)

DroPhenoPredict is a **genomic / multi-omics phenotype-prediction model for the Drosophila Genetic Reference Panel (DGRP)**. Given the genotype (and, where data exists, transcriptome-derived features) of a DGRP line, it predicts physiological/behavioral trait values (e.g. startle locomotion, sleep, lifespan, starvation resistance) **that have measured ground truth in the panel**, with rigorously cross-validated, honestly-reported accuracy and calibrated uncertainty.

Everything beyond that — predicting transgenic human-disease models, drug-rescue scoring — is a **separate, explicitly-labelled extrapolation layer** (see §6), not part of the validated core, and is framed as hypothesis-generation, never as validated prediction.

## 2. Why the DGRP makes this buildable (verified)

| Fact | Value | Source |
|---|---|---|
| Sequenced inbred lines | 205 | Mackay et al. 2012, *Nature* 482:173 |
| Variants (Freeze 2.0) | ~4.85M SNPs + ~1.3M non-SNP | Huang et al. 2014, *Genome Res* (PMC4079974) |
| Per-line transcriptome | males + females, genome-wide | Huang et al. 2015, *PNAS* 112:E6010 |
| Harmonized phenotypes | 139 studies / 838 phenotypes / 22 categories | DGRPool, Gardeux et al., *eLife* 88981 |
| Lines still live at Bloomington | 192 | DGRPool / BDSC |

These are the only data assets we treat as given. See `DATA_SOURCES.md` for access methods and verification status.

## 3. The central scientific constraint (the thing the original brief ignored)

**The DGRP is natural standing genetic variation in non-transgenic, wild-derived flies.** It contains *no* human-disease transgenes (TDP-43, FUS, SOD1, HTT). Therefore:

- A model trained on the DGRP learns **genotype → trait mapping for natural variants**.
- It has **zero training examples** of human-transgene proteotoxic overexpression.
- Predicting a transgenic disease-model phenotype from a DGRP-trained model is **extrapolation outside the training distribution** — scientifically it is hypothesis-generation, not prediction. It must be validated against held-out *transgenic* data before any predictive claim is made, and even then only within the validated domain.

This constraint defines the boundary between the validated core (§5) and the speculative extension (§6). Conflating them is the single most likely way this project becomes scientifically indefensible.

## 4. Realistic accuracy expectations (no invented targets)

The original spec stated "expected R² 0.55–0.75" for behavioral traits. **We delete that.** It is a wish, not an estimate.

Prior: genomic prediction of *behavioral* traits in the DGRP is expected to be **modest**, because (a) the panel is small (n≈205), (b) behavioral traits are highly polygenic with substantial micro-environmental variance, and (c) much trait variance is non-additive. The honest position is: **the achievable accuracy is unknown until measured by nested cross-validation on this project's actual feature matrix.** We report whatever R²/MAE we measure, including if it is low. A low-but-honest R² with calibrated uncertainty is a publishable, useful result (it quantifies how much of trait variance is genomically predictable). An inflated R² is fraud.

**MEASURED (Findings 03, 2026-06-29).** First baseline (GBLUP, startle response, n=204): genomic
h² ≈ 0.27–0.33, but **out-of-sample CV R² ≈ 0.03 (Pearson r ≈ 0.16), significant (perm p=0.02)**,
with calibrated 90% intervals. This is **consistent with the canonical literature** (Ober et al.
2012, DGRP GBLUP startle r = 0.230 ± 0.012). Root cause: the DGRP is an *unrelated* association panel
(GRM off-diagonal SD ≈ 0.03), so relatedness-based genomic prediction has a low structural ceiling
(r ≈ 0.2 within-panel; ≈ 0 for novel unrelated genotypes). The engine is verified correct on
structured-kernel simulations (`tests/test_gblup.py`). This measured ceiling — not the brief's
fabricated "R²=0.55–0.75 / 74%" — is the realistic baseline for the genotype-only approach.

## 5. The validated core (v1 scope)

1. **Targets:** DGRP traits with measured line coverage, now quantified empirically (see `FINDINGS_coverage.md`, 2026-06-29). Verified coverage: lifespan 218*, startle 201, starvation 197, climbing 192, sleep 168; all-5-jointly = **106** lines (NOT the 205 the brief assumed). **Decision: build single-trait models** (each using its full coverage) rather than a forced 106-line unified model. First targets = the high-coverage behavioral cluster (startle/starvation/climbing, ~180–193 shared lines). *(\*lifespan's 218>205 means line sets must first be intersected with the canonical genotyped panel — Findings 02.)*
2. **Features:** (a) genotype-derived — pathway-collapsed variant burden + inversion status + Wolbachia status as covariates; (b) transcriptome-derived pathway-activity scores (ssGSEA) from Huang et al. 2015 per-line expression. Raw SNP matrices are never used directly (205 × 4.85M guarantees overfit).
3. **Validation:** `GroupKFold` by DGRP line ID (no line appears in both train and test). Report R², MAE, and a calibration plot for every target. No metric without all three.
4. **Uncertainty:** a *statistically valid* method (e.g. quantile regression, conformal prediction, or a Gaussian-process / Bayesian model that yields calibrated intervals) — **not** Gaussian noise injected into inputs. Intervals must be validated for coverage (a 95% interval must contain truth ~95% of the time on held-out lines).

## 6. The extrapolation layer (explicitly separated, NOT validated core)

Disease-model prediction and drug-rescue scoring live here. Rules:
- Always rendered to the user with an explicit "extrapolation — outside training domain" flag.
- Any drug-target feature must come from a **real** queried source (ChEMBL/DrugBank API), with the source ID stored. No hand-typed affinity scores.
- Any external disease-model "measured" value used for validation must be **extracted from the actual cited paper** (figure/table), with page/figure provenance recorded. **No invented validation values, ever.** Until extracted from source, an external-validation entry does not exist.

## 7. Corrections applied to the original CLAUDE_DroPhenoPredict.md

The following were identified as fabricated or invalid and are **explicitly excluded** from implementation:

| Original location | Issue | Resolution |
|---|---|---|
| `fetch_dgrp_phenotypes` fallback (L136–150) | Hardcoded invented DGRP phenotype values | Deleted. Real fetch from DGRPool API only; failure returns empty + error, never fake data. |
| `EXTERNAL_VALIDATION_DATASETS` (L664–699) | Invented `measured_*` values with real-looking citations | Deleted. Entries created only after extraction from source paper, with provenance. |
| `DRUG_PROFILES` (L303–330) | Invented drug-target scores | Deleted. Replaced by live ChEMBL/DrugBank query with stored IDs. |
| Bootstrap CI (L636–644) | Input-noise spread mislabelled as predictive CI | Deleted. Replaced by conformal/quantile/GP intervals with measured coverage. |
| Train/predict feature mismatch (L567–597) | Inference builds different features than training | Inference must build the *identical* feature transform as training, from the same code path. |
| "Expected R²" table (L812–816) | Aspirational targets stated as expectations | Deleted. We report measured values only. |
| Citation "Huang 2015 Genetics 200:1015" | Mis-citation of the transcriptome source | Corrected to Huang et al. 2015 *PNAS* 112:E6010. |

Rules retained from the original (these were correct): GroupKFold by line ID; SHAP interpretability; pathway-collapse over raw SNPs; read-data-first / no-integer-index discipline; reproducible seed=42; sklearn Pipeline objects.

## 8. Non-negotiable rules (inherited + strengthened)

1. Every feature has a biological justification recorded in `docs/FEATURES.md`.
2. Every prediction carries a calibrated, coverage-validated interval.
3. Train/test split respects line independence (GroupKFold by `DGRP_line`).
4. No metric reported without R² + MAE + calibration.
5. No value enters the repo unless downloaded, computed, or explicitly flagged placeholder.
6. Internal CV accuracy and any external accuracy are always reported separately and never conflated.
