# DroPhenoPredict

A genotype/transcriptome → phenotype prediction model for the **Drosophila Genetic Reference Panel (DGRP)**, with honestly-reported cross-validated accuracy and calibrated uncertainty.

> **Read [`docs/SCIENTIFIC_CHARTER.md`](docs/SCIENTIFIC_CHARTER.md) before contributing.** It is the governing document. The core rule: **no fabricated data, ever** — every value is downloaded, computed, or explicitly flagged as a placeholder.

## What it does (validated core)
Predicts DGRP physiological/behavioral traits **that have measured ground truth in the panel** (e.g. startle locomotion, sleep, lifespan, starvation resistance) from genotype- and transcriptome-derived pathway features. Reports R², MAE, calibration, and coverage-validated prediction intervals using line-aware cross-validation.

## What it does NOT (yet) do
It does **not** make validated predictions for transgenic human-disease models (TDP-43/FUS/SOD1/HTT) or drug rescue. The DGRP contains no such genotypes, so those are a separate, explicitly-labelled **extrapolation layer** (Charter §6) — hypothesis-generation, not prediction.

## Status — genotype→phenotype predictability fully characterized (Findings 03–05)

Validated GBLUP pipeline + a complete, honest evaluation of what is predictable in the DGRP:

- **Genotype (GBLUP)** is the best predictor; ceiling **r ≈ 0.2** (starvation, startle, sleep♀
  significant), reproducing Ober et al. 2012. The DGRP is an *unrelated* panel → low structural
  prediction ceiling (engine verified on structured-kernel sims; `tests/` 7 pass).
- **Option 3 — predictability map** across 5 traits × sex (`Findings 04`): a real trait-specific
  gradient (starvation/sleep♀ > startle > climbing/lifespan≈0); calibrated 90% intervals throughout.
- **Option 1 — transcriptome** (real per-line expression recovered via Internet Archive): adds **no**
  out-of-sample gain over genotype, by three methods (TRM, G+T REML, supervised genes). `Findings 05`.
- **Option 2 — marker-effect models**: do **not** systematically beat GBLUP (confirms Ober 2012);
  lone lead = starvation♀. `Findings 05`.

**Model v1.1 built & packaged — MULTI-MODAL** (`scripts/13_finalize_multimodal.py`): per trait it
adaptively uses genotype (GRM) or RNA-seq expression (GSE117850), chosen by **nested CV** (unbiased).
8 trait-sex models, **6 predictable** (perm p<0.05), calibrated 90% intervals (coverage 0.90–0.92).
Headline = nested adaptive r. RNA-seq expression rescues the genotype-unpredictable traits:
**climbing♂ 0.10→0.27, climbing♀ 0.01→0.13, starvation♂ 0.27→0.34** (reproduces Morgante 2020).
Deep learning / non-linear kernels tested and do NOT help (Findings 06). Inference API:
`DroPhenoPredictor("v1.1").predict_profile("DGRP_357")`. See `docs/MODEL_CARD.md`, `docs/FINDINGS_06_*`.

(v1.0 genotype-only model retained: `scripts/10_finalize_models.py`, 10 models.)

**Web app:** `./.venv/Scripts/streamlit run app.py` — per-line phenotype prediction with calibrated
intervals + modality + predictability, benchmark view, and per-trait GWAS/QTT candidate genes.

**Genetic architecture & advanced methods** (Findings 07–09): inversion/Wolbachia-adjusted GWAS
(robust candidate genes — sleep♀ *Vmat*, climbing♂ *pdm3/shep*, starvation *Ptp61F/fz2*);
transcriptome-wide QTT association (starvation lipid genes *Lip4/Lsd-1*, FDR<0.05); and three
controlled **negative** results — deep learning, non-linear kernels, and functional stratification all
fail to beat linear GBLUP at this n. The only accuracy lever is the RNA-seq modality.

**Bottom line:** the honest product is a *calibrated genomic-prediction resource with intervals* for
the predictable DGRP traits — useful for prioritizing which lines to phenotype. Not a high-accuracy
point predictor; accuracy is structurally limited by the unrelated ~200-line panel, and neither
omics nor marker-selection nor inversion/Wolbachia covariates raised it (all tested, all documented).

Reproduce: `scripts/03_build_grm.py` → `04`…`10` → `pytest tests/` (41 tests). Predict:
`DroPhenoPredictor("v1.0").predict_profile("DGRP_021")`.

## v2 — Disease-convergence resource (separate from the DGRP predictor above)

Uses a real matched multi-disease fly panel (BCM-DMAS, Synapse syn34767207 — Farrell/Shulman,
*npj Parkinson's Disease* 2025): head RNA-seq for **AD (Aβ42, Tau), PD (αSyn), HD (HTT-128Q, HTT-200Q)**
vs. a shared driver control, days 2–57. Tests whether neurodegeneration models molecularly converge
(the fly analog of Ruffini et al. 2020, *Cells* — shared human NDD biology).

**Result:** significant convergence at both the gene level (27 genes DE in ≥3/5 models, perm p=0.0005)
and, more strongly, the pathway level (60 GO-BP pathways FDR<0.1 in ≥3/5 models via a locally-built,
verified GSEA engine). The concordant core: **nucleolar stress / rRNA biogenesis is down-regulated**
in AD/PD/HD-128Q. **HTT-200Q is a consistent directional outlier** across three independent analyses.
Convergence survives leave-one-out (excluding any single model). DE pipeline independently
cross-validated against the source paper's own DESeq2 output (7,018 genes, r=0.77).

**This is a molecular-convergence resource, not a behaviour predictor** — matched climbing/lifespan
values for this panel exist only in the source paper's figures, not as data. See
`docs/FINDINGS_12`–`15`, `docs/AUDIT.md` §7 for the full, honest gap list (ALS/SCA not yet
incorporated, human↔fly ortholog validation pending, proteomics not yet pulled in).

## Layout
```
docs/    SCIENTIFIC_CHARTER.md, DATA_SOURCES.md  (governance — read first)
src/drophenopredict/   package code
scripts/ reproducible data-fetch / training entrypoints
data/    raw/ interim/ processed/ external_validation/  (raw is download-only, never edited)
models/  trained/ configs/
tests/   correctness + leakage guards
```

## Provenance
This project supersedes the draft spec in `CLAUDE_DroPhenoPredict.md`. That draft contained fabricated placeholder values (invented phenotype/drug/validation numbers); those are catalogued and excluded in Charter §7 and `DATA_SOURCES.md`. The valid engineering rules from it (GroupKFold by line, SHAP, pathway-collapse, reproducibility) are retained.
