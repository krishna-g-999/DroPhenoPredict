# DroPhenoPredict — Feature documentation (Charter §8.1)

Every feature used by the model, with its biological justification, source, and the
**measured** verdict on whether it helped prediction. Nothing is included on faith.

**Updated 2026-07-03** — this file originally documented only the v1.0 genotype-only model.
v1.1 (shipped, `app.py`'s actual predictor) is multi-modal: RNA-seq expression is now a
**retained, deployed** feature for several traits, not merely an excluded/future resource as
the first draft of this file said. Corrected below.

## Features IN the final model (v1.1, per-trait adaptive)

### Genome-wide additive relationship (GRM) — deployed for startle (both sexes), sleep (both sexes), starvation♀
- **What:** VanRaden genomic relationship matrix over **1,493,351 autosomal QC SNPs**
  (MAF≥0.05, call-rate≥0.8) across 205 DGRP lines.
- **Biological justification:** quantitative traits are additive-polygenic; the GRM encodes
  the expected sharing of causal-allele states between lines, the substrate of GBLUP.
- **Source:** DGRP Freeze 2.0 PLINK genotypes (Zenodo 5582846). `genotypes.build_grm`.
- **Verdict:** significant out-of-sample signal for 4/8 trait-sex models (r 0.10–0.25). Retained
  as the deployed modality for those traits (`multimodal.nested_modality_cv` picks per trait).

### RNA-seq expression relationship — deployed for climbing (both sexes), starvation♂
- **What:** linear relationship matrix over RNA-seq gene-expression profiles, GSE117850
  (Everett/Huang 2020; female 15,732 genes × 200 lines, male 20,375 × 200).
- **Biological justification:** expression is a heritable molecular intermediate closer to
  phenotype than sequence; this data type (RNA-seq, not the earlier microarray — see below)
  was the one that actually delivered gain (Findings 06).
- **Source:** GEO GSE117850. `expression.load_rnaseq`, `multimodal.expression_kernel`.
- **Verdict (Findings 06/09, MODEL_CARD):** significant out-of-sample gain over genotype for
  climbing♂ (0.10→0.27), climbing♀ (0.01→0.13), starvation♂ (0.27→0.34). **Retained** as the
  deployed modality for those 3/8 trait-sex models. Requires the query line's own RNA-seq
  profile (a usage precondition, not a limitation of the method — stated in `MODEL_CARD.md`).
- **Note — corrects this file's original claim:** the first draft of this document (pre-v1.1)
  said transcriptome data was "excluded from the predictive model, retained only as a future
  resource." That was true of the Huang 2015 **microarray** (below) but is no longer true of
  **RNA-seq**, which v1.1 actively deploys. Kept both entries below for the historical record of
  what was tested and why one succeeded where the other didn't.

## Features TESTED and EXCLUDED (measured to not help — documented per charter)

### Transcriptomic relationship (TRM) — line-level **microarray** expression (superseded by RNA-seq above)
- **What:** relationship matrix over 18,140 genes (Huang 2015 microarray, sex-matched).
- **Justification tested:** expression is a heritable molecular intermediate closer to phenotype.
- **Verdict (Findings 05):** no out-of-sample gain by three methods (TRM, G+T REML, supervised
  genes). **Excluded** from the predictive model. Superseded by the RNA-seq expression modality
  above, which measurably does help (Findings 06) — the negative result here is data-type
  specific (microarray dynamic range/noise), not evidence that expression in general is useless.

### Marker-effect (sparse QTL) features
- **What:** within-fold GWAS-selected top-500 SNPs + RidgeCV.
- **Justification tested:** a few large-effect loci might beat the infinitesimal GBLUP.
- **Verdict (Findings 05):** no systematic gain over GBLUP (confirms Ober 2012); lone lead
  starvation♀. **Excluded** from the default model.

### Inversion karyotype + Wolbachia status (fixed-effect covariates)
- **What:** polymorphic inversions In(2L)t, In(2R)NS, In(3R)Mo (0/0.5/1) + Wolbachia (0/1).
- **Justification tested:** known major cofactors of DGRP trait variation (Huang 2014).
- **Verdict:** adding them as fixed effects **reduced** CV r (Δr ≈ −0.01 to −0.05): inversions are
  large LD blocks already captured by the genome-wide GRM, so they are redundant. **Excluded**
  from the predictive model (available in `covariates.py` for GWAS-style confounder control).
- **Source:** Internet Archive snapshots of NCSU `inversion.xlsx`, `wolbachia.xlsx`.

## Design principle
The model is deliberately **minimal and honest**: it uses the one feature set that demonstrably
carries out-of-sample signal, after fairly testing the richer alternatives the original brief
proposed. Every exclusion above is a measured result, not an assumption.
