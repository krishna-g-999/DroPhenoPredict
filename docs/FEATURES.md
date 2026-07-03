# DroPhenoPredict — Feature documentation (Charter §8.1)

Every feature used by the model, with its biological justification, source, and the
**measured** verdict on whether it helped prediction. Nothing is included on faith.

## Features IN the final model

### Genome-wide additive relationship (GRM) — the predictive substrate
- **What:** VanRaden genomic relationship matrix over **1,493,351 autosomal QC SNPs**
  (MAF≥0.05, call-rate≥0.8) across 205 DGRP lines.
- **Biological justification:** quantitative traits are additive-polygenic; the GRM encodes
  the expected sharing of causal-allele states between lines, the substrate of GBLUP.
- **Source:** DGRP Freeze 2.0 PLINK genotypes (Zenodo 5582846). `genotypes.build_grm`.
- **Verdict:** the ONLY feature set that yields significant out-of-sample signal (r up to ~0.24).
  Retained.

## Features TESTED and EXCLUDED (measured to not help — documented per charter)

### Transcriptomic relationship (TRM) — line-level expression
- **What:** relationship matrix over 18,140 genes (Huang 2015 microarray, sex-matched).
- **Justification tested:** expression is a heritable molecular intermediate closer to phenotype.
- **Verdict (Findings 05):** no out-of-sample gain by three methods (TRM, G+T REML, supervised
  genes). **Excluded** from the predictive model. (Retained as a resource for the future
  perturbation/disease extension, Charter §6.)

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
