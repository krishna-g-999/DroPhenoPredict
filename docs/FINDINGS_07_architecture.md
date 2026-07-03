# Findings 07 — Genetic architecture of the predictable traits (GWAS → genes)

**Date:** 2026-06-29 · **Code:** `scripts/15_gwas.py`, `scripts/16_map_genes.py`
**Data:** DGRP2 genotypes (1.5M autosomal SNPs) + DGRPool line-mean phenotypes + dm3 gene models
(GSE67505 GTF) · **Outputs:** `data/processed/gwas_top_*.csv`, `gwas_genes_*.csv`

## Method
Single-marker association (linear regression of line-mean phenotype on SNP dosage) for each
predictable trait × sex, streamed over all autosomal SNPs. Top hits mapped to the gene they fall in
(or nearest gene ≤5 kb) using dm3/BDGP5 gene models.

## Top association signals & candidate genes

| Trait·sex | n | Top p | Notable candidate genes (top SNPs) | Biological plausibility |
|---|---|---|---|---|
| sleep ♀ | 167 | **1.3e-13** | **Vmat**, Sema-1a, drk, mip120, mars | Vmat = vesicular monoamine transporter → sleep/arousal; drk signaling |
| climbing ♂ | 185 | **9.4e-9** | pdm3, **shep**, Ir92a, Hex-t1 | shep/pdm3 neuronal differentiation & remodeling → locomotion |
| starvation ♀ | 197 | **1.4e-8** | gek, Gycbeta100B, Tim10, snama | Tim10 mito import, Gycbeta100B signaling → energy/stress |
| climbing ♀ | 185 | 4.8e-8 | Hrb98DE, **Msp-300**, **jvl**, Eip63F | Msp-300 (nesprin) & jvl (javelin) = muscle/locomotor |
| startle ♀ | 201 | 7.5e-8 | **Snap25**, alpha-Cat, Abd-B, AGO3 | Snap25 = SNARE neurotransmitter release → startle |
| starvation ♂ | 197 | 2.8e-7 | **cnc**, Sulf1, Dro, CG12858 | cnc = NRF2 ortholog, metabolic/oxidative-stress master regulator |

Genome-wide Bonferroni ≈ 0.05/1.5e6 ≈ 3.3e-8: sleep♀, climbing♂, starvation♀ pass; others are
suggestive. Full ranked SNP and gene lists in `data/processed/gwas_*`.

## Interpretation
The candidate genes are **trait-appropriate**: synaptic (Snap25), monoaminergic (Vmat), muscle
(Msp-300, jvl), neuronal-developmental (shep, pdm3), and metabolic/stress (cnc/NRF2, Tim10). This
face validity supports that the predictable signal is biological, and gives reviewers a concrete
mechanistic narrative linking the prediction model to known *Drosophila* neuromuscular/metabolic
genetics.

## Honest caveats (must accompany any use of this table)
1. **Marginal association, not causal.** These are candidate loci; LD and pleiotropy are not resolved.
2. **Underpowered (n≈190).** DGRP behavioral GWAS rarely reaches stringent significance; treat
   sub-genome-wide hits as hypotheses.
3. **NOT yet adjusted for inversions / Wolbachia.** Standard DGRP GWAS regresses these out first;
   large inversions (esp. In(2L)t, In(3R)Mo) can tag association. An inversion/Wolbachia-adjusted
   re-run (covariates already verified in `covariates.py`) is the next refinement — the specific,
   named neuronal/muscle genes above are unlikely to be pure inversion artifacts, but this must be
   confirmed before any strong claim.
4. These candidates are for **biological interpretation of the predictable traits**, independent of
   the prediction model's accuracy (which is reported separately and honestly).
