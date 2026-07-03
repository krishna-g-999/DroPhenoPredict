# DroPhenoPredict — scientific roadmap (how to take this forward)

The guiding principle, stated once: **integrate matched, per-line data on the SAME panel — not a pile
of heterogeneous studies.** Supervised phenotype prediction needs (feature, phenotype) pairs on the
same animals; most disease/metabolism studies have molecular data but no matched behavioural
phenotype, and pooling them naively drowns biology in batch effects and domain shift. Every direction
below respects that principle, and each cites data we have *verified exists*.

## Tier 0 — done (this work)
Multi-modal genotype + RNA-seq prediction (v1.1), calibrated intervals, predictability flags,
adjusted GWAS, QTT gene association, GO pathway→trait map, web app. Three modern techniques tested
and found NOT to help (deep learning, non-linear kernels, functional stratification); the one lever
that helps is the RNA-seq modality. Mechanism: lipid β-oxidation ↔ starvation (genes + pathways
converge).

## Tier 1 — metabolome as a third modality  [DONE — honest negative, Findings 11]
**Executed and verified (2026-07-01).** Downloaded + processed all 673 MTBLS2060 NMR spectra → 170
lines × 475 bins (PQN+log). The metabolome is heritable (within-line replicate r=+0.20) but its
relationship matrix is **orthogonal to the cross-study DGRPool phenotypes** (h²≈0; metab CV r
negative for all traits; did not reproduce the paper's activity PA>0.4). **Not shipped** — could not
demonstrate benefit. Likely cause: the paper predicted MATCHED-sample phenotypes; cross-study transfer
of the condition-labile metabolome fails where the transcriptome succeeded. This is a cautionary,
integrity-preserving result that directly informs Tier 2. See `docs/FINDINGS_11_metabolome.md`.

### (original rationale, retained)
- **Why:** the metabolome is the layer closest to phenotype and to disease metabolism. Published
  precedent: *Prediction of complex phenotypes using the Drosophila metabolome* (Heredity 2021)
  shows **MBLUP** (metabolomic relationship matrix → GBLUP, exactly our framework) predicts
  **locomotor activity at PA>0.4 vs genomic <0.1** — it wins for traits genotype/expression fail on.
- **Data (verified):** MetaboLights **MTBLS2060** — NMR metabolome, **170 DGRP lines**, 4 reps,
  sample→line map present (`Factor Value[DGRP line]`). Also Zhou 2020 (Genome Res) 40-line MS
  metabolome (named metabolites; smaller n).
- **Work required:** MTBLS2060 ships **raw NMR only** (no processed matrix); needs spectral
  processing (read spectra → bin/align → normalize → lines×bins matrix) before building the
  metabolomic relationship matrix. Then add as a 3rd modality in the nested-CV adaptive framework.
- **Expected payoff:** new predictable traits (locomotor activity), and a metabolite-level mechanism
  complementing the lipid-pathway story.

## Tier 2 — the disease/drug aim, done correctly  [VALID design, data-limited]
The original vision (predict transgenic disease models / drug rescue) is reachable ONLY as a
**transfer-validation experiment**, not by training on DGRP natural variation:
1. Train expression→phenotype on the DGRP (have it).
2. **Curate external Drosophila datasets that have BOTH an expression signature AND a measured
   behavioural phenotype** on the same animals (e.g. neurodegeneration models with climbing + RNA-seq;
   drug-treatment cohorts). Phenotype values must be EXTRACTED from the papers with provenance —
   never fabricated.
3. Apply the DGRP-trained model to these external signatures and test whether predictions match the
   measured phenotypes (directional accuracy). This DIRECTLY tests whether natural-variation models
   transfer to perturbations — the make-or-break question for the disease aim.
- **Reality check:** the number of matched (expression + behaviour) fly datasets is small and
  heterogeneous → this is a careful curation + meta-analysis, with a real chance of a negative
  (non-transfer) result. That negative would itself be an important, publishable finding.

## Tier 3 — signature-based compound/perturbation scoring (CMap-for-flies)  [future]
Represent each perturbation by its expression signature; score it against the DGRP-derived
gene→phenotype weights (QTT) to predict phenotypic direction. Feasible only after Tier 2 establishes
that the DGRP gene→phenotype map transfers. Otherwise it is speculation.

## What we will NOT do (and why)
- Pool heterogeneous GEO studies into one training set (no common label; batch/domain confounding).
- Use the user's human/mouse ALS omics as per-fly features (cross-organism category error).
- Claim disease/drug prediction before Tier 2 transfer is demonstrated.

## Immediate next action (recommended)
Tier 1: process MTBLS2060 NMR → metabolomic relationship matrix → add as the 3rd adaptive modality.
This is concrete, the data is in hand, the method is published, and it targets exactly the
metabolism/disease-relevant traits.
