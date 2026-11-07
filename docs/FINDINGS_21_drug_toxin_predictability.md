# Findings 21 — Predictability of DGRP drug and toxin responses

**Date:** 2026-11-07 · **Scripts:** `45_collect_dgrpool_chemical.py`, `46_download_chemical_supplements.py`,
`47_download_repository_data.py`, `48_build_chemical_traits.py`, `49_chemical_predictability.py` ·
**Data:** `data/processed/chemical_traits_long.csv`, `chemical_traits_catalog.csv`,
`chemical_predictability.csv` · **Sources:** `docs/DATA_SOURCES.md` §7

## Question

Can the existing DroPhenoPredict engine predict how DGRP lines respond to drugs and toxins, and is
the line-specific drug *effect* (treated minus untreated) predictable?

## Data collected

69 trait–sex combinations for 21 compounds, 99–205 genotyped lines each (12,107 line values), from
DGRPool (15 studies), figshare (Battlay 2018, Duneau 2018), and paper supplements (Denecke 2017,
Mitchell 2017, Harrison 2024, Marconcini 2026). Seventeen paired designs gave derived drug-effect
traits (treated − control per line). 86 traits were screened in total. Sources that were found but
could not be used are listed with reasons in `DATA_SOURCES.md` §7.

Checks: no duplicate line–trait values; value ranges consistent with each assay (e.g. lead and ethanol
lower viability and lengthen development relative to their controls); deltamethrin survival
recomputed from raw records agrees with the DGRPool copy (r = 0.93).

## Method

Same engine and settings as v1.1 (Findings 06, MODEL_CARD):
- **Genotype:** GBLUP on all genotyped lines, 5-fold CV × 10 repeats, 200-permutation test.
- **Multi-modal:** nested CV choosing between the GRM and RNA-seq expression (same-sex expression for
  sex-specific traits; female and male expression both offered for traits with sexes pooled);
  permutation of the deployed kernel against the nested adaptive R².
- Benjamini–Hochberg FDR across the 86 traits, separately for each analysis. With 200 permutations
  the smallest p is 0.005, so the smallest attainable q is about 0.043.

## Results

**12 drug/toxin trait–sex combinations, covering 7 compounds, are predictable at FDR < 0.05** by
genotype, by the multi-modal procedure, or by both:

| Trait | Lines | Genotype r | Adaptive r (data used) | FDR q (geno / adaptive) |
|---|---|---|---|---|
| Methylphenidate effect on food intake (♂, treated − sucrose control) | 172 | 0.39 | 0.36 (geno) | 0.043 / 0.043 |
| Parathion survival at 48 h | 194 | 0.37 | 0.32 (geno) | 0.043 / 0.043 |
| Paraquat resistance (♀) | 167 | 0.03 | 0.34 (expression ♀) | 0.32 / 0.043 |
| Alcohol tolerance (♂) | 205 | 0.32 | 0.31 (geno) | 0.043 / 0.043 |
| Ethanol, developmental viability | 201 | 0.31 | 0.26 (expression ♀) | 0.043 / 0.076 |
| Lead, developmental viability | 199 | 0.31 | 0.25 (geno) | 0.043 / 0.043 |
| Methylphenidate food intake (♂) | 172 | 0.26 | 0.28 (geno) | 0.043 / 0.043 |
| Malathion 12 h mortality (♂) | 170 | 0.26 | 0.21 (geno) | 0.043 / 0.16 |
| Octanoic acid (noni toxin) survival (♂) | 184 | 0.25 | 0.24 (geno) | 0.043 / 0.043 |
| Methylphenidate locomotor activity (♂) | 172 | 0.25 | 0.21 (geno) | 0.043 / 0.066 |
| Alcohol tolerance (♀) | 205 | 0.24 | 0.22 (geno) | 0.043 / 0.066 |
| Malathion 12 h mortality (♀) | 170 | 0.20 | 0.23 (geno) | 0.080 / 0.043 |

Two control traits also passed (ethanol-study control viability, adaptive r = 0.42 using female
expression; menadione-study control startle ♀, r = 0.20); they are not drug responses.

**Not predictable** (q > 0.1 in both analyses): deltamethrin, azinphos-methyl, ryanoid insecticide,
imidacloprid 25 ppm, α-amanitin, rapamycin, caffeine lifespan, boric acid, methylmercury, menadione
resistance, lead development time, adult alcohol sensitivity, and low-dose paraquat climbing (99 lines).

**Drug-effect traits:** of 17 treated-minus-control traits, only the methylphenidate effect on food
intake was predictable. The other 16 (menadione climbing and startle, lead, ethanol development,
methylmercury, caffeine rescue of methylmercury, paraquat climbing, rapamycin) were not. A difference
of two noisy measurements carries the error of both, which lowers its predictability.

## Interpretation

1. The engine transfers to chemical-response phenotypes: accuracy for the best drug traits
   (r ≈ 0.25–0.39) is in the same range as the best v1.1 traits (r ≈ 0.25–0.34).
2. Data type matters per compound, as for v1.1: paraquat resistance in females is predictable only
   from RNA-seq expression (genotype r = 0.03, expression r = 0.34), matching the earlier pattern for
   climbing and male starvation.
3. Strain-specific drug *effects* are mostly not predictable with these sample sizes. The one
   exception (methylphenidate on feeding) is a drug acting on behaviour, measured with a matched
   control in the same experiment.
4. High REML heritability does not imply predictability (e.g. lead development-time effect, rapamycin
   control), the same unrelated-panel limit described in Findings 03–04.

## Caveats

- Each source measured its trait under its own protocol; values are not comparable across studies.
- Traits with sexes pooled ('B') are compared with expression from each sex separately.
- Several traits are bounded or skewed (mortality proportions, LD50); they were analysed untransformed,
  as in v1.1.
- Multiple-testing control is within this screen only (86 traits), consistent with the project's
  per-analysis FDR policy (AUDIT §3).
- REML h² at the boundary (1.0 or 0.0) reflects estimation at n ≈ 170, not exact heritability.

## Next steps (drug work)

- ~~Package the 12 predictable traits as an extension model~~ — done, see update below.
- Use the predicted resistant/susceptible lines to choose strains for a matched drug screen.
- Connect to the v2 convergence signature and to drug transcriptional signatures (CMap/LINCS via
  human orthologs; Karam 2022 amphetamine response in flies) for candidate-compound ranking.

## Update — packaged as drug-v1.0 (2026-11-07, script 50)

The 12 traits in the table above are packaged as `DroPhenoPredictor("drug-v1.0")`, in the v1.1
file format. Measured 90% interval coverage: 0.896–0.923. For parathion and malathion (both sexes)
the median interval covers 87–98% of the possible range, so those predictions carry an
`interval_note` and should be used to rank lines only. Details and per-trait figures: MODEL_CARD,
"Extension model — drug-v1.0". Tests: `tests/test_drug_model.py` (9 tests; suite total 68).
