# Findings 01 — Real DGRP trait coverage (gating question)

**Date:** 2026-06-29 · **Script:** `scripts/01_measure_coverage.py` · **Data:** `data/processed/coverage_report.csv`
**Source:** DGRPool live JSON API (VERIFIED). 903 phenotypes total, 595 numeric-continuous usable; 56 matched our 5 target families.

## Result

Best-covered phenotype per target family and per-line coverage:

| Family | Phenotype (id) | Lines |
|---|---|---|
| lifespan | `mn_Longevity` (1315) | 218* |
| locomotion_startle | `StartleRes` (2799) | 201 |
| starvation | `StarvationRes` (2798) | 197 |
| climbing | `mn_NegGeotaxis_CTRL` (1543) | 192 |
| sleep | `mn_SleepLength_Night` (1483) | 168 |
| **All 5 jointly** | — | **106** |

Pairwise overlap (lines measured in both):
```
climbing & lifespan           = 124      lifespan & sleep            = 113   <- joint bottleneck
climbing & locomotion_startle = 181      lifespan & starvation       = 128
climbing & sleep              = 167      locomotion_startle & sleep  = 163
climbing & starvation         = 177      locomotion_startle & starv. = 193   <- tightest pair
lifespan & locomotion_startle = 129      sleep & starvation          = 160
```

## Interpretation (drives Charter §5 scope)

1. **The brief's "all 205, all traits" is false.** A unified 5-trait set is **n = 106**.
2. **Lifespan is genetically/experimentally disjoint** from the behavioral traits (different studies, different line panels) — it depresses the joint set the most.
3. **Behavioral traits cluster** (startle/starvation/climbing share ~180–193 lines), suggesting a coherent behavioral sub-model with larger n.

### Scope decision
**Build single-trait models, not a forced unified model.** Each single-trait model uses that trait's full coverage (168–218 lines) instead of throwing away data to satisfy a 106-line intersection. A unified multi-output model is optional and secondary; it is not required and costs ~half the data. v1 first target = the high-coverage behavioral traits.

### Required data-hygiene step (open)
`*` `mn_Longevity` reports **218 > 205** canonical lines. DGRPool aggregates across studies and may include lines outside the genotyped DGRP Freeze 2.0. **Before modeling, intersect every trait's line set with the canonical genotyped panel** (the ~205 lines for which we have variant calls). Real modeling-n per trait is therefore ≤ the counts above. Tracked for Findings 02.
