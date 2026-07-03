# Findings 08 — Trait-panel expansion: tested, no new predictable traits

**Date:** 2026-06-29 · **Code:** `scripts/17_expand_traits.py` · **Output:** `data/processed/expanded_traits.csv`

## Question
Can the resource be broadened with additional predictable DGRP traits beyond the 4 core families?

## What was tested
Scanned DGRPool numeric-continuous phenotypes for new biologically-distinct families (oxidative
stress, chill-coma, activity, feeding, body weight, triglyceride, glucose, ethanol, lifespan) with
≥150 lines overlapping genotype + RNA-seq, then ran genotype vs expression CV + permutation test.

## Result
- Most candidate families (paraquat, chill-coma, body weight, triglyceride, glucose) **lack ≥150
  lines** overlapping genotype + RNA-seq — coverage is the binding constraint.
- The families that did qualify were **not predictable**: oxidative-MSB (geno r=0.064, p=0.075),
  locomotor activity (geno r=0.089, p=0.10), lifespan (exp r=0.053, p=0.085).
- **0 newly predictable traits.**

## Conclusion (honest)
Expanding the trait panel does **not** add predictable models with current public data. The 6
predictable trait-sex models (climbing, startle, starvation, sleep) are essentially the predictable
subset the DGRP supports at n≈190. The resource is therefore **appropriately scoped** — we do not pad
it with non-predictable traits. Broader prediction would require larger phenotyped+omics panels, not
more traits on these lines.

This is consistent with Findings 04 (most DGRP organismal traits have low genomic predictability) and
the unrelated-panel ceiling (Findings 03).
