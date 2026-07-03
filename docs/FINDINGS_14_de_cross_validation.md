# Findings 14 — Independent cross-validation of the DE pipeline against the paper's own results

**Date:** 2026-07-01 · **Code:** ad hoc verification script (folded into `scripts/33_shared_biology.py`
methodology) · **Data:** BCM-DMAS raw counts (ours) vs the paper's Supplementary Table 4 (theirs).

## Why this check matters
Findings 13's convergence result rests entirely on our own differential-expression (DE) pipeline. The
DE math was verified against `statsmodels` (identical output), but that only proves the arithmetic is
correct — not that the underlying design is a reasonable/valid genotype effect. This check asks a
stronger question: does an **independently computed** result agree with the **paper authors' own**
published DE, computed by a **different statistical method**, on the same raw data?

## Method
- **Ours:** linear model on log2-CPM, `Expr ~ 1 + disease + age_dummies` (age as a categorical fixed
  effect), aSyn vs Elav-Gal4 control, OLS by matrix algebra.
- **Theirs (Sup. Table 4):** DESeq2 negative-binomial GLM, `Expr ~ Age + Genotype + Genotype*Age`
  (age as continuous, explicit genotype×age interaction term), likelihood-ratio test.
- **Critical catch:** the paper's `Gene` column is a **gene SYMBOL** (e.g. `abd-A`), while our counts
  matrix is indexed by **FBgn ID** (from featureCounts). A naive join on the raw index only matched 153
  genes (stray FBgn-format entries in the paper's table) — an artifact, not a real comparison. Fixed by
  mapping our FBgn IDs to symbols via the existing FlyBase-derived annotation
  (`annotation.load_gene_spans`) before joining.

## Result (corrected, real comparison)
- 9,385 expressed genes (ours) → 9,091 resolved to unique symbols → **7,018 genes in common** with the
  paper's table (of their 16,602).
- **Pearson r = 0.766, Spearman r = 0.607 (both p ≈ 0)** between our log2FC and theirs.
- Median absolute discrepancy: 0.108 log2FC units.
- Regression slope (ours ~ theirs) = 0.69 — same direction, moderately attenuated (expected: two
  different statistical frameworks — OLS-on-log-CPM vs DESeq2-NB-GLM with an explicit interaction term
  — will not agree exactly, but strong concordance confirms the same underlying signal).

## Conclusion
The DE pipeline that feeds Findings 13's convergence result and the GSEA analysis (script 35) is
**independently validated**: a materially different statistical method, applied by different people
(the original authors, via DESeq2) to the same raw counts, recovers highly correlated gene-level
effects. This is strong, real evidence that our downstream convergence/GSEA results are not an artifact
of our specific DE implementation.

## Process note (transparency)
This check initially reported a spuriously strong result (153 genes, r=0.83) due to a gene-ID
format mismatch that was not caught before the first run. It was diagnosed, the join was fixed, and the
corrected, larger-n comparison (7,018 genes) is what is reported above. Recorded here per the project's
verify-every-step standard — mistakes get caught and corrected, not hidden.
