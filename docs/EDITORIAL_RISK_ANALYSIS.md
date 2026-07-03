# Editorial / reviewer risk analysis — anticipated before building the Yang 2023 integration

Written as the "chief editor + expert reviewer" pass required before adding a new, substantial piece
of science to this project. Each risk below is addressed by a specific, concrete design decision in
the pipeline (`scripts/40`–`43`), not by hoping it doesn't come up in review.

## Risks specific to the Yang 2023 DGRP×AD integration

### 1. "The phenotype is measured in F1 heterozygotes, not homozygous DGRP lines — is this the same genetic model as your core predictor?"
**No, and we say so explicitly.** The F1 progeny (DGRP mother × R32 GMR>Aβ42;tau father) carry one
allele from the DGRP parent and one from the constant R32 background at every locus. This is the
**standard, well-precedented design** for DGRP modifier screens (used by the original Mackay-lab
papers this whole project already relies on for the AD/PD DGRP-crossing literature) — the phenotype is
modeled as attributable to the DGRP parent's *additive* genetic background (general combining ability),
with the R32 background held constant across all crosses. Our GRM (built from the DGRP lines'
own genotype calls) is the correct predictor input for this design — we are not changing the encoding,
only being explicit that the **target trait's biological meaning** differs from our core DGRP
behavioral traits (measured directly on homozygous lines). Stated in the model card / findings doc,
not buried.

### 2. "Did you replicate the original authors' own trait-construction method, or invent your own?"
**Replicated theirs, verified.** The original RMD shows the exact mixed model used:
`trait ~ group + (1|line)` (REML), extracting the random-effect BLUP per DGRP line, with group (26
experimental batches) as a fixed effect to remove technical/batch variance. We reproduce this design
with `statsmodels.MixedLM` rather than inventing a simpler ad hoc line-mean. Where our BLUP values or
heritability estimates can be compared to anything reported in the source materials, we do so (same
verify-against-source discipline as Findings 14).

### 3. "Internal QC control genotypes (classes '32C'/'441C') — were they correctly excluded from the genetic mapping?"
**Yes, explicitly filtered before model fitting**, not after. These are recurring internal reference
crosses (not DGRP genetic variation) used by the original authors to help stabilize batch-effect
estimation across groups. We exclude them from the phenotype-target construction from the start — a
design choice we make explicitly rather than including them and hoping they don't contaminate a later
join. (Real risk we identified ourselves: DGRP line-32 and line-441 are also genuine DGRP line
numbers elsewhere in the panel; failing to exclude the control classes could silently substitute
control-cross data for real DGRP-32/441 genetic effects in a downstream join. Documented so this
specific failure mode is auditable.)

### 4. "Multiple testing — you have 16 correlated image features. Did you fish for a hit?"
**No — we pre-commit to the original authors' own "primary" trait designation** (the first 5 of the 16
features: nn.mean, nn.sd, cc.mean, cc.sd, mean.s.area — their own labeling, visible in the source RMD,
not ours) as the primary analysis set, with Benjamini-Hochberg correction across them. The remaining 11
"secondary" traits are reported as exploratory/supplementary, clearly labeled as such, not used to
cherry-pick a headline number.

### 5. "Batch/day confounding?"
Checked directly: all `class=='line'` rows are `day==28` (single timepoint for the main GWAS-mapping
population — the day-4-vs-28 comparison in the original study used a separate small subset, not part of
the main mapping population). No day confound to adjust for. Batch/group (26 levels) is handled via
the mixed-model fixed effect exactly as above.

### 6. "Is this genuinely new science, or just re-running someone else's GWAS?"
The original paper did **single-marker GWAS** (per-SNP association) on this phenotype. We are doing
something they did not: (a) testing whether our **already-validated GBLUP framework** (genome-wide
additive prediction, verified via structured-kernel simulations and 3 literature-benchmark
reproductions in this project) can predict the degeneration phenotype **out-of-sample**, with proper
line-aware cross-validation and calibrated uncertainty; (b) testing whether **RNA-seq expression**
(a data type the original authors did not have for this cross) adds predictive power via our
multi-modal nested-CV framework. Both are genuinely new analyses on their real, public data — the same
pattern as every other extension in this project (build new, verified methodology on real data;
never just re-describe someone else's result as our own).

### 7. "Same unrelated-panel accuracy ceiling as everywhere else in your project?"
**Expected, and we will report honestly whatever we measure.** This is a ~150-160 line subset of the
same DGRP panel whose structural prediction ceiling (r≈0.1–0.3) is already established and reproduced
three times in this project (Ober 2012, Morgante 2020 ×1). We do not expect — and will not claim — a
qualitatively different result here just because the trait is disease-relevant. If the result is
modest, we report it as modest, exactly as everywhere else.

### 8. "Inversions confound DGRP population structure — did you check?"
The original authors' own PCA on this 162-line subset found In(3R)Mo and In(2L)t as dominant structure —
**they flagged this themselves.** We re-use our existing, tested `covariates.py` infrastructure
(inversion/Wolbachia design matrix) to run an adjusted robustness check, exactly as we already did for
our own DGRP GWAS (Findings 07/AUDIT §7-adjacent work) — not a new method, applying an existing verified
one to a new phenotype.

## Broader risks to the project's reliability (not specific to this integration)

### 9. "You keep finding new data sources — is the project's scope creeping without a clear endpoint?"
We hold the line: every new dataset is evaluated against the same fixed bar (matched design, verified
accessible, doesn't require the batch-effect-risking naive pooling we already proved fails) before
being integrated, and every integration produces a scoped, documented Finding — not an open-ended
data-collection exercise. The Scientific Charter and AUDIT.md remain the source of truth for what is
in scope.

### 10. "With this many analyses (GWAS, GBLUP, GSEA, convergence, proteomics, now Yang2023), is there a multiple-comparisons problem across the WHOLE project's claims?"
Each analysis is scoped, pre-registered in its own Finding with its own correction, and none of them
individually claims significance by borrowing power from another. We do not pool p-values across
findings to manufacture significance. This is worth stating explicitly in the eventual manuscript's
statistical methods section.

### 11. "Third-party dependency / reproducibility risk (pyreadr, statsmodels, Miniconda)?"
Noted and mitigated: `pyreadr` had a transient network/certificate issue in this sandbox (unrelated to
the R data itself) and was resolved via an explicit PyPI index URL — recorded here so the fix is
reproducible, not a silent workaround.
