# Findings 18 — Comprehensive omics-data landscape search (genomic, epigenomic, metabolomic, other)

**Date:** 2026-07-02 · Systematic, verified search for additional matched neurodegenerative-disease
Drosophila omics data, per the project's standing rule: verify before recommending, never grab data
that would reintroduce the batch-effect risk already demonstrated (Findings 11).

## 1. Proteomics (BCM-DMAS) — DONE
See Findings 17. Fully integrated, real RNA↔protein concordance validated (r=0.33–0.36, p<10⁻¹⁵⁰).

## 2. Genomic layer — MAJOR FIND: Yang et al. 2023 DGRP×AD eye-degeneration screen
**Already fully downloaded, no DUA, no further user action needed** (public GitHub, EuropePMC
supplementary API).

- **Source:** Yang, Zinkgraf, Fitzgerald-Cook, Harrison, Putzier, Promislow, Wang. 2023.
  *"Using Drosophila to identify naturally occurring genetic modifiers of amyloid beta 42- and
  tau-induced toxicity."* G3 (Bethesda) 13(9):jkad132.
- **Design (verified from the paper's own Methods, quoted directly):** "Experimental F1-generation
  flies were produced by crossing unmated virgin females from 162 DGRP lines to males from the
  GMR > Aβ42;tau (R32) line" — a **double Aβ42+Tau transgenic AD model**, crossed onto DGRP genetic
  backgrounds (F1, so the phenotype reflects a **heterozygous/dominant** background effect — a real,
  important caveat, different from our homozygous-line GBLUP model).
  Eye degeneration scored at day 28 via 16 automated image-morphology features (nearest-neighbor
  distance, circularity, ommatidial area/perimeter/radius statistics).
- **Genotypes:** the paper independently downloaded and QC'd the **same DGRP2 PLINK source** we
  already have fully processed (identical thresholds: MAF≥0.05, call-rate≥0.8) — meaning our existing
  GRM (`data/processed/grm.npz`) can likely be reused directly, restricted to their 162-line subset, with
  no new genotype processing needed.
- **Data obtained:** `data/raw/yang2023_github/DGRP_rawUntransformed_16traits.rdata` — 5,256 individual
  eye images, 164 DGRP lines, 3 batches (2 internal control lines "32C"/"441C" per batch for QC).
  Not yet aggregated to per-line means/BLUP (that's the next step, using the same line-mean/GBLUP
  machinery already built in this project).
- **Why this matters:** this is the **first real, matched (DGRP genotype ↔ disease phenotype) dataset**
  found in this entire project. Everything else either lacked genotypes (BCM-DMAS engineered
  transgenics on one background), lacked matched behavior (BCM-DMAS climbing data locked in a figure),
  or lacked a disease phenotype (DGRP metabolome — natural variation only). This closes that specific
  gap, with the honest caveat that it is a heterozygous F1-cross effect, not the additive homozygous
  model our core GBLUP pipeline assumes — integration must account for that, not ignore it.

## 3. Cross-lab Drosophila multi-omics — EmoryDrosophilaTau (verified real, recommend download)
Found via systematic scan of the AD Knowledge Portal / Synapse consortium (parent syn5550383) for
Drosophila-named projects (only 2 of ~90 sibling projects are fly-based; the rest are mouse/iPSC/cell
lines — verified by name, not assumed).

- **syn7274101**, from a different lab (Emory) than BCM-DMAS (Baylor/Shulman) — genuinely independent
  data, not just another arm of the same experiment.
- **Contents (verified via folder structure):** bulk RNA-seq, **scRNA-seq**, and label-free
  proteomics, with full metadata (biospecimen, per-assay). This would give **cell-type resolution**
  BCM-DMAS's own (unlocatable) scRNA-seq couldn't provide, plus a genuine **cross-lab replication
  opportunity** for Tau-specific findings (a real strengthening of Findings 13/15's Tau results,
  which currently rest on one lab's data).
- **Access:** same DUA pattern as BCM-DMAS (403 unauthenticated) — requires you to accept terms on
  the Synapse page and download via your own session, same as before.
- **Recommendation:** worth pursuing if you want to strengthen the Tau-specific convergence findings
  with independent replication; not essential to the core project.

## 4. EmoryFly (syn52223842) — checked, currently empty
Data folder verified empty (recently modified, 2026-02-24) — likely an in-progress/embargoed project.
Not usable now; noted for awareness, not recommended for action.

## 5. Epigenomic Drosophila neurodegeneration data — searched, no usable dataset found
Checked the most promising lead (α-synuclein H3K9me2 chromatin paper, PMC5093762). **Verified
negative:** the H3K9me2 data is from targeted ChIP-qPCR at ~29 candidate loci plus western blot —
**not genome-wide** (no ChIP-seq/CUT&RUN), **not deposited** in any repository, and **no behavioral
data**. No other genome-wide epigenomic (ATAC-seq/ChIP-seq/methylation) Drosophila neurodegeneration
dataset with a public deposit was found. Honest gap, not fabricated or padded.

## 6. Disease-model metabolomics (distinct from the DGRP natural-variation metabolome already tested)
Checked the most promising lead (Drosophila HD-model 1H-MAS-NMR paper, PMC7289880: 2 HD genotypes,
real sample sizes n=8–11/group). **Verified negative:** the NMR data was **never deposited** to
MetaboLights, Metabolomics Workbench, or any other public repository — only summarized in figures/a
supplementary docx. No behavioral data on the same flies either. Not usable; another honest gap.

## 7. Conceptual clarification — "genomic data" for the disease-model panel
Unlike the DGRP (natural variation across 205 wild-derived lines, where "genomic data" means variant
calls across many diverse genotypes), the BCM-DMAS disease models are **engineered transgenics on one
controlled genetic background** — there is no natural-variation genomic layer to seek *for that
specific panel*; the disease-relevant "genomic" question only becomes meaningful when a disease
transgene is crossed onto a *diverse* background, which is exactly what the Yang et al. 2023 dataset
(§2) provides. This reframes what was being searched for, rather than continuing to look for something
that doesn't fit the BCM-DMAS panel's design.

## Summary table

| Layer | Status | Action needed |
|---|---|---|
| Proteomics (BCM-DMAS, PD) | ✅ Done, validated (Findings 17) | none |
| Genomic × disease phenotype (Yang 2023, AD) | ✅ **Downloaded**, real, verified | analysis pending (needs a design decision — see below) |
| Cross-lab multi-omics (EmoryDrosophilaTau, Tau) | Verified real, not downloaded | you download via Synapse if wanted |
| scRNA-seq (BCM-DMAS) | Referenced in metadata, data folder not locatable | none currently actionable |
| Epigenomics (any disease) | Searched, none usable found | none — honest gap |
| Disease-model metabolomics | Searched, none usable (not deposited) | none — honest gap |
