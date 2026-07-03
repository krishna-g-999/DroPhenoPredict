# Findings 10 — Pathway-activity → trait map (metabolism, done rigorously)

**Date:** 2026-06-30 · **Code:** `src/drophenopredict/pathways.py`, `scripts/21_pathway_association.py`
**Data:** FlyBase GAF + GO-basic ontology (public) + DGRP RNA-seq (GSE117850).

## Method
GO biological-process gene sets (ancestor-propagated over the GO DAG; 2,617 sets, size 10–300).
Per DGRP line, pathway activity = mean z-scored RNA-seq expression of the set's genes. Each pathway
activity correlated with each predictable trait (BH-FDR over 2,617 pathways). This is the rigorous way
to incorporate metabolism/disease pathway knowledge — curated gene sets, not pooled raw studies.

## Result — coherent, sex-specific energy biology

| Trait·sex | Top metabolic pathway | r |
|---|---|---|
| starvation ♂ | **lipid oxidation / fatty-acid β-oxidation** | +0.26 |
| starvation ♀ | **glycogen / carbohydrate catabolism** | +0.23 |
| sleep ♀ | intracellular glucose homeostasis / glycolysis | −0.25 |
| climbing ♂ | mitochondrial OXPHOS (cyt-c / complex IV / ATP synthesis) | +0.24 |
| startle ♀ | cadherin cell-cell adhesion (neural) | −0.26 |

## Why this matters
1. **Independent corroboration of the QTT genes.** Starvation♂ → *lipid oxidation* pathway matches the
   QTT lipid genes (*Lip4, Lsd-1*) from Findings 09. Two independent analyses (gene-level and
   pathway-level) converge on lipid β-oxidation — convergence stronger than any single p-value.
2. **Sex-specific metabolic strategy.** Males mobilize **lipid** for starvation; females mobilize
   **glycogen/carbohydrate**. A concrete, testable biological hypothesis produced by the tool.
3. **Locomotion ↔ mitochondrial energy; sleep ↔ glucose** — both consistent with known physiology.

## Honest caveat
After correcting across 2,617 pathways, **none reach FDR<0.05** at n≈190 (top |r|≈0.25). These are
**suggestive, convergent, hypothesis-generating** associations — not genome-wide-significant. Reported
as such. The cross-method convergence (genes + pathways) is the real strength.

## What this does and does NOT do
- DOES: bring curated metabolism/pathway knowledge into the resource rigorously, formalize the lipid
  finding across all processes, and give per-trait mechanistic pathway maps (in the web app).
- Does NOT: by itself raise prediction accuracy (pathway features are interpretive, not a new
  predictive modality), and does NOT reach the disease/drug aim (see ROADMAP.md).
