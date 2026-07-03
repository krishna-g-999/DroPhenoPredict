# Findings 09 — Advanced methods: adjusted GWAS, modern prediction techniques, web tool

**Date:** 2026-06-30 · **Scripts:** 18 (adjusted GWAS), 19 (stratified prediction),
20 (QTT), `app.py` (web UI). All results measured; nothing assumed.

## 1. Inversion/Wolbachia-adjusted GWAS (Part 2 refinement)
Frisch-Waugh-Lovell residualization on covariates X = [1, In(2L)t, In(2R)NS, In(3R)Mo, Wolbachia],
then partial association per SNP. Adjustment matters (marginal/adjusted top-200 overlap 116–181),
removing inversion-confounded hits. **Adjustment-robust candidate genes** (top, covariate-adjusted):
- **sleep♀:** Vmat, Sema-1a, CG13330/CG13331 (p=4.2e-13; 181/200 robust) — *Vmat = vesicular
  monoamine transporter, central to sleep/arousal.*
- **climbing♂:** pdm3, shep, Ir92a — neural development/differentiation.
- **starvation♂:** Ptp61F, fz2 — insulin/JAK-STAT and Wnt signaling (metabolic).
- **starvation♀:** gek, scaf6. **startle♀:** nwk, alpha-Cat (synaptic/adhesion).
Caveat: n≈190 → underpowered; candidate loci, not definitive. Files `gwas_adj_*`.

## 2. Functionally-stratified prediction (MultiBLUP-style) — TESTED, no gain
Partitioned 1.49M SNPs into **genic (1.04M)** vs **intergenic (0.45M)** relationship matrices and
fit 2-component REML (Speed & Balding LDAK; Edwards 2016). Result: **mean Δr = −0.007** (range
−0.03 to +0.02). Functional stratification does **not** improve DGRP prediction — joining deep
learning and non-linear kernels as modern techniques that do not break the unrelated-panel ceiling.
The only lever that helps is the RNA-seq modality. (`stratified_prediction.csv`)

## 3. Transcriptome-wide gene association (QTT, Morgante-style) — POSITIVE
For each predictable trait, correlated every gene's RNA-seq expression with the phenotype (BH-FDR).
- **starvation♂: 23 genes FDR<0.05** — dominated by **lipid metabolism: Lip4 (lipase), Jabba,
  Cyp18a1.**
- **starvation♀: 10 genes** incl. **Lsd-1 (Lipid storage droplet-1).**
- climbing/startle/sleep: 0 FDR-significant QTTs (consistent with their weaker expression signal).
**Mechanistic coherence:** lipid-storage/turnover genes track starvation resistance — the textbook
biology of starvation survival — recovered directly from expression. This is the mechanism behind the
expression modality's edge for starvation. Files `qtt_top_*`.

## 4. Web interface (NAR Web Server track)
`app.py` (Streamlit) wraps `DroPhenoPredictor` with: per-line phenotype-profile prediction +
calibrated intervals + modality + predictability; benchmark/accuracy view; per-trait GWAS + QTT
candidate genes; honest methods/limitations page. Verified to boot and serve (health 200, no errors).
Run: `./.venv/Scripts/streamlit run app.py`.

## Net scientific message
Two modern techniques tested: functional stratification (**no gain**) and QTT (**mechanistic genes**).
Together with the DL/kernel negative controls, this establishes that — at the DGRP's n≈190 unrelated
scale — predictive accuracy is governed by *data modality*, not algorithm sophistication, while the
RNA-seq transcriptome provides both the accuracy gain (starvation, climbing) and the mechanism (lipid
genes). That is the honest, defensible innovation.
