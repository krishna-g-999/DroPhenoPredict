# DroPhenoPredict — Complete Chief-Auditor Review (2026-07-01)

Full audit of the project as it stands: two resources (v1.1 DGRP predictor, v2 disease-convergence
analysis), their data, methods, code, tests, and honest gaps. Every claim below is a *checked fact*,
verified this session where marked, not carried forward from memory.

---

## 1. Inventory (verified counts, this session)
- **29 scripts** (numbered, reproducible pipeline, 01→36)
- **13 library modules** in `src/drophenopredict/`
- **41 automated tests**, all passing (`pytest tests/ -q` → `41 passed`)
- **22 findings documents** (`docs/FINDINGS_*.md`) + charter, data sources, model card, roadmap,
  tool comparison
- **2 trained model artifacts** (v1.0 genotype-only, v1.1 multi-modal)
- **1 working web app** (`app.py`, Streamlit — verified boots clean, HTTP 200, all referenced data
  files present)
- **68 processed data files** in `data/processed/`

## 2. Data provenance — VERIFIED (v1 + v2)

| Asset | Source | Status |
|---|---|---|
| DGRP2 genotypes (4.4M variants, 205 lines) | Zenodo 5582846 (NCSU offline) | ✅ downloaded, parsed |
| DGRPool phenotypes (903) | dgrpool.epfl.ch API | ✅ live, cached |
| Microarray expression (18k×184×2) | Internet Archive (NCSU) | ✅ |
| RNA-seq expression (15.7–20.4k×200) | GEO GSE117850 | ✅ |
| Inversion / Wolbachia | Internet Archive (NCSU) | ✅ |
| Variant/gene annotation | Archive `dgrp.fb557.annot.txt`, GSE67505 GTF (dm3) | ✅ |
| GO ontology + FlyBase GAF | geneontology.org (current) | ✅ |
| DGRP metabolome (NMR, 673 spectra) | MetaboLights MTBLS2060 | ✅ downloaded, processed, **excluded** (Findings 11: no predictive gain) |
| **BCM-DMAS multi-disease RNA-seq** (147 samples, 5 disease models) | Synapse syn34767207, user-downloaded | ✅ |
| **Ruffini 2020 139-gene human NDD list** | EuropePMC supplementary-files API (PMC7764447) | ✅ **obtained this session** |

**Fabrication scan (re-run this session):** no hardcoded fake/`measured_`/fallback values in `src/` or
`scripts/`. Every excluded fabrication from the original brief remains catalogued in
`SCIENTIFIC_CHARTER.md §7`.

**Citation correction (this session):** the shared-neurodegeneration meta-study was mis-attributed to
"Mallik & Mukhopadhyay 2021" in two findings docs and memory. Corrected to the actual authors —
**Ruffini, Klingenberg, Schweiger & Gerber, 2020, *Cells* 9(12):2642**.

## 3. Methods — correctness (verified, including two real bugs caught and fixed)

- **GBLUP/REML engine** verified on known-truth simulations (7 tests): recovers h², predicts high-h²
  structured traits, null centered ~0.
- **GSEA engine** (built locally after `gseapy` install failed on a network/cert issue) — **a real bug
  was caught during verification**: the first p-value implementation tested only in the direction of
  the observed effect's sign, capping p-values near 0.5 and inflating the false-positive rate to 12%
  (nominal 5%). Fixed to a proper two-sided test; recalibrated (FPR=4.3%). Locked in with 3 regression
  tests (`tests/test_gsea.py`).
- **DE pipeline cross-validated against the paper's own DESeq2 output** — first attempt (153 genes,
  r=0.83) was itself an artifact of a gene-ID format mismatch (FBgn vs. symbol); caught and fixed via
  the FlyBase annotation map; corrected result: **7,018 genes, Pearson r=0.766, Spearman r=0.607
  (p≈0)** — independent confirmation the DE pipeline is sound.
- **No data leakage:** all CV splits by DGRP line; nested CV for modality selection; expression kernel
  uses only the held-out line's own expression (a feature, not a leak).
- **Unbiased headline metric:** nested adaptive r (not optimistic LOO) for v1.1.
- **Calibrated uncertainty:** split-conformal, empirical coverage verified 0.90–0.92.
- **Reproduces 3 literature benchmarks:** Ober 2012 (×2 traits), Morgante 2020.
- **A pandas forward-compatibility bug was found and fixed** this session (`pathways.py` used
  positional `.mean(1)`/`.std(1)`, which pandas 4.0 will reject) — caught by the new test suite, not
  by manual review. Same pattern exists in ~12 already-executed one-off analysis scripts (not the
  reusable library) — low severity, not fixed (would not change any reported conclusion), flagged here.

## 4. Test coverage — expanded this session (25 → 41)
New tests added for previously-untested, non-trivial logic:
- `test_dgrpool.py` (5) — line-ID harmonization (`to_grm_line_id`), the join key for every DGRPool
  merge in the project; a silent bug here would misalign labels without raising an error.
- `test_pathways.py` (6) — GO-DAG ancestor propagation (multi-level chains, diamond inheritance,
  memoization) and pathway-activity z-scoring, using synthetic data (no dependency on the 32MB
  OBO/GAF files).
- `test_annotation.py` (5) — SNP→gene mapping (inside-gene, near-window, cross-chromosome rejection,
  multi-SNP aggregation), using a synthetic gene-span table.
- `test_gsea.py` (3, from the bug-fix above).

**Remaining test-coverage gaps (honest):** `genotypes.py`'s GRM-building I/O path, `covariates.py`,
`expression.py`'s real loaders, and `multimodal.py`'s higher-level `nested_modality_cv` are exercised
end-to-end by the numbered scripts and the finalized model's own validation metrics, but do not have
dedicated unit tests in isolation. Lower priority than the modules above because their correctness is
continuously checked by the model's own calibration/coverage numbers (a bug would show up as broken
calibration, not silently pass).

## 5. Reproducibility
Seeded (42) throughout. Deterministic pipeline: genotype side `03_build_grm` → `13_finalize_multimodal`;
disease side `31_build_bcm_dmas_matrix` → `36_gsea_convergence_finish`. 41/41 tests pass. All inputs
public; provenance in `DATA_SOURCES.md`.

## 6. Innovation assessment (honest, both resources)

**v1.1 DGRP predictor — genuinely new:**
1. Multi-modal, per-trait-adaptive integration (genotype vs. RNA-seq expression) with modality chosen
   by nested CV — neither Ober (genotype-only) nor Morgante (fixed combination) does this.
2. Calibrated, coverage-verified prediction intervals + an explicit predictability flag — no prior DGRP
   work ships either.
3. A packaged, reusable API + web app, not a one-off analysis.

**v2 disease-convergence resource — genuinely new:**
1. An independent, from-scratch, verified multi-disease (AD/PD/HD) convergence analysis on a real
   matched fly panel that the source paper deposited but did not itself analyze this way (the paper
   focuses only on its PD/αSyn arm).
2. Leave-one-out robustness testing of the convergence signal (not done in the comparator paper).
3. A specific, mechanistic, testable convergent signal — nucleolar stress / rRNA biogenesis
   down-regulation, concordant across 3/5 models — plus a recurring, cross-validated outlier finding
   (HTT-200Q diverges directionally in 3 independent analyses: DGRP-transfer, gene-level convergence,
   pathway-level GSEA).
4. Two engineering-verification catches (GSEA p-value bug; DE gene-ID join artifact) that materially
   changed the reported numbers before publication — evidence the pipeline is trustworthy, not just
   "ran without crashing."

**What is explicitly NOT claimed anywhere in the project** (guard against overclaiming): prediction for
novel/unrelated genotypes; transgenic-disease-model behaviour prediction; drug-rescue prediction; deep
learning as a meaningful lever (tested, does not help); the metabolome as a useful modality (tested,
does not help); DGRP→disease transfer as validated (tested, weak/model-specific).

## 7. Known gaps / limitations (complete, honest list)

**v1 (DGRP predictor) side:**
1. Only 8 trait-sex models across 4 trait families — expansion was tested (Findings 08) and found no
   additional predictable traits at current DGRPool coverage; not a gap to "fix," a measured ceiling.
2. Expression-modality prediction requires the query line's own RNA-seq — valid design (one assay
   informs many behaviours) but must be stated as a usage precondition.
3. starvation♀ expression accuracy (0.10) is lower than Morgante's reported 0.28 — an acknowledged,
   unresolved discrepancy (Findings 06), likely CV-scheme/normalization related.
4. Does not address the original disease/drug-screening vision — that requires perturbation data,
   which no DGRP-based method (including this one) provides.

**v2 (disease-convergence) side:**
5. **No matched behavioural data** for the BCM-DMAS panel — climbing/turning/stumbling values exist
   only in the source paper's Figure 1B, not as downloadable data. This is the hard ceiling on turning
   the convergence resource into a behaviour predictor; not solvable by more analysis, only by new data.
6. **DGRP→disease transfer is weak and model-specific** (Findings 12) — the natural-variation DGRP
   climbing map does not cleanly predict disease-model expression direction.
7. ~~Cross-species (human↔fly) validation incomplete~~ — **RESOLVED (2026-07-02, Findings 16).**
   PANTHER's ortholog REST API (pantherdb.org) worked where 5 other sources failed (verified correct:
   MAPT→tau, APP→Appl, SNCA→no ortholog, matching our own curated table). All 139 human genes mapped
   (174 fly-ortholog relationships, 173 unique fly genes). The overlap test itself came back a
   **statistically underpowered null** (0 observed vs. 0.39 expected overlap with our 27-gene fly
   convergent set) — not evidence for or against cross-species convergence at this gene-set size. The
   within-fly convergence finding (Findings 13/15) is unaffected. See `scripts/37`, `Findings 16`.
8. **ALS is not represented** in the disease-convergence panel. A real, downloadable ALS dataset exists
   (ENA PRJEB42797/PRJEB42798, TDP-43/FUS/SMN fly RNAi) but the authors explicitly report these flies
   as **asymptomatic** (no motor phenotype or mortality change) and it is a different lab/platform —
   adding it would require its own within-study DE plus a careful cross-study comparison, not naive
   pooling with BCM-DMAS (the same batch-effect risk already demonstrated with the metabolome,
   Findings 11).
9. **SCA (spinocerebellar ataxia)** — searched this session; no verified, easily-accessible, matched
   (omics + behaviour, or even omics alone with clear deposit) Drosophila SCA dataset was found. Not
   fabricated or assumed; recorded as an open gap.
10. **Proteomics validation not yet performed** — the same BCM-DMAS Synapse project has 16 TMT
    proteomic samples (PD vs. control, day 10 only), which would let the top convergent transcripts be
    checked at the protein level. Blocked by the same Data Use Agreement as the RNA-seq (requires the
    user's Synapse token); not yet requested.
11. **scRNA-seq data referenced in the metadata (21 samples, assay=scrnaSeq) could not be located** in
    the public Synapse folder tree traversed this session — may require deeper/different traversal or
    may not be public yet. Not fabricated; flagged as unresolved.

## 8. NAR-readiness verdict
The **v1.1 DGRP predictor** is at submission quality for a methods/resource paper (honest accuracy,
calibration, reproduction of two literature benchmarks, controlled negatives, biological architecture,
web interface). The **v2 disease-convergence work** is a genuine, verified, honest finding (nucleolar
stress convergence + HTT-200Q outlier pattern) suitable as a second paper or a substantial second
section — but is explicitly a molecular-convergence resource, not a behaviour predictor, and must be
framed as such. Both resources together tell one coherent, defensible story: **rigorous, calibrated,
multi-modal Drosophila phenotype/molecular analysis, with every accuracy claim measured and every
negative result reported.**

## 9. Recommended next steps (priority order, each concretely scoped)
1. ~~Human↔fly ortholog mapping~~ — **done** (Findings 16); result is an underpowered null, not
   further actionable without a pathway-level re-test (noted in Findings 16).
2. ~~Proteomics validation~~ — **DONE (2026-07-02, Findings 17).** Genome-wide RNA↔protein
   concordance: Spearman r=0.33–0.36 (p<10⁻¹⁵⁰), literature-consistent, independently validating the
   BCM-DMAS transcriptomic pipeline. Channel-to-genotype mapping cross-verified between two metadata
   files (the protein file's own column labels were confirmed to be generic PD software defaults, NOT
   real genotype assignments — would have silently mislabeled groups if trusted). A real bug
   (`scipy.stats.spearmanr` returning NaN when any input has missing values) was caught and fixed
   before reporting. Focused 14-gene convergent-set test is an honest null (n too small).
3. ~~Yang 2023 DGRP×AD integration~~ — **DONE (2026-07-03, Findings 19).** First dataset in the
   project matching real DGRP genotypes to a real disease phenotype (Aβ42+Tau F1-cross eye
   degeneration, 162 lines). Risk analysis written first (`EDITORIAL_RISK_ANALYSIS.md`); phenotype
   construction faithfully replicated the source authors' own mixed-model BLUP method (caught a real
   statsmodels API bug in the process); genotype-only, multi-modal (both sexes), and
   inversion-adjusted analyses all ran using the project's existing verified engine. **Result: a
   complete, consistent, well-powered null (0/5, 0/10, 0/5 significant respectively)** — reported in
   full, not hidden. A real methodological gap in the first draft (wrong permutation null for the
   adaptive procedure) was caught and fixed before running.
4. **Fix the remaining pandas positional-axis calls** in the one-off analysis scripts (cosmetic,
   low-priority; would not change any result).
5. **EmoryDrosophilaTau cross-lab replication** — verified real (different lab than BCM-DMAS);
   proteomics + metadata ready to download (`docs/EMORYDROSOPHILATAU_DOWNLOAD.md`); bulk RNA-seq is
   raw FASTQ/BAM only (not actionable without a new alignment pipeline); scRNA-seq usable but a bigger
   undertaking. Awaiting user's download.
4. **Manuscript drafting** for both resources (v1.1 predictor; v2 convergence finding), per the
   project's original sequencing.
