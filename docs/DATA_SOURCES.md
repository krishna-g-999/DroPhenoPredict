# DroPhenoPredict — Data Sources Registry

Every external data asset used by this project is registered here with: what it provides, its access method, and a **verification status**. Nothing is used in modeling until its status is `VERIFIED` (accessed and inspected by us).

Status legend: `VERIFIED` = we have fetched and inspected it · `CITED` = exists per literature, not yet fetched · `TO-VERIFY` = claimed but unconfirmed · `EXCLUDED` = will not be used (with reason).

---

## Primary training assets

### 1. DGRP genotypes (Freeze 2.0)
- **Provides:** ~4.85M SNPs + ~1.3M non-SNP variants across 205 inbred lines; inversion & Wolbachia status per line.
- **Source:** DGRP2 portal — http://dgrp2.gnets.ncsu.edu/
- **Access:** Full VCF (>1 GB). We do NOT load raw SNPs into the model (overfit). We derive pathway-collapsed variant burden + use inversion/Wolbachia as covariates.
- **Citation:** Mackay et al. 2012 *Nature* 482:173; Huang et al. 2014 *Genome Res* (PMC4079974).
- **Status:** `CITED` — portal confirmed live; VCF not yet downloaded.

### 2. DGRP harmonized phenotypes (DGRPool) — TRAINING LABELS
- **Provides:** harmonized phenotype tables across 139 studies / 838 phenotypes / 22 categories, per DGRP line.
- **Source:** https://dgrpool.epfl.ch/ · API + reproducible fetch code: https://github.com/DeplanckeLab/dgrpool
- **Access:** documented API; bulk fetch script in the GitHub repo.
- **Citation:** Gardeux et al., *eLife* 88981 (2023/2024).
- **Status:** `VERIFIED` (2026-06-28) — JSON API live and inspected. Confirmed endpoints:
  - `GET /studies.json` → 142 study metadata records.
  - `GET /phenotypes.json` → 903 phenotype definitions (595 numeric-continuous, non-obsolete).
  - `GET /phenotypes/{id}.json` → per-phenotype `original_data`: `{line: {sex: [values]}}` (the actual ground-truth labels). E.g. id 1537 `mn_Lifespan` (Days), id 2847 `Climbing_Index_CTRL`.
  - Client implemented in `src/drophenopredict/dgrpool.py`; coverage measured by `scripts/01_measure_coverage.py`.
- **⚠ Critical open question:** how many lines have *all* candidate target traits measured? The original brief assumed "all 205 for all traits" — that is false (phenotypes come from different studies on different line subsets). The usable n is an empirical unknown to be measured here first.

### 3. DGRP per-line transcriptome — TRANSCRIPTOMIC FEATURES
- **Provides:** line-level microarray expression, **18,140 genes × ~184 lines, both sexes**, log2,
  2 replicates/line. Line ids `line_NN` match the PLINK/GRM ids directly.
- **Source:** Huang et al. 2015 *PNAS* 112:E6010. Original host `dgrp2.gnets.ncsu.edu/data/website/`
  is **offline**; recovered from the **Internet Archive** snapshot 2024-07-30:
  `http://web.archive.org/web/20240730id_/http://dgrp2.gnets.ncsu.edu/data/website/dgrp.array.exp.{female,male}.txt`
- **Access:** Wayback `id_` raw download → `data/raw/dgrp_expression/`. Parser + TRM builder in
  `src/drophenopredict/expression.py`.
- **Status:** `VERIFIED` (2026-06-29) — downloaded (113 MB/sex) and parsed; structure confirmed
  (space-delimited, header `gene line_NN:rep`, FBgn rows, log2 values ~3–8). Used to build the
  transcriptomic relationship matrix (TRM) for transcriptome-based GBLUP (Option 1).
- *Channels checked and rejected (recorded so they aren't re-tried):* quantgenet.msu.edu expression
  files 404; GEO GSE67505 FPKM is **pooled** (3 cols gene/female/male), not per-line; GEO GSE60314
  is only 16 lines; ArrayExpress E-MTAB-3216 is raw CEL only.

---

## Supporting / reference assets

### 3b. DGRP per-line RNA-seq expression — STRONGER TRANSCRIPTOMIC FEATURES
- **Provides:** RNA-seq log2-FPKM line means, female **15,732 genes × 200 lines**, male
  **20,375 × 200**. Columns `NN_F_mean` → harmonized to `line_NN`.
- **Source:** GEO **GSE117850** (Everett/Huang 2020, "Gene Expression Networks in the DGRP") —
  the data used by Morgante et al. 2020. Files Table_3 (female) / Table_4 (male) Gene Line Means.
- **Access:** GEO FTP `…/GSE117nnn/GSE117850/suppl/`. Loader `expression.load_rnaseq`.
- **Status:** `VERIFIED` (2026-06-29) — downloaded, parsed, used in `scripts/12_rnaseq_prediction.py`.
- **Why it matters:** RNA-seq expression predicts several traits far better than genotype (starvation♂
  r=0.38, climbing♂ 0.29), reproducing Morgante 2020, where the **microarray** data (#3) gave ~0
  (Findings 06). This is the modality that raises the natural-variation prediction ceiling.

### 4. FlyBase — gene annotation, orthology, alleles
- **Provides:** gene IDs, function, human orthologs (DIOPT), allele/phenotype annotations.
- **Source:** https://flybase.org · API & bulk downloads.
- **Status:** `CITED`. Note: the API endpoint hardcoded in the original spec must be verified against current FlyBase API docs before use (endpoints change between releases).

### 5. GEO — disease-model & treatment transcriptomes (EXTENSION layer only)
- **Provides:** expression from transgenic/disease/drug studies, for the extrapolation layer (§6 of Charter).
- **Source:** https://www.ncbi.nlm.nih.gov/geo/ via NCBI eutils / `GEOparse`.
- **Status:** `TO-VERIFY` per accession. Every GSE referenced (e.g. those in the original spec) must be individually confirmed to exist and to contain what is claimed before use. Accessions in the original brief are **unverified** until checked.

### 6. ChEMBL / DrugBank — drug-target features (EXTENSION layer only)
- **Provides:** measured drug–target affinities for compound feature encoding.
- **Source:** ChEMBL API (https://www.ebi.ac.uk/chembl/) / DrugBank.
- **Status:** `CITED`. Replaces the invented `DRUG_PROFILES` dict. Every affinity stored with its ChEMBL/source ID.

### 7. DGRP drug / toxin response phenotypes (collected 2026-11-07, scripts 45–48)
- **Provides:** per-line responses of DGRP2 lines to 21 compounds: 69 trait–sex combinations,
  99–205 genotyped lines each, 12,107 line values. Harmonised table:
  `data/processed/chemical_traits_long.csv` + `chemical_traits_catalog.csv` (source, role, n per trait).
- **Search:** DGRPool live catalog (keyword search of 142 studies) + Europe PMC (609 DGRP papers →
  69 chemical-exposure genetics papers) + DOI lookups in Dryad, Zenodo and figshare.
- **Status:** `VERIFIED` (2026-11-07). Every file downloaded by script with URL + SHA-256 manifest
  (`data/raw/chemical_supp/download_manifest.csv`, `data/raw/chemical_repo/download_manifest.csv`).
  Cross-check: deltamethrin survival computed from the raw Duneau 2018 table agrees with the DGRPool
  copy (r = 0.93). The DGRPool climbing trait used by v1.1 (phenotype 1543) is the **untreated control
  arm** of Jordan 2012's menadione study; its treated arm is now included.

| Compound | Source (data location) | Lines in panel |
|---|---|---|
| ethanol (adult sensitivity/tolerance; developmental) | Morozova 2015, 2018 (DGRPool) | 198–205 |
| lead (development, viability, activity; + controls) | Zhou 2016 (DGRPool) | 166–200 |
| menadione sodium bisulfite (climbing, startle; + controls) | Jordan 2012 (DGRPool) | 185 |
| menadione, paraquat (resistance) | Weber 2012 (DGRPool) | 167 |
| paraquat low dose (climbing; + control) | Lovejoy 2021 (DGRPool) | 99 |
| methylmercury, caffeine (developmental survival) | Montgomery 2014 (DGRPool) | 171 |
| caffeine, boric acid (lifespan, females) | Najarro 2015, 2017 (DGRPool) | 158–163 |
| methylphenidate (food intake, activity; + sucrose control, males) | Rohde 2019 (DGRPool; raw on figshare) | 172 |
| azinphos-methyl (LD50) | Battlay 2016 (DGRPool) | 175 |
| ryanoid insecticide (LC50) | Green 2018 (DGRPool) | 152 |
| malathion (3/6/12/24 h mortality), permethrin (3 h knockdown, 24 h mortality) | Battlay 2018 (figshare 10.25387/g3.6954443) | 169–170 |
| parathion (48 h survival, from per-fly records), deltamethrin (48 h survival) | Duneau 2018 (figshare 10.25387/g3.7047734) | 194–195 |
| imidacloprid (larval wiggle index, 25 and 100 ppm) | Denecke 2017 (Sci Rep supplement) | 171 |
| α-amanitin (eclosion at 0.2 and 2 µg/g) | Mitchell 2017 (PLoS One supplement) | 145 |
| rapamycin (pupation time; + vehicle control) | Harrison 2024 (Aging Cell supplement) | 140 |
| octanoic acid / noni toxin (survival) | Marconcini 2026 (G3 supplement) | 184 |

- **Found but not usable as per-line DGRP2 data (recorded, not used):** cocaine/methamphetamine
  (Highfill 2019, 46 lines), BTEX (Adebambo 2020, 79 lines), DDT (Schmidt 2017, 38 lines) — below the
  80-line threshold; 4-methylimidazole (Collins 2024) — DGRP3 panel, not covered by our genotypes;
  cadmium/lead extreme-QTL (Zhou 2017) — pooled, no per-line values; copper (Green 2022) — European
  wild strains; nicotine (Highfill 2017), venom (Ward 2026) — DSPR panel; proteasome inhibitor
  (Tu'ifua 2025), rapamycin (Peddibhotla 2026) — 1–3 genotypes; methotrexate (Kislukhin 2013) —
  Dryad file requires login (HTTP 401); anaesthetics (Olufs 2018, Scharenbrock 2020,
  Johnson-Schlitz 2025), oxidative transcriptome (Ramnarine 2022) — figures/PDF only; peroxide
  survival (Harrison 2020) — per-fly records in DGRPool, not yet summarised; rapamycin fecundity
  (Rohde 2021) — 6 lines. Amphetamine transcriptome (Karam 2022, lab strain) kept for the drug
  signature work, not as a DGRP trait.

---

## EXCLUDED (will not be used as model features)

| Asset | Reason |
|---|---|
| Hardcoded "fallback" DGRP phenotype values (orig. spec L136–150) | Fabricated. `EXCLUDED`. |
| Hardcoded `measured_climbing_pi` external values (orig. spec L664–699) | Fabricated. `EXCLUDED` until extracted from source papers with provenance. |
| Hardcoded `DRUG_PROFILES` scores (orig. spec L303–330) | Fabricated. `EXCLUDED`; replaced by ChEMBL query. |
| User's human/mouse ALS proteomics (SIRT3/HDAC11) as a per-fly-line *feature* | Category error — different organism/context. Valid as a **mechanistic prior** (pathway selection, §6) but not a DGRP feature column. `EXCLUDED` as a feature. |

---

## Verification workflow

When promoting any asset from `CITED`/`TO-VERIFY` to `VERIFIED`:
1. Fetch programmatically (reproducible script under `scripts/`).
2. Print shape, columns, dtypes, first 5 rows, and per-column non-null counts.
3. Record n_lines actually present (not assumed).
4. Update this file's status + the date.
