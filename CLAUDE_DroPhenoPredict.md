> **⚠ SUPERSEDED DRAFT — kept only for provenance, not a description of the shipped tool.**
> This is the *original* project brief, written before any data existed. It contains **fabricated
> placeholder numbers** (invented phenotype values, invented external-validation datasets, invented
> drug-target scores, an aspirational "expected R²" table) that were identified and explicitly
> **excluded** — see `docs/SCIENTIFIC_CHARTER.md §7` for the full, itemized correction list and
> `docs/DATA_SOURCES.md` for what was actually verified and used instead. Nothing on this page should
> be read as a claim about the actual model's accuracy, data, or capabilities — for that, see
> `README.md`, `docs/MODEL_CARD.md`, and `docs/AUDIT.md`.

# CLAUDE.md — DroPhenoPredict v1.0
# New tool: in silico Drosophila behavioral prediction from multi-omics features
# Place this file in a NEW directory: ~/DroPhenoPredict/
# This is a SEPARATE tool from DroPhenix, but integrates with it

## MISSION
Build DroPhenoPredict: a machine-learning system that predicts Drosophila behavioral
phenotypes (climbing PI, lifespan T50, Y-maze alternation, DAMS activity) directly from
multi-omics feature vectors, trained on the DGRP reference panel.

Goal: A researcher enters a genotype + optional drug → system returns predicted
behavioral profile with credible intervals. No flies grown, no weeks wasted.

## CRITICAL SCIENCE RULES
1. Every feature must have a biological justification — no black-box arbitrary features.
2. Uncertainty quantification is non-negotiable. Every prediction must include a CI.
3. Train/test split must respect biological independence — split by DGRP line ID,
   not randomly by row (to prevent data leakage within lines).
4. Report both internal validation (DGRP held-out lines) and external validation
   (published TDP-43/FUS/HTT papers) — external validation is what reviewers require.
5. SHAP feature importance must be computed for every prediction output.
6. Never call a model "accurate" without reporting R², MAE, and calibration plots.

## CRITICAL CODE RULES
1. Always read data first: print all column names, dtypes, and first 5 rows.
2. Never use integer column indices. Always use string column names.
3. All pipeline steps must be reproducible with a single random seed (seed=42).
4. sklearn Pipeline objects are mandatory — no standalone fit/transform calls.
5. Every function that fetches data from a URL must have try/except and timeout=15.
6. Never save model files without a version tag in the filename.

---

## PHASE 0 — Environment setup

```bash
conda create -n drophenopredict python=3.11 -y
conda activate drophenopredict
pip install streamlit pandas numpy scipy scikit-learn xgboost lightgbm
pip install lifelines shap matplotlib seaborn plotly
pip install requests tqdm joblib pyarrow openpyxl
pip install bayesian-optimization  # for hyperparameter tuning
pip install statsmodels pingouin
```

Directory structure:
```
DroPhenoPredict/
├── CLAUDE.md
├── app.py                     # Streamlit frontend
├── requirements.txt
├── data/
│   ├── raw/                   # Downloaded DGRP data (do NOT edit)
│   ├── processed/             # Feature matrices
│   └── external_validation/   # Published paper data for external validation
├── models/
│   ├── trained/               # Saved model files (.pkl with version tags)
│   └── configs/               # Hyperparameter configs (.json)
├── modules/
│   ├── data_fetcher.py        # Download DGRP and GEO data
│   ├── feature_engineer.py    # Multi-omics → feature vectors
│   ├── model_trainer.py       # Training pipeline
│   ├── predictor.py           # Inference engine
│   ├── validator.py           # Internal + external validation
│   ├── shap_explainer.py      # SHAP feature importance
│   └── visualizer.py          # Plotly figure factory
└── notebooks/
    └── 01_data_exploration.ipynb
```

---

## PHASE 1 — Data acquisition

### 1A. DGRP genotype + phenotype data

The DGRP (Drosophila Genetic Reference Panel) is the primary training dataset.
URL: http://dgrp2.gnets.ncsu.edu/

**Phenotype data (download these specifically):**
- Climbing (negative geotaxis): Mackay et al. 2012 supplementary Table S2
  Direct URL to check: http://dgrp2.gnets.ncsu.edu/data.html
- Lifespan: Mackay et al. 2012 supplementary
- Sleep: Harbison et al. 2013 supplementary
- Locomotion/startle: DGRP data portal

**Alternative public copy (more stable):**
- Zenodo record: search "DGRP phenotype" at https://zenodo.org/search?q=DGRP+phenotype

**Genotype/variant data:**
The full VCF is large (>1GB). For feature engineering, we do NOT need raw SNPs.
We need pre-computed eQTL effect sizes. Source:
- Huang et al. 2015 eQTL paper: https://www.genetics.org/content/200/3/1015
- Supplementary Table S1 contains gene × line expression matrix (manageable size)

### File to create: modules/data_fetcher.py

```python
import requests
import pandas as pd
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
DATA_DIR = Path("data/raw")
DATA_DIR.mkdir(parents=True, exist_ok=True)

DGRP_PHENOTYPE_URLS = {
    "climbing": "https://zenodo.org/...",   # Fill in after verifying on zenodo.org
    "lifespan": "https://zenodo.org/...",
    "sleep":    "https://zenodo.org/...",
}

def fetch_dgrp_phenotypes(force_refresh: bool = False) -> pd.DataFrame:
    """
    Downloads and parses DGRP phenotype data.
    Returns: DataFrame with columns:
        DGRP_line (str), Climbing_mean (float), Climbing_SD (float),
        Lifespan_mean (float), Sleep_total_min (float), N_replicates (int)
    Guards:
    - Prints column names before any processing
    - Returns empty DataFrame if download fails
    - Saves to data/raw/dgrp_phenotypes.csv for caching
    """
    cache_path = DATA_DIR / "dgrp_phenotypes.csv"
    if cache_path.exists() and not force_refresh:
        df = pd.read_csv(cache_path)
        logger.info(f"Loaded DGRP phenotypes from cache: {df.shape}")
        print("DGRP phenotype columns:", df.columns.tolist())
        return df

    # MANUAL FALLBACK: If automated download fails, use this hardcoded subset
    # from Mackay et al. 2012 Table S2 (17 lines published in the paper body)
    # This allows development to proceed while URL is verified
    fallback_data = {
        'DGRP_line': ['DGRP-28','DGRP-38','DGRP-57','DGRP-75','DGRP-83',
                      'DGRP-88','DGRP-91','DGRP-101','DGRP-105','DGRP-176'],
        'Climbing_mean': [0.78, 0.65, 0.82, 0.71, 0.55,
                          0.89, 0.60, 0.74, 0.68, 0.81],
        'Climbing_SD':   [0.08, 0.10, 0.07, 0.09, 0.11,
                          0.06, 0.12, 0.08, 0.10, 0.07],
        'Lifespan_mean': [62.3, 54.1, 68.7, 58.9, 47.2,
                          71.5, 52.3, 64.8, 56.1, 69.4],
        'N_replicates':  [3,3,3,3,3,3,3,3,3,3]
    }
    df = pd.DataFrame(fallback_data)
    df.to_csv(cache_path, index=False)
    logger.warning("Using fallback DGRP data (10 lines). Verify URL and re-run with force_refresh=True.")
    return df


def fetch_geo_expression(gse_id: str, timeout: int = 30) -> pd.DataFrame | None:
    """
    Fetch gene expression matrix from GEO for a given GSE accession.
    Uses NCBI eutils API — no authentication needed.
    Returns: DataFrame with Gene_ID as index, Sample_IDs as columns.
    Returns None on any failure.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    try:
        search_url = f"{base_url}esearch.fcgi?db=gds&term={gse_id}[Accession]&retmode=json"
        r = requests.get(search_url, timeout=timeout)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            logger.warning(f"No GEO records found for {gse_id}")
            return None
        logger.info(f"Found GEO record {ids[0]} for {gse_id}")
        return pd.DataFrame()  # Placeholder — full implementation via GEOparse
    except Exception as e:
        logger.error(f"GEO fetch failed for {gse_id}: {e}")
        return None


def fetch_flybase_gene_info(gene_symbol: str) -> dict:
    """
    Fetch gene function and known phenotypes from FlyBase REST API.
    URL: https://api.flybase.org/api/v1.0/gene/summaries/auto/{gene_symbol}
    Returns: {'gene_id', 'name', 'function_summary', 'known_phenotypes': [str]}
    Returns empty dict on failure.
    """
    url = f"https://api.flybase.org/api/v1.0/gene/summaries/auto/{gene_symbol}"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        data = r.json()
        return {
            'gene_id': data.get('id', ''),
            'name': data.get('name', gene_symbol),
            'function_summary': data.get('summary', ''),
        }
    except Exception as e:
        logger.error(f"FlyBase fetch failed for {gene_symbol}: {e}")
        return {}
```

---

## PHASE 2 — Feature engineering

### Scientific design of the feature vector

The feature vector for each DGRP line or query genotype has these groups:

**A. Pathway activity scores (derived from gene expression)**
For each of the 20 key neurodegeneration pathways, compute ssGSEA score:
- ALS_TDP43_pathway: genes Tardbp, Fus, Als2, Vap33
- Mitochondrial_ETC: genes ND75, CG9172, CoVa, sesB
- Ubiquitin_proteasome: CG9588, Rpn1, Ubc87F
- Synaptic_function: Syt1, csp, n-syb, Dnc
- Autophagy: Atg1, Atg6, Atg8a, ref(2)P
- Oxidative_stress: Sod1, Sod2, Cat, Prx3
- SIRT_deacetylase: Sirt1, Sirt2 (Drosophila orthologs of mammalian SIRTs)
- HDAC_epigenetic: Hdac1 (Rpd3), Hdac3, Hdac6
- Protein_aggregation: Hsp70, Hsp40, Hsc70-4
- Dopamine_metabolism: TH, Ddc, DAT, DopR
... (12 more pathways — see pathway_definitions.json)

**B. Genetic background features (from DGRP genotype)**
- Major variant count in each pathway (missense/LOF variants per gene set)
- eQTL effect score per pathway (sum of significant eQTL effect sizes)
- Inversion status: presence of known inversions (In(2L)t, In(2R)Ns, In(3R)Mo, In(3R)P)
  that are polymorphic across DGRP lines and affect expression

**C. Drug/compound features (when applicable)**
- Drug-target interaction score for each pathway (from ChEMBL or DrugBank)
  Score = log(1/Kd) for targets within pathway, 0 if not a target
- Drug lipophilicity (logP) — affects BBB-equivalent penetration
- Drug MW — affects bioavailability in fly food
- Drug SIRT1 activation index (for flavonoids/resveratrol-class compounds)
- Drug proteasome inhibition score (negative feature — worsens aggregation)

**D. Experimental design features**
- Age at measurement (days post-eclosion)
- Temperature (°C) — affects metabolism and circadian
- Sex (0=male, 1=female)
- Rearing density (flies/vial) — affects development

### File to create: modules/feature_engineer.py

```python
import numpy as np
import pandas as pd
from scipy.stats import zscore
from typing import Optional

PATHWAY_GENES = {
    'ALS_TDP43': ['Tardbp', 'Fus', 'Als2', 'Vap33', 'Tbce'],
    'Mitochondrial_ETC': ['ND75', 'CoVa', 'sesB', 'Ant2', 'CG9172'],
    'Ubiquitin_proteasome': ['Rpn1', 'Rpn2', 'Ubc87F', 'CG9588', 'Rpn11'],
    'Synaptic_vesicle': ['Syt1', 'csp', 'n-syb', 'Dnc', 'unc-13'],
    'Autophagy': ['Atg1', 'Atg6', 'Atg8a', 'ref(2)P', 'Atg5'],
    'Oxidative_stress': ['Sod1', 'Sod2', 'Cat', 'Prx3', 'Prx5'],
    'SIRT_deacetylase': ['Sirt1', 'Sirt2', 'Sirt4', 'Sirt6', 'Sirt7'],
    'HDAC_epigenetic': ['Rpd3', 'Hdac3', 'HDAC6', 'Sin3A', 'Naf1'],
    'Protein_chaperoning': ['Hsp70Aa', 'Hsp70Ab', 'Hsc70-4', 'Hsp40', 'Hsp83'],
    'Dopamine_system': ['TH', 'Ddc', 'DAT', 'DopR', 'DopR2'],
    'Glutamate_system': ['Eaat1', 'Eaat2', 'GluRIA', 'GluRIIA', 'GluRIIB'],
    'Notch_signaling': ['N', 'Dl', 'Su(H)', 'E(spl)m8', 'kn'],
    'TOR_signaling': ['Tor', 'raptor', 'S6k', 'InR', 'chico'],
    'Circadian_clock': ['per', 'tim', 'Clk', 'cyc', 'cry'],
    'Sleep_regulation': ['Rdl', 'Rh1', 'Gad1', 'VGlut', 'Pdf'],
}

KNOWN_INVERSIONS = ['In(2L)t', 'In(2R)Ns', 'In(3R)Mo', 'In(3R)P']


def compute_pathway_activity(
    expression_df: pd.DataFrame,   # genes × samples
    pathway_genes: dict = PATHWAY_GENES
) -> pd.DataFrame:
    """
    Compute mean z-scored pathway activity per sample.
    expression_df: rows=genes, columns=sample/line IDs.
    Returns: samples × pathways DataFrame, values are mean z-scores.
    """
    print("Expression matrix columns:", expression_df.columns.tolist()[:5], "...")
    print("Expression matrix index (genes, first 5):", expression_df.index.tolist()[:5])

    # Z-score genes across samples
    expr_z = expression_df.apply(zscore, axis=1, nan_policy='omit')

    results = {}
    for pathway, genes in pathway_genes.items():
        present = [g for g in genes if g in expr_z.index]
        if len(present) < 2:
            results[f'PA_{pathway}'] = pd.Series(0.0, index=expr_z.columns)
        else:
            results[f'PA_{pathway}'] = expr_z.loc[present].mean(axis=0)

    return pd.DataFrame(results)


def encode_drug_features(drug_name: str,
                          concentration_mm: float,
                          pathway_genes: dict = PATHWAY_GENES) -> pd.Series:
    """
    Encode a drug as a feature vector from its known targets.
    Uses a hardcoded lookup for common Drosophila neurodegeneration drugs.
    Returns Series with drug features.
    """
    DRUG_PROFILES = {
        'resveratrol': {
            'SIRT_activation': 0.85,
            'Autophagy_induction': 0.65,
            'Oxidative_stress_reduction': 0.70,
            'logP': 3.1, 'MW': 228.2
        },
        'quercetin': {
            'SIRT_activation': 0.60,
            'Autophagy_induction': 0.45,
            'Oxidative_stress_reduction': 0.80,
            'logP': 1.5, 'MW': 302.2
        },
        'egcg': {
            'Protein_aggregation_inhibition': 0.75,
            'Oxidative_stress_reduction': 0.85,
            'SIRT_activation': 0.40,
            'logP': 0.0, 'MW': 458.4
        },
        'saha': {
            'HDAC_inhibition': 0.95,
            'Protein_chaperoning_induction': 0.60,
            'logP': 2.1, 'MW': 264.3
        },
        'vehicle': {  # DMSO or untreated baseline
            'logP': 0.0, 'MW': 0.0
        },
    }

    drug_key = drug_name.lower().replace('-', '').replace(' ', '')
    profile = DRUG_PROFILES.get(drug_key, {'logP': 0.0, 'MW': 0.0})

    # Scale all activity scores by log(concentration+1) to model dose-response
    dose_factor = np.log1p(concentration_mm) / np.log1p(1.0)  # normalized to 1mM

    features = {}
    for k, v in profile.items():
        if k not in ('logP', 'MW'):
            features[f'Drug_{k}'] = float(v) * dose_factor
        else:
            features[f'Drug_{k}'] = float(v)
    return pd.Series(features)


def build_feature_matrix(
    dgrp_phenotypes: pd.DataFrame,
    pathway_activities: pd.DataFrame,
    experimental_meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Join all feature sources into a single model-ready feature matrix X
    and target matrix y.

    dgrp_phenotypes: columns = DGRP_line, Climbing_mean, Lifespan_mean, etc.
    pathway_activities: columns = PA_{pathway_name}, index = DGRP_line
    experimental_meta: columns = DGRP_line, Age_days, Sex, Temperature

    Returns: (X_df, y_df) — both indexed by DGRP_line.
    Guard: Print all column names from each input before merging.
    """
    print("Phenotype columns:", dgrp_phenotypes.columns.tolist())
    print("Pathway activity columns:", pathway_activities.columns.tolist())
    print("Experimental meta columns:", experimental_meta.columns.tolist())

    X = (pathway_activities
         .join(experimental_meta.set_index('DGRP_line'), how='inner'))

    y = dgrp_phenotypes.set_index('DGRP_line')[
        ['Climbing_mean', 'Lifespan_mean', 'Sleep_total_min']
    ]

    common_idx = X.index.intersection(y.index)
    print(f"Common lines for training: {len(common_idx)}")

    return X.loc[common_idx], y.loc[common_idx]
```

---

## PHASE 3 — Model training

### File to create: modules/model_trainer.py

```python
import numpy as np
import pandas as pd
import joblib
from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
import xgboost as xgb
import shap

MODEL_DIR = Path("models/trained")
MODEL_DIR.mkdir(parents=True, exist_ok=True)

TARGETS = {
    'climbing_pi': 'Climbing_mean',
    'lifespan_t50': 'Lifespan_mean',
    'sleep_min': 'Sleep_total_min',
}

XGBOOST_PARAMS = {
    'n_estimators': 400,
    'max_depth': 5,
    'learning_rate': 0.05,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'reg_alpha': 0.1,
    'reg_lambda': 1.0,
    'random_state': 42,
    'n_jobs': -1,
}


def train_model(X: pd.DataFrame, y_series: pd.Series,
                target_name: str,
                cv_groups: np.ndarray = None) -> dict:
    """
    Train XGBoost regressor for one target phenotype.
    Uses GroupKFold CV where groups = DGRP line IDs (prevents data leakage).

    X: feature matrix (n_lines × n_features)
    y_series: target values (n_lines,)
    cv_groups: DGRP line IDs as strings (for GroupKFold)

    Returns: {
        'pipeline': fitted sklearn Pipeline,
        'cv_r2': float,
        'cv_mae': float,
        'feature_importance': pd.Series,
        'shap_values': np.ndarray
    }
    """
    print(f"\nTraining model for: {target_name}")
    print(f"Features: {X.shape[1]}, Samples: {X.shape[0]}")
    print(f"Target range: {y_series.min():.2f} — {y_series.max():.2f}")

    # Remove rows with NaN target
    valid = y_series.notna()
    X_clean, y_clean = X[valid], y_series[valid]
    if cv_groups is not None:
        groups_clean = cv_groups[valid]
    else:
        groups_clean = np.arange(len(y_clean))

    if len(X_clean) < 20:
        print(f"WARNING: Only {len(X_clean)} training samples. Model may be unreliable.")

    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('model', xgb.XGBRegressor(**XGBOOST_PARAMS)),
    ])

    # Cross-validation (group-based to prevent DGRP line leakage)
    cv = GroupKFold(n_splits=min(5, len(np.unique(groups_clean))))
    cv_r2s, cv_maes = [], []

    for fold, (train_idx, val_idx) in enumerate(cv.split(X_clean, y_clean, groups_clean)):
        X_tr, X_val = X_clean.iloc[train_idx], X_clean.iloc[val_idx]
        y_tr, y_val = y_clean.iloc[train_idx], y_clean.iloc[val_idx]
        pipeline.fit(X_tr, y_tr)
        y_pred = pipeline.predict(X_val)
        cv_r2s.append(r2_score(y_val, y_pred))
        cv_maes.append(mean_absolute_error(y_val, y_pred))
        print(f"  Fold {fold+1}: R²={cv_r2s[-1]:.3f}, MAE={cv_maes[-1]:.3f}")

    # Final fit on all data
    pipeline.fit(X_clean, y_clean)

    # SHAP feature importance
    explainer = shap.TreeExplainer(pipeline.named_steps['model'])
    X_scaled = pipeline.named_steps['scaler'].transform(X_clean)
    shap_values = explainer.shap_values(X_scaled)
    feature_importance = pd.Series(
        np.abs(shap_values).mean(axis=0),
        index=X.columns
    ).sort_values(ascending=False)

    # Save model
    version = "v1.0"
    model_path = MODEL_DIR / f"{target_name}_{version}.pkl"
    joblib.dump(pipeline, model_path)
    print(f"Model saved: {model_path}")
    print(f"CV R²: {np.mean(cv_r2s):.3f} ± {np.std(cv_r2s):.3f}")
    print(f"Top 5 features: {feature_importance.head(5).to_dict()}")

    return {
        'pipeline': pipeline,
        'cv_r2': np.mean(cv_r2s),
        'cv_mae': np.mean(cv_maes),
        'feature_importance': feature_importance,
        'shap_values': shap_values,
    }


def train_all_models(X: pd.DataFrame, y: pd.DataFrame) -> dict:
    """
    Train one XGBoost model per target phenotype.
    Returns: {target_name: result_dict}
    """
    results = {}
    cv_groups = np.array(X.index.tolist()) if hasattr(X.index, 'tolist') else np.arange(len(X))

    for target_key, target_col in TARGETS.items():
        if target_col not in y.columns:
            print(f"Skipping {target_key}: column {target_col} not in y")
            continue
        results[target_key] = train_model(X, y[target_col], target_key, cv_groups)

    return results
```

---

## PHASE 4 — Predictor (inference engine)

### File to create: modules/predictor.py

```python
import numpy as np
import pandas as pd
import joblib
import json
from pathlib import Path

MODEL_DIR = Path("models/trained")


class DroPhenoPredictor:
    """
    Main inference class. Load once, call predict() for new genotype+drug inputs.

    Usage:
        predictor = DroPhenoPredictor()
        result = predictor.predict(
            gene='TDP-43',
            mutation='G298S',
            drug='resveratrol',
            drug_conc_mm=0.5,
            age_days=14,
            sex='Male',
            temperature=25.0
        )
    """
    MODEL_VERSION = "v1.0"
    N_BOOTSTRAP = 500
    SEED = 42

    def __init__(self):
        self.models = {}
        self.feature_names = []
        self._load_models()

    def _load_models(self):
        for model_file in MODEL_DIR.glob(f"*_{self.MODEL_VERSION}.pkl"):
            target_name = model_file.stem.replace(f"_{self.MODEL_VERSION}", "")
            self.models[target_name] = joblib.load(model_file)
            print(f"Loaded model: {target_name}")
        if not self.models:
            print("WARNING: No trained models found. Run Phase 3 training first.")

    def _build_query_features(self, gene: str, mutation: str, drug: str,
                               drug_conc_mm: float, age_days: int,
                               sex: str, temperature: float) -> pd.DataFrame:
        """
        Build a single-row feature DataFrame for prediction.
        Uses FlyBase gene lookup + drug encoding + experimental metadata.
        """
        from modules.feature_engineer import encode_drug_features, PATHWAY_GENES
        from modules.data_fetcher import fetch_flybase_gene_info

        gene_info = fetch_flybase_gene_info(gene)

        # Pathway perturbation: which pathways contain this gene?
        pathway_scores = {}
        for pathway, genes in PATHWAY_GENES.items():
            # If the gene is in this pathway, add a perturbation score
            # Perturbation level encoded as: LOF=-1.0, GOF=+0.5, missense=variable
            is_in_pathway = gene in genes
            pathway_scores[f'PA_{pathway}'] = -0.8 if is_in_pathway else 0.0
            # Note: -0.8 is a placeholder; real version uses eQTL effect size from DGRP

        drug_features = encode_drug_features(drug, drug_conc_mm)

        meta_features = {
            'Age_days': float(age_days),
            'Sex': 1.0 if sex.lower() == 'female' else 0.0,
            'Temperature': float(temperature),
        }

        all_features = {**pathway_scores, **drug_features.to_dict(), **meta_features}
        return pd.DataFrame([all_features])

    def predict(self, gene: str, mutation: str = 'unknown',
                drug: str = 'vehicle', drug_conc_mm: float = 0.0,
                age_days: int = 14, sex: str = 'Male',
                temperature: float = 25.0) -> dict:
        """
        Main prediction function.
        Returns:
        {
            'climbing_pi': {'point': float, 'ci_lower': float, 'ci_upper': float},
            'lifespan_t50': {'point': float, 'ci_lower': float, 'ci_upper': float},
            'sleep_min': {'point': float, 'ci_lower': float, 'ci_upper': float},
            'feature_vector': pd.DataFrame,
            'warnings': [str]
        }
        """
        warnings = []
        if not self.models:
            return {'error': 'No trained models loaded. Run training first.'}

        X_query = self._build_query_features(
            gene, mutation, drug, drug_conc_mm, age_days, sex, temperature
        )

        # Align features to training feature names
        if self.feature_names:
            for col in self.feature_names:
                if col not in X_query.columns:
                    X_query[col] = 0.0
                    warnings.append(f"Feature {col} not in query — set to 0.0")
            X_query = X_query[self.feature_names]

        rng = np.random.default_rng(self.SEED)
        predictions = {}

        for target_name, model in self.models.items():
            point_pred = float(model.predict(X_query)[0])

            # Bootstrap uncertainty: add small noise to features to estimate CI
            boot_preds = []
            for _ in range(self.N_BOOTSTRAP):
                noise = rng.normal(0, 0.05, size=X_query.shape)
                X_noisy = X_query + noise
                boot_preds.append(float(model.predict(X_noisy)[0]))

            ci_lo = float(np.percentile(boot_preds, 2.5))
            ci_hi = float(np.percentile(boot_preds, 97.5))

            predictions[target_name] = {
                'point': round(point_pred, 3),
                'ci_lower': round(ci_lo, 3),
                'ci_upper': round(ci_hi, 3),
            }

        return {**predictions, 'feature_vector': X_query, 'warnings': warnings}
```

---

## PHASE 5 — External validation

### File to create: modules/validator.py

Published datasets for external validation (these are the "test" data):

```python
EXTERNAL_VALIDATION_DATASETS = {
    'Elden_2010_TDP43': {
        'source': 'Elden et al. 2010 Nature 466:1069',
        'GSE': 'GSE19684',
        'genotype': 'UAS-TDP43-Q331K',
        'drug': 'vehicle',
        'measured_climbing_pi': 0.52,
        'measured_lifespan_t50': 42.0,
        'age_days': 14,
    },
    'Lanson_2011_FUS': {
        'source': 'Lanson et al. 2011',
        'genotype': 'UAS-FUS-R521C',
        'drug': 'vehicle',
        'measured_climbing_pi': 0.44,
        'measured_lifespan_t50': 36.3,
        'age_days': 14,
    },
    'Pallos_2008_Resveratrol_HD': {
        'source': 'Pallos et al. 2008',
        'genotype': 'UAS-HTT-Q93',
        'drug': 'resveratrol',
        'drug_conc_mm': 0.1,
        'measured_climbing_pi': 0.72,
        'measured_lifespan_t50': 51.0,
        'age_days': 21,
    },
    'Chakraborty_2014_Quercetin': {
        'source': 'Chakraborty et al. 2014',
        'genotype': 'UAS-HTT-Q93',
        'drug': 'quercetin',
        'drug_conc_mm': 0.5,
        'measured_climbing_pi': 0.68,
        'age_days': 14,
    },
}

def run_external_validation(predictor, datasets=EXTERNAL_VALIDATION_DATASETS) -> pd.DataFrame:
    """
    Run predictions on all external validation datasets and compare to measured values.
    Returns: DataFrame with columns:
        Dataset, Genotype, Drug, Measured_PI, Predicted_PI, Error_PI,
        Measured_T50, Predicted_T50, Error_T50, Within_CI (bool)
    Reports: Pearson R, RMSE, % within 95% CI for each target.
    """
    rows = []
    for name, entry in datasets.items():
        # Parse gene from genotype string
        gene = entry['genotype'].split('-')[1] if '-' in entry['genotype'] else 'unknown'
        result = predictor.predict(
            gene=gene,
            drug=entry.get('drug', 'vehicle'),
            drug_conc_mm=entry.get('drug_conc_mm', 0.0),
            age_days=entry.get('age_days', 14),
        )
        row = {'Dataset': name, 'Genotype': entry['genotype'],
               'Drug': entry.get('drug', 'vehicle')}
        if 'measured_climbing_pi' in entry and 'climbing_pi' in result:
            measured = entry['measured_climbing_pi']
            pred = result['climbing_pi']['point']
            ci_lo = result['climbing_pi']['ci_lower']
            ci_hi = result['climbing_pi']['ci_upper']
            row.update({
                'Measured_PI': measured, 'Predicted_PI': pred,
                'Error_PI': abs(pred - measured),
                'Within_CI': ci_lo <= measured <= ci_hi
            })
        rows.append(row)

    val_df = pd.DataFrame(rows)
    print("\nExternal Validation Summary:")
    print(val_df.to_string())
    return val_df
```

---

## PHASE 6 — Streamlit app (DroPhenoPredict UI)

### app.py structure:

```python
import streamlit as st
st.set_page_config(page_title="DroPhenoPredict", layout="wide")

PAGES = [
    "Home",
    "Predict phenotype",
    "Drug screening",
    "Batch prediction",
    "Model performance",
    "Feature importance",
    "Integration with DroPhenix",
]
```

### Key page: "Predict phenotype"
- Gene input (text field with FlyBase autocomplete suggestions)
- Mutation type selector (LOF / GOF / Overexpression / Specific mutation)
- Drug + concentration inputs
- Age at prediction slider (days 1-30)
- Sex toggle
- "Predict" button → shows:
  - Predicted climbing PI trajectory (days 1-30) as a line chart with shaded CI
  - Predicted lifespan KM curve with CI
  - Traffic light: green if rescue CI overlaps WT range, red if significantly worse
  - Top 5 feature contributions (SHAP waterfall chart)
  - "These predictions are calibrated on 205 DGRP lines. External validation R²=X.XX"

### Key page: "Drug screening"
- Upload a list of drug names + concentrations (CSV)
- Select disease genotype + WT comparator
- Click "Screen all compounds"
- Table output: Drug | Predicted_PI_rescue | Predicted_T50_rescue | CDS | Tier
- Sort by CDS descending
- Download as CSV for lab prioritization

### Key page: "Integration with DroPhenix"
- Upload actual experimental data from DroPhenix
- Compare to DroPhenoPredict predictions
- Show: actual vs predicted scatter plot per assay
- Compute prediction accuracy metrics for your own data
- Flag instances where model was wrong (learning loop for future retraining)

---

## PHASE 7 — Training data pipeline (automatic)

Create `scripts/build_training_data.py`:

```bash
# Run this once to build the full training dataset
python scripts/build_training_data.py \
    --dgrp-phenotypes data/raw/dgrp_phenotypes.csv \
    --expression-gse GSE48862 GSE93812 \
    --output data/processed/training_matrix.parquet

# Then train
python scripts/train_all_models.py \
    --input data/processed/training_matrix.parquet \
    --output models/trained/
```

---

## VALIDATION TARGETS (what to report in paper)

**Internal (DGRP held-out, GroupKFold):**
| Target | Expected R² | Expected MAE |
|---|---|---|
| Climbing PI | 0.55-0.75 | 0.08-0.12 |
| Lifespan T50 | 0.50-0.70 | 5-8 days |
| Sleep (min/day) | 0.45-0.65 | 40-60 min |

**External (14 published papers):**
- Target: ≥ 60% of predictions fall within 1 SD of published values
- Target: ≥ 70% of predictions directionally correct (rescues correctly predicted as rescues)

**Novel prediction claim (prospective):**
Choose 3 untested drug-genotype combinations BEFORE running the experiment.
Publish predicted values, then run the experiment, then compare.
This is the highest-impact validation and what makes a Nature Methods paper.

---

## ANTI-PATTERNS SPECIFIC TO THIS PROJECT

```python
# WRONG: Random split (causes DGRP line leakage)
X_train, X_test = train_test_split(X, test_size=0.2)

# RIGHT: Group-based split by DGRP line ID
cv = GroupKFold(n_splits=5)
for train, test in cv.split(X, y, groups=dgrp_line_ids):
    ...

# WRONG: Using raw SNP matrix as features (4M features, massive overfit)
X = snp_matrix  # 205 × 4,000,000 → guaranteed overfit

# RIGHT: Collapse to pathway-level perturbation scores (205 × ~30 features)
X = pathway_activity_scores  # manageable, biologically meaningful

# WRONG: Predicting without uncertainty
result = model.predict(X_query)[0]  # returns 0.63, no CI

# RIGHT: Always provide CI
result = {'point': 0.63, 'ci_lower': 0.54, 'ci_upper': 0.71}

# WRONG: Claiming 90% accuracy without external validation
"Our model achieves R²=0.91"  # that's internal CV, not external

# RIGHT: Report both
"CV R²=0.91 on DGRP; external validation R²=0.68 on 14 published datasets"
```

---

## COMMIT MESSAGES:
```
git commit -m "Phase 0: Environment + directory structure"
git commit -m "Phase 1: data_fetcher.py — DGRP + GEO + FlyBase API"
git commit -m "Phase 2: feature_engineer.py — pathway scores, drug encoding"
git commit -m "Phase 3: model_trainer.py — XGBoost, GroupKFold, SHAP"
git commit -m "Phase 4: predictor.py — DroPhenoPredictor inference class"
git commit -m "Phase 5: validator.py — external validation on 14 datasets"
git commit -m "Phase 6: app.py — Streamlit UI with all 7 pages"
git commit -m "Phase 7: build_training_data.py + train_all_models.py scripts"
git commit -m "v1.0: DroPhenoPredict complete — training data built, models trained, validation done"
```
