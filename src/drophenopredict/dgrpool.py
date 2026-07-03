"""DGRPool API client — fetches harmonized DGRP phenotype data.

Source: https://dgrpool.epfl.ch  (Gardeux et al., eLife 88981)
Reproducible fetch reference: https://github.com/DeplanckeLab/dgrpool

Verified API shape (2026-06-28):
  GET /studies.json          -> list of study metadata dicts
  GET /phenotypes.json       -> list of phenotype-definition dicts
  GET /phenotypes/{id}.json  -> single phenotype incl. `original_data`:
        {"DGRP_021": {"F": [49.76], "M": [..]}, ...}   # per-line, per-sex values

Charter rules honored here:
  - Never fabricates values. On fetch failure raises/returns empty; never invents data.
  - Caches raw JSON under data/raw/dgrpool/ (download-only, never edited).
  - Every network call has a timeout and is wrapped for clear errors.
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)

BASE_URL = "https://dgrpool.epfl.ch"
USER_AGENT = "DroPhenoPredict/0.1 (research; +https://dgrpool.epfl.ch)"
TIMEOUT = 30
RETRIES = 3
RETRY_BACKOFF_S = 2.0

CACHE_DIR = Path("data/raw/dgrpool")


def _cache_path(name: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / name


def _get_json(path: str, *, cache_name: str | None = None, force: bool = False):
    """GET {BASE_URL}{path} as JSON, with on-disk cache and bounded retries.

    Returns parsed JSON (list/dict). Raises RuntimeError on definitive failure
    — it never returns fabricated data.
    """
    if cache_name is not None:
        cp = _cache_path(cache_name)
        if cp.exists() and not force:
            with cp.open("r", encoding="utf-8") as fh:
                return json.load(fh)

    url = f"{BASE_URL}{path}"
    last_err: Exception | None = None
    for attempt in range(1, RETRIES + 1):
        try:
            r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=TIMEOUT)
            r.raise_for_status()
            data = r.json()
            if cache_name is not None:
                with _cache_path(cache_name).open("w", encoding="utf-8") as fh:
                    json.dump(data, fh)
            return data
        except Exception as e:  # noqa: BLE001 - report and retry
            last_err = e
            logger.warning("GET %s failed (attempt %d/%d): %s", url, attempt, RETRIES, e)
            if attempt < RETRIES:
                time.sleep(RETRY_BACKOFF_S * attempt)
    raise RuntimeError(f"DGRPool fetch failed for {url}: {last_err}")


def fetch_studies(force: bool = False) -> pd.DataFrame:
    """All study metadata. Returns DataFrame indexed by study id."""
    data = _get_json("/studies.json", cache_name="studies.json", force=force)
    df = pd.DataFrame(data)
    logger.info("Fetched %d studies; columns=%s", len(df), df.columns.tolist())
    return df.set_index("id") if "id" in df.columns else df


def fetch_phenotype_catalog(force: bool = False) -> pd.DataFrame:
    """All phenotype definitions (NOT the values). Returns DataFrame indexed by phenotype id.

    Key columns: study_id, name, is_summary, is_continuous, is_numeric,
    obsolete, summary_type, unit.
    """
    data = _get_json("/phenotypes.json", cache_name="phenotypes.json", force=force)
    df = pd.DataFrame(data)
    logger.info("Fetched %d phenotype definitions", len(df))
    return df.set_index("id") if "id" in df.columns else df


def numeric_continuous_phenotypes(catalog: pd.DataFrame) -> pd.DataFrame:
    """Subset of the catalog usable as continuous regression targets/features."""
    mask = (
        catalog.get("is_numeric", False).astype(bool)
        & catalog.get("is_continuous", False).astype(bool)
        & ~catalog.get("obsolete", False).astype(bool)
    )
    return catalog[mask]


def fetch_phenotype_values(phenotype_id: int, force: bool = False) -> pd.DataFrame:
    """Per-line values for one phenotype, parsed to a tidy long DataFrame.

    Returns columns: [DGRP_line, sex, value]. One row per (line, sex, replicate value).
    `original_data` maps line -> {sex -> [values]}. Some studies report a single
    summary value per (line, sex); others report replicates.

    Raises RuntimeError on fetch failure. Returns an EMPTY DataFrame (never fake
    data) if the phenotype legitimately has no `original_data`.
    """
    rec = _get_json(
        f"/phenotypes/{phenotype_id}.json",
        cache_name=f"phenotype_{phenotype_id}.json",
        force=force,
    )
    original = rec.get("original_data") or {}
    rows: list[dict] = []
    for line, by_sex in original.items():
        if not isinstance(by_sex, dict):
            continue
        for sex, values in by_sex.items():
            if sex == "sex":  # guard: some payloads echo a 'sex' key
                continue
            if not isinstance(values, list):
                values = [values]
            for v in values:
                if isinstance(v, (int, float)):
                    rows.append({"DGRP_line": line, "sex": sex, "value": float(v)})
    df = pd.DataFrame(rows, columns=["DGRP_line", "sex", "value"])
    logger.info(
        "Phenotype %s (%s): %d lines, %d (line,sex,value) rows",
        phenotype_id, rec.get("name", "?"), df["DGRP_line"].nunique(), len(df),
    )
    return df


def to_grm_line_id(dgrpool_id: str) -> str:
    """Map a DGRPool line id to the PLINK/GRM id. 'DGRP_021' -> 'line_21'."""
    num = dgrpool_id.split("_")[-1].lstrip("0") or "0"
    return f"line_{num}"


def line_means(values: pd.DataFrame, sex: str | None = None,
               harmonize: bool = False) -> pd.Series:
    """Collapse tidy phenotype values to one mean per DGRP line.

    sex=None averages across available sexes; sex='F'/'M' filters first.
    Returns Series indexed by DGRP_line.
    """
    df = values if sex is None else values[values["sex"] == sex]
    s = df.groupby("DGRP_line")["value"].mean()
    if harmonize:
        s.index = [to_grm_line_id(i) for i in s.index]
        s = s[~s.index.duplicated()]
    return s
