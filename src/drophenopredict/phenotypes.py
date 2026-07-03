"""Load DGRP line-mean phenotypes (targets) from the quantgenet server files.

These files are headerless CSVs of `line_id,value` (line means), e.g.
startle.female.csv. Source: https://quantgenet.msu.edu/dgrp/downloads.html
(Mackay/Huang DGRP quantitative-genetics resource). Line IDs are `line_NNN`,
matching the PLINK .fam IDs exactly.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)

BASE_URL = "https://quantgenet.msu.edu/dgrp/data"
USER_AGENT = "DroPhenoPredict/0.1 (research)"
CACHE_DIR = Path("data/raw/dgrp_quantgenet")

# Verified-present files on the server (others on the page 404 / are offline).
PHENOTYPE_FILES = {
    "startle_female": "startle.female.csv",
    "startle_male": "startle.male.csv",
    "starvation_female": "starvation.female.csv",
    "starvation_male": "starvation.male.csv",
    "chillcoma_female": "chillcoma.female.csv",
    "chillcoma_male": "chillcoma.male.csv",
}


def fetch_phenotype(key: str, *, force: bool = False) -> pd.Series:
    """Return a Series indexed by DGRP line id (`line_NNN`) of line-mean values.

    Headerless `line,value`. Raises on network failure (never fabricates).
    """
    if key not in PHENOTYPE_FILES:
        raise KeyError(f"Unknown phenotype '{key}'. Known: {list(PHENOTYPE_FILES)}")
    fname = PHENOTYPE_FILES[key]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dst = CACHE_DIR / fname
    if not dst.exists() or force:
        url = f"{BASE_URL}/{fname}"
        # NOTE: host serves an incomplete TLS chain; content integrity is
        # confirmed by structure (line_NNN ids + numeric values) on load.
        r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=60, verify=False)
        r.raise_for_status()
        dst.write_bytes(r.content)
    df = pd.read_csv(dst, header=None, names=["line", "value"])
    df = df.dropna()
    s = df.set_index("line")["value"].astype(float)
    if not s.index.str.startswith("line_").all():
        raise ValueError(f"Unexpected line id format in {fname}: {s.index[:3].tolist()}")
    logger.info("Phenotype %s: %d lines, range %.3f–%.3f",
                key, s.size, s.min(), s.max())
    return s
