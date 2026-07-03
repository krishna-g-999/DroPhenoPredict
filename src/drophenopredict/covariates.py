"""DGRP fixed-effect covariates: inversion karyotype + Wolbachia status.

Known major cofactors of DGRP trait variation (Huang et al. 2014). Standard
practice is to account for them. Note: inversions are large LD blocks already
partly captured by the genome-wide GRM, so their *prediction* gain is expected
to be small — we include them per the charter and MEASURE the effect honestly.

Source (Internet Archive snapshots of the offline NCSU DGRP2 server):
  inversion.xlsx  (2020-02-23 snapshot)
  wolbachia.xlsx  (2019-12-30 snapshot)
Line ids harmonized to `line_NN` to match the GRM.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .dgrpool import to_grm_line_id

logger = logging.getLogger(__name__)

COV_DIR = Path("data/raw/dgrp_covariates")
# Common polymorphic, trait-relevant inversions in the DGRP.
CANDIDATE_INVERSIONS = ["In(2L)t", "In(2R)NS", "In(3R)P", "In(3R)Mo", "In(3R)K", "In(3R)C"]


def _read_headered_xlsx(path: Path) -> pd.DataFrame:
    """These sheets carry their real header on the row whose first cell is
    'DGRP Line' (sometimes preceded by a blank row). Locate it robustly."""
    raw = pd.read_excel(path, header=None)
    hdr = 0
    for i in range(min(6, len(raw))):
        if "dgrp line" in str(raw.iloc[i, 0]).lower():
            hdr = i
            break
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(c).strip() for c in raw.iloc[hdr]]
    return df.reset_index(drop=True)


def _code_inversion(v) -> float:
    s = str(v).upper()
    if "INV" in s and "ST" in s:
        return 0.5            # heterozygous karyotype
    if "INV" in s:
        return 1.0
    if "ST" in s:
        return 0.0
    return np.nan


def load_covariates(*, min_minor_freq: float = 0.05) -> pd.DataFrame:
    """Return covariate design matrix indexed by `line_NN`.

    Columns: polymorphic inversions (0/0.5/1) with minor frequency >= threshold,
    plus `Wolbachia` (1=infected, 0=uninfected). Missing values mean-imputed.
    """
    inv = _read_headered_xlsx(COV_DIR / "inversion.xlsx")
    inv_lines = [to_grm_line_id(str(i)) for i in inv.iloc[:, 0]]

    cols = {}
    for c in CANDIDATE_INVERSIONS:
        if c not in inv.columns:
            continue
        coded = pd.Series(inv[c].map(_code_inversion).to_numpy(), index=inv_lines)
        coded = coded.fillna(coded.mean())
        mf = min(coded.mean(), 1 - coded.mean())   # minor "allele" freq proxy
        if mf >= min_minor_freq:
            cols[c] = coded
        else:
            logger.info("drop monomorphic inversion %s (mf=%.3f)", c, mf)

    wol = _read_headered_xlsx(COV_DIR / "wolbachia.xlsx")
    wol_lines = [to_grm_line_id(str(i)) for i in wol.iloc[:, 0]]
    wstatus = wol.iloc[:, 1].astype(str).str.lower().str[0].map({"y": 1.0, "n": 0.0})
    wmap = pd.Series(wstatus.to_numpy(), index=wol_lines)
    cols["Wolbachia"] = wmap.fillna(wmap.mean())

    cov = pd.DataFrame(cols)
    cov = cov[~cov.index.duplicated()]
    logger.info("Covariates: %d lines x %d cols (%s)", cov.shape[0], cov.shape[1],
                ", ".join(cov.columns))
    return cov


def design_matrix(cov: pd.DataFrame, lines: list[str]) -> np.ndarray:
    """Intercept + covariates for the given lines (missing rows -> column means)."""
    sub = cov.reindex(lines)
    sub = sub.fillna(sub.mean())
    return np.column_stack([np.ones(len(lines)), sub.to_numpy(dtype=float)])
