"""DGRP per-line transcriptome -> transcriptomic relationship matrix (TRM).

Source: Huang et al. 2015 PNAS 112:E6010 microarray expression, line-level,
both sexes (~185 lines, 18,140 genes, log2). Recovered from the Internet Archive
snapshot of the (now-offline) NCSU DGRP2 server:
  http://web.archive.org/web/20240730id_/http://dgrp2.gnets.ncsu.edu/data/website/dgrp.array.exp.{sex}.txt

File format: space-delimited; header `gene line_NN:rep ...` (1-2 replicates per
line); rows are FBgn genes; values log2 expression. Line ids (`line_NN`) match
the PLINK/GRM ids directly.

The TRM is the expression analogue of the GRM: with Z = column-standardized
(z-scored across lines) gene expression, K_expr = (1/p) Z Z^T. Feeding it to the
same gblup engine yields transcriptome-based GBLUP ("T-BLUP"), directly
comparable to genomic GBLUP. No fabrication: all values parsed from the file.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

EXPR_DIR = Path("data/raw/dgrp_expression")
FILES = {"female": "dgrp.array.exp.female.txt", "male": "dgrp.array.exp.male.txt"}


def load_expression(sex: str) -> pd.DataFrame:
    """Return genes x lines log2 expression, replicates averaged per line.

    Columns are `line_NN`; index is FBgn gene id.
    """
    path = EXPR_DIR / FILES[sex]
    df = pd.read_csv(path, sep=r"\s+", index_col=0, engine="c")
    # columns like 'line_21:1','line_21:2' -> group by line, average replicates
    line_of = {c: c.split(":")[0] for c in df.columns}
    out = df.T.groupby(line_of).mean().T          # genes x lines
    logger.info("Expression %s: %d genes x %d lines (reps averaged)",
                sex, out.shape[0], out.shape[1])
    return out


RNASEQ_DIR = Path("data/raw/dgrp_rnaseq")
RNASEQ_FILES = {
    "female": "GSE117850_DGRP_GEO_Table_3_Gene_Female_Line_Means.txt.gz",
    "male": "GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz",
}


def load_rnaseq(sex: str) -> pd.DataFrame:
    """Per-line RNA-seq expression (log2 FPKM line means), genes x lines.

    Source: GSE117850 (Everett/Huang 2020, the data used by Morgante 2020).
    Columns like '21_F_mean' -> harmonized to 'line_21'.
    """
    df = pd.read_csv(RNASEQ_DIR / RNASEQ_FILES[sex], sep="\t", index_col=0, engine="c")
    df.columns = [f"line_{c.split('_')[0]}" for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    logger.info("RNA-seq %s: %d genes x %d lines", sex, df.shape[0], df.shape[1])
    return df


def build_trm(sex: str, *, min_var: float = 1e-6) -> tuple[np.ndarray, list[str], int]:
    """Transcriptomic relationship matrix from line-level expression.

    Returns (K_expr, line_ids, n_genes_used).
    """
    expr = load_expression(sex)               # genes x lines
    X = expr.to_numpy(dtype=float).T          # lines x genes
    lines = list(expr.columns)
    # keep finite, variable genes
    finite = np.isfinite(X).all(axis=0)
    X = X[:, finite]
    var = X.var(axis=0)
    keep = var > min_var
    X = X[:, keep]
    Z = (X - X.mean(0)) / X.std(0)            # standardize each gene across lines
    K = (Z @ Z.T) / Z.shape[1]
    logger.info("TRM %s: %d lines, %d genes used", sex, len(lines), Z.shape[1])
    return K, lines, int(Z.shape[1])


def save_trm(K, lines, n_genes, path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, K=K, iids=np.array(lines), n_genes=n_genes)


def load_trm(path) -> tuple[np.ndarray, list[str], dict]:
    d = np.load(path, allow_pickle=True)
    return d["K"], [str(x) for x in d["iids"]], {"n_genes": int(d["n_genes"])}
