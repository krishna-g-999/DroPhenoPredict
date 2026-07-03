"""DGRP genotypes -> genomic relationship matrix (GRM) for GBLUP.

Source: DGRP Freeze 2.0 PLINK genotypes (dgrp2.bed/bim/fam), 205 inbred lines,
4,438,427 variants. Obtained from Zenodo record 5582846 (canonical NCSU host
offline 2026-06). Lines are inbred -> genotypes coded 0/2 (homozygous) with
residual hets/missing.

GRM = VanRaden (2008) method 1, computed by STREAMING over SNP blocks so the
4.4M x 205 matrix is never held in memory at once. QC per Findings: autosomal
SNPs (arms 2L/2R/3L/3R), MAF >= 0.05, call-rate >= 0.8, missing mean-imputed.

No fabrication: every value is computed from the real .bed; QC params are
recorded with the output.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

logger = logging.getLogger(__name__)

AUTOSOME_PREFIXES = ("2L_", "2R_", "3L_", "3R_")


@dataclass
class GRMResult:
    K: np.ndarray              # (n_lines, n_lines) genomic relationship matrix
    iids: list[str]            # line IDs, row/col order of K
    n_snps_used: int
    n_snps_candidate: int
    maf_min: float
    callrate_min: float
    block_size: int

    def meta(self) -> dict:
        d = asdict(self)
        d.pop("K")
        d.pop("iids")
        d["n_lines"] = len(self.iids)
        return d


def _autosomal_snp_mask(bim_path: Path) -> np.ndarray:
    """Boolean mask over all variants: True for autosomal SNPs (by variant id).

    Reads only the variant-id column (col 1) of the .bim, in chunks, to avoid
    bed_reader's bim parser (which OOMs on this 148MB file).
    """
    ids = pd.read_csv(
        bim_path, sep="\t", header=None, usecols=[1], names=["vid"],
        dtype=str, engine="c",
    )["vid"]
    is_auto = ids.str.startswith(AUTOSOME_PREFIXES, na=False)
    is_snp = ids.str.endswith("_SNP", na=False)
    mask = (is_auto & is_snp).to_numpy()
    logger.info("Variant mask: %d autosomal SNPs of %d total", mask.sum(), mask.size)
    return mask


def build_grm(
    bed_path: str | Path,
    *,
    maf_min: float = 0.05,
    callrate_min: float = 0.8,
    block_size: int = 50_000,
    autosomes_only: bool = True,
) -> GRMResult:
    """Stream the .bed and accumulate the VanRaden GRM over QC-passing SNPs.

    For each SNP j passing QC, with alt-allele freq p_j, standardized genotype
    z_ij = (g_ij - 2 p_j) / sqrt(2 p_j (1 - p_j)) (missing -> 0 after centering).
    GRM = (1/m) * sum_j z_.j z_.j^T  (m = number of SNPs used).
    """
    bed_path = Path(bed_path)
    bed = open_bed(bed_path)
    n_line = bed.iid_count
    n_snp = bed.sid_count
    iids = [str(x) for x in bed.iid]

    include0 = (
        _autosomal_snp_mask(bed_path.with_suffix(".bim"))
        if autosomes_only
        else np.ones(n_snp, dtype=bool)
    )
    n_candidate = int(include0.sum())

    K = np.zeros((n_line, n_line), dtype=np.float64)
    m_used = 0
    for start in range(0, n_snp, block_size):
        stop = min(start + block_size, n_snp)
        cols = np.where(include0[start:stop])[0]
        if cols.size == 0:
            continue
        G = bed.read(index=np.s_[:, start:stop], dtype="float32")[:, cols]  # lines x snps
        # QC on this block
        n_obs = (~np.isnan(G)).sum(axis=0)
        callrate = n_obs / n_line
        p = np.nanmean(G, axis=0) / 2.0
        maf = np.minimum(p, 1.0 - p)
        keep = (callrate >= callrate_min) & (maf >= maf_min) & (p > 0) & (p < 1)
        if not keep.any():
            continue
        Gk = G[:, keep]
        pk = p[keep]
        # center + impute-missing-to-mean (NaN -> 0 after centering) + scale
        Z = Gk - 2.0 * pk
        Z[np.isnan(Gk)] = 0.0
        Z /= np.sqrt(2.0 * pk * (1.0 - pk))
        K += Z.astype(np.float64) @ Z.astype(np.float64).T
        m_used += int(keep.sum())
        if (start // block_size) % 10 == 0:
            logger.info("  block @%d: cum SNPs used=%d", start, m_used)

    if m_used == 0:
        raise RuntimeError("No SNPs passed QC — check thresholds / input.")
    K /= m_used
    logger.info("GRM built: %d lines, %d SNPs used (of %d candidate).",
                n_line, m_used, n_candidate)
    return GRMResult(K, iids, m_used, n_candidate, maf_min, callrate_min, block_size)


def build_snp_matrix(
    bed_path: str | Path,
    *,
    maf_min: float = 0.05,
    callrate_min: float = 0.8,
    target_snps: int = 100_000,
    block_size: int = 50_000,
    autosomes_only: bool = True,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Return an in-memory genotype matrix (lines x SNPs), QC'd + evenly thinned.

    For marker-effect models (Option 2). QC-passing autosomal SNPs are thinned by
    even spacing to ~target_snps to fit in memory (205 x 100k x 4B ≈ 82MB).
    Missing values mean-imputed; values left as 0/1/2 dosage (not standardized,
    so downstream models can scale within CV folds without leakage).

    Returns (X, iids, snp_index) where snp_index are the kept variant indices.
    """
    bed_path = Path(bed_path)
    bed = open_bed(bed_path)
    n_line, n_snp = bed.iid_count, bed.sid_count
    iids = [str(x) for x in bed.iid]
    include0 = (_autosomal_snp_mask(bed_path.with_suffix(".bim"))
                if autosomes_only else np.ones(n_snp, dtype=bool))

    # first pass: indices of QC-passing SNPs (stream)
    passing: list[int] = []
    for start in range(0, n_snp, block_size):
        stop = min(start + block_size, n_snp)
        cols = np.where(include0[start:stop])[0]
        if cols.size == 0:
            continue
        G = bed.read(index=np.s_[:, start:stop], dtype="float32")[:, cols]
        callrate = (~np.isnan(G)).mean(axis=0)
        p = np.nanmean(G, axis=0) / 2.0
        maf = np.minimum(p, 1 - p)
        keep = (callrate >= callrate_min) & (maf >= maf_min) & (p > 0) & (p < 1)
        passing.extend((start + cols[keep]).tolist())
    passing = np.array(passing)
    if passing.size > target_snps:                      # even thinning
        passing = passing[np.linspace(0, passing.size - 1, target_snps).astype(int)]
    logger.info("SNP matrix: %d SNPs kept (thinned from QC-passing)", passing.size)

    # second pass: read kept SNPs in blocks, mean-impute
    X = np.empty((n_line, passing.size), dtype=np.float32)
    for j0 in range(0, passing.size, block_size):
        idx = passing[j0:j0 + block_size]
        G = bed.read(index=np.s_[:, idx], dtype="float32")
        col_mean = np.nanmean(G, axis=0)
        nan = np.isnan(G)
        G[nan] = np.take(col_mean, np.where(nan)[1])
        X[:, j0:j0 + idx.size] = G
    return X, iids, passing


def save_grm(res: GRMResult, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, K=res.K, iids=np.array(res.iids),
        n_snps_used=res.n_snps_used, n_snps_candidate=res.n_snps_candidate,
        maf_min=res.maf_min, callrate_min=res.callrate_min, block_size=res.block_size,
    )


def load_grm(path: str | Path) -> tuple[np.ndarray, list[str], dict]:
    d = np.load(path, allow_pickle=True)
    iids = [str(x) for x in d["iids"]]
    meta = {k: d[k].item() for k in
            ("n_snps_used", "n_snps_candidate", "maf_min", "callrate_min", "block_size")}
    return d["K"], iids, meta
