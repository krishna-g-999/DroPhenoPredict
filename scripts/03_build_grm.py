"""Build and cache the DGRP genomic relationship matrix (GRM) for GBLUP.

Streams the 4.4M-variant PLINK .bed, applies QC (autosomal SNPs, MAF>=0.05,
call-rate>=0.8), and writes data/processed/grm.npz (205x205 K + metadata).
Run once: ./.venv/Scripts/python.exe scripts/03_build_grm.py
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import genotypes  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

BED = "data/raw/dgrp2/dgrp2.bed"
OUT = Path("data/processed/grm.npz")


def main() -> None:
    res = genotypes.build_grm(BED, maf_min=0.05, callrate_min=0.8, block_size=50_000)
    genotypes.save_grm(res, OUT)
    print("\nGRM saved ->", OUT)
    print(json.dumps(res.meta(), indent=2))
    # sanity: diagonal ~ 1 (self-relationship), symmetric
    import numpy as np
    K = res.K
    print(f"K shape={K.shape}  diag mean={np.diag(K).mean():.3f} "
          f"(expect ~1)  offdiag mean={K[~np.eye(len(K),dtype=bool)].mean():.3f}  "
          f"symmetric={np.allclose(K,K.T)}")


if __name__ == "__main__":
    main()
