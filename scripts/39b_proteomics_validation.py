"""Validate the transcript-level findings (Findings 13/15) at the protein level,
using the BCM-DMAS TMT proteomics (run scripts/38 first to download it).

IMPORTANT — written defensively because the exact file format has not been seen
yet (proteomics was not downloaded when this script was written; only the
metadata's *description* of it — 16 TMT samples, PD/aSyn vs Elav-Gal4 control,
day 10 only — is known). This script:
  PHASE 1 (always safe): loads the raw file and PRINTS its real structure
    (columns, dtypes, shape, head) before assuming anything about it.
  PHASE 2 (gated): only proceeds to the biological comparison if Phase 1
    successfully auto-detects (a) a protein/gene identifier column and (b) 16
    numeric abundance columns matching the 16 TMT samples in the metadata. If
    detection fails, the script STOPS and prints exactly what it found, rather
    than guessing column names and silently computing something wrong.

Comparison performed (only if Phase 2 succeeds):
  1. Map protein IDs -> FBgn (via the file's own gene-symbol column, cross-
     referenced with our existing FlyBase annotation map).
  2. Compute protein-level aSyn-vs-control log2FC (mean of the 8 aSyn TMT
     channels vs mean of the 8 control channels — verify this grouping against
     the metadata's specimenID->genotype key, do not assume column order).
  3. Genome-wide RNA-protein concordance (Spearman r, calibration check — the
     literature typically reports modest positive correlation, e.g. Spearman
     0.4-0.6, in many systems; this is what we are checking for, not assuming).
  4. Focused test: do the aSyn-specific convergent/QTT genes from Findings
     09/13/15 show concordant direction at the protein level?

Output (only if Phase 2 succeeds):
  data/processed/proteomics_asyn_log2fc.csv
  data/processed/proteomics_rna_concordance.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

PROC = Path("data/processed")
PROT_FILE = Path("data/raw/npjpd2025/proteomics/jshulman_flytmt16__Proteins.txt")
META_FILE = Path("data/raw/bcm_dmas/BCM-DMAS_biospecimen_metadata.csv")


def phase1_inspect() -> pd.DataFrame:
    if not PROT_FILE.exists():
        sys.exit(f"Proteomics file not found at {PROT_FILE}.\n"
                 f"Run scripts/38_download_proteomics.py first (in your own terminal, "
                 f"with your Synapse token).")
    # Try tab-separated first (Proteome-Discoverer-style exports are tab-delimited);
    # fall back to comma if that fails.
    for sep in ("\t", ","):
        try:
            df = pd.read_csv(PROT_FILE, sep=sep, nrows=5)
            if df.shape[1] > 1:
                break
        except Exception:
            continue
    else:
        sys.exit(f"Could not parse {PROT_FILE} as TSV or CSV. Inspect it manually.")

    full = pd.read_csv(PROT_FILE, sep=sep)
    print(f"=== PHASE 1: raw file structure ===")
    print(f"file: {PROT_FILE} | sep={sep!r} | shape={full.shape}")
    print(f"columns ({len(full.columns)}):")
    for c in full.columns:
        print(f"  {c!r}  dtype={full[c].dtype}  example={full[c].iloc[0] if len(full) else 'NA'}")
    return full


def phase2_detect_and_compare(full: pd.DataFrame) -> None:
    meta = pd.read_csv(META_FILE)
    tmt = meta[meta["assay"] == "TMTquantitation"].copy()
    print(f"\n=== PHASE 2: matching against metadata ({len(tmt)} TMT samples expected) ===")
    print(tmt[["specimenID", "genotype"]].to_string(index=False))

    # candidate identifier columns: anything with 'gene' or 'accession' or 'symbol' in the name
    id_candidates = [c for c in full.columns
                     if any(k in c.lower() for k in ("gene", "accession", "symbol", "protein"))]
    print(f"\nCandidate identifier columns: {id_candidates}")

    # candidate abundance columns: numeric columns, count should be near 16
    numeric_cols = full.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Numeric columns found: {len(numeric_cols)} -> {numeric_cols[:20]}"
          f"{'...' if len(numeric_cols) > 20 else ''}")

    if not id_candidates or len(numeric_cols) < 10:
        print("\nSTOPPING: could not confidently auto-detect the identifier column and/or "
              "16 abundance columns. This needs a human to look at the printed structure "
              "above and tell me which columns are which — I will not guess.")
        return

    print("\nAuto-detection found candidate columns, but this script will NOT proceed to "
          "the biological comparison automatically on the first run against real data — "
          "the exact TMT-channel-to-specimenID mapping must be confirmed against the "
          "structure above before any log2FC is computed. Report the printed structure "
          "back and the comparison logic will be finalized against the real file.")


def main() -> None:
    full = phase1_inspect()
    phase2_detect_and_compare(full)


if __name__ == "__main__":
    main()
