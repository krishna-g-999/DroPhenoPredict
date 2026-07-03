"""Build the BCM-DMAS multi-disease fly expression matrix (Synapse syn34767207).

Assembles the 148 featureCounts files (FlyBase r6.06, FBgn genes) into a
genes x samples matrix, joined to the biospecimen metadata key
(specimenID -> genotype, ageDays). This is a REAL 5-disease-model matched-omics
resource on one genetic background (Farrell/Shulman, npj PD 2025):
  controls: Elav-Gal4/+, w1118/+ ; AD: Abeta42, Tau2N4R ; PD: aSyn ;
  HD: HTT-128Q, HTT-200Q.

NOTE: the number in a sample name (A10) is an internal tube id, NOT the age.
Age comes from the metadata `ageDays` column only.

Output: data/processed/bcm_dmas_counts.parquet  (genes x samples)
        data/processed/bcm_dmas_samples.csv       (sample -> genotype, disease, age)
"""
from __future__ import annotations

import io
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ZIP = "C:/Users/gkp11/Downloads/syn36328638.zip"
META = Path("data/raw/bcm_dmas/BCM-DMAS_biospecimen_metadata.csv")
PROC = Path("data/processed")

# map the verbose genotype string -> (disease, short label)
DISEASE = {
    "Negative Control (driver)": ("control", "Elav-Gal4"),
    "Negative Control (no driver)": ("control", "w1118"),
    "AD-model-Secreted beta42": ("AD", "Abeta42"),
    "AD-model-Tauopathy": ("AD", "Tau"),
    "PD model-Synuclein": ("PD", "aSyn"),
    "HD model-NT-HTT-128Q": ("HD", "HTT128Q"),
    "HD model-FL-HTT-200Q": ("HD", "HTT200Q"),
}


def label(genotype: str):
    for key, (dis, short) in DISEASE.items():
        if key in genotype:
            return dis, short
    return "other", genotype[:20]


def read_counts(zf, name) -> pd.Series:
    with zf.open(name) as fh:
        df = pd.read_csv(io.TextIOWrapper(fh, "utf-8"), sep="\t", comment="#")
    # featureCounts: Geneid, Chr, Start, End, Strand, Length, <count col last>
    s = df.set_index("Geneid").iloc[:, -1]
    s.name = name.replace(".counts.txt", "")
    return s


def main() -> None:
    meta = pd.read_csv(META)
    meta = meta[meta["assay"] == "rnaSeq"].copy()
    meta[["disease", "model"]] = meta["genotype"].apply(lambda g: pd.Series(label(g)))
    key = meta.set_index("specimenID")

    zf = zipfile.ZipFile(ZIP)
    files = [n for n in zf.namelist() if n.endswith(".counts.txt")]
    cols = {}
    for n in files:
        sid = n.replace(".counts.txt", "")
        if sid in key.index:
            cols[sid] = read_counts(zf, n)
    mat = pd.DataFrame(cols)                       # genes x samples
    print(f"expression matrix: {mat.shape[0]} genes x {mat.shape[1]} samples")

    samp = key.loc[mat.columns, ["disease", "model", "ageDays", "genotype"]].copy()
    samp.index.name = "specimenID"
    PROC.mkdir(parents=True, exist_ok=True)
    mat.to_parquet(PROC / "bcm_dmas_counts.parquet")
    samp.to_csv(PROC / "bcm_dmas_samples.csv")

    print("\n=== samples per disease model x age ===")
    tab = samp.groupby(["disease", "model", "ageDays"]).size().rename("n").reset_index()
    print(tab.pivot_table(index=["disease", "model"], columns="ageDays",
                          values="n", fill_value=0).to_string())
    print(f"\nSaved -> {PROC}/bcm_dmas_counts.parquet + bcm_dmas_samples.csv")


if __name__ == "__main__":
    main()
