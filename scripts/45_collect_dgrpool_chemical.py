"""Collect every DGRPool phenotype from studies of drug, toxin or chemical exposure.

Study selection: DGRPool studies whose title/abstract mention a drug, toxin,
xenobiotic, insecticide, metal, alcohol or oxidant exposure (keyword search of
the live catalog, 2026-11-07), reviewed by hand. Diet-only and infection studies
are excluded: they are not chemical exposures.

Output
  data/raw/dgrpool/phenotype_<id>.json        raw API responses (cache, untouched)
  data/interim/chemical/dgrpool_manifest.csv   one row per phenotype: study, name,
                                               flags, n lines per sex
  data/interim/chemical/dgrpool_values.csv     tidy per-line values (line, sex, value)

Nothing is imputed or edited; lines and values are exactly as DGRPool serves them.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import dgrpool  # noqa: E402

# DGRPool study id -> exposure (hand-checked against each study's title/abstract)
CHEMICAL_STUDIES = {
    3: "insecticide (azinphos-methyl)",
    19: "oxidative stress (menadione sodium bisulfite, chronic)",
    23: "methylmercury (+ caffeine)",
    25: "alcohol (ethanol)",
    26: "caffeine",
    27: "boric acid",
    36: "oxidative stress (paraquat, menadione)",
    37: "lead",
    43: "drug (methylphenidate) vs sucrose",
    51: "alcohol (ethanol, developmental)",
    59: "paraquat (low dose)",
    63: "insecticide (Duneau 2018)",
    85: "hydrogen peroxide",
    95: "rapamycin",
    98: "insecticide (Green 2018)",
}

OUT = Path("data/interim/chemical")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    studies = dgrpool.fetch_studies()
    catalog = dgrpool.fetch_phenotype_catalog()

    manifest, values = [], []
    for sid, exposure in CHEMICAL_STUDIES.items():
        s = studies.loc[sid]
        for pid, ph in catalog[catalog.study_id == sid].iterrows():
            try:
                df = dgrpool.fetch_phenotype_values(int(pid))
            except RuntimeError as e:
                print(f"  FAILED phenotype {pid}: {e}")
                manifest.append({"study_id": sid, "phenotype_id": pid, "name": ph["name"],
                                 "fetch_error": str(e)})
                continue
            n_by_sex = df.groupby("sex")["DGRP_line"].nunique().to_dict()
            manifest.append({
                "study_id": sid, "first_author": s["first_author"], "year": s["year"],
                "doi": s["doi"], "exposure": exposure, "phenotype_id": pid,
                "name": ph["name"], "unit": ph.get("unit"),
                "is_numeric": ph.get("is_numeric"), "is_continuous": ph.get("is_continuous"),
                "is_summary": ph.get("is_summary"), "obsolete": ph.get("obsolete"),
                "n_lines_total": df["DGRP_line"].nunique(),
                "n_lines_F": n_by_sex.get("F", 0), "n_lines_M": n_by_sex.get("M", 0),
                "n_lines_NA": n_by_sex.get("NA", 0),
            })
            if len(df):
                df.insert(0, "phenotype_id", pid)
                values.append(df)

    man = pd.DataFrame(manifest)
    man.to_csv(OUT / "dgrpool_manifest.csv", index=False)
    if values:
        pd.concat(values).to_csv(OUT / "dgrpool_values.csv", index=False)
    print(f"{len(man)} phenotypes from {len(CHEMICAL_STUDIES)} studies")
    cols = ["study_id", "first_author", "exposure", "phenotype_id", "name",
            "is_continuous", "n_lines_F", "n_lines_M", "n_lines_NA"]
    with pd.option_context("display.width", 250, "display.max_rows", 300,
                           "display.max_colwidth", 40):
        print(man[[c for c in cols if c in man.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
