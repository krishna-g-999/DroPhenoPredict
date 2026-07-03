"""Measure REAL DGRPool line coverage for candidate target traits.

Answers the gating question from SCIENTIFIC_CHARTER.md §5.1:
  "How many DGRP lines actually share the candidate target traits?"
The original brief assumed 'all 205 for all traits' — that is false. This
script measures the truth from live data and writes a report. No assumptions.

Run:  ./.venv/Scripts/python.exe scripts/01_measure_coverage.py
Output: data/processed/coverage_report.csv  +  console summary.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import dgrpool  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

# Trait families we care about, matched by case-insensitive substring against
# the REAL phenotype names. We report everything matched — we do not pre-judge.
TRAIT_KEYWORDS = {
    "climbing": ["climb", "geotax"],
    "locomotion_startle": ["startle", "locomot", "activity"],
    "sleep": ["sleep"],
    "lifespan": ["lifespan", "longevity", "lt50", "survival"],
    "starvation": ["starv"],
}

OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)


def match_family(name: str) -> str | None:
    low = name.lower()
    for fam, kws in TRAIT_KEYWORDS.items():
        if any(k in low for k in kws):
            return fam
    return None


def main() -> None:
    catalog = dgrpool.fetch_phenotype_catalog()
    cont = dgrpool.numeric_continuous_phenotypes(catalog)
    print(f"\nPhenotype catalog: {len(catalog)} total, {len(cont)} numeric-continuous usable.\n")

    # Tag candidates by family
    cont = cont.copy()
    cont["family"] = cont["name"].map(match_family)
    candidates = cont[cont["family"].notna()]
    print(f"Candidate target phenotypes matched: {len(candidates)} across "
          f"{candidates['family'].nunique()} families.\n")

    # Fetch per-line values for each candidate; record line coverage.
    rows = []
    line_sets: dict[int, set[str]] = {}
    for pid, row in candidates.iterrows():
        try:
            vals = dgrpool.fetch_phenotype_values(int(pid))
        except Exception as e:  # noqa: BLE001
            logging.warning("skip phenotype %s: %s", pid, e)
            continue
        lines = set(vals["DGRP_line"].unique())
        line_sets[int(pid)] = lines
        rows.append({
            "phenotype_id": int(pid),
            "study_id": int(row["study_id"]),
            "family": row["family"],
            "name": row["name"],
            "unit": row.get("unit"),
            "n_lines": len(lines),
            "n_obs": len(vals),
            "sexes": ",".join(sorted(vals["sex"].unique())) if len(vals) else "",
        })

    report = pd.DataFrame(rows).sort_values(["family", "n_lines"], ascending=[True, False])
    report_path = OUT / "coverage_report.csv"
    report.to_csv(report_path, index=False)
    print(f"Wrote per-phenotype coverage -> {report_path}\n")

    # Per-family best phenotype (max line coverage) and the multi-trait overlap.
    print("=== Best-covered phenotype per family ===")
    best: dict[str, int] = {}
    for fam, grp in report.groupby("family"):
        top = grp.iloc[0]
        best[fam] = int(top["phenotype_id"])
        print(f"  {fam:<20} id={top['phenotype_id']:<5} n_lines={top['n_lines']:<4} "
              f"{top['name'][:55]}")

    print("\n=== Multi-trait line overlap (using best phenotype per family) ===")
    fams = list(best.keys())
    if fams:
        joint = set.intersection(*[line_sets[best[f]] for f in fams])
        print(f"  Families: {fams}")
        print(f"  Lines with ALL {len(fams)} families measured: {len(joint)}")
        for f in fams:
            print(f"    {f:<20} alone: {len(line_sets[best[f]])} lines")
        # pairwise so we can see which trait is the bottleneck
        print("\n  Pairwise overlap:")
        for i, a in enumerate(fams):
            for b in fams[i + 1:]:
                ov = len(line_sets[best[a]] & line_sets[best[b]])
                print(f"    {a:<20} & {b:<20} = {ov}")

    print("\nDONE. This overlap is the REAL upper bound on a unified multi-trait "
          "training set. Interpret v1 scope (Charter §5) against it.")


if __name__ == "__main__":
    main()
