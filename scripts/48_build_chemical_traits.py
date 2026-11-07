"""Build one harmonised table of DGRP drug/toxin response traits from all collected sources.

Inputs (downloaded by scripts 45-47, never edited):
  data/raw/dgrpool/phenotype_<id>.json                     DGRPool values
  data/raw/chemical_repo/battlay2018_insecticide/...      malathion, permethrin
  data/raw/chemical_repo/duneau2018_insecticides/...      parathion (per fly), deltamethrin (per line)
  data/raw/chemical_supp/denecke2017_imidacloprid/...     imidacloprid wiggle index
  data/raw/chemical_supp/mitchell2017_amanitin/...        alpha-amanitin
  data/raw/chemical_supp/harrison2024_rapamycin/...       rapamycin
  data/raw/chemical_supp/marconcini2026_toxin/...         octanoic acid (noni toxin)

Rules
  - One value per line per trait-sex: the mean of the values the source reports for that
    line. The only derived quantities are documented below (proportions from counts).
  - Line ids harmonised to the GRM form 'line_<n>'; only lines in the genotyped
    205-line DGRP2 panel are kept, and the number dropped is recorded.
  - role = 'treated' (exposure phenotype) or 'control' (same assay without the compound).

Outputs
  data/processed/chemical_traits_long.csv     trait_id, line, value
  data/processed/chemical_traits_catalog.csv  one row per trait-sex with source, role, n
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import dgrpool  # noqa: E402

RAW_REPO = Path("data/raw/chemical_repo")
RAW_SUPP = Path("data/raw/chemical_supp")
PANEL = set(np.load("data/processed/grm.npz", allow_pickle=True)["iids"].tolist())

# DGRPool phenotypes: (phenotype_id, trait label, compound, role, higher_is)
# higher_is says what a larger value means, to keep directions interpretable.
DGRPOOL = [
    (24, "azinphos_LD50", "azinphos-methyl", "treated", "resistant"),
    (1544, "climbing_MSB", "menadione sodium bisulfite", "treated", "better climbing"),
    (1542, "startle_MSB", "menadione sodium bisulfite", "treated", "longer startle"),
    (1543, "climbing_MSB_ctrl", "menadione sodium bisulfite", "control", "better climbing"),
    (1541, "startle_MSB_ctrl", "menadione sodium bisulfite", "control", "longer startle"),
    (1570, "MeHg_cum_survival_idx", "methylmercury", "treated", "tolerant"),
    (2640, "MeHg10_survival", "methylmercury", "treated", "tolerant"),
    (2643, "MeHg10_caffeine_survival", "methylmercury + caffeine", "treated", "tolerant"),
    (2642, "caffeine_dev_survival", "caffeine", "treated", "tolerant"),
    (2638, "MeHg_ctrl_survival", "methylmercury", "control", "tolerant"),
    (1589, "alcohol_sensitivity_E1", "ethanol", "treated", "resistant (longer elution)"),
    (1591, "alcohol_sensitivity_E2", "ethanol", "treated", "resistant (longer elution)"),
    (1593, "alcohol_tolerance", "ethanol", "treated", "more tolerance"),
    (2800, "caffeine_lifespan", "caffeine", "treated", "longer survival"),
    (2801, "boric_acid_lifespan", "boric acid", "treated", "longer survival"),
    (1720, "paraquat_resistance", "paraquat", "treated", "resistant"),
    (1721, "menadione_resistance", "menadione sodium bisulfite", "treated", "resistant"),
    (2843, "lead_dev_time", "lead", "treated", "slower development"),
    (2844, "lead_viability", "lead", "treated", "higher viability"),
    (2846, "lead_activity", "lead", "treated", "more active"),
    (2841, "lead_dev_time_ctrl", "lead", "control", "slower development"),
    (2842, "lead_viability_ctrl", "lead", "control", "higher viability"),
    (2845, "lead_activity_ctrl", "lead", "control", "more active"),
    (2885, "methylphenidate_food_intake", "methylphenidate", "treated", "more intake"),
    (2884, "methylphenidate_activity", "methylphenidate", "treated", "more active"),
    (2891, "methylphenidate_food_intake_ctrl", "methylphenidate", "control", "more intake"),
    (2890, "methylphenidate_activity_ctrl", "methylphenidate", "control", "more active"),
    (2956, "ethanol_dev_viability", "ethanol (developmental)", "treated", "higher viability"),
    (2958, "ethanol_dev_time", "ethanol (developmental)", "treated", "slower development"),
    (2955, "ethanol_dev_viability_ctrl", "ethanol (developmental)", "control", "higher viability"),
    (2957, "ethanol_dev_time_ctrl", "ethanol (developmental)", "control", "slower development"),
    (2410, "paraquat_climbing", "paraquat (low dose)", "treated", "better climbing"),
    (2402, "paraquat_climbing_ctrl", "paraquat (low dose)", "control", "better climbing"),
    (2502, "ryanoid_LC50", "ryanoid insecticide", "treated", "resistant"),
]
DGRPOOL_SOURCE = {24: 3, 1544: 19, 1542: 19, 1543: 19, 1541: 19, 1570: 23, 2640: 23, 2643: 23,
                  2642: 23, 2638: 23, 1589: 25, 1591: 25, 1593: 25, 2800: 26, 2801: 27, 1720: 36,
                  1721: 36, 2843: 37, 2844: 37, 2846: 37, 2841: 37, 2842: 37, 2845: 37, 2885: 43,
                  2884: 43, 2891: 43, 2890: 43, 2956: 51, 2958: 51, 2955: 51, 2957: 51, 2410: 59,
                  2402: 59, 2502: 98}


def to_line(x) -> str | None:
    """'DGRP_021', 'DGRP100', 'line_913', 'RAL-208', 100, '100.0' -> 'line_<n>'."""
    m = re.search(r"(\d{1,3})(?:\.0+)?\s*$", str(x).strip())
    return f"line_{int(m.group(1))}" if m else None


rows: list[dict] = []
catalog: list[dict] = []


def add(trait: str, sex: str, compound: str, role: str, higher_is: str, source: str,
        series: pd.Series, note: str = "") -> None:
    s = series.dropna()
    s.index = [to_line(i) for i in s.index]
    s = s[[i is not None for i in s.index]]
    s = s.groupby(level=0).mean()  # one value per line
    n_all = len(s)
    s = s[s.index.isin(PANEL)]
    tid = f"{trait}_{sex}"
    for line, v in s.items():
        rows.append({"trait_id": tid, "line": line, "value": float(v)})
    catalog.append({"trait_id": tid, "trait": trait, "sex": sex, "compound": compound,
                    "role": role, "higher_is": higher_is, "source": source,
                    "n_lines_reported": n_all, "n_lines_in_panel": len(s), "note": note})


def from_dgrpool() -> None:
    studies = dgrpool.fetch_studies()
    for pid, trait, compound, role, hi in DGRPOOL:
        vals = dgrpool.fetch_phenotype_values(pid)
        st = studies.loc[DGRPOOL_SOURCE[pid]]
        src = f"DGRPool phenotype {pid}; {st['first_author']} {st['year']} (doi {st['doi']})"
        for sex in sorted(vals["sex"].unique()):
            sx = "B" if sex == "NA" else sex  # B = sexes not separated
            add(trait, sx, compound, role, hi, src,
                vals[vals.sex == sex].groupby("DGRP_line")["value"].mean())


def from_battlay2018() -> None:
    f = RAW_REPO / "battlay2018_insecticide/SupplementaryData.xlsx"
    sheets = {"mal.3.pheno": ("malathion_3h_mortality", "malathion"),
              "mal.6.pheno": ("malathion_6h_mortality", "malathion"),
              "mal.12.pheno": ("malathion_12h_mortality", "malathion"),
              "mal.24.pheno": ("malathion_24h_mortality", "malathion"),
              "per.3KD.pheno": ("permethrin_3h_knockdown", "permethrin"),
              "per.24.pheno": ("permethrin_24h_mortality", "permethrin")}
    for sh, (trait, compound) in sheets.items():
        df = pd.read_excel(f, sheet_name=sh, header=None)
        hdr = df.index[df.iloc[:, 0].astype(str).str.contains("DGRP line", na=False)][0]
        body = df.iloc[hdr + 1:, :3]
        body.columns = ["line", "M", "F"]
        body = body.set_index("line")
        for sex in ("F", "M"):
            add(trait, sex, compound, "treated", "more susceptible",
                "Battlay et al. 2018, G3 (figshare 10.25387/g3.6954443)",
                pd.to_numeric(body[sex], errors="coerce"))


def from_duneau2018() -> None:
    d = RAW_REPO / "duneau2018_insecticides"
    p = pd.read_csv(d / "Duneau_G3_TableS8.csv", sep=";", skiprows=2, dtype=str)
    # Censor 1 = death observed, 0 = alive at 48 h. Resistance = proportion alive at 48 h.
    alive = (p["Censor"].astype(int) == 0).groupby(p["DGRP_lines"]).mean()
    add("parathion_survival_48h", "B", "parathion", "treated", "resistant",
        "Duneau et al. 2018, G3 (figshare 10.25387/g3.7047734), Table S8",
        alive, "proportion of flies alive at 48 h, pooled over experiments")
    t = pd.read_csv(d / "Duneau_G3_TableS9.csv", sep=";", skiprows=2, dtype=str)
    t = t.dropna(subset=["DGRP_lines"])
    a = pd.to_numeric(t["Delta_Alive.at.48h"], errors="coerce").groupby(t["DGRP_lines"]).sum()
    n = pd.to_numeric(t["Sample.size"], errors="coerce").groupby(t["DGRP_lines"]).sum()
    add("deltamethrin_survival_48h", "B", "deltamethrin", "treated", "resistant",
        "Duneau et al. 2018, G3 (figshare 10.25387/g3.7047734), Table S9",
        a / n, "alive / exposed at 48 h, pooled over experiments")


def from_denecke2017() -> None:
    f = RAW_SUPP / "denecke2017_imidacloprid/files/41598_2017_11092_MOESM2_ESM.xlsx"
    for sh, trait in (("TableS1", "imidacloprid_25ppm_wiggle"), ("TableS2", "imidacloprid_100ppm_wiggle")):
        df = pd.read_excel(f, sheet_name=sh, header=None)
        df = df[pd.to_numeric(df[0], errors="coerce").notna()]
        add(trait, "B", "imidacloprid", "treated", "less affected (more movement)",
            "Denecke et al. 2017, Sci Rep, Supplementary Tables S1-S2",
            pd.Series(pd.to_numeric(df[1], errors="coerce").values, index=df[0].values),
            "larval wiggle index relative to untreated")


def from_mitchell2017() -> None:
    f = RAW_SUPP / "mitchell2017_amanitin/files/pone.0173162.s001.xlsx"
    df = pd.read_excel(f, header=None)
    body = df.iloc[2:]
    # two side-by-side blocks, each with its own strain column (0 and 12) and average (5 and 17)
    for id_col, col, trait in ((0, 5, "amanitin_0.2ug_eclosion"), (12, 17, "amanitin_2ug_eclosion")):
        b = body[pd.to_numeric(body[id_col], errors="coerce").notna()]
        add(trait, "B", "alpha-amanitin", "treated", "resistant (more adults hatched)",
            "Mitchell et al. 2017, PLoS One, S1 Table",
            pd.Series(pd.to_numeric(b[col], errors="coerce").values, index=b[id_col].values),
            "mean adults eclosed per replicate")


def from_harrison2024() -> None:
    f = RAW_SUPP / "harrison2024_rapamycin/files/ACEL-23-e14292-s003.xlsx"
    df = pd.read_excel(f, header=None)
    body = df.iloc[2:]
    body = body[pd.to_numeric(body[0], errors="coerce").notna()]
    idx = body[0].values
    add("rapamycin_pupation_time", "B", "rapamycin", "treated", "slower development",
        "Harrison et al. 2024, Aging Cell, Supplementary Table 1",
        pd.Series(pd.to_numeric(body[9], errors="coerce").values, index=idx))
    add("rapamycin_pupation_time_ctrl", "B", "rapamycin", "control", "slower development",
        "Harrison et al. 2024, Aging Cell, Supplementary Table 1",
        pd.Series(pd.to_numeric(body[2], errors="coerce").values, index=idx), "ethanol vehicle")


def from_marconcini2026() -> None:
    f = RAW_SUPP / "marconcini2026_toxin/files/nested_jkag032_supplementary_data/File_S1_G3-2026-406569.xlsx"
    for sh, sex in (("DGRP_females", "F"), ("DGRP_males", "M")):
        raw = pd.read_excel(f, sheet_name=sh, header=None)
        hdr = raw.index[raw[0].astype(str).str.strip() == "Line"][0]  # header row differs by sheet
        df = raw.iloc[hdr + 1:].set_axis(raw.iloc[hdr].tolist(), axis=1).dropna(subset=["Line"])
        surv = pd.to_numeric(df["survival"], errors="coerce").groupby(df["Line"]).mean()
        add("octanoic_acid_survival", sex, "octanoic acid (noni toxin)", "treated", "resistant",
            "Marconcini et al. 2026, G3, File S1", surv, "mean replicate survival")


def main() -> None:
    from_dgrpool()
    from_battlay2018()
    from_duneau2018()
    from_denecke2017()
    from_mitchell2017()
    from_harrison2024()
    from_marconcini2026()
    long = pd.DataFrame(rows)
    cat = pd.DataFrame(catalog).sort_values(["compound", "trait", "sex"])
    Path("data/processed").mkdir(exist_ok=True)
    long.to_csv("data/processed/chemical_traits_long.csv", index=False)
    cat.to_csv("data/processed/chemical_traits_catalog.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_rows", 200, "display.max_colwidth", 34):
        print(cat[["trait_id", "compound", "role", "n_lines_reported", "n_lines_in_panel"]].to_string(index=False))
    print(f"\n{cat.trait_id.nunique()} trait-sex combinations, {cat.compound.nunique()} compounds, "
          f"{long.line.nunique()} panel lines, {len(long)} values")


if __name__ == "__main__":
    main()
