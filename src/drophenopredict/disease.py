"""Disease classifier — maps a query gene/mutation to a neurodegenerative-disease
category (honest v2 scaffold).

This is CURATION, not a learned/predictive model: it looks a gene up in a
manually curated table of well-established human disease genes (with confidently
known Drosophila orthologs where they exist). Sources are standard disease
genetics (OMIM; ALSoD; MDSGene; and review literature). Nothing here is fabricated
or fitted — it is a lookup used to route a query to the right pathway set.

Where a human gene has NO fly ortholog (SNCA, HTT, MAPT), disease modelling in
flies uses a human transgene — noted in FLY_ORTHOLOG as None.
"""
from __future__ import annotations

# Curated human disease genes per category (well-established; not exhaustive).
DISEASE_GENES: dict[str, list[str]] = {
    "ALS": ["TARDBP", "FUS", "SOD1", "C9orf72", "VCP", "OPTN", "UBQLN2",
            "TBK1", "CHCHD10", "VAPB", "ANG", "SETX", "PFN1", "MATR3"],
    "PD": ["SNCA", "PINK1", "PRKN", "PARK2", "LRRK2", "PARK7", "DJ1",
           "VPS35", "GBA", "ATP13A2"],
    "HD": ["HTT"],
    "AD": ["APP", "PSEN1", "PSEN2", "MAPT", "APOE", "TREM2", "GRN"],
    "SCA": ["ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "TBP", "ATN1"],
}

# Confidently-known human -> Drosophila ortholog (None = no ortholog; human
# transgene used in fly models). Fly symbols per FlyBase/DIOPT.
FLY_ORTHOLOG: dict[str, str | None] = {
    "TARDBP": "TBPH", "FUS": "caz", "SOD1": "Sod1", "VCP": "TER94",
    "OPTN": None, "UBQLN2": "Ubqn", "C9orf72": None, "CHCHD10": None,
    "SNCA": None, "PINK1": "Pink1", "PRKN": "park", "PARK2": "park",
    "LRRK2": "Lrrk", "PARK7": "DJ-1beta", "DJ1": "DJ-1beta", "GBA": "Gba1b",
    "HTT": None, "APP": "Appl", "PSEN1": "Psn", "PSEN2": "Psn", "MAPT": None,
    "ATXN1": None, "ATXN2": "Atx2", "ATXN3": "Atx3", "CACNA1A": "cac", "TBP": "Tbp",
}

# Common fly-model aliases -> canonical human gene (so fly-symbol queries route too).
FLY_ALIAS: dict[str, str] = {
    "TBPH": "TARDBP", "CAZ": "FUS", "SOD1": "SOD1", "TER94": "VCP",
    "PINK1": "PINK1", "PARK": "PRKN", "LRRK": "LRRK2", "DJ-1BETA": "PARK7",
    "APPL": "APP", "PSN": "PSEN1", "ATX2": "ATXN2", "ATX3": "ATXN3",
    "ALPHA-SYNUCLEIN": "SNCA", "ASYN": "SNCA", "TAU": "MAPT", "HTT": "HTT",
}

CATEGORIES = list(DISEASE_GENES) + ["OTHER"]


def _canonical(gene: str) -> str:
    g = gene.upper().strip().replace("_", "-")
    return FLY_ALIAS.get(g, g)


def classify(gene: str, mutation: str = "", disease_hint: str = "") -> dict:
    """Return a category label + evidence for a query gene.

    Deterministic curation — NOT a probabilistic ML output. Returns:
      {category, matched_gene, human_gene, fly_ortholog, basis, note}
    `category` is 'OTHER' when the gene is not in the curated table.
    """
    if disease_hint:
        h = disease_hint.upper().strip()
        if h in DISEASE_GENES:
            return {"category": h, "matched_gene": None, "human_gene": None,
                    "fly_ortholog": None, "basis": "user-provided disease hint",
                    "note": "hint trusted; gene not cross-checked"}

    g = _canonical(gene)
    for cat, genes in DISEASE_GENES.items():
        if g in genes:
            note = ""
            # HD severity signal: pathogenic polyQ length (CAG/Q > 35)
            if cat == "HD":
                import re
                m = re.search(r"Q(\d+)", mutation.upper())
                if m:
                    note = f"polyQ length {m.group(1)} " + \
                           ("(pathogenic)" if int(m.group(1)) > 35 else "(sub-threshold)")
            return {"category": cat, "matched_gene": g, "human_gene": g,
                    "fly_ortholog": FLY_ORTHOLOG.get(g, "unknown"),
                    "basis": "curated disease-gene table (OMIM/ALSoD/MDSGene)",
                    "note": note}
    return {"category": "OTHER", "matched_gene": None, "human_gene": g,
            "fly_ortholog": None, "basis": "not in curated table",
            "note": "unknown gene -> route to shared pathways only"}
