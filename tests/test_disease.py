"""Tests for the honest disease classifier (curation lookup)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import disease  # noqa: E402


def test_human_genes_route_correctly():
    assert disease.classify("TARDBP")["category"] == "ALS"
    assert disease.classify("SNCA")["category"] == "PD"
    assert disease.classify("HTT")["category"] == "HD"
    assert disease.classify("MAPT")["category"] == "AD"
    assert disease.classify("ATXN2")["category"] == "SCA"


def test_fly_aliases_route_correctly():
    assert disease.classify("TBPH")["category"] == "ALS"          # fly TDP-43
    assert disease.classify("alpha-synuclein")["category"] == "PD"
    assert disease.classify("Pink1")["category"] == "PD"


def test_unknown_gene_is_other_not_faked():
    r = disease.classify("CG12345")
    assert r["category"] == "OTHER" and r["fly_ortholog"] is None


def test_no_ortholog_flagged_honestly():
    # HTT / SNCA / MAPT have no fly ortholog (human transgene used)
    assert disease.classify("HTT")["fly_ortholog"] is None
    assert disease.classify("SNCA")["fly_ortholog"] is None


def test_polyq_severity_note():
    assert "pathogenic" in disease.classify("HTT", "Q128")["note"]
    assert "sub-threshold" in disease.classify("HTT", "Q20")["note"]


def test_disease_hint_respected():
    assert disease.classify("unknowngene", disease_hint="ALS")["category"] == "ALS"
