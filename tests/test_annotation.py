"""Tests for SNP/gene mapping (src/drophenopredict/annotation.py) using a
synthetic gene-span table, so tests don't depend on the real GTF file.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import annotation  # noqa: E402

GENES = pd.DataFrame({
    "gene_id": ["FBgn0001", "FBgn0002", "FBgn0003"],
    "chrom": ["2L", "2L", "3R"],
    "start": [1000, 5000, 2000],
    "end": [1500, 5500, 2500],
    "gene_name": ["geneA", "geneB", "geneC"],
})


def test_snp_inside_gene_maps_directly():
    top = pd.DataFrame({"variant_id": ["2L_1200_SNP"], "chrom": ["2L"],
                        "pos": [1200], "p": [1e-8]})
    out = annotation.map_snps_to_genes(top, GENES)
    assert list(out["gene_name"]) == ["geneA"]
    assert out.iloc[0]["min_dist"] == 0


def test_snp_near_gene_maps_within_window():
    top = pd.DataFrame({"variant_id": ["2L_5600_SNP"], "chrom": ["2L"],
                        "pos": [5600], "p": [1e-8]})            # 100bp downstream of geneB
    out = annotation.map_snps_to_genes(top, GENES, near_window=5000)
    assert list(out["gene_name"]) == ["geneB"]
    assert out.iloc[0]["min_dist"] == 100


def test_snp_beyond_window_is_dropped():
    top = pd.DataFrame({"variant_id": ["2L_50000_SNP"], "chrom": ["2L"],
                        "pos": [50000], "p": [1e-8]})
    out = annotation.map_snps_to_genes(top, GENES, near_window=5000)
    assert out.empty


def test_snp_on_wrong_chromosome_never_matches():
    top = pd.DataFrame({"variant_id": ["4_1200_SNP"], "chrom": ["4"],
                        "pos": [1200], "p": [1e-8]})
    out = annotation.map_snps_to_genes(top, GENES)
    assert out.empty


def test_multiple_snps_aggregate_to_gene_with_best_p():
    top = pd.DataFrame({
        "variant_id": ["2L_1100_SNP", "2L_1200_SNP"],
        "chrom": ["2L", "2L"], "pos": [1100, 1200], "p": [1e-5, 1e-8],
    })
    out = annotation.map_snps_to_genes(top, GENES)
    row = out[out.gene_name == "geneA"].iloc[0]
    assert row["n_snps"] == 2
    assert row["best_p"] == 1e-8
