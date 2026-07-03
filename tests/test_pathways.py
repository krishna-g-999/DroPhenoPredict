"""Tests for the GO-DAG ancestor propagation and pathway-activity scoring
(src/drophenopredict/pathways.py). Uses synthetic parents/expression data so
the tests don't depend on the real 32MB GO OBO / GAF files.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import pathways  # noqa: E402


def test_ancestors_walks_multi_level_dag():
    # child -> parent -> grandparent chain
    parents = {
        "GO:child": {"GO:parent"},
        "GO:parent": {"GO:grandparent"},
        "GO:grandparent": set(),
    }
    memo = {}
    anc = pathways._ancestors("GO:child", parents, memo)
    assert anc == {"GO:parent", "GO:grandparent"}


def test_ancestors_handles_diamond_inheritance():
    # child has two parents that share a common grandparent (diamond shape)
    parents = {
        "GO:child": {"GO:parentA", "GO:parentB"},
        "GO:parentA": {"GO:root"},
        "GO:parentB": {"GO:root"},
        "GO:root": set(),
    }
    memo = {}
    anc = pathways._ancestors("GO:child", parents, memo)
    assert anc == {"GO:parentA", "GO:parentB", "GO:root"}


def test_ancestors_root_term_has_no_ancestors():
    parents = {"GO:root": set()}
    assert pathways._ancestors("GO:root", {}, {}) == set()


def test_ancestors_memoization_is_consistent():
    parents = {"GO:a": {"GO:b"}, "GO:b": {"GO:c"}, "GO:c": set()}
    memo = {}
    first = pathways._ancestors("GO:a", parents, memo)
    second = pathways._ancestors("GO:a", parents, memo)   # cached path
    assert first == second == {"GO:b", "GO:c"}


def test_pathway_activity_zscores_and_averages_geneset():
    # 3 genes x 4 lines; geneset uses 2 of the 3 genes
    expr = pd.DataFrame(
        [[1.0, 2.0, 3.0, 4.0],
         [10.0, 20.0, 30.0, 40.0],
         [5.0, 5.0, 5.0, 5.0]],           # zero-variance gene, must not blow up (min_present + 1e-9 guard)
        index=["g1", "g2", "g3"], columns=["l1", "l2", "l3", "l4"],
    )
    genesets = {"GO:test": ("test pathway", frozenset(["g1", "g2"]))}
    pa = pathways.pathway_activity(expr, genesets, min_present=2)
    assert list(pa.columns) == ["GO:test"]
    assert list(pa.index) == ["l1", "l2", "l3", "l4"]
    # g1 and g2 are perfectly correlated (g2 = 10*g1) -> identical z-scores -> mean = same z-score
    g1_z = (expr.loc["g1"] - expr.loc["g1"].mean()) / expr.loc["g1"].std()
    assert np.allclose(pa["GO:test"].to_numpy(), g1_z.to_numpy(), atol=1e-6)


def test_pathway_activity_skips_geneset_below_min_present():
    expr = pd.DataFrame([[1.0, 2.0]], index=["g1"], columns=["l1", "l2"])
    genesets = {"GO:sparse": ("sparse", frozenset(["g1", "g_missing_a", "g_missing_b"]))}
    pa = pathways.pathway_activity(expr, genesets, min_present=2)
    assert pa.empty or "GO:sparse" not in pa.columns
