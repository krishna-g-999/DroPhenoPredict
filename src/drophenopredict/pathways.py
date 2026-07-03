"""GO biological-process gene sets for Drosophila + pathway-activity scoring.

Builds GO-term -> gene (FBgn) sets from the FlyBase GAF with proper ANCESTOR
PROPAGATION over the GO DAG (a gene annotated to a specific term also belongs to
all its is_a/part_of ancestors), so broad terms like 'lipid metabolic process'
get their full gene complement. Pathway activity per DGRP line = mean z-scored
expression of the set's genes (from RNA-seq). Used to map which pathways drive
which traits — formalizing the lipid->starvation finding across all processes.

Sources (public): FlyBase GAF (current.geneontology.org/annotations/fb.gaf.gz),
GO basic ontology (go-basic.obo).
"""
from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

GAF = Path("data/raw/go/fb.gaf.gz")
OBO = Path("data/raw/go/go-basic.obo")


def _parse_obo(obo: Path = OBO):
    """Return (names: id->name, parents: id->set(parent ids via is_a/part_of), aspect via namespace)."""
    names, parents, namespace = {}, {}, {}
    cur = None
    for line in open(obo, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"id": None, "parents": set()}
        elif line == "[Typedef]":
            cur = None
        elif cur is not None:
            if line.startswith("id: GO:"):
                cur["id"] = line[4:]
            elif line.startswith("name:"):
                names[cur["id"]] = line[6:]
            elif line.startswith("namespace:"):
                namespace[cur["id"]] = line[11:]
            elif line.startswith("is_a: GO:"):
                cur["parents"].add(line[6:].split(" ! ")[0])
            elif line.startswith("relationship: part_of GO:"):
                cur["parents"].add(line.split("part_of ")[1].split(" ! ")[0])
            if cur["id"]:
                parents[cur["id"]] = cur["parents"]
    return names, parents, namespace


def _ancestors(term, parents, memo):
    if term in memo:
        return memo[term]
    anc = set()
    for p in parents.get(term, ()):
        anc.add(p)
        anc |= _ancestors(p, parents, memo)
    memo[term] = anc
    return anc


def load_go_genesets(aspect: str = "P", min_genes: int = 10, max_genes: int = 300):
    """{GO_id: (name, frozenset(FBgn))} for the chosen aspect, ancestor-propagated."""
    names, parents, namespace = _parse_obo()
    asp_ns = {"P": "biological_process", "F": "molecular_function", "C": "cellular_component"}[aspect]
    memo: dict = {}
    sets: dict[str, set] = {}
    with gzip.open(GAF, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[8] != aspect:
                continue
            if "NOT" in f[3]:
                continue
            gene, go = f[1], f[4]            # FBgn, GO id
            for t in (go, *_ancestors(go, parents, memo)):
                if namespace.get(t) == asp_ns:
                    sets.setdefault(t, set()).add(gene)
    out = {t: (names.get(t, t), frozenset(g)) for t, g in sets.items()
           if min_genes <= len(g) <= max_genes}
    logger.info("GO %s gene sets (size %d-%d): %d", aspect, min_genes, max_genes, len(out))
    return out


def pathway_activity(expr_df: pd.DataFrame, genesets: dict, min_present: int = 5) -> pd.DataFrame:
    """lines x pathways: mean z-scored (across lines) expression of set genes.

    expr_df: genes(FBgn) x lines. Only genes present in expr are used.
    """
    Z = expr_df.sub(expr_df.mean(axis=1), axis=0).div(expr_df.std(axis=1) + 1e-9, axis=0)  # gene z across lines
    present = set(expr_df.index)
    cols = {}
    for go, (name, genes) in genesets.items():
        g = [x for x in genes if x in present]
        if len(g) >= min_present:
            cols[go] = Z.loc[g].mean(axis=0)
    pa = pd.DataFrame(cols).T          # pathways x lines
    return pa.T                         # lines x pathways
