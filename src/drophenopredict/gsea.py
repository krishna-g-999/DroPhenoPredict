"""Preranked GSEA (Subramanian et al. 2005, PNAS) — weighted Kolmogorov-Smirnov
running-sum enrichment statistic, with gene-permutation significance.

Implemented locally because the standard `gseapy` package could not be
installed in this environment (network/certificate failure, not a code issue).
Because we cannot cross-check against gseapy directly, this implementation is
instead VERIFIED against known ground-truth properties before use on real data
(see scripts/34_verify_gsea.py): a gene set concentrated at the top of a ranked
list must give a strong positive, significant enrichment score; a random gene
set must give an enrichment score statistically indistinguishable from the
permutation null (calibration check).

Algorithm (weighted, p=1, as in the original GSEA paper):
  Rank genes by a score (e.g. log2FC or t-statistic), descending.
  Walk down the ranked list; at each gene, increase a running-sum statistic by
  |score|/N_hit if the gene is IN the set, decrease by 1/N_miss if NOT.
  ES = the running-sum's max deviation from 0 (signed).
  Significance: permute gene labels (shuffle which genes belong to the set),
  recompute ES many times -> empirical p-value.
"""
from __future__ import annotations

import numpy as np


def enrichment_score(ranked_scores: np.ndarray, in_set: np.ndarray) -> float:
    """ranked_scores: genes already sorted descending by score.
    in_set: boolean array, same order, True if gene is a member of the gene set.
    Returns signed ES (max |running sum| deviation, sign preserved).
    """
    n = ranked_scores.size
    n_hit = in_set.sum()
    if n_hit == 0 or n_hit == n:
        return 0.0
    hit_weights = np.abs(ranked_scores) * in_set
    hit_step = hit_weights / hit_weights.sum()
    miss_step = (~in_set).astype(float) / (n - n_hit)
    running = np.cumsum(hit_step - miss_step)
    i_max = np.argmax(np.abs(running))
    return float(running[i_max])


def gsea_preranked(scores_by_gene: dict, gene_set: set, *, n_perm: int = 1000,
                   seed: int = 42) -> dict:
    """Preranked GSEA for one gene set against one ranked gene list.

    scores_by_gene: {gene_id: score} (e.g. log2FC). Sorted internally, descending.
    gene_set: set of gene_ids defining the pathway/set.
    Returns: {ES, n_hit, n_genes, p_value, null_mean, null_sd}
    """
    genes = np.array(list(scores_by_gene.keys()))
    scores = np.array([scores_by_gene[g] for g in genes], dtype=float)
    order = np.argsort(-scores)
    genes, scores = genes[order], scores[order]
    in_set = np.isin(genes, list(gene_set))
    n_hit = int(in_set.sum())
    if n_hit < 3:
        return {"ES": 0.0, "n_hit": n_hit, "n_genes": genes.size,
                "p_value": 1.0, "null_mean": 0.0, "null_sd": 0.0}

    es = enrichment_score(scores, in_set)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm_set = rng.permutation(in_set)      # shuffle set-membership labels
        null[i] = enrichment_score(scores, perm_set)
    # Two-sided: a gene set can be enriched at either end of the ranked list, so
    # significance must be judged against |null| >= |es| — NOT tested only in
    # the direction of the observed sign (that bug caps p at ~0.5 and inflates
    # the false-positive rate; caught by the calibration test in scripts/34).
    p = (1 + np.sum(np.abs(null) >= abs(es))) / (1 + n_perm)
    return {"ES": float(es), "n_hit": n_hit, "n_genes": int(genes.size),
            "p_value": float(p), "null_mean": float(null.mean()),
            "null_sd": float(null.std())}
