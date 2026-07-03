"""Map genomic positions (dm3/BDGP5) to genes via the DGRP gene-model GTF.

Source GTF: GSE67505_DGRP.gtf.gz (dm3 coordinates, matching DGRP2 variant ids).
Used to turn GWAS top SNPs into candidate genes (gene the SNP falls in, else the
nearest gene within a window).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

GTF = Path("data/raw/annotation/GSE67505_DGRP.gtf.gz")


def load_gene_spans(gtf: Path = GTF) -> pd.DataFrame:
    """gene_id, chrom, start, end, gene_name (one row per gene; min/max span)."""
    spans: dict[str, list] = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9:
                continue
            chrom, start, end, attr = f[0], int(f[3]), int(f[4]), f[8]
            if 'gene_id "' not in attr:
                continue
            gid = attr.split('gene_id "')[1].split('"')[0]
            name = attr.split('gene_name "')[1].split('"')[0] if 'gene_name "' in attr else gid
            if gid not in spans:
                spans[gid] = [chrom, start, end, name]
            else:
                spans[gid][1] = min(spans[gid][1], start)
                spans[gid][2] = max(spans[gid][2], end)
    return pd.DataFrame([(g, *v) for g, v in spans.items()],
                        columns=["gene_id", "chrom", "start", "end", "gene_name"])


def map_snps_to_genes(top: pd.DataFrame, genes: pd.DataFrame,
                      near_window: int = 5000) -> pd.DataFrame:
    """Aggregate top SNPs to candidate genes. `top` needs chrom,pos,p,variant_id.

    Returns genes with n_snps (supporting hits), best_p, min_dist (0 = inside).
    """
    by_chrom = {c: g.sort_values("start").reset_index(drop=True)
                for c, g in genes.groupby("chrom")}
    recs = []
    for _, snp in top.iterrows():
        g = by_chrom.get(snp["chrom"])
        if g is None:
            continue
        pos = snp["pos"]
        inside = g[(g["start"] <= pos) & (g["end"] >= pos)]
        if len(inside):
            for _, row in inside.iterrows():
                recs.append((row["gene_name"], row["gene_id"], 0, snp["p"]))
        else:
            d = pd.concat([(g["start"] - pos).abs(), (g["end"] - pos).abs()], axis=1).min(axis=1)
            j = int(d.idxmin())
            if d[j] <= near_window:
                recs.append((g.loc[j, "gene_name"], g.loc[j, "gene_id"], int(d[j]), snp["p"]))
    if not recs:
        return pd.DataFrame(columns=["gene_name", "gene_id", "n_snps", "best_p", "min_dist"])
    r = pd.DataFrame(recs, columns=["gene_name", "gene_id", "min_dist", "p"])
    out = (r.groupby(["gene_name", "gene_id"])
           .agg(n_snps=("p", "size"), best_p=("p", "min"), min_dist=("min_dist", "min"))
           .reset_index().sort_values("best_p"))
    return out
