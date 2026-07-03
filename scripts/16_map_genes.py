"""Part 2b — map GWAS top hits to candidate genes (dm3/BDGP5 gene models).

Uses the DGRP gene-model GTF (GSE67505, dm3 coordinates matching the DGRP2
variants). For each predictable trait, maps top GWAS SNPs to the gene they fall
in (or the nearest gene within a window) and ranks candidate genes.

Output: data/processed/gwas_genes_<trait_sex>.csv + docs/FINDINGS_07_architecture.md table
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")

GTF = "data/raw/annotation/GSE67505_DGRP.gtf.gz"
PROC = Path("data/processed")
NEAR_WINDOW = 5000   # bp: assign to nearest gene within this distance if not inside one


def load_gene_spans() -> pd.DataFrame:
    spans = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9:
                continue
            chrom, start, end, attr = f[0], int(f[3]), int(f[4]), f[8]
            gid = attr.split('gene_id "')[1].split('"')[0]
            name = attr.split('gene_name "')[1].split('"')[0] if 'gene_name "' in attr else gid
            if gid not in spans:
                spans[gid] = [chrom, start, end, name]
            else:
                spans[gid][1] = min(spans[gid][1], start)
                spans[gid][2] = max(spans[gid][2], end)
    df = pd.DataFrame([(g, *v) for g, v in spans.items()],
                      columns=["gene_id", "chrom", "start", "end", "gene_name"])
    return df


def map_trait(top: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    by_chrom = {c: g.sort_values("start").reset_index(drop=True) for c, g in genes.groupby("chrom")}
    recs = []
    for _, snp in top.iterrows():
        c, pos = snp["chrom"], snp["pos"]
        g = by_chrom.get(c)
        if g is None:
            continue
        inside = g[(g["start"] <= pos) & (g["end"] >= pos)]
        if len(inside):
            for _, row in inside.iterrows():
                recs.append((row["gene_name"], row["gene_id"], 0, snp["p"], snp["variant_id"]))
        else:
            d_start = (g["start"] - pos).abs()
            d_end = (g["end"] - pos).abs()
            dist = np.minimum(d_start, d_end)
            j = dist.idxmin()
            if dist[j] <= NEAR_WINDOW:
                recs.append((g.loc[j, "gene_name"], g.loc[j, "gene_id"], int(dist[j]),
                             snp["p"], snp["variant_id"]))
    if not recs:
        return pd.DataFrame(columns=["gene_name", "gene_id", "n_snps", "best_p", "min_dist"])
    r = pd.DataFrame(recs, columns=["gene_name", "gene_id", "dist", "p", "variant_id"])
    agg = (r.groupby(["gene_name", "gene_id"])
           .agg(n_snps=("variant_id", "nunique"), best_p=("p", "min"),
                min_dist=("dist", "min")).reset_index()
           .sort_values("best_p"))
    return agg


def main() -> None:
    genes = load_gene_spans()
    print(f"Gene models: {len(genes)} genes (dm3)")
    summary = {}
    for f in sorted(PROC.glob("gwas_top_*.csv")):
        key = f.stem.replace("gwas_top_", "")
        top = pd.read_csv(f)
        agg = map_trait(top, genes)
        agg.to_csv(PROC / f"gwas_genes_{key}.csv", index=False)
        summary[key] = agg
        names = ", ".join(agg["gene_name"].head(8))
        print(f"\n{key}: top p={top['p'].min():.2e}  candidate genes (top SNPs): {names}")
    print("\nSaved per-trait candidate gene lists -> data/processed/gwas_genes_*.csv")


if __name__ == "__main__":
    main()
