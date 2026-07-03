"""Part 2 (refined) — inversion/Wolbachia-ADJUSTED GWAS for predictable traits.

Standard DGRP practice: test each SNP for association with the trait AFTER
adjusting for major polymorphic inversions and Wolbachia status (large LD blocks /
endosymbiont that confound naive association). Implemented by Frisch-Waugh-Lovell
residualization: regress both phenotype and each SNP on the covariate design
matrix X = [1, inversions, Wolbachia], then correlate the residuals (partial r).

Outputs per predictable trait:
  data/processed/gwas_adj_top_<key>.csv     (adjusted top SNPs)
  data/processed/gwas_adj_genes_<key>.csv   (adjusted candidate genes)
  data/processed/gwas_adjustment_summary.csv (marginal vs adjusted overlap)
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from bed_reader import open_bed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, covariates, annotation  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

BED = "data/raw/dgrp2/dgrp2.bed"
TRAIT_PID = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
TOP_K, BLOCK = 200, 50_000
PROC = Path("data/processed")


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    predictable = [k for k, m in man["traits"].items() if m["predictable"]]
    cov = covariates.load_covariates()
    print("Adjusted GWAS for:", predictable, "| covariates:", list(cov.columns))

    bed = open_bed(BED)
    iids = [str(x) for x in bed.iid]
    pos = {l: i for i, l in enumerate(iids)}
    n_snp = bed.sid_count

    Y, masks, Bmats, Xs, yres = {}, {}, {}, {}, {}
    for key in predictable:
        trait, sex = key.rsplit("_", 1)
        s = dgrpool.line_means(dgrpool.fetch_phenotype_values(TRAIT_PID[trait]),
                               sex=sex, harmonize=True)
        y = np.full(len(iids), np.nan)
        for l, v in s.items():
            if l in pos:
                y[pos[l]] = v
        m = ~np.isnan(y)
        lines = [iids[i] for i in np.where(m)[0]]
        X = covariates.design_matrix(cov, lines)          # n_obs x p (incl intercept)
        B = np.linalg.solve(X.T @ X, X.T)                 # p x n_obs  (=(X'X)^-1 X')
        yo = y[m]
        yr = yo - X @ (B @ yo)                            # residualized phenotype
        Y[key], masks[key], Xs[key], Bmats[key], yres[key] = y, m, X, B, yr

    bim_ids = pd.read_csv(Path(BED).with_suffix(".bim"), sep="\t", header=None,
                          usecols=[1], names=["vid"], dtype=str)["vid"].to_numpy()
    auto = genotypes._autosomal_snp_mask(Path(BED).with_suffix(".bim"))

    top = {k: [] for k in predictable}
    for start in range(0, n_snp, BLOCK):
        stop = min(start + BLOCK, n_snp)
        cols = np.where(auto[start:stop])[0]
        if cols.size == 0:
            continue
        G = bed.read(index=np.s_[:, start:stop], dtype="float32")[:, cols]
        gidx = start + cols
        for key in predictable:
            m, X, B, yr = masks[key], Xs[key], Bmats[key], yres[key]
            Go = G[m]
            cm = np.nanmean(Go, axis=0)
            nan = np.where(np.isnan(Go)); Go[nan] = np.take(cm, nan[1])
            Gr = Go - X @ (B @ Go)                        # residualize SNPs on covariates
            ss_g = np.sqrt(np.sum(Gr * Gr, axis=0))
            ss_y = np.sqrt(np.sum(yr * yr))
            valid = ss_g > 1e-8
            r = np.zeros(Go.shape[1])
            r[valid] = (yr @ Gr[:, valid]) / (ss_y * ss_g[valid])
            r = np.clip(r, -0.999, 0.999)
            dfree = m.sum() - X.shape[1] - 1
            t = r * np.sqrt(dfree / (1 - r * r))
            k = min(TOP_K, int(valid.sum()))
            for o in np.argsort(-np.abs(t))[:k]:
                top[key].append((abs(t[o]), int(gidx[o]), float(r[o]), float(t[o]), dfree))
        if (start // BLOCK) % 5 == 0:
            for key in predictable:
                top[key] = sorted(top[key], key=lambda z: -z[0])[:TOP_K]

    genes = annotation.load_gene_spans()
    summary = []
    for key in predictable:
        rows = sorted(top[key], key=lambda z: -z[0])[:TOP_K]
        recs = []
        for _a, idx, r, t, dfree in rows:
            p = 2 * stats.t.sf(abs(t), dfree)
            arm, position, _ = bim_ids[idx].split("_")
            recs.append({"variant_id": bim_ids[idx], "chrom": arm, "pos": int(position),
                         "r": round(r, 3), "t": round(t, 2), "p": p})
        adj = pd.DataFrame(recs)
        adj.to_csv(PROC / f"gwas_adj_top_{key}.csv", index=False)
        gmap = annotation.map_snps_to_genes(adj, genes)
        gmap.to_csv(PROC / f"gwas_adj_genes_{key}.csv", index=False)

        marg = pd.read_csv(PROC / f"gwas_top_{key}.csv")
        overlap = len(set(adj["variant_id"]) & set(marg["variant_id"]))
        summary.append({"trait_sex": key, "adj_top_p": adj["p"].min(),
                        "adj_p<1e-5": int((adj["p"] < 1e-5).sum()),
                        "marginal_p<1e-5": int((marg["p"] < 1e-5).sum()),
                        "top200_overlap": overlap,
                        "top_genes": ", ".join(gmap["gene_name"].head(5))})
        print(f"{key:13s} adj top p={adj['p'].min():.2e}  adj #p<1e-5={int((adj['p']<1e-5).sum()):3d}"
              f"  overlap(marg)={overlap}/200  genes: {summary[-1]['top_genes']}")

    pd.DataFrame(summary).to_csv(PROC / "gwas_adjustment_summary.csv", index=False)
    print(f"\nSaved adjusted GWAS + genes + summary -> {PROC}/gwas_adj_*")


if __name__ == "__main__":
    main()
