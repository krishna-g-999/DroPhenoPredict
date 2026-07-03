"""Part 2a — single-marker GWAS for the predictable traits (line means).

Streams all autosomal SNPs and computes marginal association (linear regression
of line-mean phenotype on SNP dosage) per trait x sex. Reports top hits with
variant id + position for gene mapping (script 16). Honest about power: with
n~190 DGRP lines, behavioral GWAS is underpowered for genome-wide significance;
top/suggestive hits are candidate loci, not definitive.

Output: data/processed/gwas_top_<trait_sex>.csv (top SNPs per predictable trait)
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
from drophenopredict import dgrpool, genotypes  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

BED = "data/raw/dgrp2/dgrp2.bed"
TRAIT_PID = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
TOP_K = 200
BLOCK = 50_000
PROC = Path("data/processed")


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    predictable = [k for k, m in man["traits"].items() if m["predictable"]]
    print("Predictable traits for GWAS:", predictable)

    bed = open_bed(BED)
    iids = [str(x) for x in bed.iid]
    pos = {l: i for i, l in enumerate(iids)}
    n_snp = bed.sid_count

    # phenotype vectors aligned to bed line order (NaN where missing)
    Y, masks = {}, {}
    for key in predictable:
        trait, sex = key.rsplit("_", 1)
        vals = dgrpool.fetch_phenotype_values(TRAIT_PID[trait])
        s = dgrpool.line_means(vals, sex=sex, harmonize=True)
        y = np.full(len(iids), np.nan)
        for l, v in s.items():
            if l in pos:
                y[pos[l]] = v
        Y[key] = y
        masks[key] = ~np.isnan(y)

    bim_ids = pd.read_csv(Path(BED).with_suffix(".bim"), sep="\t", header=None,
                          usecols=[1], names=["vid"], dtype=str)["vid"].to_numpy()
    auto = genotypes._autosomal_snp_mask(Path(BED).with_suffix(".bim"))

    # running top-K per trait: keep arrays of (t, idx)
    top = {k: [] for k in predictable}   # list of (abs_t, snp_index, beta, r)

    for start in range(0, n_snp, BLOCK):
        stop = min(start + BLOCK, n_snp)
        cols = np.where(auto[start:stop])[0]
        if cols.size == 0:
            continue
        G = bed.read(index=np.s_[:, start:stop], dtype="float32")[:, cols]  # 205 x b
        gidx = start + cols
        for key in predictable:
            m = masks[key]
            yo = Y[key][m]
            yc = yo - yo.mean()
            ss_y = np.sqrt(np.sum(yc * yc))
            Go = G[m]                                   # n_obs x b
            colmean = np.nanmean(Go, axis=0)
            inds = np.where(np.isnan(Go))
            Go[inds] = np.take(colmean, inds[1])
            Gc = Go - Go.mean(0)
            ss_g = np.sqrt(np.sum(Gc * Gc, axis=0))
            valid = ss_g > 1e-8
            r = np.zeros(Go.shape[1])
            r[valid] = (yc @ Gc[:, valid]) / (ss_y * ss_g[valid])
            r = np.clip(r, -0.999, 0.999)
            n = m.sum()
            t = r * np.sqrt((n - 2) / (1 - r * r))
            # keep top by |t|
            k = min(TOP_K, valid.sum())
            order = np.argsort(-np.abs(t))[:k]
            for o in order:
                top[key].append((abs(t[o]), int(gidx[o]), float(r[o]), float(t[o])))
        # periodically prune to TOP_K
        if (start // BLOCK) % 5 == 0:
            for key in predictable:
                top[key] = sorted(top[key], key=lambda z: -z[0])[:TOP_K]

    for key in predictable:
        rows = sorted(top[key], key=lambda z: -z[0])[:TOP_K]
        n = masks[key].sum()
        recs = []
        for _abs_t, idx, r, t in rows:
            p = 2 * stats.t.sf(abs(t), n - 2)
            vid = bim_ids[idx]
            arm, position, _ = vid.split("_")
            recs.append({"variant_id": vid, "chrom": arm, "pos": int(position),
                         "r": round(r, 3), "t": round(t, 2), "p": p})
        df = pd.DataFrame(recs)
        df.to_csv(PROC / f"gwas_top_{key}.csv", index=False)
        print(f"{key:13s} n={n:3d}  top p={df['p'].min():.2e}  "
              f"#p<1e-5={int((df['p']<1e-5).sum())}  -> gwas_top_{key}.csv")


if __name__ == "__main__":
    main()
