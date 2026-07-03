"""Modern technique 1 — functionally-stratified genomic prediction (MultiBLUP-style).

Partitions autosomal SNPs into GENIC (inside a gene span) vs INTERGENIC relationship
matrices and fits a two-component REML model (genic + intergenic + residual). If
trait-relevant variants concentrate in genic regions, the model upweights them and
can beat a single GRM (Speed & Balding LDAK; Edwards 2016 MultiBLUP). We TEST this
honestly on the predictable traits — gains in an unrelated panel are not assumed.

Builds & caches: data/processed/grm_genic.npz, grm_intergenic.npz
Output: data/processed/stratified_prediction.csv
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, gblup, annotation  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

BED = "data/raw/dgrp2/dgrp2.bed"
TRAIT_PID = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
BLOCK = 50_000
PROC = Path("data/processed")


def genic_mask_for_snps(bim_path: Path) -> np.ndarray:
    """Boolean over ALL variants: True if the SNP position falls inside a gene."""
    bim = pd.read_csv(bim_path, sep="\t", header=None, usecols=[1], names=["vid"], dtype=str)
    vid = bim["vid"].to_numpy()
    arm = np.array([v.split("_")[0] for v in vid])
    posn = np.array([int(v.split("_")[1]) if v.count("_") >= 2 else -1 for v in vid])
    genes = annotation.load_gene_spans()
    out = np.zeros(len(vid), dtype=bool)
    for a, g in genes.groupby("chrom"):
        starts = g["start"].to_numpy(); ends = g["end"].to_numpy()
        order = np.argsort(starts); starts, ends = starts[order], ends[order]
        # merge intervals
        ms, me = [], []
        for s, e in zip(starts, ends):
            if ms and s <= me[-1]:
                me[-1] = max(me[-1], e)
            else:
                ms.append(s); me.append(e)
        ms, me = np.array(ms), np.array(me)
        idx = np.where(arm == a)[0]
        if idx.size == 0:
            continue
        p = posn[idx]
        j = np.searchsorted(ms, p, side="right") - 1
        inside = (j >= 0) & (p <= np.where(j >= 0, me[j], -1))
        out[idx] = inside
    return out


def build_two_grms(genic_mask):
    bed = open_bed(BED)
    n_line, n_snp = bed.iid_count, bed.sid_count
    iids = [str(x) for x in bed.iid]
    auto = genotypes._autosomal_snp_mask(Path(BED).with_suffix(".bim"))
    Kg = np.zeros((n_line, n_line)); Ki = np.zeros((n_line, n_line))
    mg = mi = 0
    for start in range(0, n_snp, BLOCK):
        stop = min(start + BLOCK, n_snp)
        cols = np.where(auto[start:stop])[0]
        if cols.size == 0:
            continue
        gabs = start + cols
        G = bed.read(index=np.s_[:, start:stop], dtype="float32")[:, cols]
        p = np.nanmean(G, axis=0) / 2.0
        maf = np.minimum(p, 1 - p)
        cr = (~np.isnan(G)).mean(axis=0)
        keep = (cr >= 0.8) & (maf >= 0.05) & (p > 0) & (p < 1)
        if not keep.any():
            continue
        Gk = G[:, keep]; pk = p[keep]
        Z = Gk - 2 * pk; Z[np.isnan(Gk)] = 0.0; Z /= np.sqrt(2 * pk * (1 - pk))
        is_genic = genic_mask[gabs[keep]]
        if is_genic.any():
            Zg = Z[:, is_genic]; Kg += Zg @ Zg.T; mg += int(is_genic.sum())
        if (~is_genic).any():
            Zi = Z[:, ~is_genic]; Ki += Zi @ Zi.T; mi += int((~is_genic).sum())
    Kg /= max(mg, 1); Ki /= max(mi, 1)
    np.savez_compressed(PROC / "grm_genic.npz", K=Kg, iids=np.array(iids), n=mg)
    np.savez_compressed(PROC / "grm_intergenic.npz", K=Ki, iids=np.array(iids), n=mi)
    print(f"genic SNPs={mg}, intergenic SNPs={mi}")
    return Kg, Ki, iids


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    predictable = [k for k, m in man["traits"].items() if m["predictable"]]
    Kfull, fiids, _ = genotypes.load_grm("data/processed/grm.npz")

    if (PROC / "grm_genic.npz").exists():
        dg = np.load(PROC / "grm_genic.npz", allow_pickle=True)
        di = np.load(PROC / "grm_intergenic.npz", allow_pickle=True)
        Kg, Ki, iids = dg["K"], di["K"], [str(x) for x in dg["iids"]]
        print(f"loaded cached stratified GRMs (genic n={int(dg['n'])}, intergenic n={int(di['n'])})")
    else:
        gmask = genic_mask_for_snps(Path(BED).with_suffix(".bim"))
        print(f"genic mask: {gmask.sum()} of {gmask.size} variants inside genes")
        Kg, Ki, iids = build_two_grms(gmask)

    pos = {l: i for i, l in enumerate(iids)}
    rows = []
    for key in predictable:
        trait, sex = key.rsplit("_", 1)
        s = dgrpool.line_means(dgrpool.fetch_phenotype_values(TRAIT_PID[trait]),
                               sex=sex, harmonize=True)
        common = [l for l in iids if l in set(s.index)]
        idx = [pos[l] for l in common]
        yv = s.loc[common].to_numpy(float)
        Kg_s = Kg[np.ix_(idx, idx)] / np.mean(np.diag(Kg[np.ix_(idx, idx)]))
        Ki_s = Ki[np.ix_(idx, idx)] / np.mean(np.diag(Ki[np.ix_(idx, idx)]))
        fidx = [fiids.index(l) for l in common]
        Kf = Kfull[np.ix_(fidx, fidx)] / np.mean(np.diag(Kfull[np.ix_(fidx, fidx)]))
        r_single = gblup.cross_validate(Kf, yv, n_repeats=10, seed=42).r_per_repeat.mean()
        r_strat = gblup.cross_validate_multi([Kg_s, Ki_s], yv, n_repeats=5, seed=42)["r_mean"]
        rows.append({"trait_sex": key, "n": len(common),
                     "r_single_GRM": round(float(r_single), 3),
                     "r_stratified": round(float(r_strat), 3),
                     "delta": round(float(r_strat - r_single), 3)})
        print(f"{key:13s} single={rows[-1]['r_single_GRM']:+.3f} "
              f"stratified={rows[-1]['r_stratified']:+.3f} d={rows[-1]['delta']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(PROC / "stratified_prediction.csv", index=False)
    print("\n=== FUNCTIONALLY-STRATIFIED (genic|intergenic) vs SINGLE GRM ===")
    print(df.to_string(index=False))
    print(f"mean Δr = {df['delta'].mean():+.3f}")


if __name__ == "__main__":
    main()
