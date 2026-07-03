"""Option 1 — transcriptome-based prediction vs genomic, head to head.

For each trait x sex, predict the phenotype with the SAME validated GBLUP engine
using three kernels on the SAME common lines and SAME CV folds:
  - GRM       : genomic relationship (genotypes)
  - TRM       : transcriptomic relationship (line-level expression, sex-matched)
  - GRM+TRM   : equal-weight multi-omic kernel
Reports genomic/transcriptomic/combined CV r, R2, and permutation p.

Output: data/processed/transcriptome_vs_genomic.csv, reports/figures/omics_comparison.png
Run: ./.venv/Scripts/python.exe scripts/06_transcriptome_prediction.py
"""
from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, genotypes, expression, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAITS = {"climbing": 1543, "startle": 2799, "starvation": 2798,
          "sleep": 1483, "lifespan": 1315}
MIN_LINES = 110
N_REPEATS = 10
N_PERM = 200
OUT = Path("data/processed")
FIG = Path("reports/figures")


def _norm(K):
    return K / np.mean(np.diag(K))


def _subset(K, iids, lines):
    idx = [iids.index(l) for l in lines]
    return K[np.ix_(idx, idx)]


def build_trms():
    trms = {}
    for sex in ("female", "male"):
        p = OUT / f"trm_{sex}.npz"
        if not p.exists():
            K, lines, ng = expression.build_trm(sex)
            expression.save_trm(K, lines, ng, p)
            print(f"  built TRM {sex}: {len(lines)} lines, {ng} genes")
        trms[sex] = expression.load_trm(p)
    return trms


def evaluate_kernel(K, y, label):
    h = gblup.reml(K, y)
    cv = gblup.cross_validate(K, y, n_repeats=N_REPEATS, seed=42)
    r2 = float(cv.r2_per_repeat.mean())
    perm = gblup.permutation_test(K, y, r2, n_perm=N_PERM)
    return {f"{label}_h2": round(h.h2, 3), f"{label}_r": round(float(cv.r_per_repeat.mean()), 3),
            f"{label}_R2": round(r2, 4), f"{label}_p": round(perm["p_value"], 4)}


def main() -> None:
    Kg, g_iids, gmeta = genotypes.load_grm("data/processed/grm.npz")
    trms = build_trms()
    sex_map = {"F": "female", "M": "male"}
    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        for sex_code in ("F", "M"):
            if sex_code not in set(vals["sex"].unique()):
                continue
            y = dgrpool.line_means(vals, sex=sex_code, harmonize=True)
            Kt, t_iids, tmeta = trms[sex_map[sex_code]]
            common = [l for l in g_iids if l in set(t_iids) and l in set(y.index)]
            if len(common) < MIN_LINES:
                print(f"skip {trait}/{sex_code}: {len(common)} common lines")
                continue
            yv = y.loc[common].to_numpy(float)
            Kg_s = _norm(_subset(Kg, g_iids, common))
            Kt_s = _norm(_subset(Kt, t_iids, common))
            Kc_s = 0.5 * (Kg_s + Kt_s)
            row = {"trait": trait, "sex": sex_code, "n_lines": len(common)}
            row.update(evaluate_kernel(Kg_s, yv, "geno"))
            row.update(evaluate_kernel(Kt_s, yv, "trans"))
            row.update(evaluate_kernel(Kc_s, yv, "comb"))
            rows.append(row)
            print(f"{trait:11s} {sex_code}  n={len(common):3d}  "
                  f"geno r={row['geno_r']:+.3f}(p{row['geno_p']:.3f})  "
                  f"trans r={row['trans_r']:+.3f}(p{row['trans_p']:.3f})  "
                  f"comb r={row['comb_r']:+.3f}(p{row['comb_p']:.3f})")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "transcriptome_vs_genomic.csv", index=False)
    print("\n=============== TRANSCRIPTOME vs GENOMIC ===============")
    print(df.to_string(index=False))
    print(f"\nSaved -> {OUT/'transcriptome_vs_genomic.csv'}")

    # figure: per trait/sex grouped bars of predictive r
    fig, ax = plt.subplots(figsize=(11, 5))
    labels = [f"{r.trait[:4]}.{r.sex}" for r in df.itertuples()]
    x = np.arange(len(df)); w = 0.27
    ax.bar(x - w, df["geno_r"], w, label="genomic (GRM)")
    ax.bar(x, df["trans_r"], w, label="transcriptomic (TRM)")
    ax.bar(x + w, df["comb_r"], w, label="combined")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Out-of-sample predictive r (5×10 CV)")
    ax.set_title("DGRP: transcriptome vs genome for phenotype prediction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "omics_comparison.png", dpi=130)
    print(f"Figure -> {FIG/'omics_comparison.png'}")


if __name__ == "__main__":
    main()
