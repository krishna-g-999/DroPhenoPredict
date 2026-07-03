"""Option 3 — DGRP predictability map: the validated GBLUP baseline across all
five high-coverage traits x sex.

For each (trait, sex) cell with enough lines: REML genomic h2, 5x10 CV
(R2/r/MAE vs null), permutation p, split-conformal 90% coverage. One harmonized
source (DGRPool API) + one verified engine (gblup.py). No fabrication.

Output: data/processed/predictability_map.csv, reports/figures/predictability_map.png
Run: ./.venv/Scripts/python.exe scripts/05_predictability_map.py
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
from drophenopredict import dgrpool, genotypes, gblup  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

# Best-covered phenotype per family (from Findings 01, verified DGRPool ids)
TRAITS = {
    "climbing":   1543,   # mn_NegGeotaxis_CTRL
    "startle":    2799,   # StartleRes
    "starvation": 2798,   # StarvationRes
    "sleep":      1483,   # mn_SleepLength_Night
    "lifespan":   1315,   # mn_Longevity
}
MIN_LINES = 120
N_REPEATS = 10
N_PERM = 200
OUT = Path("data/processed")
FIG = Path("reports/figures")


def evaluate(K, iids, trait, pid, sex, y) -> dict | None:
    Ksub, yv, lines = gblup.align(K, iids, y)
    if len(lines) < MIN_LINES:
        print(f"  skip {trait}/{sex}: only {len(lines)} lines")
        return None
    h = gblup.reml(Ksub, yv)
    cv = gblup.cross_validate(Ksub, yv, n_repeats=N_REPEATS, seed=42)
    r2 = float(cv.r2_per_repeat.mean())
    perm = gblup.permutation_test(Ksub, yv, r2, n_perm=N_PERM)
    conf = gblup.conformal_coverage(Ksub, yv, alpha=0.10, n_repeats=N_REPEATS)
    row = {
        "trait": trait, "phenotype_id": pid, "sex": sex, "n_lines": len(lines),
        "genomic_h2": round(h.h2, 3),
        "cv_R2": round(r2, 4), "cv_R2_sd": round(float(cv.r2_per_repeat.std()), 4),
        "pearson_r": round(float(cv.r_per_repeat.mean()), 3),
        "MAE": round(float(cv.mae_per_repeat.mean()), 3),
        "null_MAE": round(float(cv.null_mae_per_repeat.mean()), 3),
        "perm_p": round(perm["p_value"], 4),
        "cover90": round(conf["empirical_coverage"], 3),
    }
    print(f"  {trait:11s} {sex:2s} n={row['n_lines']:3d}  h2={row['genomic_h2']:.2f}  "
          f"CV_R2={row['cv_R2']:+.3f}  r={row['pearson_r']:.3f}  p={row['perm_p']:.3f}")
    return row


def main() -> None:
    K, iids, meta = genotypes.load_grm("data/processed/grm.npz")
    print(f"GRM: {meta['n_snps_used']} SNPs, {len(iids)} lines\n")
    rows = []
    for trait, pid in TRAITS.items():
        vals = dgrpool.fetch_phenotype_values(pid)
        sexes = sorted(s for s in vals["sex"].unique() if s in ("F", "M"))
        print(f"{trait} (id {pid}): sexes present = {sexes or 'none-tagged'}")
        targets = [(s, dgrpool.line_means(vals, sex=s, harmonize=True)) for s in sexes] \
            if sexes else [("all", dgrpool.line_means(vals, harmonize=True))]
        for sex, y in targets:
            r = evaluate(K, iids, trait, pid, sex, y)
            if r:
                rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "predictability_map.csv", index=False)
    print("\n================= DGRP PREDICTABILITY MAP =================")
    print(df.to_string(index=False))
    print(f"\nSaved -> {OUT/'predictability_map.csv'}")

    # figure: genomic h2 vs out-of-sample predictive r, per trait/sex
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for _, r in df.iterrows():
        ax.scatter(r["genomic_h2"], r["pearson_r"], s=90,
                   edgecolor="k", zorder=3)
        ax.annotate(f"{r['trait'][:4]}.{r['sex']}", (r["genomic_h2"], r["pearson_r"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=8)
    lim = max(0.05, df["genomic_h2"].max() + 0.05)
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5, label="r = h² (unreachable ceiling)")
    ax.axhline(0, color="gray", lw=0.7)
    ax.set_xlabel("Genomic heritability h² (REML, in-sample)")
    ax.set_ylabel("Out-of-sample predictive r (5×10 CV)")
    ax.set_title("DGRP predictability map: heritability ≠ predictability\n"
                 "(unrelated panel → large gap between h² and CV r)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "predictability_map.png", dpi=130)
    print(f"Figure -> {FIG/'predictability_map.png'}")


if __name__ == "__main__":
    main()
