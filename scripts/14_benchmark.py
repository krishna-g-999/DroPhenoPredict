"""Part 1 — head-to-head benchmark: DroPhenoPredict vs baselines & literature.

Compares, per trait x sex, prediction accuracy (CV Pearson r):
  - Genotype-only GBLUP (= v1.0 method, = Ober 2012 approach)   [this work]
  - Expression-only (RNA-seq)                                   [this work]
  - Multi-modal adaptive (v1.1, nested CV)                      [this work]
  - Published anchors: Ober 2012 (genotype), Morgante 2020 (expression)

Published values are transcribed from the papers (provenance below) — NOT our
measurements — and used only to show our pipeline reproduces them.

Outputs: data/processed/benchmark_accuracy.csv, reports/figures/benchmark.png,
         docs/TOOL_COMPARISON.md (capability matrix)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")

PROC = Path("data/processed")
FIG = Path("reports/figures")

# Verified published accuracies (CV predictive correlation r). Provenance:
#  Ober et al. 2012, PLOS Genetics 8(5):e1002685 — DGRP GBLUP: startle 0.230±0.012,
#    starvation 0.239±0.008 (not sex-split in the cited figure).
#  Morgante et al. 2020, G3 10(12):4599 — DGRP expression vs genotype:
#    starvation expression r=0.28 (F)/0.38 (M); genotype r=0.07 (F)/0.15 (M).
PUBLISHED = {
    "Ober2012_geno": {"startle": 0.230, "starvation": 0.239},
    "Morgante2020_exp": {"starvation_F": 0.28, "starvation_M": 0.38},
    "Morgante2020_geno": {"starvation_F": 0.07, "starvation_M": 0.15},
}


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    rows = []
    for key, m in man["traits"].items():
        trait = m["trait"]
        rows.append({
            "trait_sex": key, "trait": trait, "sex": m["sex"], "n": m["n_lines"],
            "geno_GBLUP": m["r_geno"], "expression": m["r_exp"],
            "combined": m["r_combined"], "multimodal_adaptive": m["adaptive_nested_r"],
            "deployed": m["deployed_modality"], "predictable": m["predictable"],
            "Ober2012_geno": PUBLISHED["Ober2012_geno"].get(trait, np.nan),
            "Morgante2020_exp": PUBLISHED["Morgante2020_exp"].get(key, np.nan),
            "Morgante2020_geno": PUBLISHED["Morgante2020_geno"].get(key, np.nan),
        })
    df = pd.DataFrame(rows)
    df.to_csv(PROC / "benchmark_accuracy.csv", index=False)
    print("=============== ACCURACY BENCHMARK (CV r) ===============")
    print(df.drop(columns=["trait", "sex"]).to_string(index=False))

    # ---- reproduction check (our vs published, matched cells) ----
    checks = []
    for key, m in man["traits"].items():
        t = m["trait"]
        if t in PUBLISHED["Ober2012_geno"]:
            checks.append(("Ober2012 geno", f"{key}", m["r_geno"], PUBLISHED["Ober2012_geno"][t]))
        if key in PUBLISHED["Morgante2020_exp"]:
            checks.append(("Morgante2020 exp", key, m["r_exp"], PUBLISHED["Morgante2020_exp"][key]))
    print("\n--- literature reproduction (ours vs published) ---")
    for label, key, ours, pub in checks:
        print(f"  {label:18s} {key:13s} ours={ours:+.3f}  published={pub:+.3f}  Δ={ours-pub:+.3f}")

    # ---- figure: grouped bars per trait-sex ----
    fig, ax = plt.subplots(figsize=(12, 5.5))
    x = np.arange(len(df)); w = 0.22
    ax.bar(x - 1.5*w, df["geno_GBLUP"], w, label="Genotype GBLUP (v1.0 / Ober-style)")
    ax.bar(x - 0.5*w, df["expression"], w, label="Expression (RNA-seq)")
    ax.bar(x + 0.5*w, df["multimodal_adaptive"], w, label="Multi-modal adaptive (v1.1)", color="C2")
    # published markers
    for i, r in df.iterrows():
        if not np.isnan(r["Ober2012_geno"]):
            ax.plot(i - 1.5*w, r["Ober2012_geno"], "k_", ms=14, mew=2)
        if not np.isnan(r["Morgante2020_exp"]):
            ax.plot(i - 0.5*w, r["Morgante2020_exp"], "k_", ms=14, mew=2)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(df["trait_sex"], rotation=45, ha="right")
    ax.set_ylabel("Cross-validated predictive r")
    ax.set_title("DroPhenoPredict benchmark — multi-modal vs genotype/expression\n"
                 "(black ticks = published values: Ober 2012 geno, Morgante 2020 exp)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "benchmark.png", dpi=130)
    print(f"\nFigure -> {FIG/'benchmark.png'}")
    print(f"Table  -> {PROC/'benchmark_accuracy.csv'}")


if __name__ == "__main__":
    main()
