"""MULTI-MODAL test: does RNA-seq expression (a data type the original authors
did not have for this cross) add predictive power for AD eye-degeneration
severity, via our validated nested-CV adaptive framework (multimodal.py)?
Genuinely new analysis — the original GWAS paper never tested this.

RNA-seq is sex-specific (GSE117850); this F1 cross's sex composition is not
specified in the metadata we have, so we test BOTH the female and male
DGRP-parent expression panels as candidate predictors (the DGRP mothers
contribute the genetic background being mapped; parental expression in either
sex is a defensible proxy, tested honestly rather than assumed to be one).

Output: data/processed/yang2023_multimodal_results.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")
from drophenopredict import genotypes, expression, multimodal  # noqa: E402

PROC = Path("data/processed")
PRIMARY_TRAITS = ["nn.mean", "nn.sd", "cc.mean", "cc.sd", "mean.s.area"]


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def adaptive_permutation_test(kernels: dict, y: np.ndarray, observed_R2: float, *,
                              n_perm: int = 100, n_repeats: int = 3, seed: int = 7) -> float:
    """Correct null for the ADAPTIVE (nested-CV modality-selection) procedure:
    permute the line labels of y, then rerun the SAME nested_modality_cv
    procedure end-to-end (including its own inner modality selection) on the
    permuted data. Using a single fixed kernel's null here would be the wrong
    test -- it wouldn't reflect the adaptive selection's own null behaviour.
    Reduced n_repeats/n_perm vs. the headline estimate for tractability
    (compute cost scales with n_perm x n_repeats x n_splits x inner-CV)."""
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        y_perm = rng.permutation(y)
        res = multimodal.nested_modality_cv(kernels, y_perm, n_repeats=n_repeats, seed=seed + 1 + i)
        null[i] = res["adaptive_R2"]
    return float((1 + np.sum(null >= observed_R2)) / (1 + n_perm))


def main() -> None:
    blup = pd.read_csv(PROC / "yang2023_blup_lines.csv", index_col=0)
    blup.index = [f"line_{n}" for n in blup.index]
    Kg_full, g_iids, _ = genotypes.load_grm("data/processed/grm.npz")

    rows = []
    for sex_label, sex in [("female_expr", "female"), ("male_expr", "male")]:
        expr = expression.load_rnaseq(sex)
        common = [l for l in g_iids if l in set(expr.columns) and l in blup.index]
        print(f"\n[{sex_label}] common lines (GRM ∩ expr ∩ Yang2023): {len(common)}")
        gi = [g_iids.index(l) for l in common]
        Kg = Kg_full[np.ix_(gi, gi)]; Kg = Kg / np.mean(np.diag(Kg))
        Ke = multimodal.expression_kernel(expr, common)

        for trait in PRIMARY_TRAITS:
            y = blup.loc[common, trait].to_numpy(float)
            kernels = {"geno": Kg, "exp": Ke}
            fixed = multimodal.fixed_modality_r(kernels, y, n_repeats=10)
            nested = multimodal.nested_modality_cv(kernels, y, n_repeats=10)
            p_adaptive = adaptive_permutation_test(kernels, y, nested["adaptive_R2"])
            row = {"modality_set": sex_label, "trait": trait, "n": len(common),
                  "r_geno": round(fixed["geno"], 3), "r_exp": round(fixed["exp"], 3),
                  "r_combined": round(fixed["combined"], 3),
                  "nested_adaptive_r": round(nested["adaptive_r"], 3),
                  "perm_p_adaptive": round(p_adaptive, 4),
                  "selections": nested["selections"]}
            rows.append(row)
            print(f"  {trait:14s} geno={row['r_geno']:+.3f} exp={row['r_exp']:+.3f} "
                  f"comb={row['r_combined']:+.3f} nested={row['nested_adaptive_r']:+.3f} "
                  f"p={row['perm_p_adaptive']:.3f}  sel={row['selections']}")

    df = pd.DataFrame(rows)
    df["fdr"] = df.groupby("modality_set")["perm_p_adaptive"].transform(
        lambda p: bh(p.to_numpy()))
    df.to_csv(PROC / "yang2023_multimodal_results.csv", index=False)
    print("\n=== MULTI-MODAL (genotype vs RNA-seq expression): Yang 2023 AD degeneration ===")
    print(df.drop(columns=["selections"]).to_string(index=False))
    n_sig = int((df.fdr < 0.05).sum())
    print(f"\n{n_sig}/{len(df)} trait x modality-set combinations significant at FDR<0.05.")


if __name__ == "__main__":
    main()
