"""Cross-validate our local GSEA implementation (gsea.py, built and verified
in Findings 15 because the standard `gseapy` package could not be installed
in this sandbox at the time) against the REAL gseapy package, now that a
network change has resolved the installation issue.

Uses the EXACT SAME input for both: the aSyn model's age-adjusted log2FC
ranking (already computed, Findings 14/17) and the same GO-BP gene sets
(pathways.load_go_genesets — identical annotation source, identical ancestor
propagation) that produced the Findings 15 results being checked here.

This is a genuine, independent verification: if the two implementations
(different code, different exact algorithm details) agree on direction and
broadly on which pathways are top hits, that is strong evidence Findings 15's
GSEA convergence result (nucleolar stress down-regulation, HTT-200Q outlier)
is not an artifact of our custom implementation.

Output: data/processed/gseapy_crossvalidation.csv
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
from drophenopredict import pathways  # noqa: E402

PROC = Path("data/processed")


def main() -> None:
    import gseapy

    rna = pd.read_csv(PROC / "asyn_log2fc_ageadjusted.csv", index_col=0)
    rnk = rna["log2FC"].sort_values(ascending=False)
    print(f"Ranked list: {len(rnk)} genes (aSyn age-adjusted log2FC)")

    gsets = pathways.load_go_genesets(aspect="P", min_genes=10, max_genes=300)
    gene_sets = {f"{go}_{name[:40]}": list(genes) for go, (name, genes) in gsets.items()}
    print(f"Gene sets: {len(gene_sets)} (same source as Findings 15)")

    print("\nRunning gseapy.prerank (this may take a few minutes for ~2600 sets)...")
    pre = gseapy.prerank(rnk=rnk, gene_sets=gene_sets, min_size=10, max_size=300,
                        permutation_num=200, seed=42, threads=4, no_plot=True,
                        outdir=None)
    gp = pre.res2d.copy()
    gp["GO_id"] = gp["Term"].str.split("_").str[0]
    gp = gp.rename(columns={"NES": "gseapy_NES", "ES": "gseapy_ES", "FDR q-val": "gseapy_fdr"})
    gp["gseapy_ES"] = gp["gseapy_ES"].astype(float)
    gp["gseapy_NES"] = gp["gseapy_NES"].astype(float)
    gp["gseapy_fdr"] = gp["gseapy_fdr"].astype(float)
    print(f"gseapy returned {len(gp)} results")

    ours = pd.read_csv(PROC / "gsea_per_model.csv")
    ours_asyn = ours[ours.model == "aSyn"][["GO_id", "pathway", "ES", "p_value", "fdr"]]
    ours_asyn = ours_asyn.rename(columns={"ES": "our_ES", "p_value": "our_p", "fdr": "our_fdr"})

    merged = ours_asyn.merge(gp[["GO_id", "gseapy_ES", "gseapy_NES", "gseapy_fdr"]],
                             on="GO_id", how="inner")
    merged.to_csv(PROC / "gseapy_crossvalidation.csv", index=False)
    print(f"\nMatched {len(merged)} pathways between our engine and gseapy")

    from scipy import stats
    r_es, p_es = stats.pearsonr(merged["our_ES"], merged["gseapy_ES"])
    r_nes, p_nes = stats.pearsonr(merged["our_ES"], merged["gseapy_NES"])
    same_sign = (np.sign(merged["our_ES"]) == np.sign(merged["gseapy_ES"])).mean()
    print(f"Our ES vs gseapy ES:  Pearson r={r_es:.3f} (p={p_es:.2e})")
    print(f"Our ES vs gseapy NES: Pearson r={r_nes:.3f} (p={p_nes:.2e})")
    print(f"Same-direction agreement: {same_sign*100:.1f}% of {len(merged)} pathways")

    # do our top-10 significant hits from Findings 15 show up as gseapy-significant too?
    our_top = ours_asyn.sort_values("our_p").head(10)
    check = our_top.merge(gp[["GO_id", "gseapy_fdr", "gseapy_NES"]], on="GO_id", how="left")
    print("\nOur top-10 aSyn hits (Findings 15) vs gseapy:")
    print(check[["pathway", "our_p", "our_fdr", "gseapy_fdr", "gseapy_NES"]].to_string(index=False))


if __name__ == "__main__":
    main()
