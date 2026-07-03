"""Final RNA<->protein concordance analysis for the aSyn/PD arm (BCM-DMAS).

Uses the REAL, VERIFIED channel-to-genotype mapping (cross-referenced between
BCM-DMAS_assay_TMTquantitation_metadata.csv [specimenID->TMT channel] and
BCM-DMAS_biospecimen_metadata.csv [specimenID->genotype] — NOT the generic
"Control"/"Sample" suffixes in the protein file's own column headers, which
were checked and found to be Proteome-Discoverer reference-channel labels,
not real genotype assignments):

  CONTROL (Elav-Gal4, n=8): channels 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N
  ASYN    (PD/aSyn,   n=8): channels 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N

Steps:
  1. Load protein abundances (normalized), log2-transform, aggregate isoform
     rows to gene (FBgn) level by mean.
  2. Compute protein-level aSyn-vs-control log2FC (two-sample t-test, n=8 v 8).
  3. Compare against BOTH transcript estimates (script 39a): age-adjusted
     (well-powered, not day-matched) and day-10-only (day-matched, underpowered
     alone) — report both for honesty.
  4. Genome-wide RNA-protein concordance (Spearman r — calibration check).
  5. Focused test: do the aSyn QTT/convergent genes (Findings 09/13/15) show
     concordant direction at the protein level?

Output: data/processed/proteomics_asyn_log2fc.csv
        data/processed/proteomics_rna_concordance.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

PROC = Path("data/processed")
PROT_FILE = Path("data/raw/npjpd2025/proteomics/jshulman_flytmt16__Proteins.txt")

# Verified channel -> group (see docstring). Protein file uses "127C" (no
# underscore); TMT metadata used "127_C" — normalized here.
CONTROL_CHANNELS = ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N"]
ASYN_CHANNELS = ["130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N"]


def bh(p):
    p = np.asarray(p); n = p.size; o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    df = pd.read_csv(PROT_FILE, sep="\t")
    # locate the normalized-abundance columns for each channel
    col_of = {}
    for ch in CONTROL_CHANNELS + ASYN_CHANNELS:
        matches = [c for c in df.columns
                  if c.startswith("Abundances (Normalized): F1: ") and c.split(": ")[2].split(",")[0] == ch]
        if len(matches) != 1:
            sys.exit(f"Could not uniquely resolve channel {ch}: matches={matches}")
        col_of[ch] = matches[0]
    print("Resolved all 16 channel columns:")
    for ch, c in col_of.items():
        grp = "CONTROL" if ch in CONTROL_CHANNELS else "ASYN"
        print(f"  {ch:>5s} [{grp}] -> {c}")

    keep = ["Ensembl Gene ID", "Gene Symbol"] + list(col_of.values())
    prot = df[keep].dropna(subset=["Ensembl Gene ID"]).copy()
    prot = prot[prot["Ensembl Gene ID"].str.startswith("FBgn")]
    logab = np.log2(prot[list(col_of.values())].astype(float) + 1)
    logab["FBgn"] = prot["Ensembl Gene ID"].to_numpy()
    gene_level = logab.groupby("FBgn").mean()      # aggregate isoform rows
    print(f"\nProtein matrix: {gene_level.shape[0]} unique FBgn genes x {gene_level.shape[1]} channels")

    ctrl_cols = [col_of[ch] for ch in CONTROL_CHANNELS]
    asyn_cols = [col_of[ch] for ch in ASYN_CHANNELS]
    Yc = gene_level[ctrl_cols].to_numpy()
    Ya = gene_level[asyn_cols].to_numpy()
    log2fc = Ya.mean(1) - Yc.mean(1)
    tstat, pval = stats.ttest_ind(Ya, Yc, axis=1, equal_var=False)
    prot_de = pd.DataFrame({"log2FC": log2fc, "p": pval, "fdr": bh(pval)}, index=gene_level.index)
    prot_de.to_csv(PROC / "proteomics_asyn_log2fc.csv")
    print(f"Protein-level aSyn-vs-control DE (n=8 v 8): "
          f"{int((prot_de.fdr < 0.05).sum())} FDR<0.05 genes")

    # ---- concordance vs both transcript estimates ----
    rna_full = pd.read_csv(PROC / "asyn_log2fc_ageadjusted.csv", index_col=0)
    rna_d10 = pd.read_csv(PROC / "asyn_log2fc_day10only.csv", index_col=0)

    n_prot_nan = prot_de["log2FC"].isna().sum()
    print(f"\nProtein genes with NaN log2FC (missing in >=1 TMT channel, excluded from "
          f"correlation): {n_prot_nan}/{len(prot_de)}")

    results = []
    for label, rna in [("age_adjusted_transcript", rna_full), ("day10_only_transcript", rna_d10)]:
        common = prot_de.index.intersection(rna.index)
        a = prot_de.loc[common, "log2FC"]
        b = rna.loc[common, "log2FC"]
        valid = a.notna() & b.notna()
        r, p = stats.spearmanr(a[valid], b[valid])
        results.append({"comparison": label, "n_genes": int(valid.sum()),
                        "n_dropped_nan": int((~valid).sum()),
                        "spearman_r": round(r, 4), "p_value": p})
        print(f"\n[{label}] n={int(valid.sum())} genes ({int((~valid).sum())} dropped for "
              f"missing data), Spearman r={r:.3f} (p={p:.2e})")

    pd.DataFrame(results).to_csv(PROC / "proteomics_rna_concordance.csv", index=False)

    # ---- focused test: aSyn QTT genes concordant at protein level? ----
    qtt = pd.read_csv(PROC / "qtt_top_starvation_M.csv") if (PROC / "qtt_top_starvation_M.csv").exists() else None
    # (QTT tables were built per predictable DGRP trait, not per BCM-DMAS disease —
    #  the relevant focused check here is the BCM-DMAS aSyn-specific convergent/GSEA hits)
    conv = pd.read_csv(PROC / "shared_biology_convergent_genes.csv")
    conv_asyn = conv.dropna(subset=["aSyn"]) if "aSyn" in conv.columns else conv
    common_conv = prot_de.index.intersection(conv["Geneid"]).intersection(rna_full.index)
    if len(common_conv):
        sub_rna = rna_full.loc[common_conv, "log2FC"]
        sub_prot = prot_de.loc[common_conv, "log2FC"]
        valid_conv = sub_rna.notna() & sub_prot.notna()
        same_dir = (np.sign(sub_rna[valid_conv]) == np.sign(sub_prot[valid_conv])).sum()
        n_valid = int(valid_conv.sum())
        # binomial test: is same-direction rate > chance (0.5)?
        binom_p = stats.binomtest(int(same_dir), n_valid, p=0.5, alternative="greater").pvalue if n_valid else float("nan")
        print(f"\n[focused: BCM-DMAS convergent genes present in proteomics] "
              f"n={len(common_conv)} ({n_valid} with complete protein data), "
              f"same-direction (transcript vs protein): {same_dir}/{n_valid}  "
              f"(binomial test vs 50% chance: p={binom_p:.3f})")

    print(f"\nSaved -> proteomics_asyn_log2fc.csv, proteomics_rna_concordance.csv")


if __name__ == "__main__":
    main()
