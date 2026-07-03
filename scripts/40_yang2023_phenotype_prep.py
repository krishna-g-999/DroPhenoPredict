"""Yang et al. 2023 (G3, jkad132) DGRP x AD eye-degeneration phenotype prep.

Replicates the ORIGINAL AUTHORS' OWN mixed-model design (verified from their
public analysis script, 02_AD_fly_eye.rmd):
    trait ~ group + (1|line), REML
extracting the BLUP (best linear unbiased predictor) random effect per DGRP
line as the batch-adjusted genetic phenotype, and the broad-sense heritability
proxy h2 = sigma_line^2 / (sigma_line^2 + sigma_resid^2) from the fitted
variance components — same formula the original script uses.

Design decisions (see docs/EDITORIAL_RISK_ANALYSIS.md for the full reasoning):
  - class=='line' rows ONLY (excludes the '32C'/'441C' internal QC control
    crosses BEFORE model fitting, not after — avoids any chance of those
    contaminating a later join with real DGRP line 32/441 genotypes).
  - Primary traits = the first 5 columns (nn.mean, nn.sd, cc.mean, cc.sd,
    mean.s.area) — the ORIGINAL AUTHORS' OWN "primary" designation (visible in
    their RMD's plotting code), not a post-hoc choice by us. Avoids fishing
    across all 16 correlated features.
  - All class=='line' rows are day==28 (verified below, not assumed) — no
    day confound in the mapping population.

Output: data/processed/yang2023_blup_lines.csv   (line x trait BLUP values)
        data/processed/yang2023_heritability.csv  (per-trait h2 + variance components)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

RDATA = Path("data/raw/yang2023_github/DGRP_rawUntransformed_16traits.rdata")
PROC = Path("data/processed")
PRIMARY_TRAITS = ["nn.mean", "nn.sd", "cc.mean", "cc.sd", "mean.s.area"]


def main() -> None:
    import pyreadr
    from statsmodels.regression.mixed_linear_model import MixedLM

    result = pyreadr.read_r(str(RDATA))
    df, info = result["df"], result["df.info"]
    print(f"Raw data: {df.shape[0]} images x {df.shape[1]} traits")
    print(f"Classes: {info['class'].value_counts().to_dict()}")

    # ---- verify day==28 for the mapping population before assuming it ----
    line_mask = info["class"] == "line"
    days = info.loc[line_mask, "day"].dropna().unique()
    print(f"Unique 'day' values for class=='line': {days}")
    if not (len(days) == 1 and days[0] == 28):
        sys.exit(f"UNEXPECTED: class=='line' rows are not all day 28 ({days}) — "
                 f"the day-confound assumption in the docstring is wrong; stopping "
                 f"rather than proceeding on a false premise.")

    sub_df = df[line_mask].reset_index(drop=True)
    sub_info = info[line_mask].reset_index(drop=True)
    n_lines = sub_info["line"].nunique()
    print(f"class=='line' only: {len(sub_df)} images, {n_lines} unique DGRP lines, "
          f"{sub_info['group'].nunique()} groups")

    group = pd.get_dummies(sub_info["group"], prefix="grp", drop_first=True).astype(float)
    exog = np.column_stack([np.ones(len(sub_info)), group.to_numpy()])
    line_ids = sub_info["line"].astype(str).to_numpy()

    blup_rows = {}
    herit_rows = []
    for trait in PRIMARY_TRAITS:
        y = sub_df[trait].to_numpy(dtype=float)
        valid = np.isfinite(y)
        model = MixedLM(y[valid], exog[valid], groups=line_ids[valid])
        fit = model.fit(reml=True)
        re = fit.random_effects   # dict: line_id -> Series with the random intercept
        blup = {k: v.iloc[0] for k, v in re.items()}
        blup_rows[trait] = pd.Series(blup)

        cov_re = np.asarray(fit.cov_re)
        sigma_line = float(cov_re[0, 0])
        sigma_resid = float(fit.scale)
        h2 = sigma_line / (sigma_line + sigma_resid)
        herit_rows.append({"trait": trait, "n_images": int(valid.sum()),
                           "n_lines": len(blup), "sigma_line": sigma_line,
                           "sigma_resid": sigma_resid, "h2_broad_sense": round(h2, 4)})
        print(f"  {trait:14s} n_img={int(valid.sum()):4d} n_lines={len(blup):3d} "
              f"h2={h2:.3f}  (sigma_line={sigma_line:.4f}, sigma_resid={sigma_resid:.4f})")

    blup_df = pd.DataFrame(blup_rows)
    blup_df.index.name = "dgrp_line_num"
    blup_df.to_csv(PROC / "yang2023_blup_lines.csv")
    pd.DataFrame(herit_rows).to_csv(PROC / "yang2023_heritability.csv", index=False)

    # ---- self-check: BLUP line ranking should correlate strongly with the naive
    # (unadjusted) per-line raw mean -- if it doesn't, something is wrong with the
    # mixed-model fit, not a subtle real effect ----
    from scipy import stats
    naive_mean = sub_df.groupby(sub_info["line"])[PRIMARY_TRAITS[0]].mean()
    common = blup_df.index.intersection(naive_mean.index.astype(str))
    r, p = stats.pearsonr(blup_df.loc[common, PRIMARY_TRAITS[0]],
                          naive_mean.loc[common.astype(naive_mean.index.dtype)])
    print(f"\n[self-check] BLUP vs naive raw line-mean ({PRIMARY_TRAITS[0]}): "
          f"Pearson r={r:.3f} (expect strongly positive, e.g. >0.8, if the mixed "
          f"model is behaving sensibly)")
    if r < 0.5:
        print("  WARNING: low correlation between BLUP and naive mean — investigate "
              "before trusting downstream results.")

    print(f"\nSaved -> {PROC}/yang2023_blup_lines.csv, yang2023_heritability.csv")


if __name__ == "__main__":
    main()
