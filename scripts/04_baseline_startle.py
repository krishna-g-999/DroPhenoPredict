"""Single-trait GBLUP baseline: startle-induced locomotion (female & male).

Complete, validated genomic-prediction baseline:
  REML genomic heritability  ->  repeated k-fold CV (R2/r/MAE/RMSE vs null)
  ->  permutation test (significance)  ->  split-conformal interval coverage.

Outputs: data/processed/baseline_startle_{sex}.json + reports/figures/*.png
Run: ./.venv/Scripts/python.exe scripts/04_baseline_startle.py
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")  # Windows console defaults to cp1252

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drophenopredict import genotypes, phenotypes, gblup  # noqa: E402

warnings.filterwarnings("ignore")  # silence urllib3 insecure-TLS notice for the data host
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

GRM_PATH = "data/processed/grm.npz"
OUT = Path("data/processed")
FIG = Path("reports/figures")
N_REPEATS = 10
N_PERM = 200


def run_one(K, iids, sex: str) -> dict:
    y = phenotypes.fetch_phenotype(f"startle_{sex}")
    Ksub, yv, lines = gblup.align(K, iids, y)
    n = len(lines)
    print(f"\n=== startle_{sex}: {n} lines (GRM ∩ phenotype) ===")

    h = gblup.reml(Ksub, yv)
    cv = gblup.cross_validate(Ksub, yv, n_repeats=N_REPEATS, seed=42)
    r2 = float(cv.r2_per_repeat.mean())
    perm = gblup.permutation_test(Ksub, yv, r2, n_perm=N_PERM)
    conf = gblup.conformal_coverage(Ksub, yv, alpha=0.10, n_repeats=N_REPEATS)

    res = {
        "trait": f"startle_{sex}",
        "n_lines": n,
        "y_mean": float(yv.mean()), "y_sd": float(yv.std(ddof=1)),
        "genomic_h2_REML": round(h.h2, 4),
        "cv": {
            "scheme": f"5-fold x {N_REPEATS} repeats, GBLUP, delta refit per fold",
            "R2_mean": round(r2, 4), "R2_sd": round(float(cv.r2_per_repeat.std()), 4),
            "pearson_r_mean": round(float(cv.r_per_repeat.mean()), 4),
            "MAE_mean": round(float(cv.mae_per_repeat.mean()), 4),
            "RMSE_mean": round(float(cv.rmse_per_repeat.mean()), 4),
            "null_MAE_mean": round(float(cv.null_mae_per_repeat.mean()), 4),
        },
        "permutation_test": {
            "n_perm": N_PERM, "p_value": round(perm["p_value"], 4),
            "null_R2_mean": round(perm["null_mean"], 4),
            "null_R2_q95": round(perm["null_q95"], 4),
        },
        "conformal_90pct_interval": {
            "target_coverage": conf["target_coverage"],
            "empirical_coverage": round(conf["empirical_coverage"], 3),
            "mean_width": round(conf["mean_interval_width"], 3),
        },
    }
    (OUT / f"baseline_startle_{sex}.json").write_text(json.dumps(res, indent=2))

    # figure: obs vs pred + permutation null
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].scatter(cv.y_true, cv.y_pred, s=18, alpha=0.7)
    lim = [min(cv.y_true.min(), cv.y_pred.min()), max(cv.y_true.max(), cv.y_pred.max())]
    ax[0].plot(lim, lim, "k--", lw=1)
    ax[0].set_xlabel("Observed startle (line mean)")
    ax[0].set_ylabel("GBLUP out-of-fold prediction")
    ax[0].set_title(f"startle_{sex}: CV R²={r2:.3f}, r={res['cv']['pearson_r_mean']:.3f}, "
                    f"h²={h.h2:.2f}")
    ax[1].hist(perm["null"], bins=30, color="gray", alpha=0.8)
    ax[1].axvline(r2, color="red", lw=2, label=f"observed R²={r2:.3f}")
    ax[1].set_xlabel("CV R² under permuted labels (null)")
    ax[1].set_ylabel("count")
    ax[1].set_title(f"Permutation test  p={perm['p_value']:.3f}")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(FIG / f"baseline_startle_{sex}.png", dpi=130)
    plt.close(fig)
    return res


def main() -> None:
    K, iids, meta = genotypes.load_grm(GRM_PATH)
    print("GRM:", {k: meta[k] for k in ("n_snps_used", "n_lines") if k in meta},
          "| n_lines", len(iids))
    summary = [run_one(K, iids, s) for s in ("female", "male")]
    print("\n================ BASELINE SUMMARY ================")
    for r in summary:
        c = r["cv"]
        print(f"{r['trait']:16s} n={r['n_lines']:3d}  h²={r['genomic_h2_REML']:.2f}  "
              f"CV_R²={c['R2_mean']:+.3f}±{c['R2_sd']:.3f}  r={c['pearson_r_mean']:.3f}  "
              f"MAE={c['MAE_mean']:.2f}(null {c['null_MAE_mean']:.2f})  "
              f"perm_p={r['permutation_test']['p_value']:.3f}  "
              f"cover90={r['conformal_90pct_interval']['empirical_coverage']:.2f}")
    print("\nArtifacts: data/processed/baseline_startle_*.json, reports/figures/*.png")


if __name__ == "__main__":
    main()
