"""Modern technique 2 — transcriptome-wide gene association (QTT analysis).

For each predictable trait, correlate every gene's RNA-seq expression (line means)
with the trait across lines -> "quantitative trait transcripts" (Morgante 2020).
This identifies the genes whose EXPRESSION state tracks the phenotype — the
mechanism behind the expression-modality prediction — complementing SNP GWAS.
Benjamini-Hochberg FDR control. Correlational (both line-level), not causal.

Output: data/processed/qtt_top_<key>.csv
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
from drophenopredict import dgrpool, expression, annotation  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAIT_PID = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
SEXMAP = {"F": "female", "M": "male"}
PROC = Path("data/processed")


def bh_fdr(p):
    p = np.asarray(p); n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = (p[order] * n) / (np.arange(n) + 1)
    # enforce monotonicity
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    predictable = [k for k, m in man["traits"].items() if m["predictable"]]
    fbgn2sym = annotation.load_gene_spans().set_index("gene_id")["gene_name"].to_dict()

    for key in predictable:
        trait, sex = key.rsplit("_", 1)
        y = dgrpool.line_means(dgrpool.fetch_phenotype_values(TRAIT_PID[trait]),
                               sex=sex, harmonize=True)
        expr = expression.load_rnaseq(SEXMAP[sex])           # genes x lines
        common = [l for l in expr.columns if l in set(y.index)]
        yv = y.loc[common].to_numpy(float)
        X = expr[common].to_numpy(float)                     # genes x n
        keep = np.isfinite(X).all(1) & (X.std(1) > 1e-9)
        X = X[keep]; gene_ids = expr.index[keep].to_numpy()
        n = len(common)
        yc = yv - yv.mean()
        Xc = X - X.mean(1, keepdims=True)
        r = (Xc @ yc) / (np.sqrt((Xc ** 2).sum(1)) * np.sqrt((yc ** 2).sum()))
        r = np.clip(r, -0.999, 0.999)
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
        q = bh_fdr(p)
        df = pd.DataFrame({"gene_id": gene_ids,
                           "gene_name": [fbgn2sym.get(g, g) for g in gene_ids],
                           "r": np.round(r, 3), "p": p, "fdr": np.round(q, 4)})
        df = df.sort_values("p")
        df.head(200).to_csv(PROC / f"qtt_top_{key}.csv", index=False)
        n_sig = int((q < 0.05).sum())
        top = ", ".join(df["gene_name"].head(6))
        print(f"{key:13s} n={n:3d} genes={len(gene_ids)} FDR<0.05={n_sig:4d}  top: {top}")

    print(f"\nSaved QTT tables -> {PROC}/qtt_top_*.csv")


if __name__ == "__main__":
    main()
