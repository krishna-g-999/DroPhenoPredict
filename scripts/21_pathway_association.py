"""Modern technique 3 — pathway-activity -> trait association map.

Formalizes the lipid->starvation finding across ALL GO biological processes:
per DGRP line, computes pathway-activity scores (mean z-scored RNA-seq expression
of each GO-BP gene set, ancestor-propagated), then associates each pathway with
each predictable trait (BH-FDR). Answers "which metabolic/biological pathways
drive which behaviours" using curated knowledge — the rigorous way to bring
metabolism/disease pathway information in (vs dumping heterogeneous raw studies).

Output: data/processed/pathway_assoc_<key>.csv  (+ console summary)
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
from drophenopredict import dgrpool, expression, pathways  # noqa: E402

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

TRAIT_PID = {"climbing": 1543, "startle": 2799, "starvation": 2798, "sleep": 1483}
SEXMAP = {"F": "female", "M": "male"}
PROC = Path("data/processed")
METAB_KEYWORDS = ("metabol", "lipid", "fatty", "carbohydrate", "glucose", "tca",
                  "oxidation", "mitochond", "atp", "glyco", "insulin", "energy")


def bh_fdr(p):
    p = np.asarray(p); n = p.size
    order = np.argsort(p); q = np.empty(n)
    q[order] = p[order] * n / (np.arange(n) + 1)
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    man = json.loads(Path("models/trained/drophenopredict_v1.1.json").read_text())
    predictable = [k for k, m in man["traits"].items() if m["predictable"]]
    gsets = pathways.load_go_genesets(aspect="P", min_genes=10, max_genes=300)
    print(f"GO-BP gene sets: {len(gsets)}")

    pa_cache = {}
    for key in predictable:
        trait, sex = key.rsplit("_", 1)
        y = dgrpool.line_means(dgrpool.fetch_phenotype_values(TRAIT_PID[trait]),
                               sex=sex, harmonize=True)
        sx = SEXMAP[sex]
        if sx not in pa_cache:
            pa_cache[sx] = pathways.pathway_activity(expression.load_rnaseq(sx), gsets)
        PA = pa_cache[sx]                                # lines x pathways
        common = [l for l in PA.index if l in set(y.index)]
        yv = y.loc[common].to_numpy(float)
        M = PA.loc[common].to_numpy(float)              # n x pathways
        yc = yv - yv.mean()
        Mc = M - M.mean(0)
        r = (yc @ Mc) / (np.sqrt((yc ** 2).sum()) * np.sqrt((Mc ** 2).sum(0)) + 1e-12)
        r = np.clip(r, -0.999, 0.999)
        n = len(common)
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
        q = bh_fdr(p)
        go_ids = PA.columns.to_numpy()
        names = [gsets[g][0] for g in go_ids]
        df = pd.DataFrame({"GO_id": go_ids, "pathway": names, "r": np.round(r, 3),
                           "p": p, "fdr": np.round(q, 4)})
        df["is_metabolic"] = df["pathway"].str.lower().str.contains("|".join(METAB_KEYWORDS))
        df = df.sort_values("p")
        df.head(60).to_csv(PROC / f"pathway_assoc_{key}.csv", index=False)
        nsig = int((q < 0.05).sum())
        topmet = df[df["is_metabolic"]].head(3)["pathway"].tolist()
        print(f"{key:13s} n={n:3d} paths={len(go_ids)} FDR<0.05={nsig:3d}  "
              f"top: {df.iloc[0]['pathway'][:40]} (r={df.iloc[0]['r']:+.2f}) | "
              f"metabolic: {topmet}")

    print(f"\nSaved pathway-trait association maps -> {PROC}/pathway_assoc_*.csv")


if __name__ == "__main__":
    main()
