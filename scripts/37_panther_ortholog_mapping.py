"""Map the Ruffini et al. 2020 139-gene human shared-NDD list to Drosophila
orthologs via the PANTHER classification system ortholog API (pantherdb.org),
then test whether they overlap our BCM-DMAS fly-convergent genes more than
expected by chance.

PANTHER is a peer-reviewed, curated ortholog resource (used by the GO
Consortium). Verified reachable and correct on spot checks: MAPT->tau,
APP->Appl, SNCA->no fly ortholog (matches our existing curated table in
disease.py, which encodes SNCA as a human-transgene-only model).

Output: data/processed/ruffini139_fly_orthologs.csv
        data/processed/ortholog_overlap_test.csv
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

PROC = Path("data/processed")
H = {"User-Agent": "DroPhenoPredict/0.1 (research)"}
PANTHER_URL = "http://www.pantherdb.org/services/oai/pantherdb/ortholog/matchortho"
CHUNK = 10   # verified limit: 10 works, 15 fails with "more items than supported"


def query_panther(genes: list[str]) -> list[dict]:
    gene_str = ",".join(genes)
    params = {"geneInputList": gene_str, "organism": "9606",
             "targetOrganism": "7227", "targetGeneList": "", "orthologType": "all"}
    for attempt in range(5):
        try:
            r = requests.get(PANTHER_URL, params=params, headers=H, timeout=60, verify=False)
            r.raise_for_status()
            j = r.json()
            err = j.get("search", {}).get("error")
            if err:
                raise RuntimeError(f"PANTHER API error: {err}")
            mapped = j.get("search", {}).get("mapping", {}).get("mapped", [])
            if isinstance(mapped, dict):    # single-result responses aren't wrapped in a list
                mapped = [mapped]
            return mapped
        except Exception as e:
            print(f"  retry ({e}); attempt {attempt+1}/5")
            time.sleep(8)
    raise RuntimeError(f"PANTHER query failed after 5 attempts for chunk: {genes}")


def main() -> None:
    genes_df = pd.read_csv(PROC / "ruffini2020_139_shared_ndd_genes.csv")
    symbols = genes_df["geneSymbol"].dropna().unique().tolist()
    print(f"Querying PANTHER for {len(symbols)} human genes -> fly orthologs "
          f"(chunks of {CHUNK})...")

    checkpoint = PROC / "ruffini139_fly_orthologs_partial.csv"
    rows = []
    if checkpoint.exists():
        prior = pd.read_csv(checkpoint)
        rows = prior.to_dict("records")
        done_genes = set(prior["human_gene"])
        print(f"  resuming from checkpoint: {len(done_genes)} genes already queried")
    else:
        done_genes = set()

    for i in range(0, len(symbols), CHUNK):
        chunk = [g for g in symbols[i:i + CHUNK] if g not in done_genes]
        if not chunk:
            continue
        mapped = query_panther(chunk)
        queried_ids = {m.get("id") for m in mapped}
        for m in mapped:
            rows.append({
                "human_gene": m.get("id"),
                "fly_gene_symbol": m.get("target_gene_symbol"),
                "fly_gene_full": m.get("target_gene"),
                "ortholog_type": m.get("ortholog"),
            })
        # PANTHER omits genes it fully fails to recognize; record them as no-match too
        for g in chunk:
            if g not in queried_ids:
                rows.append({"human_gene": g, "fly_gene_symbol": None,
                            "fly_gene_full": None, "ortholog_type": None})
        pd.DataFrame(rows).to_csv(checkpoint, index=False)   # checkpoint after every chunk
        print(f"  [{i+len(chunk)}/{len(symbols)}] queried (checkpointed)")
        time.sleep(3)

    out = pd.DataFrame(rows)
    out["FBgn"] = out["fly_gene_full"].str.extract(r"FlyBase=(FBgn\d+)")
    out.to_csv(PROC / "ruffini139_fly_orthologs.csv", index=False)

    genes_with_ortholog = out.loc[out["FBgn"].notna(), "human_gene"].nunique()
    print(f"\n{len(symbols)} human genes queried; {len(out)} human-fly relationship rows "
          f"(some genes have multiple fly paralogs); {out['FBgn'].nunique()} unique fly orthologs; "
          f"{genes_with_ortholog}/{len(symbols)} human genes ({genes_with_ortholog/len(symbols)*100:.0f}%) "
          f"resolved to >=1 fly ortholog.")

    # ---- overlap test vs our BCM-DMAS convergent genes ----
    conv = pd.read_csv(PROC / "shared_biology_convergent_genes.csv")
    our_convergent_fbgn = set(conv["Geneid"])
    ruffini_fbgn = set(out["FBgn"].dropna())

    overlap = our_convergent_fbgn & ruffini_fbgn
    print(f"\nOur BCM-DMAS convergent genes (DE in >=3/5 models): {len(our_convergent_fbgn)}")
    print(f"Ruffini 139-gene orthologs resolved to fly: {len(ruffini_fbgn)}")
    print(f"Overlap: {len(overlap)} gene(s): {sorted(overlap)}")

    # background = all genes expressed in BCM-DMAS (from the counts matrix)
    counts = pd.read_parquet(PROC / "bcm_dmas_counts.parquet")
    cpm = counts / counts.sum(0) * 1e6
    background = set(counts.index[(cpm >= 1).sum(1) >= (0.5 * counts.shape[1])])
    ruffini_in_bg = ruffini_fbgn & background
    conv_in_bg = our_convergent_fbgn & background
    overlap_in_bg = conv_in_bg & ruffini_in_bg

    # hypergeometric: is overlap more than expected given background/set sizes?
    M = len(background)
    n = len(ruffini_in_bg)      # "successes" in population (Ruffini orthologs present in our data)
    N = len(conv_in_bg)         # draws (our convergent gene set, restricted to background)
    k = len(overlap_in_bg)
    p = stats.hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
    expected = N * n / M if M else 0

    result = pd.DataFrame([{
        "background_genes": M, "ruffini_orthologs_in_background": n,
        "our_convergent_genes_in_background": N, "observed_overlap": k,
        "expected_overlap": round(expected, 3), "hypergeom_p": p,
    }])
    result.to_csv(PROC / "ortholog_overlap_test.csv", index=False)
    print(f"\n=== Cross-species overlap test ===")
    print(f"Background (expressed genes, BCM-DMAS): {M}")
    print(f"Ruffini 139-gene orthologs present in background: {n}")
    print(f"Our convergent genes present in background: {N}")
    print(f"Observed overlap: {k} (expected by chance: {expected:.3f})")
    print(f"Hypergeometric p-value: {p:.4g}")


if __name__ == "__main__":
    main()
