"""Tier 1 — build the DGRP metabolome matrix (3rd modality) from MTBLS2060.

Downloads the 673 Bruker NMR sample zips, extracts each AUTHORS'-processed real
spectrum (pdata/1/1r + procs — already FT'd/phased/baseline-corrected, so we use
their processing, not a re-derivation), bins to 0.04-ppm bins over 0.5-10 ppm
(water 4.5-5.1 excluded), constant-sum normalizes, and averages replicates per
line -> lines x bins matrix. Line ids harmonized to line_NN.

Output: data/processed/metabolome_matrix.parquet  (170 lines x ~237 bins)
        data/processed/metabolome_bins.csv         (bin ppm centers)
Run: ./.venv/Scripts/python.exe scripts/22_build_metabolome.py
"""
from __future__ import annotations

import io
import os
import re
import sys
import zipfile
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

BASE = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS2060/FILES/RAW_FILES"
RAW = Path("data/raw/dgrp_metabolome/raw")
PROC = Path("data/processed")
H = {"User-Agent": "DroPhenoPredict/0.1 (research)"}
EDGES = np.arange(0.5, 10.0001, 0.02)          # 0.02 ppm bins
WATER = (4.5, 5.1)


def fetch(src):
    dst = RAW / f"{src}.zip"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    for _ in range(3):
        try:
            r = requests.get(f"{BASE}/{src}.zip", headers=H, timeout=120, verify=False)
            r.raise_for_status()
            dst.write_bytes(r.content)
            return dst
        except Exception:
            continue
    return None


def parse_procs(txt):
    d = {}
    for line in txt.splitlines():
        m = re.match(r"##\$(\w+)= ?(.*)", line)
        if m:
            d[m.group(1)] = m.group(2).strip()
    return d


def spectrum_to_bins(zip_path):
    z = zipfile.ZipFile(zip_path)
    names = z.namelist()
    pr = next(n for n in names if n.endswith("pdata/1/procs"))
    rr = next(n for n in names if n.endswith("pdata/1/1r"))
    procs = parse_procs(z.read(pr).decode("latin-1"))
    SI = int(procs["SI"]); NC = float(procs["NC_proc"]); SF = float(procs["SF"])
    SW = float(procs["SW_p"]); OFF = float(procs["OFFSET"])
    dt = "<i4" if int(procs.get("BYTORDP", "0")) == 0 else ">i4"
    raw = np.frombuffer(z.read(rr), dtype=dt).astype(float) * (2.0 ** NC)
    if raw.size != SI:
        raw = raw[:SI]
    ppm = OFF - np.arange(SI) * (SW / SF) / SI
    mask = (ppm >= 0.5) & (ppm <= 10) & ~((ppm >= WATER[0]) & (ppm <= WATER[1]))
    which = np.digitize(ppm[mask], EDGES)
    vals = np.clip(raw[mask], 0, None)
    return np.array([vals[which == i].sum() for i in range(1, len(EDGES))])  # raw bins


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    smap = pd.read_csv("data/raw/dgrp_metabolome/sample_map.csv")
    sources = smap["source"].astype(int).tolist()

    print(f"Downloading {len(sources)} NMR zips (threaded)...")
    with ThreadPoolExecutor(max_workers=10) as ex:
        paths = list(ex.map(fetch, sources))
    ok = sum(p is not None for p in paths)
    print(f"  downloaded/cached: {ok}/{len(sources)}")

    rows, used = [], []
    for src, p in zip(sources, paths):
        if p is None:
            continue
        try:
            rows.append(spectrum_to_bins(p)); used.append(src)
        except Exception as e:
            print(f"  skip {src}: {e}")
    M = np.vstack(rows).astype(float)
    # sample-level normalization: constant-sum -> PQN -> log1p (recovers heritable signal)
    M = M / (M.sum(1, keepdims=True) + 1e-12)
    ref = np.median(M, axis=0) + 1e-12
    quot = np.median(M / ref, axis=1, keepdims=True)
    M = M / quot                                         # PQN
    M = np.log1p(M * 1e4)                                # compress dynamic range
    df = pd.DataFrame(M, index=used)
    lines = smap.set_index("source").loc[used, "dgrp"].values
    df["line"] = [f"line_{int(s.split('-')[1])}" for s in lines]   # DGRP-021 -> line_21
    per_line = df.groupby("line").mean()                # average replicates
    centers = (EDGES[:-1] + EDGES[1:]) / 2
    per_line.columns = [f"ppm_{c:.2f}" for c in centers]
    per_line.to_parquet(PROC / "metabolome_matrix.parquet")
    pd.Series(centers, name="ppm").to_csv(PROC / "metabolome_bins.csv", index=False)
    print(f"\nMetabolome matrix: {per_line.shape[0]} lines x {per_line.shape[1]} bins")
    print("  saved -> data/processed/metabolome_matrix.parquet")


if __name__ == "__main__":
    main()
