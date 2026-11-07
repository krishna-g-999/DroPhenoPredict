"""Download DGRP drug/toxin datasets deposited in figshare and Dryad.

Found by searching each candidate paper's DOI in the Dryad, Zenodo and figshare
APIs (2026-11-07). Files land in their own folder per dataset:
    data/raw/chemical_repo/<label>/<file>
with a manifest of URL, size and SHA-256. Files are data only; nothing is executed.
"""
from __future__ import annotations

import hashlib
import io
import urllib.parse
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd
import requests

FIGSHARE = {  # label -> figshare article id
    "duneau2018_insecticides": 7047734,
    "battlay2018_insecticide": 6954443,
    "rohde2019_drug_response": 9368636,
}
DRYAD = {  # label -> dataset DOI
    "kislukhin2013_methotrexate": "doi:10.5061/dryad.d20qc",
}

ROOT = Path("data/raw/chemical_repo")
H = {"User-Agent": "DroPhenoPredict/0.1 (academic research)"}


def save(label: str, name: str, content: bytes, url: str, log: list) -> None:
    d = ROOT / label
    d.mkdir(parents=True, exist_ok=True)
    safe = Path(name).name
    (d / safe).write_bytes(content)
    log.append({"label": label, "file": safe, "url": url, "bytes": len(content),
                "sha256": hashlib.sha256(content).hexdigest(), "downloaded": str(date.today())})
    print(f"  {label}: {safe} ({len(content):,} B)")
    if Path(safe).suffix.lower() == ".zip":  # unpack archives next to the original, regular files only
        out = d / ("unzipped_" + Path(safe).stem)
        out.mkdir(exist_ok=True)
        with zipfile.ZipFile(io.BytesIO(content)) as z:
            for info in z.infolist():
                p = Path(info.filename)
                if info.is_dir() or p.is_absolute() or ".." in p.parts or p.name.startswith("._"):
                    continue
                (out / p.name).write_bytes(z.read(info))


def main() -> None:
    log: list = []
    for label, aid in FIGSHARE.items():
        meta = requests.get(f"https://api.figshare.com/v2/articles/{aid}", headers=H, timeout=60).json()
        for f in meta.get("files", []):
            r = requests.get(f["download_url"], headers=H, timeout=300)
            if r.ok:
                save(label, f["name"], r.content, f["download_url"], log)
            else:
                print(f"  {label}: {f['name']} HTTP {r.status_code}")
    for label, doi in DRYAD.items():
        enc = urllib.parse.quote(doi, safe="")
        vers = requests.get(f"https://datadryad.org/api/v2/datasets/{enc}/versions", headers=H, timeout=60).json()
        latest = vers["_embedded"]["stash:versions"][-1]["_links"]["self"]["href"]
        files = requests.get(f"https://datadryad.org{latest}/files", headers=H, timeout=60).json()
        for f in files["_embedded"]["stash:files"]:
            url = "https://datadryad.org" + f["_links"]["stash:download"]["href"]
            r = requests.get(url, headers=H, timeout=300)
            if r.ok:
                save(label, f["path"], r.content, url, log)
            else:
                print(f"  {label}: {f['path']} HTTP {r.status_code}")
    pd.DataFrame(log).to_csv(ROOT / "download_manifest.csv", index=False)


if __name__ == "__main__":
    main()
