"""Download supplementary files for DGRP drug/toxin studies not covered by DGRPool values.

Candidates come from a Europe PMC search (2026-11-07) for DGRP papers on drug,
toxin, insecticide, metal, anaesthetic or oxidant exposure
(data/interim/chemical/europepmc_search.csv), checked by title. Whether each
paper actually used DGRP lines is decided later from the downloaded tables
themselves (script 47), not assumed here.

Each paper's files go into their own folder:
    data/raw/chemical_supp/<label>/supplementary.zip   (as served, untouched)
    data/raw/chemical_supp/<label>/files/...           (extracted copy)
A manifest records URL, HTTP status, size and SHA-256 for every download.

Downloaded files are treated as data only: nothing in them is executed.
"""
from __future__ import annotations

import hashlib
import io
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd
import requests

# label -> (PMCID, DOI, exposure)
CANDIDATES = {
    "kislukhin2013_methotrexate": ("PMC3737169", "10.1534/g3.113.006619", "methotrexate (chemotherapy)"),
    "denecke2017_imidacloprid": ("PMC5596007", "10.1038/s41598-017-11092-5", "imidacloprid (insecticide)"),
    "mitchell2017_amanitin": ("PMC5328632", "10.1371/journal.pone.0173162", "alpha-amanitin (toxin)"),
    "schmidt2017_ddt": ("PMC5676240", "10.1534/genetics.117.300310", "DDT (insecticide)"),
    "highfill2017_nicotine": ("PMC5586381", "10.1534/genetics.117.300058", "nicotine"),
    "zhou2017_heavymetals": ("PMC5544243", "10.1371/journal.pgen.1006907", "cadmium, lead"),
    "duneau2018_insecticides": ("PMC6222576", "10.1534/g3.118.200537", "insecticides"),
    "battlay2018_insecticide": ("PMC6222580", "10.1534/g3.118.200619", "insecticides"),
    "olufs2018_anesthetics": ("PMC5799260", "10.1038/s41598-018-20720-7", "isoflurane, sevoflurane"),
    "highfill2019_cocaine_meth": ("PMC6527214", "10.1371/journal.pgen.1007834", "cocaine, methamphetamine"),
    "scharenbrock2020_anesthetic_tbi": ("PMC7555144", "10.3390/ijms21186926", "anaesthetics + brain injury"),
    "adebambo2020_btex": ("PMC7593870", "10.3389/fgene.2020.594179", "BTEX solvents"),
    "ramnarine2022_oxidative_tx": ("PMC8727983", "10.1093/g3journal/jkab366", "oxidative stress transcriptome"),
    "karam2022_amphetamine": ("PMC8894854", "10.3389/fpsyt.2022.831597", "amphetamine"),
    "green2022_copper": ("PMC9733279", "10.1186/s12915-022-01479-w", "copper"),
    "harrison2024_rapamycin": ("PMC11561674", "10.1111/acel.14292", "rapamycin"),
    "collins2024_4mei": ("PMC11533318", "10.1186/s40246-024-00689-3", "4-methylimidazole"),
    "tuifua2025_proteasome": ("PMC12342305", "10.1371/journal.pgen.1011823", "proteasome inhibitor"),
    "johnsonschlitz2025_anesthetic_tbi": ("PMC12677224", "10.1089/neu.2024.0314", "anaesthetic preconditioning"),
    "marconcini2026_toxin": ("PMC13042304", "10.1093/g3journal/jkag032", "natural toxin"),
    "ward2026_venom": ("PMC12978544", "10.1093/molbev/msag048", "venom"),
    "peddibhotla2026_rapamycin": ("PMC12896840", "10.3390/cells15030236", "rapamycin"),
}

ROOT = Path("data/raw/chemical_supp")
UA = {"User-Agent": "DroPhenoPredict/0.1 (academic research)"}


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    log = []
    for label, (pmcid, doi, exposure) in CANDIDATES.items():
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        rec = {"label": label, "pmcid": pmcid, "doi": doi, "exposure": exposure, "url": url,
               "downloaded": str(date.today())}
        try:
            r = requests.get(url, headers=UA, timeout=180)
            rec["http_status"] = r.status_code
            if r.status_code == 200 and r.content[:2] == b"PK":
                d = ROOT / label
                (d / "files").mkdir(parents=True, exist_ok=True)
                (d / "supplementary.zip").write_bytes(r.content)
                rec["bytes"] = len(r.content)
                rec["sha256"] = hashlib.sha256(r.content).hexdigest()
                with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                    names = []
                    for info in z.infolist():
                        p = Path(info.filename)
                        # extract regular files only, flattened, no path traversal
                        if info.is_dir() or p.is_absolute() or ".." in p.parts:
                            continue
                        (d / "files" / p.name).write_bytes(z.read(info))
                        names.append(p.name)
                rec["files"] = "; ".join(names)
            else:
                rec["note"] = "no supplementary zip served"
        except Exception as e:  # noqa: BLE001 - record and continue
            rec["error"] = str(e)
        print(f"{label:38s} {rec.get('http_status')} {rec.get('bytes', '')} {rec.get('files', rec.get('note', rec.get('error', '')))[:150]}")
        log.append(rec)
    pd.DataFrame(log).to_csv(ROOT / "download_manifest.csv", index=False)


if __name__ == "__main__":
    main()
