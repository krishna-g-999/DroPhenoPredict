"""Download the npj PD 2025 alpha-synuclein Drosophila OMICS from Synapse.

Data: Baylor "BCM-DMAS" project, syn34767207 (Farrell/Shulman et al., npj
Parkinson's Disease 2025). RNA-seq counts + proteomics + metadata.
NOTE: the behavioural data (climbing/turning/stumbling) is NOT in Synapse — it
is in the paper's open-access supplementary (fetched separately).

>>> ONE-TIME ACCESS SETUP (you must do this — it's your consent, not mine) <<<
1. Log in at https://www.synapse.org and open  https://www.synapse.org/Synapse:syn34767207
2. If it shows an access requirement / Data Use agreement, click through & ACCEPT it.
   (Without this, downloads return HTTP 403 "unmet access requirements".)
3. Create a Personal Access Token: Synapse -> Account Settings ->
   Personal Access Tokens -> generate one with 'Download' scope.
4. In THIS terminal (do not paste the token into chat):
       set SYNAPSE_AUTH_TOKEN=<your token>        (Windows CMD)
       $env:SYNAPSE_AUTH_TOKEN="<your token>"     (PowerShell)
5. pip install synapseclient   (into the project venv)
6. ./.venv/Scripts/python.exe scripts/30_download_npjpd_synapse.py

Output: data/raw/npjpd2025/{metadata,counts,proteomics}/
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("data/raw/npjpd2025")
# exact Synapse entities (verified from the public folder tree, 2026-07)
TARGETS = {
    "metadata": "syn34767252",     # biospecimen + assay metadata CSVs (the sample->genotype/age key)
    "counts": "syn36328638",       # RNA-seq per-sample counts (genotypes A/B/C/D x days 2-57 x 3 reps)
    "proteomics": "syn34767275",   # TMT proteomics (optional)
}


def main() -> None:
    try:
        import synapseclient
        import synapseutils
    except ImportError:
        sys.exit("Install first:  ./.venv/Scripts/python.exe -m pip install synapseclient")

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("Set SYNAPSE_AUTH_TOKEN (see the setup instructions at the top of this file).")

    syn = synapseclient.login(authToken=token)
    for name, syn_id in TARGETS.items():
        dest = OUT / name
        dest.mkdir(parents=True, exist_ok=True)
        print(f"\n=== downloading {name} ({syn_id}) -> {dest}")
        try:
            synapseutils.syncFromSynapse(syn, syn_id, path=str(dest))
        except Exception as e:
            print(f"  FAILED for {name}: {e}")
            print("  If this is HTTP 403, accept the Data Use agreement on the "
                  "syn34767207 page first (step 2 above).")
    print("\nDone. Next: scripts/31_build_npjpd_matrix.py will assemble the counts + "
          "metadata into a genes x (genotype,age) matrix for the transfer experiment.")


if __name__ == "__main__":
    main()
