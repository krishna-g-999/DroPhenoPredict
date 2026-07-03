"""Download the BCM-DMAS TMT proteomics file from Synapse (16 samples: PD/aSyn
vs Elav-Gal4 control, day 10 only). RNA-seq + metadata are already in hand
(user-downloaded); this script fetches ONLY the remaining proteomics file so it
can validate the top convergent transcripts (Findings 13/15) at the protein level.

>>> RUN THIS IN YOUR OWN TERMINAL — do not paste your token into chat <<<

1. You already accepted the syn34767207 Data Use Agreement (done, since the
   RNA-seq download worked). No new consent step needed.
2. Create a Personal Access Token if you don't have one:
   https://www.synapse.org -> profile icon -> Account Settings ->
   Personal Access Tokens -> Create New Token (scope: "Download").
3. Open a terminal in D:\\DroPhenoPredict and run (PowerShell):
       $env:SYNAPSE_AUTH_TOKEN = "paste-your-token-here"
       .\\.venv\\Scripts\\python.exe -m pip install synapseclient
       .\\.venv\\Scripts\\python.exe scripts\\38_download_proteomics.py
   (or in cmd.exe:  set SYNAPSE_AUTH_TOKEN=paste-your-token-here)
4. The token only lives in that terminal's environment for that session — it is
   never written to any file this script creates, and never appears in chat.

Output: data/raw/npjpd2025/proteomics/jshulman_flytmt16__Proteins.txt
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("data/raw/npjpd2025/proteomics")
PROTEOMICS_FILE_ID = "syn36911668"   # jshulman_flytmt16__Proteins.txt (verified, 2026-07)


def main() -> None:
    try:
        import synapseclient
    except ImportError:
        sys.exit("Install first:  ./.venv/Scripts/python.exe -m pip install synapseclient")

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("Set SYNAPSE_AUTH_TOKEN in THIS terminal first (see instructions at the "
                 "top of this file). Do not paste your token into chat.")

    OUT.mkdir(parents=True, exist_ok=True)
    syn = synapseclient.login(authToken=token)
    entity = syn.get(PROTEOMICS_FILE_ID, downloadLocation=str(OUT))
    print(f"\nDownloaded: {entity.path}")
    print("Next: scripts/39_proteomics_validation.py will check whether the top "
          "convergent transcripts (Findings 13/15) show concordant protein-level changes "
          "in the PD/aSyn arm.")


if __name__ == "__main__":
    main()
