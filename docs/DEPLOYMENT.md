# Deployment — Streamlit Community Cloud

DroPhenoPredict follows the same deployment pattern as this account's other Streamlit app
(SAI-BBB-Predictor): a public/private GitHub repo + Streamlit Community Cloud (free, connects
directly to GitHub, auto-redeploys on push).

## Why this is the right choice here
- **Free**, zero server maintenance, matches an already-proven pattern on this account.
- **Data footprint is small enough to commit directly** — verified: the exact files `app.py` needs
  (GRM, model artifacts, benchmark/GWAS/GSEA result tables) total **~2.4 MB**. No Git LFS, no external
  object storage needed.
- **DUA-gated raw data is excluded** (`.gitignore`) — the large/restricted Synapse downloads
  (BCM-DMAS RNA-seq, proteomics) are never committed, both because redistribution would likely violate
  their data-use terms and because one of them (103 MB) exceeds GitHub's 100 MB hard file-size limit
  anyway. The deployed app only needs our own *derived, small* result files, not the raw inputs.

## What's set up (this repo)
- `requirements.txt` — **lean, deploy-only** dependencies (verified: exactly what `app.py`'s import
  chain needs — streamlit, pandas, numpy, scipy, plotly, bed-reader, requests). The full research
  stack (xgboost, lightgbm, shap, gseapy, ...) lives in `requirements-dev.txt` and is NOT installed on
  the deployed app, to keep cloud builds fast and reduce failure surface.
- `runtime.txt` / `.python-version` — pinned to Python 3.12 (matches the already-working
  SAI-BBB-Predictor deployment; safer than assuming this dev environment's 3.13 is supported).
- `.streamlit/config.toml` — basic theme/server config.
- `.gitignore` — excludes all raw/DUA-gated/large data by default, with explicit exceptions for the
  ~18 small files the app actually reads (see comments in the file itself).

## To go live (steps requiring your action)
1. **Create an empty GitHub repository** (github.com/new) — decide the name (e.g.
   `DroPhenoPredict`) and visibility (public/private — Streamlit Community Cloud's free tier deploys
   both, but public is simpler for reviewers/collaborators to browse the code).
2. Tell me the repo URL and I'll add the remote and push the existing local commit (I won't push
   without your explicit go-ahead — that's a publishing action).
3. Go to **share.streamlit.io**, sign in with GitHub, click "New app", select the repo/branch
   (`main`)/file (`app.py`).
4. Streamlit Cloud installs `requirements.txt` (the lean one) and boots the app — should take under a
   minute given the small footprint.

## Verified locally before any of this
- `app.py` boots clean (HTTP 200, no console errors) with the exact lean `requirements.txt` import set.
- All ~18 files it reads exist and are correctly *not* gitignored; all large/DUA-gated files are
  correctly gitignored (verified with `git check-ignore` on every relevant path, not assumed).
