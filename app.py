"""DroPhenoPredict — web interface (NAR Web Server style).

Wraps the validated DroPhenoPredictor (v1.1 multi-modal) and the benchmark /
genetic-architecture results into an honest, calibrated prediction tool.

Run:  ./.venv/Scripts/streamlit run app.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from drophenopredict.predictor import DroPhenoPredictor  # noqa: E402

PROC = Path("data/processed")
st.set_page_config(page_title="DroPhenoPredict", layout="wide", page_icon="🪰")


@st.cache_resource
def get_predictor():
    return DroPhenoPredictor("v1.1")


@st.cache_data
def load_csv(name):
    p = PROC / name
    return pd.read_csv(p) if p.exists() else None


pred = get_predictor()
man = pred.manifest
panel_lines = sorted({l for m in man["traits"].values() for l in m["lines"]},
                     key=lambda s: int(s.split("_")[1]))

st.sidebar.title("🪰 DroPhenoPredict")
st.sidebar.caption(f"v{man['version']} · multi-modal DGRP phenotype prediction")
page = st.sidebar.radio("Page", ["Predict", "Benchmark & predictability",
                                 "Genetic architecture", "Disease convergence (v2)",
                                 "About / methods"])
st.sidebar.info("Calibrated **prior with uncertainty** for prioritizing DGRP "
                "experiments — not a precise point predictor. Read the intervals.")

# ----------------------------------------------------------------- Predict ----
if page == "Predict":
    st.title("Predict a DGRP line's phenotype profile")
    col1, col2 = st.columns([1, 2])
    with col1:
        line = st.selectbox("DGRP line (in panel)", panel_lines,
                            index=panel_lines.index("line_357") if "line_357" in panel_lines else 0)
        st.caption("Accepts line_NN; novel/unrelated genotypes are out of scope.")
    prof = pred.predict_profile(line)
    rows = []
    for r in prof["predictions"]:
        if "error" in r:
            rows.append({"trait": r["trait"], "note": r["error"]})
            continue
        rows.append({
            "trait": r["trait"], "modality": r["modality"], "predicted": r["predicted"],
            "PI90 low": r["pi90_lower"], "PI90 high": r["pi90_upper"],
            "observed": r.get("observed", None), "unit": r["unit"],
            "predictable": "✓" if r["predictable"] else "—", "cv_r": r["model_cv_r"],
        })
    df = pd.DataFrame(rows)
    with col2:
        st.dataframe(df, use_container_width=True, hide_index=True)
    # plot predicted vs observed with intervals
    plot = df[df["predictable"] == "✓"].dropna(subset=["predicted"])
    if len(plot):
        fig = go.Figure()
        fig.add_bar(x=plot["trait"], y=plot["predicted"],
                    error_y=dict(type="data",
                                 array=plot["PI90 high"] - plot["predicted"],
                                 arrayminus=plot["predicted"] - plot["PI90 low"]),
                    name="predicted ±90% PI")
        if plot["observed"].notna().any():
            fig.add_scatter(x=plot["trait"], y=plot["observed"], mode="markers",
                            marker=dict(size=11, color="black", symbol="diamond"),
                            name="observed")
        fig.update_layout(title=f"{line}: predicted (bars, 90% PI) vs observed (♦)",
                          height=420, yaxis_title="trait value (native units)")
        st.plotly_chart(fig, use_container_width=True)
    st.caption("Modality is chosen per trait by nested CV (geno = genotype GBLUP, "
               "exp = RNA-seq expression). '—' = not better than the trait mean.")

# -------------------------------------------------- Benchmark & predictability
elif page == "Benchmark & predictability":
    st.title("Benchmark & predictability")
    bench = load_csv("benchmark_accuracy.csv")
    if bench is not None:
        st.subheader("Cross-validated accuracy (Pearson r) by trait")
        fig = go.Figure()
        for col, name in [("geno_GBLUP", "Genotype GBLUP"),
                          ("expression", "Expression (RNA-seq)"),
                          ("multimodal_adaptive", "Multi-modal adaptive")]:
            fig.add_bar(x=bench["trait_sex"], y=bench[col], name=name)
        fig.update_layout(barmode="group", height=460, yaxis_title="CV predictive r")
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(bench, use_container_width=True, hide_index=True)
        st.markdown("**Reproduces published benchmarks** — Ober 2012 (genotype "
                    "starvation r≈0.24) and Morgante 2020 (expression starvation♂ r≈0.38). "
                    "Multi-modal wins where genotype fails (climbing, starvation♂).")

# ----------------------------------------------------- Genetic architecture ---
elif page == "Genetic architecture":
    st.title("Genetic architecture of the predictable traits")
    keys = [k for k, m in man["traits"].items() if m["predictable"]]
    key = st.selectbox("Trait", keys)
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("GWAS candidate genes (inversion/Wolbachia-adjusted)")
        g = load_csv(f"gwas_adj_genes_{key}.csv")
        st.dataframe(g.head(20) if g is not None else pd.DataFrame(),
                     use_container_width=True, hide_index=True)
    with c2:
        st.subheader("Expression-associated genes (QTT, FDR)")
        q = load_csv(f"qtt_top_{key}.csv")
        st.dataframe(q.head(20) if q is not None else pd.DataFrame(),
                     use_container_width=True, hide_index=True)
    st.subheader("Pathway activity → trait associations (GO biological process)")
    pa = load_csv(f"pathway_assoc_{key}.csv")
    if pa is not None:
        st.dataframe(pa.head(15), use_container_width=True, hide_index=True)
        st.caption("Metabolic pathways flagged in `is_metabolic`. e.g. starvation♂ → lipid "
                   "β-oxidation, starvation♀ → glycogen catabolism (sex-specific energy strategy).")
    st.caption("GWAS = DNA-variant association (single-marker, covariate-adjusted). "
               "QTT = gene-expression association. Pathway = GO-set activity association. "
               "All correlational, underpowered at n≈190 (none FDR<0.05) — candidate/convergent "
               "signals, not definitive. Convergence across methods (e.g. lipid genes + lipid "
               "pathway for starvation) is the strength.")

# ---------------------------------------------- Disease convergence (v2) ------
elif page == "Disease convergence (v2)":
    st.title("Cross-disease molecular convergence (BCM-DMAS, real matched data)")
    st.markdown("""
    **Separate resource from the DGRP predictor above.** Uses the BCM-DMAS multi-disease fly panel
    (Synapse syn34767207 — Farrell/Shulman, *npj Parkinson's Disease* 2025): head RNA-seq for
    **AD (Aβ42, Tau), PD (αSyn), HD (HTT-128Q, HTT-200Q)** vs a shared Elav-Gal4 driver control,
    timecourse days 2–57. Tests whether neurodegeneration models molecularly converge
    (Ruffini et al. 2020, *Cells* — the human multi-omics analog of this question).
    """)
    convp = load_csv("gsea_convergent_pathways.csv")
    loo = load_csv("convergence_leave_one_out.csv")
    conv_genes = load_csv("shared_biology_convergent_genes.csv")

    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Leave-one-out robustness (gene-level convergence)")
        if loo is not None:
            st.dataframe(loo, use_container_width=True, hide_index=True)
            st.caption("Genes DE in ≥3/4 remaining models after excluding each model in turn — "
                      "all perm p=0.002, so convergence is not driven by any single model.")
    with c2:
        st.subheader("Convergent genes (DE in ≥3 of 5 models)")
        if conv_genes is not None:
            st.dataframe(conv_genes.head(15), use_container_width=True, hide_index=True)

    st.subheader("Genome-wide GSEA: pathways significant (FDR<0.1) in ≥3 of 5 models")
    if convp is not None and len(convp):
        show = convp[["pathway", "n_models_sig", "concordant_direction",
                      "Abeta42", "Tau", "aSyn", "HTT128Q", "HTT200Q"]]
        st.dataframe(show.head(20), use_container_width=True, hide_index=True)
        st.markdown("""
        **Key finding:** nucleolar stress / rRNA biogenesis is robustly **down-regulated** in Aβ42,
        αSyn, and HTT-128Q (concordant direction) — a documented hallmark of proteotoxic stress and
        aging. **HTT-200Q is a consistent directional outlier** across three independent analyses in
        this project (DGRP transfer test, gene-level convergence, and this pathway-level GSEA) —
        likely reflecting a distinct, more acute disease process for the severe full-length model.
        **Tau shows 0 significant pathways genome-wide** — the weakest molecular signal of the five.
        """)
    st.caption("GSEA engine is locally implemented (standard `gseapy` package could not be installed "
              "in this environment) and independently verified against ground-truth calibration "
              "properties before use (see `docs/FINDINGS_15_gsea_convergence.md`). DE pipeline "
              "independently cross-validated against the paper's own published DESeq2 results "
              "(7,018 genes, Pearson r=0.77 — `docs/FINDINGS_14_de_cross_validation.md`).")
    st.warning("This is a **molecular convergence resource**, not a behaviour predictor. Matched "
              "behavioural values (climbing/turning/stumbling) for this panel exist only in the "
              "paper's Figure 1B, not as downloadable data — so this analysis does not, and cannot "
              "yet, predict climbing/lifespan for these disease models.")

# ------------------------------------------------------------- About ----------
else:
    st.title("About / methods")
    st.markdown(f"""
**DroPhenoPredict v{man['version']}** predicts DGRP organismal phenotypes from genotype and/or
RNA-seq expression, choosing the modality **per trait by nested cross-validation**, and reports
**calibrated 90% prediction intervals** with an honest **predictability flag**.

**Method.** {man['method']}

**What it is for.** Prioritizing *which DGRP lines to phenotype* before doing fly work — a calibrated
prior with uncertainty, not a precise point predictor (accuracy r ≈ 0.09–0.34).

**What it is NOT.** It does not predict novel/unrelated genotypes, transgenic disease models, or drug
responses — those need perturbation data. Deep learning and functional-stratification were tested and
do **not** improve accuracy at this scale; the gain comes from the RNA-seq modality.

**Data (all public, provenance tracked):** DGRP Freeze 2.0 genotypes; DGRPool phenotypes;
GSE117850 RNA-seq expression; FlyBase/dm3 annotation. See `docs/DATA_SOURCES.md`.

**Validation:** nested CV (unbiased), permutation significance, conformal coverage verified,
engine checked on known-truth simulations (41 automated tests). Reproduces Ober 2012 & Morgante 2020.

**Also in this tool:** a separate multi-disease molecular-convergence resource (BCM-DMAS panel;
see "Disease convergence (v2)") — real matched AD/PD/HD fly omics, not merged with the DGRP predictor.
    """)
    st.dataframe(pd.DataFrame([
        {"trait": k, "deployed": m["deployed_modality"], "n": m["n_lines"],
         "adaptive_r": m["adaptive_nested_r"], "predictable": m["predictable"],
         "90%_coverage": m.get("interval_coverage_90")}
        for k, m in man["traits"].items()]), use_container_width=True, hide_index=True)
