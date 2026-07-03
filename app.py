"""DroPhenoPredict — web interface (NAR Web Server style).

Wraps the validated DroPhenoPredictor (v1.1 multi-modal) and the benchmark /
genetic-architecture results into an honest, calibrated prediction tool.

Visual design shares the navy/orange system used by this lab's other public tool,
DroPhenix (github.com/krishna-g-999/DroPhenix), for brand consistency across both.

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
st.set_page_config(page_title="DroPhenoPredict", layout="wide", page_icon="🪰",
                    initial_sidebar_state="expanded")

# --------------------------------------------------------------------- design ----

def inject_css() -> None:
    st.markdown("""
    <style>
    .dpp-hero {
        background: linear-gradient(135deg, #0D2137 0%, #09192A 100%);
        border-radius: 18px; padding: 2.1rem 2.4rem; margin-bottom: 1.4rem;
        box-shadow: 0 10px 40px rgba(0,0,0,0.30);
    }
    .dpp-hero h1 { color: #FFFFFF; font-weight: 900; margin: 0 0 .35rem 0; font-size: 2.1rem; }
    .dpp-hero p { color: #C5D7F0; font-size: 1.02rem; margin: 0; line-height: 1.6; max-width: 62ch; }
    .dpp-hero .dpp-accent { color: #F0A500; }

    .dpp-card {
        background: #FFFFFF; border: 1px solid #dde3ec; border-radius: 12px;
        padding: 1.2rem 1.35rem; margin-bottom: 1rem; height: 100%;
    }
    .dpp-card h3 { color: #0D2137; font-weight: 800; font-size: 1.05rem; margin: 0 0 .4rem 0; }
    .dpp-card p { color: #475569; font-size: 0.88rem; line-height: 1.55; margin: 0 0 .5rem 0; }

    .dpp-section-hdr {
        color: #0D2137; font-weight: 800; font-size: 1.3rem;
        border-left: 5px solid #F0A500; padding-left: .6rem; margin: .2rem 0 1rem 0;
    }

    .dpp-badge {
        display: inline-block; border-radius: 999px; padding: 3px 13px;
        font-size: 0.76rem; font-weight: 700; margin: 2px 4px 2px 0; letter-spacing: .01em;
    }
    .dpp-badge-ok    { background: rgba(34,197,94,.13);  color: #15803d; border: 1px solid rgba(34,197,94,.35); }
    .dpp-badge-warn  { background: rgba(240,165,0,.13);  color: #92600a; border: 1px solid rgba(240,165,0,.35); }
    .dpp-badge-off   { background: rgba(148,163,184,.18); color: #475569; border: 1px solid rgba(148,163,184,.4); }
    .dpp-badge-info  { background: rgba(59,130,246,.13); color: #1d4ed8; border: 1px solid rgba(59,130,246,.35); }

    .dpp-footer { color: #64748B; font-size: 0.78rem; margin-top: 2rem; padding-top: 1rem;
                  border-top: 1px solid #dde3ec; }
    section[data-testid="stSidebar"] { background: #0D2137; }
    section[data-testid="stSidebar"] * { color: #EEF2F9 !important; }
    section[data-testid="stSidebar"] .dpp-side-title { color: #F0A500 !important; font-weight: 900; }
    </style>
    """, unsafe_allow_html=True)


def hero(title: str, subtitle: str) -> None:
    st.markdown(f"""
    <div class="dpp-hero"><h1>{title}</h1><p>{subtitle}</p></div>
    """, unsafe_allow_html=True)


def section_hdr(text: str) -> None:
    st.markdown(f'<div class="dpp-section-hdr">{text}</div>', unsafe_allow_html=True)


def badge(text: str, kind: str = "info") -> str:
    return f'<span class="dpp-badge dpp-badge-{kind}">{text}</span>'


def feature_card(title: str, body: str, badges_html: str = "") -> str:
    return f'<div class="dpp-card"><h3>{title}</h3><p>{body}</p>{badges_html}</div>'


inject_css()


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
n_predictable = sum(1 for m in man["traits"].values() if m["predictable"])

# --------------------------------------------------------------------- sidebar ---
st.sidebar.markdown('<p class="dpp-side-title" style="font-size:1.3rem;margin-bottom:0;">🪰 DroPhenoPredict</p>',
                     unsafe_allow_html=True)
st.sidebar.caption(f"{man['version']} · multi-modal DGRP phenotype prediction")
page = st.sidebar.radio("Page", ["Home", "Predict", "Benchmark & predictability",
                                 "Genetic architecture", "Disease convergence (v2)",
                                 "About / methods"], label_visibility="collapsed")
st.sidebar.markdown("---")
st.sidebar.markdown(f"**{len(panel_lines)}** DGRP lines in panel &nbsp;·&nbsp; "
                     f"**{n_predictable}/8** trait models predictable")
st.sidebar.info("Calibrated **prior with uncertainty** for prioritizing DGRP "
                "experiments — not a precise point predictor. Read the intervals.")
st.sidebar.markdown(
    '<div style="font-size:0.75rem;opacity:.75;margin-top:1rem;">Companion tool: '
    '<a href="https://drophenix.streamlit.app" style="color:#F0A500;">DroPhenix</a> '
    '— multi-assay Drosophila phenotyping & compound ranking.</div>',
    unsafe_allow_html=True)

# ------------------------------------------------------------------------ Home ---
if page == "Home":
    hero("DroPhenoPredict",
         "Calibrated, multi-modal phenotype prediction for the <span class='dpp-accent'>"
         "Drosophila Genetic Reference Panel</span> — genotype and RNA-seq expression, "
         "chosen per trait by nested cross-validation, with honest uncertainty and an "
         "explicit predictability flag on every prediction.")

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown(feature_card(
            "Predict", "Look up any of the 205 DGRP lines and get a full trait profile with "
            "90% prediction intervals and the modality (genotype or expression) used.",
            badge("8 trait models", "info") + badge(f"{n_predictable}/8 predictable", "ok")),
            unsafe_allow_html=True)
    with c2:
        st.markdown(feature_card(
            "Benchmark & architecture", "Cross-validated accuracy per trait, reproducing two "
            "published benchmarks, plus GWAS/expression/pathway signals behind each prediction.",
            badge("Ober 2012 ✓", "ok") + badge("Morgante 2020 ✓", "ok")),
            unsafe_allow_html=True)
    with c3:
        st.markdown(feature_card(
            "Disease convergence (v2)", "A separate, independently-verified resource: molecular "
            "convergence across 5 real AD/PD/HD Drosophila models on the same matched panel.",
            badge("Molecular only", "warn") + badge("Not a behaviour predictor", "off")),
            unsafe_allow_html=True)

    section_hdr("What this tool is — and isn't")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(feature_card(
            "Validated core", "Genotype (GBLUP) and RNA-seq expression prediction for climbing, "
            "startle, starvation and sleep, with calibrated 90% intervals (coverage-verified "
            "0.90–0.92), a permutation-based predictability flag, and reproduction of 3 literature "
            "benchmarks. Engine checked on known-truth simulations (41 automated tests)."),
            unsafe_allow_html=True)
    with c2:
        st.markdown(feature_card(
            "Explicitly out of scope", "Novel/unrelated genotypes, transgenic disease models, and "
            "drug responses need perturbation data this project doesn't have. A direct test of "
            "disease-phenotype transfer (Yang et al. 2023, DGRP×AD) was run and came back an "
            "honest, well-powered <b>null</b> — see About / methods."),
            unsafe_allow_html=True)

    st.markdown('<div class="dpp-footer">All data public and provenance-tracked '
                '(<code>docs/DATA_SOURCES.md</code>). No fabricated or placeholder values anywhere '
                'in this tool — every number is downloaded, computed, or explicitly flagged as a '
                'hypothesis. See the sibling tool <a href="https://drophenix.streamlit.app">'
                'DroPhenix</a> for assay-level Drosophila phenotyping and compound ranking.</div>',
                unsafe_allow_html=True)

# ----------------------------------------------------------------- Predict ----
elif page == "Predict":
    hero("Predict a DGRP line's phenotype profile",
         "Select an in-panel line to see its full calibrated trait profile — predicted value, "
         "90% interval, modality used, and whether the trait is predictable at all.")
    col1, col2 = st.columns([1, 2])
    with col1:
        st.markdown('<div class="dpp-card">', unsafe_allow_html=True)
        line = st.selectbox("DGRP line (in panel)", panel_lines,
                            index=panel_lines.index("line_357") if "line_357" in panel_lines else 0)
        st.caption("Accepts line_NN; novel/unrelated genotypes are out of scope.")
        st.markdown('</div>', unsafe_allow_html=True)
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
        fig.add_bar(x=plot["trait"], y=plot["predicted"], marker_color="#F0A500",
                    error_y=dict(type="data", color="#0D2137",
                                 array=plot["PI90 high"] - plot["predicted"],
                                 arrayminus=plot["predicted"] - plot["PI90 low"]),
                    name="predicted ±90% PI")
        if plot["observed"].notna().any():
            fig.add_scatter(x=plot["trait"], y=plot["observed"], mode="markers",
                            marker=dict(size=11, color="#0D2137", symbol="diamond"),
                            name="observed")
        fig.update_layout(title=f"{line}: predicted (bars, 90% PI) vs observed (♦)",
                          height=420, yaxis_title="trait value (native units)",
                          plot_bgcolor="#FFFFFF", paper_bgcolor="#FFFFFF")
        st.plotly_chart(fig, use_container_width=True)
    st.caption("Modality is chosen per trait by nested CV (geno = genotype GBLUP, "
               "exp = RNA-seq expression). '—' = not better than the trait mean.")

# -------------------------------------------------- Benchmark & predictability
elif page == "Benchmark & predictability":
    hero("Benchmark & predictability",
         "Cross-validated accuracy per trait, and where our numbers reproduce two prior "
         "published DGRP prediction studies.")
    bench = load_csv("benchmark_accuracy.csv")
    if bench is not None:
        section_hdr("Cross-validated accuracy (Pearson r) by trait")
        fig = go.Figure()
        for col, name, color in [("geno_GBLUP", "Genotype GBLUP", "#0D2137"),
                                  ("expression", "Expression (RNA-seq)", "#3B82F6"),
                                  ("multimodal_adaptive", "Multi-modal adaptive", "#F0A500")]:
            fig.add_bar(x=bench["trait_sex"], y=bench[col], name=name, marker_color=color)
        fig.update_layout(barmode="group", height=460, yaxis_title="CV predictive r",
                          plot_bgcolor="#FFFFFF", paper_bgcolor="#FFFFFF")
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(bench, use_container_width=True, hide_index=True)
        st.markdown(feature_card("Reproduces published benchmarks",
            "Ober 2012 (genotype starvation r≈0.24) and Morgante 2020 (expression starvation♂ "
            "r≈0.38). Multi-modal wins where genotype fails (climbing, starvation♂).",
            badge("Ober 2012 reproduced", "ok") + badge("Morgante 2020 reproduced", "ok")),
            unsafe_allow_html=True)

# ----------------------------------------------------- Genetic architecture ---
elif page == "Genetic architecture":
    hero("Genetic architecture of the predictable traits",
         "GWAS, expression-association (QTT), and pathway-activity signals behind each "
         "predictable trait — correlational and exploratory, not part of the predictor itself.")
    keys = [k for k, m in man["traits"].items() if m["predictable"]]
    key = st.selectbox("Trait", keys)
    c1, c2 = st.columns(2)
    with c1:
        section_hdr("GWAS candidate genes")
        st.caption("inversion/Wolbachia-adjusted, single-marker")
        g = load_csv(f"gwas_adj_genes_{key}.csv")
        st.dataframe(g.head(20) if g is not None else pd.DataFrame(),
                     use_container_width=True, hide_index=True)
    with c2:
        section_hdr("Expression-associated genes")
        st.caption("QTT, FDR-controlled")
        q = load_csv(f"qtt_top_{key}.csv")
        st.dataframe(q.head(20) if q is not None else pd.DataFrame(),
                     use_container_width=True, hide_index=True)
    section_hdr("Pathway activity → trait associations (GO biological process)")
    pa = load_csv(f"pathway_assoc_{key}.csv")
    if pa is not None:
        st.dataframe(pa.head(15), use_container_width=True, hide_index=True)
        st.caption("Metabolic pathways flagged in `is_metabolic`. e.g. starvation♂ → lipid "
                   "β-oxidation, starvation♀ → glycogen catabolism (sex-specific energy strategy).")
    st.markdown(badge("all correlational", "warn") + badge("underpowered at n≈190", "warn") +
                badge("none FDR<0.05 individually", "off"), unsafe_allow_html=True)
    st.caption("GWAS = DNA-variant association (single-marker, covariate-adjusted). "
               "QTT = gene-expression association. Pathway = GO-set activity association. "
               "Candidate/convergent signals, not definitive. Convergence across methods (e.g. "
               "lipid genes + lipid pathway for starvation) is the strength.")

# ---------------------------------------------- Disease convergence (v2) ------
elif page == "Disease convergence (v2)":
    hero("Cross-disease molecular convergence",
         "A separate resource from the DGRP predictor: real matched AD/PD/HD Drosophila "
         "head RNA-seq (BCM-DMAS panel), independently DE-tested and GSEA-verified.")
    st.markdown(feature_card("Data source", "BCM-DMAS multi-disease fly panel (Synapse "
        "syn34767207 — Farrell/Shulman, <i>npj Parkinson's Disease</i> 2025): head RNA-seq for "
        "<b>AD (Aβ42, Tau), PD (αSyn), HD (HTT-128Q, HTT-200Q)</b> vs a shared Elav-Gal4 driver "
        "control, timecourse days 2–57. Tests whether neurodegeneration models molecularly "
        "converge (Ruffini et al. 2020, <i>Cells</i> — the human multi-omics analog of this "
        "question)."), unsafe_allow_html=True)

    convp = load_csv("gsea_convergent_pathways.csv")
    loo = load_csv("convergence_leave_one_out.csv")
    conv_genes = load_csv("shared_biology_convergent_genes.csv")

    c1, c2 = st.columns(2)
    with c1:
        section_hdr("Leave-one-out robustness")
        st.caption("gene-level convergence")
        if loo is not None:
            st.dataframe(loo, use_container_width=True, hide_index=True)
            st.caption("Genes DE in ≥3/4 remaining models after excluding each model in turn — "
                      "all perm p=0.002, so convergence is not driven by any single model.")
    with c2:
        section_hdr("Convergent genes")
        st.caption("DE in ≥3 of 5 models")
        if conv_genes is not None:
            st.dataframe(conv_genes.head(15), use_container_width=True, hide_index=True)

    section_hdr("Genome-wide GSEA: pathways significant (FDR<0.1) in ≥3 of 5 models")
    if convp is not None and len(convp):
        show = convp[["pathway", "n_models_sig", "concordant_direction",
                      "Abeta42", "Tau", "aSyn", "HTT128Q", "HTT200Q"]]
        st.dataframe(show.head(20), use_container_width=True, hide_index=True)
        st.markdown(feature_card("Key finding",
            "Nucleolar stress / rRNA biogenesis is robustly <b>down-regulated</b> in Aβ42, αSyn, "
            "and HTT-128Q (concordant direction) — a documented hallmark of proteotoxic stress and "
            "aging. <b>HTT-200Q is a consistent directional outlier</b> across three independent "
            "analyses in this project. <b>Tau shows 0 significant pathways genome-wide</b> — the "
            "weakest molecular signal of the five.",
            badge("independently verified vs gseapy (r=1.000)", "ok")),
            unsafe_allow_html=True)
    st.caption("GSEA engine is locally implemented (standard `gseapy` package could not be "
              "installed originally) and independently verified — against ground-truth "
              "calibration properties and, separately, against real `gseapy` on identical input "
              "(ES Pearson r=1.000, `docs/FINDINGS_20`). DE pipeline cross-validated against the "
              "paper's own published DESeq2 results (7,018 genes, r=0.77, `docs/FINDINGS_14`).")
    st.warning("This is a **molecular convergence resource**, not a behaviour predictor. Matched "
              "behavioural values (climbing/turning/stumbling) for this panel exist only in the "
              "paper's Figure 1B, not as downloadable data — so this analysis does not, and cannot "
              "yet, predict climbing/lifespan for these disease models.")

# ------------------------------------------------------------- About ----------
else:
    hero("About / methods", "Full method summary, validation evidence, and honest scope "
         "boundaries for both resources shipped in this tool.")
    st.markdown(f"""
**DroPhenoPredict {man['version']}** predicts DGRP organismal phenotypes from genotype and/or
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

**An honest null result, reported for completeness.** This predictor's own GBLUP/multi-modal engine
was applied to a real, matched DGRP-genotype × AD-transgene phenotype (Yang et al. 2023, eye
degeneration from Aβ42+Tau crossed onto 162 DGRP backgrounds) — a direct test of whether this
approach extends beyond the 8 trait-sex models shipped above. **Result: a complete, well-powered
null** (0/5 genotype-only, 0/10 multi-modal, 0/5 inversion-adjusted models significant at FDR<0.05).
Pre-registered risk analysis, replicated phenotype-construction method, and full results in
`docs/FINDINGS_19_yang2023_ad_degeneration.md`. Reported honestly rather than omitted, per this
project's standing rule to report every measured result, positive or negative.
    """)
    st.dataframe(pd.DataFrame([
        {"trait": k, "deployed": m["deployed_modality"], "n": m["n_lines"],
         "adaptive_r": m["adaptive_nested_r"], "predictable": m["predictable"],
         "90%_coverage": m.get("interval_coverage_90")}
        for k, m in man["traits"].items()]), use_container_width=True, hide_index=True)
    st.markdown('<div class="dpp-footer">DroPhenoPredict and DroPhenix are independent, '
                'separately-scoped tools by the same author. Neither shares code, data, or model '
                'artifacts with the other — cross-referenced here only for discoverability.</div>',
                unsafe_allow_html=True)
