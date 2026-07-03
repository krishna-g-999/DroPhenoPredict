"""DroPhenoPredict inference API (supports v1.0 genomic and v1.1 multi-modal).

Returns, for a DGRP line, the predicted trait profile with calibrated 90%
prediction intervals, the data modality used, and an honest predictability flag.

Scope (honest): predictions are leave-one-out GBLUP for lines IN the DGRP panel.
For a novel, unrelated genotype the relatedness-based signal is ~0 (Findings 03),
so the API serves the validated in-panel use case (prioritizing which DGRP lines
to phenotype) rather than over-promising novel-genotype accuracy.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import genotypes, gblup

MODEL_DIR = Path("models/trained")


class DroPhenoPredictor:
    def __init__(self, version: str = "v1.1"):
        self.version = version
        self.manifest = json.loads((MODEL_DIR / f"drophenopredict_{version}.json").read_text())
        self.models = self.manifest["traits"]
        kpath = MODEL_DIR / f"drophenopredict_{version}_kernels.npz"
        if kpath.exists():                       # v1.1 multi-modal: per-model kernels
            self.mode = "multimodal"
            d = np.load(kpath, allow_pickle=True)
            self.arrays = {k: d[k] for k in d.files}
        else:                                    # v1.0 genomic: shared GRM + target vectors
            self.mode = "genomic"
            self.K, self.iids, _ = genotypes.load_grm(self.manifest["grm"])
            self.pos = {l: i for i, l in enumerate(self.iids)}
            td = np.load(MODEL_DIR / f"drophenopredict_{version}_targets.npz", allow_pickle=True)
            self.targets = {k: td[k] for k in td.files if k != "iids"}

    def available_traits(self) -> list[str]:
        return list(self.models.keys())

    @staticmethod
    def normalize_line(line_id: str) -> str:
        """Accept 'line_21', 'DGRP_021', 'DGRP-21', '21' -> 'line_21'."""
        s = str(line_id).strip()
        if s.startswith("line_"):
            return s
        digits = "".join(ch for ch in s if ch.isdigit())
        return f"line_{int(digits)}" if digits else s

    def predict(self, line_id: str, trait_sex: str) -> dict:
        line = self.normalize_line(line_id)
        meta = self.models.get(trait_sex)
        if meta is None:
            raise KeyError(f"Unknown trait '{trait_sex}'. Have: {self.available_traits()}")
        hw = meta["interval_halfwidth_90"]

        if self.mode == "multimodal":
            lines = list(meta["lines"])
            if line not in lines:
                return {"trait": trait_sex, "line": line,
                        "error": "line not covered by this model's data (genotype+expression panel)"}
            y = self.arrays[f"{trait_sex}__y"]
            K = self.arrays[f"{trait_sex}__K"]
            li = lines.index(line)
            train = np.array([i for i in range(len(lines)) if i != li])
            point = float(gblup.gblup_predict(K, y, train, np.array([li]), delta=meta["delta"])[0])
            observed = round(float(y[li]), 3)
            modality = meta["deployed_modality"]
        else:  # genomic v1.0
            if line not in self.pos:
                return {"trait": trait_sex, "line": line, "error": "line not in DGRP panel"}
            y = self.targets[trait_sex]
            obs = ~np.isnan(y); li = self.pos[line]
            train = np.where(obs & (np.arange(len(y)) != li))[0]
            point = float(gblup.gblup_predict(self.K, np.nan_to_num(y), train,
                                              np.array([li]), delta=meta["delta"])[0])
            observed = round(float(y[li]), 3) if obs[li] else None
            modality = "geno"

        # honest accuracy: nested adaptive r (v1.1) or LOO r (v1.0)
        model_r = meta.get("adaptive_nested_r", meta.get("loo_r", meta.get("loo_r_inpanel")))
        out = {
            "trait": trait_sex, "line": line, "unit": meta["unit"], "modality": modality,
            "predicted": round(point, 3),
            "pi90_lower": round(point - hw, 3), "pi90_upper": round(point + hw, 3),
            "predictable": meta["predictable"], "model_cv_r": model_r,
            "confidence": "calibrated genomic/transcriptomic signal" if meta["predictable"]
                          else "not better than trait mean (report as prior only)",
        }
        if observed is not None:
            out["observed"] = observed
        return out

    def predict_profile(self, line_id: str) -> dict:
        return {"line": self.normalize_line(line_id), "model_version": self.version,
                "predictions": [self.predict(line_id, t) for t in self.available_traits()]}
