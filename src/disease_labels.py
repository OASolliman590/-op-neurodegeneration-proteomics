"""
disease_labels.py
-----------------
Shared disease/control label detection across Expression Atlas, E-PROT, and PRIDE.

Import from here instead of duplicating keyword logic per script.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Canonical disease keyword mapping — single source of truth
DISEASE_KEYWORDS: dict[str, list[str]] = {
    "Alzheimers": [
        "alzheimer", "alzheimer's", "alzheimers",
        "amyloid", "tau pathology",
    ],
    "Parkinsons": [
        "parkinson", "parkinson's", "parkinsons",
        "lewy body",
    ],
    "MS": [
        "multiple sclerosis", "demyelinating",
    ],
    "FTD": [
        "frontotemporal", "frontotemporal dementia", "ftd",
        "grn", "mapt", "tau pathology", "semantic dementia",
        "progressive supranuclear",
    ],
}

# Control/normal label keywords
CONTROL_KEYWORDS = ["normal", "control", "ctrl", "healthy", "hc"]

# PRIDE column-name heuristics for label inference
PRIDE_CASE_KEYWORDS = {
    "ad", "alzheimer", "parkinson", "pd", "ms", "multiple sclerosis",
    "disease", "affected", "patient",
}
PRIDE_CTRL_KEYWORDS = {"control", "ctrl", "healthy", "normal", "hc"}


def is_case(text: str, disease: str) -> bool:
    """Return True if text matches any keyword for the given disease."""
    t = text.strip().lower()
    return any(kw in t for kw in DISEASE_KEYWORDS.get(disease, []))


def is_control(text: str) -> bool:
    """Return True if text matches any control/normal keyword."""
    t = text.strip().lower()
    return any(kw in t for kw in CONTROL_KEYWORDS)


def assign_diseases(text: str) -> list[str]:
    """Return all disease names whose keywords appear in text (for experiment assignment)."""
    t = text.strip().lower()
    return [d for d, kws in DISEASE_KEYWORDS.items() if any(kw in t for kw in kws)]


def log2fc(mean_case: float, mean_ctrl: float) -> float:
    """Compute log2 fold change; fallback to linear difference if values not both positive."""
    if pd.isna(mean_case) or pd.isna(mean_ctrl):
        return float("nan")
    if mean_case > 0 and mean_ctrl > 0:
        return float(np.log2(mean_case / mean_ctrl))
    return float(mean_case - mean_ctrl)


def direction_from_fc(fc: float, threshold: float = 0.3) -> str:
    if pd.isna(fc):
        return "flat"
    if fc < -threshold:
        return "down"
    if fc > threshold:
        return "up"
    return "flat"


def concordance(direction: str, expected: str) -> str:
    if direction in ("up", "down") and expected not in ("none", ""):
        return "yes" if direction == expected else "no"
    return "na"


def pride_classify_text(text: str, case_lbl: str, ctrl_lbl: str) -> str:
    """
    Classify a PRIDE sample column/label as 'case', 'control', or 'unknown'.
    Prefers explicit cohort labels; falls back to PRIDE_CASE/CTRL_KEYWORDS.
    """
    t = text.lower()
    cl = case_lbl.lower()
    ctl = ctrl_lbl.lower()

    is_ctrl = ctl in t or any(k in t for k in PRIDE_CTRL_KEYWORDS)
    is_dis = cl in t or any(k in t for k in PRIDE_CASE_KEYWORDS if len(k) > 2)

    if is_ctrl:
        return "control"
    if is_dis:
        return "case"
    return "unknown"
