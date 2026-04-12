"""
guardrail_check.py
------------------
Automated QA gate for external validation outputs.

Usage:
  python src/guardrail_check.py
  python src/guardrail_check.py --strict
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
QA = RESULTS / "qa"
CONFIG = ROOT / "config"

PANEL_SIZE = 10
REQUIRED_COLUMNS = [
    "accession", "disease", "biospecimen", "modality", "exposure_type",
    "gene", "present", "n_case", "n_control", "mean_case", "mean_ctrl",
    "direction", "expected", "concordant",
]

SOURCE_FILES = {
    "geo": "geo_query.csv",
    "pride": "pride_query.csv",
    "hpa": "hpa_query.csv",
    "expression_atlas": "expression_atlas_query.csv",
    "proteomicsdb": "proteomicsdb_query.csv",
    "opentargets": "opentargets_query.csv",
    "diseases": "diseases_query.csv",
    "master": "master_query.csv",
}


@dataclass
class Check:
    name: str
    status: str  # PASS/WARN/FAIL
    detail: str


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def source_expected_rows(cohorts: pd.DataFrame, repo: str) -> int:
    sub = cohorts[(cohorts["repository"] == repo) & (cohorts["access_type"] == "open")]
    return len(sub) * PANEL_SIZE


def run_checks(strict: bool) -> tuple[list[Check], bool]:
    checks: list[Check] = []
    hard_fail = False

    cohorts = pd.read_csv(CONFIG / "cohorts.csv")

    data = {k: read_csv(RESULTS / v) for k, v in SOURCE_FILES.items()}

    # 1) Required files
    for key, fname in SOURCE_FILES.items():
        p = RESULTS / fname
        if not p.exists():
            checks.append(Check(f"file_exists::{fname}", "FAIL", "missing"))
            hard_fail = True
        elif p.stat().st_size == 0:
            checks.append(Check(f"file_exists::{fname}", "FAIL", "empty file"))
            hard_fail = True
        else:
            checks.append(Check(f"file_exists::{fname}", "PASS", f"{p.stat().st_size} bytes"))

    # 2) Schema checks
    for key in ["geo", "pride", "hpa", "expression_atlas", "proteomicsdb", "opentargets", "diseases"]:
        df = data[key]
        if df.empty:
            checks.append(Check(f"schema::{key}", "FAIL", "cannot parse or empty dataframe"))
            hard_fail = True
            continue
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            checks.append(Check(f"schema::{key}", "FAIL", f"missing columns: {missing}"))
            hard_fail = True
        else:
            checks.append(Check(f"schema::{key}", "PASS", "required columns present"))

    # 3) Row counts (hard for GEO/PRIDE)
    geo_rows = len(data["geo"])
    pride_rows = len(data["pride"])
    exp_geo = source_expected_rows(cohorts, "GEO")
    exp_pride = source_expected_rows(cohorts, "PRIDE")

    if geo_rows == exp_geo:
        checks.append(Check("rows::geo", "PASS", f"{geo_rows}/{exp_geo}"))
    else:
        checks.append(Check("rows::geo", "FAIL", f"{geo_rows}/{exp_geo}"))
        hard_fail = True

    if pride_rows == exp_pride:
        checks.append(Check("rows::pride", "PASS", f"{pride_rows}/{exp_pride}"))
    else:
        checks.append(Check("rows::pride", "FAIL", f"{pride_rows}/{exp_pride}"))
        hard_fail = True

    # 4) Multiple-of-panel guardrail
    for key in ["geo", "pride", "hpa", "expression_atlas", "proteomicsdb", "opentargets", "diseases"]:
        n = len(data[key])
        if n == 0 or n % PANEL_SIZE != 0:
            checks.append(Check(f"panel_multiple::{key}", "FAIL", f"rows={n}"))
            hard_fail = True
        else:
            checks.append(Check(f"panel_multiple::{key}", "PASS", f"rows={n}"))

    # 5) Critical label guardrails
    geo = data["geo"]
    if not geo.empty:
        critical_geo = cohorts[
            (cohorts["repository"] == "GEO")
            & (cohorts["access_type"] == "open")
            & (cohorts["disease"] == "OP_exposure")
        ]["accession"].drop_duplicates().tolist()
        for acc in critical_geo:
            sub = geo[geo["accession"] == acc]
            if sub.empty:
                checks.append(Check(f"labels::{acc}", "FAIL", "no rows"))
                hard_fail = True
                continue
            n_case = pd.to_numeric(sub["n_case"], errors="coerce").fillna(0).max()
            n_ctrl = pd.to_numeric(sub["n_control"], errors="coerce").fillna(0).max()
            if n_case > 0 and n_ctrl > 0:
                checks.append(Check(f"labels::{acc}", "PASS", f"case={int(n_case)}, ctrl={int(n_ctrl)}"))
            else:
                checks.append(Check(f"labels::{acc}", "FAIL", f"case={int(n_case)}, ctrl={int(n_ctrl)}"))
                hard_fail = True

    # 6) Soft checks
    # PRIDE useful cohorts
    pride = data["pride"]
    if not pride.empty:
        coh = pride[["accession", "n_case", "n_control"]].drop_duplicates().copy()
        coh["n_case"] = pd.to_numeric(coh["n_case"], errors="coerce").fillna(0)
        coh["n_control"] = pd.to_numeric(coh["n_control"], errors="coerce").fillna(0)
        good = int(((coh["n_case"] > 0) & (coh["n_control"] > 0)).sum())
        if good >= 2:
            checks.append(Check("soft::pride_useful_cohorts", "PASS", f"{good} cohorts"))
        else:
            checks.append(Check("soft::pride_useful_cohorts", "WARN", f"{good} cohorts"))

    # Expression atlas all-absent
    ea = data["expression_atlas"]
    if not ea.empty:
        present = int(pd.to_numeric(ea["present"], errors="coerce").fillna(0).astype(bool).sum())
        if present > 0:
            checks.append(Check("soft::expression_atlas_present", "PASS", f"present_rows={present}"))
        else:
            checks.append(Check("soft::expression_atlas_present", "WARN", "all absent"))

    # ProteomicsDB/PeptideAtlas all-absent
    pdb = data["proteomicsdb"]
    if not pdb.empty:
        present = int(pd.to_numeric(pdb["present"], errors="coerce").fillna(0).astype(bool).sum())
        if present > 0:
            checks.append(Check("soft::proteomicsdb_present", "PASS", f"present_rows={present}"))
        else:
            checks.append(Check("soft::proteomicsdb_present", "WARN", "all absent"))

    # Scoreable volume trend sanity
    master = data["master"]
    if not master.empty and "concordant" in master.columns:
        scoreable = int(master["concordant"].isin(["yes", "no"]).sum())
        if scoreable >= 80:
            checks.append(Check("soft::scoreable_volume", "PASS", f"scoreable_rows={scoreable}"))
        else:
            checks.append(Check("soft::scoreable_volume", "WARN", f"scoreable_rows={scoreable}"))

    if strict:
        # In strict mode, WARN does not fail automatically, only hard FAILs.
        pass

    return checks, hard_fail


def write_report(checks: list[Check], failed: bool) -> Path:
    QA.mkdir(parents=True, exist_ok=True)
    out = QA / "guardrail_report.md"
    lines = []
    lines.append("# Guardrail Report")
    lines.append("")
    lines.append(f"Overall: {'FAIL' if failed else 'PASS'}")
    lines.append("")
    lines.append("| Check | Status | Detail |")
    lines.append("|---|---|---|")
    for c in checks:
        lines.append(f"| {c.name} | {c.status} | {c.detail} |")
    out.write_text("\n".join(lines))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", action="store_true", help="Fail on hard guardrail failures")
    args = parser.parse_args()

    checks, hard_fail = run_checks(strict=args.strict)
    report = write_report(checks, hard_fail)
    print(f"Guardrail report: {report}")

    n_pass = sum(1 for c in checks if c.status == "PASS")
    n_warn = sum(1 for c in checks if c.status == "WARN")
    n_fail = sum(1 for c in checks if c.status == "FAIL")
    print(f"PASS={n_pass} WARN={n_warn} FAIL={n_fail}")

    if args.strict and hard_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
