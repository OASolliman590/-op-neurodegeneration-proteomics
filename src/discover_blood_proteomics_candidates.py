#!/usr/bin/env python3
"""
discover_blood_proteomics_candidates.py
--------------------------------------
Build a ranked intake list of external blood proteomics datasets for AD/PD/MS
from curated public repository records (PRIDE / ProteomeXchange-hosted).

Outputs:
- config/blood_proteomics_candidates_ranked.csv
- results/analysis/blood_proteomics_candidates_ranked.csv
- results/analysis/blood_proteomics_candidates_summary.md
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    return p.parse_args()


def rows() -> List[Dict[str, str]]:
    return [
        {
            "rank": "1",
            "priority_tier": "A",
            "disease": "Alzheimers",
            "accession": "PXD022265",
            "repository": "PRIDE",
            "hosting_repository": "PRIDE",
            "biospecimen": "blood_plasma",
            "design": "case_control",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "113",
            "quant_ready_likelihood": "medium",
            "notes": "AD plasma biomarker discovery cohort (AA and non-Hispanic-White groups).",
            "source_url": "https://www.omicsdi.org/dataset/pride/PXD022265",
        },
        {
            "rank": "2",
            "priority_tier": "A",
            "disease": "Alzheimers",
            "accession": "PXD011482",
            "repository": "PRIDE",
            "hosting_repository": "PRIDE",
            "biospecimen": "blood_serum",
            "design": "case_control",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "medium",
            "notes": "Deep undepleted human serum proteome profiling for AD biomarker discovery.",
            "source_url": "https://www.omicsdi.org/dataset/pride/PXD011482",
        },
        {
            "rank": "3",
            "priority_tier": "A",
            "disease": "MS",
            "accession": "PXD040101",
            "repository": "PRIDE",
            "hosting_repository": "PRIDE",
            "biospecimen": "blood_plasma",
            "design": "case_control",
            "human_only": "yes",
            "estimated_case_n": "22",
            "estimated_control_n": "22",
            "estimated_total_n": "44",
            "quant_ready_likelihood": "medium",
            "notes": "MS plasma TMT proteomics with explicit MS vs healthy control groups.",
            "source_url": "https://www.omicsdi.org/dataset/pride/PXD040101",
        },
        {
            "rank": "4",
            "priority_tier": "A",
            "disease": "Parkinsons",
            "accession": "PXD026439",
            "repository": "ProteomeXchange",
            "hosting_repository": "iProX",
            "biospecimen": "blood_plasma",
            "design": "case_control",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "PD plasma TMT dataset hosted via iProX (not PRIDE); may need non-PRIDE ingestion route.",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD026439",
        },
        {
            "rank": "5",
            "priority_tier": "B",
            "disease": "Alzheimers",
            "accession": "PXD028392",
            "repository": "PRIDE",
            "hosting_repository": "PRIDE",
            "biospecimen": "blood_serum",
            "design": "case_control",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "128",
            "quant_ready_likelihood": "low",
            "notes": "Serum autoantibody-focused proteomics workflow; supportive but not plain global plasma proteome.",
            "source_url": "https://www.omicsdi.org/dataset/pride/PXD028392",
        },
        {
            "rank": "6",
            "priority_tier": "B",
            "disease": "MS",
            "accession": "PXD035422",
            "repository": "ProteomeXchange",
            "hosting_repository": "PRIDE/jPOST_unspecified",
            "biospecimen": "blood_plasma",
            "design": "longitudinal_relapse_activity",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "RRMS temporal plasma biomarker profiling (activity-linked; not classic simple case-control).",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD035422",
        },
        {
            "rank": "7",
            "priority_tier": "B",
            "disease": "MS",
            "accession": "PXD043337-1",
            "repository": "ProteomeXchange",
            "hosting_repository": "jPOST",
            "biospecimen": "blood_plasma_EV",
            "design": "relapse_vs_remission",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "Plasma EV + CSF EV proteomics in RRMS disease activity context.",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD043337-1&test=no",
        },
        {
            "rank": "8",
            "priority_tier": "B",
            "disease": "Parkinsons",
            "accession": "PXD051421-1",
            "repository": "ProteomeXchange",
            "hosting_repository": "jPOST",
            "biospecimen": "blood_plasma_EV",
            "design": "parkinsonism_subtype",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "Plasma EV proteomics across PD cognitive subtypes, MSA, and healthy controls.",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD051421-1&test=no",
        },
        {
            "rank": "9",
            "priority_tier": "B",
            "disease": "Parkinsons",
            "accession": "PXD051111",
            "repository": "ProteomeXchange",
            "hosting_repository": "jPOST",
            "biospecimen": "blood_plasma_EV",
            "design": "parkinsonism_differential",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "Plasma EV profiling in PD vs MSA differential diagnosis setting.",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD051111",
        },
        {
            "rank": "10",
            "priority_tier": "B",
            "disease": "Alzheimers",
            "accession": "PXD059280",
            "repository": "ProteomeXchange",
            "hosting_repository": "PRIDE/unspecified",
            "biospecimen": "blood_serum_plus_CSF",
            "design": "stage_progression",
            "human_only": "yes",
            "estimated_case_n": "",
            "estimated_control_n": "",
            "estimated_total_n": "",
            "quant_ready_likelihood": "low",
            "notes": "Paired serum/CSF AD progression structural+glycoproteomic profiling; useful for blood layer if serum extractable.",
            "source_url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD059280",
        },
    ]


def write_csv(path: Path, rows_in: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank",
        "priority_tier",
        "disease",
        "accession",
        "repository",
        "hosting_repository",
        "biospecimen",
        "design",
        "human_only",
        "estimated_case_n",
        "estimated_control_n",
        "estimated_total_n",
        "quant_ready_likelihood",
        "notes",
        "source_url",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows_in:
            w.writerow(r)


def write_summary_md(path: Path, rows_in: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_dis: Dict[str, List[Dict[str, str]]] = {}
    for r in rows_in:
        by_dis.setdefault(r["disease"], []).append(r)

    lines: List[str] = []
    lines.append("# Ranked Blood Proteomics Candidate Intake")
    lines.append("")
    lines.append("Curated external candidates for AD/PD/MS blood proteomics expansion.")
    lines.append("")
    lines.append("## Counts")
    for d in sorted(by_dis):
        lines.append(f"- {d}: **{len(by_dis[d])}** candidates")
    lines.append("")
    lines.append("## Ranked List")
    lines.append("| rank | tier | disease | accession | biospecimen | design | host |")
    lines.append("|---:|---|---|---|---|---|---|")
    for r in sorted(rows_in, key=lambda x: int(x["rank"])):
        lines.append(
            f"| {r['rank']} | {r['priority_tier']} | {r['disease']} | {r['accession']} | {r['biospecimen']} | {r['design']} | {r['hosting_repository']} |"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("- Tier A: direct human blood/plasma/serum disease cohorts with likely case-control utility.")
    lines.append("- Tier B: supportive blood layers (EV-focused, progression/relapse designs, mixed-matrix).")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()

    data = rows()

    write_csv(root / "config" / "blood_proteomics_candidates_ranked.csv", data)
    write_csv(root / "results" / "analysis" / "blood_proteomics_candidates_ranked.csv", data)
    write_summary_md(root / "results" / "analysis" / "blood_proteomics_candidates_summary.md", data)

    print(f"Wrote: {root / 'config' / 'blood_proteomics_candidates_ranked.csv'}")
    print(f"Wrote: {root / 'results' / 'analysis' / 'blood_proteomics_candidates_ranked.csv'}")
    print(f"Wrote: {root / 'results' / 'analysis' / 'blood_proteomics_candidates_summary.md'}")


if __name__ == "__main__":
    main()
