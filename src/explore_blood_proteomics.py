#!/usr/bin/env python3
"""
explore_blood_proteomics.py
---------------------------
Inventory blood proteomics datasets per disease and per dataset from local
project manifests/results.

Outputs (results/analysis):
- blood_proteomics_manifest_inventory.csv
- blood_proteomics_results_inventory.csv
- blood_proteomics_by_disease_summary.csv
- blood_proteomics_exploration.md
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


BLOOD_TERMS = (
    "blood",
    "whole_blood",
    "plasma",
    "serum",
    "pbmc",
    "leukocyte",
    "lymphocyte",
    "immune",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--out-dir", default=None)
    return p.parse_args()


def is_truthy(v: str) -> bool:
    return str(v).strip().lower() in {"1", "true", "t", "yes", "y"}


def has_blood_term(text: str) -> bool:
    t = (text or "").strip().lower()
    return any(term in t for term in BLOOD_TERMS)


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fieldnames))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def inventory_manifest(root: Path, blood_only: bool = True) -> List[Dict[str, object]]:
    cohorts = read_csv_rows(root / "config" / "cohorts.csv")
    out: List[Dict[str, object]] = []

    for r in cohorts:
        modality = (r.get("modality") or "").strip().lower()
        biospec = (r.get("biospecimen") or "").strip()
        notes = (r.get("notes") or "").strip()

        if "proteomics" not in modality:
            continue

        blood_flag = has_blood_term(biospec) or has_blood_term(notes)
        if blood_only and not blood_flag:
            continue

        out.append(
            {
                "accession": r.get("accession", ""),
                "disease": r.get("disease", ""),
                "repository": r.get("repository", ""),
                "biospecimen": biospec,
                "modality": modality,
                "access_type": r.get("access_type", ""),
                "priority": r.get("priority", ""),
                "case_label": r.get("case_label", ""),
                "control_label": r.get("control_label", ""),
                "notes": notes,
                "source": "config/cohorts.csv",
                "blood_flag": str(bool(blood_flag)),
            }
        )
    return out


def summarize_result_file(path: Path, label: str, blood_only: bool = True) -> List[Dict[str, object]]:
    rows = read_csv_rows(path)
    if not rows:
        return []

    groups: Dict[Tuple[str, str, str], Dict[str, object]] = {}

    for r in rows:
        modality = (r.get("modality") or "").strip().lower()
        biospec = (r.get("biospecimen") or "").strip()

        if "proteomics" not in modality:
            continue
        blood_flag = has_blood_term(biospec)
        if blood_only and not blood_flag:
            continue

        acc = (r.get("accession") or "").strip()
        dis = (r.get("disease") or "").strip()
        key = (acc, dis, biospec)

        if key not in groups:
            groups[key] = {
                "accession": acc,
                "disease": dis,
                "biospecimen": biospec,
                "modality": modality,
                "source_file": label,
                "blood_flag": str(bool(blood_flag)),
                "n_rows": 0,
                "genes": set(),
                "present_genes": set(),
                "n_case_max": None,
                "n_control_max": None,
            }

        g = groups[key]
        g["n_rows"] += 1
        gene = (r.get("gene") or "").strip()
        if gene:
            g["genes"].add(gene)
            if is_truthy(r.get("present", "")):
                g["present_genes"].add(gene)

        for fld, out_fld in (("n_case", "n_case_max"), ("n_case_final", "n_case_max"),
                             ("n_control", "n_control_max"), ("n_control_final", "n_control_max")):
            v = (r.get(fld) or "").strip()
            if not v:
                continue
            try:
                val = float(v)
            except ValueError:
                continue
            if g[out_fld] is None or val > g[out_fld]:
                g[out_fld] = val

    out: List[Dict[str, object]] = []
    for g in groups.values():
        out.append(
            {
                "accession": g["accession"],
                "disease": g["disease"],
                "biospecimen": g["biospecimen"],
                "modality": g["modality"],
                "source_file": g["source_file"],
                "blood_flag": g["blood_flag"],
                "n_rows": g["n_rows"],
                "n_genes": len(g["genes"]),
                "n_present_genes": len(g["present_genes"]),
                "n_case": "" if g["n_case_max"] is None else int(g["n_case_max"]),
                "n_control": "" if g["n_control_max"] is None else int(g["n_control_max"]),
            }
        )

    out.sort(key=lambda x: (x["disease"], x["accession"], x["source_file"]))
    return out


def build_disease_summary(manifest_rows: List[Dict[str, object]], result_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    diseases = set([r.get("disease", "") for r in manifest_rows] + [r.get("disease", "") for r in result_rows])
    diseases = sorted(d for d in diseases if d)

    man_by_dis = defaultdict(list)
    for r in manifest_rows:
        man_by_dis[r["disease"]].append(r)

    res_by_dis = defaultdict(list)
    for r in result_rows:
        res_by_dis[r["disease"]].append(r)

    out: List[Dict[str, object]] = []
    for d in diseases:
        m = man_by_dis[d]
        rs = res_by_dis[d]
        out.append(
            {
                "disease": d,
                "n_manifest_blood_proteomics_datasets": len({x["accession"] for x in m}),
                "n_result_blood_proteomics_datasets": len({x["accession"] for x in rs}),
                "n_result_rows": sum(int(x["n_rows"]) for x in rs) if rs else 0,
                "n_result_present_genes_total": sum(int(x["n_present_genes"]) for x in rs) if rs else 0,
                "datasets_manifest": ";".join(sorted({x["accession"] for x in m})),
                "datasets_results": ";".join(sorted({x["accession"] for x in rs})),
            }
        )

    if not out:
        out.append(
            {
                "disease": "NONE",
                "n_manifest_blood_proteomics_datasets": 0,
                "n_result_blood_proteomics_datasets": 0,
                "n_result_rows": 0,
                "n_result_present_genes_total": 0,
                "datasets_manifest": "",
                "datasets_results": "",
            }
        )
    return out


def build_proteomics_context_summary(manifest_all_rows: List[Dict[str, object]], manifest_blood_rows: List[Dict[str, object]], result_all_rows: List[Dict[str, object]], result_blood_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    diseases = sorted(
        set([r.get("disease", "") for r in manifest_all_rows] +
            [r.get("disease", "") for r in result_all_rows] +
            [r.get("disease", "") for r in result_blood_rows])
    )
    diseases = [d for d in diseases if d]

    out: List[Dict[str, object]] = []
    for d in diseases:
        man_all = {r["accession"] for r in manifest_all_rows if r.get("disease") == d}
        man_blood = {r["accession"] for r in manifest_blood_rows if r.get("disease") == d}
        res_all = {r["accession"] for r in result_all_rows if r.get("disease") == d}
        res_blood = {r["accession"] for r in result_blood_rows if r.get("disease") == d}
        out.append(
            {
                "disease": d,
                "manifest_all_proteomics": len(man_all),
                "manifest_blood_proteomics": len(man_blood),
                "results_all_proteomics": len(res_all),
                "results_blood_proteomics": len(res_blood),
            }
        )
    return out


def write_report(path: Path, manifest_rows: List[Dict[str, object]], result_rows: List[Dict[str, object]], summary_rows: List[Dict[str, object]], context_rows: List[Dict[str, object]]) -> None:
    lines: List[str] = []
    lines.append("# Blood Proteomics Exploration")
    lines.append("")
    lines.append("## Scope")
    lines.append("- Inventory blood/plasma/serum/PBMC/immune proteomics datasets by disease and dataset.")
    lines.append("- Sources: `config/cohorts.csv` and local result files.")
    lines.append("")

    lines.append("## Quick Findings")
    lines.append(f"- Manifest blood-proteomics datasets: **{len({r['accession'] for r in manifest_rows})}**")
    lines.append(f"- Results blood-proteomics datasets currently parsed: **{len({r['accession'] for r in result_rows})}**")
    if result_rows:
        diseases = sorted({r.get('disease','') for r in result_rows if r.get('disease')})
        lines.append(f"- Diseases with blood-proteomics result rows: **{', '.join(diseases)}**")
    else:
        lines.append("- Diseases with blood-proteomics result rows: **none yet**")
    lines.append("")

    lines.append("## Per-Disease Summary")
    lines.append("| disease | manifest_datasets | result_datasets | result_rows | present_genes_total |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in summary_rows:
        lines.append(
            f"| {r['disease']} | {r['n_manifest_blood_proteomics_datasets']} | {r['n_result_blood_proteomics_datasets']} | {r['n_result_rows']} | {r['n_result_present_genes_total']} |"
        )

    lines.append("")
    if context_rows:
        lines.append("## Proteomics Context (All vs Blood)")
        lines.append("| disease | manifest_all | manifest_blood | results_all | results_blood |")
        lines.append("|---|---:|---:|---:|---:|")
        for r in context_rows:
            lines.append(
                f"| {r['disease']} | {r['manifest_all_proteomics']} | {r['manifest_blood_proteomics']} | {r['results_all_proteomics']} | {r['results_blood_proteomics']} |"
            )
        lines.append("")

    lines.append("## Output Files")
    lines.append("- `results/analysis/blood_proteomics_manifest_inventory.csv`")
    lines.append("- `results/analysis/blood_proteomics_results_inventory.csv`")
    lines.append("- `results/analysis/blood_proteomics_by_disease_summary.csv`")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (root / "results" / "analysis")
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = inventory_manifest(root, blood_only=True)
    manifest_all_rows = inventory_manifest(root, blood_only=False)

    result_sources = [
        (root / "results" / "eprot_query.csv", "results/eprot_query.csv"),
        (root / "results" / "pride_query.csv", "results/pride_query.csv"),
        (root / "results" / "pride_discovery" / "pride_ad_pd_ms_query.csv", "results/pride_discovery/pride_ad_pd_ms_query.csv"),
        (root / "results" / "analysis" / "publication_effects_with_uncertainty.csv", "results/analysis/publication_effects_with_uncertainty.csv"),
        (root / "results" / "hpa_query.csv", "results/hpa_query.csv"),
    ]

    result_rows: List[Dict[str, object]] = []
    result_all_rows: List[Dict[str, object]] = []
    for p, label in result_sources:
        result_rows.extend(summarize_result_file(p, label, blood_only=True))
        result_all_rows.extend(summarize_result_file(p, label, blood_only=False))

    # de-duplicate exact repeated rows across sources
    dedup_key = set()
    deduped: List[Dict[str, object]] = []
    for r in result_rows:
        k = (r["accession"], r["disease"], r["biospecimen"], r["source_file"])
        if k in dedup_key:
            continue
        dedup_key.add(k)
        deduped.append(r)
    result_rows = deduped

    # de-duplicate all-proteomics summary rows as well
    dedup_all_key = set()
    deduped_all: List[Dict[str, object]] = []
    for r in result_all_rows:
        k = (r["accession"], r["disease"], r["biospecimen"], r["source_file"])
        if k in dedup_all_key:
            continue
        dedup_all_key.add(k)
        deduped_all.append(r)
    result_all_rows = deduped_all

    summary_rows = build_disease_summary(manifest_rows, result_rows)
    context_rows = build_proteomics_context_summary(
        manifest_all_rows=manifest_all_rows,
        manifest_blood_rows=manifest_rows,
        result_all_rows=result_all_rows,
        result_blood_rows=result_rows,
    )

    manifest_fields = [
        "accession", "disease", "repository", "biospecimen", "modality", "access_type",
        "priority", "case_label", "control_label", "notes", "source", "blood_flag",
    ]
    result_fields = [
        "accession", "disease", "biospecimen", "modality", "source_file", "blood_flag", "n_rows",
        "n_genes", "n_present_genes", "n_case", "n_control",
    ]
    summary_fields = [
        "disease", "n_manifest_blood_proteomics_datasets", "n_result_blood_proteomics_datasets",
        "n_result_rows", "n_result_present_genes_total", "datasets_manifest", "datasets_results",
    ]

    write_csv(out_dir / "blood_proteomics_manifest_inventory.csv", manifest_fields, manifest_rows)
    write_csv(out_dir / "blood_proteomics_results_inventory.csv", result_fields, result_rows)
    write_csv(out_dir / "blood_proteomics_by_disease_summary.csv", summary_fields, summary_rows)
    write_csv(out_dir / "blood_proteomics_manifest_all_proteomics.csv", manifest_fields, manifest_all_rows)
    write_csv(out_dir / "blood_proteomics_results_all_proteomics.csv", result_fields, result_all_rows)
    write_report(out_dir / "blood_proteomics_exploration.md", manifest_rows, result_rows, summary_rows, context_rows)

    print(f"Wrote: {out_dir / 'blood_proteomics_manifest_inventory.csv'}")
    print(f"Wrote: {out_dir / 'blood_proteomics_results_inventory.csv'}")
    print(f"Wrote: {out_dir / 'blood_proteomics_by_disease_summary.csv'}")
    print(f"Wrote: {out_dir / 'blood_proteomics_manifest_all_proteomics.csv'}")
    print(f"Wrote: {out_dir / 'blood_proteomics_results_all_proteomics.csv'}")
    print(f"Wrote: {out_dir / 'blood_proteomics_exploration.md'}")


if __name__ == "__main__":
    main()
