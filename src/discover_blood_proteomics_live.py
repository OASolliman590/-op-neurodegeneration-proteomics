#!/usr/bin/env python3
"""
Live discovery of blood proteomics datasets for AD/PD/MS from OmicsDI + PRIDE metadata.

Outputs:
- results/analysis/blood_proteomics_discovery_live.csv
- results/analysis/blood_proteomics_discovery_live_summary.csv
- results/analysis/blood_proteomics_discovery_live.md

Optional copy to T7 via --t7-out-dir.
"""

from __future__ import annotations

import argparse
import csv
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import requests

DISEASE_QUERIES: Dict[str, List[str]] = {
    "Alzheimers": [
        "alzheimer blood proteomics human",
        "alzheimer plasma proteomics human",
        "alzheimer serum proteomics human",
        "alzheimer peripheral blood proteomics human",
        "alzheimer pbmc proteomics human",
    ],
    "Parkinsons": [
        "parkinson blood proteomics human",
        "parkinson plasma proteomics human",
        "parkinson serum proteomics human",
        "parkinson peripheral blood proteomics human",
        "parkinson pbmc proteomics human",
    ],
    "MS": [
        "multiple sclerosis blood proteomics human",
        "multiple sclerosis plasma proteomics human",
        "multiple sclerosis serum proteomics human",
        "multiple sclerosis peripheral blood proteomics human",
        "multiple sclerosis pbmc proteomics human",
    ],
}

BLOOD_TERMS = (
    "blood",
    "whole blood",
    "plasma",
    "serum",
    "pbmc",
    "peripheral blood",
    "mononuclear",
    "leukocyte",
    "lymphocyte",
    "immune cell",
)

CSF_TERMS = ("cerebrospinal fluid", "csf")
BRAIN_TERMS = ("brain", "cortex", "hippocampus", "substantia nigra", "frontal lobe")

PROTEOMICS_TERMS = (
    "proteom",
    "protein",
    "peptide",
    "tmt",
    "lfq",
    "itraq",
    "swath",
    "maxquant",
    "protein group",
    "proteome",
)

NON_PROTEOMICS_TERMS = (
    "metabolom",
    "lipidom",
    "transcriptom",
    "genom",
)

MATRIX_TERMS: List[Tuple[str, Tuple[str, ...]]] = [
    ("plasma", ("plasma",)),
    ("serum", ("serum",)),
    ("pbmc", ("pbmc", "peripheral blood mononuclear", "mononuclear")),
    ("whole_blood", ("whole blood", "blood")),
    ("immune_cells", ("leukocyte", "lymphocyte", "immune cell")),
]

DISEASE_TERMS = {
    "Alzheimers": ("alzheimer", "alzheimers", "dementia"),
    "Parkinsons": ("parkinson", "parkinson's"),
    "MS": ("multiple sclerosis", "rrms", "spms", "ppms"),
}

CASE_TERMS = ("case", "patient", "disease", "ms", "pd", "ad")
CONTROL_TERMS = ("control", "healthy", "normal", "hc", "ctrl")

OMICSDI_SEARCH = "https://www.omicsdi.org/ws/dataset/search"
OMICSDI_DETAIL = "https://www.omicsdi.org/ws/dataset/{source}/{accession}"

# Dataset sources we accept as repository records (exclude literature wrappers).
REPOSITORY_SOURCES = {"pride", "iprox", "jpost", "massive", "panorama"}

# Curated blood-focused seeds to guarantee inclusion even if search ranking/paging misses them.
SEED_DATASETS = [
    {"target_disease": "Alzheimers", "accession": "PXD022265", "source": "pride"},
    {"target_disease": "Alzheimers", "accession": "PXD011482", "source": "pride"},
    {"target_disease": "Alzheimers", "accession": "PXD028392", "source": "pride"},
    {"target_disease": "Alzheimers", "accession": "PXD059280", "source": "pride"},
    {"target_disease": "MS", "accession": "PXD040101", "source": "pride"},
    {"target_disease": "MS", "accession": "PXD035422", "source": "pride"},
    {"target_disease": "MS", "accession": "PXD043337-1", "source": "jpost"},
    {"target_disease": "Parkinsons", "accession": "PXD026439", "source": "iprox"},
    {"target_disease": "Parkinsons", "accession": "PXD051111", "source": "jpost"},
    {"target_disease": "Parkinsons", "accession": "PXD051421-1", "source": "jpost"},
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--out-dir", default=None)
    p.add_argument("--t7-out-dir", default=None)
    p.add_argument("--size", type=int, default=100, help="Results per search page")
    p.add_argument("--max-pages", type=int, default=4, help="Max pages per query")
    p.add_argument("--sleep-seconds", type=float, default=0.2)
    return p.parse_args()


def _norm_text(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, (list, tuple)):
        return " ".join(_norm_text(x) for x in v)
    if isinstance(v, dict):
        return " ".join(_norm_text(x) for x in v.values())
    return str(v)


def _contains_any(text: str, terms: Iterable[str]) -> bool:
    t = text.lower()
    return any(term in t for term in terms)


def _infer_matrix(text: str) -> str:
    t = text.lower()
    for label, terms in MATRIX_TERMS:
        if any(term in t for term in terms):
            return label
    return "unknown"


def _is_human(text: str) -> bool:
    t = text.lower()
    return ("homo sapiens" in t) or ("human" in t)


def _is_disease_match(target_disease: str, text: str) -> bool:
    t = text.lower()
    terms = DISEASE_TERMS[target_disease]
    return any(term in t for term in terms)


def _case_control_hint(text: str) -> bool:
    t = text.lower()
    return _contains_any(t, CASE_TERMS) and _contains_any(t, CONTROL_TERMS)


def _is_repository_source(source: str) -> bool:
    return source.lower() in REPOSITORY_SOURCES


def _is_proteomics_record(source: str, text: str) -> bool:
    t = text.lower()
    has_proteomics = _contains_any(t, PROTEOMICS_TERMS)
    has_non_proteomics = _contains_any(t, NON_PROTEOMICS_TERMS)
    if has_non_proteomics and not has_proteomics:
        return False
    # Repository source + proteomics vocabulary is the strongest signal.
    if _is_repository_source(source) and has_proteomics:
        return True
    # Allow non-standard source labels only if proteomics terms are explicit.
    return has_proteomics


def _score(
    source: str,
    disease_match: bool,
    blood_name_desc: bool,
    blood_meta: bool,
    human_match: bool,
    case_control: bool,
    quant_hint: bool,
    proteomics_hint: bool,
) -> int:
    s = 0
    if _is_repository_source(source):
        s += 2
    if source.lower() == "pride":
        s += 1
    if proteomics_hint:
        s += 2
    if disease_match:
        s += 3
    if blood_meta:
        s += 4
    elif blood_name_desc:
        s += 2
    if human_match:
        s += 2
    if case_control:
        s += 2
    if quant_hint:
        s += 1
    return s


def _search_query(session: requests.Session, query: str, size: int, start: int) -> Dict[str, object]:
    r = session.get(OMICSDI_SEARCH, params={"query": query, "size": size, "start": start}, timeout=45)
    r.raise_for_status()
    return r.json()


def _detail(session: requests.Session, source: str, accession: str) -> Dict[str, object]:
    url = OMICSDI_DETAIL.format(source=source, accession=accession)
    r = session.get(url, params={"detail": "full"}, timeout=45)
    if r.status_code != 200:
        return {}
    try:
        return r.json()
    except Exception:
        return {}


def discover(args: argparse.Namespace) -> List[Dict[str, object]]:
    session = requests.Session()
    session.headers.update({"User-Agent": "op-external-validation/1.0"})

    seen: Dict[Tuple[str, str], Dict[str, object]] = {}

    for disease, queries in DISEASE_QUERIES.items():
        for q in queries:
            for page in range(args.max_pages):
                start = page * args.size
                payload = _search_query(session, q, args.size, start)
                datasets = payload.get("datasets", []) or []
                if not datasets:
                    break

                for ds in datasets:
                    accession = _norm_text(ds.get("id")).strip()
                    source = _norm_text(ds.get("source")).strip().lower() or "unknown"
                    if not accession:
                        continue

                    key = (source, accession)
                    if key not in seen:
                        seen[key] = {
                            "target_disease": disease,
                            "accession": accession,
                            "source": source,
                            "search_title": _norm_text(ds.get("title", "")).strip(),
                            "search_description": _norm_text(ds.get("description", "")).strip(),
                            "matched_queries": set([q]),
                        }
                    else:
                        seen[key]["matched_queries"].add(q)
                        # Keep the first target disease but include all matched diseases in report later.

                if len(datasets) < args.size:
                    break
                time.sleep(args.sleep_seconds)

    # Inject curated seed datasets so known blood cohorts are always evaluated.
    for seed in SEED_DATASETS:
        source = seed["source"].lower()
        accession = seed["accession"]
        key = (source, accession)
        if key in seen:
            seen[key]["matched_queries"].add("seed_list")
            continue
        seen[key] = {
            "target_disease": seed["target_disease"],
            "accession": accession,
            "source": source,
            "search_title": "",
            "search_description": "",
            "matched_queries": set(["seed_list"]),
        }

    out: List[Dict[str, object]] = []

    for (source, accession), base in sorted(seen.items()):
        detail = _detail(session, source, accession)

        name = _norm_text(detail.get("name", "")).strip()
        description = _norm_text(detail.get("description", "")).strip()
        additional = detail.get("additional", {}) if isinstance(detail.get("additional", {}), dict) else {}

        add_tissue = _norm_text(additional.get("tissue", "")).strip()
        add_sample_protocol = _norm_text(additional.get("sample_protocol", "")).strip()
        add_species = _norm_text(additional.get("species", "")).strip()
        add_keywords = _norm_text(additional.get("submitter_keywords", "")).strip()
        add_quant = _norm_text(additional.get("quantification_method", "")).strip()
        add_tech = _norm_text(additional.get("technology_type", "")).strip()
        add_repo = _norm_text(additional.get("repository", "")).strip()
        add_omics_type = _norm_text(additional.get("omics_type", "")).strip()

        text_name_desc = " ".join([base["search_title"], base["search_description"], name, description]).strip()
        text_meta = " ".join([add_tissue, add_sample_protocol, add_keywords, add_quant, add_tech]).strip()
        all_text = " ".join([text_name_desc, text_meta, add_species, add_omics_type]).strip()

        disease_match = _is_disease_match(base["target_disease"], all_text)
        blood_name_desc = _contains_any(text_name_desc, BLOOD_TERMS)
        blood_meta = _contains_any(text_meta, BLOOD_TERMS)
        blood_any = blood_name_desc or blood_meta

        human_match = _is_human(" ".join([add_species, text_name_desc, text_meta]))
        case_control = _case_control_hint(all_text)
        quant_hint = _contains_any(" ".join([add_quant, add_tech, text_name_desc]), ("tmt", "label-free", "lfq", "dia", "dda", "quant", "abundance", "intensity"))
        proteomics_hint = _is_proteomics_record(source, " ".join([text_name_desc, text_meta, add_omics_type]))
        matrix = _infer_matrix(" ".join([add_tissue, add_sample_protocol, text_name_desc, text_meta]))
        csf_mention = _contains_any(all_text, CSF_TERMS)
        brain_mention = _contains_any(all_text, BRAIN_TERMS)
        blood_primary = matrix in {"plasma", "serum", "pbmc", "whole_blood", "immune_cells"} and blood_any
        compartment_mixed = blood_primary and (csf_mention or brain_mention)

        score = _score(
            source=source,
            disease_match=disease_match,
            blood_name_desc=blood_name_desc,
            blood_meta=blood_meta,
            human_match=human_match,
            case_control=case_control,
            quant_hint=quant_hint,
            proteomics_hint=proteomics_hint,
        )

        out.append(
            {
                "target_disease": base["target_disease"],
                "accession": accession,
                "source": source,
                "name": name,
                "title": base["search_title"],
                "description": (description or base["search_description"])[:500],
                "blood_match_name_desc": str(bool(blood_name_desc)).lower(),
                "blood_match_meta": str(bool(blood_meta)).lower(),
                "blood_any": str(bool(blood_any)).lower(),
                "disease_match": str(bool(disease_match)).lower(),
                "human_match": str(bool(human_match)).lower(),
                "case_control_hint": str(bool(case_control)).lower(),
                "quant_hint": str(bool(quant_hint)).lower(),
                "proteomics_hint": str(bool(proteomics_hint)).lower(),
                "repository_source": str(bool(_is_repository_source(source))).lower(),
                "matrix_inferred": matrix,
                "blood_primary": str(bool(blood_primary)).lower(),
                "csf_mention": str(bool(csf_mention)).lower(),
                "brain_mention": str(bool(brain_mention)).lower(),
                "compartment_mixed": str(bool(compartment_mixed)).lower(),
                "species": add_species,
                "tissue": add_tissue,
                "sample_protocol": add_sample_protocol[:500],
                "omics_type": add_omics_type,
                "quantification_method": add_quant,
                "technology_type": add_tech,
                "repository": add_repo,
                "score": score,
                "matched_queries": " | ".join(sorted(base["matched_queries"])),
                "omicsdi_url": f"https://www.omicsdi.org/dataset/{source}/{accession}",
            }
        )

        time.sleep(args.sleep_seconds)

    out.sort(key=lambda r: (r["target_disease"], -int(r["score"]), r["accession"]))
    return out


def summarize(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    by_dis: Dict[str, List[Dict[str, object]]] = {}
    for r in rows:
        by_dis.setdefault(str(r["target_disease"]), []).append(r)

    summary: List[Dict[str, object]] = []
    for disease in sorted(by_dis):
        rs = by_dis[disease]
        def count(field: str, val: str = "true") -> int:
            return sum(1 for x in rs if str(x.get(field, "")).lower() == val)

        summary.append(
            {
                "target_disease": disease,
                "n_datasets_total": len(rs),
                "n_blood_any": count("blood_any"),
                "n_disease_match": count("disease_match"),
                "n_human_match": count("human_match"),
                "n_case_control_hint": count("case_control_hint"),
                "n_proteomics_hint": count("proteomics_hint"),
                "n_repository_source": count("repository_source"),
                "n_blood_primary": count("blood_primary"),
                "n_compartment_mixed": count("compartment_mixed"),
                "n_high_priority": sum(
                    1
                    for x in rs
                    if x.get("blood_any") == "true"
                    and x.get("disease_match") == "true"
                    and x.get("human_match") == "true"
                    and x.get("proteomics_hint") == "true"
                    and x.get("repository_source") == "true"
                    and x.get("blood_primary") == "true"
                    and int(x.get("score", 0)) >= 9
                ),
            }
        )
    return summary


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_md(path: Path, rows: List[Dict[str, object]], summary: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    lines.append("# Live Blood Proteomics Discovery (OmicsDI/PRIDE)")
    lines.append("")
    lines.append("This report inventories blood-related proteomics datasets for AD/PD/MS discovered programmatically.")
    lines.append("")

    lines.append("## Summary")
    lines.append("| disease | total | blood_any | blood_primary | disease_match | human_match | proteomics_hint | repo_source | case_control_hint | mixed_compartment | high_priority |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for s in summary:
        lines.append(
            f"| {s['target_disease']} | {s['n_datasets_total']} | {s['n_blood_any']} | {s['n_blood_primary']} | {s['n_disease_match']} | {s['n_human_match']} | {s['n_proteomics_hint']} | {s['n_repository_source']} | {s['n_case_control_hint']} | {s['n_compartment_mixed']} | {s['n_high_priority']} |"
        )

    lines.append("")
    lines.append("## Top Datasets Per Disease")
    for disease in ["Alzheimers", "Parkinsons", "MS"]:
        top = [
            r for r in rows
            if r["target_disease"] == disease
            and r["blood_any"] == "true"
            and r["blood_primary"] == "true"
            and r["disease_match"] == "true"
            and r["human_match"] == "true"
            and r["proteomics_hint"] == "true"
            and r["repository_source"] == "true"
        ]
        top = sorted(top, key=lambda x: (-int(x["score"]), x["accession"]))[:15]
        lines.append(f"### {disease}")
        if not top:
            lines.append("- none found under current filters")
            continue
        lines.append("| accession | source | matrix | score | case_control | title |")
        lines.append("|---|---|---|---:|---|---|")
        for r in top:
            lines.append(
                f"| {r['accession']} | {r['source']} | {r['matrix_inferred']} | {r['score']} | {r['case_control_hint']} | {r['title'][:80]} |"
            )
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (root / "results" / "analysis")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = discover(args)
    summary = summarize(rows)
    strict_rows = [
        r for r in rows
        if r["blood_any"] == "true"
        and r["blood_primary"] == "true"
        and r["disease_match"] == "true"
        and r["human_match"] == "true"
        and r["proteomics_hint"] == "true"
        and r["repository_source"] == "true"
    ]
    strict_summary = summarize(strict_rows)
    primary_rows = [r for r in strict_rows if r["compartment_mixed"] == "false"]
    primary_summary = summarize(primary_rows)

    data_path = out_dir / "blood_proteomics_discovery_live.csv"
    summary_path = out_dir / "blood_proteomics_discovery_live_summary.csv"
    md_path = out_dir / "blood_proteomics_discovery_live.md"
    strict_data_path = out_dir / "blood_proteomics_discovery_live_strict.csv"
    strict_summary_path = out_dir / "blood_proteomics_discovery_live_strict_summary.csv"
    primary_data_path = out_dir / "blood_proteomics_discovery_live_primary.csv"
    primary_summary_path = out_dir / "blood_proteomics_discovery_live_primary_summary.csv"

    write_csv(data_path, rows)
    write_csv(summary_path, summary)
    write_md(md_path, rows, summary)
    write_csv(strict_data_path, strict_rows)
    write_csv(strict_summary_path, strict_summary)
    write_csv(primary_data_path, primary_rows)
    write_csv(primary_summary_path, primary_summary)

    print(f"Wrote: {data_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {md_path}")
    print(f"Wrote: {strict_data_path}")
    print(f"Wrote: {strict_summary_path}")
    print(f"Wrote: {primary_data_path}")
    print(f"Wrote: {primary_summary_path}")

    if args.t7_out_dir:
        t7_dir = Path(args.t7_out_dir).resolve()
        t7_dir.mkdir(parents=True, exist_ok=True)
        write_csv(t7_dir / data_path.name, rows)
        write_csv(t7_dir / summary_path.name, summary)
        write_md(t7_dir / md_path.name, rows, summary)
        write_csv(t7_dir / strict_data_path.name, strict_rows)
        write_csv(t7_dir / strict_summary_path.name, strict_summary)
        write_csv(t7_dir / primary_data_path.name, primary_rows)
        write_csv(t7_dir / primary_summary_path.name, primary_summary)
        print(f"Copied outputs to: {t7_dir}")


if __name__ == "__main__":
    main()
