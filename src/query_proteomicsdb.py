"""
query_proteomicsdb.py
---------------------
Query the live ProteomicsDB public REST API for a panel of genes and extract
expression evidence in neurologic and biofluid tissues. If ProteomicsDB is not
usable, fall back to UniProt tissue-specificity annotations.

Output: results/proteomicsdb_query.csv
Columns:
  gene, uniprot_id, tissue, disease_context, n_samples,
  mean_abundance, detected, source, notes
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]

UNIPROT_IDS = {
    "ACTG1": "P63261",
    "DNAH9": "Q9ULW0",
    "GPX3": "P22352",
    "VWF": "P04275",
    "C4B": "P0C0L4",
    "CD44": "P16070",
    "CFHR2": "P36980",
    "ITIH3": "Q06033",
    "LRG1": "P02750",
    "MYH7B": "Q8WZ82",
}

OUTPUT_COLUMNS = [
    "gene",
    "uniprot_id",
    "tissue",
    "disease_context",
    "n_samples",
    "mean_abundance",
    "detected",
    "source",
    "notes",
]

BASE = "https://www.proteomicsdb.org/api/v2/proteomicsdb/biology"

RELEVANT_TISSUE_TERMS = [
    "brain",
    "cortex",
    "hippocampus",
    "substantia nigra",
    "frontal",
    "temporal",
    "cerebellum",
    "cerebrospinal",
    "csf",
    "plasma",
    "serum",
    "blood",
]

DISEASE_TERMS = ["alzheimer", "parkinson", "multiple sclerosis", "ftd"]


def shorten(text: str, max_len: int = 200) -> str:
    clean = re.sub(r"\s+", " ", text or "").strip()
    return clean[:max_len]


def build_url(url: str, params: dict[str, Any] | None = None) -> str:
    req = requests.Request("GET", url, params=params).prepare()
    return req.url or url


def norm_key(key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(key).lower())


def get_ci(record: dict[str, Any], candidates: list[str]) -> Any:
    norm_map = {norm_key(k): v for k, v in record.items()}
    for key in candidates:
        nk = norm_key(key)
        if nk in norm_map:
            value = norm_map[nk]
            if value not in (None, "", [], {}):
                return value
    return None


def is_numeric(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        return re.fullmatch(r"[-+]?\d+(?:\.\d+)?", value.strip()) is not None
    return False


def to_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def to_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def collect_dicts(obj: Any) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def _walk(node: Any) -> None:
        if isinstance(node, dict):
            out.append(node)
            for value in node.values():
                _walk(value)
        elif isinstance(node, list):
            for item in node:
                _walk(item)

    _walk(obj)
    return out


def safe_get(
    session: requests.Session,
    url: str,
    params: dict[str, Any] | None,
    timeout: int,
    logs: list[dict[str, Any]],
    gene: str,
    endpoint: str,
) -> dict[str, Any]:
    full_url = build_url(url, params)
    try:
        response = session.get(url, params=params, timeout=timeout)
        preview = shorten(response.text, max_len=200)
        logs.append({
            "gene": gene,
            "endpoint": endpoint,
            "url": full_url,
            "status_code": response.status_code,
            "preview": preview,
        })
        print(
            f"URL: {full_url}\n"
            f"status: {response.status_code}\n"
            f"preview: {preview}\n"
        )
        payload = None
        try:
            payload = response.json()
        except json.JSONDecodeError:
            payload = None
        return {
            "ok": True,
            "status_code": response.status_code,
            "url": full_url,
            "text": response.text,
            "json": payload,
        }
    except requests.RequestException as exc:
        err_preview = shorten(str(exc), max_len=200)
        logs.append({
            "gene": gene,
            "endpoint": endpoint,
            "url": full_url,
            "status_code": None,
            "preview": err_preview,
        })
        print(
            f"URL: {full_url}\n"
            f"status: ERROR\n"
            f"preview: {err_preview}\n"
        )
        return {
            "ok": False,
            "status_code": None,
            "url": full_url,
            "text": "",
            "json": None,
            "error": str(exc),
        }


def parse_service_endpoints(payload: Any) -> list[str]:
    endpoints: list[str] = []
    if isinstance(payload, dict):
        value = payload.get("value")
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str):
                    endpoints.append(item)
                elif isinstance(item, dict):
                    candidate = item.get("url") or item.get("name") or item.get("endpoint")
                    if candidate:
                        endpoints.append(str(candidate))
        d_obj = payload.get("d")
        if isinstance(d_obj, dict):
            entity_sets = d_obj.get("EntitySets")
            if isinstance(entity_sets, list):
                endpoints.extend(str(x) for x in entity_sets)

    deduped: list[str] = []
    seen: set[str] = set()
    for endpoint in endpoints:
        cleaned = endpoint.strip().strip("/")
        if cleaned and cleaned not in seen:
            deduped.append(cleaned)
            seen.add(cleaned)
    return deduped


def is_relevant_tissue(text: str) -> bool:
    text_l = text.lower()
    return any(term in text_l for term in RELEVANT_TISSUE_TERMS)


def extract_relevant_tissues(text: str) -> list[str]:
    text_l = text.lower()
    matches = [term for term in RELEVANT_TISSUE_TERMS if term in text_l]
    # Preserve order while deduping
    return list(dict.fromkeys(matches))


def has_disease_context(text: str) -> bool:
    text_l = text.lower()
    return any(term in text_l for term in DISEASE_TERMS)


def parse_protein_info(payload: Any, gene: str) -> tuple[str | None, str | None]:
    if payload is None:
        return None, None
    records = collect_dicts(payload)
    for record in records:
        rec_gene = str(get_ci(record, ["gene_name", "gene", "symbol"]) or "").upper()
        if rec_gene and rec_gene != gene.upper():
            continue
        protein_id = get_ci(record, ["protein_id", "proteinid", "id"])
        uniprot_id = get_ci(
            record,
            ["uniprot_id", "uniprotid", "uniprot", "accession", "uniprot_accession"],
        )
        protein_id_s = str(protein_id).strip() if protein_id not in (None, "") else None
        uniprot_id_s = str(uniprot_id).strip() if uniprot_id not in (None, "") else None
        if protein_id_s or uniprot_id_s:
            return protein_id_s, uniprot_id_s
    return None, None


def guess_mean_abundance(record: dict[str, Any]) -> float | None:
    direct = get_ci(
        record,
        [
            "mean_abundance",
            "abundance",
            "mean",
            "expression",
            "expression_value",
            "normalized_abundance",
            "intensity",
            "value",
        ],
    )
    value = to_float(direct)
    if value is not None:
        return value

    for key, raw in record.items():
        key_l = key.lower()
        if any(tag in key_l for tag in ("abundance", "expression", "intensity", "mean", "value")):
            if is_numeric(raw):
                numeric = to_float(raw)
                if numeric is not None:
                    return numeric
    return None


def extract_expression_rows(
    payload: Any,
    gene: str,
    uniprot_id: str | None,
    endpoint: str,
    protein_id: str | None,
) -> list[dict[str, Any]]:
    if payload is None:
        return []

    rows: list[dict[str, Any]] = []
    for record in collect_dicts(payload):
        tissue = get_ci(
            record,
            [
                "tissue",
                "tissue_name",
                "tissue_name_long",
                "organ",
                "organ_name",
                "sample_type",
                "bio_source",
                "biosource",
            ],
        )
        description = get_ci(
            record,
            [
                "description",
                "details",
                "comment",
                "disease",
                "condition",
                "experiment_name",
                "study_name",
                "title",
                "context",
            ],
        )
        tissue_s = str(tissue).strip() if tissue is not None else ""
        description_s = str(description).strip() if description is not None else ""
        combined = f"{tissue_s} {description_s}".strip()
        if not combined or not is_relevant_tissue(combined):
            continue

        n_samples = to_int(get_ci(record, ["n_samples", "sample_count", "number_of_samples", "nsamples"]))
        mean_abundance = guess_mean_abundance(record)
        disease_context = has_disease_context(description_s)
        notes = [f"endpoint={endpoint}"]
        if protein_id:
            notes.append(f"protein_id={protein_id}")
        if description_s:
            notes.append(f"description={shorten(description_s, 120)}")

        output_tissue = tissue_s if tissue_s and is_relevant_tissue(tissue_s) else extract_relevant_tissues(combined)[0]
        rows.append({
            "gene": gene,
            "uniprot_id": uniprot_id or UNIPROT_IDS.get(gene, ""),
            "tissue": output_tissue,
            "disease_context": disease_context,
            "n_samples": n_samples,
            "mean_abundance": mean_abundance,
            "detected": True,
            "source": "ProteomicsDB",
            "notes": "; ".join(notes),
        })
    return rows


def extract_texts_from_comment(comment: dict[str, Any]) -> list[str]:
    texts: list[str] = []
    for text_obj in comment.get("texts", []):
        if isinstance(text_obj, dict):
            value = text_obj.get("value")
            if value:
                texts.append(str(value))
    note = comment.get("note")
    if isinstance(note, dict):
        for text_obj in note.get("texts", []):
            if isinstance(text_obj, dict):
                value = text_obj.get("value")
                if value:
                    texts.append(str(value))
    return texts


def query_uniprot_fallback(
    session: requests.Session,
    gene: str,
    uniprot_id: str,
    logs: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    resp = safe_get(
        session=session,
        url=url,
        params=None,
        timeout=30,
        logs=logs,
        gene=gene,
        endpoint="uniprot_fallback",
    )
    if not resp["ok"] or resp["status_code"] is None or resp["status_code"] >= 400:
        notes = f"UniProt request failed: {shorten(resp.get('error', 'HTTP error'), 150)}"
        return [{
            "gene": gene,
            "uniprot_id": uniprot_id,
            "tissue": "",
            "disease_context": False,
            "n_samples": None,
            "mean_abundance": None,
            "detected": False,
            "source": "UniProt_fallback",
            "notes": notes,
        }]

    payload = resp.get("json")
    if not isinstance(payload, dict):
        return [{
            "gene": gene,
            "uniprot_id": uniprot_id,
            "tissue": "",
            "disease_context": False,
            "n_samples": None,
            "mean_abundance": None,
            "detected": False,
            "source": "UniProt_fallback",
            "notes": "UniProt response was not JSON object.",
        }]

    tissue_texts: list[str] = []
    subcellular_locations: list[str] = []
    for comment in payload.get("comments", []):
        if not isinstance(comment, dict):
            continue
        ctype = str(comment.get("commentType", "")).upper()
        if ctype == "TISSUE SPECIFICITY":
            tissue_texts.extend(extract_texts_from_comment(comment))
        if ctype == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                if not isinstance(loc, dict):
                    continue
                location = loc.get("location")
                if isinstance(location, dict):
                    value = location.get("value")
                    if value:
                        subcellular_locations.append(str(value))

    tissue_blob = " ".join(tissue_texts)
    matched_tissues = extract_relevant_tissues(tissue_blob)
    disease_context = has_disease_context(tissue_blob)
    subcell = ", ".join(sorted(set(subcellular_locations)))

    if matched_tissues:
        rows = []
        for tissue in matched_tissues:
            notes = "UniProt tissue specificity match"
            if subcell:
                notes += f"; subcellular={shorten(subcell, 120)}"
            rows.append({
                "gene": gene,
                "uniprot_id": uniprot_id,
                "tissue": tissue,
                "disease_context": disease_context,
                "n_samples": None,
                "mean_abundance": None,
                "detected": True,
                "source": "UniProt_fallback",
                "notes": notes,
            })
        return rows

    notes = "No relevant tissues found in UniProt tissue specificity."
    if subcell:
        notes += f" subcellular={shorten(subcell, 120)}"
    return [{
        "gene": gene,
        "uniprot_id": uniprot_id,
        "tissue": "",
        "disease_context": disease_context,
        "n_samples": None,
        "mean_abundance": None,
        "detected": False,
        "source": "UniProt_fallback",
        "notes": notes,
    }]


def dedupe_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[Any, ...]] = set()
    deduped: list[dict[str, Any]] = []
    for row in rows:
        key = (
            row.get("gene"),
            row.get("uniprot_id"),
            row.get("tissue"),
            row.get("disease_context"),
            row.get("n_samples"),
            row.get("mean_abundance"),
            row.get("detected"),
            row.get("source"),
            row.get("notes"),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(OUT / "proteomicsdb_query.csv"))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    api_logs: list[dict[str, Any]] = []

    # Step 1: Explore ProteomicsDB API service root.
    print("=== STEP 1: Explore ProteomicsDB API ===")
    root_resp = safe_get(
        session=session,
        url=f"{BASE}/",
        params=None,
        timeout=30,
        logs=api_logs,
        gene="GLOBAL",
        endpoint="biology_root",
    )
    endpoints = parse_service_endpoints(root_resp.get("json"))
    print("Available endpoints from /biology/:")
    if endpoints:
        for endpoint in endpoints:
            print(endpoint)
    else:
        print("No endpoints parsed from response.")
    print()

    all_rows: list[dict[str, Any]] = []

    # Step 2: Query candidate endpoints per gene.
    for gene in PANEL:
        print(f"=== Gene: {gene} ===")
        protein_resp = safe_get(
            session=session,
            url=f"{BASE}/protein/",
            params={"gene_name": gene, "no_isoforms": 1},
            timeout=30,
            logs=api_logs,
            gene=gene,
            endpoint="protein",
        )
        protein_id, uniprot_from_pdb = parse_protein_info(protein_resp.get("json"), gene)
        uniprot_id = uniprot_from_pdb or UNIPROT_IDS.get(gene, "")

        expression_calls = [
            (
                "proteinexpression",
                f"{BASE}/proteinexpression/",
                {
                    "protein_id": protein_id or "",
                    "experiment_type": "MS",
                    "scope": "proteomicsdb",
                },
            ),
            (
                "humanproteinatlas",
                f"{BASE}/humanproteinatlas/",
                {"gene_name": gene},
            ),
            (
                "expressionprofile",
                f"{BASE}/expressionprofile/",
                {"protein_id": protein_id or ""},
            ),
            (
                "tissueexpression",
                f"{BASE}/tissueexpression/",
                {"gene_name": gene},
            ),
        ]

        gene_rows: list[dict[str, Any]] = []
        for endpoint_name, url, params in expression_calls:
            resp = safe_get(
                session=session,
                url=url,
                params=params,
                timeout=30,
                logs=api_logs,
                gene=gene,
                endpoint=endpoint_name,
            )
            if resp.get("status_code") is not None and 200 <= resp["status_code"] < 400:
                gene_rows.extend(
                    extract_expression_rows(
                        payload=resp.get("json"),
                        gene=gene,
                        uniprot_id=uniprot_id,
                        endpoint=endpoint_name,
                        protein_id=protein_id,
                    )
                )

        if gene_rows:
            all_rows.extend(dedupe_rows(gene_rows))
        else:
            # Temporary placeholder, replaced below if global fallback is triggered.
            all_rows.append({
                "gene": gene,
                "uniprot_id": uniprot_id or UNIPROT_IDS.get(gene, ""),
                "tissue": "",
                "disease_context": False,
                "n_samples": None,
                "mean_abundance": None,
                "detected": False,
                "source": "ProteomicsDB",
                "notes": f"No relevant ProteomicsDB expression rows. protein_id={protein_id or 'NA'}",
            })
        print()

    # Step 3: Fallback logic.
    pdb_calls = [
        log for log in api_logs
        if log.get("endpoint") in {
            "biology_root",
            "protein",
            "proteinexpression",
            "humanproteinatlas",
            "expressionprofile",
            "tissueexpression",
        }
    ]
    pdb_any_2xx = any(
        isinstance(log.get("status_code"), int) and 200 <= int(log["status_code"]) < 300
        for log in pdb_calls
    )

    if not pdb_any_2xx:
        print("=== STEP 3: ProteomicsDB unavailable -> UniProt fallback for all genes ===")
        all_rows = []
        for gene in PANEL:
            uniprot_id = UNIPROT_IDS[gene]
            fallback_rows = query_uniprot_fallback(
                session=session,
                gene=gene,
                uniprot_id=uniprot_id,
                logs=api_logs,
            )
            all_rows.extend(fallback_rows)
    else:
        print("=== STEP 3: ProteomicsDB provided at least one 2xx response ===")
        # Optionally fill per-gene no-data rows with UniProt context.
        genes_with_detected = {row["gene"] for row in all_rows if row.get("source") == "ProteomicsDB" and row.get("detected")}
        for gene in PANEL:
            if gene in genes_with_detected:
                continue
            fallback_rows = query_uniprot_fallback(
                session=session,
                gene=gene,
                uniprot_id=UNIPROT_IDS[gene],
                logs=api_logs,
            )
            all_rows.extend(fallback_rows)

    all_rows = dedupe_rows(all_rows)
    df = pd.DataFrame(all_rows, columns=OUTPUT_COLUMNS)
    df.to_csv(output_path, index=False)

    logs_df = pd.DataFrame(api_logs, columns=["gene", "endpoint", "url", "status_code", "preview"])
    logs_df.to_csv(OUT / "proteomicsdb_api_calls.csv", index=False)

    print(f"\nSaved {len(df)} rows to {output_path}")
    print(f"Saved API call log to {OUT / 'proteomicsdb_api_calls.csv'}")

    detected_df = df[df["detected"] == True]
    if detected_df.empty:
        print("\nPanel genes detected per tissue: none")
    else:
        print("\nPanel genes detected per tissue:")
        by_tissue = (
            detected_df.groupby("tissue")["gene"]
            .apply(lambda x: ", ".join(sorted(set(x))))
            .sort_index()
        )
        for tissue, genes in by_tissue.items():
            print(f"- {tissue}: {genes}")

    print("\nFull output CSV:")
    print(df.to_csv(index=False))


if __name__ == "__main__":
    main()
