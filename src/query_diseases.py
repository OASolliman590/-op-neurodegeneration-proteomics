"""
query_diseases.py
-----------------
Query Jensen Lab DISEASES bulk downloads for disease-gene evidence
(used as an open DisGeNET-like evidence layer).

Sources (filtered channels, human):
  - human_disease_knowledge_filtered.tsv
  - human_disease_experiments_filtered.tsv
  - optional: human_disease_textmining_filtered.tsv

This layer is association evidence (not directional fold-change), so:
  - `present` means the gene has disease association evidence in channel
  - `mean_case` stores mean confidence score (0-1 where available)
  - `direction` = "no_data" when present
  - `expected` = "none", `concordant` = "na"

Output: results/diseases_query.csv
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

PREFERRED_EXTERNAL_ROOT = Path("/Volumes/T7/5-Alzhimers_Parkisons_MS_External_Valida/op_external_validation_data")


def choose_external_root() -> Path:
    env_root = os.getenv("OP_EXTERNAL_DATA_ROOT")
    if env_root:
        return Path(env_root).expanduser()
    if PREFERRED_EXTERNAL_ROOT.exists():
        return PREFERRED_EXTERNAL_ROOT
    return ROOT / "data" / "external"


EXTERNAL_DATA_ROOT = choose_external_root()
RAW = EXTERNAL_DATA_ROOT / "diseases_raw"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]

OUTPUT_COLUMNS = [
    "accession", "disease", "biospecimen", "modality", "exposure_type",
    "gene", "present", "n_case", "n_control", "mean_case", "mean_ctrl",
    "direction", "expected", "concordant",
]

BASE_URL = "https://download.jensenlab.org"

CHANNELS = {
    "knowledge": {
        "filename": "human_disease_knowledge_filtered.tsv",
        "accession": "DISEASES_KNOWLEDGE",
    },
    "experiments": {
        "filename": "human_disease_experiments_filtered.tsv",
        "accession": "DISEASES_EXPERIMENTS",
    },
    "textmining": {
        "filename": "human_disease_textmining_filtered.tsv",
        "accession": "DISEASES_TEXTMINING",
    },
}

DISEASE_KEYWORDS = {
    "Alzheimers": ["alzheimer", "alzheimers", "alzheimer's"],
    "Parkinsons": ["parkinson", "parkinsons", "parkinson's"],
    "MS": ["multiple sclerosis"],
}


def ensure_file(filename: str, raw_dir: Path) -> Optional[Path]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    path = raw_dir / filename
    if path.exists() and path.stat().st_size > 1024:
        return path

    url = f"{BASE_URL}/{filename}"
    try:
        with requests.get(url, stream=True, timeout=180) as r:
            r.raise_for_status()
            with open(path, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    if chunk:
                        fh.write(chunk)
        log.info(f"Downloaded {filename}")
        return path
    except Exception as e:
        log.warning(f"Failed downloading {filename}: {e}")
        return None


def _norm_cols(df: pd.DataFrame) -> dict[str, str]:
    return {c.lower().strip().replace(" ", "_"): c for c in df.columns}


def _resolve_columns(df: pd.DataFrame) -> tuple[Optional[str], Optional[str], Optional[str]]:
    cols = _norm_cols(df)

    gene_col = None
    for key in ("gene_name", "gene", "symbol"):
        if key in cols:
            gene_col = cols[key]
            break
    if gene_col is None and len(df.columns) >= 2:
        gene_col = df.columns[1]

    disease_col = None
    for key in ("disease_name", "disease"):
        if key in cols:
            disease_col = cols[key]
            break
    if disease_col is None and len(df.columns) >= 4:
        disease_col = df.columns[3]

    conf_col = None
    for key in ("confidence", "confidence_score", "score"):
        if key in cols:
            conf_col = cols[key]
            break
    if conf_col is None:
        # fallback: choose the last mostly numeric column
        best_col = None
        best_ratio = 0.0
        for c in df.columns:
            ratio = pd.to_numeric(df[c], errors="coerce").notna().mean()
            if ratio > best_ratio:
                best_ratio, best_col = ratio, c
        if best_ratio > 0.4:
            conf_col = best_col

    return gene_col, disease_col, conf_col


def _disease_match_mask(series: pd.Series, disease_name: str) -> pd.Series:
    kws = DISEASE_KEYWORDS.get(disease_name, [])
    text = series.astype(str).str.lower()
    if not kws:
        return pd.Series([False] * len(series), index=series.index)
    mask = pd.Series(False, index=series.index)
    for kw in kws:
        mask = mask | text.str.contains(kw, na=False, regex=False)
    return mask


def rows_for_channel(df: pd.DataFrame, accession: str) -> list[dict]:
    gene_col, disease_col, conf_col = _resolve_columns(df)
    if gene_col is None or disease_col is None:
        log.warning(f"{accession}: unable to resolve gene/disease columns")
        return []

    work = df.copy()
    work[gene_col] = work[gene_col].astype(str).str.upper().str.strip()
    work = work[work[gene_col].isin(set(PANEL))]

    if conf_col is not None:
        work["_confidence"] = pd.to_numeric(work[conf_col], errors="coerce")
    else:
        work["_confidence"] = np.nan

    records = []
    for disease in DISEASE_KEYWORDS:
        dmask = _disease_match_mask(work[disease_col], disease)
        dsub = work[dmask]

        for gene in PANEL:
            gsub = dsub[dsub[gene_col] == gene]
            present = not gsub.empty
            n_hits = int(len(gsub))
            mean_conf = float(gsub["_confidence"].mean()) if present and gsub["_confidence"].notna().any() else np.nan

            records.append(dict(
                accession=accession,
                disease=disease,
                biospecimen="multi",
                modality="disease_association",
                exposure_type="not_applicable",
                gene=gene,
                present=present,
                n_case=n_hits,
                n_control=0,
                mean_case=round(mean_conf, 4) if not np.isnan(mean_conf) else None,
                mean_ctrl=None,
                direction="no_data" if present else "absent",
                expected="none",
                concordant="na",
            ))

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(OUT / "diseases_query.csv"))
    parser.add_argument("--raw-dir", default=str(RAW))
    parser.add_argument("--channels", nargs="*", choices=list(CHANNELS.keys()),
                        default=["knowledge", "experiments"],
                        help="DISEASES channels to include")
    parser.add_argument("--include-textmining", action="store_true",
                        help="Include textmining channel (larger file)")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = Path(args.raw_dir)

    channels = list(args.channels)
    if args.include_textmining and "textmining" not in channels:
        channels.append("textmining")

    all_records = []
    for ch in channels:
        cfg = CHANNELS[ch]
        p = ensure_file(cfg["filename"], raw_dir)
        if p is None:
            log.warning(f"{ch}: source unavailable, writing absent rows")
            for disease in DISEASE_KEYWORDS:
                for gene in PANEL:
                    all_records.append(dict(
                        accession=cfg["accession"],
                        disease=disease,
                        biospecimen="multi",
                        modality="disease_association",
                        exposure_type="not_applicable",
                        gene=gene,
                        present=False,
                        n_case=0,
                        n_control=0,
                        mean_case=None,
                        mean_ctrl=None,
                        direction="absent",
                        expected="none",
                        concordant="na",
                    ))
            continue

        try:
            df = pd.read_csv(p, sep="\t", low_memory=False)
        except Exception as e:
            log.warning(f"{ch}: failed to parse {p.name}: {e}")
            continue

        rows = rows_for_channel(df, cfg["accession"])
        log.info(f"{cfg['accession']}: generated {len(rows)} rows")
        all_records.extend(rows)

    out_df = pd.DataFrame(all_records, columns=OUTPUT_COLUMNS)
    out_df.to_csv(args.out, index=False)
    log.info(f"Saved {len(out_df)} rows to {args.out}")


if __name__ == "__main__":
    main()
