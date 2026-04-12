"""
query_eprot.py
--------------
Expression Atlas E-PROT proteomics parser for AD/PD brain datasets.

Reads expected chronic directions from config/targets.yaml and writes:
  results/eprot_query.csv
"""

from __future__ import annotations

import argparse
import io
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

from disease_labels import is_case, is_control, log2fc as compute_log2fc, direction_from_fc, concordance as compute_concordance

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config"
OUT = ROOT / "results"

FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

PANEL = [
    "ACTG1",
    "DNAH9",
    "GPX3",
    "VWF",
    "C4B",
    "CD44",
    "CFHR2",
    "ITIH3",
    "LRG1",
    "MYH7B",
]

OUTPUT_COLUMNS = [
    "accession",
    "disease",
    "biospecimen",
    "modality",
    "exposure_type",
    "gene",
    "present",
    "n_case",
    "n_control",
    "mean_case",
    "mean_ctrl",
    "log2fc",
    "direction",
    "expected",
    "concordant",
]

BASELINE_EXPERIMENTS = [
    {"accession": "E-PROT-53", "disease": "Alzheimers", "biospecimen": "brain_DLPFC"},
    {"accession": "E-PROT-56", "disease": "Alzheimers", "biospecimen": "brain_temporal"},
    {"accession": "E-PROT-31", "disease": "Alzheimers", "biospecimen": "brain_atlas"},
    {"accession": "E-PROT-61", "disease": "Alzheimers", "biospecimen": "brain"},
    {"accession": "E-PROT-57", "disease": "Alzheimers", "biospecimen": "brain_MtSinai"},
    {"accession": "E-PROT-65",  "disease": "Parkinsons", "biospecimen": "brain_prefrontal"},
    {"accession": "E-PROT-32",  "disease": "Alzheimers", "biospecimen": "brain_tau_braak_stages"},
    {"accession": "E-PROT-137", "disease": "FTD", "biospecimen": "brain_GRN_MAPT"},
    {"accession": "E-PROT-147", "disease": "FTD", "biospecimen": "brain_celltype"},
    {"accession": "E-PROT-148", "disease": "FTD", "biospecimen": "brain_dentate_gyrus"},
]

DIFF_EXPERIMENT = {
    "accession": "E-PROT-39",
    "disease": "Alzheimers",
    "biospecimen": "brain",
}


def load_targets(path: Path) -> dict[str, dict]:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg.get("targets", {})


def expected_direction(gene: str, targets: dict[str, dict]) -> str:
    t = targets.get(gene, {})
    return t.get("chronic_direction") or "none"


def normalize_col(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def fetch_text(url: str, timeout: int = 90) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


def parse_tsv_without_hash_comments(text: str) -> pd.DataFrame:
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    if not lines:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", low_memory=False)


# direction_from_effect, concordance, is_case, is_control imported from disease_labels


def extract_gene_name_col(df: pd.DataFrame) -> str | None:
    cols = list(df.columns)
    # First try standard gene name columns — prefer ones with non-Ensembl content
    for candidate in ("Gene Name", "Gene.Name"):
        if candidate in cols:
            sample = df[candidate].dropna().astype(str)
            # If values look like gene symbols (not ENSG...), use this column
            if not sample.empty and not sample.str.startswith("ENSG").all():
                return candidate
    # Fall back: look for a column with gene symbol-like values
    for c in cols:
        nc = normalize_col(c)
        if nc in ("genename", "genesymbol", "symbol", "geneid", "gene_id"):
            sample = df[c].dropna().astype(str)
            if not sample.empty and not sample.str.startswith("ENSG").all():
                return c
    # Last resort: Gene Name even if Ensembl
    for candidate in ("Gene Name", "Gene.Name"):
        if candidate in cols:
            return candidate
    return None


def flatten_numeric(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    if df.empty or not cols:
        return np.array([], dtype=float)
    vals = pd.to_numeric(df[cols].to_numpy().ravel(), errors="coerce")
    vals = vals[~np.isnan(vals)]
    return vals


def parse_config_xml(text: str, disease: str) -> tuple[list[str], list[str]]:
    """
    Parse Expression Atlas configuration XML to map group IDs to disease/normal.
    Returns (disease_group_ids, normal_group_ids) as TSV column prefixes (e.g. ['g1', 'g3']).
    """
    disease_ids: list[str] = []
    normal_ids: list[str] = []
    try:
        root = ET.fromstring(text)
        for ag in root.iter("assay_group"):
            gid = ag.get("id", "")
            label = ag.get("label", "").lower()
            if is_case(label, disease):
                disease_ids.append(gid)
            elif is_control(label):
                normal_ids.append(gid)
    except Exception as e:
        log.warning("Failed to parse config XML: %s", e)
    return disease_ids, normal_ids


def process_baseline_experiment(
    exp: dict, targets: dict[str, dict]
) -> list[dict]:
    accession = exp["accession"]
    disease = exp["disease"]
    biospec = exp["biospecimen"]
    rows: list[dict] = []

    config_url = f"{FTP_BASE}/{accession}/{accession}-configuration.xml"
    table_url = f"{FTP_BASE}/{accession}/{accession}.tsv"

    abundance_df = pd.DataFrame()
    disease_cols: list[str] = []
    normal_cols: list[str] = []
    gene_name_col: str | None = None

    # Step 1: parse configuration XML for group-to-disease mapping
    disease_group_ids: list[str] = []
    normal_group_ids: list[str] = []
    try:
        config_text = fetch_text(config_url)
        disease_group_ids, normal_group_ids = parse_config_xml(config_text, disease)
        log.info(
            "%s: config XML parsed — disease groups=%s, normal groups=%s",
            accession, disease_group_ids, normal_group_ids,
        )
    except Exception as e:
        log.warning("%s: failed to fetch/parse config XML (%s)", accession, e)

    # Step 2: load abundance TSV and find matching columns
    try:
        table_text = fetch_text(table_url)
        abundance_df = parse_tsv_without_hash_comments(table_text)
        gene_name_col = extract_gene_name_col(abundance_df)
        if gene_name_col is None and abundance_df.shape[1] >= 2:
            gene_name_col = abundance_df.columns[1]

        # TSV columns are like "g1.WithInSampleAbundance" — match by group id prefix
        for c in abundance_df.columns:
            col_prefix = c.split(".")[0]  # e.g. "g1" from "g1.WithInSampleAbundance"
            if col_prefix in disease_group_ids:
                disease_cols.append(c)
            elif col_prefix in normal_group_ids:
                normal_cols.append(c)

        log.info(
            "%s: table loaded (rows=%d), case_cols=%d, control_cols=%d",
            accession,
            len(abundance_df),
            len(disease_cols),
            len(normal_cols),
        )
    except Exception as e:
        log.warning("%s: failed to parse abundance table (%s)", accession, e)

    for gene in PANEL:
        expected = expected_direction(gene, targets)
        present = False
        mean_case = np.nan
        mean_ctrl = np.nan
        direction = "flat"

        gene_rows = pd.DataFrame()
        if not abundance_df.empty and gene_name_col in abundance_df.columns:
            name_series = abundance_df[gene_name_col].astype(str).str.strip().str.upper()
            gene_rows = abundance_df[name_series == gene.upper()]

        fc = np.nan
        if not gene_rows.empty:
            case_vals = flatten_numeric(gene_rows, disease_cols)
            ctrl_vals = flatten_numeric(gene_rows, normal_cols)

            mean_case = float(np.mean(case_vals)) if case_vals.size else np.nan
            mean_ctrl = float(np.mean(ctrl_vals)) if ctrl_vals.size else np.nan

            if pd.notna(mean_case) and pd.notna(mean_ctrl):
                fc = compute_log2fc(mean_case, mean_ctrl)
                direction = direction_from_fc(fc)

            present = bool(
                pd.notna(mean_case)
                and pd.notna(mean_ctrl)
                and mean_case > 0
                and mean_ctrl > 0
            )

        rows.append(
            {
                "accession": accession,
                "disease": disease,
                "biospecimen": biospec,
                "modality": "proteomics",
                "exposure_type": "chronic",
                "gene": gene,
                "present": present,
                "n_case": len(disease_cols),
                "n_control": len(normal_cols),
                "mean_case": mean_case,
                "mean_ctrl": mean_ctrl,
                "log2fc": round(fc, 4) if pd.notna(fc) else np.nan,
                "direction": direction,
                "expected": expected,
                "concordant": compute_concordance(direction, expected),
            }
        )

    return rows


def process_differential_experiment(
    exp: dict, targets: dict[str, dict]
) -> list[dict]:
    accession = exp["accession"]
    disease = exp["disease"]
    biospec = exp["biospecimen"]
    rows: list[dict] = []

    analytics_url = f"{FTP_BASE}/{accession}/{accession}-analytics.tsv"
    df = pd.DataFrame()

    try:
        text = fetch_text(analytics_url)
        df = parse_tsv_without_hash_comments(text)
        log.info("%s: analytics table loaded (rows=%d)", accession, len(df))
    except Exception as e:
        log.warning("%s: failed to parse analytics TSV (%s)", accession, e)

    gene_col = extract_gene_name_col(df)
    fc_cols = [c for c in df.columns if "log2foldchange" in normalize_col(c)]

    for gene in PANEL:
        expected = expected_direction(gene, targets)
        direction = "flat"
        present = False
        mean_fc = np.nan

        gene_rows = pd.DataFrame()
        if not df.empty and gene_col in df.columns:
            names = df[gene_col].astype(str).str.strip().str.upper()
            gene_rows = df[names == gene.upper()]

        if not gene_rows.empty and fc_cols:
            fc_vals = flatten_numeric(gene_rows, fc_cols)
            if fc_vals.size:
                mean_fc = float(np.mean(fc_vals))
                direction = direction_from_fc(mean_fc)
                present = True

        rows.append(
            {
                "accession": accession,
                "disease": disease,
                "biospecimen": biospec,
                "modality": "proteomics",
                "exposure_type": "chronic",
                "gene": gene,
                "present": present,
                "n_case": len(fc_cols),
                "n_control": len(fc_cols),
                "mean_case": mean_fc,
                "mean_ctrl": 0.0 if pd.notna(mean_fc) else np.nan,
                "log2fc": round(mean_fc, 4) if pd.notna(mean_fc) else np.nan,
                "direction": direction,
                "expected": expected,
                "concordant": compute_concordance(direction, expected),
            }
        )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default=str(CFG / "targets.yaml"))
    parser.add_argument("--out", default=str(OUT / "eprot_query.csv"))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    targets = load_targets(Path(args.targets))

    records: list[dict] = []
    for exp in BASELINE_EXPERIMENTS:
        log.info("Processing baseline %s (%s)", exp["accession"], exp["disease"])
        records.extend(process_baseline_experiment(exp, targets))

    log.info("Processing differential %s", DIFF_EXPERIMENT["accession"])
    records.extend(process_differential_experiment(DIFF_EXPERIMENT, targets))

    out_df = pd.DataFrame(records, columns=OUTPUT_COLUMNS)
    out_df = out_df.reindex(columns=OUTPUT_COLUMNS)
    out_df.to_csv(args.out, index=False)

    log.info("Wrote %d rows to %s", len(out_df), args.out)
    if "present" in out_df.columns:
        log.info("Rows with present=True: %d", int(out_df["present"].fillna(False).sum()))


if __name__ == "__main__":
    main()
