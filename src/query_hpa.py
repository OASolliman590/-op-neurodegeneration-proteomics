"""
query_hpa.py
------------
Integrate Human Protein Atlas (HPA) reference datasets into the same schema used
by GEO/PRIDE/Census first-pass validation.

This is a reference layer (not case-control), so:
  - `present` indicates whether the panel gene is detectable in dataset slice
  - `mean_case` stores mean reference expression/concentration
  - `direction` is "no_data" for present rows
  - `expected` is "none", `concordant` is "na"

Output: results/hpa_query.csv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ROOT.parent
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
RAW = EXTERNAL_DATA_ROOT / "hpa_raw"
TARGETS_YAML = ROOT / "config" / "targets.yaml"

DEFAULT_PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]

HPA_BASES = [
    "https://www.proteinatlas.org/download/tsv",
    "https://www.proteinatlas.org/download",
    "https://v25.proteinatlas.org/download/tsv",
    "https://v25.proteinatlas.org/download",
    "https://v24.proteinatlas.org/download/tsv",
    "https://v24.proteinatlas.org/download",
]

# Files correspond to HPA v25 downloadable datasets.
# If a file is not available in a given release, script gracefully emits absent rows.
DATASETS = {
    "RNA_CONSENSUS": {
        "filenames": ["rna_tissue_consensus.tsv.zip"],
        "modality": "transcriptomics",
        "buckets": ["BRAIN", "BLOOD_IMMUNE"],
    },
    "RNA_HPA": {
        "filenames": ["rna_tissue_hpa.tsv.zip"],
        "modality": "transcriptomics",
        "buckets": ["BRAIN", "BLOOD_IMMUNE"],
    },
    "RNA_GTEX": {
        "filenames": ["rna_tissue_gtex.tsv.zip"],
        "modality": "transcriptomics",
        "buckets": ["BRAIN", "BLOOD_IMMUNE"],
    },
    "RNA_BRAIN_HPA": {
        # rna_brain_hpa_subregions.tsv.zip is a brain region hierarchy lookup
        # (only Subregion/Brain region columns — no gene expression data).
        # Use rna_brain_hpa.tsv.zip which contains actual subregion-level expression.
        "filenames": ["rna_brain_hpa.tsv.zip"],
        "modality": "transcriptomics",
        "buckets": ["BRAIN"],
    },
    "RNA_BRAIN_GTEX": {
        "filenames": ["rna_brain_gtex.tsv.zip"],
        "modality": "transcriptomics",
        "buckets": ["BRAIN"],
    },
    "PROTEIN_IHC": {
        "filenames": ["normal_ihc_data.tsv.zip", "normal_tissue.tsv.zip"],
        "modality": "proteomics",
        "buckets": ["BRAIN", "BLOOD_IMMUNE"],
    },
    "BLOOD_MS_CONC": {
        "filenames": ["blood_ms_concentration.tsv.zip"],
        "modality": "proteomics",
        "buckets": ["BLOOD_IMMUNE"],
        "fallback_tissue": "blood",
    },
    "BLOOD_IMMUNOASSAY_CONC": {
        "filenames": ["blood_concentration_immunoassay.tsv.zip"],
        "modality": "proteomics",
        "buckets": ["BLOOD_IMMUNE"],
        "fallback_tissue": "blood",
    },
    "RNA_SINGLE_CELL_TYPE": {
        "filenames": ["rna_single_cell_type.tsv.zip"],
        "modality": "single_cell",
        "buckets": ["BRAIN", "BLOOD_IMMUNE"],
    },
    "RNA_SINGLE_NUCLEI_BRAIN": {
        "filenames": ["rna_single_nuclei_cluster_type.tsv.zip", "rna_single_nuclei_brain.tsv.zip"],
        "modality": "single_cell",
        "buckets": ["BRAIN"],
    },
}

TISSUE_BUCKETS = {
    "BRAIN": [
        "brain", "cortex", "hippocamp", "cerebell", "substantia", "amygdala",
        "thalam", "putamen", "caudate", "hypothalam", "midbrain", "spinal cord",
        "oligodendro", "astrocyte", "microglia", "neuron", "retina",
    ],
    "BLOOD_IMMUNE": [
        "blood", "pbmc", "lymphocyte", "monocyte", "neutrophil", "eosinophil",
        "basophil", "t cell", "b cell", "nk cell", "immune", "plasma",
    ],
}

GENE_COL_CANDS = [
    "Gene name", "gene_name", "Gene symbol", "gene_symbol", "Symbol", "symbol",
    "Gene", "gene",
]

TEXT_COL_CANDS = [
    "Tissue", "tissue", "Blood cell", "blood_cell", "Cell type", "cell_type",
    "Cluster type", "cluster_type", "Cell type class", "cell_type_class",
    "Cluster", "cluster", "Main region", "main_region", "Brain region", "brain_region",
    "Subregion", "subregion", "Region", "region",
]

EXPR_COL_CANDS = [
    "nTPM", "NX", "pTPM", "TPM", "nCPM", "CPM",
    "Expression", "expression", "Level", "level",
    "Concentration", "concentration",
]

LEVEL_MAP = {
    "high": 3.0,
    "medium": 2.0,
    "low": 1.0,
    "not detected": 0.0,
}

GENE_TSV_FIELDS = {
    "RNA_TISSUE_SPECIFIC": {
        "column": "RNA tissue specific nTPM",
        "modality": "transcriptomics",
        "bucket": "BRAIN",
        "accession": "HPA_GENE_TSV_RNA_TISSUE_BRAIN",
        "biospecimen": "brain",
    },
    "RNA_BRAIN_SPECIFIC": {
        "column": "RNA brain regional specific nTPM",
        "modality": "transcriptomics",
        "bucket": "BRAIN",
        "accession": "HPA_GENE_TSV_RNA_BRAIN",
        "biospecimen": "brain",
    },
    "RNA_BLOOD_CELL_SPECIFIC": {
        "column": "RNA blood cell specific nTPM",
        "modality": "transcriptomics",
        "bucket": "BLOOD_IMMUNE",
        "accession": "HPA_GENE_TSV_RNA_BLOOD_IMMUNE",
        "biospecimen": "blood_immune",
    },
    "RNA_SINGLE_CELL_SPECIFIC": {
        "column": "RNA single cell type specific nCPM",
        "modality": "single_cell",
        "bucket": "BRAIN",
        "accession": "HPA_GENE_TSV_SC_CELLTYPE",
        "biospecimen": "brain",
    },
    "RNA_SINGLE_NUCLEI_BRAIN_SPECIFIC": {
        "column": "RNA single nuclei brain specific nCPM",
        "modality": "single_cell",
        "bucket": "BRAIN",
        "accession": "HPA_GENE_TSV_SC_NUCLEI_BRAIN",
        "biospecimen": "brain",
    },
    "BLOOD_MS_CONCENTRATION": {
        "column": "Blood concentration - Conc. blood MS [pg/L]",
        "modality": "proteomics",
        "bucket": "BLOOD_IMMUNE",
        "accession": "HPA_GENE_TSV_BLOOD_MS_CONC",
        "biospecimen": "blood_immune",
    },
    "BLOOD_IM_CONCENTRATION": {
        "column": "Blood concentration - Conc. blood IM [pg/L]",
        "modality": "proteomics",
        "bucket": "BLOOD_IMMUNE",
        "accession": "HPA_GENE_TSV_BLOOD_IMMUNOASSAY_CONC",
        "biospecimen": "blood_immune",
    },
}


def load_panel(targets_yaml: Path) -> list[str]:
    if not targets_yaml.exists():
        return DEFAULT_PANEL

    try:
        with open(targets_yaml, "r", encoding="utf-8") as fh:
            cfg = yaml.safe_load(fh) or {}
        targets = cfg.get("targets", {})
        panel = [str(g).upper() for g in targets.keys()]
        return panel if panel else DEFAULT_PANEL
    except Exception as e:
        log.warning(f"Failed reading targets.yaml ({e}); using default panel")
        return DEFAULT_PANEL


def _first(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _best_gene_col(df: pd.DataFrame, panel: list[str]) -> Optional[str]:
    existing = [c for c in GENE_COL_CANDS if c in df.columns]
    if not existing:
        return None

    panel_set = set(panel)
    best_col, best_score = None, -1
    for c in existing:
        vals = (
            df[c]
            .astype(str)
            .str.strip()
            .str.upper()
            .replace({"NAN": "", "NONE": ""})
        )
        score = int(vals.isin(panel_set).sum())
        if score > best_score:
            best_score, best_col = score, c

    # If nothing matches panel genes, prefer symbol-like columns over Ensembl IDs.
    if best_score <= 0:
        for c in existing:
            if "name" in c.lower() or "symbol" in c.lower():
                return c

    return best_col


def ensure_hpa_file(filenames: list[str], raw_dir: Path) -> Optional[Path]:
    raw_dir.mkdir(parents=True, exist_ok=True)

    def _looks_valid_download(path: Path) -> bool:
        if not path.exists() or path.stat().st_size < 1024:
            return False
        if "".join(path.suffixes).endswith(".zip"):
            return zipfile.is_zipfile(path)
        return True

    for fn in filenames:
        local = raw_dir / fn
        if _looks_valid_download(local):
            return local
        if local.exists():
            try:
                local.unlink()
            except Exception:
                pass

    for fn in filenames:
        out = raw_dir / fn
        for base in HPA_BASES:
            url = f"{base}/{fn}"
            try:
                with requests.get(
                    url,
                    stream=True,
                    timeout=90,
                    headers={
                        "User-Agent": (
                            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                            "AppleWebKit/537.36 (KHTML, like Gecko) "
                            "Chrome/124.0.0.0 Safari/537.36"
                        ),
                        "Referer": "https://www.proteinatlas.org/about/download",
                    },
                ) as r:
                    if r.status_code != 200:
                        continue
                    with open(out, "wb") as fh:
                        for chunk in r.iter_content(1 << 20):
                            if chunk:
                                fh.write(chunk)
                log.info(f"Downloaded HPA file: {url}")
                return out
            except Exception:
                continue
    return None


def read_hpa_table(path: Path) -> Optional[pd.DataFrame]:
    try:
        if "".join(path.suffixes).endswith(".tsv.zip"):
            return pd.read_csv(path, sep="\t", compression="zip", low_memory=False)
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:
        log.warning(f"Cannot parse HPA table {path.name}: {e}")
        return None


def parse_expr(v) -> float:
    if pd.isna(v):
        return np.nan
    if isinstance(v, (int, float, np.number)):
        return float(v)

    s = str(v).strip()
    sl = s.lower()
    if sl in LEVEL_MAP:
        return LEVEL_MAP[sl]

    # Handle values like "12.3 pg/L".
    m = re.search(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", s)
    if m:
        try:
            return float(m.group(0))
        except Exception:
            pass

    return np.nan


def _best_expr_col(df: pd.DataFrame) -> Optional[str]:
    for c in EXPR_COL_CANDS:
        if c in df.columns:
            return c

    skip = set(GENE_COL_CANDS + TEXT_COL_CANDS)
    candidates = [c for c in df.columns if c not in skip]
    best_col, best_ratio = None, 0.0
    for c in candidates:
        parsed = df[c].map(parse_expr)
        ratio = float(parsed.notna().mean())
        if ratio > best_ratio:
            best_ratio, best_col = ratio, c
    return best_col if best_ratio > 0.2 else None


def to_long_table(df: pd.DataFrame,
                  panel: list[str],
                  fallback_tissue: Optional[str] = None) -> pd.DataFrame:
    gene_col = _best_gene_col(df, panel)
    if gene_col is None:
        return pd.DataFrame(columns=["gene", "tissue_text", "expr"])

    text_cols = [c for c in TEXT_COL_CANDS if c in df.columns]
    expr_col = _best_expr_col(df)

    out = pd.DataFrame()
    out["gene"] = df[gene_col].astype(str).str.upper().str.strip()

    if text_cols:
        parts = [df[c].astype(str).str.lower().str.strip() for c in text_cols]
        merged = parts[0]
        for p in parts[1:]:
            merged = merged + " " + p
        out["tissue_text"] = merged
    else:
        out["tissue_text"] = (fallback_tissue or "all").lower()

    out["expr"] = df[expr_col].map(parse_expr) if expr_col else np.nan
    out = out.dropna(subset=["gene"])
    out = out[out["gene"].isin(set(panel))]
    return out


def bucket_mask(tissue_series: pd.Series, keywords: list[str]) -> pd.Series:
    pat = "|".join([re.escape(k).replace(r"\ ", r"\s+") for k in keywords])
    return tissue_series.str.contains(pat, case=False, na=False, regex=True)


def empty_rows_for_dataset(dataset_key: str,
                           modality: str,
                           buckets: list[str],
                           panel: list[str]) -> list[dict]:
    rows = []
    for bucket in buckets:
        acc = f"HPA_{dataset_key}_{bucket}"
        for gene in panel:
            rows.append(dict(
                accession=acc,
                disease="HPA_reference",
                biospecimen=bucket.lower(),
                modality=modality,
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
    return rows


def build_bucket_rows(dataset_key: str,
                      modality: str,
                      long_df: pd.DataFrame,
                      buckets: list[str],
                      panel: list[str]) -> list[dict]:
    if long_df.empty:
        return empty_rows_for_dataset(dataset_key, modality, buckets, panel)

    rows = []
    for bucket in buckets:
        acc = f"HPA_{dataset_key}_{bucket}"
        keywords = TISSUE_BUCKETS.get(bucket, [])
        if keywords:
            bdf = long_df[bucket_mask(long_df["tissue_text"], keywords)].copy()
        else:
            bdf = long_df.copy()

        for gene in panel:
            gsub = bdf[bdf["gene"] == gene]
            present = not gsub.empty
            mean_expr = float(gsub["expr"].mean()) if present and gsub["expr"].notna().any() else np.nan
            rows.append(dict(
                accession=acc,
                disease="HPA_reference",
                biospecimen=bucket.lower(),
                modality=modality,
                exposure_type="not_applicable",
                gene=gene,
                present=present,
                n_case=0,
                n_control=0,
                mean_case=round(mean_expr, 4) if not np.isnan(mean_expr) else None,
                mean_ctrl=None,
                direction="no_data" if present else "absent",
                expected="none",
                concordant="na",
            ))

    return rows


def _mean_from_hpa_specific_field(v) -> float:
    if pd.isna(v):
        return np.nan
    s = str(v).strip()
    if not s:
        return np.nan

    # Formats like: "bone marrow: 49.9;lymphoid tissue: 70.9"
    vals = []
    for tok in s.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            tok = tok.split(":", 1)[1].strip()
        x = parse_expr(tok)
        if not np.isnan(x):
            vals.append(x)

    if not vals:
        # Could still be a single scalar string
        x = parse_expr(s)
        return x if not np.isnan(x) else np.nan

    return float(np.mean(vals))


def _collect_gene_tsv_paths(gene_tsv_dir: Optional[Path],
                            gene_tsv_files: list[str]) -> list[Path]:
    paths: list[Path] = []

    for fp in gene_tsv_files:
        p = Path(fp).expanduser().resolve()
        if p.exists() and p.is_file():
            paths.append(p)

    if gene_tsv_dir and gene_tsv_dir.exists():
        paths.extend(sorted(gene_tsv_dir.glob("ENSG*.tsv")))

    # de-duplicate while keeping order
    seen = set()
    dedup = []
    for p in paths:
        if p not in seen:
            dedup.append(p)
            seen.add(p)
    return dedup


def build_rows_from_gene_tsvs(panel: list[str],
                              gene_tsv_dir: Optional[Path],
                              gene_tsv_files: list[str]) -> list[dict]:
    paths = _collect_gene_tsv_paths(gene_tsv_dir, gene_tsv_files)
    if not paths:
        return []

    gene_values: dict[str, dict[str, float]] = {}

    for p in paths:
        try:
            gdf = pd.read_csv(p, sep="\t", low_memory=False)
            if gdf.empty:
                continue
            row = gdf.iloc[0]
            gene = str(row.get("Gene", "")).upper().strip()
            if gene not in set(panel):
                continue

            if gene not in gene_values:
                gene_values[gene] = {}

            for key, meta in GENE_TSV_FIELDS.items():
                col = meta["column"]
                if col in gdf.columns:
                    gene_values[gene][key] = _mean_from_hpa_specific_field(row.get(col))
        except Exception as e:
            log.warning(f"Cannot parse local gene TSV {p}: {e}")

    rows = []
    panel_set = set(panel)
    for key, meta in GENE_TSV_FIELDS.items():
        acc = meta["accession"]
        modality = meta["modality"]
        biospecimen = meta["biospecimen"]

        for gene in panel:
            val = gene_values.get(gene, {}).get(key, np.nan)
            present = bool(not np.isnan(val))
            rows.append(dict(
                accession=acc,
                disease="HPA_reference",
                biospecimen=biospecimen,
                modality=modality,
                exposure_type="not_applicable",
                gene=gene,
                present=present,
                n_case=0,
                n_control=0,
                mean_case=round(float(val), 4) if present else None,
                mean_ctrl=None,
                direction="no_data" if present else "absent",
                expected="none",
                concordant="na",
            ))

    detected = len([g for g in panel_set if g in gene_values])
    if detected:
        log.info(f"Local ENSG TSV support: found panel data for {detected}/{len(panel)} genes")

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(OUT / "hpa_query.csv"))
    parser.add_argument("--raw-dir", default=str(RAW))
    parser.add_argument("--datasets", nargs="*", choices=list(DATASETS.keys()),
                        help="Optional subset of HPA dataset keys")
    parser.add_argument("--gene-tsv-dir", default=str(PROJECT_ROOT),
                        help="Directory containing optional ENSG*.tsv gene files")
    parser.add_argument("--gene-tsv-files", nargs="*", default=[],
                        help="Optional explicit paths to ENSG gene TSV files")
    parser.add_argument("--no-gene-tsv", action="store_true",
                        help="Disable local ENSG*.tsv integration")
    parser.add_argument("--download-only", action="store_true",
                        help="Only download HPA files into --raw-dir and exit")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = Path(args.raw_dir)
    panel = load_panel(TARGETS_YAML)

    dataset_keys = args.datasets if args.datasets else list(DATASETS.keys())
    all_rows = []
    download_rows = []

    for key in dataset_keys:
        cfg = DATASETS[key]
        p = ensure_hpa_file(cfg["filenames"], raw_dir)
        download_rows.append({
            "dataset_key": key,
            "requested_files": ";".join(cfg["filenames"]),
            "status": "ok" if p is not None else "missing",
            "resolved_path": str(p) if p is not None else "",
        })

        if args.download_only:
            continue

        if p is None:
            log.warning(f"{key}: file unavailable, writing absent rows")
            all_rows.extend(empty_rows_for_dataset(key, cfg["modality"], cfg.get("buckets", []), panel))
            continue

        df = read_hpa_table(p)
        if df is None or df.empty:
            log.warning(f"{key}: unreadable/empty table, writing absent rows")
            all_rows.extend(empty_rows_for_dataset(key, cfg["modality"], cfg.get("buckets", []), panel))
            continue

        long_df = to_long_table(
            df,
            panel=panel,
            fallback_tissue=cfg.get("fallback_tissue"),
        )
        log.info(f"{key}: parsed {len(long_df)} panel-matched rows from {p.name}")
        all_rows.extend(
            build_bucket_rows(
                dataset_key=key,
                modality=cfg["modality"],
                long_df=long_df,
                buckets=cfg.get("buckets", []),
                panel=panel,
            )
        )

    if not args.no_gene_tsv:
        gene_tsv_rows = build_rows_from_gene_tsvs(
            panel=panel,
            gene_tsv_dir=Path(args.gene_tsv_dir),
            gene_tsv_files=args.gene_tsv_files,
        )
        all_rows.extend(gene_tsv_rows)

    if args.download_only:
        manifest = pd.DataFrame(download_rows)
        manifest_path = raw_dir / "hpa_download_manifest.csv"
        manifest.to_csv(manifest_path, index=False)
        ok = int((manifest["status"] == "ok").sum())
        log.info(f"Download-only complete: {ok}/{len(manifest)} datasets resolved")
        log.info(f"Download manifest: {manifest_path}")
        return

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(args.out, index=False)
    log.info(f"Saved {len(out_df)} rows to {args.out}")


if __name__ == "__main__":
    main()
