"""
download_all.py
---------------
Single entry point to download ALL open-access datasets to the T7 SSD.

Downloads to:
  /Volumes/T7/.../geo_raw/<accession>/   — GEO expression + metadata
  /Volumes/T7/.../pride_raw/<accession>/ — PRIDE protein group tables
  /Volumes/T7/.../hpa_raw/               — HPA reference files (already done)

Skips anything already fully cached.
Writes a manifest CSV at the end: data/download_manifest.csv

Usage:
  python3 src/download_all.py                  # all sources
  python3 src/download_all.py --only geo
  python3 src/download_all.py --only pride
  python3 src/download_all.py --status         # show manifest only
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("/tmp/op_download.log"),
    ],
)
log = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).resolve().parents[1]
CFG       = ROOT / "config"
T7        = Path("/Volumes/T7/5-Alzhimers_Parkisons_MS_External_Valida/op_external_validation_data")
GEO_RAW   = T7 / "geo_raw"
PRIDE_RAW = T7 / "pride_raw"
HPA_RAW   = T7 / "hpa_raw"
MANIFEST  = ROOT / "data" / "download_manifest.csv"

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2"

GENE_COL_CANDIDATES = [
    "Gene Symbol", "gene_symbol", "GENE_SYMBOL", "Symbol",
    "Gene_Symbol", "gene symbol", "Gene symbol", "Gene Symbol (HGNC)",
]
PROBE_COL_CANDIDATES = ["ID", "ID_REF", "Probe ID", "ProbeID", "probe_id"]
INTENSITY_PREFIXES   = [
    "LFQ intensity ", "Intensity ", "iBAQ ",
    "Abundance: ", "Reporter intensity ", "MS2 ",
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_t7() -> bool:
    if not T7.exists():
        log.error(f"T7 SSD not mounted at {T7}. Please connect it and retry.")
        return False
    for d in [GEO_RAW, PRIDE_RAW, HPA_RAW]:
        d.mkdir(parents=True, exist_ok=True)
    return True


def load_cohorts(modality_filter: list[str] | None = None) -> pd.DataFrame:
    df = pd.read_csv(CFG / "cohorts.csv")
    df = df[df["access_type"] == "open"]
    if modality_filter:
        df = df[df["repository"].isin(modality_filter)]
    return df


def _first_col(cols: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _parse_gene(val) -> str | None:
    if pd.isna(val):
        return None
    s = re.split(r"\s*///\s*|[;|]+", str(val).strip())[0].strip()
    return s or None


def record(acc: str, status: str, path: str = "", notes: str = "") -> dict:
    return dict(
        accession=acc,
        status=status,
        path=path,
        notes=notes,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )


# ── GEO download + process ────────────────────────────────────────────────────

def geo_is_cached(acc: str) -> bool:
    d = GEO_RAW / acc
    return (d / f"{acc}_expression.parquet").exists() and \
           (d / f"{acc}_metadata.csv").exists()


def process_geo(acc: str, case_label: str, ctrl_label: str) -> dict:
    out_dir = GEO_RAW / acc
    out_dir.mkdir(parents=True, exist_ok=True)

    parquet_path  = out_dir / f"{acc}_expression.parquet"
    metadata_path = out_dir / f"{acc}_metadata.csv"

    if parquet_path.exists() and metadata_path.exists():
        log.info(f"  {acc}: already cached — skipping")
        return record(acc, "cached", str(parquet_path))

    # Expression parquet already computed — just redo metadata
    parquet_only = parquet_path.exists() and not metadata_path.exists()

    log.info(f"  {acc}: {'loading cached parquet, rebuilding metadata' if parquet_only else 'downloading from GEO'} …")
    try:
        gse = GEOparse.get_GEO(geo=acc, destdir=str(out_dir), silent=True)
    except Exception as e:
        log.error(f"  {acc}: GEO download failed — {e}")
        return record(acc, "failed", notes=str(e))

    if not parquet_only:
        # Build expression matrix
        frames = []
        for gsm_name, gsm in gse.gsms.items():
            if gsm.table is None or gsm.table.empty:
                continue
            tbl = gsm.table.copy()
            id_col  = _first_col(tbl.columns.tolist(),
                                 ["ID_REF", "ID", "ProbeID"]) or tbl.columns[0]
            val_col = _first_col(tbl.columns.tolist(),
                                 ["VALUE", "value", "Signal"]) or (
                                 tbl.columns[1] if len(tbl.columns) > 1 else None)
            if val_col is None:
                continue
            s = tbl[[id_col, val_col]].copy()
            s.columns = ["probe_id", gsm_name]
            s[gsm_name] = pd.to_numeric(s[gsm_name], errors="coerce")
            frames.append(s.set_index("probe_id"))

        if not frames:
            # RNA-seq fallback: look for a series-level count matrix in supplementary files
            suppl_urls = gse.metadata.get("supplementary_file", [])
            if isinstance(suppl_urls, str):
                suppl_urls = [suppl_urls]
            expr = None
            for url in suppl_urls:
                url = url.strip()
                # GEO FTP URLs need to be converted to HTTPS
                if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
                    url = "https://ftp.ncbi.nlm.nih.gov/" + url[len("ftp://ftp.ncbi.nlm.nih.gov/"):]
                fname = url.split("/")[-1]
                # Only handle single flat count/expression tables (not tar archives)
                if any(k in fname.lower() for k in ["count", "expression", "tpm", "fpkm"]) \
                   and any(fname.endswith(e) for e in [".tsv.gz", ".csv.gz", ".txt.gz", ".tsv", ".csv"]) \
                   and "raw.tar" not in fname.lower():
                    dest = out_dir / fname
                    if not dest.exists():
                        log.info(f"  {acc}: downloading supplementary count matrix {fname} …")
                        try:
                            with requests.get(url, stream=True, timeout=180) as resp:
                                resp.raise_for_status()
                                with open(dest, "wb") as fh:
                                    for chunk in resp.iter_content(1 << 20):
                                        if chunk:
                                            fh.write(chunk)
                        except Exception as e:
                            log.warning(f"  {acc}: supplementary download failed — {e}")
                            continue
                    try:
                        sep = "\t" if fname.endswith(".tsv.gz") or fname.endswith(".tsv") else ","
                        mat = pd.read_csv(dest, sep=sep, index_col=0, low_memory=False)
                        # Drop non-numeric columns
                        mat = mat.select_dtypes(include="number")
                        if not mat.empty:
                            expr = mat
                            expr.index = expr.index.astype(str)
                            log.info(f"  {acc}: loaded supplementary matrix {mat.shape}")
                            break
                    except Exception as e:
                        log.warning(f"  {acc}: failed to parse {fname} — {e}")

            if expr is None:
                log.warning(f"  {acc}: no expression tables in GSMs or supplementary files")
                return record(acc, "failed", notes="no expression tables")

        if not frames and expr is not None:
            # Skip probe→gene mapping (supplementary matrices are already gene-level)
            expr.to_parquet(parquet_path)
            # Jump straight to metadata
            n_genes = expr.shape[0]
        else:
            expr = pd.concat(frames, axis=1)
            expr = expr.loc[~expr.index.duplicated(keep="first")]
            expr.index = expr.index.astype(str)

        # Probe → gene map (only for array data)
        probe_map = {}
        for gpl in gse.gpls.values():
            tbl = gpl.table
            if tbl is None or tbl.empty:
                continue
            gene_col  = _first_col(tbl.columns.tolist(), GENE_COL_CANDIDATES)
            probe_col = _first_col(tbl.columns.tolist(), PROBE_COL_CANDIDATES) or tbl.columns[0]
            if gene_col is None:
                continue
            for _, row in tbl[[probe_col, gene_col]].dropna().iterrows():
                sym = _parse_gene(row[gene_col])
                if sym:
                    probe_map[str(row[probe_col])] = sym

        # Collapse probes → genes
        if probe_map:
            expr["gene_symbol"] = expr.index.map(probe_map)
            expr = expr.dropna(subset=["gene_symbol"])
            expr = expr.groupby("gene_symbol").mean(numeric_only=True)
        else:
            log.warning(f"  {acc}: no probe map — saving at probe level")

        expr.to_parquet(parquet_path)
    else:
        expr = pd.read_parquet(parquet_path)

    # Build metadata (infer case/control from GSM metadata)
    case_terms = {case_label.lower(), case_label.lower().replace("_", " ")}
    ctrl_terms = {ctrl_label.lower(), ctrl_label.lower().replace("_", " ")}
    if case_label.upper() == "AD":
        case_terms |= {"alzheimer", "alzheimers"}
    if case_label.upper() == "PD":
        case_terms |= {"parkinson", "parkinsons"}
    if case_label.upper() == "MS":
        case_terms |= {"multiple sclerosis"}
    if "chlorpyrifos" in case_label.lower() or "op" in case_label.lower():
        case_terms |= {"chlorpyrifos", "exposed", "organophosphate",
                       # Common OP ester compound names in GEO metadata
                       "bpdp", "ippp", "tboep", "tmpp", "tphp", "tcpp", "tdcpp",
                       "tbep", "tcep", "tris", "tipp", "triphenyl phosphate"}
    if "unexposed" in ctrl_label.lower() or ctrl_label.lower() in ("control", "controls"):
        ctrl_terms |= {"vehicle", "untreated", "unexposed", "healthy", "normal",
                       "ctl", "ctrl", "hc", "nc", "hcs", "dmso", "solvent", "control"}

    rows = []
    for gsm_name, gsm in gse.gsms.items():
        meta = gsm.metadata or {}

        # 1. Try to extract explicit "status: X" value from characteristics_ch1
        chars = meta.get("characteristics_ch1", [])
        if isinstance(chars, str):
            chars = [chars]
        status_val = None
        for c in (chars or []):
            # Accept both "status: X" and "treatment: X" as explicit discriminators
            m = re.match(r"(?:status|treatment):\s*(.+)", str(c).strip(), re.IGNORECASE)
            if m:
                status_val = m.group(1).strip().lower()
                break

        if status_val is not None:
            # Primary: match within status/treatment field (word boundary, not full-string)
            has_case = any(re.search(rf"\b{re.escape(t)}\b", status_val) for t in case_terms)
            has_ctrl = any(re.search(rf"\b{re.escape(t)}\b", status_val) for t in ctrl_terms)
        else:
            # Fallback: full metadata text search (exclude generic characteristics)
            fields = ["title", "source_name_ch1", "description"]
            text = " | ".join(
                (v if isinstance(v, str) else " ".join(str(x) for x in (v or [])))
                for k in fields for v in [meta.get(k, "")]
            ).lower()
            # Also include characteristics but strip "case-control" to avoid false positives
            chars_text = re.sub(r"case.?control", "", " ".join(str(c) for c in (chars or []))).lower()
            full_text = text + " | " + chars_text
            has_case = any(re.search(rf"\b{re.escape(t)}\b", full_text) for t in case_terms)
            has_ctrl = any(re.search(rf"\b{re.escape(t)}\b", full_text) for t in ctrl_terms)

        if has_case and not has_ctrl:
            cond = case_label
        elif has_ctrl and not has_case:
            cond = ctrl_label
        else:
            cond = "unknown"
        rows.append({"sample_id": gsm_name, "condition": cond})

    meta_df = pd.DataFrame(rows).set_index("sample_id")
    meta_df.to_csv(metadata_path)

    n_case = int((meta_df["condition"] == case_label).sum())
    n_ctrl = int((meta_df["condition"] == ctrl_label).sum())
    log.info(f"  {acc}: saved {expr.shape[0]} genes x {expr.shape[1]} samples "
             f"| case={n_case} ctrl={n_ctrl}")

    return record(acc, "ok", str(parquet_path),
                  f"{expr.shape[0]} genes, case={n_case}, ctrl={n_ctrl}")


# ── PRIDE download ────────────────────────────────────────────────────────────

def pride_is_cached(acc: str) -> bool:
    d = PRIDE_RAW / acc
    if not d.exists():
        return False
    files = list(d.glob("*.txt")) + list(d.glob("*.tsv")) + list(d.glob("*.csv"))
    return len(files) > 0


def _pride_file_score(name: str) -> int:
    n = name.lower()
    score = 0
    if any(k in n for k in ["proteingroup", "protein_group", "protein-group",
                             "proteingroupsreport", "proteinreport"]):
        score += 10
    if any(k in n for k in ["lfq", "ibaq", "abundance", "quant", "msstats"]):
        score += 5
    if "sdrf" in n:
        score += 4
    if any(n.endswith(e) for e in [".txt", ".tsv", ".csv", ".zip"]):
        score += 2
    if any(k in n for k in ["raw", ".mzml", ".wiff", ".mzxml", "peak",
                              "fraction", "library"]):
        score -= 20
    # Penalise very large archives
    return score


def _pride_pick_url(f: dict) -> str:
    """Return a downloadable HTTPS URL from a PRIDE file record."""
    link = f.get("downloadLink", "")
    if link and link.startswith("http"):
        return link
    for loc in f.get("publicFileLocations", []):
        val = loc.get("value", "")
        if val.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            # EBI PRIDE FTP → HTTPS mirror at ftp.ebi.ac.uk
            val = "https://ftp.ebi.ac.uk/" + val[len("ftp://ftp.pride.ebi.ac.uk/"):]
            return val
        if val.startswith("ftp://"):
            val = "https://" + val[len("ftp://"):]
            return val
        if val.startswith("http"):
            return val
    return ""


def process_pride(acc: str) -> dict:
    out_dir = PRIDE_RAW / acc
    out_dir.mkdir(parents=True, exist_ok=True)

    if pride_is_cached(acc):
        existing = list(out_dir.glob("*.txt")) + list(out_dir.glob("*.tsv"))
        log.info(f"  {acc}: already cached ({len(existing)} files) — skipping")
        return record(acc, "cached", str(out_dir))

    log.info(f"  {acc}: fetching PRIDE file list …")
    files = []
    try:
        # v2 projects endpoint (returns a JSON array directly)
        url = f"{PRIDE_API}/projects/{acc}/files?pageSize=500"
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        body = r.text.strip()
        if body:
            data = r.json()
            files = data if isinstance(data, list) else data.get("_embedded", {}).get("files", [])
    except Exception as e:
        log.error(f"  {acc}: PRIDE API failed — {e}")
        return record(acc, "failed", notes=str(e))

    if not files:
        # Fallback: old endpoint
        try:
            url2 = f"{PRIDE_API}/files/byProject?accession={acc}&pageSize=200"
            r2 = requests.get(url2, timeout=30)
            r2.raise_for_status()
            files = r2.json().get("_embedded", {}).get("files", [])
        except Exception as e2:
            log.warning(f"  {acc}: fallback API also failed — {e2}")

    if not files:
        log.error(f"  {acc}: no files found from PRIDE API")
        return record(acc, "failed", notes="no files from API")

    # Filter and rank
    result_files = [f for f in files
                    if _pride_file_score(f.get("fileName", "")) > 0]
    result_files.sort(key=lambda f: _pride_file_score(f.get("fileName", "")),
                      reverse=True)

    log.info(f"  {acc}: {len(result_files)} candidate files (from {len(files)} total)")

    downloaded = []
    for f in result_files[:8]:
        name = f.get("fileName", "unknown")
        url  = _pride_pick_url(f)
        if not url:
            continue
        dest = out_dir / name
        if dest.exists():
            downloaded.append(dest)
            continue
        try:
            log.info(f"    downloading {name} …")
            with requests.get(url, stream=True, timeout=180) as resp:
                resp.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in resp.iter_content(1 << 20):
                        if chunk:
                            fh.write(chunk)
            downloaded.append(dest)
            time.sleep(0.5)
        except Exception as e:
            log.warning(f"    {name}: failed — {e}")

    if not downloaded:
        return record(acc, "failed", notes="no files downloaded")

    log.info(f"  {acc}: {len(downloaded)} files saved to {out_dir}")
    return record(acc, "ok", str(out_dir),
                  f"{len(downloaded)} files: {[f.name for f in downloaded[:3]]}")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def show_status():
    if not MANIFEST.exists():
        log.info("No manifest yet — run a download first.")
        return
    df = pd.read_csv(MANIFEST)
    print("\n" + "=" * 65)
    print("DOWNLOAD MANIFEST")
    print("=" * 65)
    print(df[["accession", "status", "timestamp", "notes"]].to_string(index=False))
    ok    = (df["status"].isin(["ok", "cached"])).sum()
    fail  = (df["status"] == "failed").sum()
    print(f"\nTotal: {len(df)}  OK/cached: {ok}  Failed: {fail}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only",   choices=["geo", "pride"],
                        help="Download only this source")
    parser.add_argument("--accessions", nargs="*",
                        help="Specific accessions to download")
    parser.add_argument("--status", action="store_true",
                        help="Show download manifest and exit")
    parser.add_argument("--force",  action="store_true",
                        help="Re-download even if cached")
    args = parser.parse_args()

    if args.status:
        show_status()
        return

    if not check_t7():
        sys.exit(1)

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    # ── GEO ───────────────────────────────────────────────────────────────
    if args.only in (None, "geo"):
        geo_cohorts = load_cohorts(["GEO"])
        if args.accessions:
            geo_cohorts = geo_cohorts[geo_cohorts["accession"].isin(args.accessions)]

        log.info(f"\n{'='*60}\nGEO — {len(geo_cohorts)} cohorts\n{'='*60}")
        for _, row in geo_cohorts.iterrows():
            acc = row["accession"]
            if not args.force and geo_is_cached(acc):
                log.info(f"  {acc}: cached — skipping")
                manifest_rows.append(record(acc, "cached",
                                            str(GEO_RAW / acc)))
                continue
            result = process_geo(acc, row["case_label"], row["control_label"])
            manifest_rows.append(result)

    # ── PRIDE ─────────────────────────────────────────────────────────────
    if args.only in (None, "pride"):
        pride_cohorts = load_cohorts(["PRIDE"])
        if args.accessions:
            pride_cohorts = pride_cohorts[pride_cohorts["accession"].isin(args.accessions)]

        log.info(f"\n{'='*60}\nPRIDE — {len(pride_cohorts)} cohorts\n{'='*60}")
        for _, row in pride_cohorts.iterrows():
            acc = row["accession"]
            if not args.force and pride_is_cached(acc):
                log.info(f"  {acc}: cached — skipping")
                manifest_rows.append(record(acc, "cached",
                                            str(PRIDE_RAW / acc)))
                continue
            result = process_pride(acc)
            manifest_rows.append(result)

    # ── Save manifest ──────────────────────────────────────────────────────
    if manifest_rows:
        manifest_df = pd.DataFrame(manifest_rows)
        # Merge with existing manifest if present
        if MANIFEST.exists():
            existing = pd.read_csv(MANIFEST)
            updated = existing[~existing["accession"].isin(manifest_df["accession"])]
            manifest_df = pd.concat([updated, manifest_df], ignore_index=True)
        manifest_df.to_csv(MANIFEST, index=False)

    show_status()
    log.info(f"\nFull log: /tmp/op_download.log")


if __name__ == "__main__":
    main()
