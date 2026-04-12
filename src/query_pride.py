"""
query_pride.py
--------------
For each PRIDE cohort in cohorts.csv:
  - Download protein group table via PRIDE REST API
  - Check whether each panel gene/protein is present
  - Compute mean intensity in cases vs controls
  - Record direction

Output: results/pride_query.csv
"""

import argparse
import datetime as dt
import json
import logging
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

from disease_labels import pride_classify_text as _pride_classify_text, PRIDE_CASE_KEYWORDS, PRIDE_CTRL_KEYWORDS

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
PREFERRED_EXTERNAL_ROOT = Path("/Volumes/T7/5-Alzhimers_Parkisons_MS_External_Valida/op_external_validation_data")


def choose_external_root() -> Path:
    env_root = os.getenv("OP_EXTERNAL_DATA_ROOT")
    if env_root:
        return Path(env_root).expanduser()
    if PREFERRED_EXTERNAL_ROOT.exists():
        return PREFERRED_EXTERNAL_ROOT
    return ROOT / "data" / "external"


EXTERNAL_DATA_ROOT = choose_external_root()
RAW  = EXTERNAL_DATA_ROOT / "pride_raw"
OUT  = ROOT / "results"
CFG  = ROOT / "config"
PRIDE_OVERRIDES = CFG / "pride_overrides.csv"

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]

OUTPUT_COLUMNS = [
    "accession", "disease", "biospecimen", "modality", "exposure_type",
    "gene", "present", "n_case", "n_control", "mean_case", "mean_ctrl",
    "direction", "expected", "concordant",
]

GENE_COL_CANDIDATES = [
    "Gene Names", "Gene names", "gene_names", "Gene Name",
    "Gene Symbol", "gene_symbol", "Genes",
    "PG.Genes", "Protein.names", "gene",
]
INTENSITY_PREFIXES = [
    "LFQ intensity ", "Intensity ", "iBAQ ", "Abundance: ",
    "Reporter intensity ", "MS2 ",
]

REQ_HEADERS = {
    "Accept": "application/json",
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
}


# ── helpers ──────────────────────────────────────────────────────────────────

def load_expected(targets_yaml: Path) -> dict:
    with open(targets_yaml) as f:
        cfg = yaml.safe_load(f)
    return {g: v for g, v in cfg["targets"].items()}


def load_overrides(path: Path) -> dict[str, list[dict]]:
    """
    Optional deterministic override list for PRIDE files.
    Expected columns: accession,file_name,download_url,file_role,notes
    """
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path)
    except Exception:
        return {}

    required = {"accession", "download_url"}
    if not required.issubset(df.columns):
        return {}

    out: dict[str, list[dict]] = {}
    for _, row in df.iterrows():
        acc = str(row.get("accession", "")).strip()
        url = str(row.get("download_url", "")).strip()
        if not acc or not url:
            continue
        name = str(row.get("file_name", "")).strip() or Path(url).name
        out.setdefault(acc, []).append({"fileName": name, "downloadLink": url})
    return out


def absent_records_for_cohort(row: pd.Series, targets: dict) -> list[dict]:
    records = []
    for gene in PANEL:
        exp = expected_direction(gene, row["exposure_type"], targets)
        records.append(dict(
            accession=row["accession"],
            disease=row["disease"],
            biospecimen=row["biospecimen"],
            modality="proteomics",
            exposure_type=row["exposure_type"],
            gene=gene,
            present=False,
            n_case=0,
            n_control=0,
            mean_case=None,
            mean_ctrl=None,
            direction="absent",
            expected=exp,
            concordant="na",
        ))
    return records


def expected_direction(gene: str, exposure_type: str, targets: dict) -> str:
    t = targets.get(gene, {})
    if exposure_type == "acute":
        return t.get("acute_direction") or "none"
    if exposure_type in ("chronic", "chronic_specific", "prenatal_chronic"):
        return t.get("chronic_direction") or "none"
    chronic = t.get("chronic_direction") or "none"
    return chronic if chronic != "none" else "none"


def direction(mean_case: float, mean_ctrl: float, threshold: float = 0.1) -> str:
    diff = mean_case - mean_ctrl
    if diff < -threshold:
        return "down"
    if diff > threshold:
        return "up"
    return "flat"


def _request_json(url: str, params: dict | None = None) -> dict | list:
    r = requests.get(url, params=params, headers=REQ_HEADERS, timeout=45)
    r.raise_for_status()
    text = r.text.strip()
    if not text:
        raise ValueError(f"Empty response from {url}")
    try:
        return r.json()
    except json.JSONDecodeError:
        preview = text[:180].replace("\n", " ")
        raise ValueError(f"Non-JSON response from {url}: {preview}")


def _extract_files(payload) -> list[dict]:
    if isinstance(payload, list):
        return payload
    if not isinstance(payload, dict):
        return []

    # Common REST envelope patterns.
    if "_embedded" in payload and isinstance(payload["_embedded"], dict):
        emb = payload["_embedded"]
        for key in ("files", "compactFiles", "compactfiles", "fileList", "items", "content"):
            v = emb.get(key)
            if isinstance(v, list):
                return v
        for v in emb.values():
            if isinstance(v, list):
                return v

    for key in ("files", "content", "items", "list"):
        v = payload.get(key)
        if isinstance(v, list):
            return v

    return []


def _normalize_file_items(items: list[dict]) -> list[dict]:
    out = []
    for it in items:
        if not isinstance(it, dict):
            continue
        name = (it.get("fileName") or it.get("filename") or it.get("name") or "").strip()
        link = (
            it.get("downloadLink")
            or it.get("downloadlink")
            or it.get("ftpLink")
            or it.get("httpLink")
            or next((loc.get("value", "") for loc in it.get("publicFileLocations", []) if isinstance(loc, dict)), "")
            or ""
        ).strip()
        if name or link:
            out.append({
                **it,
                "fileName": name or Path(link).name,
                "downloadLink": link,
            })
    return out


def _files_from_pride_api(accession: str) -> list[dict]:
    candidates = [
        (f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files", {"page": 0, "pageSize": 500}),
        (f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files/all", None),
        (f"{PRIDE_API}/files/byProject", {"accession": accession, "pageSize": 500}),
    ]
    for url, params in candidates:
        try:
            payload = _request_json(url, params=params)
            items = _normalize_file_items(_extract_files(payload))
            if items:
                return items
        except Exception:
            continue
    return []


def _files_from_ftp_listing(accession: str) -> list[dict]:
    """
    Fallback: discover project directory from ProteomeCentral page and parse
    the PRIDE FTP HTTP listing.
    """
    px_url = f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={accession}"
    try:
        r = requests.get(px_url, headers=REQ_HEADERS, timeout=45)
        r.raise_for_status()
        html = r.text
    except Exception:
        return []

    m = re.search(
        rf"ftp://ftp\.pride\.ebi\.ac\.uk/pride/data/archive/\d{{4}}/\d{{2}}/{re.escape(accession)}",
        html,
        flags=re.IGNORECASE,
    )
    ftp_http = None
    if m:
        ftp_http = m.group(0).replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
        if not ftp_http.endswith("/"):
            ftp_http += "/"

    if not ftp_http:
        ftp_http = _discover_archive_dir_bruteforce(accession)
        if ftp_http:
            log.info(f"  {accession}: discovered FTP archive path by scan: {ftp_http}")
        else:
            return []

    listing = _download_listing(ftp_http)
    if not listing:
        # One more try via brute-force discovery (covers stale ProteomeCentral links).
        scanned = _discover_archive_dir_bruteforce(accession)
        if scanned and scanned != ftp_http:
            ftp_http = scanned
            listing = _download_listing(ftp_http)
    if not listing:
        return []

    items = []
    for href, name in re.findall(r'href="([^"]+)">([^<]+)</a>', listing, flags=re.IGNORECASE):
        nm = name.strip()
        if not nm or nm in ("Parent Directory", "../"):
            continue
        if nm.endswith("/"):
            continue
        link = href if href.startswith("http") else ftp_http + href.lstrip("/")
        items.append({"fileName": nm, "downloadLink": link})

    return items


def _download_listing(ftp_http: str) -> str | None:
    try:
        r = requests.get(ftp_http, headers=REQ_HEADERS, timeout=45)
        r.raise_for_status()
        return r.text
    except Exception:
        return None


def _discover_archive_dir_bruteforce(accession: str) -> str | None:
    """
    Last-resort resolver for PRIDE archive directories:
    https://ftp.pride.ebi.ac.uk/pride/data/archive/YYYY/MM/PXDxxxxxx/
    """
    current_year = dt.date.today().year
    for year in range(current_year, 2011, -1):
        for month in range(12, 0, -1):
            url = f"https://ftp.pride.ebi.ac.uk/pride/data/archive/{year:04d}/{month:02d}/{accession}/"
            try:
                r = requests.get(url, headers=REQ_HEADERS, timeout=15)
                if r.status_code == 200 and ("href=" in r.text.lower() or "<html" in r.text.lower()):
                    return url
            except Exception:
                continue
    return None


def get_files(accession: str) -> list[dict]:
    items = _files_from_pride_api(accession)
    if items:
        return items

    items = _files_from_ftp_listing(accession)
    if items:
        log.info(f"  {accession}: using FTP listing fallback ({len(items)} files)")
        return items

    raise RuntimeError(f"No file list could be resolved for {accession}")


def score_file(name: str) -> int:
    n = name.lower()
    score = 0
    if any(k in n for k in ["proteingroup", "protein_group", "protein-group",
                             "proteingroupsreport", "proteinreport"]):
        score += 10
    if any(k in n for k in ["lfq", "ibaq", "abundance", "quant", "msstats"]):
        score += 5
    if "sdrf" in n:
        score += 3
    if any(n.endswith(e) for e in [".txt", ".tsv", ".csv"]):
        score += 2
    if n.endswith(".zip") and any(k in n for k in ["proteingroup", "proteingroupsreport",
                                                    "result", "protein"]):
        score += 8  # Spectronaut/DIA-NN report zips are high priority
    if any(k in n for k in [".raw", ".mzml", ".wiff", ".mzxml", "peak"]):
        score -= 20
    if any(k in n for k in ["fraction", "library", "dda", "acquisition"]):
        score -= 5
    return score


def pick_url(f: dict) -> str:
    return (f.get("downloadLink")
            or next((loc.get("value", "") for loc in f.get("publicFileLocations", [])), "")
            or "")


def download_file(url: str, dest: Path) -> Path | None:
    if dest.exists():
        return dest
    try:
        with requests.get(url, stream=True, timeout=120, headers=REQ_HEADERS) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    if chunk:
                        fh.write(chunk)
        return dest
    except Exception as e:
        log.warning(f"  download failed {dest.name}: {e}")
        return None


def load_protein_table(path: Path) -> pd.DataFrame | None:
    try:
        if path.suffix == ".zip":
            import zipfile, io
            with zipfile.ZipFile(path) as zf:
                # Find the best CSV/TSV inside the zip
                members = sorted(
                    zf.namelist(),
                    key=lambda n: (
                        0 if any(k in n.lower() for k in ["proteingroup", "protein group", "protein_group"]) else 1,
                        -zf.getinfo(n).file_size,
                    )
                )
                SKIP_SUFFIXES = {"log", "analysislog", "parameter", "params", "summary", "readme"}
                csv_members = [n for n in members
                               if any(n.lower().endswith(e) for e in [".csv", ".tsv", ".txt"])
                               and not any(n.lower().endswith(f"_{s}.txt") or n.lower().endswith(f"{s}.txt")
                                           for s in SKIP_SUFFIXES)
                               and "analysislog" not in n.lower()
                               and "log.txt" not in n.lower()
                               and "summary" not in n.lower()
                               and "parameter" not in n.lower()]
                if not csv_members:
                    log.debug(f"  no CSV/TSV inside {path.name}")
                    return None
                with zf.open(csv_members[0]) as f:
                    buf = io.BytesIO(f.read())
                nm = csv_members[0].lower()
                sep = "\t" if nm.endswith(".tsv") or nm.endswith(".txt") else ","
                df = pd.read_csv(buf, sep=sep, low_memory=False)
                log.debug(f"  loaded {csv_members[0]} from {path.name}: {df.shape}")
                return df
        sep = "\t" if path.suffix in (".txt", ".tsv") else None
        df = pd.read_csv(path, sep=sep, engine="python", low_memory=False)
        return df
    except Exception as e:
        log.debug(f"  cannot parse {path.name}: {e}")
        return None


def find_gene_col(df: pd.DataFrame) -> str | None:
    for c in GENE_COL_CANDIDATES:
        if c in df.columns:
            return c
    for c in df.columns:
        if "gene" in c.lower():
            return c
    return None


def intensity_cols(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns
            if any(c.startswith(p) for p in INTENSITY_PREFIXES)]
    # Spectronaut format: columns end in .PG.Quantity
    if not cols:
        cols = [c for c in df.columns if c.endswith(".PG.Quantity")]
    if not cols:
        # fallback: all numeric columns that aren't gene/protein identifiers
        skip = {"gene", "protein", "majority", "id", "sequence", "mass",
                "score", "count", "peptide", "unique", "razor", "pg.genes",
                "pg.proteinaccessions", "pg.proteingroups"}
        cols = [c for c in df.select_dtypes(include=[np.number]).columns
                if not any(k in c.lower() for k in skip)]
    return cols


def parse_first_gene(val) -> str | None:
    if pd.isna(val):
        return None
    s = re.split(r"[;,|\s]+", str(val).strip())[0].strip()
    return s.upper() if s else None


def build_gene_intensity(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    """Gene x sample log2-intensity matrix."""
    int_cols = intensity_cols(df)
    if not int_cols:
        return pd.DataFrame()

    sub = df[[gene_col] + int_cols].copy()
    sub[int_cols] = sub[int_cols].apply(pd.to_numeric, errors="coerce")
    sub["gene"] = sub[gene_col].map(parse_first_gene)
    sub = sub.dropna(subset=["gene"])
    sub = sub.groupby("gene")[int_cols].mean()

    # Replace 0 → NaN (not detected), then log2
    sub = sub.replace(0, np.nan)
    numeric = sub.select_dtypes(include=[np.number])
    log2_sub = np.log2(numeric + 1)
    return log2_sub


def infer_sample_labels(cols: list[str], case_label: str, ctrl_label: str) -> dict[str, str]:
    """Guess case/control from intensity column names."""
    labels = {}
    cl = case_label.lower().replace("_", " ")
    ctl = ctrl_label.lower().replace("_", " ")
    for c in cols:
        cn = c.lower()
        if cl in cn or case_label.lower() in cn:
            labels[c] = "case"
        elif ctl in cn or ctrl_label.lower() in cn:
            labels[c] = "control"
        else:
            labels[c] = "unknown"
    return labels


def load_sdrf_labels(acc_dir: Path, int_cols: list[str],
                     case_label: str, ctrl_label: str) -> dict[str, str] | None:
    """Try to read sample labels from an SDRF file."""
    for sdrf in acc_dir.glob("*.sdrf.tsv"):
        try:
            sdrf_df = pd.read_csv(sdrf, sep="\t", low_memory=False)
            sample_col = next((c for c in sdrf_df.columns
                               if "source name" in c.lower()), None)
            char_cols  = [c for c in sdrf_df.columns
                          if "characteristics" in c.lower() or "factor value" in c.lower()]
            if not sample_col or not char_cols:
                continue

            labels = {}
            for _, row in sdrf_df.iterrows():
                sample = str(row[sample_col])
                text = " ".join(str(row[c]) for c in char_cols
                                if pd.notna(row[c])).lower()
                if case_label.lower() in text:
                    labels[sample] = "case"
                elif ctrl_label.lower() in text:
                    labels[sample] = "control"
                else:
                    labels[sample] = "unknown"

            # Map SDRF sample names to intensity column names (partial match)
            col_labels = {}
            for col in int_cols:
                for sample, lbl in labels.items():
                    if sample in col or col in sample:
                        col_labels[col] = lbl
                        break
            if col_labels:
                return col_labels
        except Exception:
            pass

    # Fallback: look for Excel annotation files (e.g. Spectronaut annotation xlsx inside zips)
    # _classify_text -> _pride_classify_text from disease_labels (shared module)
    for xlsx in list(acc_dir.glob("*.xlsx")) + list(acc_dir.glob("*.xls")):
        try:
            ann = pd.read_excel(xlsx, engine="openpyxl")
            sample_col = next((c for c in ann.columns if "sample" in c.lower() and "name" in c.lower()), None)
            if sample_col is None:
                continue
            label_col = next((c for c in ann.columns
                              if any(k in c.lower() for k in ["classification", "diagnosis", "disease",
                                                               "status", "condition", "group"])), None)
            if label_col is None:
                continue
            col_labels = {}
            for col in int_cols:
                for _, row in ann.iterrows():
                    sname = str(row[sample_col])
                    if sname in col or col in sname:
                        lbl = _pride_classify_text(str(row[label_col]), case_label, ctrl_label)
                        col_labels[col] = lbl
                        break
            if col_labels:
                return col_labels
        except Exception:
            pass

    # Also check xlsx files inside any zip in acc_dir
    for zpath in acc_dir.glob("*.zip"):
        try:
            import zipfile, io
            with zipfile.ZipFile(zpath) as zf:
                xlsx_members = [n for n in zf.namelist() if n.lower().endswith((".xlsx", ".xls"))]
                for nm in xlsx_members:
                    with zf.open(nm) as f:
                        buf = io.BytesIO(f.read())
                    ann = pd.read_excel(buf, engine="openpyxl")
                    sample_col = next((c for c in ann.columns if "sample" in c.lower() and "name" in c.lower()), None)
                    if sample_col is None:
                        continue
                    label_col = next((c for c in ann.columns
                                      if any(k in c.lower() for k in ["classification", "diagnosis", "disease",
                                                                       "status", "condition", "group"])), None)
                    if label_col is None:
                        continue
                    col_labels = {}
                    for col in int_cols:
                        for _, row in ann.iterrows():
                            sname = str(row[sample_col])
                            if sname in col or col in sname:
                                lbl = _pride_classify_text(str(row[label_col]), case_label, ctrl_label)
                                col_labels[col] = lbl
                                break
                    if col_labels:
                        log.info(f"  Labels from {zpath.name}/{nm}: {sum(1 for l in col_labels.values() if l=='case')} case, {sum(1 for l in col_labels.values() if l=='control')} control")
                        return col_labels
        except Exception:
            pass

    return None


# ── main query ───────────────────────────────────────────────────────────────

def query_cohort(row: pd.Series, targets: dict, raw_dir: Path,
                 overrides: dict[str, list[dict]] | None = None) -> list[dict]:
    acc = row["accession"]
    acc_dir = raw_dir / acc
    acc_dir.mkdir(parents=True, exist_ok=True)

    if overrides and acc in overrides:
        files = overrides[acc]
        log.info(f"  {acc}: using {len(files)} override file(s)")
    else:
        try:
            files = get_files(acc)
        except Exception as e:
            log.error(f"  {acc}: PRIDE API failed — {e}")
            cached = [p for p in acc_dir.glob("*") if p.is_file()]
            if not cached:
                return absent_records_for_cohort(row, targets)
            log.warning(f"  {acc}: using {len(cached)} cached local files")
            files = [{"fileName": p.name, "downloadLink": ""} for p in cached]

    # Sort by score, download top candidates
    files_sorted = sorted(files, key=lambda f: score_file(f.get("fileName", "")), reverse=True)
    downloaded = []
    for f in files_sorted[:10]:
        name = f.get("fileName", "unknown")
        url  = pick_url(f)
        dest = acc_dir / name
        if not url:
            # API returned no URL — use cached file if present
            if dest.exists():
                downloaded.append(dest)
            continue
        p = download_file(url, dest)
        if p:
            downloaded.append(p)
        elif dest.exists():
            # Download failed (e.g. FTP/ASCP URL) but file is cached — use it
            downloaded.append(dest)
        time.sleep(0.3)

    # Find first parseable protein group file
    gene_matrix = None
    for path in sorted(downloaded, key=lambda p: score_file(p.name), reverse=True):
        df = load_protein_table(path)
        if df is None or df.empty:
            continue
        gene_col = find_gene_col(df)
        if gene_col is None:
            continue
        gene_matrix = build_gene_intensity(df, gene_col)
        if not gene_matrix.empty:
            log.info(f"  {acc}: using {path.name} — {gene_matrix.shape[0]} proteins")
            break

    if gene_matrix is None or gene_matrix.empty:
        log.warning(f"  {acc}: no usable protein table found")
        gene_matrix = pd.DataFrame()

    # Sample labels
    int_cols_list = gene_matrix.columns.tolist() if not gene_matrix.empty else []
    labels = load_sdrf_labels(acc_dir, int_cols_list,
                              row["case_label"], row["control_label"])
    if labels is None:
        labels = infer_sample_labels(int_cols_list, row["case_label"], row["control_label"])

    cases = [c for c, l in labels.items() if l == "case"]
    ctrls = [c for c, l in labels.items() if l == "control"]
    log.info(f"  {acc}: {len(cases)} case cols, {len(ctrls)} control cols")

    records = []
    for gene in PANEL:
        present = (not gene_matrix.empty) and (gene in gene_matrix.index)
        if present and cases and ctrls:
            case_cols = [c for c in cases if c in gene_matrix.columns]
            ctrl_cols = [c for c in ctrls if c in gene_matrix.columns]
            mc = float(gene_matrix.loc[gene, case_cols].mean()) if case_cols else np.nan
            mk = float(gene_matrix.loc[gene, ctrl_cols].mean()) if ctrl_cols else np.nan
            obs_dir = direction(mc, mk) if not (np.isnan(mc) or np.isnan(mk)) else "no_labels"
        else:
            mc = mk = np.nan
            obs_dir = "absent" if not present else "no_labels"

        exp = expected_direction(gene, row["exposure_type"], targets)
        concordant = ("yes" if obs_dir == exp
                      else "no" if obs_dir in ("up", "down") and exp != "none"
                      else "na")

        records.append(dict(
            accession=acc,
            disease=row["disease"],
            biospecimen=row["biospecimen"],
            modality="proteomics",
            exposure_type=row["exposure_type"],
            gene=gene,
            present=present,
            n_case=len(cases),
            n_control=len(ctrls),
            mean_case=round(mc, 4) if not np.isnan(mc) else None,
            mean_ctrl=round(mk, 4) if not np.isnan(mk) else None,
            direction=obs_dir,
            expected=exp,
            concordant=concordant,
        ))

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohorts",    default=str(CFG / "cohorts.csv"))
    parser.add_argument("--targets",    default=str(CFG / "targets.yaml"))
    parser.add_argument("--out",        default=str(OUT / "pride_query.csv"))
    parser.add_argument("--raw-dir",    default=str(RAW))
    parser.add_argument("--accessions", nargs="*")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    targets = load_expected(Path(args.targets))
    overrides = load_overrides(PRIDE_OVERRIDES)

    cohorts = pd.read_csv(args.cohorts)
    cohorts = cohorts[(cohorts["repository"] == "PRIDE") & (cohorts["access_type"] == "open")]
    if args.accessions:
        cohorts = cohorts[cohorts["accession"].isin(args.accessions)]

    all_records = []
    for _, row in cohorts.iterrows():
        log.info(f"Querying {row['accession']} — {row['disease']}")
        records = query_cohort(row, targets, Path(args.raw_dir), overrides=overrides)
        all_records.extend(records)

    df = pd.DataFrame(all_records, columns=OUTPUT_COLUMNS)
    df.to_csv(args.out, index=False)
    log.info(f"Saved {len(df)} rows to {args.out}")


if __name__ == "__main__":
    main()
