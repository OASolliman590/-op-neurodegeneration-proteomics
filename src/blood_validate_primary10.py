#!/usr/bin/env python3
"""
blood_validate_primary10.py
---------------------------
Primary blood proteomics validation runner for the frozen 10-target OP panel.

Implements dataset-level stopping rules:
1) blocked_no_quant_table
2) blocked_no_id_mapping
3) blocked_no_group_labels
4) low_coverage (panel coverage < 2 genes)

Outputs (results/blood_validation):
- per_dataset_effects.csv
- per_dataset_qc.csv
- meta_plasma.csv
- meta_serum.csv
- summary_primary10.md
- selected_file_inventory.csv (intermediate)
"""

from __future__ import annotations

import argparse
import csv
import ftplib
import gzip
import io
import math
import os
import re
import time
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests
import yaml

from disease_labels import direction_from_fc, concordance, pride_classify_text

ROOT = Path(__file__).resolve().parents[1]
PREFERRED_EXTERNAL_ROOT = Path("/Volumes/T7/5-Alzhimers_Parkisons_MS_External_Valida/op_external_validation_data")

PANEL = ["ACTG1", "DNAH9", "GPX3", "VWF", "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B"]

# Candidate quant-like files
ACCEPT_EXT = (
    ".txt", ".tsv", ".csv", ".txt.gz", ".tsv.gz", ".csv.gz", ".zip",
    ".mztab", ".mztab.gz", ".xlsx", ".xls",
)
BLOCK_EXT = (".raw", ".wiff", ".d", ".sne", ".mzml", ".mgf")
NAME_KEYWORDS = (
    "proteingroup", "protein_group", "protein groups", "proteins", "proteinreport",
    "quant", "abundance", "lfq", "intensity", "result", "report", "mztab",
)
ANTI_QUANT_KEYWORDS = (
    "rawfile", "rawfiles", "spectronautsession", "session", "library", "fractionraw",
    "singleshot", "wiff", "mzml", "ms2", "dda", "dia",
)
NON_QUANT_HINTS = (
    "readme", "license", "fasta", "fa.", ".fa", "metadata", "parameter", "params",
    "peptide", "peptides", "psm", "spectra", "scan", "chromatogram", "feature",
    "evidence", "identification", "mzid", "msms", "qc_", "qualitycontrol", "manifest",
    "sample_annotation", "design", "protocol", "checksum", "md5", "sha1", "sha256",
)
FALLBACK_MAX_FILES = 5
FALLBACK_MAX_SIZE = 700 * 1024 * 1024  # allow moderately large processed archives
ZIP_MEMBER_MAX_COUNT = 250
ZIP_MEMBER_MAX_BYTES = 150 * 1024 * 1024

CASE_KEYS = ("case", "disease", "ms", "pd", "parkinson", "alzheimer", "patient", "affected")
CTRL_KEYS = ("control", "ctrl", "healthy", "normal", "hc")

STATUS_PRIORITY = {
    "ok": 5,
    "low_coverage": 4,
    "blocked_no_group_labels": 3,
    "blocked_no_id_mapping": 2,
    "blocked_no_quant_table": 1,
    "parse_failed": 0,
}


def choose_external_root() -> Path:
    env_root = os.getenv("OP_EXTERNAL_DATA_ROOT")
    if env_root:
        return Path(env_root).expanduser()
    if PREFERRED_EXTERNAL_ROOT.exists():
        return PREFERRED_EXTERNAL_ROOT
    return ROOT / "data" / "external"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default=str(ROOT / "config" / "blood_primary10.csv"))
    p.add_argument("--targets", default=str(ROOT / "config" / "targets.yaml"))
    p.add_argument("--out-dir", default=str(ROOT / "results" / "blood_validation"))
    p.add_argument("--raw-dir", default=None, help="Cache directory for downloaded PRIDE candidate files")
    p.add_argument("--sample-map-dir", default=str(ROOT / "config" / "pride_sample_maps"), help="Optional accession-level sample label maps: <ACCESSION>.csv with columns sample,label")
    p.add_argument("--max-file-bytes", type=int, default=700 * 1024 * 1024)
    p.add_argument("--max-files-per-accession", type=int, default=20)
    p.add_argument("--threshold", type=float, default=0.3)
    p.add_argument("--sleep-seconds", type=float, default=0.05)
    p.add_argument("--accessions", nargs="*", default=None, help="Optional subset accessions")
    return p.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"accession", "disease", "biospecimen", "source", "case_label", "control_label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    df["accession"] = df["accession"].astype(str).str.strip().str.upper()
    return df


def read_expected_chronic(path: Path) -> Dict[str, str]:
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    expected: Dict[str, str] = {}
    for gene in PANEL:
        g = cfg.get("targets", {}).get(gene, {})
        direction = str(g.get("chronic_direction", "none")).strip().lower()
        if direction not in {"up", "down", "none"}:
            direction = "none"
        expected[gene] = direction
    return expected


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _ftp_dir_entries(ftp: ftplib.FTP, path: str) -> List[Tuple[str, bool, int]]:
    lines: List[str] = []
    ftp.dir(path, lines.append)
    out: List[Tuple[str, bool, int]] = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < 9:
            continue
        name = " ".join(parts[8:]).strip()
        if not name or name in (".", ".."):
            continue
        is_dir = parts[0].startswith("d")
        try:
            size = int(parts[4])
        except Exception:
            size = -1
        out.append((name, is_dir, size))
    return out


def _ftp_collect_files(
    ftp: ftplib.FTP,
    root: str,
    max_depth: int = 5,
    max_files: int = 8000,
) -> List[Tuple[str, int, str]]:
    files: List[Tuple[str, int, str]] = []
    stack: List[Tuple[str, int]] = [(root, 0)]
    visited = set()
    while stack and len(files) < max_files:
        curr, depth = stack.pop()
        if curr in visited:
            continue
        visited.add(curr)
        try:
            entries = _ftp_dir_entries(ftp, curr)
        except Exception:
            continue
        for name, is_dir, size in entries:
            full = f"{curr}/{name}"
            if is_dir:
                if depth < max_depth:
                    stack.append((full, depth + 1))
                continue
            rel = full[len(root) + 1:] if full.startswith(root + "/") else name
            url = f"https://ftp.pride.ebi.ac.uk{full}"
            files.append((rel, size, url))
            if len(files) >= max_files:
                break
    return files


def list_pride_files(accession: str, timeout: int = 30) -> Tuple[Optional[str], List[Tuple[str, int, str]]]:
    try:
        ftp = ftplib.FTP("ftp.pride.ebi.ac.uk", timeout=timeout)
        ftp.login()
        base = "/pride/data/archive"
        years = sorted(ftp.nlst(base), reverse=True)
        first_seen_path: Optional[str] = None
        for y in years:
            try:
                months = sorted(ftp.nlst(y), reverse=True)
            except Exception:
                continue
            for m in months:
                path = f"{m}/{accession}"
                lines: List[str] = []
                try:
                    ftp.dir(path, lines.append)
                except Exception:
                    continue
                if first_seen_path is None:
                    first_seen_path = path
                # Shallow scan (top-level + one nested directory level) for stability.
                out: List[Tuple[str, int, str]] = []
                try:
                    entries = _ftp_dir_entries(ftp, path)
                except Exception:
                    entries = []
                for name, is_dir, size in entries:
                    full = f"{path}/{name}"
                    if is_dir:
                        try:
                            sub_entries = _ftp_dir_entries(ftp, full)
                        except Exception:
                            sub_entries = []
                        for sname, sis_dir, ssize in sub_entries:
                            if sis_dir:
                                continue
                            rel = f"{name}/{sname}"
                            url = f"https://ftp.pride.ebi.ac.uk{full}/{sname}"
                            out.append((rel, ssize, url))
                    else:
                        rel = name
                        url = f"https://ftp.pride.ebi.ac.uk{full}"
                        out.append((rel, size, url))
                if out:
                    ftp.quit()
                    return path, out
        ftp.quit()
        if first_seen_path is not None:
            return first_seen_path, []
    except Exception:
        return None, []
    return None, []


def list_pride_files_deep(root_path: str, timeout: int = 30, max_depth: int = 4, max_files: int = 5000) -> List[Tuple[str, int, str]]:
    try:
        ftp = ftplib.FTP("ftp.pride.ebi.ac.uk", timeout=timeout)
        ftp.login()
        files = _ftp_collect_files(ftp, root_path, max_depth=max_depth, max_files=max_files)
        ftp.quit()
        return files
    except Exception:
        return []


def list_iprox_files(accession: str, timeout: int = 30) -> List[Tuple[str, int, str]]:
    url = f"https://www.iprox.cn/proxi/datasets/{accession}"
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        payload = r.json()
    except Exception:
        return []

    out: List[Tuple[str, int, str]] = []
    for item in payload.get("dataFiles", []) or []:
        if not isinstance(item, dict):
            continue
        v = str(item.get("value", "")).strip()
        if not v or not v.startswith(("http://", "https://")):
            continue
        name = v.rsplit("/", 1)[-1] or f"{accession}_file"
        out.append((name, -1, v))
    return out


def list_jpost_files(accession: str, timeout: int = 30) -> List[Tuple[str, int, str]]:
    # jPOST ProXI endpoints can be intermittently unavailable; best-effort only.
    probes = [
        f"https://repository.jpostdb.org/proxi/datasets/{accession}",
        f"https://repository.jpostdb.org/proxi/v0.1/datasets/{accession}",
    ]
    for url in probes:
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code != 200:
                continue
            payload = r.json()
        except Exception:
            continue

        out: List[Tuple[str, int, str]] = []
        for item in payload.get("dataFiles", []) or []:
            if not isinstance(item, dict):
                continue
            v = str(item.get("value", "")).strip()
            if not v or not v.startswith(("http://", "https://")):
                continue
            name = v.rsplit("/", 1)[-1] or f"{accession}_file"
            out.append((name, -1, v))
        if out:
            return out
    return []


def merge_file_lists(primary: List[Tuple[str, int, str]], extra: List[Tuple[str, int, str]]) -> List[Tuple[str, int, str]]:
    merged: Dict[str, Tuple[str, int, str]] = {}
    for name, size, url in primary + extra:
        key = name.lower().strip()
        if key not in merged:
            merged[key] = (name, size, url)
            continue
        old_name, old_size, old_url = merged[key]
        # Prefer entries with known file size and https URLs.
        if (old_size <= 0 and size > 0) or (old_url.startswith("http://") and url.startswith("https://")):
            merged[key] = (name, size, url)
    return list(merged.values())


def pick_quant_candidates(files: List[Tuple[str, int, str]], max_file_bytes: int, max_files_per_accession: int) -> List[Tuple[str, int, str]]:
    quant_files = [(n, s, u) for (n, s, u) in files if is_quant_candidate(n, s, max_file_bytes)]
    dedup: Dict[str, Tuple[str, int, str]] = {}
    for n, s, u in quant_files:
        k = n.lower()
        if k not in dedup:
            dedup[k] = (n, s, u)
    quant_files = list(dedup.values())
    quant_files.sort(key=quant_candidate_sort_key)
    return quant_files[:max_files_per_accession]


def is_quant_candidate(name: str, size: int, max_bytes: int) -> bool:
    n = name.lower()
    if any(n.endswith(ext) for ext in BLOCK_EXT):
        return False
    if any(k in n for k in ANTI_QUANT_KEYWORDS):
        return False
    if any(k in n for k in NON_QUANT_HINTS):
        return False
    if size > 0 and size > max_bytes:
        return False
    # For primary candidates, require both parsable extension and quant-like naming.
    if not any(n.endswith(ext) for ext in ACCEPT_EXT) and not n.endswith((".tab", ".dat")):
        return False
    if any(k in n for k in NAME_KEYWORDS):
        return True
    if n.endswith((".mztab", ".mztab.gz")):
        return True
    return False


def is_fallback_parse_candidate(name: str, size: int) -> bool:
    n = name.lower()
    if any(n.endswith(ext) for ext in BLOCK_EXT):
        return False
    if any(k in n for k in NON_QUANT_HINTS):
        return False
    if size > 0 and size > FALLBACK_MAX_SIZE:
        return False
    if n.endswith(".zip"):
        return True
    if size > 0 and size < 2048:
        return False
    if any(k in n for k in ("matrix", "table", "report", "results", "summary", "processed", "protein", "quant", "intensity", "abundance", "sample")):
        return any(n.endswith(ext) for ext in ACCEPT_EXT) or n.endswith((".tab", ".dat"))
    if n.endswith((".mztab", ".mztab.gz", ".xlsx", ".xls")):
        return True
    return False


def quant_candidate_sort_key(item: Tuple[str, int, str]) -> Tuple[int, int, int, int]:
    name, size, _ = item
    n = name.lower()
    strong = (
        "proteingroup" in n
        or "protein_group" in n
        or "protein groups" in n
        or "lfq" in n
        or "intensity" in n
        or "abundance" in n
        or n.endswith(".mztab")
        or n.endswith(".mztab.gz")
    )
    kw_hits = sum(1 for k in NAME_KEYWORDS if k in n)
    size_key = size if size >= 0 else 10**18
    return (0 if strong else 1, -kw_hits, size_key, len(n))


def fallback_candidate_sort_key(item: Tuple[str, int, str]) -> Tuple[int, int]:
    name, size, _ = item
    n = name.lower()
    tabular = n.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz", ".mztab", ".mztab.gz", ".xlsx", ".xls"))
    size_key = size if size >= 0 else 10**18
    return (0 if tabular else 1, size_key)


def safe_cache_name(rel_name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", rel_name)


def download_bytes(url: str, cache_path: Optional[Path] = None, timeout: int = 120) -> bytes:
    if cache_path is not None and cache_path.exists() and cache_path.stat().st_size > 0:
        return cache_path.read_bytes()
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    b = r.content
    if cache_path is not None:
        ensure_dir(cache_path.parent)
        cache_path.write_bytes(b)
    return b


def load_table_from_bytes(name: str, b: bytes) -> Optional[pd.DataFrame]:
    name_l = name.lower()

    def strip_leading_hash_comments(raw_bytes: bytes) -> bytes:
        try:
            txt = raw_bytes.decode("utf-8", errors="ignore")
        except Exception:
            return raw_bytes
        lines = txt.splitlines()
        i = 0
        while i < len(lines) and lines[i].lstrip().startswith("#"):
            i += 1
        if i == 0:
            return raw_bytes
        return ("\n".join(lines[i:]) + "\n").encode("utf-8")

    if name_l.endswith(".zip"):
        try:
            with zipfile.ZipFile(io.BytesIO(b)) as zf:
                infos = [zi for zi in zf.infolist() if not zi.is_dir()]

                def score_member(zi: zipfile.ZipInfo) -> tuple:
                    ml = zi.filename.lower()
                    is_table = ml.endswith((".csv", ".tsv", ".txt", ".mztab"))
                    bad = any(x in ml for x in ["analysislog", "parameter", "params", "summary", "readme", "metadata", "peptide", "psm", "fasta"])
                    quantish = any(k in ml for k in NAME_KEYWORDS)
                    too_big = zi.file_size > ZIP_MEMBER_MAX_BYTES
                    return (0 if is_table else 1, 0 if quantish else 1, 1 if bad else 0, 1 if too_big else 0, len(ml))

                for zi in sorted(infos, key=score_member)[:ZIP_MEMBER_MAX_COUNT]:
                    m = zi.filename
                    ml = m.lower()
                    if not ml.endswith((".csv", ".tsv", ".txt", ".mztab")):
                        continue
                    if any(x in ml for x in ["analysislog", "parameter", "params", "summary", "readme", "metadata", "peptide", "psm", "fasta"]):
                        continue
                    if zi.file_size > ZIP_MEMBER_MAX_BYTES:
                        continue
                    raw = strip_leading_hash_comments(zf.read(m))
                    for sep in ["\t", ",", ";"]:
                        try:
                            df = pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)
                            if df.shape[1] >= 2:
                                return df
                        except Exception:
                            continue
        except Exception:
            return None
        return None

    if name_l.endswith((".xlsx", ".xls")):
        try:
            xls = pd.ExcelFile(io.BytesIO(b))
            for sh in xls.sheet_names:
                try:
                    df = pd.read_excel(io.BytesIO(b), sheet_name=sh)
                    if df.shape[1] >= 2:
                        return df
                except Exception:
                    continue
        except Exception:
            return None
        return None

    raw = b
    if name_l.endswith(".gz"):
        try:
            raw = gzip.decompress(b)
        except Exception:
            return None

    raw = strip_leading_hash_comments(raw)
    for sep in ["\t", ",", ";"]:
        try:
            df = pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    return None


def pick_gene_col(df: pd.DataFrame) -> Optional[str]:
    cands = [
        "Gene names", "Gene Names", "Gene.names", "gene_name", "Gene name",
        "Gene", "Genes", "Gene Symbol", "gene_symbol",
        "Description", "description",
        "Accession", "accession", "Protein ID", "Protein IDs", "ProteinID",
        "Protein", "ProteinName", "Majority protein IDs",
        "PG.Genes", "Gene Symbol", "gene_symbol",
    ]
    for c in cands:
        if c in df.columns:
            return c
    for c in df.columns:
        cl = c.lower()
        if "protein group" in cl:
            continue
        if "gene" in cl:
            return c
        if "accession" in cl or "description" in cl:
            return c
        if "protein" in cl and "group" not in cl:
            return c
    return None


def is_probable_uniprot_accession(token: str) -> bool:
    t = token.strip().upper()
    if re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9](?:-\d+)?", t):
        return True
    if re.fullmatch(r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:-\d+)?", t):
        return True
    return False


def split_genes(v: object) -> List[str]:
    if pd.isna(v):
        return []
    s = str(v)

    out: List[str] = []

    # UniProt-style description strings often include `GN=GENE`.
    for m in re.finditer(r"\bGN=([A-Za-z0-9\-]+)", s, flags=re.IGNORECASE):
        g = m.group(1).strip().upper()
        if g and re.fullmatch(r"[A-Z][A-Z0-9\-]{1,20}", g):
            out.append(g)

    toks = re.split(r"[;|,/\\]\s*", s)
    for t in toks:
        tt = t.strip()
        if not tt:
            continue

        # e.g. APOB_HUMAN or APOB_MOUSE -> APOB
        if "_" in tt:
            stem = tt.split("_", 1)[0].strip().upper()
            if stem and re.fullmatch(r"[A-Z][A-Z0-9\-]{1,20}", stem) and not is_probable_uniprot_accession(stem):
                out.append(stem)

        g = tt.upper()
        g = re.sub(r"_[0-9]+$", "", g)
        if g and re.fullmatch(r"[A-Z0-9\-]+", g) and not is_probable_uniprot_accession(g):
            out.append(g)

    # preserve order, dedupe
    out = list(dict.fromkeys(out))
    return out


def pick_intensity_cols(df: pd.DataFrame) -> List[str]:
    cols: List[str] = []
    for c in df.columns:
        cl = c.lower()
        if any(k in cl for k in ["ibaq", "lfq intensity", "intensity", "abundance", "quantity", ".pg.quantity", "pg.quantity"]):
            cols.append(c)
    if cols:
        return cols

    skip = ("gene", "protein", "accession", "description", "group", "id", "score", "qvalue", "fdr")
    for c in df.columns:
        cl = c.lower()
        if any(s in cl for s in skip):
            continue
        n = pd.to_numeric(df[c], errors="coerce")
        if n.notna().mean() > 0.5:
            cols.append(c)
    return cols


def _label_from_text(text: str, disease_hint: str, case_label: str, control_label: str) -> Optional[str]:
    t = text.lower()
    cl = case_label.lower()
    ctl = control_label.lower()

    lbl = pride_classify_text(text, case_lbl=case_label, ctrl_lbl=control_label)
    if lbl in {"case", "control"}:
        return lbl

    has_case = (cl and cl in t) or any(k in t for k in CASE_KEYS)
    has_ctrl = (ctl and ctl in t) or any(k in t for k in CTRL_KEYS)

    if disease_hint == "MS" and re.search(r"(^|[^a-z0-9])ms([^a-z0-9]|$)", t):
        has_case = True
    if disease_hint == "AD" and ("alzheimer" in t or re.search(r"(^|[^a-z0-9])ad([^a-z0-9]|$)", t)):
        has_case = True
    if disease_hint == "PD" and ("parkinson" in t or re.search(r"(^|[^a-z0-9])pd([^a-z0-9]|$)", t)):
        has_case = True

    if has_case and not has_ctrl:
        return "case"
    if has_ctrl and not has_case:
        return "control"
    return None


def _extract_annotation_labels(
    files: List[Tuple[str, int, str]],
    disease_hint: str,
    case_label: str,
    control_label: str,
    cache_dir: Optional[Path],
) -> Dict[str, str]:
    ann_files: List[Tuple[str, int, str]] = []
    for name, size, url in files:
        nl = name.lower()
        if any(k in nl for k in ["annotation", "sample", "metadata", "clinical", "sdrf"]) and \
           nl.endswith((".csv", ".tsv", ".txt", ".xlsx", ".xls", ".zip")) and \
           not any(nl.endswith(ext) for ext in BLOCK_EXT):
            ann_files.append((name, size, url))

    out: Dict[str, str] = {}

    def parse_df(df: pd.DataFrame) -> None:
        if df is None or df.empty:
            return
        cols = list(df.columns)
        sample_col = None
        for c in cols:
            cl = c.lower()
            if "sample" in cl and any(x in cl for x in ["name", "id", "file", "raw"]):
                sample_col = c
                break
        if sample_col is None:
            for c in cols:
                if "sample" in c.lower():
                    sample_col = c
                    break
        if sample_col is None:
            return

        label_cols = [
            c for c in cols
            if any(k in c.lower() for k in ["diagnos", "classification", "group", "condition", "status", "disease", "phenotype", "class"])
        ]
        if not label_cols:
            return

        for _, r in df.iterrows():
            sname = str(r.get(sample_col, "")).strip().lower()
            if not sname or sname == "nan":
                continue
            label_text = " | ".join(str(r.get(c, "")) for c in label_cols)
            lab = _label_from_text(label_text, disease_hint=disease_hint, case_label=case_label, control_label=control_label)
            if lab:
                out[sname] = lab

    for name, _, url in ann_files:
        cache_path = None
        if cache_dir is not None:
            cache_path = cache_dir / safe_cache_name(name)
        try:
            b = download_bytes(url, cache_path=cache_path)
        except Exception:
            continue

        nl = name.lower()
        try:
            if nl.endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(b)) as zf:
                    for m in zf.namelist():
                        ml = m.lower()
                        if m.endswith("/"):
                            continue
                        if ml.endswith((".csv", ".tsv", ".txt")):
                            raw = zf.read(m)
                            for sep in [",", "\t", ";"]:
                                try:
                                    df = pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)
                                    if df.shape[1] > 1:
                                        parse_df(df)
                                        break
                                except Exception:
                                    continue
                        elif ml.endswith((".xlsx", ".xls")):
                            raw = zf.read(m)
                            try:
                                xls = pd.ExcelFile(io.BytesIO(raw))
                                for sh in xls.sheet_names:
                                    df = pd.read_excel(io.BytesIO(raw), sheet_name=sh)
                                    parse_df(df)
                            except Exception:
                                pass
            elif nl.endswith((".xlsx", ".xls")):
                xls = pd.ExcelFile(io.BytesIO(b))
                for sh in xls.sheet_names:
                    df = pd.read_excel(io.BytesIO(b), sheet_name=sh)
                    parse_df(df)
            elif nl.endswith((".csv", ".tsv", ".txt")):
                for sep in [",", "\t", ";"]:
                    try:
                        df = pd.read_csv(io.BytesIO(b), sep=sep, low_memory=False)
                        if df.shape[1] > 1:
                            parse_df(df)
                            break
                    except Exception:
                        continue
        except Exception:
            continue

    return out


def load_manual_sample_labels(accession: str, sample_map_dir: Optional[Path]) -> Dict[str, str]:
    if sample_map_dir is None:
        return {}
    path = sample_map_dir / f"{accession}.csv"
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path)
    except Exception:
        return {}
    if df.empty:
        return {}

    cols_lower = {c.lower(): c for c in df.columns}
    sample_col = cols_lower.get("sample")
    label_col = cols_lower.get("label")
    if not sample_col or not label_col:
        return {}

    out: Dict[str, str] = {}
    for _, r in df.iterrows():
        s = str(r.get(sample_col, "")).strip().lower()
        l = str(r.get(label_col, "")).strip().lower()
        if not s:
            continue
        if l in {"case", "disease", "patient"}:
            out[s] = "case"
        elif l in {"control", "healthy", "ctrl"}:
            out[s] = "control"
    return out


def split_case_control(
    cols: List[str],
    disease_hint: str,
    case_label: str,
    control_label: str,
    sample_labels: Optional[Dict[str, str]] = None,
) -> Tuple[List[str], List[str]]:
    case: List[str] = []
    ctrl: List[str] = []
    for c in cols:
        cl = c.lower()
        lbl = pride_classify_text(c, case_lbl=case_label, ctrl_lbl=control_label)
        if lbl == "case":
            case.append(c)
            continue
        if lbl == "control":
            ctrl.append(c)
            continue

        has_case = any(k in cl for k in CASE_KEYS)
        has_ctrl = any(k in cl for k in CTRL_KEYS)

        if disease_hint == "MS" and re.search(r"(^|[^a-z0-9])ms([^a-z0-9]|$)", cl):
            has_case = True
        if disease_hint == "AD" and ("alzheimer" in cl or re.search(r"(^|[^a-z0-9])ad([^a-z0-9]|$)", cl)):
            has_case = True
        if disease_hint == "PD" and ("parkinson" in cl or re.search(r"(^|[^a-z0-9])pd([^a-z0-9]|$)", cl)):
            has_case = True

        if has_case and not has_ctrl:
            case.append(c)
        elif has_ctrl and not has_case:
            ctrl.append(c)

    if sample_labels and (not case or not ctrl):
        for c in cols:
            cl = c.lower()
            lab = None
            for sname, l in sample_labels.items():
                if sname in cl or cl in sname:
                    lab = l
                    break
            if lab == "case" and c not in case:
                case.append(c)
            elif lab == "control" and c not in ctrl:
                ctrl.append(c)

    # preserve order and dedupe
    case = list(dict.fromkeys(case))
    ctrl = list(dict.fromkeys(ctrl))
    return case, ctrl


def infer_disease_hint(disease: str) -> str:
    d = disease.strip().lower()
    if "alz" in d or d == "ad":
        return "AD"
    if "multiple" in d or d == "ms":
        return "MS"
    return "PD"


def normalize_sample_values(values: List[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    # Most proteomics intensities are non-negative; allow zeros and stabilize with +1.
    if np.all(arr >= 0):
        return np.log2(arr + 1.0)
    return arr


def sample_level_group_values(sub_df: pd.DataFrame, cols: List[str]) -> List[float]:
    vals: List[float] = []
    for c in cols:
        v = pd.to_numeric(sub_df[c], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if len(v) == 0:
            continue
        vals.append(float(np.nanmedian(v.values.astype(float))))
    return vals


def ci_from_effect_se(effect: float, se: float) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if pd.isna(effect) or pd.isna(se) or se <= 0:
        return None, None, None
    z = float(effect) / float(se)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    lo = float(effect) - 1.96 * float(se)
    hi = float(effect) + 1.96 * float(se)
    return lo, hi, p


def summarize_gene(
    gene: str,
    sub_df: pd.DataFrame,
    case_cols: List[str],
    ctrl_cols: List[str],
    threshold: float,
    expected: str,
) -> Dict[str, object]:
    case_vals_raw = sample_level_group_values(sub_df, case_cols)
    ctrl_vals_raw = sample_level_group_values(sub_df, ctrl_cols)

    case_vals = normalize_sample_values(case_vals_raw)
    ctrl_vals = normalize_sample_values(ctrl_vals_raw)

    n_case = int(case_vals.size)
    n_control = int(ctrl_vals.size)
    present = n_case > 0 and n_control > 0

    if not present:
        return {
            "gene": gene,
            "present": False,
            "n_case": n_case,
            "n_control": n_control,
            "mean_case": None,
            "mean_ctrl": None,
            "log2fc": None,
            "se": None,
            "ci95_low": None,
            "ci95_high": None,
            "p_value": None,
            "direction": "absent",
            "expected": expected,
            "concordant": "na",
        }

    mean_case = float(np.mean(case_vals))
    mean_ctrl = float(np.mean(ctrl_vals))
    effect = float(mean_case - mean_ctrl)

    se = np.nan
    if n_case >= 2 and n_control >= 2:
        var_case = float(np.var(case_vals, ddof=1))
        var_ctrl = float(np.var(ctrl_vals, ddof=1))
        if var_case >= 0 and var_ctrl >= 0:
            se = math.sqrt(var_case / n_case + var_ctrl / n_control)

    lo, hi, p = ci_from_effect_se(effect, se)
    obs_dir = direction_from_fc(effect, threshold=threshold)
    conc = concordance(obs_dir, expected)

    return {
        "gene": gene,
        "present": True,
        "n_case": n_case,
        "n_control": n_control,
        "mean_case": round(mean_case, 6),
        "mean_ctrl": round(mean_ctrl, 6),
        "log2fc": round(effect, 6),
        "se": None if pd.isna(se) else round(float(se), 6),
        "ci95_low": None if lo is None else round(float(lo), 6),
        "ci95_high": None if hi is None else round(float(hi), 6),
        "p_value": None if p is None else float(p),
        "direction": obs_dir,
        "expected": expected,
        "concordant": conc,
    }


def parse_quant_file(
    accession: str,
    disease: str,
    biospecimen: str,
    case_label: str,
    control_label: str,
    expected_map: Dict[str, str],
    threshold: float,
    file_name: str,
    file_url: str,
    file_size: int,
    all_files: List[Tuple[str, int, str]],
    cache_dir: Optional[Path],
    manual_sample_labels: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    cache_path = None
    if cache_dir is not None:
        cache_path = cache_dir / accession / safe_cache_name(file_name)

    try:
        b = download_bytes(file_url, cache_path=cache_path)
        df = load_table_from_bytes(file_name, b)
        if df is None or df.empty:
            return [], {
                "status": "parse_failed",
                "selected_file": file_name,
                "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
            }

        gene_col = pick_gene_col(df)
        if not gene_col:
            return [], {
                "status": "blocked_no_id_mapping",
                "selected_file": file_name,
                "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
            }

        intensity_cols = pick_intensity_cols(df)
        if not intensity_cols:
            return [], {
                "status": "blocked_no_quant_table",
                "selected_file": file_name,
                "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
            }

        # Map rows to genes once, then aggregate by gene.
        gene_to_rows: Dict[str, List[int]] = {g: [] for g in PANEL}
        for idx, row in df.iterrows():
            genes = split_genes(row.get(gene_col))
            if not genes:
                continue
            gset = set(genes)
            for g in PANEL:
                if g in gset:
                    gene_to_rows[g].append(idx)

        if not any(gene_to_rows[g] for g in PANEL):
            return [], {
                "status": "blocked_no_id_mapping",
                "selected_file": file_name,
                "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
            }

        disease_hint = infer_disease_hint(disease)
        sample_labels = _extract_annotation_labels(
            files=all_files,
            disease_hint=disease_hint,
            case_label=case_label,
            control_label=control_label,
            cache_dir=(cache_dir / accession) if cache_dir is not None else None,
        )
        merged_labels: Dict[str, str] = {}
        if sample_labels:
            merged_labels.update({str(k).lower(): v for k, v in sample_labels.items()})
        if manual_sample_labels:
            merged_labels.update({str(k).lower(): v for k, v in manual_sample_labels.items()})

        case_cols, ctrl_cols = split_case_control(
            cols=intensity_cols,
            disease_hint=disease_hint,
            case_label=case_label,
            control_label=control_label,
            sample_labels=merged_labels if merged_labels else None,
        )

        if not case_cols or not ctrl_cols:
            return [], {
                "status": "blocked_no_group_labels",
                "selected_file": file_name,
                "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
            }

        rows_out: List[Dict[str, object]] = []
        n_case_max = 0
        n_ctrl_max = 0
        n_present = 0

        for gene in PANEL:
            idxs = gene_to_rows.get(gene, [])
            expected = expected_map.get(gene, "none")
            if not idxs:
                stats = {
                    "gene": gene,
                    "present": False,
                    "n_case": 0,
                    "n_control": 0,
                    "mean_case": None,
                    "mean_ctrl": None,
                    "log2fc": None,
                    "se": None,
                    "ci95_low": None,
                    "ci95_high": None,
                    "p_value": None,
                    "direction": "absent",
                    "expected": expected,
                    "concordant": "na",
                }
            else:
                sub_df = df.loc[idxs, intensity_cols]
                stats = summarize_gene(
                    gene=gene,
                    sub_df=sub_df,
                    case_cols=case_cols,
                    ctrl_cols=ctrl_cols,
                    threshold=threshold,
                    expected=expected,
                )

            if stats["present"]:
                n_present += 1
                n_case_max = max(n_case_max, int(stats["n_case"]))
                n_ctrl_max = max(n_ctrl_max, int(stats["n_control"]))

            rows_out.append(
                {
                    "accession": accession,
                    "disease": disease,
                    "biospecimen": biospecimen,
                    "source": "pride",
                    "file_name": file_name,
                    **stats,
                }
            )

        coverage = n_present / float(len(PANEL))
        status = "ok" if n_present >= 2 else "low_coverage"

        return rows_out, {
            "status": status,
            "selected_file": file_name,
            "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
            "n_case": n_case_max,
            "n_control": n_ctrl_max,
            "n_panel_genes_present": n_present,
            "coverage_rate": round(coverage, 4),
        }

    except Exception:
        return [], {
            "status": "parse_failed",
            "selected_file": file_name,
            "file_size_mb": round(file_size / (1024 * 1024), 2) if file_size > 0 else None,
            "n_case": 0,
            "n_control": 0,
            "n_panel_genes_present": 0,
            "coverage_rate": 0.0,
        }


def pick_best_parse(candidates: List[Tuple[str, int, str]], parsed_items: List[Tuple[List[Dict[str, object]], Dict[str, object]]]) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    best_rows: List[Dict[str, object]] = []
    best_qc: Dict[str, object] = {
        "status": "blocked_no_quant_table",
        "selected_file": "",
        "file_size_mb": None,
        "n_case": 0,
        "n_control": 0,
        "n_panel_genes_present": 0,
        "coverage_rate": 0.0,
    }

    for rows, qc in parsed_items:
        cur_status = str(qc.get("status", "parse_failed"))
        cur_priority = STATUS_PRIORITY.get(cur_status, -1)
        best_priority = STATUS_PRIORITY.get(str(best_qc.get("status", "parse_failed")), -1)

        better = False
        if cur_priority > best_priority:
            better = True
        elif cur_priority == best_priority:
            cur_found = int(qc.get("n_panel_genes_present", 0))
            best_found = int(best_qc.get("n_panel_genes_present", 0))
            if cur_found > best_found:
                better = True

        if better:
            best_rows = rows
            best_qc = qc

        if cur_status == "ok":
            break

    return best_rows, best_qc


def random_effects_meta(yi: np.ndarray, sei: np.ndarray) -> Dict[str, float]:
    vi = sei ** 2
    k = yi.size
    if k == 1:
        mu = float(yi[0])
        se = float(sei[0])
        lo, hi, p = ci_from_effect_se(mu, se)
        return {
            "k": 1,
            "meta_log2fc": mu,
            "meta_se": se,
            "ci95_low": lo,
            "ci95_high": hi,
            "tau2": 0.0,
            "I2": 0.0,
            "Q": 0.0,
            "p_value": p,
        }

    wi = 1.0 / vi
    mu_fixed = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * ((yi - mu_fixed) ** 2)))
    df = k - 1
    c = float(np.sum(wi) - (np.sum(wi ** 2) / np.sum(wi)))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0

    wi_star = 1.0 / (vi + tau2)
    mu = float(np.sum(wi_star * yi) / np.sum(wi_star))
    se = math.sqrt(1.0 / float(np.sum(wi_star)))
    lo, hi, p = ci_from_effect_se(mu, se)
    i2 = max(0.0, ((q - df) / q) * 100.0) if q > 0 else 0.0

    return {
        "k": int(k),
        "meta_log2fc": mu,
        "meta_se": se,
        "ci95_low": lo,
        "ci95_high": hi,
        "tau2": float(tau2),
        "I2": float(i2),
        "Q": float(q),
        "p_value": p,
    }


def build_meta(
    effects_df: pd.DataFrame,
    expected_map: Dict[str, str],
    threshold: float,
    compartment: str,
) -> pd.DataFrame:
    meta_cols = [
        "compartment", "disease_scope", "gene", "k", "meta_log2fc", "meta_se",
        "ci95_low", "ci95_high", "tau2", "I2", "Q", "p_value", "direction", "expected", "concordant",
    ]
    required_cols = {"biospecimen", "analysis_tier", "present", "se", "disease", "gene", "log2fc"}
    if effects_df.empty or not required_cols.issubset(set(effects_df.columns)):
        return pd.DataFrame(columns=meta_cols)

    sub = effects_df.copy()
    sub = sub[sub["biospecimen"].str.lower() == compartment.lower()]
    sub = sub[sub["analysis_tier"] == "primary"]
    sub = sub[sub["present"] == True]  # noqa: E712
    sub = sub[pd.to_numeric(sub["se"], errors="coerce").notna()]
    sub = sub[pd.to_numeric(sub["se"], errors="coerce") > 0]

    rows: List[Dict[str, object]] = []
    if sub.empty:
        return pd.DataFrame(columns=meta_cols)

    diseases = sorted(sub["disease"].dropna().astype(str).unique().tolist())
    scopes = diseases + ["ALL"]

    for scope in scopes:
        for gene in PANEL:
            if scope == "ALL":
                g = sub[sub["gene"] == gene]
            else:
                g = sub[(sub["gene"] == gene) & (sub["disease"] == scope)]
            if g.empty:
                continue

            yi = pd.to_numeric(g["log2fc"], errors="coerce").values.astype(float)
            sei = pd.to_numeric(g["se"], errors="coerce").values.astype(float)
            keep = np.isfinite(yi) & np.isfinite(sei) & (sei > 0)
            yi = yi[keep]
            sei = sei[keep]
            if yi.size == 0:
                continue

            meta = random_effects_meta(yi, sei)
            exp = expected_map.get(gene, "none")
            direction = direction_from_fc(meta["meta_log2fc"], threshold=threshold)
            conc = concordance(direction, exp)

            rows.append(
                {
                    "compartment": compartment,
                    "disease_scope": scope,
                    "gene": gene,
                    "k": meta["k"],
                    "meta_log2fc": round(float(meta["meta_log2fc"]), 6),
                    "meta_se": round(float(meta["meta_se"]), 6),
                    "ci95_low": None if meta["ci95_low"] is None else round(float(meta["ci95_low"]), 6),
                    "ci95_high": None if meta["ci95_high"] is None else round(float(meta["ci95_high"]), 6),
                    "tau2": round(float(meta["tau2"]), 6),
                    "I2": round(float(meta["I2"]), 4),
                    "Q": round(float(meta["Q"]), 6),
                    "p_value": None if meta["p_value"] is None else float(meta["p_value"]),
                    "direction": direction,
                    "expected": exp,
                    "concordant": conc,
                }
            )

    return pd.DataFrame(rows)


def write_summary_md(
    path: Path,
    qc_df: pd.DataFrame,
    effects_df: pd.DataFrame,
    meta_plasma_df: pd.DataFrame,
    meta_serum_df: pd.DataFrame,
) -> None:
    ensure_dir(path.parent)
    lines: List[str] = []
    lines.append("# Blood Primary10 Validation Summary")
    lines.append("")
    lines.append("Primary blood proteomics external-validation run across 10 preselected AD/PD/MS datasets.")
    lines.append("")

    lines.append("## Dataset QC")
    status_counts = qc_df.groupby("status").size().sort_values(ascending=False)
    for status, n in status_counts.items():
        lines.append(f"- `{status}`: **{int(n)}**")
    lines.append("")

    lines.append("## Dataset Table")
    lines.append("| accession | disease | biospecimen | status | genes_present | n_case | n_control | selected_file |")
    lines.append("|---|---|---|---|---:|---:|---:|---|")
    for _, r in qc_df.sort_values(["disease", "accession"]).iterrows():
        lines.append(
            f"| {r['accession']} | {r['disease']} | {r['biospecimen']} | {r['status']} | {int(r['n_panel_genes_present'])} | {int(r['n_case']) if pd.notna(r['n_case']) else 0} | {int(r['n_control']) if pd.notna(r['n_control']) else 0} | {r['selected_file']} |"
        )
    lines.append("")

    if not effects_df.empty:
        p = effects_df[effects_df["present"] == True]  # noqa: E712
        conc = p[p["concordant"].isin(["yes", "no"])].groupby("concordant").size()
        lines.append("## Effect Summary")
        lines.append(f"- Total gene-dataset rows: **{len(effects_df)}**")
        lines.append(f"- Present rows: **{len(p)}**")
        lines.append(f"- Concordant `yes`: **{int(conc.get('yes', 0))}**")
        lines.append(f"- Concordant `no`: **{int(conc.get('no', 0))}**")
        lines.append("")

    lines.append("## Meta-analysis")
    lines.append(f"- Plasma meta rows: **{len(meta_plasma_df)}**")
    lines.append(f"- Serum meta rows: **{len(meta_serum_df)}**")
    lines.append("")

    lines.append("## Output Files")
    lines.append("- `per_dataset_effects.csv`")
    lines.append("- `per_dataset_qc.csv`")
    lines.append("- `meta_plasma.csv`")
    lines.append("- `meta_serum.csv`")
    lines.append("- `summary_primary10.md`")
    lines.append("- `selected_file_inventory.csv`")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()

    manifest_path = Path(args.manifest).resolve()
    targets_path = Path(args.targets).resolve()
    out_dir = Path(args.out_dir).resolve()
    sample_map_dir = Path(args.sample_map_dir).resolve() if args.sample_map_dir else None

    external_root = choose_external_root()
    raw_dir = Path(args.raw_dir).resolve() if args.raw_dir else (external_root / "pride_blood_primary10_raw")

    ensure_dir(out_dir)
    ensure_dir(raw_dir)

    manifest_df = read_manifest(manifest_path)
    if args.accessions:
        keep = {a.strip().upper() for a in args.accessions if str(a).strip()}
        manifest_df = manifest_df[manifest_df["accession"].isin(keep)].copy()

    expected_map = read_expected_chronic(targets_path)

    qc_rows: List[Dict[str, object]] = []
    effects_rows: List[Dict[str, object]] = []
    inventory_rows: List[Dict[str, object]] = []

    total = len(manifest_df)
    for i, row in manifest_df.iterrows():
        accession = str(row["accession"]).strip().upper()
        disease = str(row["disease"]).strip()
        biospecimen = str(row["biospecimen"]).strip().lower()
        source = str(row["source"]).strip().lower()
        case_label = str(row["case_label"]).strip()
        control_label = str(row["control_label"]).strip()
        manual_labels = load_manual_sample_labels(accession, sample_map_dir)

        print(f"[{len(qc_rows)+1}/{total}] {accession} ({disease}, {biospecimen})")

        if source not in {"pride", "iprox", "jpost"}:
            qc_rows.append({
                "accession": accession,
                "disease": disease,
                "biospecimen": biospecimen,
                "source": source,
                "status": "blocked_no_quant_table",
                "selected_file": "",
                "file_size_mb": None,
                "n_files_total": 0,
                "n_candidate_files": 0,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
                "analysis_tier": "blocked",
                "stop_rule": "blocked_no_quant_table",
            })
            continue

        pride_root, files = list_pride_files(accession)
        used_deep_probe = False
        used_source_probe = ""
        quant_files = pick_quant_candidates(
            files=files,
            max_file_bytes=args.max_file_bytes,
            max_files_per_accession=args.max_files_per_accession,
        )

        probe_files: List[Tuple[str, int, str]] = []
        if not quant_files and files:
            probe_pool = [(n, s, u) for (n, s, u) in files if is_fallback_parse_candidate(n, s)]
            probe_pool.sort(key=fallback_candidate_sort_key)
            probe_files = probe_pool[:FALLBACK_MAX_FILES]
            quant_files = probe_files

        if not quant_files and pride_root and (len(files) == 0 or len(files) <= 120):
            deep_files = list_pride_files_deep(pride_root, timeout=40)
            if deep_files:
                files = merge_file_lists(files, deep_files)
                used_deep_probe = True
                quant_files = pick_quant_candidates(
                    files=files,
                    max_file_bytes=args.max_file_bytes,
                    max_files_per_accession=args.max_files_per_accession,
                )
                if not quant_files:
                    probe_pool = [(n, s, u) for (n, s, u) in files if is_fallback_parse_candidate(n, s)]
                    probe_pool.sort(key=fallback_candidate_sort_key)
                    probe_files = probe_pool[:FALLBACK_MAX_FILES]
                    quant_files = probe_files

        if not quant_files and source == "iprox":
            iprox_files = list_iprox_files(accession, timeout=40)
            if iprox_files:
                files = merge_file_lists(files, iprox_files)
                used_source_probe = "iprox"
                quant_files = pick_quant_candidates(
                    files=files,
                    max_file_bytes=args.max_file_bytes,
                    max_files_per_accession=args.max_files_per_accession,
                )
                if not quant_files:
                    probe_pool = [(n, s, u) for (n, s, u) in files if is_fallback_parse_candidate(n, s)]
                    probe_pool.sort(key=fallback_candidate_sort_key)
                    probe_files = probe_pool[:FALLBACK_MAX_FILES]
                    quant_files = probe_files

        if not quant_files and source == "jpost":
            jpost_files = list_jpost_files(accession, timeout=40)
            if jpost_files:
                files = merge_file_lists(files, jpost_files)
                used_source_probe = "jpost"
                quant_files = pick_quant_candidates(
                    files=files,
                    max_file_bytes=args.max_file_bytes,
                    max_files_per_accession=args.max_files_per_accession,
                )
                if not quant_files:
                    probe_pool = [(n, s, u) for (n, s, u) in files if is_fallback_parse_candidate(n, s)]
                    probe_pool.sort(key=fallback_candidate_sort_key)
                    probe_files = probe_pool[:FALLBACK_MAX_FILES]
                    quant_files = probe_files

        inventory_rows.append({
            "accession": accession,
            "disease": disease,
            "biospecimen": biospecimen,
            "n_files_total": len(files),
            "n_candidate_files": len(quant_files),
            "used_fallback_probe": bool(probe_files),
            "used_deep_probe": used_deep_probe,
            "used_source_probe": used_source_probe,
            "candidate_files": " | ".join([f[0] for f in quant_files[:10]]),
        })

        if not quant_files:
            for g in PANEL:
                effects_rows.append({
                    "accession": accession,
                    "disease": disease,
                    "biospecimen": biospecimen,
                    "source": source,
                    "file_name": "",
                    "gene": g,
                    "present": False,
                    "n_case": 0,
                    "n_control": 0,
                    "mean_case": None,
                    "mean_ctrl": None,
                    "log2fc": None,
                    "se": None,
                    "ci95_low": None,
                    "ci95_high": None,
                    "p_value": None,
                    "direction": "absent",
                    "expected": expected_map.get(g, "none"),
                    "concordant": "na",
                    "dataset_status": "blocked_no_quant_table",
                    "analysis_tier": "blocked",
                })
            qc_rows.append({
                "accession": accession,
                "disease": disease,
                "biospecimen": biospecimen,
                "source": source,
                "status": "blocked_no_quant_table",
                "selected_file": "",
                "file_size_mb": None,
                "n_files_total": len(files),
                "n_candidate_files": 0,
                "n_case": 0,
                "n_control": 0,
                "n_panel_genes_present": 0,
                "coverage_rate": 0.0,
                "analysis_tier": "blocked",
                "stop_rule": "blocked_no_quant_table",
            })
            continue

        parsed_items: List[Tuple[List[Dict[str, object]], Dict[str, object]]] = []
        for n, s, u in quant_files:
            rows_out, qc = parse_quant_file(
                accession=accession,
                disease=disease,
                biospecimen=biospecimen,
                case_label=case_label,
                control_label=control_label,
                expected_map=expected_map,
                threshold=args.threshold,
                file_name=n,
                file_url=u,
                file_size=s,
                all_files=files,
                cache_dir=raw_dir,
                manual_sample_labels=manual_labels,
            )
            parsed_items.append((rows_out, qc))
            if qc.get("status") == "ok":
                break
            time.sleep(args.sleep_seconds)

        best_rows, best_qc = pick_best_parse(quant_files, parsed_items)
        status = str(best_qc.get("status", "parse_failed"))

        if status == "ok":
            analysis_tier = "primary"
            stop_rule = "passed"
        elif status == "low_coverage":
            analysis_tier = "supplementary"
            stop_rule = "low_coverage"
        elif status == "blocked_no_id_mapping":
            analysis_tier = "blocked"
            stop_rule = "blocked_no_id_mapping"
        elif status == "blocked_no_group_labels":
            analysis_tier = "blocked"
            stop_rule = "blocked_no_group_labels"
        elif status == "blocked_no_quant_table":
            analysis_tier = "blocked"
            stop_rule = "blocked_no_quant_table"
        else:
            analysis_tier = "blocked"
            stop_rule = status

        if best_rows:
            for r in best_rows:
                effects_rows.append({
                    **r,
                    "dataset_status": status,
                    "analysis_tier": analysis_tier,
                })
        else:
            for g in PANEL:
                effects_rows.append({
                    "accession": accession,
                    "disease": disease,
                    "biospecimen": biospecimen,
                    "source": source,
                    "file_name": str(best_qc.get("selected_file", "")),
                    "gene": g,
                    "present": False,
                    "n_case": 0,
                    "n_control": 0,
                    "mean_case": None,
                    "mean_ctrl": None,
                    "log2fc": None,
                    "se": None,
                    "ci95_low": None,
                    "ci95_high": None,
                    "p_value": None,
                    "direction": "absent",
                    "expected": expected_map.get(g, "none"),
                    "concordant": "na",
                    "dataset_status": status,
                    "analysis_tier": analysis_tier,
                })

        qc_rows.append({
            "accession": accession,
            "disease": disease,
            "biospecimen": biospecimen,
            "source": source,
            "status": status,
            "selected_file": best_qc.get("selected_file", ""),
            "file_size_mb": best_qc.get("file_size_mb", None),
            "n_files_total": len(files),
            "n_candidate_files": len(quant_files),
            "n_case": int(best_qc.get("n_case", 0)),
            "n_control": int(best_qc.get("n_control", 0)),
            "n_panel_genes_present": int(best_qc.get("n_panel_genes_present", 0)),
            "coverage_rate": float(best_qc.get("coverage_rate", 0.0)),
            "analysis_tier": analysis_tier,
            "stop_rule": stop_rule,
        })

    effects_df = pd.DataFrame(effects_rows)
    qc_df = pd.DataFrame(qc_rows)
    inv_df = pd.DataFrame(inventory_rows)

    # Plasma and serum random-effects meta-analysis, primary tier only.
    meta_plasma_df = build_meta(effects_df, expected_map=expected_map, threshold=args.threshold, compartment="plasma")
    meta_serum_df = build_meta(effects_df, expected_map=expected_map, threshold=args.threshold, compartment="serum")

    effects_fields = [
        "accession", "disease", "biospecimen", "source", "file_name", "gene", "present", "n_case", "n_control",
        "mean_case", "mean_ctrl", "log2fc", "se", "ci95_low", "ci95_high", "p_value", "direction", "expected",
        "concordant", "dataset_status", "analysis_tier",
    ]
    qc_fields = [
        "accession", "disease", "biospecimen", "source", "status", "selected_file", "file_size_mb", "n_files_total",
        "n_candidate_files", "n_case", "n_control", "n_panel_genes_present", "coverage_rate", "analysis_tier", "stop_rule",
    ]
    meta_fields = [
        "compartment", "disease_scope", "gene", "k", "meta_log2fc", "meta_se", "ci95_low", "ci95_high",
        "tau2", "I2", "Q", "p_value", "direction", "expected", "concordant",
    ]
    inv_fields = [
        "accession", "disease", "biospecimen", "n_files_total", "n_candidate_files",
        "used_fallback_probe", "used_deep_probe", "used_source_probe", "candidate_files",
    ]

    write_csv(out_dir / "per_dataset_effects.csv", effects_df.to_dict("records"), effects_fields)
    write_csv(out_dir / "per_dataset_qc.csv", qc_df.to_dict("records"), qc_fields)
    write_csv(out_dir / "meta_plasma.csv", meta_plasma_df.to_dict("records"), meta_fields)
    write_csv(out_dir / "meta_serum.csv", meta_serum_df.to_dict("records"), meta_fields)
    write_csv(out_dir / "selected_file_inventory.csv", inv_df.to_dict("records"), inv_fields)

    write_summary_md(
        path=out_dir / "summary_primary10.md",
        qc_df=qc_df,
        effects_df=effects_df,
        meta_plasma_df=meta_plasma_df,
        meta_serum_df=meta_serum_df,
    )

    print(f"Wrote: {out_dir / 'per_dataset_effects.csv'}")
    print(f"Wrote: {out_dir / 'per_dataset_qc.csv'}")
    print(f"Wrote: {out_dir / 'meta_plasma.csv'}")
    print(f"Wrote: {out_dir / 'meta_serum.csv'}")
    print(f"Wrote: {out_dir / 'summary_primary10.md'}")
    print(f"Wrote: {out_dir / 'selected_file_inventory.csv'}")


if __name__ == "__main__":
    main()
