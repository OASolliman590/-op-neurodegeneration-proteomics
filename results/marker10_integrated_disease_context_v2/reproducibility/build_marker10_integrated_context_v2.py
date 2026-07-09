#!/usr/bin/env python3
"""Build Marker10 integrated disease-context v2 figures.

V2 changes the figure logic:
- Panel A is a cohort-level Marker10 signature cluster, not a CKD reference panel.
- Panel B is a readable standalone concordance heatmap with disease/accession/n.
- Panel C saves per-disease marker pathway-proximity maps and an overlap summary.
- Panel E is an enhanced pathway/database-role map with definitions.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

import pandas as pd
from PIL import Image, ImageChops, ImageDraw, ImageFont


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[5]
OP = ROOT / "op_external_validation"
OUT = OP / "results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709"
FIG = OUT / "figures"
TABLES = OUT / "tables"
REPRO = OUT / "reproducibility"
PANEL_A_DOMAINS = FIG / "panel_A_domain_signature_clusters"
PANEL_C_MAPS = FIG / "panel_C_disease_marker_proximity_maps"

MARKER_XLSX = OP / "log fold change 10 Markers.xlsx"
DE_TABLE = OP / "results/all_cohorts_moderated_tiered_20260709/tables/all_harmonized_de_rows.tsv"
SUMMARY_TABLE = OP / "results/all_cohorts_moderated_tiered_20260709/tables/all_contrasts_tier_summary.tsv"
STRING_DIR = OP / "results/de_visualization_hub_network/string_marker10_disease_roles/tables"
STRING_EDGES = STRING_DIR / "string_marker10_network_edges.tsv"
STRING_ENRICHMENT = STRING_DIR / "string_marker10_enrichment_raw.tsv"


for path in [OUT, FIG, TABLES, REPRO, PANEL_A_DOMAINS, PANEL_C_MAPS]:
    path.mkdir(parents=True, exist_ok=True)


COLORS = {
    "ink": "#1F2933",
    "muted": "#5B6770",
    "grid": "#D8DEE6",
    "concordant_dep": "#009E73",
    "discordant_dep": "#D55E00",
    "concordant_exploratory": "#8FD1B3",
    "discordant_exploratory": "#E9A37C",
    "weak_same_direction": "#D9F0E4",
    "weak_opposite_direction": "#F6D9CA",
    "directionless": "#E5E7EB",
    "not_detected_or_unmapped": "#FFFFFF",
    "ref_up": "#C2410C",
    "ref_down": "#2563EB",
}

STATUS_ORDER = [
    "concordant_dep",
    "discordant_dep",
    "concordant_exploratory",
    "discordant_exploratory",
    "weak_same_direction",
    "weak_opposite_direction",
    "directionless",
    "not_detected_or_unmapped",
]

STATUS_LABELS = {
    "concordant_dep": "Concordant DEP",
    "discordant_dep": "Discordant DEP",
    "concordant_exploratory": "Concordant exploratory",
    "discordant_exploratory": "Discordant exploratory",
    "weak_same_direction": "Same direction, weak",
    "weak_opposite_direction": "Opposite direction, weak",
    "directionless": "No direction",
    "not_detected_or_unmapped": "Not detected",
}

DOMAIN_COLORS = {
    "Nervous system": "#4C78A8",
    "Metabolic/endocrine": "#F58518",
    "Kidney/urinary": "#54A24B",
    "Hepatic": "#B279A2",
    "Cardiovascular/inflammatory": "#E45756",
    "Other": "#8C8C8C",
}


def font(size: int, bold: bool = False):
    candidates = []
    if bold:
        candidates += [
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "/Library/Fonts/Arial Bold.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        ]
    candidates += [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


FONTS = {
    "title": font(34, True),
    "subtitle": font(24, True),
    "body": font(20),
    "small": font(16),
    "tiny": font(13),
    "panel": font(22, True),
}

STATUS_CODES = {
    "concordant_dep": "C",
    "discordant_dep": "D",
    "concordant_exploratory": "c",
    "discordant_exploratory": "d",
    "weak_same_direction": "w+",
    "weak_opposite_direction": "w-",
    "directionless": "0",
    "not_detected_or_unmapped": ".",
}

DOMAIN_ORDER = [
    "Nervous system",
    "Metabolic/endocrine",
    "Kidney/urinary",
    "Hepatic",
    "Cardiovascular/inflammatory",
    "Other",
]

ROLE_CATEGORY_COLORS = {
    "NetworkNeighborAL": "#6A4C93",
    "KEGG": "#1982C4",
    "Component": "#8AC926",
    "COMPARTMENTS": "#52B788",
    "TISSUES": "#FFCA3A",
    "PMID": "#FF595E",
}

MARKER_POSITIONS = {
    "ACTG1": (0.14, 0.40),
    "VWF": (0.50, 0.24),
    "CD44": (0.22, 0.58),
    "C4B": (0.70, 0.34),
    "CFHR2": (0.82, 0.52),
    "DNAH9": (0.50, 0.56),
    "GPX3": (0.16, 0.76),
    "LRG1": (0.42, 0.80),
    "ITIH3": (0.64, 0.82),
    "MYH7B": (0.88, 0.78),
}


def save_with_pdf(img: Image.Image, path: Path) -> None:
    img.save(path, dpi=(300, 300))
    img.convert("RGB").save(path.with_suffix(".pdf"), "PDF", resolution=300)


def crop_content(im: Image.Image, pad: int = 28) -> Image.Image:
    bg = Image.new(im.mode, im.size, "white")
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox is None:
        return im
    x0, y0, x1, y1 = bbox
    x0 = max(0, x0 - pad)
    y0 = max(0, y0 - pad)
    x1 = min(im.size[0], x1 + pad)
    y1 = min(im.size[1], y1 + pad)
    return im.crop((x0, y0, x1, y1))


def code_color(status: str) -> str:
    return COLORS.get(status, "#FFFFFF")


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def blend(color_a: str, color_b: str, alpha: float) -> str:
    a = hex_to_rgb(color_a)
    b = hex_to_rgb(color_b)
    mixed = tuple(int(round((1 - alpha) * x + alpha * y)) for x, y in zip(a, b))
    return "#%02x%02x%02x" % mixed


def text_size(draw: ImageDraw.ImageDraw, text: str, font) -> tuple[int, int]:
    bbox = draw.textbbox((0, 0), str(text), font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def wrapped(draw: ImageDraw.ImageDraw, xy, text, width, font, fill=COLORS["muted"], max_lines=None) -> int:
    words = str(text).split()
    lines = []
    current = ""
    for word in words:
        proposal = word if not current else f"{current} {word}"
        if text_size(draw, proposal, font)[0] <= width:
            current = proposal
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if max_lines is not None and len(lines) > max_lines:
        lines = lines[:max_lines]
        while lines and text_size(draw, lines[-1] + "...", font)[0] > width:
            lines[-1] = lines[-1][:-1]
        if lines:
            lines[-1] += "..."
    x, y = xy
    line_h = text_size(draw, "Ag", font)[1] + 5
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += line_h
    return y


def panel_header(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, title: str, subtitle: str, letter: str) -> int:
    draw.rounded_rectangle([x, y, x + 34, y + 30], radius=5, fill=COLORS["ink"])
    draw.text((x + 10, y + 3), letter, font=FONTS["panel"], fill="white")
    draw.text((x + 44, y), title, font=FONTS["panel"], fill=COLORS["ink"])
    yy = wrapped(draw, (x + 44, y + 32), subtitle, w - 44, FONTS["small"], fill=COLORS["muted"], max_lines=2)
    return yy + 8


def safe_name(text: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text))).strip("_")[:170]


def fmt(value, digits: int = 2) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def fmt_p(value) -> str:
    if value is None or pd.isna(value):
        return "NA"
    value = float(value)
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def sign(value) -> int:
    if value is None or pd.isna(value):
        return 0
    value = float(value)
    return 1 if value > 0 else -1 if value < 0 else 0


def scaled_color(value, max_abs: float = 2.0) -> str:
    if value is None or pd.isna(value):
        return "#FFFFFF"
    value = max(-max_abs, min(max_abs, float(value)))
    if value < 0:
        alpha = abs(value) / max_abs
        return blend("#FFFFFF", COLORS["ref_down"], alpha)
    if value > 0:
        alpha = abs(value) / max_abs
        return blend("#FFFFFF", COLORS["ref_up"], alpha)
    return "#F3F4F6"


TIER_PRIORITY = {"strict": 0, "suggestive": 1, "exploratory": 2, "not_tiered": 3, "not_detected": 4}

ONTOLOGY = {
    "alzheimer_disease": ("Nervous system", "neurodegenerative disease", "Alzheimer disease; dementia-related neurodegeneration", 1, 1),
    "parkinson_disease": ("Nervous system", "neurodegenerative disease", "Parkinson disease; movement-disorder neurodegeneration", 1, 2),
    "multiple_sclerosis": ("Nervous system", "demyelinating autoimmune disease", "Multiple sclerosis; neuroimmune demyelination", 1, 3),
    "metabolic_dysfunction": ("Metabolic/endocrine", "metabolic disease / diabetes mellitus", "Metabolic dysfunction including diabetes mellitus", 2, 1),
    "primary_serum_diabetes": ("Metabolic/endocrine", "diabetes mellitus", "Serum diabetes contrasts in transplant recipients", 2, 2),
    "conditional_urine_dkd_prognosis": ("Kidney/urinary", "diabetic kidney disease / renal prognosis", "Urine T2D renal dysfunction/prognosis contrast", 3, 1),
    "kidney_dysfunction": ("Kidney/urinary", "renal dysfunction", "Kidney dysfunction case-control proteomics", 3, 2),
    "chronic_kidney_disease": ("Kidney/urinary", "chronic kidney disease", "CKD all-stage and stage-specific urine contrasts", 3, 3),
    "conditional_plasma_aki_subphenotype": ("Kidney/urinary", "acute kidney injury", "AKI plasma subphenotype 2 versus subphenotype 1", 3, 4),
    "tissue_context_proximal_tubule": ("Kidney/urinary", "proximal tubule tissue context", "Human FFPE proximal tubule diabetes-related tissue context", 3, 5),
    "liver_dysfunction": ("Hepatic", "liver disease / hepatic dysfunction", "Liver dysfunction case-control proteomics", 4, 1),
    "cardiac_disease": ("Cardiovascular/inflammatory", "cardiovascular disease / inflammatory cardiac involvement", "Cardiac disease and MIS-C cardiac involvement", 5, 1),
}


def prepare_marker_reference() -> pd.DataFrame:
    try:
        marker = pd.read_excel(MARKER_XLSX)
    except PermissionError:
        marker = pd.DataFrame(
            [
                ("ACTG1", 25.36, 19.38, -5.98),
                ("VWF", 23.70, 27.98, 4.28),
                ("DNAH9", 21.56, 27.40, 5.84),
                ("GPX3", 25.43, 22.60, -2.83),
                ("C4B", 25.71, 21.56, -4.15),
                ("CD44", 24.76, 26.34, 1.58),
                ("CFHR2", 25.30, 26.77, 1.47),
                ("ITIH3", 27.18, 25.97, -1.21),
                ("LRG1", 27.22, 25.50, -1.72),
                ("MYH7B", 26.26, 20.94, -5.32),
            ],
            columns=["Marker", "log2(abundance control)", "Log2(abundance chronic)", "log2FC (Chronic/control)"],
        )
    marker = marker.rename(
        columns={
            "Marker": "marker",
            "log2(abundance control)": "reference_log2_control",
            "Log2(abundance chronic)": "reference_log2_chronic",
            "log2FC (Chronic/control)": "reference_log2FC",
        }
    )
    marker["marker"] = marker["marker"].astype(str).str.strip().str.upper()
    for col in ["reference_log2_control", "reference_log2_chronic", "reference_log2FC"]:
        marker[col] = pd.to_numeric(marker[col], errors="coerce")
    marker["reference_direction"] = marker["reference_log2FC"].apply(
        lambda v: "up_in_chronic" if sign(v) > 0 else "down_in_chronic" if sign(v) < 0 else "flat"
    )
    return marker


def add_ontology(summary: pd.DataFrame) -> pd.DataFrame:
    summary = summary.copy()
    rows = []
    for _, row in summary.iterrows():
        domain, parent, detail, domain_rank, context_rank = ONTOLOGY.get(
            str(row["context"]),
            ("Other", "human disease context", "Human disease contrast", 9, 99),
        )
        r = row.to_dict()
        r.update(
            {
                "ontology_domain": domain,
                "ontology_parent": parent,
                "ontology_detail": detail,
                "domain_rank": domain_rank,
                "context_rank": context_rank,
                "human_disease_scope": "human",
            }
        )
        rows.append(r)
    out = pd.DataFrame(rows)
    return out.sort_values(["domain_rank", "context_rank", "disease", "accession", "case_label"], kind="mergesort").reset_index(drop=True)


def token_set(values) -> set[str]:
    tokens = set()
    for value in values:
        if value is None or pd.isna(value):
            continue
        for token in re.split(r"[;,\s|/]+", str(value).upper()):
            token = token.strip()
            if token:
                tokens.add(token)
    return tokens


def tier_for_row(row: pd.Series) -> str:
    tier = str(row.get("tier", "")).strip().lower()
    if tier in {"strict", "suggestive", "exploratory"}:
        return tier
    for col in ["strict", "suggestive", "exploratory"]:
        if bool(row.get(col, False)):
            return col
    return "not_tiered"


def build_concordance(de: pd.DataFrame, summary: pd.DataFrame, markers: pd.DataFrame):
    de = de.copy()
    for col in ["log2FC", "p_value", "FDR", "abs_log2FC"]:
        if col in de.columns:
            de[col] = pd.to_numeric(de[col], errors="coerce")
    de["_tokens"] = de.apply(lambda r: token_set([r.get("gene_symbol"), r.get("display_gene"), r.get("feature_id")]), axis=1)
    de["_tier_clean"] = de.apply(tier_for_row, axis=1)
    de["_tier_priority"] = de["_tier_clean"].map(TIER_PRIORITY).fillna(9).astype(int)
    de["_sort_fdr"] = de["FDR"].fillna(float("inf"))
    de["_sort_p"] = de["p_value"].fillna(float("inf"))
    de["_sort_abs"] = de["abs_log2FC"].fillna(de["log2FC"].abs()).fillna(0)
    ref = markers.set_index("marker")["reference_log2FC"].to_dict()
    records = []
    for _, meta in summary.iterrows():
        subset = de[de["contrast_id"] == meta["contrast_id"]]
        for marker in markers["marker"]:
            base = {
                "marker": marker,
                "reference_log2FC": ref[marker],
                "contrast_id": meta["contrast_id"],
                "context": meta["context"],
                "disease": meta["disease"],
                "accession": meta["accession"],
                "case_label": meta["case_label"],
                "control_label": meta["control_label"],
                "n_case": meta["n_case"],
                "n_control": meta["n_control"],
                "ontology_domain": meta["ontology_domain"],
                "ontology_parent": meta["ontology_parent"],
                "ontology_detail": meta["ontology_detail"],
            }
            hits = subset[subset["_tokens"].apply(lambda tokens: marker in tokens)]
            if hits.empty:
                records.append({**base, "detected": False, "log2FC": pd.NA, "p_value": pd.NA, "FDR": pd.NA, "tier": "not_detected", "concordance_status": "not_detected_or_unmapped"})
                continue
            best = hits.sort_values(["_tier_priority", "_sort_fdr", "_sort_p", "_sort_abs"], ascending=[True, True, True, False]).iloc[0]
            lfc = best["log2FC"]
            same = sign(lfc) != 0 and sign(ref[marker]) != 0 and sign(lfc) == sign(ref[marker])
            tier = tier_for_row(best)
            if sign(lfc) == 0:
                status = "directionless"
            elif tier in {"strict", "suggestive"}:
                status = "concordant_dep" if same else "discordant_dep"
            elif tier == "exploratory":
                status = "concordant_exploratory" if same else "discordant_exploratory"
            else:
                status = "weak_same_direction" if same else "weak_opposite_direction"
            records.append(
                {
                    **base,
                    "detected": True,
                    "display_gene": best.get("display_gene", ""),
                    "feature_id": best.get("feature_id", ""),
                    "uniprot_id": best.get("uniprot_id", ""),
                    "protein_name": best.get("protein_name", ""),
                    "log2FC": lfc,
                    "p_value": best.get("p_value", pd.NA),
                    "FDR": best.get("FDR", pd.NA),
                    "tier": tier,
                    "concordance_status": status,
                }
            )
    concordance = pd.DataFrame(records)
    disease_counts = concordance.groupby("contrast_id").agg(markers_detected=("detected", "sum")).reset_index()
    marker_summary = concordance.groupby("marker").agg(detected_contrasts=("detected", "sum")).reset_index()
    return concordance, disease_counts, marker_summary


def build_role_matrix(enrichment: pd.DataFrame, markers: list[str]):
    enrichment = enrichment.copy()
    enrichment["fdr"] = pd.to_numeric(enrichment["fdr"], errors="coerce")
    enrichment["role_label"] = enrichment.apply(
        lambda r: f"{r['category']}: {str(r['description'])[:70]} (FDR {fmt_p(r['fdr'])})",
        axis=1,
    )
    records = []
    for _, role in enrichment.iterrows():
        members = {x.strip().upper() for x in str(role.get("preferredNames", "")).split(",") if x.strip()}
        for marker in markers:
            records.append(
                {
                    "marker": marker,
                    "role_label": role["role_label"],
                    "category": role["category"],
                    "term": role["term"],
                    "description": role["description"],
                    "fdr": role["fdr"],
                    "is_member": marker in members,
                }
            )
    return pd.DataFrame(records), enrichment


def load_data():
    markers = prepare_marker_reference()
    edges = pd.read_csv(STRING_EDGES, sep="\t")
    enrichment = pd.read_csv(STRING_ENRICHMENT, sep="\t")
    enrichment["fdr"] = pd.to_numeric(enrichment["fdr"], errors="coerce")

    try:
        de = pd.read_csv(DE_TABLE, sep="\t")
        summary = pd.read_csv(SUMMARY_TABLE, sep="\t")
        summary = add_ontology(summary)
        for col in [
            "n_case",
            "n_control",
            "n_features_tested",
            "strict_count",
            "suggestive_count",
            "exploratory_count",
            "tiered_count",
        ]:
            summary[col] = pd.to_numeric(summary[col], errors="coerce").fillna(0).astype(int)
        concordance, disease_counts, marker_summary = build_concordance(de, summary, markers)
        role_matrix, enrichment_labeled = build_role_matrix(enrichment, list(markers["marker"]))
        role_defs = enrichment_labeled.copy()
        role_defs["fdr"] = pd.to_numeric(role_defs["fdr"], errors="coerce")
        role_defs = role_defs.sort_values(["category", "fdr"], kind="mergesort").reset_index(drop=True)
        role_defs["role_code"] = [f"R{i + 1}" for i in range(len(role_defs))]
        role_matrix = role_matrix.merge(role_defs[["role_label", "role_code"]], on="role_label", how="left")
        summary = add_cluster_order(summary, concordance, list(markers["marker"]))
    except PermissionError:
        de = pd.DataFrame()
        summary = pd.read_csv(TABLES / "contrast_ontology_cluster_metadata_v2.tsv", sep="\t")
        concordance = pd.read_csv(TABLES / "marker10_concordance_by_contrast_v2.tsv", sep="\t")
        role_defs = pd.read_csv(TABLES / "panel_E_pathway_database_role_definitions.tsv", sep="\t")
        role_matrix = pd.read_csv(TABLES / "panel_E_marker10_pathway_role_membership.tsv", sep="\t")
        concordance["detected"] = concordance["detected"].astype(str).str.lower().isin(["true", "1", "yes"])
        role_matrix["is_member"] = role_matrix["is_member"].astype(str).str.lower().isin(["true", "1", "yes"])
        disease_counts = concordance.groupby("contrast_id").agg(markers_detected=("detected", "sum")).reset_index()
        marker_summary = concordance.groupby("marker").agg(detected_contrasts=("detected", "sum")).reset_index()
    return markers, de, summary, concordance, disease_counts, marker_summary, edges, enrichment, role_matrix, role_defs


def contrast_vector(contrast_id: str, concordance: pd.DataFrame, marker_order: list[str]) -> list[float]:
    sub = concordance[concordance["contrast_id"] == contrast_id].set_index("marker")
    vals = []
    for marker in marker_order:
        if marker not in sub.index or pd.isna(sub.loc[marker, "log2FC"]):
            vals.append(0.0)
        else:
            vals.append(max(-2.0, min(2.0, float(sub.loc[marker, "log2FC"]))))
    return vals


def euclid(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def greedy_cluster(ids: list[str], vectors: dict[str, list[float]]) -> list[str]:
    if not ids:
        return []
    remaining = set(ids)
    start = max(ids, key=lambda cid: sum(abs(v) for v in vectors[cid]))
    ordered = [start]
    remaining.remove(start)
    while remaining:
        last = ordered[-1]
        nxt = min(remaining, key=lambda cid: (euclid(vectors[last], vectors[cid]), cid))
        ordered.append(nxt)
        remaining.remove(nxt)
    return ordered


def add_cluster_order(summary: pd.DataFrame, concordance: pd.DataFrame, marker_order: list[str]) -> pd.DataFrame:
    summary = summary.copy()
    vectors = {cid: contrast_vector(cid, concordance, marker_order) for cid in summary["contrast_id"]}
    ordered: list[str] = []
    for domain in DOMAIN_ORDER:
        domain_ids = summary.loc[summary["ontology_domain"] == domain, "contrast_id"].tolist()
        ordered.extend(greedy_cluster(domain_ids, vectors))
    leftover = [cid for cid in summary["contrast_id"] if cid not in ordered]
    ordered.extend(leftover)
    order_map = {cid: i + 1 for i, cid in enumerate(ordered)}
    summary["ontology_cluster_order"] = summary["contrast_id"].map(order_map)
    summary["cohort_number"] = [f"C{i:02d}" for i in summary["ontology_cluster_order"]]
    return summary.sort_values("ontology_cluster_order").reset_index(drop=True)


def build_signature_tables(markers, summary, concordance, role_matrix, role_defs):
    marker_order = list(markers["marker"])
    lfc = concordance.pivot(index="contrast_id", columns="marker", values="log2FC").reindex(summary["contrast_id"])[marker_order]
    status = concordance.pivot(index="contrast_id", columns="marker", values="concordance_status").reindex(summary["contrast_id"])[marker_order]
    detected = concordance.pivot(index="contrast_id", columns="marker", values="detected").reindex(summary["contrast_id"])[marker_order]
    lfc.insert(0, "cohort_number", summary["cohort_number"].values)
    lfc.insert(1, "disease", summary["disease"].values)
    lfc.insert(2, "accession", summary["accession"].values)
    lfc.insert(3, "case_control_n", [f"{r.n_case}/{r.n_control}" for r in summary.itertuples()])
    status.insert(0, "cohort_number", summary["cohort_number"].values)
    status.insert(1, "disease", summary["disease"].values)
    status.insert(2, "accession", summary["accession"].values)
    lfc.to_csv(TABLES / "panel_A_marker10_signature_matrix_log2fc.tsv", sep="\t")
    status.to_csv(TABLES / "panel_B_marker10_concordance_status_matrix.tsv", sep="\t")
    summary.to_csv(TABLES / "contrast_ontology_cluster_metadata_v2.tsv", sep="\t", index=False)
    concordance.to_csv(TABLES / "marker10_concordance_by_contrast_v2.tsv", sep="\t", index=False)

    # Disease pathway proximity score: for each cohort and STRING/database role,
    # score concordant members as +1, discordant members as -1, weak signs half.
    role_lookup = role_matrix[role_matrix["is_member"]].groupby("role_code")["marker"].apply(list).to_dict()
    conc_lookup = concordance.set_index(["contrast_id", "marker"])
    records = []
    for _, meta in summary.iterrows():
        for _, role in role_defs.iterrows():
            members = role_lookup.get(role["role_code"], [])
            score = 0.0
            detected_members = 0
            dep_members = 0
            for marker in members:
                if (meta["contrast_id"], marker) not in conc_lookup.index:
                    continue
                row = conc_lookup.loc[(meta["contrast_id"], marker)]
                if not bool(row["detected"]):
                    continue
                detected_members += 1
                status = row["concordance_status"]
                if status == "concordant_dep":
                    score += 1.0
                    dep_members += 1
                elif status == "discordant_dep":
                    score -= 1.0
                    dep_members += 1
                elif status == "concordant_exploratory":
                    score += 0.75
                elif status == "discordant_exploratory":
                    score -= 0.75
                elif status == "weak_same_direction":
                    score += 0.35
                elif status == "weak_opposite_direction":
                    score -= 0.35
            records.append(
                {
                    "contrast_id": meta["contrast_id"],
                    "cohort_number": meta["cohort_number"],
                    "disease": meta["disease"],
                    "accession": meta["accession"],
                    "ontology_domain": meta["ontology_domain"],
                    "role_code": role["role_code"],
                    "category": role["category"],
                    "term": role["term"],
                    "description": role["description"],
                    "fdr": role["fdr"],
                    "role_marker_count": len(members),
                    "detected_marker_count": detected_members,
                    "dep_marker_count": dep_members,
                    "proximity_score": score,
                }
            )
    prox = pd.DataFrame(records)
    prox.to_csv(TABLES / "panel_D_disease_pathway_proximity_scores.tsv", sep="\t", index=False)
    role_defs.to_csv(TABLES / "panel_E_pathway_database_role_definitions.tsv", sep="\t", index=False)
    role_matrix.to_csv(TABLES / "panel_E_marker10_pathway_role_membership.tsv", sep="\t", index=False)
    return lfc, status, detected, prox


def draw_signature_cluster(path: Path, markers, summary, concordance, title: str, domain_filter: str | None = None, letter: str = "A"):
    marker_order = list(markers["marker"])
    rows = summary if domain_filter is None else summary[summary["ontology_domain"] == domain_filter].copy()
    rows = rows.sort_values("ontology_cluster_order").reset_index(drop=True)
    n = len(rows)
    width = 2600
    height = max(900, 250 + n * 44)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        48,
        width - 116,
        title,
        "Rows are human disease cohorts ordered by ontology/domain and clustered within domain by Marker10 log2FC similarity. This panel is cohort-centered, not CKD-centered.",
        letter,
    )
    left_label = 520
    heat_x = left_label
    heat_y = yy + 70
    cell_w = 88
    cell_h = 36
    max_abs = 2.0
    lookup = concordance.set_index(["contrast_id", "marker"])
    for j, marker in enumerate(marker_order):
        x = heat_x + j * cell_w
        draw.text((x + cell_w // 2, heat_y - 34), marker, font=FONTS["small"], fill=COLORS["ink"], anchor="ma")
    prev_domain = None
    for i, row in rows.iterrows():
        y = heat_y + i * cell_h
        domain = row["ontology_domain"]
        dcol = DOMAIN_COLORS.get(domain, "#8C8C8C")
        if domain != prev_domain:
            draw.line([40, y - 4, width - 50, y - 4], fill=dcol, width=3)
            prev_domain = domain
        draw.rectangle([58, y + 5, 68, y + cell_h - 5], fill=dcol)
        label = f"{row['cohort_number']} | {row['accession']} | {row['disease']} | {row['case_label']} vs {row['control_label']} | n={row['n_case']}/{row['n_control']}"
        draw.text((78, y + 8), label[:78], font=FONTS["tiny"], fill=COLORS["ink"])
        for j, marker in enumerate(marker_order):
            x = heat_x + j * cell_w
            r = lookup.loc[(row["contrast_id"], marker)]
            value = r["log2FC"]
            fill = scaled_color(value, max_abs=max_abs)
            draw.rectangle([x, y, x + cell_w - 3, y + cell_h - 3], fill=fill, outline="#CBD5E1")
            if not bool(r["detected"]):
                draw.text((x + cell_w // 2, y + cell_h // 2 - 1), ".", font=FONTS["tiny"], fill="#94A3B8", anchor="mm")
            else:
                tier = str(r["tier"])
                if tier == "strict":
                    glyph = "S"
                elif tier == "suggestive":
                    glyph = "G"
                elif tier == "exploratory":
                    glyph = "E"
                else:
                    glyph = ""
                if glyph:
                    draw.text((x + cell_w // 2, y + cell_h // 2 - 1), glyph, font=FONTS["tiny"], fill=COLORS["ink"], anchor="mm")
    legend_y = height - 72
    draw.rectangle([heat_x, legend_y, heat_x + 40, legend_y + 18], fill=scaled_color(-2))
    draw.text((heat_x + 48, legend_y), "negative log2FC", font=FONTS["tiny"], fill=COLORS["muted"])
    draw.rectangle([heat_x + 210, legend_y, heat_x + 250, legend_y + 18], fill=scaled_color(2))
    draw.text((heat_x + 258, legend_y), "positive log2FC", font=FONTS["tiny"], fill=COLORS["muted"])
    draw.text((heat_x + 460, legend_y), "S/G/E = strict/suggestive/exploratory marker tier; . = not detected/unmapped", font=FONTS["tiny"], fill=COLORS["muted"])
    save_with_pdf(img, path)
    return img


def draw_signature_multibar(path: Path, markers, summary, concordance, title: str, domain_filter: str | None = None, letter: str = "A") -> Image.Image:
    marker_order = list(markers["marker"])
    rows = summary if domain_filter is None else summary[summary["ontology_domain"] == domain_filter].copy()
    rows = rows.sort_values("ontology_cluster_order").reset_index(drop=True)
    n = len(rows)
    width = 3200
    height = max(980, 260 + n * 48)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        48,
        width - 116,
        title,
        "Multi-bar view: each disease/cohort row contains one mini bar per Marker10 protein. Bar direction and length show cohort log2FC; S/G/E marks the DE tier.",
        letter,
    )
    label_w = 690
    plot_x = label_w
    plot_y = yy + 78
    cell_w = 116
    row_h = 40
    max_abs = 2.0
    lookup = concordance.set_index(["contrast_id", "marker"])
    for j, marker in enumerate(marker_order):
        x = plot_x + j * cell_w
        draw.text((x + cell_w // 2, plot_y - 38), marker, font=FONTS["small"], fill=COLORS["ink"], anchor="ma")
        draw.line([x + cell_w // 2, plot_y - 8, x + cell_w // 2, plot_y + n * row_h], fill="#E2E8F0", width=1)

    prev_domain = None
    for i, row in rows.iterrows():
        y = plot_y + i * row_h
        domain = row["ontology_domain"]
        dcol = DOMAIN_COLORS.get(domain, "#8C8C8C")
        if domain != prev_domain:
            draw.line([40, y - 4, width - 55, y - 4], fill=dcol, width=3)
            prev_domain = domain
        draw.rectangle([58, y + 6, 70, y + row_h - 6], fill=dcol)
        label = f"{row['cohort_number']} | {row['accession']} | {row['disease']} | {row['case_label']} vs {row['control_label']} | n={row['n_case']}/{row['n_control']}"
        draw.text((80, y + 9), label[:90], font=FONTS["tiny"], fill=COLORS["ink"])
        for j, marker in enumerate(marker_order):
            x = plot_x + j * cell_w
            center = x + cell_w // 2
            r = lookup.loc[(row["contrast_id"], marker)]
            draw.rounded_rectangle([x + 8, y + 7, x + cell_w - 8, y + row_h - 7], radius=3, fill="#F8FAFC", outline="#E2E8F0")
            draw.line([center, y + 8, center, y + row_h - 8], fill="#94A3B8", width=1)
            if not bool(r["detected"]) or pd.isna(r["log2FC"]):
                draw.text((center, y + row_h // 2 - 1), ".", font=FONTS["tiny"], fill="#94A3B8", anchor="mm")
                continue
            value = max(-max_abs, min(max_abs, float(r["log2FC"])))
            half = (cell_w - 24) / 2
            end = center + int((value / max_abs) * half)
            color = COLORS["ref_up"] if value > 0 else COLORS["ref_down"]
            draw.rounded_rectangle(
                [min(center, end), y + 12, max(center, end), y + row_h - 12],
                radius=3,
                fill=color,
            )
            tier = str(r["tier"])
            glyph = "S" if tier == "strict" else "G" if tier == "suggestive" else "E" if tier == "exploratory" else ""
            if glyph:
                draw.text((center, y + row_h // 2 - 1), glyph, font=FONTS["tiny"], fill="white" if abs(value) > 0.6 else COLORS["ink"], anchor="mm")

    legend_y = height - 72
    draw.rounded_rectangle([plot_x, legend_y, plot_x + 60, legend_y + 16], radius=3, fill=COLORS["ref_down"])
    draw.text((plot_x + 70, legend_y - 2), "negative log2FC", font=FONTS["tiny"], fill=COLORS["muted"])
    draw.rounded_rectangle([plot_x + 250, legend_y, plot_x + 310, legend_y + 16], radius=3, fill=COLORS["ref_up"])
    draw.text((plot_x + 320, legend_y - 2), "positive log2FC", font=FONTS["tiny"], fill=COLORS["muted"])
    draw.text((plot_x + 520, legend_y - 2), "S/G/E = strict/suggestive/exploratory marker tier; . = not detected/unmapped", font=FONTS["tiny"], fill=COLORS["muted"])
    save_with_pdf(img, path)
    return img


def draw_concordance_panel(path: Path, markers, summary, concordance, title: str, compact: bool = False, letter: str = "B") -> Image.Image:
    marker_order = list(markers["marker"])
    rows = summary.sort_values("ontology_cluster_order").reset_index(drop=True)
    width = 3200 if not compact else 1800
    cell_w = 72 if not compact else 42
    cell_h = 38 if not compact else 28
    label_w = 960 if not compact else 330
    height = max(900, 250 + len(rows) * cell_h)
    if compact:
        height = 1050
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        48,
        width - 116,
        title,
        "Each row includes disease name, accession/cohort number, case/control labels, and n case/control. Concordance is against the supplied Marker10 reference direction.",
        letter,
    )
    x0 = label_w
    y0 = yy + 62
    lookup = concordance.set_index(["contrast_id", "marker"])
    if not compact:
        draw.text((58, y0 - 36), "cohort | disease | accession | contrast | n", font=FONTS["small"], fill=COLORS["muted"])
    for j, marker in enumerate(marker_order):
        x = x0 + j * cell_w
        draw.text((x + cell_w // 2, y0 - 34), marker, font=FONTS["small"] if not compact else FONTS["tiny"], fill=COLORS["ink"], anchor="ma")
    for i, row in rows.iterrows():
        y = y0 + i * cell_h
        dcol = DOMAIN_COLORS.get(row["ontology_domain"], "#8C8C8C")
        draw.rectangle([58, y + 5, 70, y + cell_h - 5], fill=dcol)
        if compact:
            label = f"{row['cohort_number']} {row['accession']}"
            font = FONTS["tiny"]
        else:
            label = f"{row['cohort_number']} | {row['disease']} | {row['accession']} | {row['case_label']} vs {row['control_label']} | n={row['n_case']}/{row['n_control']}"
            font = FONTS["tiny"]
        draw.text((78, y + 8), label[:120], font=font, fill=COLORS["ink"])
        for j, marker in enumerate(marker_order):
            x = x0 + j * cell_w
            r = lookup.loc[(row["contrast_id"], marker)]
            status = r["concordance_status"]
            fill = code_color(status)
            draw.rectangle([x, y, x + cell_w - 3, y + cell_h - 3], fill=fill, outline="#CBD5E1")
            code = STATUS_CODES.get(status, "")
            draw.text((x + cell_w // 2, y + cell_h // 2 - 1), code, font=FONTS["tiny"], fill=COLORS["ink"], anchor="mm")
    legend_y = height - 82
    lx = 58
    for status in STATUS_ORDER[:6]:
        draw.rectangle([lx, legend_y, lx + 18, legend_y + 14], fill=code_color(status), outline="#94A3B8")
        draw.text((lx + 24, legend_y - 1), STATUS_LABELS[status], font=FONTS["tiny"], fill=COLORS["muted"])
        lx += 230
        if lx > width - 360:
            lx = 58
            legend_y += 24
    save_with_pdf(img, path)
    return img


def role_zone_members(role_defs: pd.DataFrame, role_matrix: pd.DataFrame) -> dict[str, set[str]]:
    zones = {
        "Complement/coagulation": set(),
        "Extracellular/vesicle": set(),
        "Tissue/literature roles": set(),
    }
    member = role_matrix[role_matrix["is_member"]].merge(role_defs[["role_label", "category", "description"]], on="role_label", how="left")
    for _, row in member.iterrows():
        desc = str(row["description"]).lower()
        cat = str(row["category"])
        marker = row["marker"]
        if "complement" in desc or "coagulation" in desc or "opsonization" in desc:
            zones["Complement/coagulation"].add(marker)
        if "extracellular" in desc or "vesicle" in desc or "exosome" in desc:
            zones["Extracellular/vesicle"].add(marker)
        if cat in {"TISSUES", "PMID"}:
            zones["Tissue/literature roles"].add(marker)
    return zones


def draw_marker_network_consensus(path: Path, markers, summary, concordance, edges, role_defs, role_matrix, title: str, letter: str = "C") -> Image.Image:
    marker_order = list(markers["marker"])
    img = Image.new("RGB", (2400, 1550), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        48,
        2284,
        title,
        "Fixed Marker10 layout using STRING direct edges and enriched role zones. Node size reflects detected cohorts; node color reflects aggregate concordance across human disease cohorts.",
        letter,
    )
    px0, py0, px1, py1 = 90, yy + 70, 1600, 1260
    draw.rounded_rectangle([px0, py0, px1, py1], radius=12, fill="#FBFCFD", outline="#CBD5E1", width=2)
    zones = [
        ("Complement/coagulation", [px0 + 710, py0 + 70, px1 - 85, py0 + 410], "#FDF2E9"),
        ("Extracellular/vesicle", [px0 + 90, py1 - 430, px0 + 1040, py1 - 80], "#ECFDF3"),
        ("Tissue/literature roles", [px0 + 90, py0 + 90, px0 + 590, py0 + 420], "#FFF7ED"),
    ]
    for label, rect, fill in zones:
        draw.rounded_rectangle(rect, radius=16, fill=fill, outline="#CBD5E1")
        draw.text((rect[0] + 16, rect[1] + 14), label, font=FONTS["small"], fill=COLORS["muted"])

    def xy(marker):
        rx, ry = MARKER_POSITIONS.get(marker, (0.5, 0.5))
        return int(px0 + rx * (px1 - px0)), int(py0 + ry * (py1 - py0))

    for _, edge in edges.iterrows():
        a = str(edge["preferredName_A"]).upper()
        b = str(edge["preferredName_B"]).upper()
        if a not in MARKER_POSITIONS or b not in MARKER_POSITIONS:
            continue
        ax, ay = xy(a)
        bx, by = xy(b)
        score = float(edge.get("score", 0.4))
        draw.line([ax, ay, bx, by], fill="#64748B", width=max(2, int(score * 8)))
        draw.text(((ax + bx) // 2 + 8, (ay + by) // 2 - 12), f"{score:.2f}", font=FONTS["tiny"], fill=COLORS["muted"])

    agg_records = []
    for marker in marker_order:
        sub = concordance[concordance["marker"] == marker]
        detected = int(sub["detected"].sum())
        c_dep = int((sub["concordance_status"] == "concordant_dep").sum())
        d_dep = int((sub["concordance_status"] == "discordant_dep").sum())
        c_exp = int((sub["concordance_status"] == "concordant_exploratory").sum())
        d_exp = int((sub["concordance_status"] == "discordant_exploratory").sum())
        w_same = int((sub["concordance_status"] == "weak_same_direction").sum())
        w_opp = int((sub["concordance_status"] == "weak_opposite_direction").sum())
        score = (c_dep - d_dep) + 0.75 * (c_exp - d_exp) + 0.25 * (w_same - w_opp)
        agg_records.append(
            {
                "marker": marker,
                "detected": detected,
                "concordant_dep": c_dep,
                "discordant_dep": d_dep,
                "weak_same": w_same,
                "weak_opposite": w_opp,
                "aggregate_score": score,
            }
        )
    agg = pd.DataFrame(agg_records).set_index("marker")
    max_detected = max(1, int(agg["detected"].max()))

    for marker in marker_order:
        mx, my = xy(marker)
        row = agg.loc[marker]
        score = float(row["aggregate_score"])
        if score > 0:
            fill = blend("#FFFFFF", COLORS["concordant_dep"], min(1.0, score / 6.0))
        elif score < 0:
            fill = blend("#FFFFFF", COLORS["discordant_dep"], min(1.0, abs(score) / 6.0))
        else:
            fill = "#F8FAFC"
        radius = 25 + int(26 * float(row["detected"]) / max_detected)
        draw.ellipse([mx - radius, my - radius, mx + radius, my + radius], fill=fill, outline="#111827", width=3)
        draw.text((mx, my - 5), str(int(row["detected"])), font=FONTS["subtitle"], fill=COLORS["ink"], anchor="mm")
        draw.text((mx, my + 15), "det.", font=FONTS["tiny"], fill=COLORS["muted"], anchor="mm")
        tw, th = text_size(draw, marker, FONTS["small"])
        draw.rounded_rectangle([mx - tw // 2 - 7, my + radius + 8, mx + tw // 2 + 7, my + radius + th + 12], radius=4, fill="white", outline="#CBD5E1")
        draw.text((mx - tw // 2, my + radius + 10), marker, font=FONTS["small"], fill=COLORS["ink"])

    side_x = 1665
    draw.rounded_rectangle([side_x, py0, 2315, py1], radius=10, fill="#F8FAFC", outline="#CBD5E1")
    draw.text((side_x + 24, py0 + 26), "Aggregate marker context", font=FONTS["subtitle"], fill=COLORS["ink"])
    draw.text((side_x + 24, py0 + 64), "marker | detected | Cdep | Ddep | weak +/-", font=FONTS["tiny"], fill=COLORS["muted"])
    y = py0 + 95
    for marker in marker_order:
        row = agg.loc[marker]
        score = float(row["aggregate_score"])
        color = COLORS["concordant_dep"] if score > 0 else COLORS["discordant_dep"] if score < 0 else "#CBD5E1"
        draw.rectangle([side_x + 24, y + 4, side_x + 42, y + 22], fill=color, outline="#94A3B8")
        label = f"{marker:5s} | {int(row['detected']):2d} | {int(row['concordant_dep']):2d} | {int(row['discordant_dep']):2d} | {int(row['weak_same'])}/{int(row['weak_opposite'])}"
        draw.text((side_x + 52, y + 1), label, font=FONTS["small"], fill=COLORS["ink"])
        y += 48
    wrapped(
        draw,
        (side_x + 24, py1 - 170),
        "Panel C is the global Marker10 network/proximity context. The per-cohort versions are saved separately in the Panel C disease map folder.",
        560,
        FONTS["small"],
        fill=COLORS["muted"],
        max_lines=5,
    )
    save_with_pdf(img, path)
    agg.reset_index().to_csv(TABLES / "panel_C_marker10_network_consensus_summary.tsv", sep="\t", index=False)
    return img


def draw_disease_proximity_map(path: Path, meta: pd.Series, markers, concordance, edges, role_defs, role_matrix) -> dict:
    marker_order = list(markers["marker"])
    sub = concordance[concordance["contrast_id"] == meta["contrast_id"]].set_index("marker")
    img = Image.new("RGB", (1800, 1280), "white")
    draw = ImageDraw.Draw(img)
    domain_color = DOMAIN_COLORS.get(meta["ontology_domain"], "#8C8C8C")
    draw.rectangle([0, 0, 1800, 18], fill=domain_color)
    title = f"{meta['cohort_number']} | {meta['accession']} | {meta['disease']}"
    draw.text((56, 44), title, font=FONTS["title"], fill=COLORS["ink"])
    subtitle = f"{meta['ontology_domain']} > {meta['ontology_parent']} | {meta['case_label']} vs {meta['control_label']} | n={meta['n_case']}/{meta['n_control']}"
    wrapped(draw, (58, 90), subtitle, 1200, FONTS["body"], fill=COLORS["muted"], max_lines=2)
    px0, py0, px1, py1 = 80, 190, 1260, 1060
    draw.rounded_rectangle([px0, py0, px1, py1], radius=10, fill="#FBFCFD", outline="#CBD5E1", width=2)
    zones = [
        ("Complement/coagulation", [px0 + 565, py0 + 50, px1 - 70, py0 + 310], "#FDF2E9"),
        ("Extracellular/vesicle", [px0 + 70, py1 - 335, px0 + 805, py1 - 65], "#ECFDF3"),
        ("Tissue/literature roles", [px0 + 70, py0 + 70, px0 + 460, py0 + 330], "#FFF7ED"),
    ]
    for label, rect, fill in zones:
        draw.rounded_rectangle(rect, radius=14, fill=fill, outline="#CBD5E1")
        draw.text((rect[0] + 14, rect[1] + 12), label, font=FONTS["tiny"], fill=COLORS["muted"])

    def xy(marker):
        rx, ry = MARKER_POSITIONS.get(marker, (0.5, 0.5))
        return int(px0 + rx * (px1 - px0)), int(py0 + ry * (py1 - py0))

    for _, edge in edges.iterrows():
        a = str(edge["preferredName_A"]).upper()
        b = str(edge["preferredName_B"]).upper()
        if a not in MARKER_POSITIONS or b not in MARKER_POSITIONS:
            continue
        ax, ay = xy(a)
        bx, by = xy(b)
        score = float(edge.get("score", 0.4))
        draw.line([ax, ay, bx, by], fill="#64748B", width=max(2, int(score * 7)))

    for marker in marker_order:
        mx, my = xy(marker)
        row = sub.loc[marker]
        status = row["concordance_status"]
        fill = code_color(status)
        detected = bool(row["detected"])
        tier = str(row["tier"])
        radius = 28
        if tier == "strict":
            radius = 36
        elif tier == "suggestive":
            radius = 32
        elif tier == "exploratory":
            radius = 30
        outline = "#111827" if detected else "#94A3B8"
        draw.ellipse([mx - radius, my - radius, mx + radius, my + radius], fill=fill, outline=outline, width=3)
        code = STATUS_CODES.get(status, ".")
        draw.text((mx, my - 2), code, font=FONTS["small"], fill=COLORS["ink"], anchor="mm")
        tw, th = text_size(draw, marker, FONTS["small"])
        draw.rounded_rectangle([mx - tw // 2 - 6, my + radius + 6, mx + tw // 2 + 6, my + radius + th + 10], radius=4, fill="white", outline="#CBD5E1")
        draw.text((mx - tw // 2, my + radius + 8), marker, font=FONTS["small"], fill=COLORS["ink"])
    side_x = 1320
    draw.rounded_rectangle([side_x, 190, 1738, 1060], radius=10, fill="#F8FAFC", outline="#CBD5E1")
    draw.text((side_x + 24, 216), "Marker10 status", font=FONTS["subtitle"], fill=COLORS["ink"])
    y = 262
    for marker in marker_order:
        row = sub.loc[marker]
        draw.rectangle([side_x + 24, y + 4, side_x + 44, y + 24], fill=code_color(row["concordance_status"]), outline="#94A3B8")
        lfc = "NA" if pd.isna(row["log2FC"]) else fmt(row["log2FC"])
        label = f"{marker} | {lfc} | {row['tier']}"
        draw.text((side_x + 54, y + 2), label[:36], font=FONTS["small"], fill=COLORS["ink"])
        y += 48
    wrapped(
        draw,
        (side_x + 24, 930),
        "Fixed marker layout uses STRING direct edges and enriched role zones. Node color is the cohort-specific Marker10 concordance status.",
        360,
        FONTS["small"],
        fill=COLORS["muted"],
        max_lines=5,
    )
    save_with_pdf(img, path)
    return {
        "contrast_id": meta["contrast_id"],
        "cohort_number": meta["cohort_number"],
        "disease": meta["disease"],
        "accession": meta["accession"],
        "ontology_domain": meta["ontology_domain"],
        "figure_png": str(path.relative_to(ROOT)),
        "figure_pdf": str(path.with_suffix(".pdf").relative_to(ROOT)),
    }


def draw_pathway_proximity_overlap(path: Path, summary, prox, title: str, letter: str = "D") -> Image.Image:
    roles = prox.groupby(["role_code", "category", "description"], dropna=False)["detected_marker_count"].sum().reset_index()
    roles = roles.sort_values("detected_marker_count", ascending=False).head(14)
    role_order = roles["role_code"].tolist()
    rows = summary.sort_values("ontology_cluster_order").reset_index(drop=True)
    lookup = prox.set_index(["contrast_id", "role_code"])
    width = 3000
    height = max(1000, 250 + len(rows) * 36)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        46,
        width - 116,
        title,
        "Rows are disease cohorts; columns are STRING/database roles. Color shows role-level Marker10 concordance proximity score.",
        letter,
    )
    label_w = 780
    x0 = label_w
    y0 = yy + 76
    cell_w = 72
    cell_h = 30
    for j, code in enumerate(role_order):
        x = x0 + j * cell_w
        cat = roles.loc[roles["role_code"] == code, "category"].iloc[0]
        draw.rounded_rectangle([x + 6, y0 - 48, x + cell_w - 7, y0 - 18], radius=4, fill=ROLE_CATEGORY_COLORS.get(cat, "#8C8C8C"))
        draw.text((x + cell_w // 2, y0 - 40), code, font=FONTS["small"], fill="white", anchor="ma")
    for i, row in rows.iterrows():
        y = y0 + i * cell_h
        dcol = DOMAIN_COLORS.get(row["ontology_domain"], "#8C8C8C")
        draw.rectangle([58, y + 5, 70, y + cell_h - 5], fill=dcol)
        label = f"{row['cohort_number']} | {row['accession']} | {row['disease']} | n={row['n_case']}/{row['n_control']}"
        draw.text((80, y + 6), label[:95], font=FONTS["tiny"], fill=COLORS["ink"])
        for j, code in enumerate(role_order):
            x = x0 + j * cell_w
            r = lookup.loc[(row["contrast_id"], code)]
            score = float(r["proximity_score"])
            fill = "#FFFFFF"
            if score > 0:
                fill = blend("#FFFFFF", COLORS["concordant_dep"], min(1.0, score / 3.0))
            elif score < 0:
                fill = blend("#FFFFFF", COLORS["discordant_dep"], min(1.0, abs(score) / 3.0))
            draw.rectangle([x, y, x + cell_w - 3, y + cell_h - 3], fill=fill, outline="#CBD5E1")
            if abs(score) >= 1:
                draw.text((x + cell_w // 2, y + cell_h // 2 - 1), f"{score:.0f}", font=FONTS["tiny"], fill=COLORS["ink"], anchor="mm")
    role_x = x0 + len(role_order) * cell_w + 50
    draw.rounded_rectangle([role_x, y0 - 58, width - 58, min(height - 120, y0 + 560)], radius=8, fill="#F8FAFC", outline="#CBD5E1")
    draw.text((role_x + 20, y0 - 38), "Role key", font=FONTS["subtitle"], fill=COLORS["ink"])
    y = y0 + 5
    for _, row in roles.iterrows():
        code = row["role_code"]
        cat = row["category"]
        draw.rounded_rectangle([role_x + 20, y, role_x + 68, y + 24], radius=4, fill=ROLE_CATEGORY_COLORS.get(cat, "#8C8C8C"))
        draw.text((role_x + 31, y + 3), code, font=FONTS["small"], fill="white")
        wrapped(draw, (role_x + 80, y - 2), f"{cat}: {row['description']}", width - role_x - 160, FONTS["tiny"], fill=COLORS["muted"], max_lines=2)
        y += 48
    save_with_pdf(img, path)
    return img


def draw_enhanced_role_map(path: Path, markers, role_matrix, role_defs, title: str, letter: str = "E") -> Image.Image:
    marker_order = list(markers["marker"])
    roles = role_defs.sort_values(["category", "fdr"], kind="mergesort").reset_index(drop=True)
    lookup = role_matrix.set_index(["marker", "role_code"])
    width, height = 3200, 1900
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yy = panel_header(
        draw,
        58,
        46,
        width - 116,
        title,
        "Expanded readable role matrix with database category, term ID, FDR, and marker membership.",
        letter,
    )
    matrix_x, matrix_y = 180, yy + 160
    cell_w, cell_h = 62, 56
    row_label_w = 95
    for j, role in roles.iterrows():
        x = matrix_x + row_label_w + j * cell_w
        color = ROLE_CATEGORY_COLORS.get(role["category"], "#8C8C8C")
        draw.rounded_rectangle([x + 5, matrix_y - 52, x + cell_w - 6, matrix_y - 16], radius=4, fill=color)
        draw.text((x + cell_w // 2, matrix_y - 42), role["role_code"], font=FONTS["small"], fill="white", anchor="ma")
    for i, marker in enumerate(marker_order):
        y = matrix_y + i * cell_h
        draw.text((matrix_x, y + 16), marker, font=FONTS["body"], fill=COLORS["ink"])
        for j, role in roles.iterrows():
            x = matrix_x + row_label_w + j * cell_w
            member = False
            if (marker, role["role_code"]) in lookup.index:
                member = bool(lookup.loc[(marker, role["role_code"]), "is_member"])
            fill = ROLE_CATEGORY_COLORS.get(role["category"], "#8C8C8C") if member else "#F8FAFC"
            draw.rectangle([x, y, x + cell_w - 4, y + cell_h - 4], fill=fill, outline="#CBD5E1")
            if member:
                draw.text((x + cell_w // 2 - 1, y + cell_h // 2 - 2), "+", font=FONTS["subtitle"], fill="white", anchor="mm")
    key_x = 1750
    draw.rounded_rectangle([key_x - 24, matrix_y - 90, width - 80, height - 165], radius=10, fill="#F8FAFC", outline="#CBD5E1")
    draw.text((key_x, matrix_y - 62), "Role definitions", font=FONTS["subtitle"], fill=COLORS["ink"])
    y = matrix_y - 20
    for _, role in roles.iterrows():
        color = ROLE_CATEGORY_COLORS.get(role["category"], "#8C8C8C")
        draw.rounded_rectangle([key_x, y + 4, key_x + 50, y + 30], radius=4, fill=color)
        draw.text((key_x + 10, y + 6), role["role_code"], font=FONTS["small"], fill="white")
        header = f"{role['category']} | {role['term']} | FDR {fmt_p(role['fdr'])}"
        draw.text((key_x + 64, y + 1), header, font=FONTS["small"], fill=COLORS["ink"])
        y = wrapped(draw, (key_x + 64, y + 25), role["description"], width - key_x - 180, FONTS["tiny"], fill=COLORS["muted"], max_lines=2)
        y += 8
        if y > height - 220:
            break
    legend_y = height - 100
    lx = 180
    for cat, color in ROLE_CATEGORY_COLORS.items():
        draw.rectangle([lx, legend_y, lx + 20, legend_y + 16], fill=color, outline="#64748B")
        draw.text((lx + 28, legend_y - 1), cat, font=FONTS["small"], fill=COLORS["muted"])
        lx += 265
    save_with_pdf(img, path)
    return img


def draw_composite(path: Path, markers, summary, concordance, prox, role_matrix, role_defs, edges) -> Image.Image:
    img = Image.new("RGB", (4600, 4100), "white")
    draw = ImageDraw.Draw(img)
    draw.text((72, 38), "Marker10 v2: disease-centered concordance, network, and pathway-proximity context", font=FONTS["title"], fill=COLORS["ink"])
    draw.text((72, 84), "All panels use human disease cohorts. CKD is only one domain among the ontology-ordered disease contexts.", font=FONTS["body"], fill=COLORS["muted"])

    # Draw compact versions directly into temporary images, then paste.
    tmp_a = FIG / "_tmp_panel_A_compact.png"
    tmp_b = FIG / "_tmp_panel_B_compact.png"
    tmp_c = FIG / "_tmp_panel_C_compact.png"
    tmp_d = FIG / "_tmp_panel_D_compact.png"
    tmp_e = FIG / "_tmp_panel_E_compact.png"
    draw_signature_multibar(tmp_a, markers, summary, concordance, "A. Cohort-level Marker10 multi-bar signatures", letter="A")
    draw_concordance_panel(tmp_b, markers, summary, concordance, "B. Marker concordance by human disease cohort", compact=True, letter="B")
    draw_marker_network_consensus(tmp_c, markers, summary, concordance, edges, role_defs, role_matrix, "C. Marker10 STRING/pathway proximity network", letter="C")
    draw_pathway_proximity_overlap(tmp_d, summary, prox, "D. Disease pathway-proximity overlap", letter="D")
    for p, box in [
        (tmp_a, (70, 130, 2200, 1340)),
        (tmp_b, (2380, 130, 2100, 1340)),
        (tmp_c, (70, 1560, 2120, 1160)),
        (tmp_d, (2300, 1560, 2220, 1160)),
    ]:
        im = crop_content(Image.open(p).convert("RGB"))
        bw, bh = box[2], box[3]
        im.thumbnail((bw, bh), Image.Resampling.LANCZOS)
        img.paste(im, (box[0], box[1]))
    # Panel E compact role matrix.
    draw_enhanced_role_map(tmp_e, markers, role_matrix, role_defs, "E. Marker pathway/database-role matrix", letter="E")
    im = crop_content(Image.open(tmp_e).convert("RGB"))
    im.thumbnail((4200, 1060), Image.Resampling.LANCZOS)
    img.paste(im, (190, 2840))

    draw.text((72, 3990), "Standalone versions include full labels, disease names, accession/cohort numbers, sample sizes, role definitions, per-disease Panel C maps, and Panel D overlap scores.", font=FONTS["small"], fill=COLORS["muted"])
    save_with_pdf(img, path)
    for p in [
        tmp_a,
        tmp_a.with_suffix(".pdf"),
        tmp_b,
        tmp_b.with_suffix(".pdf"),
        tmp_c,
        tmp_c.with_suffix(".pdf"),
        tmp_d,
        tmp_d.with_suffix(".pdf"),
        tmp_e,
        tmp_e.with_suffix(".pdf"),
    ]:
        try:
            p.unlink()
        except FileNotFoundError:
            pass
    return img


def write_report(summary, concordance, prox, role_defs, disease_manifest):
    lines = [
        "# Marker10 Integrated Disease Context v2 Report",
        "",
        "## What Changed From v1",
        "- Panel A is now a cohort/disease Marker10 multi-bar chart, not a CKD reference lollipop.",
        "- Panel B now has a standalone readable concordance heatmap with disease, accession/cohort number, case/control labels, and n.",
        "- Panel C is now the Marker10 STRING/pathway proximity network and also saves one disease-specific map per cohort.",
        "- Panel D is now the cross-disease pathway-proximity overlap heatmap.",
        "- Panel E is expanded with role definitions, database category, term ID, FDR, and marker membership.",
        "",
        "## Scope",
        f"- Human disease contrasts: {len(summary)}",
        f"- Marker-concordance cells: {len(concordance)}",
        f"- Disease-specific pathway-proximity maps: {len(disease_manifest)}",
        f"- STRING/database role terms used: {len(role_defs)}",
        "",
        "## Main Outputs",
        f"- Composite v2: `{(FIG / 'marker10_integrated_disease_context_v2_composite.png').relative_to(ROOT)}`",
        f"- Panel A standalone: `{(FIG / 'panel_A_marker10_multibar_by_cohort.png').relative_to(ROOT)}`",
        f"- Panel B standalone: `{(FIG / 'panel_B_marker10_concordance_by_human_disease_cohort.png').relative_to(ROOT)}`",
        f"- Panel C network: `{(FIG / 'panel_C_marker10_network_proximity_consensus.png').relative_to(ROOT)}`",
        f"- Panel C individual maps: `{PANEL_C_MAPS.relative_to(ROOT)}`",
        f"- Panel D overlap: `{(FIG / 'panel_D_disease_pathway_proximity_overlap.png').relative_to(ROOT)}`",
        f"- Panel E enhanced: `{(FIG / 'panel_E_marker10_pathway_database_roles_enhanced.png').relative_to(ROOT)}`",
        "",
        "## Status Counts",
    ]
    for status, count in concordance["concordance_status"].value_counts().reindex(STATUS_ORDER).fillna(0).astype(int).items():
        lines.append(f"- {STATUS_LABELS.get(status, status)}: {count}")
    (OUT / "marker10_integrated_disease_context_v2_report.md").write_text("\n".join(lines) + "\n")


def checksums(paths: list[Path]) -> None:
    rows = []
    for path in paths:
        if path.exists() and path.is_file():
            rows.append({"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(REPRO / "output_checksums.tsv", sep="\t", index=False)


def main() -> None:
    markers, de, summary, concordance, disease_counts, marker_summary, edges, enrichment, role_matrix, role_defs = load_data()
    lfc, status, detected, prox = build_signature_tables(markers, summary, concordance, role_matrix, role_defs)

    outputs: list[Path] = []
    outputs.append(FIG / "panel_A_marker10_multibar_by_cohort.png")
    draw_signature_multibar(outputs[-1], markers, summary, concordance, "Panel A. Marker10 multi-bar signatures by disease/cohort", letter="A")
    for domain in summary["ontology_domain"].drop_duplicates():
        path = PANEL_A_DOMAINS / f"{safe_name(domain)}_marker10_multibar.png"
        draw_signature_multibar(path, markers, summary, concordance, f"Panel A domain subset: {domain}", domain_filter=domain, letter="A")
        outputs.append(path)

    outputs.append(FIG / "panel_B_marker10_concordance_by_human_disease_cohort.png")
    draw_concordance_panel(outputs[-1], markers, summary, concordance, "Panel B. Marker10 concordance across human disease cohorts", compact=False, letter="B")

    disease_manifest = []
    for _, meta in summary.iterrows():
        path = PANEL_C_MAPS / f"{safe_name(meta['cohort_number'] + '_' + meta['contrast_id'])}_marker_pathway_proximity_map.png"
        disease_manifest.append(draw_disease_proximity_map(path, meta, markers, concordance, edges, role_defs, role_matrix))
        outputs.append(path)
    disease_manifest_df = pd.DataFrame(disease_manifest)
    disease_manifest_df.to_csv(TABLES / "panel_C_disease_marker_proximity_map_manifest.tsv", sep="\t", index=False)

    outputs.append(FIG / "panel_C_marker10_network_proximity_consensus.png")
    draw_marker_network_consensus(outputs[-1], markers, summary, concordance, edges, role_defs, role_matrix, "Panel C. Marker10 STRING/pathway proximity network", letter="C")

    outputs.append(FIG / "panel_D_disease_pathway_proximity_overlap.png")
    draw_pathway_proximity_overlap(outputs[-1], summary, prox, "Panel D. Disease-specific Marker10 pathway-proximity overlap", letter="D")

    outputs.append(FIG / "panel_E_marker10_pathway_database_roles_enhanced.png")
    draw_enhanced_role_map(outputs[-1], markers, role_matrix, role_defs, "Panel E. Enhanced Marker10 pathway and database-role map", letter="E")

    outputs.append(FIG / "marker10_integrated_disease_context_v2_composite.png")
    draw_composite(outputs[-1], markers, summary, concordance, prox, role_matrix, role_defs, edges)

    write_report(summary, concordance, prox, role_defs, disease_manifest_df)
    outputs.extend(
        [
            TABLES / "panel_A_marker10_signature_matrix_log2fc.tsv",
            TABLES / "panel_B_marker10_concordance_status_matrix.tsv",
            TABLES / "contrast_ontology_cluster_metadata_v2.tsv",
            TABLES / "marker10_concordance_by_contrast_v2.tsv",
            TABLES / "panel_C_marker10_network_consensus_summary.tsv",
            TABLES / "panel_D_disease_pathway_proximity_scores.tsv",
            TABLES / "panel_C_disease_marker_proximity_map_manifest.tsv",
            TABLES / "panel_E_pathway_database_role_definitions.tsv",
            TABLES / "panel_E_marker10_pathway_role_membership.tsv",
            OUT / "marker10_integrated_disease_context_v2_report.md",
        ]
    )
    checksums(outputs + [p.with_suffix(".pdf") for p in outputs if p.suffix == ".png"])
    print(f"Composite: {FIG / 'marker10_integrated_disease_context_v2_composite.png'}")
    print(f"Panel A: {FIG / 'panel_A_marker10_multibar_by_cohort.png'}")
    print(f"Panel B: {FIG / 'panel_B_marker10_concordance_by_human_disease_cohort.png'}")
    print(f"Panel C: {FIG / 'panel_C_marker10_network_proximity_consensus.png'}")
    print(f"Panel C maps: {len(disease_manifest_df)} in {PANEL_C_MAPS}")
    print(f"Panel D: {FIG / 'panel_D_disease_pathway_proximity_overlap.png'}")
    print(f"Panel E: {FIG / 'panel_E_marker10_pathway_database_roles_enhanced.png'}")
    print(f"Report: {OUT / 'marker10_integrated_disease_context_v2_report.md'}")


if __name__ == "__main__":
    main()
