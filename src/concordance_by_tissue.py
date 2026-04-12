#!/usr/bin/env python3
"""Validate panel-gene concordance against OP chronic proteomics directions, stratified by tissue."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from statsmodels.stats.proportion import proportion_confint


ROOT = Path(__file__).resolve().parents[1]
TARGETS_YAML = ROOT / "config" / "targets.yaml"
COHORTS_CSV = ROOT / "config" / "cohorts.csv"
DEA_DIR = ROOT / "results" / "geo_dea"
OUT_DIR = ROOT / "results" / "rna" / "concordance"

PANEL_GENES = [
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

SIGNATURE_ORDER = [
    "shared_op_signature",
    "acute_specific_signature",
    "chronic_specific_signature",
]
SIGNATURE_TAG = {
    "shared_op_signature": "shared\u2193",
    "acute_specific_signature": "acute-none",
    "chronic_specific_signature": "chronic\u2193",
}

ACTIVE_ACCESSIONS: dict[str, list[str]] = {
    "OP": ["GSE30335"],
    "AD": ["GSE63060", "GSE63061", "GSE4226", "GSE18309"],
    "PD": ["GSE99039", "GSE6613", "GSE72267", "GSE7621", "GSE20292"],
    "MS": ["GSE17048", "GSE21942", "GSE41890", "GSE43591"],
}
ACCESSION_TO_DISEASE = {
    accession: disease for disease, accessions in ACTIVE_ACCESSIONS.items() for accession in accessions
}

ACCESSION_GROUPS = {
    "OP_blood": ["GSE30335"],
    "AD_blood": ["GSE63060", "GSE63061", "GSE4226", "GSE18309"],
    "PD_blood": ["GSE99039", "GSE6613", "GSE72267"],
    "PD_brain": ["GSE7621", "GSE20292"],
    "MS_blood": ["GSE17048", "GSE21942", "GSE41890", "GSE43591"],
}
HEATMAP_GROUP_ORDER = ["OP_blood", "AD_blood", "PD_blood", "PD_brain", "MS_blood"]
BARPLOT_GROUP_ORDER = ["OP_blood", "AD_blood", "PD_blood", "PD_brain", "MS_blood"]
ALL_ACTIVE_ACCESSIONS = [acc for grp in HEATMAP_GROUP_ORDER for acc in ACCESSION_GROUPS[grp]]

BIOSPECIMEN_TO_TISSUE = {
    "whole_blood": "blood",
    "pbmc": "blood",
    "peripheral_blood_leukocytes": "blood",
    "peripheral_blood_t_cells": "blood",
    "brain_substantia_nigra": "brain",
    "brain_frontal_temporal_cortex": "brain",
    "brain_dlpfc": "brain",
}

GROUP_COLORS = {
    "OP_blood": "#7f7f7f",
    "AD_blood": "#1f77b4",
    "PD_blood": "#ff7f0e",
    "PD_brain": "#ff7f0e",
    "MS_blood": "#2ca02c",
}

DEA_REQUIRED_COLUMNS = {"gene", "log2fc", "padj", "n_case", "n_control"}


@dataclass
class GeneMeta:
    gene: str
    signature_group: str
    expected_direction: str


def warn(message: str) -> None:
    print(f"[WARN] {message}")


def info(message: str) -> None:
    print(f"[INFO] {message}")


def observed_direction(log2fc: float) -> str:
    if pd.isna(log2fc):
        return "flat"
    if log2fc < -0.3:
        return "down"
    if log2fc > 0.3:
        return "up"
    return "flat"


def map_tissue_label(biospecimen: Any) -> str:
    key = str(biospecimen).strip().lower()
    if not key:
        return "unmapped"
    return BIOSPECIMEN_TO_TISSUE.get(key, "unmapped")


def load_targets(path: Path) -> tuple[list[GeneMeta], list[str]]:
    warnings: list[str] = []
    if not path.exists():
        raise FileNotFoundError(f"Missing targets file: {path}")

    with path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    targets = cfg.get("targets", {}) or {}
    if not isinstance(targets, dict):
        raise ValueError("targets.yaml has invalid 'targets' section")

    metas: list[GeneMeta] = []
    for gene in PANEL_GENES:
        raw = targets.get(gene, {})
        if not isinstance(raw, dict):
            raw = {}
            warnings.append(f"{gene} missing metadata in targets.yaml; defaulting expected_direction='none'")

        signature_group = str(raw.get("signature", "unknown"))
        expected_direction = str(raw.get("chronic_direction", "none")).strip().lower() or "none"
        if expected_direction not in {"down", "none"}:
            warnings.append(
                f"{gene} has unexpected chronic_direction='{expected_direction}' in targets.yaml; treating as 'none'"
            )
            expected_direction = "none"
        metas.append(
            GeneMeta(
                gene=gene,
                signature_group=signature_group,
                expected_direction=expected_direction,
            )
        )

    metas.sort(
        key=lambda m: (
            SIGNATURE_ORDER.index(m.signature_group) if m.signature_group in SIGNATURE_ORDER else 99,
            m.gene,
        )
    )
    return metas, warnings


def load_cohort_metadata(path: Path) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    if not path.exists():
        raise FileNotFoundError(f"Missing cohorts metadata file: {path}")

    cohorts = pd.read_csv(path, dtype=str).fillna("")
    cohorts["accession"] = cohorts["accession"].astype(str)
    cohorts = cohorts[cohorts["accession"].isin(ALL_ACTIVE_ACCESSIONS)].copy()

    missing_in_cohorts = [acc for acc in ALL_ACTIVE_ACCESSIONS if acc not in set(cohorts["accession"])]
    for accession in missing_in_cohorts:
        warnings.append(f"{accession} not found in cohorts.csv; tissue_label set to 'unmapped'")

    rows: list[dict[str, str]] = []
    for accession in ALL_ACTIVE_ACCESSIONS:
        match = cohorts[cohorts["accession"] == accession]
        if match.empty:
            biospecimen = ""
        else:
            if len(match) > 1:
                warnings.append(f"{accession} appears multiple times in cohorts.csv; using first row")
            biospecimen = str(match.iloc[0].get("biospecimen", ""))

        tissue_label = map_tissue_label(biospecimen)
        if tissue_label == "unmapped":
            warnings.append(
                f"{accession} has unmapped biospecimen '{biospecimen}'; expected blood/brain mapping"
            )

        rows.append(
            {
                "accession": accession,
                "disease": ACCESSION_TO_DISEASE.get(accession, "UNK"),
                "biospecimen": biospecimen,
                "tissue_label": tissue_label,
            }
        )

    return pd.DataFrame(rows), warnings


def load_dea_file(accession: str, path: Path) -> tuple[pd.DataFrame | None, list[str]]:
    warnings: list[str] = []
    if not path.exists():
        warnings.append(f"Missing DEA file for {accession}: {path}")
        return None, warnings

    try:
        dea = pd.read_csv(path)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"Failed reading {path.name}: {exc}")
        return None, warnings

    if not DEA_REQUIRED_COLUMNS.issubset(dea.columns):
        warnings.append(
            f"{path.name} missing required columns {sorted(DEA_REQUIRED_COLUMNS)}; cohort skipped"
        )
        return None, warnings

    use = dea[list(DEA_REQUIRED_COLUMNS)].copy()
    use["gene"] = use["gene"].astype(str).str.strip().str.upper()
    for col in ["log2fc", "padj", "n_case", "n_control"]:
        use[col] = pd.to_numeric(use[col], errors="coerce")
    use = use.drop_duplicates(subset=["gene"], keep="first")
    return use, warnings


def concordance_value(expected: str, observed_dir: str) -> str:
    if expected == "none":
        return "na"
    if observed_dir == "flat":
        return "flat"
    if observed_dir == expected:
        return "yes"
    return "no"


def build_concordance_table(gene_meta: list[GeneMeta], cohort_meta: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    rows: list[dict[str, Any]] = []

    meta_by_accession = cohort_meta.set_index("accession").to_dict(orient="index")

    for accession in ALL_ACTIVE_ACCESSIONS:
        dea_path = DEA_DIR / f"{accession}_dea.csv"
        dea_df, w = load_dea_file(accession=accession, path=dea_path)
        warnings.extend(w)

        if dea_df is None:
            continue

        cohort_n_case = float(dea_df["n_case"].median()) if dea_df["n_case"].notna().any() else np.nan
        cohort_n_control = (
            float(dea_df["n_control"].median()) if dea_df["n_control"].notna().any() else np.nan
        )

        by_gene = dea_df.set_index("gene")
        cohort_info = meta_by_accession.get(accession, {})
        disease = str(cohort_info.get("disease", ACCESSION_TO_DISEASE.get(accession, "UNK")))
        tissue_label = str(cohort_info.get("tissue_label", "unmapped"))

        for meta in gene_meta:
            if meta.gene in by_gene.index:
                row = by_gene.loc[meta.gene]
                log2fc = float(row["log2fc"]) if pd.notna(row["log2fc"]) else np.nan
                padj = float(row["padj"]) if pd.notna(row["padj"]) else np.nan
                n_case = float(row["n_case"]) if pd.notna(row["n_case"]) else cohort_n_case
                n_control = (
                    float(row["n_control"]) if pd.notna(row["n_control"]) else cohort_n_control
                )
                observed_dir = observed_direction(log2fc)
                concordant = concordance_value(meta.expected_direction, observed_dir)
            else:
                log2fc = np.nan
                padj = np.nan
                n_case = cohort_n_case
                n_control = cohort_n_control
                observed_dir = "absent"
                concordant = "absent"

            rows.append(
                {
                    "gene": meta.gene,
                    "signature_group": meta.signature_group,
                    "expected_direction": meta.expected_direction,
                    "accession": accession,
                    "disease": disease,
                    "tissue_label": tissue_label,
                    "n_case": n_case,
                    "n_control": n_control,
                    "log2fc": log2fc,
                    "padj": padj,
                    "observed_dir": observed_dir,
                    "concordant": concordant,
                }
            )

    table = pd.DataFrame(rows)
    if table.empty:
        return table, warnings

    table["gene"] = pd.Categorical(table["gene"], categories=[m.gene for m in gene_meta], ordered=True)
    table["accession"] = pd.Categorical(table["accession"], categories=ALL_ACTIVE_ACCESSIONS, ordered=True)
    table = table.sort_values(["gene", "accession"], ignore_index=True)
    table["gene"] = table["gene"].astype(str)
    table["accession"] = table["accession"].astype(str)
    return table, warnings


def compute_group_rates(table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group in BARPLOT_GROUP_ORDER:
        accessions = ACCESSION_GROUPS[group]
        sub = table[table["accession"].isin(accessions)]
        yes = int((sub["concordant"] == "yes").sum())
        no = int((sub["concordant"] == "no").sum())
        scoreable = yes + no

        if scoreable > 0:
            rate = yes / scoreable
            ci_low, ci_high = proportion_confint(
                count=yes, nobs=scoreable, alpha=0.05, method="wilson"
            )
        else:
            warn(f"No scoreable gene-cohort pairs in {group}; rate set to NaN")
            rate = np.nan
            ci_low = np.nan
            ci_high = np.nan

        rows.append(
            {
                "group": group,
                "yes": yes,
                "no": no,
                "scoreable_n": scoreable,
                "rate": rate,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "rate_pct": rate * 100 if pd.notna(rate) else np.nan,
                "ci_low_pct": ci_low * 100 if pd.notna(ci_low) else np.nan,
                "ci_high_pct": ci_high * 100 if pd.notna(ci_high) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def write_per_gene_summary(table: pd.DataFrame, gene_meta: list[GeneMeta], out_csv: Path) -> pd.DataFrame:
    groups = ["AD_blood", "PD_blood", "PD_brain", "MS_blood"]
    rows: list[dict[str, Any]] = []

    meta_lookup = {m.gene: m for m in gene_meta}
    for gene in [m.gene for m in gene_meta]:
        meta = meta_lookup[gene]
        row: dict[str, Any] = {
            "gene": gene,
            "signature_group": meta.signature_group,
            "expected_direction": meta.expected_direction,
        }
        sub_gene = table[table["gene"] == gene]

        for group in groups:
            sub = sub_gene[sub_gene["accession"].isin(ACCESSION_GROUPS[group])]
            row[f"{group}_yes"] = int((sub["concordant"] == "yes").sum())
            row[f"{group}_no"] = int((sub["concordant"] == "no").sum())
            row[f"{group}_flat"] = int((sub["concordant"] == "flat").sum())

        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(out_csv, index=False)
    return summary


def plot_gene_cohort_heatmap(table: pd.DataFrame, gene_meta: list[GeneMeta], out_png: Path) -> None:
    row_order = [m.gene for m in gene_meta]
    col_order = ALL_ACTIVE_ACCESSIONS

    fc = table.pivot(index="gene", columns="accession", values="log2fc").reindex(
        index=row_order, columns=col_order
    )
    conc = table.pivot(index="gene", columns="accession", values="concordant").reindex(
        index=row_order, columns=col_order
    )

    disease_by_acc = ACCESSION_TO_DISEASE.copy()
    tissue_by_acc = (
        table.drop_duplicates(subset=["accession"])[["accession", "tissue_label"]]
        .set_index("accession")["tissue_label"]
        .to_dict()
    )
    xlabels = [
        f"{acc}\n{disease_by_acc.get(acc, 'UNK')} {tissue_by_acc.get(acc, 'unmapped')}"
        for acc in col_order
    ]

    row_labels = []
    for m in gene_meta:
        tag = SIGNATURE_TAG.get(m.signature_group, m.signature_group)
        row_labels.append(f"{m.gene} [{tag}]")

    sns.set_theme(style="white", context="paper")
    fig_w = max(12, 0.95 * len(col_order) + 4)
    fig_h = max(6, 0.5 * len(row_order) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)

    sns.heatmap(
        fc,
        cmap="RdBu_r",
        vmin=-2,
        vmax=2,
        linewidths=0.5,
        linecolor="#d9d9d9",
        cbar_kws={"label": "log2FC"},
        ax=ax,
    )

    for i, gene in enumerate(row_order):
        for j, accession in enumerate(col_order):
            val = conc.loc[gene, accession]
            if val == "yes":
                ax.scatter(
                    j + 0.5,
                    i + 0.5,
                    marker="*",
                    s=80,
                    c="white",
                    edgecolors="black",
                    linewidths=0.4,
                    zorder=3,
                )
            elif val == "no":
                ax.scatter(
                    j + 0.5,
                    i + 0.5,
                    marker="x",
                    s=55,
                    c="black",
                    linewidths=1.2,
                    zorder=3,
                )

    boundaries: list[int] = []
    running = 0
    for group in HEATMAP_GROUP_ORDER:
        running += len(ACCESSION_GROUPS[group])
        boundaries.append(running)

    for boundary in boundaries[:-1]:
        ax.axvline(boundary, color="black", linewidth=1.2)

    ax.set_xticklabels(xlabels, rotation=0, fontsize=8)
    ax.set_yticklabels(row_labels, rotation=0, fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Panel Gene Concordance with OP Chronic Signature \u2014 GEO Disease Cohorts",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_concordance_bars(rates: pd.DataFrame, out_png: Path) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(9.5, 5.0), dpi=150)

    x = np.arange(len(BARPLOT_GROUP_ORDER))
    y = np.array(
        [
            float(rates.loc[rates["group"] == grp, "rate_pct"].iloc[0])
            if not rates.loc[rates["group"] == grp, "rate_pct"].empty
            else np.nan
            for grp in BARPLOT_GROUP_ORDER
        ]
    )

    ci_low = np.array(
        [
            float(rates.loc[rates["group"] == grp, "ci_low_pct"].iloc[0])
            if not rates.loc[rates["group"] == grp, "ci_low_pct"].empty
            else np.nan
            for grp in BARPLOT_GROUP_ORDER
        ]
    )
    ci_high = np.array(
        [
            float(rates.loc[rates["group"] == grp, "ci_high_pct"].iloc[0])
            if not rates.loc[rates["group"] == grp, "ci_high_pct"].empty
            else np.nan
            for grp in BARPLOT_GROUP_ORDER
        ]
    )
    yerr = np.vstack(
        [
            np.nan_to_num(y - ci_low, nan=0.0),
            np.nan_to_num(ci_high - y, nan=0.0),
        ]
    )

    bar_colors = [GROUP_COLORS.get(grp, "#999999") for grp in BARPLOT_GROUP_ORDER]
    bars = ax.bar(
        x,
        np.nan_to_num(y, nan=0.0),
        color=bar_colors,
        yerr=yerr,
        capsize=4,
        edgecolor="black",
        linewidth=0.6,
    )

    for idx, grp in enumerate(BARPLOT_GROUP_ORDER):
        n = int(rates.loc[rates["group"] == grp, "scoreable_n"].iloc[0])
        y_val = y[idx] if np.isfinite(y[idx]) else 0.0
        ax.text(
            idx,
            y_val + 2.5,
            f"n={n} scoreable gene-cohort pairs",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=0,
        )
        if not np.isfinite(y[idx]):
            bars[idx].set_alpha(0.35)

    ax.axhline(50, color="black", linestyle="--", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(BARPLOT_GROUP_ORDER, rotation=0)
    ax.set_ylabel("Concordance rate (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Concordance Rate with OP Chronic Proteomics Signature by Tissue")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_warnings: list[str] = []

    try:
        gene_meta, warnings = load_targets(TARGETS_YAML)
        all_warnings.extend(warnings)
    except Exception as exc:  # noqa: BLE001
        warn(f"Cannot continue without targets metadata: {exc}")
        return 0

    try:
        cohort_meta, warnings = load_cohort_metadata(COHORTS_CSV)
        all_warnings.extend(warnings)
    except Exception as exc:  # noqa: BLE001
        warn(f"Cannot continue without cohorts metadata: {exc}")
        return 0

    try:
        table, warnings = build_concordance_table(gene_meta=gene_meta, cohort_meta=cohort_meta)
        all_warnings.extend(warnings)
    except Exception as exc:  # noqa: BLE001
        warn(f"Failed while building concordance table: {exc}")
        return 0

    if table.empty:
        warn("Concordance table is empty; no outputs produced")
        return 0

    out_table = OUT_DIR / "concordance_table.csv"
    table.to_csv(out_table, index=False)
    info(f"Wrote {out_table}")

    try:
        out_heatmap = OUT_DIR / "gene_cohort_heatmap.png"
        plot_gene_cohort_heatmap(table=table, gene_meta=gene_meta, out_png=out_heatmap)
        info(f"Wrote {out_heatmap}")
    except Exception as exc:  # noqa: BLE001
        warn(f"Failed to generate heatmap: {exc}")

    try:
        rates = compute_group_rates(table)
        out_bar = OUT_DIR / "concordance_rate_barplot.png"
        plot_concordance_bars(rates=rates, out_png=out_bar)
        info(f"Wrote {out_bar}")

        print("\nConcordance rates (Wilson 95% CI):")
        printable = rates.copy()
        printable["rate_pct"] = printable["rate_pct"].round(2)
        printable["ci_low_pct"] = printable["ci_low_pct"].round(2)
        printable["ci_high_pct"] = printable["ci_high_pct"].round(2)
        print(
            printable[
                ["group", "yes", "no", "scoreable_n", "rate_pct", "ci_low_pct", "ci_high_pct"]
            ].to_string(index=False)
        )
    except Exception as exc:  # noqa: BLE001
        warn(f"Failed to compute/plot concordance rates: {exc}")
        rates = pd.DataFrame()

    try:
        out_gene_summary = OUT_DIR / "per_gene_concordance_summary.csv"
        write_per_gene_summary(table=table, gene_meta=gene_meta, out_csv=out_gene_summary)
        info(f"Wrote {out_gene_summary}")
    except Exception as exc:  # noqa: BLE001
        warn(f"Failed to write per-gene concordance summary: {exc}")

    for message in all_warnings:
        warn(message)
    info(f"Completed with {len(all_warnings)} warnings")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
