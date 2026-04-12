#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
import yaml


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

E_PROT_AD_ITEMS = [
    ("E-PROT-31", "AD(atlas)"),
    ("E-PROT-61", "AD(multi)"),
    ("E-PROT-57", "AD(Sinai)"),
    ("E-PROT-32", "AD(Braak)"),
]

E_PROT_PD_ITEMS = [("E-PROT-65", "PD brain")]

E_PROT_FTD_ITEMS = [
    ("E-PROT-137", "FTD(frontal)"),
    ("E-PROT-147", "FTD(GRN/MAPT)"),
    ("E-PROT-148", "FTD(SD)"),
]

PRIDE_AD_ITEMS = [("PXD016278", "AD CSF")]
PRIDE_PD_ITEMS = [("PXD026491", "PD CSF")]

HPA_ACCESSION = "HPA_BLOOD_MS_CONC_BLOOD_IMMUNE"
HPA_LABEL = "HPA plasma"


def is_present(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"true", "1", "yes", "y", "t"}


def normalize_concordance(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"yes", "true", "1", "y"}:
        return "yes"
    if text in {"no", "false", "0", "n"}:
        return "no"
    return None


def reduce_concordance(values: Iterable[object]) -> str | None:
    normalized = [c for c in (normalize_concordance(v) for v in values) if c]
    if not normalized:
        return None
    yes_n = sum(c == "yes" for c in normalized)
    no_n = sum(c == "no" for c in normalized)
    return "yes" if yes_n >= no_n else "no"


def concordance_marker(value: object) -> str:
    concord = normalize_concordance(value)
    if concord == "yes":
        return " ●"
    if concord == "no":
        return " ○"
    return ""


def build_source_matrix(
    df: pd.DataFrame,
    *,
    source_col: str,
    mapping_items: list[tuple[str, str]],
    value_col: str,
    panel: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    label_map = dict(mapping_items)
    label_order = [label for _, label in mapping_items]

    if "present" in df.columns:
        df = df[df["present"].map(is_present)]
    df = df[df["gene"].isin(panel)].copy()
    df["label"] = df[source_col].map(label_map)
    df = df[df["label"].notna()].copy()

    if value_col not in df.columns:
        values = pd.DataFrame(index=panel, columns=label_order, dtype=float)
        concord = pd.DataFrame(index=panel, columns=label_order, dtype=object)
        return values, concord

    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    values = (
        df.groupby(["gene", "label"], dropna=False)[value_col]
        .mean()
        .unstack("label")
        .reindex(index=panel, columns=label_order)
    )

    if "concordant" in df.columns:
        concord = (
            df.groupby(["gene", "label"], dropna=False)["concordant"]
            .apply(reduce_concordance)
            .unstack("label")
            .reindex(index=panel, columns=label_order)
        )
    else:
        concord = pd.DataFrame(index=panel, columns=label_order, dtype=object)

    return values, concord


def make_numeric_annotations(
    values: pd.DataFrame,
    *,
    concordance: pd.DataFrame | None = None,
    decimals: int = 1,
) -> pd.DataFrame:
    ann = pd.DataFrame("", index=values.index, columns=values.columns, dtype=object)
    for row in values.index:
        for col in values.columns:
            val = values.loc[row, col]
            if pd.isna(val):
                continue
            base = f"{val:.{decimals}f}"
            if concordance is not None:
                base += concordance_marker(concordance.loc[row, col])
            ann.loc[row, col] = base
    return ann


def make_op_reference(config_path: Path, panel: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    with config_path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    targets = config.get("targets", {})
    direction_to_num = {"up": 1.0, "down": -1.0, "none": 0.0}

    values = pd.DataFrame(index=panel, columns=["OP Ref"], dtype=float)
    ann = pd.DataFrame("", index=panel, columns=["OP Ref"], dtype=object)

    for gene in panel:
        direction = str(targets.get(gene, {}).get("chronic_direction", "none")).strip().lower()
        if direction not in direction_to_num:
            direction = "none"
        values.loc[gene, "OP Ref"] = direction_to_num[direction]

    return values, ann


def plot_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    annotations: pd.DataFrame,
    *,
    cmap_name: str,
    vmin: float | None,
    vmax: float | None,
    center: float | None,
    separator: bool,
    annotate: bool = True,
) -> None:
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("white")

    sns.heatmap(
        data,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        mask=data.isna(),
        cbar=False,
        linewidths=0.8,
        linecolor="lightgray",
        annot=annotations if annotate else False,
        fmt="",
        annot_kws={"fontsize": 8, "color": "black"},
        yticklabels=False,
        xticklabels=list(data.columns),
    )

    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    if separator:
        ax.axvline(data.shape[1], color="black", linewidth=1.2)


def add_group_strip(ax: plt.Axes, label: str, color: str) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(patches.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="none"))
    ax.text(
        0.5,
        0.5,
        label,
        ha="center",
        va="center",
        fontsize=10,
        color="white",
        fontweight="bold",
    )
    ax.axis("off")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "config" / "targets.yaml"
    eprot_path = root / "results" / "eprot_query.csv"
    pride_path = root / "results" / "pride_discovery" / "pride_ad_pd_ms_query.csv"
    hpa_path = root / "results" / "hpa_query.csv"

    op_values, op_ann = make_op_reference(config_path, PANEL)

    eprot_df = pd.read_csv(eprot_path)
    ad_values, ad_conc = build_source_matrix(
        eprot_df,
        source_col="accession",
        mapping_items=E_PROT_AD_ITEMS,
        value_col="log2fc",
        panel=PANEL,
    )
    ad_ann = make_numeric_annotations(ad_values, concordance=ad_conc, decimals=1)

    pd_brain_values, pd_brain_conc = build_source_matrix(
        eprot_df,
        source_col="accession",
        mapping_items=E_PROT_PD_ITEMS,
        value_col="log2fc",
        panel=PANEL,
    )
    pd_brain_ann = make_numeric_annotations(pd_brain_values, concordance=pd_brain_conc, decimals=1)

    ftd_values, ftd_conc = build_source_matrix(
        eprot_df,
        source_col="accession",
        mapping_items=E_PROT_FTD_ITEMS,
        value_col="log2fc",
        panel=PANEL,
    )
    ftd_ann = make_numeric_annotations(ftd_values, concordance=ftd_conc, decimals=1)

    pride_df = pd.read_csv(pride_path)
    ad_csf_values, ad_csf_conc = build_source_matrix(
        pride_df,
        source_col="accession",
        mapping_items=PRIDE_AD_ITEMS,
        panel=PANEL,
        value_col="log2fc",
    )
    ad_csf_ann = make_numeric_annotations(ad_csf_values, concordance=ad_csf_conc, decimals=1)

    pd_csf_values, pd_csf_conc = build_source_matrix(
        pride_df,
        source_col="accession",
        mapping_items=PRIDE_PD_ITEMS,
        panel=PANEL,
        value_col="log2fc",
    )
    pd_csf_ann = make_numeric_annotations(pd_csf_values, concordance=pd_csf_conc, decimals=1)

    hpa_df = pd.read_csv(hpa_path)
    if "present" in hpa_df.columns:
        hpa_df = hpa_df[hpa_df["present"].map(is_present)]
    hpa_df = hpa_df[
        (hpa_df["accession"] == HPA_ACCESSION) & (hpa_df["gene"].isin(PANEL))
    ].copy()
    hpa_df["mean_case"] = pd.to_numeric(hpa_df["mean_case"], errors="coerce")
    hpa_df["hpa_log10"] = np.where(hpa_df["mean_case"] > 0, np.log10(hpa_df["mean_case"]), np.nan)

    hpa_values = (
        hpa_df.groupby("gene", dropna=False)["hpa_log10"]
        .mean()
        .reindex(PANEL)
        .to_frame(HPA_LABEL)
    )
    hpa_ann = make_numeric_annotations(hpa_values, concordance=None, decimals=1)

    hpa_min = np.nanmin(hpa_values.to_numpy()) if np.isfinite(hpa_values.to_numpy()).any() else 0.0
    hpa_max = np.nanmax(hpa_values.to_numpy()) if np.isfinite(hpa_values.to_numpy()).any() else 1.0
    if np.isclose(hpa_min, hpa_max):
        hpa_min -= 0.5
        hpa_max += 0.5

    sns.set_theme(style="white")
    fig = plt.figure(figsize=(20, 9), dpi=150)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=7,
        height_ratios=[0.45, 9.55],
        width_ratios=[1, 1, 4, 1, 3, 1, 1],
        hspace=0.08,
        wspace=0.06,
    )

    strip_axes = [fig.add_subplot(gs[0, i]) for i in range(7)]
    ax0 = fig.add_subplot(gs[1, 0])
    heat_axes = [
        ax0,
        fig.add_subplot(gs[1, 1], sharey=ax0),
        fig.add_subplot(gs[1, 2], sharey=ax0),
        fig.add_subplot(gs[1, 3], sharey=ax0),
        fig.add_subplot(gs[1, 4], sharey=ax0),
        fig.add_subplot(gs[1, 5], sharey=ax0),
        fig.add_subplot(gs[1, 6], sharey=ax0),
    ]

    plot_panel(
        heat_axes[0],
        op_values,
        op_ann,
        cmap_name="RdBu_r",
        vmin=-1,
        vmax=1,
        center=0,
        separator=True,
        annotate=False,
    )
    plot_panel(
        heat_axes[1],
        hpa_values,
        hpa_ann,
        cmap_name="YlOrRd",
        vmin=float(hpa_min),
        vmax=float(hpa_max),
        center=None,
        separator=True,
    )
    plot_panel(
        heat_axes[2],
        ad_values,
        ad_ann,
        cmap_name="RdBu_r",
        vmin=-3,
        vmax=3,
        center=0,
        separator=True,
    )
    plot_panel(
        heat_axes[3],
        pd_brain_values,
        pd_brain_ann,
        cmap_name="RdBu_r",
        vmin=-3,
        vmax=3,
        center=0,
        separator=True,
    )
    plot_panel(
        heat_axes[4],
        ftd_values,
        ftd_ann,
        cmap_name="RdBu_r",
        vmin=-3,
        vmax=3,
        center=0,
        separator=True,
    )
    plot_panel(
        heat_axes[5],
        ad_csf_values,
        ad_csf_ann,
        cmap_name="RdBu_r",
        vmin=-3,
        vmax=3,
        center=0,
        separator=True,
    )
    plot_panel(
        heat_axes[6],
        pd_csf_values,
        pd_csf_ann,
        cmap_name="RdBu_r",
        vmin=-3,
        vmax=3,
        center=0,
        separator=False,
    )

    # Force row labels on left-most shared axis after all panels are built.
    heat_axes[0].set_yticks(np.arange(len(PANEL)) + 0.5)
    heat_axes[0].set_yticklabels(PANEL, rotation=0, fontsize=10)
    heat_axes[0].tick_params(axis="y", which="both", left=True, labelleft=True)

    strip_colors = ["#1a3a5c", "#D98E04", "#B22222", "#8B0000", "#6B0000", "#4682B4", "#4169E1"]
    strip_labels = [
        "OP Ref",
        "HPA Normal plasma",
        "E-PROT AD brain",
        "E-PROT PD brain",
        "E-PROT FTD brain",
        "PRIDE AD CSF",
        "PRIDE PD CSF",
    ]
    for sax, label, color in zip(strip_axes, strip_labels, strip_colors):
        add_group_strip(sax, label, color)

    fig.suptitle("OP Chronic Signature vs Neurodegenerative Disease Proteomics", fontsize=16, y=0.98)
    plt.subplots_adjust(left=0.10, right=0.99, bottom=0.20, top=0.90)

    out_dir = root / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "cross_modal_panel_heatmap.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    zero_present_columns: list[str] = []
    for df in [hpa_values, ad_values, pd_brain_values, ftd_values, ad_csf_values, pd_csf_values]:
        for col in df.columns:
            if int(df[col].notna().sum()) == 0:
                zero_present_columns.append(col)

    print("Saved: results/figures/cross_modal_panel_heatmap.png")
    if zero_present_columns:
        print("Columns with 0 present genes:", ", ".join(zero_present_columns))
    else:
        print("Columns with 0 present genes: none")


if __name__ == "__main__":
    main()
