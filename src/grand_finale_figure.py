#!/usr/bin/env python3
"""
Grand Finale Figure
-------------------
Chronic Organophosphate Exposure Protein Signature Validated Across
Neuroinflammatory Disease Proteomics (AD / PD / MS).

Layout (left → right):
  [OP Chronic Signal] | [AD Brain x2] [PD Brain] | [AD CSF] [PD CSF] [MS CSF]

Concordance markers on disease columns reflect agreement with OP exposure direction.
"""

from __future__ import annotations
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ["ACTG1", "DNAH9", "GPX3", "VWF", "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B"]

# OP chronic log2FC from experimental proteomics data (chronic vs control)
OP_FC: dict[str, float] = {
    "ACTG1": -5.98, "DNAH9":  5.84, "GPX3": -2.83, "VWF":  4.28,
    "C4B":   -4.15, "CD44":   1.58, "CFHR2": 1.47, "ITIH3": -1.21,
    "LRG1":  -1.72, "MYH7B": -5.32,
}

# OP direction per gene (threshold 0.5 — all 10 markers exceed this)
OP_DIR: dict[str, str] = {
    g: ("up" if v > 0.5 else "down" if v < -0.5 else "flat")
    for g, v in OP_FC.items()
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def direction_from_fc(fc: float | None, thr: float = 0.3) -> str:
    if fc is None or (isinstance(fc, float) and np.isnan(fc)):
        return "flat"
    return "down" if fc < -thr else ("up" if fc > thr else "flat")


def concordance_vs_op(gene: str, fc: float | None) -> str:
    """Does disease log2FC direction agree with OP direction?"""
    op_d = OP_DIR.get(gene, "flat")
    dis_d = direction_from_fc(fc)
    if op_d == "flat" or dis_d == "flat":
        return "na"
    return "yes" if dis_d == op_d else "no"


def is_present(val: object) -> bool:
    return str(val).strip().lower() in ("true", "1", "yes")


def build_eprot_col(eprot_df: pd.DataFrame, accession: str) -> pd.Series:
    sub = eprot_df[eprot_df["accession"] == accession]
    s = pd.Series(dtype=float, index=PANEL)
    for gene in PANEL:
        row = sub[sub["gene"] == gene]
        if not row.empty and bool(row.iloc[0]["present"]) and pd.notna(row.iloc[0]["log2fc"]):
            s[gene] = float(row.iloc[0]["log2fc"])
    return s


def build_pride_col(df: pd.DataFrame, accession: str, log2fc_col: str = "log2fc") -> pd.Series:
    sub = df[df["accession"] == accession]
    s = pd.Series(dtype=float, index=PANEL)
    for gene in PANEL:
        row = sub[sub["gene"] == gene]
        if row.empty:
            continue
        r = row.iloc[0]
        if not is_present(r.get("present", True)):
            continue
        if log2fc_col in r.index and pd.notna(r[log2fc_col]):
            s[gene] = float(r[log2fc_col])
        elif "mean_case" in r.index and "mean_ctrl" in r.index:
            mc, mctl = r["mean_case"], r["mean_ctrl"]
            if pd.notna(mc) and pd.notna(mctl):
                s[gene] = float(mc) - float(mctl)
    return s


def make_annotations(mat: pd.DataFrame) -> pd.DataFrame:
    ann = pd.DataFrame("", index=mat.index, columns=mat.columns, dtype=object)
    for gene in mat.index:
        for col in mat.columns:
            val = mat.loc[gene, col]
            if pd.isna(val):
                continue
            c = concordance_vs_op(gene, float(val))
            marker = " ●" if c == "yes" else (" ○" if c == "no" else "")
            ann.loc[gene, col] = f"{val:.1f}{marker}"
    return ann


def make_op_annotations(mat: pd.DataFrame) -> pd.DataFrame:
    ann = pd.DataFrame("", index=mat.index, columns=mat.columns, dtype=object)
    for gene in mat.index:
        val = mat.loc[gene, mat.columns[0]]
        if pd.notna(val):
            ann.loc[gene, mat.columns[0]] = f"{val:.1f}"
    return ann


def plot_heat(ax: plt.Axes, data: pd.DataFrame, ann: pd.DataFrame,
              cmap: str, vmin: float, vmax: float, center: float | None,
              show_yticks: bool = False, annot_fs: int = 8) -> None:
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#e8e8e8")
    sns.heatmap(
        data, ax=ax, cmap=cm, vmin=vmin, vmax=vmax, center=center,
        mask=data.isna(), cbar=False,
        linewidths=0.5, linecolor="#d0d0d0",
        annot=ann, fmt="",
        annot_kws={"fontsize": annot_fs, "color": "black"},
        yticklabels=PANEL if show_yticks else False,
        xticklabels=list(data.columns),
    )
    if show_yticks:
        ax.set_yticklabels(PANEL, rotation=0, fontsize=11, fontweight="bold")
        ax.tick_params(axis="y", left=True, labelleft=True, length=0)
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelrotation=40, labelsize=8, length=0)
    ax.set_xlabel(""); ax.set_ylabel("")
    for sp in ax.spines.values():
        sp.set_visible(False)


def add_strip(ax: plt.Axes, label: str, color: str, fs: int = 9) -> None:
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(patches.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="none"))
    ax.text(0.5, 0.5, label, ha="center", va="center",
            fontsize=fs, color="white", fontweight="bold", linespacing=1.3)
    ax.axis("off")


def add_section_brace(fig: plt.Figure, ax_left: plt.Axes, ax_right: plt.Axes,
                      label: str, color: str, y_fig: float) -> None:
    """Draw a colored section bracket above a group of axes."""
    x0 = ax_left.get_position().x0
    x1 = ax_right.get_position().x1
    fig.patches.append(patches.FancyArrowPatch(
        (x0, y_fig), (x1, y_fig),
        arrowstyle="-", color=color, linewidth=2.5,
        transform=fig.transFigure, clip_on=False,
    ))
    fig.text((x0 + x1) / 2, y_fig + 0.008, label,
             ha="center", va="bottom", fontsize=10,
             color=color, fontweight="bold", transform=fig.transFigure)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    eprot_df = pd.read_csv(ROOT / "results" / "eprot_query.csv")
    pride_df  = pd.read_csv(ROOT / "results" / "pride_discovery" / "pride_ad_pd_ms_query.csv")
    pride_old = pd.read_csv(ROOT / "results" / "pride_query.csv")

    # --- OP column ---
    op_mat = pd.DataFrame(
        {"OP\nChronic": [OP_FC[g] for g in PANEL]}, index=PANEL
    )
    op_ann = make_op_annotations(op_mat)

    # --- E-PROT Brain columns ---
    eprot_mat = pd.DataFrame(index=PANEL)
    eprot_mat["AD Brain\n(E-PROT-31)"] = build_eprot_col(eprot_df, "E-PROT-31")
    eprot_mat["AD Brain\n(E-PROT-61)"] = build_eprot_col(eprot_df, "E-PROT-61")
    eprot_mat["PD Brain\n(E-PROT-65)"] = build_eprot_col(eprot_df, "E-PROT-65")
    eprot_ann = make_annotations(eprot_mat)

    # --- CSF Proteomics columns ---
    csf_mat = pd.DataFrame(index=PANEL)
    csf_mat["AD CSF\n(PXD016278)"] = build_pride_col(pride_old, "PXD016278")
    csf_mat["PD CSF\n(PXD026491)"] = build_pride_col(pride_df, "PXD026491")
    csf_mat["MS CSF\n(PXD045058)"] = build_pride_col(pride_df, "PXD045058")
    csf_ann = make_annotations(csf_mat)

    # --- Figure layout ---
    # 7 columns: [OP] [AD-31] [AD-61] [PD-65] [AD-CSF] [PD-CSF] [MS-CSF]
    sns.set_theme(style="white")
    fig = plt.figure(figsize=(22, 10), dpi=150)
    gs = fig.add_gridspec(
        nrows=2, ncols=7,
        height_ratios=[0.35, 9.65],
        width_ratios=[2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
        hspace=0.06, wspace=0.06,
    )

    # Header color strips
    STRIP_META = [
        (0, "Chronic OP\nSignal",       "#0d1b4b"),   # deep navy
        (1, "AD Brain\nE-PROT-31",      "#7b0000"),   # dark crimson
        (2, "AD Brain\nE-PROT-61",      "#9b1414"),   # crimson
        (3, "PD Brain\nE-PROT-65",      "#8B4513"),   # saddle brown
        (4, "AD CSF\nPXD016278",        "#1a4a7a"),   # dark steel blue
        (5, "PD CSF\nPXD026491",        "#2c6e9e"),   # medium blue
        (6, "MS CSF\nPXD045058",        "#1b6ca8"),   # royal blue
    ]
    for col_i, label, color in STRIP_META:
        sax = fig.add_subplot(gs[0, col_i])
        add_strip(sax, label, color, fs=8)

    # Heatmap axes (shared y)
    ax0 = fig.add_subplot(gs[1, 0])
    axes = [ax0] + [fig.add_subplot(gs[1, i], sharey=ax0) for i in range(1, 7)]

    # OP column — wider diverging scale (±7)
    plot_heat(axes[0], op_mat, op_ann, "RdBu_r", -7, 7, 0,
              show_yticks=True, annot_fs=9)

    # E-PROT columns — disease scale (±3)
    for i, col in enumerate(eprot_mat.columns):
        plot_heat(axes[1 + i],
                  eprot_mat[[col]], eprot_ann[[col]],
                  "RdBu_r", -3, 3, 0, annot_fs=8)

    # CSF columns — disease scale (±3)
    for i, col in enumerate(csf_mat.columns):
        plot_heat(axes[4 + i],
                  csf_mat[[col]], csf_ann[[col]],
                  "RdBu_r", -3, 3, 0, annot_fs=8)

    # Vertical separators between sections
    for ax_idx in [1, 4]:
        axes[ax_idx].axvline(0, color="#333333", linewidth=2.5, clip_on=False)

    # Section brace labels (drawn after tight_layout so positions are final)
    plt.subplots_adjust(left=0.11, right=0.99, bottom=0.20, top=0.88)

    # Section brackets
    add_section_brace(fig, axes[1], axes[3],
                      "Brain Proteomics  (Expression Atlas E-PROT)",
                      "#7b0000", 0.905)
    add_section_brace(fig, axes[4], axes[6],
                      "CSF Proteomics  (PRIDE)",
                      "#1a4a7a", 0.905)

    # OP label
    x_op = (axes[0].get_position().x0 + axes[0].get_position().x1) / 2
    fig.text(x_op, 0.905, "OP Exposure", ha="center", va="bottom",
             fontsize=10, color="#0d1b4b", fontweight="bold",
             transform=fig.transFigure)

    # Colorbar legend (manual)
    # Disease scale bar
    sm_dis = plt.cm.ScalarMappable(cmap="RdBu_r",
                                    norm=plt.Normalize(vmin=-3, vmax=3))
    sm_dis.set_array([])
    cax_dis = fig.add_axes([0.60, 0.035, 0.18, 0.022])
    cb_dis = fig.colorbar(sm_dis, cax=cax_dis, orientation="horizontal")
    cb_dis.set_label("Disease log₂FC  (±3)", fontsize=8, labelpad=2)
    cb_dis.ax.tick_params(labelsize=7)

    # OP scale bar
    sm_op = plt.cm.ScalarMappable(cmap="RdBu_r",
                                   norm=plt.Normalize(vmin=-7, vmax=7))
    sm_op.set_array([])
    cax_op = fig.add_axes([0.11, 0.035, 0.14, 0.022])
    cb_op = fig.colorbar(sm_op, cax=cax_op, orientation="horizontal")
    cb_op.set_label("OP log₂FC  (±7)", fontsize=8, labelpad=2)
    cb_op.ax.tick_params(labelsize=7)

    # Concordance legend
    fig.text(0.82, 0.01,
             "● concordant with OP direction   ○ discordant",
             ha="left", va="bottom", fontsize=9, style="italic", color="#333333")

    # Title
    fig.suptitle(
        "Chronic Organophosphate Exposure Proteomics Signature Validated Across"
        " Neuroinflammatory Disease Datasets",
        fontsize=14, fontweight="bold", y=0.975,
    )

    # Subtitle
    fig.text(0.5, 0.945,
             "88–336 cases vs controls per dataset  |  10-protein OP panel  |"
             "  AD, PD, MS brain & CSF proteomics",
             ha="center", va="top", fontsize=9, color="#555555")

    out_path = OUT / "grand_finale_op_neurodegeneration.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_path}")
    print(f"Size:  {out_path.stat().st_size // 1024} KB")

    # Print summary
    print("\n=== Concordance summary across all disease datasets ===")
    all_conc = {"yes": 0, "no": 0, "na": 0}
    for mat, label in [(eprot_mat, "E-PROT"), (csf_mat, "CSF")]:
        for col in mat.columns:
            for gene in PANEL:
                fc = mat.loc[gene, col]
                c = concordance_vs_op(gene, float(fc) if pd.notna(fc) else None)
                all_conc[c] = all_conc.get(c, 0) + 1
    print(f"  Concordant: {all_conc['yes']}")
    print(f"  Discordant: {all_conc['no']}")
    print(f"  NA (flat/absent): {all_conc['na']}")


if __name__ == "__main__":
    main()
