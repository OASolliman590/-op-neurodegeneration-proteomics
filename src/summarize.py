"""
summarize.py
------------
Aggregate all query CSVs into one master table, then produce:

  1. results/master_query.csv  — all gene × dataset rows
  2. results/presence_heatmap.png  — gene × dataset presence (green = present)
  3. results/direction_heatmap.png — gene × dataset direction
     (blue = down as expected, red = up/discordant, grey = absent/no data)
  4. results/summary_table.csv — per-gene: n_datasets_present, n_concordant,
                                  n_discordant, n_expected_none, concordance_rate

Run after all query scripts have completed.
"""

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "results"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",      # shared stable-down
    "C4B", "CD44",                           # acute-specific / reversal
    "CFHR2",                                 # shared stable-down
    "ITIH3", "LRG1", "MYH7B",              # chronic-specific
]

SIGNATURE = {
    "ACTG1": "shared↓",  "DNAH9": "shared↓", "GPX3": "shared↓",
    "VWF":   "shared↓",  "CFHR2": "shared↓",
    "C4B":   "reversal", "CD44":  "reversal",
    "ITIH3": "chronic↓", "LRG1":  "chronic↓", "MYH7B": "chronic↓",
}

SIG_COLOR = {"shared↓": "#1f77b4", "reversal": "#ff7f0e", "chronic↓": "#2ca02c"}

DIR_COLOR = {
    "down_concordant":   "#2166ac",   # blue — down, as expected
    "up_discordant":     "#d73027",   # red — up, wrong direction
    "down_discordant":   "#f4a582",   # salmon — down but expected up
    "up_concordant":     "#4dac26",   # green — up, as expected
    "flat":              "#ffffbf",   # pale yellow
    "absent":            "#d9d9d9",   # grey
    "no_data":           "#d9d9d9",
    "expected_none":     "#e0e0e0",   # light grey — direction is N/A for this gene/cohort
}


def load_all(results_dir: Path) -> pd.DataFrame:
    frames = []
    excluded = {"master_query.csv"}
    for csv in results_dir.glob("*_query.csv"):
        if csv.name in excluded:
            continue
        try:
            df = pd.read_csv(csv)
            df["source_file"] = csv.name
            frames.append(df)
            log.info(f"  loaded {csv.name}: {len(df)} rows")
        except Exception as e:
            log.warning(f"  cannot load {csv.name}: {e}")
    if not frames:
        log.error("No query CSVs found. Run query scripts first.")
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def cohort_label(row: pd.Series) -> str:
    return f"{row['accession']}\n({row['disease']})"


def dir_cell_color(row: pd.Series) -> str:
    d   = str(row.get("direction", "absent")).lower()
    con = str(row.get("concordant", "na")).lower()
    exp = str(row.get("expected",   "none")).lower()

    if d == "absent":
        return "absent"
    if d == "no_data" or d == "no_labels":
        return "no_data"
    if exp == "none":
        return "expected_none"
    if d == "down":
        return "down_concordant" if con == "yes" else "down_discordant"
    if d == "up":
        return "up_concordant"   if con == "yes" else "up_discordant"
    return "flat"


def make_presence_heatmap(df: pd.DataFrame, out_path: Path):
    cohorts = df["accession"].unique().tolist()
    pivot = pd.DataFrame(index=PANEL, columns=cohorts, dtype=object)

    for _, row in df.iterrows():
        gene = row.get("gene")
        acc  = row.get("accession")
        if gene in PANEL and acc in cohorts:
            pivot.at[gene, acc] = bool(row.get("present", False))

    # Numerical for coloring
    numeric = pivot.map(lambda v: 1 if v is True else (0 if v is False else -1))

    fig, ax = plt.subplots(figsize=(max(8, len(cohorts) * 0.9), 5))
    cmap = plt.cm.RdYlGn

    im = ax.imshow(numeric.values.astype(float), cmap=cmap, vmin=-1, vmax=1, aspect="auto")

    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(PANEL)))

    # Colour gene labels by signature
    ax.set_yticklabels(PANEL, fontsize=9)
    for tick, gene in zip(ax.get_yticklabels(), PANEL):
        tick.set_color(SIG_COLOR.get(SIGNATURE.get(gene, ""), "black"))

    ax.set_title("Panel gene presence across validation datasets\n"
                 "(green=present, red=absent, grey=not queried)", fontsize=10)

    # Annotate cells
    for i, gene in enumerate(PANEL):
        for j, acc in enumerate(cohorts):
            val = numeric.iat[i, j]
            txt = "✓" if val == 1 else ("✗" if val == 0 else "?")
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if val == 1 else "black")

    # Signature legend
    patches = [mpatches.Patch(color=c, label=s) for s, c in SIG_COLOR.items()]
    ax.legend(handles=patches, loc="lower right", fontsize=7,
              title="Signature", title_fontsize=7)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Saved presence heatmap: {out_path}")


def make_direction_heatmap(df: pd.DataFrame, out_path: Path):
    cohorts = df["accession"].unique().tolist()
    pivot = pd.DataFrame(index=PANEL, columns=cohorts, data="absent", dtype=object)

    for _, row in df.iterrows():
        gene = row.get("gene")
        acc  = row.get("accession")
        if gene in PANEL and acc in cohorts:
            pivot.at[gene, acc] = dir_cell_color(row)

    fig, ax = plt.subplots(figsize=(max(8, len(cohorts) * 0.9), 5))
    ax.set_xlim(-0.5, len(cohorts) - 0.5)
    ax.set_ylim(-0.5, len(PANEL) - 0.5)
    ax.set_aspect("auto")

    for i, gene in enumerate(PANEL):
        for j, acc in enumerate(cohorts):
            color_key = pivot.at[gene, acc]
            color = DIR_COLOR.get(color_key, "#d9d9d9")
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                   facecolor=color, edgecolor="white", linewidth=0.5)
            ax.add_patch(rect)

            # Direction label inside cell
            d_map = {
                "down_concordant": "↓✓", "up_concordant":   "↑✓",
                "down_discordant": "↓✗", "up_discordant":   "↑✗",
                "flat": "—", "absent": "", "no_data": "?", "expected_none": "·",
            }
            txt = d_map.get(color_key, "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                    color="white" if "concordant" in color_key else "black")

    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(PANEL)))
    ax.set_yticklabels(PANEL, fontsize=9)
    for tick, gene in zip(ax.get_yticklabels(), PANEL):
        tick.set_color(SIG_COLOR.get(SIGNATURE.get(gene, ""), "black"))

    ax.set_title("Panel gene direction across validation datasets\n"
                 "↓✓ down-concordant  ↑✗ up-discordant  · expected N/A", fontsize=9)

    legend_items = [
        mpatches.Patch(color=DIR_COLOR["down_concordant"], label="down (concordant)"),
        mpatches.Patch(color=DIR_COLOR["up_discordant"],   label="up (discordant)"),
        mpatches.Patch(color=DIR_COLOR["down_discordant"], label="down (discordant)"),
        mpatches.Patch(color=DIR_COLOR["up_concordant"],   label="up (concordant)"),
        mpatches.Patch(color=DIR_COLOR["flat"],            label="flat"),
        mpatches.Patch(color=DIR_COLOR["absent"],          label="absent / not found"),
        mpatches.Patch(color=DIR_COLOR["expected_none"],   label="no expected direction"),
    ]
    ax.legend(handles=legend_items, loc="lower right", fontsize=6,
              title="Direction key", title_fontsize=7)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Saved direction heatmap: {out_path}")


def make_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in PANEL:
        sub = df[df["gene"] == gene]
        n_present    = int(sub["present"].sum()) if "present" in sub else 0
        n_concordant = int((sub["concordant"] == "yes").sum())
        n_discordant = int((sub["concordant"] == "no").sum())
        n_exp_none   = int((sub["concordant"] == "na").sum())
        n_scoreable  = n_concordant + n_discordant
        concordance  = round(n_concordant / n_scoreable, 2) if n_scoreable > 0 else None

        rows.append(dict(
            gene=gene,
            signature=SIGNATURE.get(gene, ""),
            n_datasets_queried=len(sub),
            n_datasets_present=n_present,
            n_concordant=n_concordant,
            n_discordant=n_discordant,
            n_expected_none=n_exp_none,
            concordance_rate=concordance,
        ))

    return pd.DataFrame(rows).sort_values("concordance_rate", ascending=False, na_position="last")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default=str(OUT))
    parser.add_argument("--out-dir",     default=str(OUT))
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_dir     = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_all(results_dir)
    if df.empty:
        return

    # Save master table
    master_path = out_dir / "master_query.csv"
    df.to_csv(master_path, index=False)
    log.info(f"Master table: {master_path} ({len(df)} rows)")

    # Figures
    make_presence_heatmap(df, out_dir / "presence_heatmap.png")
    make_direction_heatmap(df, out_dir / "direction_heatmap.png")

    # Summary table
    summary = make_summary_table(df)
    summary_path = out_dir / "summary_table.csv"
    summary.to_csv(summary_path, index=False)

    print("\n" + "=" * 60)
    print("EXTERNAL VALIDATION SUMMARY")
    print("=" * 60)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
