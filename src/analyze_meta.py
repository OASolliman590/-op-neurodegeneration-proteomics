"""
analyze_meta.py
---------------
Post-query analysis and meta summary for the fixed 10-gene panel.

Inputs:
  - results/master_query.csv (preferred), or all results/*_query.csv

Outputs:
  - results/analysis/cohort_effects.csv
  - results/analysis/disease_meta.csv
  - results/analysis/cross_disease_meta.csv
  - results/analysis/weighted_concordance_heatmap.png

Notes:
  - This is a transparent directional meta layer designed for heterogeneous
    open datasets where full variance terms are not always available.
  - Scoreable rows are those with concordant in {"yes", "no"}.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
ANALYSIS = RESULTS / "analysis"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]


def load_master(results_dir: Path) -> pd.DataFrame:
    master = results_dir / "master_query.csv"
    if master.exists():
        return pd.read_csv(master)

    frames = []
    for p in sorted(results_dir.glob("*_query.csv")):
        if p.name == "master_query.csv":
            continue
        try:
            df = pd.read_csv(p)
            df["source_file"] = p.name
            frames.append(df)
        except Exception:
            continue
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def safe_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def directional_effect(mean_case: float, mean_ctrl: float, threshold: float = 0.1) -> str:
    if np.isnan(mean_case) or np.isnan(mean_ctrl):
        return "no_data"
    d = mean_case - mean_ctrl
    if d < -threshold:
        return "down"
    if d > threshold:
        return "up"
    return "flat"


def binom_two_sided_p(k: int, n: int, p0: float = 0.5) -> float | None:
    if n <= 0:
        return None
    # exact two-sided p-value (small n here, so direct combinatorics is fine)
    probs = [math.comb(n, i) * (p0**i) * ((1 - p0) ** (n - i)) for i in range(n + 1)]
    p_obs = probs[k]
    p = sum(pi for pi in probs if pi <= p_obs + 1e-15)
    return min(1.0, float(p))


def weighted_mean(values: pd.Series, weights: pd.Series) -> float | None:
    m = values.notna() & weights.notna() & (weights > 0)
    if not m.any():
        return None
    v = values[m].astype(float)
    w = weights[m].astype(float)
    return float((v * w).sum() / w.sum())


def build_cohort_effects(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["n_case"] = safe_num(out.get("n_case"))
    out["n_control"] = safe_num(out.get("n_control"))
    out["mean_case"] = safe_num(out.get("mean_case"))
    out["mean_ctrl"] = safe_num(out.get("mean_ctrl"))
    out["effect_raw"] = out["mean_case"] - out["mean_ctrl"]
    # Cross-cohort robust effect proxy:
    # - use log2 fold-change when both means are non-negative
    # - fallback to raw difference otherwise (already log/intensity-like scales)
    out["effect_meta"] = np.where(
        (out["mean_case"] >= 0) & (out["mean_ctrl"] >= 0),
        np.log2((out["mean_case"] + 1.0) / (out["mean_ctrl"] + 1.0)),
        out["effect_raw"],
    )
    out["effect_dir"] = [
        directional_effect(mc, mk) for mc, mk in zip(out["mean_case"], out["mean_ctrl"])
    ]
    out["scoreable"] = out.get("concordant", pd.Series(index=out.index, dtype=object)).isin(["yes", "no"])
    out["is_concordant"] = out.get("concordant", pd.Series(index=out.index, dtype=object)).eq("yes")
    # Conservative weight proxy when variance is unavailable.
    out["weight_proxy"] = (out["n_case"].fillna(0) * out["n_control"].fillna(0)) / (
        out["n_case"].fillna(0) + out["n_control"].fillna(0) + 1e-9
    )
    return out


def disease_meta(cohort: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scoreable = cohort[cohort["scoreable"] == True].copy()
    for disease in sorted(scoreable["disease"].dropna().unique()):
        dsub = scoreable[scoreable["disease"] == disease]
        for gene in PANEL:
            g = dsub[dsub["gene"] == gene]
            n = len(g)
            k = int((g["is_concordant"] == True).sum())
            no = int((g["is_concordant"] == False).sum())
            p = binom_two_sided_p(k, n) if n > 0 else None
            w_yes = float(g.loc[g["is_concordant"] == True, "weight_proxy"].sum())
            w_all = float(g["weight_proxy"].sum())
            w_rate = (w_yes / w_all) if w_all > 0 else None
            eff = weighted_mean(g["effect_meta"], g["weight_proxy"])
            eff_dir = directional_effect(eff, 0.0) if eff is not None else "no_data"

            rows.append(dict(
                disease=disease,
                gene=gene,
                n_scoreable=n,
                n_concordant=k,
                n_discordant=no,
                concordance_rate=(k / n) if n > 0 else None,
                weighted_concordance=w_rate,
                binom_p_two_sided=p,
                weighted_effect=eff,
                weighted_effect_direction=eff_dir,
            ))
    return pd.DataFrame(rows)


def cross_disease_meta(dm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in PANEL:
        g = dm[dm["gene"] == gene].copy()
        n = int(g["n_scoreable"].fillna(0).sum())
        k = int(g["n_concordant"].fillna(0).sum())
        no = int(g["n_discordant"].fillna(0).sum())
        p = binom_two_sided_p(k, n) if n > 0 else None

        # Weight diseases by their scoreable study counts.
        w = g["n_scoreable"].fillna(0).astype(float)
        wc = weighted_mean(g["weighted_concordance"], w)
        we = weighted_mean(g["weighted_effect"], w)
        eff_dir = directional_effect(we, 0.0) if we is not None else "no_data"

        rows.append(dict(
            gene=gene,
            n_scoreable=n,
            n_concordant=k,
            n_discordant=no,
            concordance_rate=(k / n) if n > 0 else None,
            weighted_concordance=wc,
            binom_p_two_sided=p,
            weighted_effect=we,
            weighted_effect_direction=eff_dir,
        ))
    return pd.DataFrame(rows)


def plot_weighted_concordance(dm: pd.DataFrame, out_png: Path) -> None:
    if dm.empty:
        return
    pivot = dm.pivot(index="gene", columns="disease", values="weighted_concordance")
    pivot = pivot.reindex(PANEL)
    arr = pivot.values.astype(float)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    im = ax.imshow(np.nan_to_num(arr, nan=-0.01), aspect="auto", vmin=0, vmax=1, cmap="RdYlGn")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Weighted Concordance by Disease and Gene")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Weighted concordance")

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            txt = "NA" if np.isnan(v) else f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="black" if np.isnan(v) or v < 0.7 else "white")
    plt.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default=str(RESULTS))
    parser.add_argument("--out-dir", default=str(ANALYSIS))
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    master = load_master(results_dir)
    if master.empty:
        print("No data found in results directory.")
        return

    cohort = build_cohort_effects(master)
    cohort.to_csv(out_dir / "cohort_effects.csv", index=False)

    dm = disease_meta(cohort)
    dm.to_csv(out_dir / "disease_meta.csv", index=False)

    cm = cross_disease_meta(dm)
    cm.to_csv(out_dir / "cross_disease_meta.csv", index=False)

    plot_weighted_concordance(dm, out_dir / "weighted_concordance_heatmap.png")

    print("\nDisease meta (top rows):")
    print(dm.sort_values(["disease", "weighted_concordance"], ascending=[True, False])
          .head(20).to_string(index=False))
    print("\nCross-disease meta:")
    print(cm.sort_values("weighted_concordance", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
