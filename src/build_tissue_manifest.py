"""
build_tissue_manifest.py
------------------------
Create a cohort-level tissue manifest for cross-tissue / cross-sample analysis.

Outputs:
  - results/qa/tissue_group_manifest.csv
  - results/qa/tissue_group_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config"
RES = ROOT / "results"
QA = RES / "qa"


def normalize_tissue(biospecimen: str) -> tuple[str, str, str]:
    """
    Returns:
      tissue_supergroup: peripheral | central | other
      tissue_group:      blood | brain | csf | other
      tissue_subgroup:   blood_whole | blood_pbmc | blood_leukocyte | blood_t_cell |
                         blood_plasma | blood_serum | brain_region | brain_unspecified |
                         csf | other
    """
    b = str(biospecimen).strip().lower()

    if "cerebrospinal" in b or "csf" in b:
        return ("central", "csf", "csf")

    if "brain" in b:
        if "substantia_nigra" in b:
            return ("central", "brain", "brain_substantia_nigra")
        if "_" in b:
            # Keep explicit brain region labels when provided.
            return ("central", "brain", b)
        return ("central", "brain", "brain_unspecified")

    if "blood" in b or "pbmc" in b or "leukocyte" in b or "t_cell" in b or "plasma" in b or "serum" in b:
        if "whole_blood" in b:
            return ("peripheral", "blood", "blood_whole")
        if "pbmc" in b:
            return ("peripheral", "blood", "blood_pbmc")
        if "leukocyte" in b:
            return ("peripheral", "blood", "blood_leukocyte")
        if "t_cell" in b:
            return ("peripheral", "blood", "blood_t_cell")
        if "plasma" in b:
            return ("peripheral", "blood", "blood_plasma")
        if "serum" in b:
            return ("peripheral", "blood", "blood_serum")
        return ("peripheral", "blood", "blood_other")

    return ("other", "other", "other")


def build_manifest(
    cohorts_csv: Path,
    geo_query_csv: Path,
    out_manifest_csv: Path,
    out_summary_csv: Path,
) -> None:
    cohorts = pd.read_csv(cohorts_csv)
    geo = pd.read_csv(geo_query_csv)

    # Cohort-level counts from query output.
    geo_counts = (
        geo.groupby("accession", as_index=False)[["n_case", "n_control"]]
        .max()
        .rename(columns={"n_case": "observed_n_case", "n_control": "observed_n_control"})
    )

    # Keep only open GEO transcriptomics cohorts, then attach observed counts.
    m = cohorts[
        (cohorts["repository"] == "GEO")
        & (cohorts["access_type"] == "open")
        & (cohorts["modality"] == "transcriptomics")
    ].copy()
    m = m.merge(geo_counts, on="accession", how="left")

    # Mark whether cohort is numerically usable for case-control analysis.
    m["observed_n_case"] = pd.to_numeric(m["observed_n_case"], errors="coerce").fillna(0).astype(int)
    m["observed_n_control"] = pd.to_numeric(m["observed_n_control"], errors="coerce").fillna(0).astype(int)
    m["usable_case_control"] = (m["observed_n_case"] > 0) & (m["observed_n_control"] > 0)
    m["total_n"] = m["observed_n_case"] + m["observed_n_control"]

    # Tissue harmonization columns.
    tissue_cols = m["biospecimen"].apply(normalize_tissue).apply(pd.Series)
    tissue_cols.columns = ["tissue_supergroup", "tissue_group", "tissue_subgroup"]
    m = pd.concat([m, tissue_cols], axis=1)

    # Cross-tissue analysis strata.
    m["cross_tissue_stratum"] = m["disease"].astype(str) + "__" + m["tissue_group"].astype(str)
    m["cross_sample_stratum"] = m["disease"].astype(str) + "__" + m["tissue_subgroup"].astype(str)

    out_manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(out_manifest_csv, index=False)

    summary = (
        m.groupby(["disease", "tissue_supergroup", "tissue_group", "tissue_subgroup"], as_index=False)
        .agg(
            n_cohorts=("accession", "nunique"),
            n_usable_cohorts=("usable_case_control", "sum"),
            total_case=("observed_n_case", "sum"),
            total_control=("observed_n_control", "sum"),
        )
        .sort_values(["disease", "tissue_supergroup", "tissue_group", "tissue_subgroup"])
    )
    summary.to_csv(out_summary_csv, index=False)

    print(f"Wrote: {out_manifest_csv}")
    print(f"Wrote: {out_summary_csv}")
    print("\nTissue summary:")
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohorts", default=str(CFG / "cohorts.csv"))
    parser.add_argument("--geo-query", default=str(RES / "geo_query.csv"))
    parser.add_argument("--out-manifest", default=str(QA / "tissue_group_manifest.csv"))
    parser.add_argument("--out-summary", default=str(QA / "tissue_group_summary.csv"))
    args = parser.parse_args()

    build_manifest(
        cohorts_csv=Path(args.cohorts),
        geo_query_csv=Path(args.geo_query),
        out_manifest_csv=Path(args.out_manifest),
        out_summary_csv=Path(args.out_summary),
    )


if __name__ == "__main__":
    main()

