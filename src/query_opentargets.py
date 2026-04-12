"""
query_opentargets.py
--------------------
Query Open Targets Platform for disease-specific evidence for panel genes.

Open Targets aggregates evidence from:
  - Differential expression (RNA-seq, microarray)
  - Proteomics (pQTL, protein abundance studies)
  - GWAS / genetic association
  - Literature / text-mining

For each gene × disease pair it returns an overall association score (0–1)
and a breakdown by evidence type including expression and proteomics.

Uses the Open Targets GraphQL API — free, no registration.

Output: results/opentargets_query.csv
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "results"
CFG  = ROOT / "config"

PANEL = [
    "ACTG1", "DNAH9", "GPX3", "VWF",
    "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B",
]
OUTPUT_COLUMNS = [
    "accession", "disease", "biospecimen", "modality", "exposure_type",
    "gene", "present", "n_case", "n_control", "mean_case", "mean_ctrl",
    "direction", "expected", "concordant",
    "ot_overall_score", "ot_expression_score", "ot_proteomics_score",
]

# Open Targets disease IDs (MONDO, as of OT Platform v4)
DISEASE_IDS = {
    "Alzheimers": "MONDO_0004975",   # Alzheimer disease
    "Parkinsons": "MONDO_0005180",   # Parkinson disease
    "MS":         "MONDO_0005301",   # multiple sclerosis
}

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"


def load_expected(targets_yaml: Path) -> dict:
    with open(targets_yaml) as f:
        cfg = yaml.safe_load(f)
    return {g: v for g, v in cfg["targets"].items()}


def expected_direction(gene: str, targets: dict) -> str:
    t = targets.get(gene, {})
    chronic = t.get("chronic_direction") or "none"
    return chronic if chronic != "none" else "none"


def get_ensembl_id(gene_symbol: str) -> str | None:
    """Resolve HGNC symbol to Ensembl gene ID via Open Targets search."""
    query = """
    query SearchGene($q: String!) {
      search(queryString: $q, entityNames: ["target"]) {
        hits {
          id
          entity
          object { ... on Target { approvedSymbol } }
        }
      }
    }
    """
    try:
        r = requests.post(OT_API,
                          json={"query": query, "variables": {"q": gene_symbol}},
                          timeout=30)
        r.raise_for_status()
        hits = r.json()["data"]["search"]["hits"]
        for hit in hits:
            sym = hit.get("object", {}).get("approvedSymbol", "")
            if sym.upper() == gene_symbol.upper():
                return hit["id"]
    except Exception as e:
        log.debug(f"  Ensembl lookup failed for {gene_symbol}: {e}")
    return None


def get_association(ensembl_id: str, disease_id: str) -> dict:
    """
    Fetch association score and evidence breakdown for a target-disease pair.
    Uses target.associatedDiseases and matches by disease ID client-side.
    Returns dict with overall score and per-datatype scores.
    """
    query = """
    query TargetAssoc($target: String!, $size: Int!) {
      target(ensemblId: $target) {
        associatedDiseases(
          enableIndirect: true
          page: { index: 0, size: $size }
        ) {
          rows {
            disease { id }
            score
            datatypeScores { id score }
          }
        }
      }
    }
    """
    try:
        # Try up to 500 to find the specific disease
        for page_size in [50, 200, 500]:
            r = requests.post(OT_API,
                              json={"query": query,
                                    "variables": {"target": ensembl_id,
                                                  "size": page_size}},
                              timeout=30)
            r.raise_for_status()
            rows = (r.json().get("data", {})
                    .get("target", {})
                    .get("associatedDiseases", {})
                    .get("rows", []))
            for row in rows:
                if row.get("disease", {}).get("id") == disease_id:
                    return {"score": row["score"],
                            "datatypeScores": row.get("datatypeScores", [])}
            if len(rows) < page_size:
                break  # No more results
        return {}
    except Exception as e:
        log.debug(f"  Association query failed: {e}")
        return {}


def get_expression_evidence(ensembl_id: str, disease_efo: str) -> list[dict]:
    """
    Fetch differential expression evidence for a target-disease pair.
    Returns list of {log2fc, pvalue, tissue, experiment, type} dicts.
    """
    query = """
    query ExpressionEvidence($target: String!, $disease: String!) {
      disease(efoId: $disease) {
        evidences(ensemblIds: [$target],
                  datasourceIds: ["expression_atlas", "gtex"]) {
          rows {
            score
            datasourceId
            literature
            studyOverview
            log2FoldChangeValue
            log2FoldChangePercentileRank
            resourceScore
            pValueExponent
            pValueMantissa
            biologicalModelAllelicComposition
          }
        }
      }
    }
    """
    try:
        r = requests.post(OT_API,
                          json={"query": query,
                                "variables": {"target": ensembl_id,
                                              "disease": disease_efo}},
                          timeout=30)
        r.raise_for_status()
        rows = (r.json().get("data", {})
                .get("disease", {})
                .get("evidences", {})
                .get("rows", []))
        results = []
        for row in rows:
            fc = row.get("log2FoldChangeValue")
            if fc is not None:
                exp_val = row.get("pValueMantissa")
                exp_exp = row.get("pValueExponent")
                pval = (float(exp_val) * 10 ** float(exp_exp)
                        if exp_val is not None and exp_exp is not None
                        else np.nan)
                results.append({
                    "log2fc": float(fc),
                    "pvalue": pval,
                    "datasource": row.get("datasourceId", ""),
                    "study": row.get("studyOverview", ""),
                })
        return results
    except Exception as e:
        log.debug(f"  Expression evidence query failed: {e}")
        return []


def get_proteomics_evidence(ensembl_id: str, disease_efo: str) -> list[dict]:
    """
    Fetch proteomics-specific evidence (pQTL, protein expression).
    """
    query = """
    query ProteomicsEvidence($target: String!, $disease: String!) {
      disease(efoId: $disease) {
        evidences(ensemblIds: [$target],
                  datasourceIds: ["impc", "ot_genetics_portal",
                                   "proteomics_hpa", "uniprot_variants"]) {
          rows {
            score
            datasourceId
            studyOverview
          }
        }
      }
    }
    """
    try:
        r = requests.post(OT_API,
                          json={"query": query,
                                "variables": {"target": ensembl_id,
                                              "disease": disease_efo}},
                          timeout=30)
        r.raise_for_status()
        rows = (r.json().get("data", {})
                .get("disease", {})
                .get("evidences", {})
                .get("rows", []))
        return [{"score": row.get("score", 0),
                 "datasource": row.get("datasourceId", ""),
                 "study": row.get("studyOverview", "")}
                for row in rows]
    except Exception as e:
        log.debug(f"  Proteomics evidence query failed: {e}")
        return []


def score_to_direction(ot_score: float) -> str:
    """
    Open Targets association score (0-1) cannot directly give direction.
    Use it as a presence/strength indicator only.
    A score > 0 means evidence exists; we return 'no_data' (present, no FC).
    """
    if ot_score > 0:
        return "no_data"
    return "absent"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default=str(CFG / "targets.yaml"))
    parser.add_argument("--out",     default=str(OUT / "opentargets_query.csv"))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    targets = load_expected(Path(args.targets))

    # Cache Ensembl IDs (one lookup per gene)
    log.info("Resolving Ensembl IDs for panel genes …")
    ensembl_map = {}
    for gene in PANEL:
        eid = get_ensembl_id(gene)
        ensembl_map[gene] = eid
        status = eid if eid else "NOT FOUND"
        log.info(f"  {gene}: {status}")
        time.sleep(0.3)

    all_records = []

    for disease, efo_id in DISEASE_IDS.items():
        log.info(f"\nQuerying Open Targets — {disease} ({efo_id})")

        for gene in PANEL:
            eid = ensembl_map.get(gene)
            exp_dir = expected_direction(gene, targets)

            if not eid:
                all_records.append(dict(
                    accession="OpenTargets_RNA",
                    disease=disease, biospecimen="multi",
                    modality="transcriptomics",
                    exposure_type="not_applicable",
                    gene=gene, present=False,
                    n_case=0, n_control=0,
                    mean_case=None, mean_ctrl=None,
                    direction="absent", expected=exp_dir, concordant="na",
                    ot_overall_score=None,
                ))
                continue

            # Overall association score
            assoc = get_association(eid, efo_id)
            ot_score = assoc.get("score", 0) or 0

            # Per-datatype scores
            dtype_scores = {d["id"]: d["score"]
                            for d in assoc.get("datatypeScores", [])}
            expr_score  = dtype_scores.get("rna_expression", 0) or 0
            prot_score  = dtype_scores.get("affected_pathway", 0) or 0
            genetic_score = dtype_scores.get("genetic_association", 0) or 0

            # Expression evidence with fold changes
            expr_evidence = get_expression_evidence(eid, efo_id)
            time.sleep(0.3)
            prot_evidence = get_proteomics_evidence(eid, efo_id)
            time.sleep(0.3)

            # ── RNA / expression row ──────────────────────────────────────
            fcs = [e["log2fc"] for e in expr_evidence if not np.isnan(e["log2fc"])]
            if fcs:
                mean_fc = float(np.mean(fcs))
                obs_dir = ("down" if mean_fc < -0.3
                           else "up" if mean_fc > 0.3 else "flat")
                concordant = ("yes" if obs_dir == exp_dir
                              else "no" if exp_dir != "none" else "na")
                present = True
                mc, mk = round(mean_fc, 4), 0.0
            else:
                obs_dir = score_to_direction(expr_score)
                concordant = "na"
                present = ot_score > 0
                mc = mk = None

            all_records.append(dict(
                accession="OpenTargets_RNA",
                disease=disease,
                biospecimen="multi",
                modality="transcriptomics",
                exposure_type="not_applicable",
                gene=gene,
                present=present,
                n_case=len(fcs),
                n_control=len(fcs),
                mean_case=mc,
                mean_ctrl=mk,
                direction=obs_dir,
                expected=exp_dir,
                concordant=concordant,
                ot_overall_score=round(ot_score, 4),
                ot_expression_score=round(expr_score, 4),
            ))

            # ── Proteomics row ────────────────────────────────────────────
            prot_present = len(prot_evidence) > 0 or prot_score > 0
            all_records.append(dict(
                accession="OpenTargets_PROT",
                disease=disease,
                biospecimen="multi",
                modality="proteomics",
                exposure_type="not_applicable",
                gene=gene,
                present=prot_present,
                n_case=len(prot_evidence),
                n_control=len(prot_evidence),
                mean_case=None,
                mean_ctrl=None,
                direction="no_data" if prot_present else "absent",
                expected=exp_dir,
                concordant="na",
                ot_overall_score=round(ot_score, 4),
                ot_proteomics_score=round(prot_score, 4),
            ))

            log.info(f"  {gene}: OT_score={ot_score:.3f}, "
                     f"expr={expr_score:.3f}, "
                     f"n_fc_evidence={len(fcs)}, "
                     f"direction={obs_dir}")
            time.sleep(0.3)

    df = pd.DataFrame(all_records, columns=OUTPUT_COLUMNS)
    df.to_csv(args.out, index=False)
    log.info(f"\nSaved {len(df)} rows to {args.out}")

    # Quick summary
    rna_df = df[(df["accession"] == "OpenTargets_RNA") &
                (df["concordant"].isin(["yes", "no"]))]
    if not rna_df.empty:
        print("\nOpen Targets — RNA concordance:")
        print(rna_df.groupby(["disease", "gene", "concordant"]).size()
              .unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
