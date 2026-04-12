# Proteomics External Validation of Organophosphate Exposure Biomarkers in Neuroinflammatory Disease

## Overview

This repository contains the proteomics external validation pipeline for a 10-protein panel identified under chronic organophosphate (OP) exposure. The analysis cross-references these OP-responsive proteins against independent brain and CSF proteomics datasets from Alzheimer's disease (AD), Parkinson's disease (PD), Multiple Sclerosis (MS), and Frontotemporal Dementia (FTD).

**Central hypothesis:** Proteins dysregulated by chronic OP exposure show overlapping abundance changes in neuroinflammatory disease proteomics, supporting shared molecular mechanisms between pesticide exposure and neurodegeneration.

---

## The 10-Protein OP Panel

| Gene | OP Chronic log2FC | Direction |
|------|-------------------|-----------|
| ACTG1 | −5.98 | Down |
| DNAH9 | +5.84 | Up |
| GPX3 | −2.83 | Down |
| VWF | +4.28 | Up |
| C4B | −4.15 | Down |
| CD44 | +1.58 | Up |
| CFHR2 | +1.47 | Up |
| ITIH3 | −1.21 | Down |
| LRG1 | −1.72 | Down |
| MYH7B | −5.32 | Down |

Source: `log fold change 10 Markers.xlsx` — original OP exposure proteomics study.

---

## Datasets

### Brain Proteomics — Expression Atlas E-PROT

| Accession | Disease | Biospecimen | Panel Coverage |
|-----------|---------|-------------|----------------|
| E-PROT-31 | Alzheimer's | Brain (multi-region atlas) | 7/10 |
| E-PROT-53 | Alzheimer's | Brain DLPFC | 2/10 |
| E-PROT-56 | Alzheimer's | Brain temporal | 3/10 |
| E-PROT-57 | Alzheimer's | Brain MtSinai | 5/10 |
| E-PROT-61 | Alzheimer's | Brain multi-region | 7/10 |
| E-PROT-32 | Alzheimer's (Braak stages) | Brain tau staging | 1/10 |
| E-PROT-65 | Parkinson's | Brain prefrontal | 5/10 |
| E-PROT-137 | FTD | Brain GRN/MAPT | 1/10 |
| E-PROT-147 | FTD | Brain celltype | 1/10 |
| E-PROT-148 | FTD | Brain dentate gyrus | 2/10 |
| E-PROT-39 | Alzheimer's | Brain (differential) | — |

### CSF Proteomics — PRIDE

| Accession | Disease | n case | n control | Panel Coverage | Method |
|-----------|---------|--------|-----------|----------------|--------|
| PXD016278 | Alzheimer's | 88 | 109 | 7/10 | Label-free MS |
| PXD026491 | Parkinson's | ~96 | ~130 | 8/10 | Spectronaut DIA |
| PXD045058 | Multiple Sclerosis | 265 | 336 | 7/10 | DIA-NN Astral |

> **Note:** PXD045058 (MS CSF) gene matrix (12.3 MB) and protein intensities (28 MB) are excluded from the repository due to file size. Re-extract from PRIDE FTP using `src/query_pride.py` or follow instructions in `results/pride_discovery/`.

### Protein Reference — Human Protein Atlas

- Plasma MS concentration (ng/mL) for normal human plasma
- Used as detectability baseline for CSF/brain proteomics layers

---

## Key Results

### Grand Finale Figure

`results/figures/grand_finale_op_neurodegeneration.png`

Shows the 10-protein OP signature (log2FC ±7) alongside 6 independent disease proteomics datasets (log2FC ±3, brain and CSF). Concordance markers (●) indicate proteins changing in the same direction under OP exposure and in disease.

**Summary across 6 datasets:**
- Concordant detections: **10**
- Discordant detections: **7**
- Absent / flat: **43**

**Strongest concordant signals:**
- ITIH3 −0.86● in MS CSF (PXD045058)
- LRG1 −0.60● in MS CSF; −0.30● in AD CSF
- GPX3 −0.70● in PD brain (E-PROT-65)
- MYH7B −0.60● in AD brain (E-PROT-31)
- CFHR2 −0.33● in MS CSF

### Publication-tier Meta-analysis

Stratified by compartment (brain / CSF) and disease, with uncertainty quantification (SE proxy, 95% CI, BH-FDR). See `results/analysis/publication_effects_with_uncertainty.csv`.

Dataset tiers:
- **Tier A** — Reported case/control metadata + case-control design (primary)
- **Tier B** — Metadata inferred from file headers (sensitivity)
- **Tier C** — Non-case-control or missing critical metadata (excluded from primary)

---

## Repository Structure

```
op_external_validation/
├── config/
│   ├── targets.yaml                  # Panel gene definitions & OP directions
│   ├── cohorts.csv                   # Dataset cohort metadata
│   └── pride_sample_maps/            # Manual sample group assignments
├── src/
│   ├── disease_labels.py             # Shared disease/control keyword library
│   ├── query_eprot.py                # Expression Atlas E-PROT parser
│   ├── query_pride.py                # PRIDE proteomics querying
│   ├── query_hpa.py                  # Human Protein Atlas querying
│   ├── query_proteomicsdb.py         # ProteomicsDB / UniProt fallback
│   ├── query_opentargets.py          # OpenTargets disease evidence
│   ├── blood_validate_primary10.py   # Blood proteomics validation
│   ├── cross_modal_figure.py         # Cross-modal proteomics heatmap
│   ├── grand_finale_figure.py        # Grand finale visualization
│   ├── concordance_by_tissue.py      # Concordance stratified by tissue
│   ├── build_manuscript_package.R    # Manuscript assembly
│   ├── manuscript_consistency_check.R# QA gate & consistency validation
│   ├── mature_grand_finale_ggplot.R  # ggplot figure maturation
│   ├── null_benchmark_random_panels.R# Null benchmark (random panels)
│   └── guardrail_check.py            # Automated QA guardrails
├── results/
│   ├── eprot_query.csv               # E-PROT results (110 rows)
│   ├── pride_query.csv               # PRIDE query results
│   ├── hpa_query.csv                 # HPA protein expression results
│   ├── proteomicsdb_query.csv        # ProteomicsDB / UniProt results
│   ├── pride_discovery/              # PRIDE CSF proteomics parsed data
│   ├── blood_validation/             # Blood proteomics meta-analysis
│   ├── analysis/                     # Publication-ready statistics & manifests
│   ├── figures/                      # All figures (PNG/PDF)
│   └── qa/                           # QA audit reports
├── manuscript/
│   └── current/                      # Active manuscript package
│       ├── figures/                  # Publication figures
│       ├── tables/                   # Data tables for submission
│       ├── text/                     # Methods & claims documentation
│       └── qa/                       # Manuscript QA reports
├── spec_kit/
│   ├── gap_register.csv              # Dataset status tracker
│   ├── guardrails.md                 # QA requirements
│   └── runbook.md                    # Operational procedures
├── log fold change 10 Markers.xlsx   # Source-of-truth OP fold changes
├── VALIDATION_UPDATES.md             # Full session change log
├── requirements.txt                  # Python dependencies
└── .gitignore
```

---

## Setup & Usage

### Requirements

```bash
pip install -r requirements.txt
```

Key dependencies: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `requests`, `pyyaml`, `GEOparse`, `openpyxl`

R dependencies (for manuscript figures): `ggplot2`, `metafor`, `dplyr`, `readr`, `pheatmap`

### Running the pipeline

**Step 1 — E-PROT brain proteomics (AD/PD/FTD):**
```bash
python src/query_eprot.py
# Output: results/eprot_query.csv (110 rows)
```

**Step 2 — HPA plasma protein reference:**
```bash
python src/query_hpa.py
# Output: results/hpa_query.csv
```

**Step 3 — Grand finale visualization:**
```bash
python src/grand_finale_figure.py
# Output: results/figures/grand_finale_op_neurodegeneration.png
```

**Step 4 — Cross-modal heatmap:**
```bash
python src/cross_modal_figure.py
# Output: results/figures/cross_modal_panel_heatmap.png
```

**Step 5 — Manuscript package (R):**
```bash
Rscript src/build_manuscript_package.R
# Output: manuscript/current/
```

### PRIDE CSF data (PXD045058 MS CSF)

The gene matrix is too large for GitHub. To regenerate:
1. The data is at `ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD045058/`
2. Use byte-offset extraction from `PILOTstudy_searchresults_sampleannotations.zip`
3. See extraction offsets documented in `spec_kit/gap_register.csv`

---

## Claim Framing

This analysis demonstrates **cross-disease proteomic overlap** between the OP exposure signature and neuroinflammatory disease datasets. This is **not** a claim of definitive external validation of OP biomarker performance in independent OP-exposed human cohorts. The appropriate framing is:

> Proteins responsive to chronic organophosphate exposure show partial directional concordance with abundance changes observed in AD, PD, and MS brain and CSF proteomics, consistent with shared neuroinflammatory mechanisms.

See `results/analysis/scientific_validity_appraisal.md` and `manuscript/current/text/publication_methods_and_claims.md` for full details.

---

## Citation & Data Sources

- Expression Atlas E-PROT: [ebi.ac.uk/gxa](https://www.ebi.ac.uk/gxa)
- PRIDE proteomics: [ebi.ac.uk/pride](https://www.ebi.ac.uk/pride)
- Human Protein Atlas: [proteinatlas.org](https://www.proteinatlas.org)
- PXD016278: AD CSF proteomics (88 vs 109)
- PXD026491: PD CSF Spectronaut DIA (Hansson et al. 2021, Brain)
- PXD045058: MS CSF large-scale DIA-NN PILOT (265 MS vs 336 nonMS neurological controls)
