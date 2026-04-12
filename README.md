# OP Neurodegeneration Proteomics

**Cross-disease proteomics validation of a chronic organophosphate exposure protein signature against Alzheimer's disease, Parkinson's disease, Multiple Sclerosis, and Frontotemporal Dementia.**

---

## Background

Chronic organophosphate (OP) pesticide exposure is associated with neurological dysfunction, yet the protein-level mechanisms linking OP exposure to neurodegeneration remain poorly characterised. This repository provides the analysis pipeline used to test whether a ten-protein panel identified in a primary OP exposure proteomics study shows directional concordance in independent post-mortem brain and cerebrospinal fluid (CSF) proteomics datasets from four neuroinflammatory diseases.

---

## The 10-Protein Panel

Derived from primary OP exposure proteomics (chronic vs. control). All proteins are fixed prior to external comparison.

| Gene | OP Chronic log₂FC | Direction |
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

---

## External Datasets

### Brain Proteomics — Expression Atlas E-PROT

| Accession | Disease | Tissue | n (case / control) |
|-----------|---------|--------|--------------------|
| E-PROT-31 | Alzheimer's disease | Brain multi-region atlas | 71 / 32 |
| E-PROT-53 | Alzheimer's disease | DLPFC | — |
| E-PROT-56 | Alzheimer's disease | Temporal cortex | — |
| E-PROT-57 | Alzheimer's disease | Brain (Mount Sinai) | — |
| E-PROT-61 | Alzheimer's disease | Brain multi-region | — |
| E-PROT-32 | Alzheimer's disease (Braak staging) | Brain tau regions | 71 / 32 |
| E-PROT-39 | Alzheimer's disease | Brain (differential) | — |
| E-PROT-65 | Parkinson's disease | Prefrontal cortex | — |
| E-PROT-137 | FTD (GRN/MAPT) | Brain frontal/temporal | — |
| E-PROT-147 | FTD | Brain cell-type resolved | — |
| E-PROT-148 | FTD (semantic dementia) | Dentate gyrus | — |

### CSF Proteomics — PRIDE

| Accession | Disease | Method | n (case / control) |
|-----------|---------|--------|--------------------|
| PXD016278 | Alzheimer's disease | Label-free LC-MS/MS | 88 / 109 |
| PXD026491 | Parkinson's disease | Spectronaut DIA | ~96 / ~130 |
| PXD045058 | Multiple Sclerosis | DIA-NN Astral | 265 / 336 |

---

## Repository Structure

```
op-neurodegeneration-proteomics/
├── config/
│   ├── targets.yaml                  # Panel gene definitions & OP directions
│   ├── cohorts.csv                   # Cohort metadata
│   ├── pride_overrides.csv           # Manual PRIDE metadata corrections
│   └── pride_sample_maps/            # Sample group assignments per PXD accession
├── src/
│   ├── disease_labels.py             # Shared disease/control keyword library
│   ├── query_eprot.py                # Expression Atlas E-PROT brain proteomics parser
│   ├── query_pride.py                # PRIDE CSF proteomics querying
│   ├── query_hpa.py                  # Human Protein Atlas protein expression
│   ├── query_proteomicsdb.py         # ProteomicsDB / UniProt fallback
│   ├── query_opentargets.py          # OpenTargets disease evidence
│   ├── grand_finale_figure.py        # Main validation figure (OP vs disease proteomics)
│   ├── cross_modal_figure.py         # Cross-modal proteomics heatmap
│   ├── concordance_by_tissue.py      # Concordance stratified by compartment
│   ├── analyze_meta.py               # Random-effects meta-analysis
│   ├── guardrail_check.py            # Automated QA gate
│   ├── build_manuscript_package.R    # Manuscript figure & table assembly
│   ├── manuscript_consistency_check.R# Reproducibility QA
│   ├── mature_grand_finale_ggplot.R  # Publication ggplot figure
│   ├── null_benchmark_random_panels.R# Permutation null benchmark
│   └── publication_overlap_analysis.R# Tier-stratified overlap analysis
├── log fold change 10 Markers.xlsx   # Source-of-truth OP fold changes
├── requirements.txt
└── README.md
```

---

## Setup

### Python

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.11. Key dependencies: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `requests`, `pyyaml`, `openpyxl`.

### R

```r
install.packages(c("ggplot2", "metafor", "dplyr", "readr", "pheatmap", "officer"))
```

---

## Running the Pipeline

**Step 1 — E-PROT brain proteomics (AD / PD / FTD)**
```bash
python src/query_eprot.py
# → results/eprot_query.csv  (110 rows across 11 experiments)
```

**Step 2 — HPA plasma protein reference**
```bash
python src/query_hpa.py
# → results/hpa_query.csv
```

**Step 3 — Grand finale figure**
```bash
python src/grand_finale_figure.py
# → results/figures/grand_finale_op_neurodegeneration.png
```

**Step 4 — Cross-modal heatmap**
```bash
python src/cross_modal_figure.py
# → results/figures/cross_modal_panel_heatmap.png
```

**Step 5 — Meta-analysis & manuscript package (R)**
```bash
Rscript src/publication_overlap_analysis.R
Rscript src/build_manuscript_package.R
```

**Step 6 — QA gate**
```bash
python src/guardrail_check.py
Rscript src/manuscript_consistency_check.R
```

### PRIDE CSF data (PXD045058 — MS CSF, 2,828 proteins)

The DIA-NN gene matrix (12 MB) is not versioned here due to file size. To regenerate:

1. Connect to PRIDE FTP: `ftp.pride.ebi.ac.uk`
2. Path: `/pride/data/archive/2025/12/PXD045058/PILOTstudy_searchresults_sampleannotations.zip`
3. Extract `unique_genes_matrix.tsv` at byte offset `3,370,047,009` (compressed size `4,201,403` bytes)
4. Extract sample annotations at byte offset `5,271,864,849`
5. Place in `results/pride_discovery/`

---

## Key Results

Across 6 independent disease datasets (3 brain, 3 CSF):

| Finding | Count |
|---------|-------|
| Panel proteins concordant with OP direction | **10** |
| Panel proteins discordant | **7** |
| Absent or flat | **43** |

**Strongest concordant signals:**

- **ITIH3** −0.86 in MS CSF (PXD045058) ●
- **LRG1** −0.60 in MS CSF; −0.30 in AD CSF ●
- **GPX3** −0.70 in PD brain (E-PROT-65) ●
- **MYH7B** −0.60 in AD brain (E-PROT-31) ●
- **CFHR2** −0.33 in MS CSF (PXD045058) ●

---

## Claim Framing

This analysis demonstrates **cross-disease proteomic concordance** between the OP exposure signature and independent neuroinflammatory disease datasets. It is not a claim of definitive external validation in independent OP-exposed human cohorts. The appropriate framing is:

> Proteins responsive to chronic organophosphate exposure show partial directional concordance with abundance changes in AD, PD, and MS brain and CSF proteomics, consistent with shared neuroinflammatory mechanisms.

---

## Data Sources

| Resource | URL |
|----------|-----|
| Expression Atlas E-PROT | https://www.ebi.ac.uk/gxa/experiments |
| PRIDE Archive | https://www.ebi.ac.uk/pride |
| Human Protein Atlas | https://www.proteinatlas.org |
| ProteomicsDB | https://www.proteomicsdb.org |
| OpenTargets | https://www.opentargets.org |

---

## License

Code: MIT. Data accessed from public repositories under their respective terms of use (EMBL-EBI, PRIDE, HPA).
