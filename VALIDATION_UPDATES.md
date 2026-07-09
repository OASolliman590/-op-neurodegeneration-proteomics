# Proteomics External Validation — Session Updates

**Project:** OP Chronic Exposure Signature vs Neuroinflammatory Disease Proteomics
**Panel:** ACTG1, DNAH9, GPX3, VWF, C4B, CD44, CFHR2, ITIH3, LRG1, MYH7B
**Date:** 2026-04-11

---

## July 2026 v2 Disease-Context Expansion

The repository has been updated from a neurodegeneration-only external comparison to a broader human disease-context validation analysis.

New v2 outputs are saved in `results/marker10_integrated_disease_context_v2/`.

### v2 scope

| Metric | Count |
| --- | ---: |
| Human disease contrasts | 28 |
| Marker-concordance cells | 280 |
| Disease-specific pathway/proximity maps | 28 |
| STRING/database role terms used | 19 |

### Disease ontology coverage

- Nervous system: Alzheimer disease, Parkinson disease, multiple sclerosis.
- Metabolic/endocrine: diabetes mellitus and metabolic dysfunction.
- Kidney/urinary: CKD, diabetic kidney disease/prognosis, AKI subphenotype, kidney dysfunction, proximal tubule tissue context.
- Hepatic: liver dysfunction.
- Cardiovascular/inflammatory: cardiac disease and MIS-C/cardiac involvement.

### Claim boundary

The v2 analysis supports external disease-context plausibility for the fixed OP marker panel. It does not by itself prove OP-induced organ injury, clinical diagnostic utility, or causal mechanism in the primary OP cohort.

---

## Overview

Goal: Validate 10 protein markers identified under chronic organophosphate (OP) exposure against independent proteomics datasets from neuroinflammatory diseases (Alzheimer's, Parkinson's, Multiple Sclerosis, FTD).

---

## Updates This Session

### 1. FTD E-PROT Experiments Added (`src/query_eprot.py`, `src/disease_labels.py`)

- Added **FTD** to `DISEASE_KEYWORDS` in `disease_labels.py`
- Added 4 new experiments to `BASELINE_EXPERIMENTS`:

| Accession | Disease | Biospecimen | Panel Hits |
|-----------|---------|-------------|------------|
| E-PROT-32 | Alzheimers (reclassified from FTD) | brain_tau_braak_stages | VWF −1.29● |
| E-PROT-137 | FTD | brain_GRN_MAPT | VWF flat |
| E-PROT-147 | FTD | brain_celltype | VWF −0.40● |
| E-PROT-148 | FTD | brain_dentate_gyrus | LRG1 up○, CD44 up |

- Fixed `extract_gene_name_col()` to prefer gene symbol columns over Ensembl ID columns (was causing 0 hits in FTD experiments)
- **eprot_query.csv**: 70 → 110 rows, 38 present

---

### 2. MS CSF Proteomics — PXD045058 (NEW)

**Dataset:** Large-scale MS CSF DIA-NN proteomics (PILOT study)
**FTP path:** `/pride/data/archive/2025/12/PXD045058/`
**Method:** Selective byte-offset extraction from 5.27 GB ZIP (no full download needed)

**Files extracted:**
- `PXD045058_unique_genes_matrix.tsv` — 12.3 MB, 2,828 proteins × 645 samples
- `PXD045058_sample_annotations.xlsx` — 115 KB, full sample metadata
- `PXD045058_protein_intensities.tsv` — 28.1 MB, directLFQ intensities

**Cohort:** 265 MS vs 336 non-MS neurological controls (CSF, DIA-NN quantification)

| Gene | log2FC | Direction | Concordant (vs targets) |
|------|--------|-----------|------------------------|
| CFHR2 | −0.33 | down | yes ● |
| ITIH3 | −0.86 | down | yes ● |
| LRG1 | −0.60 | down | yes ● |
| VWF | −0.26 | flat | na |
| GPX3 | +0.04 | flat | na |
| C4B | −0.23 | flat | na |
| CD44 | +0.18 | flat | na |
| ACTG1 | — | not detected | na |
| DNAH9 | — | not detected | na |
| MYH7B | — | not detected | na |

**Result:** 7/10 detected, 3/3 concordant, 0 discordant
Appended to `results/pride_discovery/pride_ad_pd_ms_query.csv`

---

### 3. PXD064570 — Definitively Blocked

**Reason:** `Astral_proteome_search_results.zip` = **73.9 GB**; all other files are raw instrument ZIPs (25–116 GB each). No accessible protein-level quant file without institutional HPC.
**Status:** `blocked` in gap_register.csv — root_cause: `quant_zip_is_73GB_not_downloadable`

---

### 4. Gap Register Updated (`spec_kit/gap_register.csv`)

- PXD045058: added as `working` — 265 MS vs 336 nonMS, 3 concordant
- PXD064570: updated root_cause to `quant_zip_is_73GB_not_downloadable`

---

### 5. Grand Finale Figure (`src/grand_finale_figure.py`)

**Output:** `results/figures/grand_finale_op_neurodegeneration.png` (180 KB)

**Layout:**

| Section | Datasets | Description |
|---------|----------|-------------|
| OP Chronic Signal | Original OP study | log2FC ±7 scale, actual exposure data |
| Brain Proteomics | E-PROT-31 (AD atlas), E-PROT-61 (AD multi), E-PROT-65 (PD frontal) | Post-mortem brain, ±3 scale |
| CSF Proteomics | PXD016278 (AD CSF), PXD026491 (PD CSF), PXD045058 (MS CSF) | 88–336 cases vs controls |

**Concordance markers:** ● = same direction as OP, ○ = opposite direction
**Overall:** 10 concordant, 7 discordant, 43 absent/flat across 6 independent disease datasets

**Key concordant hits:**
- GPX3, ITIH3, LRG1 — down in PD brain (E-PROT-65) matching OP direction
- ITIH3 −0.9, LRG1 −0.6 — down in MS CSF (PXD045058)
- LRG1 −0.3 — down in AD CSF (PXD016278)
- MYH7B −0.6 — down in AD brain (E-PROT-31)

---

## File Inventory

| File | Description |
|------|-------------|
| `src/query_eprot.py` | E-PROT parser (updated with FTD experiments + gene col fix) |
| `src/disease_labels.py` | Shared disease keywords (FTD added) |
| `src/grand_finale_figure.py` | Grand finale visualization script |
| `src/cross_modal_figure.py` | Proteomics-only cross-modal heatmap |
| `results/eprot_query.csv` | 110 rows, 38 present |
| `results/pride_discovery/pride_ad_pd_ms_query.csv` | 30 rows (PXD026491 + PXD045058) |
| `results/pride_discovery/PXD045058_unique_genes_matrix.tsv` | Raw DIA-NN gene matrix |
| `results/pride_discovery/PXD045058_sample_annotations.xlsx` | Sample metadata (265 MS / 336 nonMS) |
| `results/figures/grand_finale_op_neurodegeneration.png` | **Grand finale figure** |
| `results/figures/cross_modal_panel_heatmap.png` | Proteomics cross-modal heatmap |
| `spec_kit/gap_register.csv` | Updated status for all datasets |
| `spec_kit/codex_cycle.md` | 6-task Codex cycle (TASK-01 through TASK-06) |

---

## Datasets Used in Final Analysis

| Accession | Disease | Biospecimen | n_case | n_ctrl | Panel hits | Source |
|-----------|---------|-------------|--------|--------|------------|--------|
| E-PROT-31 | AD | Brain (atlas) | 71 | 32 | 7/10 | Expression Atlas |
| E-PROT-53 | AD | Brain DLPFC | — | — | 2/10 | Expression Atlas |
| E-PROT-56 | AD | Brain temporal | — | — | 3/10 | Expression Atlas |
| E-PROT-57 | AD | Brain MtSinai | — | — | 5/10 | Expression Atlas |
| E-PROT-61 | AD | Brain multi | — | — | 7/10 | Expression Atlas |
| E-PROT-32 | AD (Braak) | Brain tau stages | 71 | 32 | 1/10 | Expression Atlas |
| E-PROT-65 | PD | Brain prefrontal | — | — | 5/10 | Expression Atlas |
| E-PROT-137 | FTD | Brain GRN/MAPT | — | — | 1/10 | Expression Atlas |
| E-PROT-147 | FTD | Brain celltype | — | — | 1/10 | Expression Atlas |
| E-PROT-148 | FTD | Brain dentate gyrus | — | — | 2/10 | Expression Atlas |
| PXD016278 | AD | CSF | 88 | 109 | 7/10 | PRIDE |
| PXD026491 | PD | CSF | — | — | 8/10 | PRIDE |
| PXD045058 | MS | CSF | 265 | 336 | 7/10 | PRIDE |

---

## Blocked Datasets

| Accession | Disease | Reason |
|-----------|---------|--------|
| PXD064570 | MS CSF | 73.9 GB quant ZIP — not downloadable |
| PXD034840 | MS CSF | Only .mzid.gz peptide files, no protein quant |
| PXD011216 | PD CSF | Incomplete submission |

---

## Narrative Summary

The 10-protein panel identified under chronic organophosphate (OP) exposure was compared against 13 independent proteomics datasets spanning AD, PD, MS, and FTD. Across brain and CSF proteomics:

- **ITIH3** and **LRG1** show consistent downregulation in PD brain and MS CSF, concordant with their downward trend under chronic OP exposure.
- **GPX3** shows downregulation in PD brain tissue (E-PROT-65, −0.7●), consistent with OP-induced oxidative stress signature.
- **MYH7B** is significantly reduced in AD brain (E-PROT-31, −0.6●) and under chronic OP exposure (−5.3).
- **CFHR2** is reduced in MS CSF (−0.33●), consistent with complement pathway dysregulation under OP exposure.

These signals indicate **partial cross-disease proteomic concordance** with the OP signature. This supports hypothesis generation around shared biology, but is not by itself definitive external validation of OP exposure biomarkers in independent OP-exposed human cohorts.

---

## Figure Maturation Pass (ggplot + validity appraisal)

### New script

- `src/mature_grand_finale_ggplot.R`

### New outputs

- `results/figures/grand_finale_op_neurodegeneration_ggplot.png`
- `results/figures/grand_finale_op_neurodegeneration_ggplot.pdf`
- `results/analysis/grand_finale_ggplot_data.csv`
- `results/analysis/grand_finale_dataset_validity.csv`
- `results/analysis/scientific_validity_appraisal.md`

### What improved

- **Scientific validity appraisal added** with explicit tiering of dataset comparability.
- **ggplot figure added** with:
  - full **gene names** on y-axis (`SYMBOL — Full Name`)
  - two panels: effect-size view + concordance view
  - dataset labels with sample size where available.
- **Direction-consistency check added**:
  - highlights OP direction conflicts between spreadsheet OP fold-change and `targets.yaml`.

---

## Publication Rigor Implementation (2026-04-11)

This section supersedes the earlier "validated" framing for manuscript-facing outputs.

### Implemented checklist

- OP direction source-of-truth fixed to spreadsheet:
  - `log fold change 10 Markers.xlsx` is now the sole source for expected OP direction.
  - Exported to `results/analysis/op_direction_source_of_truth.csv`.
- Uncertainty-aware per-dataset statistics added:
  - `SE_proxy = sqrt(1/n_case + 1/n_control)` with per-row 95% CI, p-value, and BH-FDR.
  - Exported in `results/analysis/publication_effects_with_uncertainty.csv`.
- Random-effects meta-analysis added:
  - Stratified by compartment (`brain`, `csf`) and by disease+compartment.
  - Model: REML + Hartung-Knapp (`metafor::rma.uni(..., test='knha')`).
- Brain and CSF kept stratified:
  - No pooled primary effect-size claims across compartments.
- Tier policy operationalized:
  - Tier A: reported case/control metadata + case-control design (primary).
  - Tier B: metadata inferred from machine-readable headers (sensitivity only).
  - Tier C: non-case-control or missing critical metadata (excluded from primary).
- Tier C sensitivity reporting added:
  - Primary outputs: Tier A only.
  - Sensitivity outputs: Tier A + Tier B, excluding Tier C.

### Tier-C repair outcome

- `PXD026491` now has header-derived inferred counts (`n_case=96`, `n_control=130`).
- It remains Tier B (not Tier A) because group-mapping uncertainty is retained.

### Publication-facing outputs generated

- `results/analysis/publication_dataset_manifest.csv`
- `results/analysis/publication_meta_primary_by_compartment.csv`
- `results/analysis/publication_meta_primary_by_disease_compartment.csv`
- `results/analysis/publication_meta_sensitivity_no_tier_c_by_compartment.csv`
- `results/analysis/publication_meta_sensitivity_no_tier_c_by_disease_compartment.csv`
- `results/figures/publication_main_tierA_heatmap.png`
- `results/figures/publication_sensitivity_no_tierC_heatmap.png`
- `results/analysis/publication_methods_and_claims.md`

### Claim framing now used

- Preferred claim: **cross-disease proteomic overlap / partial concordance** with OP signature.
- Not claimed: definitive external validation of OP biomarker performance in independent OP-exposed human cohorts.

---

## Manuscript Lockdown Build (2026-04-11)

### New build/QA scripts

- `src/build_manuscript_package.R`
- `src/manuscript_consistency_check.R`
- `src/null_benchmark_random_panels.R`

### What this build now guarantees

- Single manuscript source-of-truth: `results/analysis/publication_effects_with_uncertainty.csv`.
- Regenerated manuscript package in `manuscript/current/` with:
  - main figure (Tier A),
  - sensitivity figure (Tier A + Tier B),
  - supplementary two-panel ggplot figure,
  - null-benchmark figure/table,
  - methods + validity text,
  - QA reports.
- Automated consistency checks pass:
  - sample-size consistency,
  - concordance recomputation,
  - OP direction source lock,
  - cross-artifact count consistency,
  - effect-state taxonomy and rule logic,
  - language guardrails.
- Legacy conflicting figures moved to `results/figures/archive_legacy/`.
