# External Validation v2 Technical Note

## Purpose

This analysis documents the external-validation component for a fixed 10-protein chronic organophosphate exposure signature. The goal is to test whether the marker panel shows disease-context concordance or pathway/network proximity across independent public human proteomics cohorts.

The analysis is intentionally scoped as external disease-context support. It does not claim that non-OP disease cohorts independently validate OP exposure, organ injury, diagnostic utility, or in vivo DNA damage.

## Fixed Marker Set

The marker set was fixed before cross-disease querying:

| Marker | OP reference direction |
| --- | --- |
| ACTG1 | Down in chronic OP exposure |
| DNAH9 | Up in chronic OP exposure |
| GPX3 | Down in chronic OP exposure |
| VWF | Up in chronic OP exposure |
| C4B | Down in chronic OP exposure |
| CD44 | Up in chronic OP exposure |
| CFHR2 | Up in chronic OP exposure |
| ITIH3 | Down in chronic OP exposure |
| LRG1 | Down in chronic OP exposure |
| MYH7B | Down in chronic OP exposure |

## Cohort Eligibility Principles

Human disease cohorts were prioritized when they had:

- Human samples or human tissue context.
- Patient-level disease/comparator labels.
- Processed protein-level abundance data.
- A case/control, case/comparator, or disease-subphenotype contrast.
- Protein identifiers that could be harmonized to gene symbols and/or UniProt accessions.

Non-human models, cell lines, raw-only datasets, peptide-only matrices without protein rollup, PTM-only datasets, and opaque sample labels without a usable key are not treated as primary external-validation cohorts.

## Included Human Disease Contrasts

| Cohort | Domain | Disease/context | Accession | Case label | Control/comparator | n case | n control |
| --- | --- | --- | --- | --- | --- | ---: | ---: |
| C01 | Nervous system | Alzheimer disease | PXD021718 | case | control | 10 | 5 |
| C02 | Nervous system | Alzheimer disease | PXD016278 | case | control | 88 | 109 |
| C03 | Nervous system | Parkinson disease | PXD055996 | case | control | 40 | 40 |
| C04 | Nervous system | Multiple sclerosis | PXD050479 | case | control | 30 | 19 |
| C05 | Nervous system | Multiple sclerosis | PXD045058 | case | control | 267 | 321 |
| C06 | Nervous system | Parkinson disease | PXD038555 | case | control | 10 | 20 |
| C07 | Nervous system | Multiple sclerosis | PXD032287 | case | control | 5 | 5 |
| C08 | Metabolic/endocrine | T2DM vs PTDM | PXD054961 | T2DM | PTDM | 8 | 8 |
| C09 | Metabolic/endocrine | T2DM vs normoglycemic transplant recipients | PXD054961 | T2DM | NG | 8 | 8 |
| C10 | Metabolic/endocrine | Diabetes mellitus / metabolic dysfunction | PXD054961 | case | control | 16 | 8 |
| C11 | Metabolic/endocrine | Metabolic dysfunction | PAD000040 | case | control | 47 | 69 |
| C12 | Metabolic/endocrine | PTDM vs normoglycemic transplant recipients | PXD054961 | PTDM | NG | 8 | 8 |
| C13 | Kidney/urinary | CKD stage 5 vs control | PXD016433 | ckd_stage_5 | control | 10 | 6 |
| C14 | Kidney/urinary | All CKD vs control | PXD016433 | all_ckd | control | 30 | 6 |
| C15 | Kidney/urinary | CKD stage 3 vs control | PXD016433 | ckd_stage_3 | control | 11 | 6 |
| C16 | Kidney/urinary | CKD stage 1 vs control | PXD016433 | ckd_stage_1 | control | 9 | 6 |
| C17 | Kidney/urinary | T2D renal dysfunction/prognosis | PXD016571 | Bad | Good | 19 | 45 |
| C18 | Kidney/urinary | Proximal tubule tissue PTDM vs NG | PXD054937 | PTDM | NG | 6 | 5 |
| C19 | Kidney/urinary | Proximal tubule tissue T2DM vs NG | PXD054937 | T2DM | NG | 6 | 5 |
| C20 | Kidney/urinary | Proximal tubule tissue T2DM vs PTDM | PXD054937 | T2DM | PTDM | 6 | 6 |
| C21 | Kidney/urinary | AKI plasma subphenotype 2 vs 1 | PXD044264 | AKI_subphenotype_2 | AKI_subphenotype_1 | 42 | 14 |
| C22 | Kidney/urinary | Kidney dysfunction | PXD046550 | case | control | 56 | 18 |
| C23 | Hepatic | Liver dysfunction | PXD011839 | case | control | 48 | 24 |
| C24 | Hepatic | Liver dysfunction | GSE251855 | case | control | 37 | 20 |
| C25 | Cardiovascular/inflammatory | Cardiac disease | PXD009356 | case | control | 13 | 14 |
| C26 | Cardiovascular/inflammatory | Cardiac disease | GSE95368 | case | control | 27 | 6 |
| C27 | Cardiovascular/inflammatory | MIS-C / cardiac involvement | PXD029375 | case | control | 34 | 25 |
| C28 | Cardiovascular/inflammatory | Cardiac disease | GSE181091 | case | control | 8 | 4 |

PXD054937 is included as human proximal tubule FFPE tissue context. It should be described as tissue-context support, not as serum/plasma/urine external validation.

## Differential Expression and Harmonization

For each cohort/contrast, protein-level abundance matrices were harmonized to gene symbols and/or UniProt accessions before marker lookup. Differential expression was represented as case/comparator log2 fold change with p values and Benjamini-Hochberg style FDR values where available from the harmonized DE outputs.

The v2 tables preserve the analysis method and caveat flags per contrast in `tables/contrast_ontology_cluster_metadata_v2.tsv`. Method labels include:

- `limma_trend_robust_row_median_imputed_after_strict_filter`
- `limma_like_empirical_bayes_moderated_t_20260709`

Marker-level rows are in `tables/marker10_concordance_by_contrast_v2.tsv`.

For writing and downstream summarization, the most complete merged export is:

`results/marker10_integrated_disease_context_v2/tables/marker10_merged_de_log2fc_across_studies.tsv`

This TSV has one row per marker per cohort/contrast, with 280 rows total. It joins marker-level DE/concordance results to cohort metadata so that each row contains:

- Cohort number, disease ontology domain, disease name/context, accession, case/comparator labels, and sample counts.
- Source collection, DE method label, display contrast, and caveat flags.
- OP reference log2FC and OP reference direction.
- External cohort marker detection status, display gene, feature ID, UniProt ID, and protein name.
- External disease/cohort log2FC, p value, FDR, absolute log2FC, DE tier, and disease-cohort direction.
- Direction agreement and concordance status against the OP reference direction.
- A short `summary_sentence_fragment` field to help draft result summaries while preserving missing/discordant markers.

## Concordance Definitions

Each marker-by-cohort cell was classified against the OP reference direction:

| Status | Meaning |
| --- | --- |
| `concordant_dep` | Marker is differentially expressed and changes in the same direction as the OP reference. |
| `discordant_dep` | Marker is differentially expressed and changes in the opposite direction. |
| `concordant_exploratory` | Marker is exploratory-tier and changes in the same direction. |
| `discordant_exploratory` | Marker is exploratory-tier and changes in the opposite direction. |
| `weak_same_direction` | Marker is detected and has the same direction but does not meet a DE/exploratory tier. |
| `weak_opposite_direction` | Marker is detected and has the opposite direction but does not meet a DE/exploratory tier. |
| `not_detected_or_unmapped` | Marker was not detected or could not be mapped in that cohort. |

## Output Panels

| Panel | File | Purpose |
| --- | --- | --- |
| A | `figures/panel_A_marker10_multibar_by_cohort.png` | Multi-bar marker log2FC view across ontology-ordered cohorts. |
| B | `figures/panel_B_marker10_concordance_by_human_disease_cohort.png` | Concordance heatmap with disease name, cohort/accession, and case/control n. |
| C | `figures/panel_C_marker10_network_proximity_consensus.png` | Marker10 network/proximity consensus across disease contexts. |
| C individual | `figures/panel_C_individual_disease_maps/` | One detailed disease/proximity map per cohort, saved as PNG and PDF. |
| D | `figures/panel_D_disease_pathway_proximity_overlap.png` | Cross-disease pathway-proximity overlap. |
| E | `figures/panel_E_marker10_pathway_database_roles_enhanced.png` | STRING/database role matrix with term IDs, FDR values, and marker membership. |

## Summary Counts

| Metric | Count |
| --- | ---: |
| Human disease contrasts | 28 |
| Marker-concordance cells | 280 |
| Disease-specific pathway/proximity maps | 28 |
| STRING/database role terms used | 19 |
| Concordant DEP cells | 15 |
| Discordant DEP cells | 10 |
| Concordant exploratory cells | 3 |
| Discordant exploratory cells | 0 |
| Same-direction weak cells | 56 |
| Opposite-direction weak cells | 66 |
| Not detected or unmapped cells | 130 |

## Claim Boundary

This analysis can support the following type of statement:

> The fixed OP-associated 10-protein panel shows partial concordance and pathway proximity across independent human disease proteomics cohorts, supporting external disease-context plausibility for the marker set.

It should not be used to state that:

- Public non-OP disease cohorts prove OP-induced organ injury.
- The 10 markers are clinically diagnostic without independent OP-exposed validation.
- Pathway proximity proves a causal mechanism in the primary OP cohort.
- In vitro mechanistic assays prove in vivo DNA damage in participants.

## Reproducibility

The v2 figure/table builder is saved at:

`results/marker10_integrated_disease_context_v2/reproducibility/build_marker10_integrated_context_v2.py`

Primary result tables and figures are saved under:

`results/marker10_integrated_disease_context_v2/`
