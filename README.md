# OP Proteomics External Validation

External disease-context validation of a fixed 10-protein chronic organophosphate exposure signature across independent public human proteomics cohorts.

This repository supports the external-validation component of the OP proteomics manuscript. The central question is whether a pre-specified 10-marker OP-associated protein panel shows concordant differential-expression behavior or pathway/network proximity in independent human disease proteomics datasets.

This is not a claim that the public disease cohorts prove OP-induced organ injury. The disease cohorts are used as external biological context for marker plausibility, directionality, and pathway overlap.

## Marker Panel

The 10 proteins are fixed before external comparison:

ACTG1, DNAH9, GPX3, VWF, C4B, CD44, CFHR2, ITIH3, LRG1, MYH7B

## Current External-Validation Scope

The v2 analysis expands beyond the older neurodegeneration-only repository framing. It includes 28 human disease contrasts ordered by disease ontology:

| Domain | Included disease contexts |
| --- | --- |
| Nervous system | Alzheimer disease, Parkinson disease, multiple sclerosis |
| Metabolic/endocrine | Diabetes mellitus and metabolic dysfunction |
| Kidney/urinary | Chronic kidney disease, diabetic kidney disease/prognosis, AKI subphenotype, kidney dysfunction, proximal tubule tissue context |
| Hepatic | Liver dysfunction |
| Cardiovascular/inflammatory | Cardiac disease and MIS-C/cardiac involvement |

## Main Result Bundle

The current publication-facing outputs are in:

`results/marker10_integrated_disease_context_v2/`

Important files:

| Path | Description |
| --- | --- |
| `figures/marker10_integrated_disease_context_v2_composite.png` | Main composite figure with Panels A-E |
| `figures/panel_A_marker10_multibar_by_cohort.png` | Marker10 log2FC multi-bar chart across human disease cohorts |
| `figures/panel_B_marker10_concordance_by_human_disease_cohort.png` | Marker concordance heatmap with disease/accession/case-control n |
| `figures/panel_C_marker10_network_proximity_consensus.png` | Marker10 STRING/pathway proximity consensus |
| `figures/panel_C_individual_disease_maps/` | 28 cohort-specific disease/proximity maps, PNG and PDF |
| `figures/panel_D_disease_pathway_proximity_overlap.png` | Cross-disease pathway-proximity overlap heatmap |
| `figures/panel_E_marker10_pathway_database_roles_enhanced.png` | Expanded database role/pathway membership matrix |
| `tables/contrast_ontology_cluster_metadata_v2.tsv` | Cohort metadata, disease ontology, sample counts, DE method, caveats |
| `tables/marker10_merged_de_log2fc_across_studies.tsv` | ChatGPT/summary-ready merged DE table for all 10 markers across all 28 contrasts |
| `tables/marker10_concordance_by_contrast_v2.tsv` | Full 280-cell marker-by-cohort concordance table |
| `tables/panel_E_pathway_database_role_definitions.tsv` | STRING/database role definitions and FDR values |
| `reproducibility/build_marker10_integrated_context_v2.py` | Figure/table builder used for v2 outputs |

## Technical Method Summary

The v2 analysis uses a fixed-marker external-validation design:

1. The 10 OP-associated proteins were locked before external disease-cohort testing.
2. Eligible public human proteomics cohorts were curated for protein-level abundance, usable case/comparator labels, and interpretable human disease context.
3. Protein identifiers were harmonized to gene symbols and/or UniProt accessions before marker lookup.
4. Each cohort was analyzed as a disease/comparator contrast and represented as case-vs-control or case-vs-comparator log2FC.
5. Per-contrast DE outputs preserve p values, FDR values, DE tier, source collection, method label, and caveat flags.
6. Each marker/cohort row was classified against the OP reference direction as concordant DEP, discordant DEP, exploratory, weak same/opposite direction, or not detected/unmapped.
7. Disease contexts were ordered by ontology/domain for cross-disease interpretation.
8. STRING/database role and pathway proximity were used as contextual biology, not as causal proof.

The most useful file for result writing is:

`results/marker10_integrated_disease_context_v2/tables/marker10_merged_de_log2fc_across_studies.tsv`

This table has one row per marker per cohort/contrast: 10 markers x 28 contrasts = 280 rows. It includes disease metadata, accession, sample counts, OP reference direction/log2FC, external log2FC, p value, FDR, DE tier, mapping fields, caveats, and a short summary sentence fragment for drafting.

## v2 Summary Counts

| Metric | Count |
| --- | ---: |
| Human disease contrasts | 28 |
| Marker-concordance cells | 280 |
| Disease-specific pathway/proximity maps | 28 |
| STRING/database role terms used | 19 |
| Concordant DEP cells | 15 |
| Discordant DEP cells | 10 |
| Concordant exploratory cells | 3 |
| Same-direction weak cells | 56 |
| Opposite-direction weak cells | 66 |
| Not detected or unmapped cells | 130 |

## Claim-Safe Interpretation

Recommended wording:

> The fixed OP-associated 10-protein panel shows partial concordance and pathway proximity across independent human disease proteomics cohorts, supporting external disease-context plausibility for the marker set.

Avoid wording that says the external disease cohorts prove OP-induced organ injury, clinical diagnostic utility, or in vivo pathology in the original OP participants.

## Documentation

- `docs/external_validation_v2.md`: technical methods, cohort scope, concordance definitions, and claim boundaries.
- `VALIDATION_UPDATES.md`: historical session notes from the earlier neurodegeneration-focused analysis.

## Original Repository Components

The older pipeline scripts and configuration files are retained for provenance:

- `src/query_eprot.py`
- `src/query_pride.py`
- `src/grand_finale_figure.py`
- `src/cross_modal_figure.py`
- `src/concordance_by_tissue.py`
- `src/analyze_meta.py`
- `config/targets.yaml`
- `config/cohorts.csv`

## Data Sources

| Resource | URL |
| --- | --- |
| PRIDE Archive | https://www.ebi.ac.uk/pride |
| ProteomeXchange | https://www.proteomexchange.org |
| Expression Atlas | https://www.ebi.ac.uk/gxa/experiments |
| Human Protein Atlas | https://www.proteinatlas.org |
| ProteomicsDB | https://www.proteomicsdb.org |
| STRING | https://string-db.org |

## License

Code: MIT. Public proteomics data remain under the terms of their source repositories and publications.
