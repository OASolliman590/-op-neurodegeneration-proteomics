# Marker10 Integrated Disease Context v2

This folder contains the July 2026 external disease-context validation outputs for the fixed OP-associated 10-protein panel.

## Contents

- `figures/`: publication-facing composite and standalone Panels A-E as PNG/PDF.
- `figures/panel_C_individual_disease_maps/`: 28 cohort-specific disease/proximity maps as PNG/PDF.
- `tables/`: cohort metadata, concordance tables, pathway/proximity scores, and STRING/database role definitions.
- `reproducibility/`: script used to generate the v2 figure/table bundle.
- `marker10_integrated_disease_context_v2_report.md`: compact run report.

## Primary Summary TSV

Use this file for result writing:

`tables/marker10_merged_de_log2fc_across_studies.tsv`

It contains 280 rows, one for each Marker10 protein in each of the 28 human disease contrasts. Columns include disease ontology, accession, case/comparator labels, n, DE method, caveat flags, OP reference log2FC/direction, external cohort log2FC, p value, FDR, DE tier, marker mapping fields, concordance status, and a short summary sentence fragment.

## Scope

| Metric | Count |
| --- | ---: |
| Human disease contrasts | 28 |
| Marker-concordance cells | 280 |
| Disease-specific pathway/proximity maps | 28 |
| STRING/database role terms used | 19 |

## Interpretation

These outputs document disease-context concordance and pathway proximity for the fixed marker panel. They should be cited as external biological context, not as direct proof of OP-induced organ injury or clinical diagnostic performance.
