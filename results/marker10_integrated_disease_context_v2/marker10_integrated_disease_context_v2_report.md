# Marker10 Integrated Disease Context v2 Report

## What Changed From v1
- Panel A is now a cohort/disease Marker10 multi-bar chart, not a CKD reference lollipop.
- Panel B now has a standalone readable concordance heatmap with disease, accession/cohort number, case/control labels, and n.
- Panel C is now the Marker10 STRING/pathway proximity network and also saves one disease-specific map per cohort.
- Panel D is now the cross-disease pathway-proximity overlap heatmap.
- Panel E is expanded with role definitions, database category, term ID, FDR, and marker membership.

## Scope
- Human disease contrasts: 28
- Marker-concordance cells: 280
- Disease-specific pathway-proximity maps: 28
- STRING/database role terms used: 19

## Main Outputs
- Composite v2: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/marker10_integrated_disease_context_v2_composite.png`
- Panel A standalone: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_A_marker10_multibar_by_cohort.png`
- Panel B standalone: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_B_marker10_concordance_by_human_disease_cohort.png`
- Panel C network: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_C_marker10_network_proximity_consensus.png`
- Panel C individual maps: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_C_disease_marker_proximity_maps`
- Panel D overlap: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_D_disease_pathway_proximity_overlap.png`
- Panel E enhanced: `op_external_validation/results/de_visualization_hub_network/marker10_integrated_disease_context_v2_20260709/figures/panel_E_marker10_pathway_database_roles_enhanced.png`

## Status Counts
- Concordant DEP: 15
- Discordant DEP: 10
- Concordant exploratory: 3
- Discordant exploratory: 0
- Same direction, weak: 56
- Opposite direction, weak: 66
- No direction: 0
- Not detected: 130
