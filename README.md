# OP Proteomics External Validation

Methods-only repository for the external disease-context validation workflow used to evaluate a fixed 10-protein chronic organophosphate exposure marker panel against independent public human proteomics cohorts.

Generated result tables, figures, merged DE exports, and manuscript-ready output bundles are intentionally not versioned in this public repository.

## Scope

This repository documents and scripts the external-validation workflow only. It is designed to:

1. Curate eligible public human proteomics cohorts.
2. Verify that cohorts contain processed protein-level abundance data and interpretable sample labels.
3. Harmonize protein identifiers to gene symbols and UniProt accessions where possible.
4. Run or ingest cohort-level differential-expression contrasts.
5. Compare external disease-cohort directionality against a fixed OP marker panel.
6. Classify marker-by-cohort evidence as concordant, discordant, exploratory, weak, or not detected.
7. Use pathway/network databases as contextual biological support.

The workflow does not claim that non-OP disease cohorts prove OP-induced organ injury, clinical diagnostic utility, or causal mechanism in the original OP-exposed participants.

## Marker Panel

The external-validation workflow uses a fixed 10-protein panel:

ACTG1, DNAH9, GPX3, VWF, C4B, CD44, CFHR2, ITIH3, LRG1, MYH7B

The public repository keeps marker identities and analysis code. Discovery fold changes, generated DE tables, and figures should be regenerated or supplied locally and should not be committed.

## Repository Contents

| Path | Purpose |
| --- | --- |
| `config/targets.yaml` | Neutral marker manifest with gene/protein identifiers only |
| `config/cohorts.csv` | Cohort/source manifest used by the legacy external-validation workflow |
| `config/external_data_manifest_v2.csv` | Candidate public-data source manifest |
| `config/pride_overrides.csv` | Manual PRIDE metadata corrections |
| `config/sciphera_queries.yaml` | Search/query configuration |
| `src/` | Query, harmonization, DE, visualization, and QA scripts |
| `docs/external_validation_v2.md` | Technical method note for the v2 workflow |
| `requirements.txt` | Python dependencies |

## What Is Not Versioned

The following are intentionally excluded from the public repository:

- Generated result tables.
- Merged marker-by-cohort DE/log2FC TSV exports.
- Volcano plots, heatmaps, pathway/network figures, and composite figures.
- Local manuscript figure bundles.
- Discovery fold-change spreadsheets or other primary result files.

Keep these artifacts in local analysis folders, private storage, or release archives only when appropriate.

## Technical Method Summary

### Cohort Curation

Candidate cohorts are screened for:

- Human samples or human tissue context.
- Patient disease/comparator labels.
- Processed protein-level abundance matrices.
- A case/control, case/comparator, disease-stage, prognosis, or subphenotype contrast that can be interpreted.
- Mappable protein identifiers, preferably gene symbols and/or UniProt accessions.

Datasets are excluded or deferred when they are raw-only, non-human model-only, cell-line-only, peptide-only without protein rollup, PTM-only, or lack a usable sample-label key.

### Identifier Harmonization

Protein rows are normalized to stable marker lookup fields:

- Gene symbol.
- UniProt accession.
- Source feature identifier.
- Protein name, when available.

When multiple source identifiers map to the same marker, the selected row should be documented by the harmonization script or downstream DE table.

### Differential Expression

External cohorts are analyzed as disease/comparator contrasts. The expected output schema for each contrast is:

- `marker` or harmonized gene symbol.
- External cohort `log2FC` for case versus comparator.
- `p_value`.
- `FDR`.
- Optional DE tier.
- Method label and caveat flags.

The previous v2 implementation used limma-style moderated testing, robust trend options where available, strict filtering, and row-median imputation after filtering. Cohort-specific caveats should remain attached to each contrast.

### Concordance Classification

Each marker-by-cohort result is compared against the fixed OP reference direction from the local/private discovery source. The recommended categories are:

| Status | Meaning |
| --- | --- |
| `concordant_dep` | Differentially expressed and same direction as the OP reference |
| `discordant_dep` | Differentially expressed and opposite direction |
| `concordant_exploratory` | Exploratory-tier signal in the same direction |
| `discordant_exploratory` | Exploratory-tier signal in the opposite direction |
| `weak_same_direction` | Detected, same direction, but not tiered as DE/exploratory |
| `weak_opposite_direction` | Detected, opposite direction, but not tiered as DE/exploratory |
| `not_detected_or_unmapped` | Marker absent from the cohort table or not mappable |

Missing and discordant markers should remain visible in all summary tables.

### Pathway and Network Context

STRING or related pathway/network resources can be used to contextualize the marker panel and disease-cohort proteins. These analyses should be described as pathway proximity or biological context, not as causal proof.

## Suggested Local Output Layout

Generated artifacts can be written locally using a non-versioned structure such as:

```text
local_outputs/
  tables/
  figures/
  reproducibility/
```

Keep these folders out of public Git commits unless there is a deliberate release decision.

## Claim Boundary

Appropriate framing:

> The fixed OP-associated marker panel can be compared with independent human disease proteomics cohorts to assess disease-context concordance and pathway proximity.

Avoid:

- Claiming that non-OP disease cohorts prove OP-induced organ injury.
- Claiming clinical diagnostic utility without independent OP-exposed validation.
- Treating pathway proximity as causal mechanism.
- Treating in vitro mechanistic assays as proof of in vivo participant pathology.

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
