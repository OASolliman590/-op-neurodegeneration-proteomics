# External Validation v2 Method Note

This document describes the methods used by the external disease-context validation workflow. It intentionally avoids storing generated results, figures, or merged DE tables.

## Objective

The workflow evaluates whether a fixed OP-associated protein marker panel shows concordance or pathway proximity in independent human disease proteomics cohorts.

The analysis is external disease-context support. It is not direct validation of OP exposure in independent OP-exposed cohorts, and it does not establish organ injury, clinical diagnostic utility, or causal mechanism.

## Inputs

Expected local/private inputs:

- Fixed marker list.
- Discovery reference direction for each marker.
- Curated public cohort metadata.
- Processed protein-level abundance matrices or harmonized DE tables.
- Sample annotation files mapping samples to case/comparator groups.

The public repository contains method code and neutral marker identifiers only. Generated outputs and discovery result spreadsheets should remain outside public Git history.

## Cohort Eligibility

Prioritize cohorts that have:

- Human samples or human tissue context.
- Patient-level disease/comparator labels.
- Processed protein-level abundance data.
- A case/control, case/comparator, disease-stage, prognosis, or subphenotype contrast.
- Protein identifiers that can be harmonized to gene symbols and/or UniProt accessions.

Exclude or defer:

- Non-human model-only studies.
- Cell-line-only or in vitro-only studies.
- Raw-only deposits without processed protein quantification.
- Peptide-only data without protein rollup.
- PTM-only, glycoproteomics-only, phosphoproteomics-only, or spectral-count-only datasets when incompatible with the main analysis.
- Datasets with opaque sample IDs and no usable label key.

## Harmonization

For each cohort, build a harmonized protein table with:

- Source feature identifier.
- Gene symbol.
- UniProt accession, when available.
- Protein name, when available.
- Protein-level abundance or DE statistics.

When a marker maps to multiple source rows, document the row-selection rule in the analysis script or output metadata. Typical choices are strongest DE tier, lowest FDR, lowest p value, then largest absolute log2FC.

## Differential Expression

External disease cohorts are modeled as case versus comparator contrasts. Each contrast should preserve:

- Case label and comparator label.
- Case and comparator sample sizes.
- DE method label.
- Filtering rule.
- Missingness/imputation rule.
- log2FC direction convention.
- p value and FDR calculation.
- Cohort caveat flags.

The prior v2 implementation used limma-style moderated testing, robust/trend settings where available, strict filtering, and row-median imputation after filtering. Equivalent methods may be used if documented per contrast.

## Marker Concordance

For each marker and cohort contrast:

1. Determine whether the marker is detected and mappable.
2. Extract external cohort log2FC, p value, FDR, and tier.
3. Compare the external log2FC direction with the local OP reference direction.
4. Assign a concordance status.

Recommended status labels:

| Status | Meaning |
| --- | --- |
| `concordant_dep` | Differentially expressed and same direction as OP reference |
| `discordant_dep` | Differentially expressed and opposite direction |
| `concordant_exploratory` | Exploratory-tier signal in same direction |
| `discordant_exploratory` | Exploratory-tier signal in opposite direction |
| `weak_same_direction` | Detected, same direction, but weak/non-tiered |
| `weak_opposite_direction` | Detected, opposite direction, but weak/non-tiered |
| `not_detected_or_unmapped` | Marker absent or not mappable |

Missing and discordant markers should be retained in tables and plots to avoid over-selection.

## Disease Ontology Ordering

Disease cohorts can be ordered by broad ontology/domain for interpretation. Suggested domains include:

- Nervous system.
- Metabolic/endocrine.
- Kidney/urinary.
- Hepatic.
- Cardiovascular/inflammatory.

Within each domain, cohorts may be ordered by disease context, accession, or marker-signature similarity.

## Pathway and Network Context

STRING or related databases can be used to evaluate:

- Marker-marker interaction context.
- Enriched compartments/pathways.
- Disease-protein proximity to the fixed marker panel.
- Cross-disease pathway overlap.

Report database name, access date, species/taxid, confidence cutoff, enrichment background, and multiple-testing correction where applicable.

These analyses provide biological context only. They should not be framed as causal proof.

## Recommended Non-Versioned Outputs

When running the analysis locally, a useful output layout is:

```text
local_outputs/
  tables/
  figures/
  reproducibility/
```

These files should stay local or be shared through a deliberate release package, not committed by default.

## Claim Boundary

Acceptable:

> The fixed OP marker panel was evaluated for external disease-context concordance and pathway proximity across eligible public human proteomics cohorts.

Avoid:

- Directly claiming independent OP validation from non-OP disease cohorts.
- Claiming organ injury based only on disease-context similarity.
- Claiming clinical diagnostic utility from marker concordance alone.
- Treating pathway proximity or STRING enrichment as causal mechanism.
