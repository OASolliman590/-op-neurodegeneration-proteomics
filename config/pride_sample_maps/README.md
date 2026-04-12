# Manual PRIDE Sample Maps

Optional override maps for cohorts that have quant tables but no machine-readable case/control labels.

Create one CSV per accession named:
- `PXDxxxxxx.csv`

Required columns:
- `sample`: token to match inside intensity column names (case-insensitive)
- `label`: `case` or `control`

Example:

```csv
sample,label
c2,control
c3,control
c4,control
c8,control
cp3,control
f5,case
f9,case
f13,case
f15,case
f16,case
s6,case
s8,case
s10,case
s11,case
s17,case
```

The pipeline loads these files via `--sample-map-dir` (default: `config/pride_sample_maps`).

Notes:
- `PXD041336.csv` is intentionally empty. The selected workbook in this accession is an interaction-level table and does not expose sample-level case/control intensity columns.
