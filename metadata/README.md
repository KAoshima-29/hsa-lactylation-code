# metadata/

Place sample metadata required by the analysis scripts here.

**Required columns** (CSV):
- `sample`       — must match column names in the count matrices
- `group`        — biological grouping (e.g., HSA2, HSA3, EC_A, EC_B, etc.)
- `study`        — e.g., HU or Public
- `instrument`   — sequencing instrument model
- `lib_prep`     — library preparation protocol

Typical file:
- `Completed_pheno_table.csv`

This folder is git-ignored in `.gitignore`.
