# legacy_v1 — scripts of the first version

These scripts reproduce the tables and figure of **Preprints 202608.1930.v1**.
They are kept because that version is public and is cited by the revision; they are
**not** part of the revised analysis.

| script | produced |
|---|---|
| `01_table1_relative_expression.R` | old Table 1 — library-level lower-relative-expression counts |
| `02_patient_level_analysis.R` | old Tables 2, 3, 4 — standardized expression, gene loci, ranks |
| `03_table5_control_sets.R` | old Table 5 — expression-matched controls, 20 replicates |
| `04_clinical_association.R` | old Table 7 — PASI correlation and multivariate regression |
| `05_figure2.R` | old Figure 2 |
| `07_validate_manuscript_counts.R` | consistency checks on the counts quoted in v1 |
| `run_analysis.R`, `run_all.sh` | the v1 runners |

## Why they were superseded

The revision withdraws the analysis these scripts implement, for reasons the
reviewers raised and we confirmed:

- The per-library and per-patient nominal p values, and the Benjamini–Hochberg
  adjustment built on them, were replaced by a resampling reference that respects
  the correlation between transcripts.
- The claim of patient-level directional heterogeneity does not survive a negative
  control that breaks the within-patient pairing while keeping the condition
  contrast. The same proportions of significant and positive coefficients are
  obtained from arbitrary pairings.
- The control transcripts in old Table 5 were matched on the mean of the patient's
  own paired normal and disease values (`pair_mean`). On synthetic data in which
  low-abundance transcripts move together, that matching produces false positives in
  every dataset tested. The revision matches on other patients' normal mean and on
  cross-patient detection frequency instead.
- GSE40419 (LAC) is excluded from the revised primary analysis: its normal and
  disease libraries belong to disjoint ENA submissions and differ systematically in
  depth, pseudoalignment rate and zero rate, so condition cannot be separated from
  batch.

These scripts expect the repository root as the working directory, as before. They
share `R/config.R` and `R/00_functions.R` with the current analysis; if those files
change in a way that breaks them, that is expected, and the v1 state is recoverable
from the git history.
