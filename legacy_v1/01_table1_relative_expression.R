## ---------------------------------------------------------------------------
## 01_table1_relative_expression.R
## Manuscript Table 1: library-level positive control.
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
out <- list()

for (cohort in names(COHORTS)) {
  cat("\n== ", cohort, " ==\n", sep = "")
  libraries <- read_library_sheet(cohort)
  dat <- load_library_tpm(libraries, KALLISTO_ROOT,
                          abundance_file = "abundance.tsv",
                          expected_features = N_TRANSCRIPTS)

  is_orphan <- orphan_index(dat$id, ORPHAN_ID_FILE,
                            expected_n = N_ORPHAN_TRANSCRIPT)
  cat("  libraries: ", ncol(dat$tpm),
      " transcripts: ", nrow(dat$tpm),
      " orphan: ", sum(is_orphan), "\n", sep = "")

  inc <- orphan_lower_expression(dat$tpm, is_orphan,
                                 zero = "included", alpha = ALPHA_TABLE1)
  exc <- orphan_lower_expression(dat$tpm, is_orphan,
                                 zero = "excluded", alpha = ALPHA_TABLE1)

  out[[cohort]] <- data.frame(
    GEO_ID = cohort,
    n_libraries = ncol(dat$tpm),
    included_ge_alpha = inc$counts$ge,
    included_lt_alpha = inc$counts$lt,
    excluded_ge_alpha = exc$counts$ge,
    excluded_lt_alpha = exc$counts$lt,
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(
      library_id = libraries$library_id,
      p_zero_included = inc$p,
      q_zero_included = inc$q,
      p_zero_excluded = exc$p,
      q_zero_excluded = exc$q
    ),
    file.path(RESULTS_DIR, paste0("table1_pvalues_", cohort, ".csv")),
    row.names = FALSE
  )
  print(out[[cohort]])
  rm(dat); gc()
}

table1 <- do.call(rbind, out)
write.csv(table1, file.path(RESULTS_DIR, "table1.csv"), row.names = FALSE)
cat("\nTable 1 written to ", file.path(RESULTS_DIR, "table1.csv"), "\n", sep = "")
print(table1)
