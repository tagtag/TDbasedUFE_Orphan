## ---------------------------------------------------------------------------
## 02_patient_level_analysis.R
## Manuscript Tables 2-4 and per-patient statistics.
##
##   Table 2: transcript level, standardized expression
##   Table 3: gene-locus level, standardized expression
##   Table 4: transcript level, within-sample ranks
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

SPECS <- list(
  list(table = "table2", unit = "transcript", transform = "scale"),
  list(table = "table3", unit = "locus",      transform = "scale"),
  list(table = "table4", unit = "transcript", transform = "rank")
)
ZERO_RULES <- c(included = "none", excluded = "or")
all_counts <- list()

for (spec in SPECS) {
  cat("\n########  ", spec$table, ": ", spec$unit, " / ", spec$transform,
      "  ########\n", sep = "")

  abundance_file <- if (spec$unit == "transcript") "abundance.tsv" else "abundance_gene.tsv"
  id_file <- if (spec$unit == "transcript") ORPHAN_ID_FILE else ORPHAN_LOCUS_ID_FILE
  expected_orphan <- if (spec$unit == "transcript") N_ORPHAN_TRANSCRIPT else N_ORPHAN_LOCUS
  expected_features <- if (spec$unit == "transcript") N_TRANSCRIPTS else N_GENE_LOCI

  for (cohort in names(COHORTS)) {
    cat("\n== ", cohort, " ==\n", sep = "")
    sheet <- read_sheet(cohort)
    dat <- load_cohort_tpm(sheet, KALLISTO_ROOT,
                           abundance_file = abundance_file,
                           expected_features = expected_features)
    is_orphan <- orphan_index(dat$id, id_file, expected_n = expected_orphan)

    for (zr in names(ZERO_RULES)) {
      res <- run_cohort(
        dat$tpm, sheet, is_orphan,
        transform = spec$transform,
        zero_rule = ZERO_RULES[[zr]],
        alpha = ALPHA
      )
      cnt <- count_directions(res)
      cat("  zero ", zr, ": up=", cnt$up, " down=", cnt$down,
          " none=", cnt$none, "\n", sep = "")

      all_counts[[length(all_counts) + 1L]] <- data.frame(
        table = spec$table,
        unit = spec$unit,
        transform = spec$transform,
        cohort = cohort,
        zero = zr,
        cnt,
        stringsAsFactors = FALSE
      )

      write.csv(
        res,
        file.path(RESULTS_DIR,
                  sprintf("perpatient_%s_%s_zero_%s.csv", spec$table, cohort, zr)),
        row.names = FALSE
      )
    }
    rm(dat); gc()
  }
}

counts <- do.call(rbind, all_counts)
write.csv(counts, file.path(RESULTS_DIR, "tables_2_3_4_counts_long.csv"), row.names = FALSE)

for (tb in c("table2", "table3", "table4")) {
  pretty <- format_direction_table(counts, tb)
  write.csv(pretty, file.path(RESULTS_DIR, paste0(tb, ".csv")), row.names = FALSE)
}

cat("\n---------------- summary ----------------\n")
print(counts)
