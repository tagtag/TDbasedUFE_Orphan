## ---------------------------------------------------------------------------
## 03_table5_control_sets.R
## Manuscript Table 5: expression-matched non-orphan controls.
##
## This is a prespecified sensitivity analysis for the explanation that the
## rank-based orphan pattern arises solely because orphan transcripts are low
## abundance. It is not presented as a proof that the pattern is unique to
## orphan genes under every possible control-matching definition.
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
rows <- list()

for (cohort in names(COHORTS)) {
  cat("\n== ", cohort, " ==\n", sep = "")
  sheet <- read_sheet(cohort)
  dat <- load_cohort_tpm(sheet, KALLISTO_ROOT,
                         abundance_file = "abundance.tsv",
                         expected_features = N_TRANSCRIPTS)
  is_orphan <- orphan_index(dat$id, ORPHAN_ID_FILE,
                            expected_n = N_ORPHAN_TRANSCRIPT)

  ## The manuscript control analysis uses ranks and the paired zero rule:
  ## retain a transcript when it is nonzero in either member of the pair.
  orph <- run_cohort(dat$tpm, sheet, is_orphan,
                     transform = "rank", zero_rule = "or", alpha = ALPHA)
  orph_cnt <- count_directions(orph)

  ctl <- run_control_sets(
    dat$tpm, sheet, is_orphan,
    n_rep = N_REP,
    n_strata = N_STRATA,
    transform = "rank",
    zero_rule = "or",
    alpha = ALPHA,
    base_seed = BASE_SEED,
    keep_patient_details = TRUE
  )

  if (!all(ctl$n_control_transcripts == N_ORPHAN_TRANSCRIPT)) {
    stop("control sets do not all contain ", N_ORPHAN_TRANSCRIPT,
         " transcripts; range was ",
         paste(range(ctl$n_control_transcripts), collapse = "-"))
  }

  ctl_sum <- summarise_control(ctl$counts)
  rows[[length(rows) + 1L]] <- data.frame(
    GEO_ID = cohort,
    transcript_set = "Orphan",
    disease_gt_normal = as.character(orph_cnt$up),
    disease_lt_normal = as.character(orph_cnt$down),
    no_direction = as.character(orph_cnt$none),
    stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1L]] <- data.frame(
    GEO_ID = cohort,
    transcript_set = "Control",
    disease_gt_normal = ctl_sum$up,
    disease_lt_normal = ctl_sum$down,
    no_direction = ctl_sum$none,
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(replicate = seq_len(N_REP), ctl$counts, check.names = FALSE),
    file.path(RESULTS_DIR, paste0("table5_control_replicates_", cohort, ".csv")),
    row.names = FALSE
  )
  write.csv(
    ctl$details,
    file.path(RESULTS_DIR, paste0("table5_control_perpatient_", cohort, ".csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(replicate = seq_len(N_REP), ctl$n_control_transcripts,
               check.names = FALSE),
    file.path(RESULTS_DIR, paste0("table5_control_sizes_", cohort, ".csv")),
    row.names = FALSE
  )

  cat("  orphan: up=", orph_cnt$up, " down=", orph_cnt$down,
      " none=", orph_cnt$none, "\n", sep = "")
  cat("  control: up=", ctl_sum$up, " down=", ctl_sum$down,
      " none=", ctl_sum$none, "\n", sep = "")
  rm(dat); gc()
}

table5 <- do.call(rbind, rows)
write.csv(table5, file.path(RESULTS_DIR, "table5.csv"), row.names = FALSE)
cat("\nTable 5 written to ", file.path(RESULTS_DIR, "table5.csv"), "\n", sep = "")
print(table5)
