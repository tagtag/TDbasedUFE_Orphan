## ---------------------------------------------------------------------------
## 04_clinical_association.R
## Manuscript clinical-variable analyses (Methods 3.5 / Table 7).
##
## IMPORTANT: this script does NOT use separate clinical_GSE*.csv files.
## All clinical variables come from the same table used in the original local
## analysis:
##
##   metadata_output/GEO_patient_metadata_combined.tsv
##
## Analysed patients are aligned to that table through normal_run_accession,
## reproducing the original pattern
##
##   index <- match(normal_run_accession, y$normal_run_accession)
##
## rather than relying on row order or patient labels.
## ---------------------------------------------------------------------------

source("R/config.R")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

normalize_run <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    s <- sub("_kallisto$", "", basename(s), ignore.case = TRUE)
    m <- regexpr("(?:SRR|ERR|DRR)[0-9]+", s, perl = TRUE, ignore.case = TRUE)
    if (m[1] >= 0) return(toupper(regmatches(s, m)))
    toupper(s)
  }, character(1))
}


if (!file.exists(GEO_PATIENT_METADATA_FILE)) {
  stop(
    "Combined GEO metadata file not found: ", GEO_PATIENT_METADATA_FILE, "\n",
    "Create it first with:\n",
    "  Rscript scripts/00_build_geo_patient_metadata.R\n",
    "or place the exact GEO_patient_metadata_combined.tsv used in the original ",
    "analysis at that path."
  )
}

y <- read.delim(
  GEO_PATIENT_METADATA_FILE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "")
)
if (!"cohort" %in% colnames(y) && "GSE" %in% colnames(y)) y$cohort <- y$GSE

required_meta <- c(
  "cohort", "normal_run_accession", "PASI", "age", "stage", "alcohol", "smoking"
)
missing_meta <- setdiff(required_meta, colnames(y))
if (length(missing_meta)) {
  stop(
    "GEO_patient_metadata_combined.tsv is missing columns: ",
    paste(missing_meta, collapse = ", ")
  )
}
y$normal_run_accession <- normalize_run(y$normal_run_accession)
if ("disease_run_accession" %in% colnames(y)) {
  y$disease_run_accession <- normalize_run(y$disease_run_accession)
}

## ---------------------------------------------------------------------------
## Patient score and accession matching
## ---------------------------------------------------------------------------
load_score <- function(cohort,
                       table = CLINICAL_SCORE_TABLE,
                       zero = CLINICAL_SCORE_ZERO,
                       p_column = CLINICAL_P_COLUMN) {
  score_file <- file.path(
    RESULTS_DIR,
    sprintf("perpatient_%s_%s_zero_%s.csv", table, cohort, zero)
  )
  if (!file.exists(score_file)) {
    stop("Patient-level result not found. Run R/02_patient_level_analysis.R first: ",
         score_file)
  }

  res <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("patient_id", p_column) %in% colnames(res))) {
    stop(score_file, " must contain patient_id and ", p_column)
  }

  ## The sample sheet identifies the normal RNA-seq run belonging to each
  ## patient.  Matching clinical metadata by this accession reproduces the
  ## original use of y$normal_run_accession.
  sh <- read_sheet(cohort)
  normal <- sh[sh$condition == "normal",
               c("patient_id", "library_id", "directory"), drop = FALSE]
  idx <- match(res$patient_id, normal$patient_id)
  if (anyNA(idx)) {
    stop("Could not find a normal library for all analysed patients in ", cohort)
  }

  run <- normalize_run(normal$library_id[idx])
  bad <- is.na(run) | !grepl("^(SRR|ERR|DRR)[0-9]+$", run)
  if (any(bad)) {
    run2 <- normalize_run(normal$directory[idx][bad])
    run[bad] <- run2
  }
  if (any(is.na(run) | !grepl("^(SRR|ERR|DRR)[0-9]+$", run))) {
    stop(
      "Could not derive normal SRA/ENA run accession(s) for ", cohort,
      ". Ensure library_id or directory in the sample sheet contains SRR/ERR/DRR accessions."
    )
  }

  p <- suppressWarnings(as.numeric(res[[p_column]]))
  data.frame(
    patient_id = res$patient_id,
    normal_run_accession = run,
    p_value = p,
    score = rank(p, ties.method = "average", na.last = "keep"),
    stringsAsFactors = FALSE
  )
}

attach_geo_metadata <- function(cohort) {
  sc <- load_score(cohort)
  md <- y[y$cohort == cohort, , drop = FALSE]
  if (!nrow(md)) {
    stop("No rows for ", cohort, " in ", GEO_PATIENT_METADATA_FILE)
  }

  dup <- duplicated(md$normal_run_accession) & !is.na(md$normal_run_accession)
  if (any(dup)) {
    stop(
      "Duplicate normal_run_accession values in combined metadata for ", cohort,
      ": ", paste(unique(md$normal_run_accession[dup]), collapse = ", ")
    )
  }

  index <- match(sc$normal_run_accession, md$normal_run_accession)
  if (anyNA(index)) {
    bad <- sc$normal_run_accession[is.na(index)]
    stop(
      "Analysed normal runs were not found in ", GEO_PATIENT_METADATA_FILE,
      " for ", cohort, ":\n  ", paste(bad, collapse = "\n  ")
    )
  }

  md2 <- md[index, , drop = FALSE]
  ## Avoid duplicating columns that already exist in sc.
  add <- setdiff(colnames(md2), c("patient_id", "normal_run_accession"))
  out <- cbind(sc, md2[, add, drop = FALSE])

  audit_file <- file.path(RESULTS_DIR, paste0("clinical_metadata_matched_", cohort, ".tsv"))
  write.table(out, audit_file, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  out
}

require_clinical <- function(df, cohort, fields, min_n = 3L) {
  absent <- fields[!fields %in% colnames(df)]
  if (length(absent)) {
    stop(cohort, ": missing clinical column(s): ", paste(absent, collapse = ", "))
  }

  all_missing <- fields[vapply(fields, function(v) all(is.na(df[[v]])), logical(1))]
  if (length(all_missing)) {
    stop(
      cohort, ": clinical field(s) are entirely missing: ",
      paste(all_missing, collapse = ", "), "\n",
      "Do not fabricate them. Rebuild the public-data metadata with:\n",
      "  Rscript scripts/00_build_geo_patient_metadata.R\n",
      "For GSE127165, patient age is retrieved automatically from the published ",
      "Additional file 1, Table S1; inspect metadata_output/GSE127165_supplement_clinical.tsv ",
      "if reconstruction fails."
    )
  }

  keep <- complete.cases(df[, c("score", fields), drop = FALSE])
  z <- df[keep, , drop = FALSE]
  if (nrow(z) < min_n) {
    stop(cohort, ": only ", nrow(z), " complete cases for score + ",
         paste(fields, collapse = ", "))
  }
  z
}

report_lm <- function(fit, label, file) {
  cat("\n--- ", label, " (n = ", stats::nobs(fit), ") ---\n", sep = "")
  s <- summary(fit)
  print(s)
  capture.output(s, file = file)
  invisible(s)
}

## ---------------------------------------------------------------------------
## GSE244679: Bayesian correlation between rank(P) and PASI
## ---------------------------------------------------------------------------
cohort <- "GSE244679"
cat("\n########  ", cohort, " : PASI  ########\n", sep = "")
sc <- require_clinical(attach_geo_metadata(cohort), cohort, "PASI")
cat("Complete cases: ", nrow(sc), "\n", sep = "")

if (!requireNamespace("BayesFactor", quietly = TRUE)) {
  stop(
    "BayesFactor is required for the manuscript PASI analysis. Install it with:\n",
    "  install.packages('BayesFactor')"
  )
}

bf <- BayesFactor::correlationBF(sc$score, sc$PASI, rscale = 1/3)
set.seed(BASE_SEED)
post <- BayesFactor::posterior(bf, iterations = 10000, progress = FALSE)
if (!"rho" %in% colnames(post)) {
  stop("BayesFactor posterior did not contain a rho column")
}
rho <- as.numeric(post[, "rho"])
qs <- stats::quantile(rho, c(0.025, 0.5, 0.975), na.rm = TRUE)

pasi_out <- data.frame(
  lower_2.5 = unname(qs[1]),
  median_50 = unname(qs[2]),
  upper_97.5 = unname(qs[3]),
  posterior_Pr_rho_gt_0 = mean(rho > 0, na.rm = TRUE),
  n = nrow(sc)
)
write.csv(pasi_out, file.path(RESULTS_DIR, "table7_pasi_correlation.csv"),
          row.names = FALSE)
write.table(sc, file.path(RESULTS_DIR, "clinical_GSE244679_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(pasi_out)

## ---------------------------------------------------------------------------
## GSE127165: rank(P) ~ age + alcohol + stage
## ---------------------------------------------------------------------------
cohort <- "GSE127165"
cat("\n########  ", cohort, " : age + alcohol + stage  ########\n", sep = "")
sc <- require_clinical(
  attach_geo_metadata(cohort), cohort, c("age", "alcohol", "stage")
)
sc$alcohol <- droplevels(factor(sc$alcohol))
sc$stage <- droplevels(factor(sc$stage))
fit <- lm(score ~ age + alcohol + stage, data = sc)
report_lm(
  fit,
  "rank(P) ~ age + alcohol + stage",
  file.path(RESULTS_DIR, "clinical_lm_GSE127165.txt")
)
write.table(sc, file.path(RESULTS_DIR, "clinical_GSE127165_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ---------------------------------------------------------------------------
## GSE40419: rank(P) ~ age + smoking + stage
## ---------------------------------------------------------------------------
cohort <- "GSE40419"
cat("\n########  ", cohort, " : age + smoking + stage  ########\n", sep = "")
sc <- require_clinical(
  attach_geo_metadata(cohort), cohort, c("age", "smoking", "stage")
)
sc$smoking <- droplevels(factor(sc$smoking))
sc$stage <- droplevels(factor(sc$stage))
fit <- lm(score ~ age + smoking + stage, data = sc)
report_lm(
  fit,
  "rank(P) ~ age + smoking + stage",
  file.path(RESULTS_DIR, "clinical_lm_GSE40419.txt")
)
write.table(sc, file.path(RESULTS_DIR, "clinical_GSE40419_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\nGSE144269: no clinical-label association analysis was performed.\n")
