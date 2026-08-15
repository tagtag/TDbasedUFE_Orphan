## ---------------------------------------------------------------------------
## config.R -- paths, cohort definitions, and manuscript analysis settings
## ---------------------------------------------------------------------------

KALLISTO_ROOT <- Sys.getenv("KALLISTO_ROOT", unset = "data/kallisto")
RESULTS_DIR   <- Sys.getenv("RESULTS_DIR",   unset = "results")
FIGURE_DIR    <- Sys.getenv("FIGURE_DIR",    unset = "figures")
METADATA_DIR         <- Sys.getenv("METADATA_DIR",         unset = "metadata")
METADATA_OUTPUT_DIR  <- Sys.getenv("METADATA_OUTPUT_DIR",  unset = "metadata_output")
GEO_PATIENT_METADATA_FILE <- Sys.getenv(
  "GEO_PATIENT_METADATA_FILE",
  unset = file.path(METADATA_OUTPUT_DIR, "GEO_patient_metadata_combined.tsv")
)

ORPHAN_ID_FILE       <- file.path(METADATA_DIR, "orphan_transcript_ids.txt")
ORPHAN_LOCUS_ID_FILE <- file.path(METADATA_DIR, "orphan_gene_ids.txt")
SRA_RUNS_FILE        <- file.path(METADATA_DIR, "sra_runs.csv")

ALPHA        <- 0.01
ALPHA_TABLE1 <- 0.05
N_STRATA     <- 20
N_REP        <- 20
BASE_SEED    <- 1L

## Expected sizes of the reference used in the manuscript.
N_TRANSCRIPTS       <- 198507L
N_ORPHAN_TRANSCRIPT <- 2190L
N_GENE_LOCI         <- 58962L
N_ORPHAN_LOCUS      <- 1226L

## Clinical analysis in the manuscript was based on the patient-level
## standardized-expression analysis with zero-abundance transcripts excluded.
## Change these two values only if a different score is intentionally analysed.
CLINICAL_SCORE_TABLE <- "table2"
CLINICAL_SCORE_ZERO  <- "excluded"
CLINICAL_P_COLUMN    <- "p_up"

COHORTS <- list(
  GSE244679 = list(
    disease   = "Psoriasis",
    n_pairs   = 24L,
    sheet     = file.path(METADATA_DIR, "sample_sheet_GSE244679.csv"),
    clin_type = "bayes_pasi"
  ),
  GSE127165 = list(
    disease   = "Laryngeal squamous cell carcinoma",
    n_pairs   = 57L,
    sheet     = file.path(METADATA_DIR, "sample_sheet_GSE127165.csv"),
    clin_type = "lm_stage_age_alcohol"
  ),
  GSE144269 = list(
    disease   = "Hepatocellular carcinoma",
    n_pairs   = 70L,
    sheet     = file.path(METADATA_DIR, "sample_sheet_GSE144269.csv"),
    clin_type = "none"
  ),
  GSE40419 = list(
    disease   = "Lung adenocarcinoma",
    n_pairs   = 69L,
    sheet     = file.path(METADATA_DIR, "sample_sheet_GSE40419.csv"),
    clin_type = "lm_stage_age_smoking"
  )
)

read_sheet <- function(cohort) {
  cfg <- COHORTS[[cohort]]
  if (is.null(cfg)) stop("unknown cohort: ", cohort)
  if (!file.exists(cfg$sheet)) stop("sample sheet not found: ", cfg$sheet)

  s <- read.csv(cfg$sheet, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("library_id", "patient_id", "condition", "directory")
  missing <- setdiff(required, colnames(s))
  if (length(missing)) {
    stop("sample sheet is missing columns: ", paste(missing, collapse = ", "))
  }
  if (!all(s$condition %in% c("normal", "disease"))) {
    stop("condition must be 'normal' or 'disease' in ", cfg$sheet)
  }

  tab <- table(s$patient_id, s$condition)
  if (ncol(tab) != 2L || !all(c("normal", "disease") %in% colnames(tab)) ||
      !all(tab[, c("normal", "disease"), drop = FALSE] == 1L)) {
    stop("sample sheet is not a one-normal/one-disease matched-pair design: ", cfg$sheet)
  }
  if (nrow(tab) != cfg$n_pairs) {
    stop("expected ", cfg$n_pairs, " matched pairs in ", cohort,
         " but found ", nrow(tab))
  }
  s[order(s$patient_id, match(s$condition, c("normal", "disease"))), , drop = FALSE]
}

## Table 1 uses every screened RNA-seq library, not only matched pairs.
## The full library list is taken from metadata/sra_runs.csv. If that file is
## absent, the matched-pair sample sheet is used as a fallback with a warning.
read_library_sheet <- function(cohort) {
  if (file.exists(SRA_RUNS_FILE)) {
    z <- read.csv(SRA_RUNS_FILE, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("run_accession", "cohort", "library_layout")
    missing <- setdiff(required, colnames(z))
    if (length(missing)) {
      stop("sra_runs.csv is missing columns: ", paste(missing, collapse = ", "))
    }
    z <- z[z$cohort == cohort, , drop = FALSE]
    if (nrow(z) > 0L) {
      return(data.frame(
        library_id = z$run_accession,
        directory  = paste0(cohort, "_", z$run_accession),
        stringsAsFactors = FALSE
      ))
    }
  }

  warning("no full run list found for ", cohort,
          "; Table 1 will use only the matched-pair libraries")
  s <- read_sheet(cohort)
  unique(s[, c("library_id", "directory")])
}
