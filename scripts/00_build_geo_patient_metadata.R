#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## 00_build_geo_patient_metadata.R
##
## Fully automatic reconstruction of the public-data manifest used by the
## manuscript.  No hand-written run list, normal/disease table, or clinical
## override file is required.
##
## Public sources used by the code
##   * GEO series/sample metadata (GEOquery)
##   * ENA Portal read_run metadata and direct FASTQ URLs
##   * GSE127165 primary-paper Supplementary Table S1 for patient age
##
## Outputs
##   metadata/sra_runs.csv
##   metadata/sample_sheet_GSE244679.csv
##   metadata/sample_sheet_GSE127165.csv
##   metadata/sample_sheet_GSE144269.csv
##   metadata/sample_sheet_GSE40419.csv
##   metadata_output/GEO_patient_metadata_combined.tsv
##   metadata_output/GEO_sample_metadata_long.tsv
##   metadata_output/GEO_SRA_run_mapping.tsv
##   metadata_output/FASTQ_selection_audit.tsv
##   metadata_output/ERP001058_read_run.tsv
##   metadata_output/GSE40419_excluded_libraries.tsv
##   metadata_output/GSE127165_supplement_clinical.tsv
##   metadata_output/public_source_provenance.tsv
##
## Cohort-specific public-metadata rules encoded below
##   GSE244679: [patientN] = adjacent normal, [patientD] = psoriasis lesion
##   GSE127165: ANM_<subject> = normal, LSCC_<subject> = disease
##   GSE144269: "pat <n> non-tumor" = normal, "pat <n> tumor" = disease
##   GSE40419: LC_*_nor = normal and matching LC_* = disease.  The manuscript
##               screened the 162 ERP001058 runs having direct ENA FASTQ URLs;
##               69 disease + 76 normal runs mapped to GEO titles and 17 runs
##               were unclassified, yielding 69 complete pairs and 24 excluded
##               unmatched/unclassified libraries.
##
## The script deliberately stops if current public metadata no longer reproduces
## the manuscript counts.  This is safer than silently choosing a different
## cohort after an archive metadata change.
## ---------------------------------------------------------------------------

required_bioc <- c("GEOquery", "Biobase")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  stop(
    "Missing Bioconductor package(s): ", paste(missing_bioc, collapse = ", "), "\n",
    "Install with BiocManager::install(c('GEOquery','Biobase'))."
  )
}
required_cran <- c("jsonlite", "readxl", "xml2")
missing_cran <- required_cran[!vapply(required_cran, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  stop(
    "Missing CRAN package(s): ", paste(missing_cran, collapse = ", "), "\n",
    "Install with install.packages(c('jsonlite','readxl','xml2'))."
  )
}

OUTPUT_DIR <- Sys.getenv("METADATA_OUTPUT_DIR", unset = "metadata_output")
METADATA_DIR <- Sys.getenv("METADATA_DIR", unset = "metadata")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADATA_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FILE <- file.path(OUTPUT_DIR, "GEO_patient_metadata_combined.tsv")
SAMPLE_LONG_FILE <- file.path(OUTPUT_DIR, "GEO_sample_metadata_long.tsv")
SRA_MAPPING_FILE <- file.path(OUTPUT_DIR, "GEO_SRA_run_mapping.tsv")
FASTQ_AUDIT_FILE <- file.path(OUTPUT_DIR, "FASTQ_selection_audit.tsv")
SRA_RUNS_FILE <- file.path(METADATA_DIR, "sra_runs.csv")
GSE40419_ENA_RUN_TABLE <- file.path(OUTPUT_DIR, "ERP001058_read_run.tsv")
GSE40419_EXCLUDED_FILE <- file.path(OUTPUT_DIR, "GSE40419_excluded_libraries.tsv")
GSE127165_SUPP_FILE <- file.path(OUTPUT_DIR, "GSE127165_supplement_clinical.tsv")
PUBLIC_SOURCE_FILE <- file.path(OUTPUT_DIR, "public_source_provenance.tsv")
ORPHAN_DELTA_FILE <- Sys.getenv(
  "ORPHAN_DELTA_INPUT_FILE",
  unset = file.path(METADATA_DIR, "orphan_delta_input.tsv")
)

STRICT_MANUSCRIPT_COUNTS <- !identical(
  Sys.getenv("STRICT_MANUSCRIPT_COUNTS", unset = "1"), "0"
)

GSE_IDS <- c("GSE244679", "GSE127165", "GSE144269", "GSE40419")
EXPECTED_SCREENED_LIBRARIES <- c(
  GSE244679 = 48L,
  GSE127165 = 114L,
  GSE144269 = 140L,
  GSE40419  = 162L
)
EXPECTED_MATCHED_PAIRS <- c(
  GSE244679 = 24L,
  GSE127165 = 57L,
  GSE144269 = 70L,
  GSE40419  = 69L
)
EXPECTED_GSE40419_CLASS <- c(disease = 69L, normal = 76L, unknown = 17L)
EXPECTED_GEO_CLASSES <- list(
  GSE244679 = c(disease = 24L, normal = 24L),
  GSE127165 = c(disease = 57L, normal = 57L),
  GSE144269 = c(disease = 70L, normal = 70L),
  GSE40419  = c(disease = 87L, normal = 77L)
)
GSE40419_ENA_STUDY <- "ERP001058"
GSE127165_FIGSHARE_ARTICLE <- "12414305"
GSE127165_SPRINGER_ZIP <- paste0(
  "https://media.springernature.com/original/springer-static/esm/",
  "art%3A10.1186%2Fs12943-020-01215-4/MediaObjects/",
  "12943_2020_1215_MOESM1_ESM.zip"
)

DISEASE_NAMES <- c(
  GSE244679 = "Psoriasis",
  GSE127165 = "Laryngeal squamous cell carcinoma",
  GSE144269 = "Hepatocellular carcinoma",
  GSE40419  = "Lung adenocarcinoma"
)

## ---------------------------------------------------------------------------
## Generic helpers
## ---------------------------------------------------------------------------
trim_na <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | tolower(x) %in% c(
    "", "na", "n/a", "not available", "not applicable", "unknown", "none", "-"
  )
  x[bad] <- NA_character_
  x
}

safe_numeric <- function(x) {
  x <- trim_na(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    m <- regexpr("-?[0-9]+(?:\\.[0-9]+)?", s, perl = TRUE)
    if (m[1] < 0) return(NA_real_)
    suppressWarnings(as.numeric(regmatches(s, m)))
  }, numeric(1))
}

sanitize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

normalize_run <- function(x) {
  x <- trim_na(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    m <- regexpr("(?:SRR|ERR|DRR)[0-9]+", s, perl = TRUE, ignore.case = TRUE)
    if (m[1] < 0) return(NA_character_)
    toupper(regmatches(s, m))
  }, character(1))
}

extract_first <- function(x, pattern) {
  x <- trim_na(x)
  if (!length(x) || all(is.na(x))) return(NA_character_)
  s <- paste(x[!is.na(x)], collapse = " ; ")
  m <- regexpr(pattern, s, perl = TRUE, ignore.case = TRUE)
  if (m[1] < 0) return(NA_character_)
  toupper(regmatches(s, m))
}

split_semicolon <- function(x) {
  x <- trim_na(x)
  if (length(x) == 0L || is.na(x[1])) return(character())
  z <- trimws(strsplit(x[1], ";", fixed = TRUE)[[1]])
  z[nzchar(z)]
}

httpsify_ena <- function(x) {
  x <- trim_na(x)
  ifelse(
    is.na(x), NA_character_,
    ifelse(grepl("^https?://", x, ignore.case = TRUE), x, paste0("https://", x))
  )
}

flatten_characteristics <- function(pd) {
  char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE,
                    ignore.case = TRUE)
  if (!length(char_cols)) return(pd)

  kv_per_row <- vector("list", nrow(pd))
  keys <- character()
  for (i in seq_len(nrow(pd))) {
    vals <- unlist(pd[i, char_cols, drop = FALSE], use.names = FALSE)
    vals <- trim_na(vals)
    vals <- vals[!is.na(vals)]
    kv <- list()
    for (v in vals) {
      if (!grepl(":", v, fixed = TRUE)) next
      key <- sanitize_key(sub(":.*$", "", v))
      val <- trim_na(sub("^[^:]*:\\s*", "", v))
      if (!nzchar(key) || is.na(val)) next
      if (is.null(kv[[key]]) || is.na(kv[[key]])) kv[[key]] <- val
      keys <- union(keys, key)
    }
    kv_per_row[[i]] <- kv
  }
  for (key in keys) {
    pd[[key]] <- vapply(kv_per_row, function(kv) {
      val <- kv[[key]]
      if (is.null(val)) NA_character_ else as.character(val)
    }, character(1))
  }
  pd
}

row_first <- function(df, names) {
  names <- intersect(names, colnames(df))
  out <- rep(NA_character_, nrow(df))
  for (nm in names) {
    z <- trim_na(df[[nm]])
    use <- is.na(out) & !is.na(z)
    out[use] <- z[use]
  }
  out
}

write_tsv <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

assert_equal_count <- function(label, observed, expected) {
  if (STRICT_MANUSCRIPT_COUNTS && as.integer(observed) != as.integer(expected)) {
    stop(label, ": expected ", expected, " but reconstructed ", observed,
         ". Public archive metadata may have changed; inspect the audit files rather ",
         "than silently changing the manuscript cohort.")
  }
}

## ---------------------------------------------------------------------------
## GEO metadata
## ---------------------------------------------------------------------------
fetch_geo_samples <- function(gse) {
  message("Fetching GEO metadata: ", gse)
  obj <- GEOquery::getGEO(gse, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
  if (is.list(obj)) {
    if (length(obj) != 1L) {
      stop(gse, ": expected one expression platform but GEOquery returned ", length(obj))
    }
    obj <- obj[[1]]
  }
  pd <- Biobase::pData(obj)
  pd <- flatten_characteristics(pd)
  if (!"geo_accession" %in% colnames(pd)) stop(gse, ": GEO metadata lacks geo_accession")
  if (!"title" %in% colnames(pd)) stop(gse, ": GEO metadata lacks title")

  rel_cols <- grep("^relation", colnames(pd), value = TRUE, ignore.case = TRUE)
  pd$sra_experiment <- vapply(seq_len(nrow(pd)), function(i) {
    vals <- if (length(rel_cols)) unlist(pd[i, rel_cols, drop = FALSE], use.names = FALSE) else character()
    extract_first(vals, "(?:SRX|ERX|DRX)[0-9]+")
  }, character(1))
  pd$biosample <- vapply(seq_len(nrow(pd)), function(i) {
    vals <- if (length(rel_cols)) unlist(pd[i, rel_cols, drop = FALSE], use.names = FALSE) else character()
    extract_first(vals, "(?:SAMN|SAMEA|SAMD)[0-9]+")
  }, character(1))
  pd$direct_run <- vapply(seq_len(nrow(pd)), function(i) {
    vals <- if (length(rel_cols)) unlist(pd[i, rel_cols, drop = FALSE], use.names = FALSE) else character()
    extract_first(vals, "(?:SRR|ERR|DRR)[0-9]+")
  }, character(1))

  pd$cohort <- gse
  pd$gsm <- as.character(pd$geo_accession)
  pd$source_name <- if ("source_name_ch1" %in% colnames(pd)) {
    as.character(pd$source_name_ch1)
  } else {
    rep(NA_character_, nrow(pd))
  }
  front <- c("cohort", "gsm", "title", "source_name")
  pd[, c(front, setdiff(colnames(pd), front)), drop = FALSE]
}

annotate_cohort <- function(pd, gse) {
  title <- trimws(as.character(pd$title))
  patient <- rep(NA_character_, nrow(pd))
  condition <- rep(NA_character_, nrow(pd))

  if (gse == "GSE244679") {
    m <- regexec("\\[([0-9]+)([ND])\\]", title, perl = TRUE, ignore.case = TRUE)
    mm <- regmatches(title, m)
    ok <- lengths(mm) == 3L
    patient[ok] <- vapply(mm[ok], `[[`, character(1), 2L)
    code <- rep(NA_character_, length(mm))
    code[ok] <- toupper(vapply(mm[ok], `[[`, character(1), 3L))
    condition[code == "N"] <- "normal"
    condition[code == "D"] <- "disease"
  } else if (gse == "GSE127165") {
    m <- regexec("^(LSCC|ANM)_([0-9]+)$", title, perl = TRUE, ignore.case = TRUE)
    mm <- regmatches(title, m)
    ok <- lengths(mm) == 3L
    typ <- rep(NA_character_, length(mm))
    typ[ok] <- toupper(vapply(mm[ok], `[[`, character(1), 2L))
    patient[ok] <- vapply(mm[ok], `[[`, character(1), 3L)
    condition[typ == "LSCC"] <- "disease"
    condition[typ == "ANM"] <- "normal"
    if ("subject_number" %in% colnames(pd)) {
      s <- trim_na(pd$subject_number)
      use <- !is.na(s)
      patient[use] <- sub(".*?([0-9]+).*", "\\1", s[use])
    }
  } else if (gse == "GSE144269") {
    m <- regexec("^pat\\s+([0-9]+)\\s+(non-tumor|tumor)", title,
                 perl = TRUE, ignore.case = TRUE)
    mm <- regmatches(title, m)
    ok <- lengths(mm) == 3L
    patient[ok] <- vapply(mm[ok], `[[`, character(1), 2L)
    typ <- rep(NA_character_, length(mm))
    typ[ok] <- tolower(vapply(mm[ok], `[[`, character(1), 3L))
    condition[typ == "tumor"] <- "disease"
    condition[typ == "non-tumor"] <- "normal"
  } else if (gse == "GSE40419") {
    is_lc <- grepl("^LC_[A-Za-z]+[0-9]+(?:_nor)?$", title, perl = TRUE)
    patient[is_lc] <- sub("_nor$", "", title[is_lc], ignore.case = TRUE)
    condition[is_lc & grepl("_nor$", title, ignore.case = TRUE)] <- "normal"
    condition[is_lc & !grepl("_nor$", title, ignore.case = TRUE)] <- "disease"
  } else {
    stop("No cohort annotation rule for ", gse)
  }

  pd$patient_id <- patient
  pd$condition <- condition

  expected <- EXPECTED_GEO_CLASSES[[gse]]
  if (!is.null(expected)) {
    observed <- c(
      disease = sum(pd$condition == "disease", na.rm = TRUE),
      normal = sum(pd$condition == "normal", na.rm = TRUE)
    )
    if (STRICT_MANUSCRIPT_COUNTS &&
        (any(observed != expected) || anyNA(pd$condition) || anyNA(pd$patient_id))) {
      stop(
        gse, ": public GEO title/characteristic parsing failed. Expected ",
        paste(names(expected), expected, sep = "=", collapse = ", "),
        "; observed ", paste(names(observed), observed, sep = "=", collapse = ", "),
        "; unclassified GEO samples=", sum(is.na(pd$condition) | is.na(pd$patient_id)),
        ". Refusing to invent a condition or patient pairing."
      )
    }
  }
  pd
}

## ---------------------------------------------------------------------------
## ENA read_run metadata
## ---------------------------------------------------------------------------
ENA_FIELDS <- paste(
  c("run_accession", "experiment_accession", "experiment_alias",
    "sample_accession", "secondary_sample_accession", "sample_alias", "sample_title",
    "library_layout", "fastq_ftp", "fastq_md5", "fastq_bytes"),
  collapse = ","
)

fetch_ena_read_run <- function(accession) {
  url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
    utils::URLencode(accession, reserved = TRUE),
    "&result=read_run&fields=", ENA_FIELDS,
    "&format=tsv&download=false"
  )
  tab <- tryCatch(
    utils::read.delim(url, stringsAsFactors = FALSE, check.names = FALSE,
                      quote = "", comment.char = ""),
    error = function(e) stop("ENA query failed for ", accession, ": ", conditionMessage(e))
  )
  if (!nrow(tab)) return(tab)
  tab$run_accession <- normalize_run(tab$run_accession)
  tab$library_layout <- toupper(trim_na(tab$library_layout))
  tab
}

expand_fastq_fields <- function(tab) {
  if (!nrow(tab)) return(tab)
  tab$fastq_1 <- NA_character_
  tab$fastq_2 <- NA_character_
  tab$md5_1 <- NA_character_
  tab$md5_2 <- NA_character_
  tab$bytes_1 <- NA_character_
  tab$bytes_2 <- NA_character_
  tab$usable_direct_fastq <- FALSE

  for (i in seq_len(nrow(tab))) {
    fq <- httpsify_ena(split_semicolon(tab$fastq_ftp[i]))
    md <- split_semicolon(tab$fastq_md5[i])
    by <- split_semicolon(tab$fastq_bytes[i])
    layout <- toupper(trim_na(tab$library_layout[i]))

    if (!length(fq)) next
    if (identical(layout, "PAIRED") && length(fq) == 2L) {
      ## ENA normally returns mate 1 then mate 2.  If suffixes are explicit,
      ## enforce that order and apply the same permutation to MD5/byte fields.
      b <- basename(fq)
      mate <- ifelse(grepl("_1\\.f(ast)?q\\.gz$", b, ignore.case = TRUE), 1L,
              ifelse(grepl("_2\\.f(ast)?q\\.gz$", b, ignore.case = TRUE), 2L, NA_integer_))
      if (setequal(mate, c(1L, 2L))) {
        ord <- order(mate)
        fq <- fq[ord]
        if (length(md) == 2L) md <- md[ord]
        if (length(by) == 2L) by <- by[ord]
      }
      tab$fastq_1[i] <- fq[1]
      tab$fastq_2[i] <- fq[2]
      if (length(md) >= 1L) tab$md5_1[i] <- md[1]
      if (length(md) >= 2L) tab$md5_2[i] <- md[2]
      if (length(by) >= 1L) tab$bytes_1[i] <- by[1]
      if (length(by) >= 2L) tab$bytes_2[i] <- by[2]
      tab$usable_direct_fastq[i] <- TRUE
    } else if (identical(layout, "SINGLE") && length(fq) == 1L) {
      tab$fastq_1[i] <- fq[1]
      if (length(md) >= 1L) tab$md5_1[i] <- md[1]
      if (length(by) >= 1L) tab$bytes_1[i] <- by[1]
      tab$usable_direct_fastq[i] <- TRUE
    }
  }
  tab
}

resolve_geo_library <- function(row) {
  query <- trim_na(row$sra_experiment)
  if (is.na(query)) query <- trim_na(row$direct_run)
  if (is.na(query)) query <- trim_na(row$biosample)
  if (is.na(query)) {
    stop(row$cohort, " ", row$gsm, " (", row$title, "): no SRA/ENA relation in GEO")
  }

  tab <- expand_fastq_fields(fetch_ena_read_run(query))
  if (!nrow(tab)) stop("ENA returned no read_run row for ", row$gsm, " via ", query)
  tab <- tab[tab$usable_direct_fastq, , drop = FALSE]
  if (!nrow(tab)) {
    stop(row$cohort, " ", row$gsm, " (", row$title,
         "): no direct FASTQ URL in the GEO-linked ENA record")
  }

  ## Normally one GEO RNA-seq library corresponds to one run in these three
  ## cohorts.  If ENA ever exposes several runs for the same experiment, do not
  ## make a hidden arbitrary choice: current public metadata must still identify
  ## one direct-FASTQ run for manuscript reproduction.
  if (nrow(tab) != 1L) {
    ## Exact sample-alias/title agreement is a defensible archive-level tie break.
    alias_match <- (!is.na(trim_na(tab$sample_alias)) & trim_na(tab$sample_alias) == row$title) |
                   (!is.na(trim_na(tab$sample_title)) & trim_na(tab$sample_title) == row$title)
    tab2 <- tab[alias_match, , drop = FALSE]
    if (nrow(tab2) == 1L) tab <- tab2
  }
  if (nrow(tab) != 1L) {
    stop(row$cohort, " ", row$gsm, " (", row$title, ") maps to ", nrow(tab),
         " usable ENA runs: ", paste(tab$run_accession, collapse = ", "),
         ". The public metadata no longer define the manuscript library uniquely.")
  }

  data.frame(
    run_accession = tab$run_accession[1],
    cohort = row$cohort,
    gsm = row$gsm,
    patient_id = row$patient_id,
    condition = row$condition,
    title = row$title,
    sra_experiment = query,
    experiment_accession = tab$experiment_accession[1],
    sample_alias = trim_na(tab$sample_alias[1]),
    library_layout = tab$library_layout[1],
    fastq_1 = tab$fastq_1[1],
    fastq_2 = tab$fastq_2[1],
    md5_1 = tab$md5_1[1],
    md5_2 = tab$md5_2[1],
    bytes_1 = tab$bytes_1[1],
    bytes_2 = tab$bytes_2[1],
    mapping_status = "GEO-linked experiment",
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Special reconstruction of the 162 screened GSE40419 libraries.
##
## The original public-data workflow started from the complete ERP001058
## read_run report and retained rows having direct ENA FASTQ URLs.  This gives
## the manuscript's 162 screened libraries.  sample_alias is then matched to
## the 164 GEO titles.  The expected mapping is 69 disease, 76 normal, 17
## unclassified; the 69 disease runs all have a matched normal run.
## ---------------------------------------------------------------------------
resolve_gse40419 <- function(geo40419) {
  message("Fetching complete ENA study: ", GSE40419_ENA_STUDY)
  all <- expand_fastq_fields(fetch_ena_read_run(GSE40419_ENA_STUDY))
  write_tsv(all, GSE40419_ENA_RUN_TABLE)

  selected <- all[all$usable_direct_fastq, , drop = FALSE]
  assert_equal_count("GSE40419 direct-FASTQ screened libraries",
                     nrow(selected), EXPECTED_SCREENED_LIBRARIES[["GSE40419"]])

  alias <- trim_na(selected$sample_alias)
  title <- trimws(as.character(geo40419$title))
  idx <- match(alias, title)

  selected$cohort <- "GSE40419"
  selected$gsm <- geo40419$gsm[idx]
  selected$patient_id <- geo40419$patient_id[idx]
  selected$condition <- geo40419$condition[idx]
  selected$title <- geo40419$title[idx]
  selected$mapping_status <- ifelse(is.na(idx), "unclassified ENA alias", "matched GEO title")
  selected$sra_experiment <- selected$experiment_accession

  cls <- c(
    disease = sum(selected$condition == "disease", na.rm = TRUE),
    normal = sum(selected$condition == "normal", na.rm = TRUE),
    unknown = sum(is.na(selected$condition))
  )
  if (STRICT_MANUSCRIPT_COUNTS && any(cls != EXPECTED_GSE40419_CLASS)) {
    stop(
      "GSE40419 classification no longer reproduces the manuscript: expected ",
      paste(names(EXPECTED_GSE40419_CLASS), EXPECTED_GSE40419_CLASS,
            sep = "=", collapse = ", "), "; observed ",
      paste(names(cls), cls, sep = "=", collapse = ", ")
    )
  }

  out <- data.frame(
    run_accession = selected$run_accession,
    cohort = selected$cohort,
    gsm = selected$gsm,
    patient_id = selected$patient_id,
    condition = selected$condition,
    title = selected$title,
    sra_experiment = selected$sra_experiment,
    experiment_accession = selected$experiment_accession,
    sample_alias = trim_na(selected$sample_alias),
    library_layout = selected$library_layout,
    fastq_1 = selected$fastq_1,
    fastq_2 = selected$fastq_2,
    md5_1 = selected$md5_1,
    md5_2 = selected$md5_2,
    bytes_1 = selected$bytes_1,
    bytes_2 = selected$bytes_2,
    mapping_status = selected$mapping_status,
    stringsAsFactors = FALSE
  )

  ## The manuscript excludes all unclassified rows and the seven mapped normal
  ## libraries whose patient has no selected disease library.
  disease_ids <- unique(out$patient_id[out$condition == "disease" & !is.na(out$patient_id)])
  normal_ids <- unique(out$patient_id[out$condition == "normal" & !is.na(out$patient_id)])
  paired_ids <- intersect(disease_ids, normal_ids)
  assert_equal_count("GSE40419 matched pairs", length(paired_ids),
                     EXPECTED_MATCHED_PAIRS[["GSE40419"]])

  bad <- is.na(out$condition) | is.na(out$patient_id) | !out$patient_id %in% paired_ids
  excluded <- out[bad, , drop = FALSE]
  excluded$exclusion_reason <- ifelse(
    is.na(excluded$condition), "ENA run alias not mapped to a GEO disease/normal title",
    "mapped GEO library without a selected matched counterpart"
  )
  assert_equal_count("GSE40419 excluded unmatched/unclassified libraries",
                     nrow(excluded), 24L)
  write_tsv(excluded, GSE40419_EXCLUDED_FILE)

  list(runs = out, paired_ids = paired_ids)
}

## ---------------------------------------------------------------------------
## Clinical metadata helpers
## ---------------------------------------------------------------------------
clinical_from_geo <- function(pd) {
  data.frame(
    PASI = safe_numeric(row_first(pd, c("pasi_score", "pasi"))),
    age = safe_numeric(row_first(pd, c("age", "age_at_diagnosis", "age_at_diagnosis_years"))),
    age_of_onset = safe_numeric(row_first(pd, c("age_of_onset", "age_onset"))),
    stage = trim_na(row_first(pd, c("tumor_stage", "tumour_stage", "stage"))),
    alcohol = trim_na(row_first(pd, c("alcohol_consumption", "alcohol"))),
    smoking = trim_na(row_first(pd, c("tobacco_smoking", "smoking_status", "smoking"))),
    sex = trim_na(row_first(pd, c("sex", "gender"))),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## GSE127165 age is not present in the GEO sample characteristics.  It is part
## of Additional file 1, Table S1 of the primary Molecular Cancer paper.  The
## code downloads the public Figshare/Springer supplement and automatically
## identifies the Table S1-like table by matching its subject-number column to
## the 57 GEO subject numbers and requiring an age column.
## ---------------------------------------------------------------------------
extract_integer_token <- function(x) {
  x <- trim_na(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    m <- regexpr("[0-9]+", s, perl = TRUE)
    if (m[1] < 0) return(NA_character_)
    regmatches(s, m)
  }, character(1))
}

rectangularize <- function(x) {
  if (!length(x)) return(NULL)
  n <- max(lengths(x))
  m <- matrix(NA_character_, nrow = length(x), ncol = n)
  for (i in seq_along(x)) {
    if (length(x[[i]])) m[i, seq_along(x[[i]])] <- as.character(x[[i]])
  }
  as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
}

read_docx_tables <- function(file) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("The GSE127165 supplement is DOCX; install.packages('xml2') is required to parse it.")
  }
  td <- tempfile("docx_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(file, files = "word/document.xml", exdir = td)
  xmlfile <- file.path(td, "word", "document.xml")
  if (!file.exists(xmlfile)) return(list())
  doc <- xml2::read_xml(xmlfile)
  ## Use local-name() rather than assuming the DOCX namespace prefix is
  ## literally "w"; Word-processing XML can bind the same namespace under a
  ## different prefix without changing the document.
  tbls <- xml2::xml_find_all(doc, ".//*[local-name()='tbl']")
  lapply(tbls, function(tbl) {
    trs <- xml2::xml_find_all(tbl, "./*[local-name()='tr']")
    rows <- lapply(trs, function(tr) {
      cells <- xml2::xml_find_all(tr, "./*[local-name()='tc']")
      vapply(cells, function(tc) {
        tx <- xml2::xml_text(xml2::xml_find_all(tc, ".//*[local-name()='t']"))
        trimws(paste(tx, collapse = " "))
      }, character(1))
    })
    rectangularize(rows)
  })
}

read_supplement_tables <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(file)
    return(lapply(sheets, function(s) {
      as.data.frame(readxl::read_excel(file, sheet = s, col_names = FALSE),
                    stringsAsFactors = FALSE, check.names = FALSE)
    }))
  }
  if (ext == "csv") {
    return(list(read.csv(file, header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE)))
  }
  if (ext %in% c("tsv", "txt")) {
    return(list(read.delim(file, header = FALSE, stringsAsFactors = FALSE,
                           check.names = FALSE)))
  }
  if (ext == "docx") return(read_docx_tables(file))
  list()
}

download_gse127165_supplement <- function() {
  td <- tempfile("gse127165_supp_")
  dir.create(td)
  api <- paste0("https://api.figshare.com/v2/articles/", GSE127165_FIGSHARE_ARTICLE)
  meta <- tryCatch(jsonlite::fromJSON(api), error = function(e) NULL)
  files <- character()

  if (!is.null(meta) && !is.null(meta$files) && nrow(meta$files)) {
    for (i in seq_len(nrow(meta$files))) {
      nm <- as.character(meta$files$name[i])
      url <- as.character(meta$files$download_url[i])
      if (is.na(url) || !nzchar(url)) next
      dst <- file.path(td, nm)
      ok <- tryCatch({
        utils::download.file(url, dst, mode = "wb", quiet = TRUE)
        file.exists(dst) && file.info(dst)$size > 0
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (ok) files <- c(files, dst)
    }
  }

  if (!length(files)) {
    dst <- file.path(td, "Additional_file_1.zip")
    ok <- tryCatch({
      utils::download.file(GSE127165_SPRINGER_ZIP, dst, mode = "wb", quiet = TRUE)
      file.exists(dst) && file.info(dst)$size > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) stop("Could not download GSE127165 Additional file 1 from Figshare or Springer")
    files <- dst
  }

  expanded <- character()
  for (f in files) {
    if (tolower(tools::file_ext(f)) == "zip") {
      ex <- file.path(td, paste0("unzipped_", basename(f)))
      dir.create(ex)
      utils::unzip(f, exdir = ex)
      expanded <- c(expanded, list.files(ex, recursive = TRUE, full.names = TRUE))
    } else {
      expanded <- c(expanded, f)
    }
  }
  expanded[file.info(expanded)$isdir %in% FALSE]
}

find_lscc_age_table <- function(files, known_ids) {
  candidates <- list()
  k <- 0L
  for (f in files) {
    tabs <- tryCatch(read_supplement_tables(f), error = function(e) list())
    if (!length(tabs)) next
    for (tab_i in seq_along(tabs)) {
      tab <- tabs[[tab_i]]
      if (is.null(tab) || nrow(tab) < 5L || ncol(tab) < 2L) next
      tab[] <- lapply(tab, as.character)

      ## Locate a header row containing an age label.  Additional file 1 contains
      ## both Table S1 (57 RNA-sequenced patients) and Table S2 (107 qPCR
      ## patients), so candidate ranking below also considers the total number of
      ## patient-like rows; this prevents Table S2 from being selected merely
      ## because some of its subject numbers overlap the 57 GEO IDs.
      for (hr in seq_len(min(20L, nrow(tab)))) {
        headers <- trimws(as.character(unlist(tab[hr, , drop = FALSE], use.names = FALSE)))
        hnorm <- sanitize_key(headers)
        age_cols <- which(grepl("(^|_)age($|_)", hnorm))
        if (!length(age_cols)) next

        body <- tab[(hr + 1L):nrow(tab), , drop = FALSE]
        if (!nrow(body)) next
        for (idcol in seq_len(ncol(body))) {
          ids <- extract_integer_token(body[[idcol]])
          nmatch <- length(unique(ids[!is.na(ids) & ids %in% known_ids]))
          if (nmatch < 20L) next

          for (agecol in age_cols) {
            ## Never allow the age column itself to masquerade as the subject-ID
            ## column simply because some ages equal GEO subject numbers.
            if (idcol == agecol) next
            ages <- safe_numeric(body[[agecol]])
            valid_age <- is.finite(ages) & ages >= 18 & ages <= 100
            all_patient_like <- !is.na(ids) & valid_age
            ok <- all_patient_like & ids %in% known_ids
            score <- length(unique(ids[ok]))
            if (score < 20L) next

            total_subjects <- length(unique(ids[all_patient_like]))
            k <- k + 1L
            candidates[[k]] <- data.frame(
              patient_id = ids[ok],
              age = ages[ok],
              source_file = basename(f),
              source_table = tab_i,
              stringsAsFactors = FALSE
            )
            attr(candidates[[k]], "score") <- score
            attr(candidates[[k]], "total_subjects") <- total_subjects
            attr(candidates[[k]], "distance_from_57") <- abs(total_subjects - length(unique(known_ids)))
          }
        }
      }
    }
  }

  if (!length(candidates)) {
    stop("Could not identify GSE127165 Table S1 patient ages in Additional file 1")
  }

  ## Primary key: coverage of the known 57 GEO subjects (higher is better).
  ## Secondary key: total table size closest to 57 (Table S1 rather than S2).
  ## Tertiary key: earlier table in the supplementary file.
  score <- vapply(candidates, function(x) attr(x, "score"), numeric(1))
  dist57 <- vapply(candidates, function(x) attr(x, "distance_from_57"), numeric(1))
  tabno <- vapply(candidates, function(x) unique(x$source_table)[1], numeric(1))
  ord <- order(-score, dist57, tabno)
  out <- candidates[[ord[1]]]
  out <- out[!duplicated(out$patient_id), , drop = FALSE]
  out <- out[out$patient_id %in% known_ids, , drop = FALSE]

  if (STRICT_MANUSCRIPT_COUNTS && length(unique(out$patient_id)) != length(unique(known_ids))) {
    stop("GSE127165 supplement parser recovered ", length(unique(out$patient_id)),
         " of ", length(unique(known_ids)),
         " RNA-sequenced patient ages; refusing to guess any missing value")
  }
  out
}

## ---------------------------------------------------------------------------
## Pairing and final patient table
## ---------------------------------------------------------------------------
pair_run_manifest <- function(runmap, geo, gse, paired_ids = NULL) {
  z <- runmap[runmap$cohort == gse & !is.na(runmap$condition) &
              !is.na(runmap$patient_id), , drop = FALSE]
  if (!is.null(paired_ids)) z <- z[z$patient_id %in% paired_ids, , drop = FALSE]

  patients <- intersect(unique(z$patient_id[z$condition == "normal"]),
                        unique(z$patient_id[z$condition == "disease"]))
  patients <- sort(patients)
  assert_equal_count(paste0(gse, " matched pairs"), length(patients),
                     EXPECTED_MATCHED_PAIRS[[gse]])

  rows <- lapply(patients, function(pid) {
    nrowz <- z[z$patient_id == pid & z$condition == "normal", , drop = FALSE]
    drowz <- z[z$patient_id == pid & z$condition == "disease", , drop = FALSE]
    if (nrow(nrowz) != 1L || nrow(drowz) != 1L) {
      stop(gse, " patient ", pid, ": expected exactly one selected normal and disease run; got ",
           nrow(nrowz), " normal and ", nrow(drowz), " disease")
    }

    ng <- geo[match(nrowz$gsm[1], geo$gsm), , drop = FALSE]
    dg <- geo[match(drowz$gsm[1], geo$gsm), , drop = FALSE]
    if (!nrow(ng) || !nrow(dg)) stop(gse, " patient ", pid, ": GEO row missing")
    cn <- clinical_from_geo(ng)
    cd <- clinical_from_geo(dg)
    ## Patient covariates are taken from the disease sample first, with the
    ## matched normal metadata as fallback when the GEO record stores them there.
    pick <- function(a, b) ifelse(!is.na(a), a, b)

    data.frame(
      GSE = gse,
      cohort = gse,
      disease_name = DISEASE_NAMES[[gse]],
      patient_id = pid,
      normal_gsm = nrowz$gsm[1],
      disease_gsm = drowz$gsm[1],
      normal_title = nrowz$title[1],
      disease_title = drowz$title[1],
      normal_run_accession = nrowz$run_accession[1],
      disease_run_accession = drowz$run_accession[1],
      PASI = pick(cd$PASI, cn$PASI),
      age = pick(cd$age, cn$age),
      age_of_onset = pick(cd$age_of_onset, cn$age_of_onset),
      stage = pick(cd$stage, cn$stage),
      alcohol = pick(cd$alcohol, cn$alcohol),
      smoking = pick(cd$smoking, cn$smoking),
      sex = pick(cd$sex, cn$sex),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

join_orphan_delta <- function(combined) {
  combined$orphan_delta <- NA_real_
  if (!file.exists(ORPHAN_DELTA_FILE) || file.info(ORPHAN_DELTA_FILE)$size == 0) {
    return(combined)
  }
  d <- read.delim(ORPHAN_DELTA_FILE, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = c("NA", ""))
  if (!"orphan_delta" %in% colnames(d)) return(combined)
  if (!"cohort" %in% colnames(d) && "GSE" %in% colnames(d)) d$cohort <- d$GSE
  if (!"cohort" %in% colnames(d)) return(combined)

  if ("normal_run_accession" %in% colnames(d)) {
    key_d <- paste(d$cohort, normalize_run(d$normal_run_accession), sep = "::")
    key_c <- paste(combined$cohort, normalize_run(combined$normal_run_accession), sep = "::")
  } else if ("patient_id" %in% colnames(d)) {
    key_d <- paste(d$cohort, d$patient_id, sep = "::")
    key_c <- paste(combined$cohort, combined$patient_id, sep = "::")
  } else {
    return(combined)
  }
  idx <- match(key_c, key_d)
  combined$orphan_delta <- suppressWarnings(as.numeric(d$orphan_delta[idx]))
  combined
}

write_sample_sheets <- function(combined) {
  for (gse in GSE_IDS) {
    z <- combined[combined$cohort == gse, , drop = FALSE]
    assert_equal_count(paste0(gse, " sample-sheet pairs"), nrow(z),
                       EXPECTED_MATCHED_PAIRS[[gse]])
    ss <- do.call(rbind, lapply(seq_len(nrow(z)), function(i) {
      data.frame(
        library_id = c(z$normal_run_accession[i], z$disease_run_accession[i]),
        patient_id = c(z$patient_id[i], z$patient_id[i]),
        condition = c("normal", "disease"),
        directory = c(
          paste0(gse, "_", z$normal_run_accession[i]),
          paste0(gse, "_", z$disease_run_accession[i])
        ),
        stringsAsFactors = FALSE
      )
    }))
    write.csv(ss, file.path(METADATA_DIR, paste0("sample_sheet_", gse, ".csv")),
              row.names = FALSE, quote = FALSE)
  }
}

## ---------------------------------------------------------------------------
## Execute reconstruction
## ---------------------------------------------------------------------------
provenance <- data.frame(
  cohort = c(
    "GSE244679", "GSE127165", "GSE127165", "GSE127165",
    "GSE144269", "GSE40419", "GSE40419"
  ),
  source = c(
    "GEO series/sample metadata",
    "GEO series/sample metadata",
    "Primary publication",
    "Published supplementary dataset",
    "GEO series/sample metadata",
    "GEO series/sample metadata",
    "ENA read_run study metadata"
  ),
  accession_or_doi = c(
    "GSE244679", "GSE127165", "10.1186/s12943-020-01215-4",
    "10.6084/m9.figshare.12414305", "GSE144269", "GSE40419", "ERP001058"
  ),
  url = c(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244679",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127165",
    "https://link.springer.com/article/10.1186/s12943-020-01215-4",
    paste0("https://api.figshare.com/v2/articles/", GSE127165_FIGSHARE_ARTICLE),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419",
    paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
           GSE40419_ENA_STUDY, "&result=read_run")
  ),
  use = c(
    "24 N/D title-coded psoriasis pairs and SRA relations",
    "57 LSCC/ANM title-coded pairs; stage/smoking/alcohol and SRA relations",
    "Confirms 57 RNA-sequenced LSCC/ANM pairs and Additional file 1 Table S1",
    "Patient age for the 57 RNA-sequenced LSCC subjects",
    "70 tumor/non-tumor title-coded HCC pairs and SRA relations",
    "GEO LC_* / LC_*_nor titles and linkage to ERP001058",
    "Exact run accession, library layout, fastq_ftp and MD5 selection"
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, PUBLIC_SOURCE_FILE)

geo <- setNames(lapply(GSE_IDS, function(gse) annotate_cohort(fetch_geo_samples(gse), gse)),
                GSE_IDS)

## Save GEO-only audit table before run selection.
geo_long <- do.call(rbind, lapply(geo, function(z) {
  keep <- unique(c(
    "cohort", "gsm", "title", "source_name", "patient_id", "condition",
    "sra_experiment", "biosample", "direct_run",
    intersect(c("pasi_score", "age_of_onset", "age", "age_at_diagnosis",
                "tumor_stage", "tumour_stage", "stage", "alcohol_consumption",
                "tobacco_smoking", "smoking_status", "sex", "gender"), colnames(z))
  ))
  z[, keep, drop = FALSE]
}))
write_tsv(geo_long, SAMPLE_LONG_FILE)

## First three cohorts: every GEO library is linked to one direct-FASTQ run.
run_parts <- list()
for (gse in c("GSE244679", "GSE127165", "GSE144269")) {
  message("Resolving GEO-linked ENA FASTQs: ", gse)
  rows <- lapply(seq_len(nrow(geo[[gse]])), function(i) {
    if (i %% 20L == 0L || i == nrow(geo[[gse]])) {
      message("  ", gse, " ", i, "/", nrow(geo[[gse]]))
    }
    resolve_geo_library(geo[[gse]][i, , drop = FALSE])
  })
  z <- do.call(rbind, rows)
  assert_equal_count(paste0(gse, " screened libraries"), nrow(z),
                     EXPECTED_SCREENED_LIBRARIES[[gse]])
  run_parts[[gse]] <- z
}

## GSE40419: reproduce the historical ERP001058 direct-FASTQ screening.
r404 <- resolve_gse40419(geo[["GSE40419"]])
run_parts[["GSE40419"]] <- r404$runs

runmap <- do.call(rbind, run_parts)
rownames(runmap) <- NULL
runmap <- runmap[order(match(runmap$cohort, GSE_IDS)), , drop = FALSE]

## Validate total screened-library counts and write exact download manifest.
observed <- table(factor(runmap$cohort, levels = GSE_IDS))
for (gse in GSE_IDS) {
  assert_equal_count(paste0(gse, " screened libraries"), observed[[gse]],
                     EXPECTED_SCREENED_LIBRARIES[[gse]])
}

sra_runs <- runmap[, c(
  "run_accession", "cohort", "gsm", "patient_id", "condition",
  "library_layout", "fastq_1", "fastq_2", "md5_1", "md5_2"
), drop = FALSE]
write.csv(sra_runs, SRA_RUNS_FILE, row.names = FALSE, quote = FALSE, na = "")
write_tsv(runmap, SRA_MAPPING_FILE)
write_tsv(runmap[, c(
  "cohort", "run_accession", "gsm", "title", "sample_alias", "patient_id",
  "condition", "mapping_status", "library_layout", "fastq_1", "fastq_2",
  "md5_1", "md5_2"
), drop = FALSE], FASTQ_AUDIT_FILE)

## Build matched-patient table.
combined_parts <- list()
for (gse in c("GSE244679", "GSE127165", "GSE144269")) {
  combined_parts[[gse]] <- pair_run_manifest(runmap, geo[[gse]], gse)
}
combined_parts[["GSE40419"]] <- pair_run_manifest(
  runmap, geo[["GSE40419"]], "GSE40419", paired_ids = r404$paired_ids
)
combined <- do.call(rbind, combined_parts)
rownames(combined) <- NULL

## Add GSE127165 age from the published Additional file 1, Table S1.
message("Downloading/parsing GSE127165 Supplementary Table S1 for age")
known_lscc_ids <- unique(combined$patient_id[combined$cohort == "GSE127165"])
supp_files <- download_gse127165_supplement()
lscc_age <- find_lscc_age_table(supp_files, known_lscc_ids)
write_tsv(lscc_age, GSE127165_SUPP_FILE)
idx_lscc <- which(combined$cohort == "GSE127165")
age_idx <- match(combined$patient_id[idx_lscc], lscc_age$patient_id)
combined$age[idx_lscc] <- lscc_age$age[age_idx]
if (STRICT_MANUSCRIPT_COUNTS && anyNA(combined$age[idx_lscc])) {
  stop("GSE127165: at least one of the 57 matched patients lacks age after parsing Table S1")
}

combined <- join_orphan_delta(combined)

col_order <- c(
  "GSE", "cohort", "disease_name", "patient_id",
  "normal_gsm", "disease_gsm", "normal_title", "disease_title",
  "normal_run_accession", "disease_run_accession",
  "PASI", "age", "age_of_onset", "stage", "alcohol", "smoking", "sex",
  "orphan_delta"
)
combined <- combined[, col_order, drop = FALSE]
write_tsv(combined, OUTPUT_FILE)
write_sample_sheets(combined)

## Final clinical-coverage audit.
message("\nReconstructed matched pairs and clinical coverage:")
for (gse in GSE_IDS) {
  z <- combined[combined$cohort == gse, , drop = FALSE]
  cov <- vapply(c("PASI", "age", "age_of_onset", "stage", "alcohol", "smoking", "sex"),
                function(v) sum(!is.na(z[[v]])), integer(1))
  message("  ", gse, ": pairs=", nrow(z), "; ",
          paste(names(cov), cov, sep = "=", collapse = ", "))
}
message("\nWrote exact FASTQ manifest: ", SRA_RUNS_FILE)
message("Wrote combined clinical metadata: ", OUTPUT_FILE)
message("Wrote run/GEO audit: ", SRA_MAPPING_FILE)
message("No manual run-selection or clinical override file is used.")
