## ---------------------------------------------------------------------------
## 00_functions.R
## Core functions for the patient-level analysis of human orphan transcripts.
##
## Manuscript notation:
##   y_ij : raw TPM of transcript/locus i in sample j
##   x_ij : within-sample standardized expression
##   r_ij : within-sample expression rank
##   d_ij : normal - disease for the chosen expression representation
##
## Sign convention:
##   mean(d) < 0 -> higher expression in disease  ("Disease > normal")
##   mean(d) > 0 -> lower  expression in disease  ("Disease < normal")
## ---------------------------------------------------------------------------

## --- utilities --------------------------------------------------------------

assert_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x)) stop(name, " must be a non-missing scalar")
  invisible(TRUE)
}

safe_t_test <- function(x, alternative) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || length(unique(x)) < 2L || sd(x) == 0) return(NA_real_)
  t.test(x, mu = 0, alternative = alternative)$p.value
}

## --- reading abundance files ------------------------------------------------

read_abundance <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    z <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    z <- read.delim(file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
  required <- c("target_id", "tpm")
  missing <- setdiff(required, colnames(z))
  if (length(missing)) stop("missing columns in ", file, ": ", paste(missing, collapse = ", "))
  if (anyDuplicated(z$target_id)) stop("duplicated target_id values in ", file)
  if (any(!is.finite(z$tpm)) || any(z$tpm < 0)) stop("invalid TPM values in ", file)
  list(id = as.character(z$target_id), tpm = as.numeric(z$tpm))
}

## Resolve TSV/CSV abundance files transparently.  The manuscript analysis was
## originally run on kallisto abundance.tsv files, but the public code also
## accepts abundance.csv (and the corresponding gene-level CSV) without any
## change to the statistical analysis.
resolve_abundance_path <- function(base_dir, directory, abundance_file) {
  f <- file.path(base_dir, directory, abundance_file)
  if (file.exists(f)) return(f)

  stem <- sub("\\.(tsv|csv)$", "", abundance_file, ignore.case = TRUE)
  candidates <- unique(c(
    paste0(stem, ".tsv"),
    paste0(stem, ".csv"),
    ## tolerate the abbreviated filename used in some local working copies
    if (stem == "abundance") c("abundan.tsv", "abundan.csv") else character(0),
    if (stem == "abundance_gene") c("abundan_gene.tsv", "abndan_gene.tsv",
                                     "abundan_gene.csv", "abndan_gene.csv") else character(0)
  ))
  paths <- file.path(base_dir, directory, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 1L) return(hit)
  if (length(hit) > 1L) {
    stop("multiple candidate abundance files found in ", file.path(base_dir, directory),
         ": ", paste(basename(hit), collapse = ", "),
         ". Remove duplicates or specify a single canonical file.")
  }
  stop("missing abundance file: ", f,
       " (also tried: ", paste(candidates, collapse = ", "), ")")
}

## Load a collection of libraries into one TPM matrix. The caller supplies a
## data.frame with at least library_id and directory columns.
load_library_tpm <- function(library_sheet, base_dir, abundance_file = "abundance.tsv",
                             expected_features = NULL) {
  required <- c("library_id", "directory")
  missing <- setdiff(required, colnames(library_sheet))
  if (length(missing)) stop("library sheet is missing columns: ", paste(missing, collapse = ", "))

  n <- nrow(library_sheet)
  ref_id <- NULL
  mat <- NULL

  for (k in seq_len(n)) {
    f <- resolve_abundance_path(base_dir, library_sheet$directory[k], abundance_file)
    a <- read_abundance(f)

    if (is.null(ref_id)) {
      ref_id <- a$id
      if (!is.null(expected_features) && length(ref_id) != expected_features) {
        stop("reference contains ", length(ref_id), " features but ",
             expected_features, " were expected")
      }
      mat <- matrix(NA_real_, nrow = length(ref_id), ncol = n,
                    dimnames = list(ref_id, library_sheet$library_id))
    } else if (!identical(a$id, ref_id)) {
      stop("target_id order differs in ", f,
           ". All kallisto outputs must use the same reference/order.")
    }

    mat[, k] <- a$tpm
    if (k %% 20L == 0L || k == n) cat("  read ", k, "/", n, "\n", sep = "")
  }

  list(id = ref_id, tpm = mat)
}

load_cohort_tpm <- function(sample_sheet, base_dir, abundance_file = "abundance.tsv",
                            expected_features = NULL) {
  load_library_tpm(sample_sheet, base_dir, abundance_file, expected_features)
}

## --- orphan feature identification -----------------------------------------

## The original analysis used identifiers such as ORPHAN_hsa_00255579 and
## ORPHAN_XLOC_022250, while the raw Ruiz-Orera GTF may contain the same IDs
## without the ORPHAN_ prefix. Accept either representation, but require an
## explicit ID list; silent regex-based feature selection is not used.
orphan_index <- function(ids, orphan_id_file, expected_n = NULL) {
  if (!file.exists(orphan_id_file)) {
    stop("orphan ID file not found: ", orphan_id_file,
         "\nRun scripts/03_extract_orphan_ids.R first.")
  }
  orph <- unique(readLines(orphan_id_file, warn = FALSE))
  orph <- orph[nzchar(orph)]

  candidates <- unique(c(
    orph,
    ifelse(startsWith(orph, "ORPHAN_"), sub("^ORPHAN_", "", orph), paste0("ORPHAN_", orph))
  ))
  idx <- ids %in% candidates

  if (!is.null(expected_n) && sum(idx) != expected_n) {
    not_found <- setdiff(orph, c(ids, sub("^ORPHAN_", "", ids)))
    stop("orphan count is ", sum(idx), " but ", expected_n, " was expected.",
         if (length(not_found)) paste0(" Example unmatched ID: ", not_found[1]) else "")
  }
  idx
}

## --- within-sample transformations -----------------------------------------

transform_expression <- function(v, transform = c("rank", "scale")) {
  transform <- match.arg(transform)
  if (transform == "rank") {
    return(rank(v, ties.method = "average"))
  }
  s <- sd(v)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(v)))
  (v - mean(v)) / s
}

## --- one-patient paired gene-set test --------------------------------------

## zero_rule = "or" retains a feature when TPM is nonzero in either member of
## the matched pair. zero_rule = "none" keeps all features, including zeros.
## Selection occurs before transformation, reproducing the manuscript rank
## analysis in which the rank universe depends on the retained feature set.
paired_test <- function(tpm_n, tpm_d, target,
                        transform = c("rank", "scale"),
                        zero_rule = c("or", "none")) {
  transform <- match.arg(transform)
  zero_rule <- match.arg(zero_rule)
  stopifnot(length(tpm_n) == length(tpm_d), length(target) == length(tpm_n))

  keep <- if (zero_rule == "or") (tpm_n != 0 | tpm_d != 0) else rep(TRUE, length(tpm_n))
  if (sum(keep) < 3L) {
    return(data.frame(p_up = NA_real_, p_down = NA_real_, mean_d = NA_real_,
                      sd_d = NA_real_, n_tested = 0L, n_kept = sum(keep)))
  }

  xn <- transform_expression(tpm_n[keep], transform)
  xd <- transform_expression(tpm_d[keep], transform)
  sel <- target[keep]
  d <- xn[sel] - xd[sel]
  d <- d[is.finite(d)]

  n_tested <- length(d)
  if (n_tested < 3L) {
    return(data.frame(p_up = NA_real_, p_down = NA_real_, mean_d = NA_real_,
                      sd_d = NA_real_, n_tested = n_tested, n_kept = sum(keep)))
  }

  p_up   <- safe_t_test(d, alternative = "less")     # disease > normal
  p_down <- safe_t_test(d, alternative = "greater")  # disease < normal

  data.frame(
    p_up = p_up,
    p_down = p_down,
    mean_d = mean(d),
    sd_d = sd(d),
    n_tested = n_tested,
    n_kept = sum(keep)
  )
}

## --- cohort-level patient analysis -----------------------------------------

run_cohort <- function(tpm, sample_sheet, target,
                       transform = c("rank", "scale"),
                       zero_rule = c("or", "none"), alpha = 0.01) {
  transform <- match.arg(transform)
  zero_rule <- match.arg(zero_rule)
  patients <- unique(sample_sheet$patient_id)
  out <- vector("list", length(patients))

  for (k in seq_along(patients)) {
    p <- patients[k]
    jN <- which(sample_sheet$patient_id == p & sample_sheet$condition == "normal")
    jD <- which(sample_sheet$patient_id == p & sample_sheet$condition == "disease")
    if (length(jN) != 1L || length(jD) != 1L) {
      stop("patient ", p, " does not have exactly one normal and one disease library")
    }

    res <- paired_test(tpm[, jN], tpm[, jD], target,
                       transform = transform, zero_rule = zero_rule)
    out[[k]] <- cbind(patient_id = p, res, stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, out)
  res$q_up   <- p.adjust(res$p_up,   method = "BH")
  res$q_down <- p.adjust(res$p_down, method = "BH")
  res$direction <- ifelse(!is.na(res$q_up) & res$q_up < alpha, "up",
                   ifelse(!is.na(res$q_down) & res$q_down < alpha, "down", "none"))
  res$transform <- transform
  res$zero_rule <- zero_rule
  res
}

count_directions <- function(res) {
  data.frame(
    up   = sum(res$direction == "up",   na.rm = TRUE),
    down = sum(res$direction == "down", na.rm = TRUE),
    none = sum(res$direction == "none", na.rm = TRUE)
  )
}

## Convert counts to the column order used in manuscript Tables 2-4.
format_direction_table <- function(counts, table_name) {
  z <- counts[counts$table == table_name, , drop = FALSE]
  cohorts <- unique(z$cohort)
  out <- vector("list", length(cohorts))

  for (k in seq_along(cohorts)) {
    c0 <- cohorts[k]
    inc <- z[z$cohort == c0 & z$zero == "included", , drop = FALSE]
    exc <- z[z$cohort == c0 & z$zero == "excluded", , drop = FALSE]
    if (nrow(inc) != 1L || nrow(exc) != 1L) stop("incomplete results for ", table_name, "/", c0)

    n <- inc$up + inc$down + inc$none
    out[[k]] <- data.frame(
      GEO_ID = c0,
      disease_gt_normal_included_ge_alpha = n - inc$up,
      disease_gt_normal_included_lt_alpha = inc$up,
      disease_gt_normal_excluded_ge_alpha = n - exc$up,
      disease_gt_normal_excluded_lt_alpha = exc$up,
      disease_lt_normal_included_ge_alpha = n - inc$down,
      disease_lt_normal_included_lt_alpha = inc$down,
      disease_lt_normal_excluded_ge_alpha = n - exc$down,
      disease_lt_normal_excluded_lt_alpha = exc$down,
      no_direction_included = inc$none,
      no_direction_excluded = exc$none,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

## --- expression-matched non-orphan controls --------------------------------

## One prespecified low-abundance control: match the orphan set by the pair
## mean m=(normal+disease)/2 using a zero stratum plus n_strata quantile strata.
## This does NOT prove uniqueness of orphan genes; it tests whether low abundance
## alone is sufficient to explain the observed rank-based directional pattern.
draw_control_set <- function(tpm_n, tpm_d, is_orphan, n_strata = 20) {
  stopifnot(length(tpm_n) == length(tpm_d), length(is_orphan) == length(tpm_n))
  m <- (tpm_n + tpm_d) / 2
  stratum <- rep(1L, length(m))   # 1 = exact zero pair mean
  pos <- m > 0

  if (any(pos)) {
    br <- unique(quantile(m[pos], probs = (0:n_strata) / n_strata,
                          na.rm = TRUE, names = FALSE, type = 7))
    if (length(br) < 2L) {
      stratum[pos] <- 2L
    } else {
      stratum[pos] <- 1L + as.integer(cut(m[pos], breaks = br,
                                          include.lowest = TRUE, labels = FALSE))
    }
  }

  ctrl <- logical(length(m))
  orphan_strata <- sort(unique(stratum[is_orphan]))

  for (st in orphan_strata) {
    need <- sum(is_orphan & stratum == st)
    pool <- which(!is_orphan & stratum == st)
    if (length(pool) < need) {
      stop("insufficient non-orphan controls in stratum ", st,
           ": need ", need, ", available ", length(pool),
           ". Consider fewer strata, but record that choice explicitly.")
    }
    ctrl[sample(pool, need, replace = FALSE)] <- TRUE
  }

  if (sum(ctrl) != sum(is_orphan)) {
    stop("control set contains ", sum(ctrl), " transcripts; expected ", sum(is_orphan))
  }
  ctrl
}

run_control_sets <- function(tpm, sample_sheet, is_orphan,
                             n_rep = 20, n_strata = 20,
                             transform = "rank", zero_rule = "or",
                             alpha = 0.01, base_seed = 1L,
                             keep_patient_details = TRUE) {
  patients <- unique(sample_sheet$patient_id)
  counts <- matrix(NA_integer_, n_rep, 3,
                   dimnames = list(NULL, c("up", "down", "none")))
  n_ctrl <- matrix(NA_integer_, n_rep, length(patients),
                   dimnames = list(NULL, patients))
  detail <- if (keep_patient_details) vector("list", n_rep) else NULL

  for (rep in seq_len(n_rep)) {
    set.seed(base_seed + rep - 1L)
    cat("  control replicate ", rep, "/", n_rep, "\n", sep = "")
    per_patient <- vector("list", length(patients))

    for (k in seq_along(patients)) {
      p <- patients[k]
      jN <- which(sample_sheet$patient_id == p & sample_sheet$condition == "normal")
      jD <- which(sample_sheet$patient_id == p & sample_sheet$condition == "disease")

      ctrl <- draw_control_set(tpm[, jN], tpm[, jD], is_orphan, n_strata = n_strata)
      n_ctrl[rep, k] <- sum(ctrl)
      ans <- paired_test(tpm[, jN], tpm[, jD], ctrl,
                         transform = transform, zero_rule = zero_rule)
      per_patient[[k]] <- cbind(patient_id = p, ans, stringsAsFactors = FALSE)
    }

    res <- do.call(rbind, per_patient)
    res$q_up   <- p.adjust(res$p_up,   "BH")
    res$q_down <- p.adjust(res$p_down, "BH")
    res$direction <- ifelse(!is.na(res$q_up) & res$q_up < alpha, "up",
                     ifelse(!is.na(res$q_down) & res$q_down < alpha, "down", "none"))
    res$replicate <- rep

    counts[rep, ] <- c(
      sum(res$direction == "up"),
      sum(res$direction == "down"),
      sum(res$direction == "none")
    )
    if (keep_patient_details) detail[[rep]] <- res
  }

  list(
    counts = counts,
    n_control_transcripts = n_ctrl,
    details = if (keep_patient_details) do.call(rbind, detail) else NULL
  )
}

summarise_control <- function(counts) {
  fmt <- function(v) sprintf("%s (%s-%s)", format(median(v)), min(v), max(v))
  data.frame(
    up = fmt(counts[, "up"]),
    down = fmt(counts[, "down"]),
    none = fmt(counts[, "none"]),
    stringsAsFactors = FALSE
  )
}

## --- positive control: low relative expression -----------------------------

orphan_lower_expression <- function(tpm, is_orphan,
                                    zero = c("included", "excluded"),
                                    alpha = 0.05) {
  zero <- match.arg(zero)
  n <- ncol(tpm)
  p <- rep(NA_real_, n)

  for (j in seq_len(n)) {
    keep <- if (zero == "included") rep(TRUE, nrow(tpm)) else tpm[, j] != 0
    if (sum(keep) < 3L) next
    x <- transform_expression(tpm[keep, j], "scale")
    sel <- is_orphan[keep]
    if (sum(sel) < 2L || sum(!sel) < 2L) next
    p[j] <- t.test(x[sel], x[!sel], alternative = "less")$p.value
  }

  q <- p.adjust(p, "BH")
  list(
    p = p,
    q = q,
    counts = data.frame(
      ge = sum(is.na(q) | q >= alpha),
      lt = sum(!is.na(q) & q < alpha)
    )
  )
}
