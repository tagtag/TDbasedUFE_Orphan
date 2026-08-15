## ---------------------------------------------------------------------------
## 04_aggregate_gene_level.R
##
## Aggregate kallisto transcript-level abundance to annotated gene loci.
## This script supports both the standard kallisto TSV and CSV copies:
##
##   abundance.tsv  -> abundance_gene.tsv
##   abundance.csv  -> abundance_gene.csv
##
## Arbitrary filenames are also accepted, e.g.:
##   Rscript scripts/04_aggregate_gene_level.R \
##     data/reference/combined_GRCh37.gtf sample/abundan.csv sample/abndan_gene.csv
##
## Orphan locus IDs retain an explicit ORPHAN_ prefix, e.g.
## ORPHAN_XLOC_022250, so orphan status remains visible in the identifier.
##
## Modes
## -----
## 1) One file:
##    Rscript scripts/04_aggregate_gene_level.R GTF INPUT [OUTPUT]
##
## 2) Every immediate subdirectory under a kallisto root:
##    Rscript scripts/04_aggregate_gene_level.R GTF KALLISTO_ROOT [overwrite]
##
## In batch mode the script auto-detects abundance.tsv or abundance.csv in each
## library directory and writes the same format with _gene inserted before the
## extension.
## ---------------------------------------------------------------------------

source("scripts/gtf_attr.R")

read_abundance_table <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    z <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    z <- read.delim(file, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "")
  }
  required <- c("target_id", "tpm")
  missing <- setdiff(required, colnames(z))
  if (length(missing)) {
    stop("missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(z$target_id)) stop("duplicated target_id values in ", file)
  if (any(!is.finite(z$tpm)) || any(z$tpm < 0)) stop("invalid TPM values in ", file)
  z
}

write_abundance_table <- function(z, file) {
  ext <- tolower(tools::file_ext(file))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  if (ext == "csv") {
    write.csv(z, file, row.names = FALSE, quote = FALSE)
  } else {
    write.table(z, file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

derive_gene_filename <- function(file) {
  ext <- tools::file_ext(file)
  if (!nzchar(ext)) return(paste0(file, "_gene"))
  stem <- sub(paste0("\\.", ext, "$"), "", file)
  paste0(stem, "_gene.", ext)
}

read_tx2gene <- function(gtf_file) {
  if (!file.exists(gtf_file)) stop("combined GTF not found: ", gtf_file)
  cat("reading ", gtf_file, "\n", sep = "")
  gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#",
                    stringsAsFactors = FALSE, quote = "", fill = TRUE)
  if (ncol(gtf) < 9L) stop("GTF has fewer than 9 columns: ", gtf_file)

  attr <- gtf[[9]]
  map <- unique(data.frame(
    transcript_id = gtf_attribute(attr, "transcript_id"),
    gene_id       = gtf_attribute(attr, "gene_id"),
    stringsAsFactors = FALSE
  ))
  map <- map[!is.na(map$transcript_id) & !is.na(map$gene_id), , drop = FALSE]

  ## A transcript may occur on many exon rows, but after unique() it must map
  ## to exactly one gene_id.
  n_gene <- tapply(map$gene_id, map$transcript_id, function(x) length(unique(x)))
  bad <- names(n_gene)[n_gene > 1L]
  if (length(bad)) stop("transcript maps to multiple gene IDs: ", bad[1])

  cat("transcripts in map: ", length(unique(map$transcript_id)), "\n", sep = "")
  cat("gene loci in map:   ", length(unique(map$gene_id)),
      " (manuscript: 58962)\n", sep = "")
  map
}

map_targets_to_loci <- function(target_id, map) {
  lookup <- setNames(map$gene_id, map$transcript_id)
  gene <- unname(lookup[target_id])

  ## Compatibility with a raw Ruiz-Orera GTF lacking ORPHAN_ prefixes, or with
  ## abundance files produced from a prefixed reference.
  miss <- is.na(gene)
  if (any(miss)) {
    q <- target_id[miss]
    alt <- ifelse(startsWith(q, "ORPHAN_"),
                  sub("^ORPHAN_", "", q), paste0("ORPHAN_", q))
    gene[miss] <- unname(lookup[alt])
  }

  if (anyNA(gene)) {
    i <- which(is.na(gene))[1]
    stop("transcript is absent from the combined GTF: ", target_id[i])
  }

  ## Even if the GTF itself was not prefixed, keep orphan status auditable in
  ## the gene-level output identifier.
  need <- startsWith(target_id, "ORPHAN_") & !startsWith(gene, "ORPHAN_")
  gene[need] <- paste0("ORPHAN_", gene[need])
  gene
}

aggregate_one_abundance <- function(input_file, output_file, map) {
  z <- read_abundance_table(input_file)
  gene <- map_targets_to_loci(as.character(z$target_id), map)
  levels_this <- sort(unique(gene))

  ## TPM is additive across transcript isoforms belonging to the same locus.
  tpm <- tapply(z$tpm, factor(gene, levels = levels_this), sum)
  out <- data.frame(
    target_id = levels_this,
    tpm = as.numeric(tpm),
    stringsAsFactors = FALSE
  )

  ## est_counts is also additive, and is retained when present for auditability.
  if ("est_counts" %in% colnames(z)) {
    ec <- tapply(z$est_counts, factor(gene, levels = levels_this), sum)
    out$est_counts <- as.numeric(ec)
  }

  ## Record how many transcript rows contributed to each gene locus.  The
  ## downstream manuscript analysis uses target_id and tpm only.
  nt <- table(factor(gene, levels = levels_this))
  out$n_transcripts <- as.integer(nt)

  ## TPM total should be preserved apart from floating-point rounding.
  if (!isTRUE(all.equal(sum(z$tpm), sum(out$tpm), tolerance = 1e-8))) {
    stop("TPM total was not preserved while aggregating ", input_file)
  }

  write_abundance_table(out, output_file)
  cat("[ok] ", input_file, " -> ", output_file,
      "  (", nrow(z), " transcripts -> ", nrow(out), " loci; ",
      sum(startsWith(out$target_id, "ORPHAN_")), " orphan loci)\n", sep = "")
  invisible(out)
}

find_transcript_abundance <- function(d) {
  preferred <- c("abundance.tsv", "abundance.csv", "abundan.tsv", "abundan.csv")
  hit <- file.path(d, preferred)
  hit <- hit[file.exists(hit)]
  if (!length(hit)) return(NA_character_)
  if (length(hit) > 1L) {
    ## Prefer the canonical kallisto TSV if both a TSV and a CSV copy exist.
    canonical <- file.path(d, "abundance.tsv")
    if (canonical %in% hit) return(canonical)
    stop("multiple transcript abundance files found in ", d, ": ",
         paste(basename(hit), collapse = ", "))
  }
  hit[1]
}

args <- commandArgs(trailingOnly = TRUE)
gtf_file <- if (length(args) >= 1L) args[1] else "data/reference/combined_GRCh37.gtf"
second   <- if (length(args) >= 2L) args[2] else "data/kallisto"
third    <- if (length(args) >= 3L) args[3] else NA_character_

map <- read_tx2gene(gtf_file)

if (file.exists(second) && !dir.exists(second)) {
  ## Single-file mode.
  input_file <- second
  output_file <- if (!is.na(third) && nzchar(third)) third else derive_gene_filename(input_file)
  aggregate_one_abundance(input_file, output_file, map)
} else {
  ## Batch mode.
  kal_dir <- second
  if (!dir.exists(kal_dir)) stop("kallisto directory not found: ", kal_dir)
  overwrite <- if (!is.na(third)) as.logical(third) else FALSE
  if (is.na(overwrite)) stop("third argument in batch mode must be TRUE or FALSE")

  dirs <- list.dirs(kal_dir, recursive = FALSE, full.names = TRUE)
  cat("library directories found: ", length(dirs), "\n", sep = "")
  n_done <- 0L
  for (d in dirs) {
    fin <- find_transcript_abundance(d)
    if (is.na(fin)) next
    fout <- derive_gene_filename(fin)
    if (file.exists(fout) && !overwrite) {
      cat("[skip, exists] ", fout, "\n", sep = "")
      next
    }
    aggregate_one_abundance(fin, fout, map)
    n_done <- n_done + 1L
  }
  cat("gene-level abundance files created: ", n_done, "\n", sep = "")
}
