## Run the reproducible analysis from transcript-level abundance files onward.
## Gene-locus abundance files are generated first from abundance.tsv/csv.
## Metadata/run selection is reconstructed from public GEO/ENA plus the
## published GSE127165 supplementary table when required.  Set
## REBUILD_GEO_METADATA=1 to force a fresh reconstruction.

source("R/config.R")

rscript <- file.path(R.home("bin"), "Rscript")
kroot <- Sys.getenv("KALLISTO_ROOT", unset = "data/kallisto")
gtf <- Sys.getenv("COMBINED_GTF", unset = "data/reference/combined_GRCh37.gtf")
rebuild_meta <- identical(Sys.getenv("REBUILD_GEO_METADATA", unset = "0"), "1")

required_meta <- c(
  GEO_PATIENT_METADATA_FILE,
  SRA_RUNS_FILE,
  vapply(COHORTS, function(x) x$sheet, character(1))
)
file_ready <- vapply(required_meta, function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(FALSE)
  x <- tryCatch(readLines(f, n = 20L, warn = FALSE), error = function(e) character())
  if (any(grepl("xxxx|placeholder", x, ignore.case = TRUE))) return(FALSE)
  if (grepl("^(sra_runs|sample_sheet_GSE)", basename(f)) && length(x) <= 1L) {
    return(FALSE)
  }
  TRUE
}, logical(1))

if (rebuild_meta || !all(file_ready)) {
  cat("\n=============== GEO/SRA/ENA metadata: build ===============\n")
  status <- system2(rscript, "scripts/00_build_geo_patient_metadata.R")
  if (!identical(status, 0L)) stop("GEO/SRA/ENA metadata construction failed")
} else {
  cat("\n=============== GEO/SRA/ENA metadata: reuse ===============\n")
  cat("Using existing validated combined metadata, sra_runs.csv, and sample sheets.\n")
  cat("Set REBUILD_GEO_METADATA=1 to reconstruct all metadata from public GEO/ENA sources.\n")
}

cat("\n=============== transcript -> gene loci ===============\n")
status <- system2(rscript, c("scripts/04_aggregate_gene_level.R", gtf, kroot))
if (!identical(status, 0L)) stop("gene-locus aggregation failed")

steps <- c(
  "R/01_table1_relative_expression.R",
  "R/02_patient_level_analysis.R",
  "R/03_table5_control_sets.R",
  "R/04_clinical_association.R",
  "R/05_figure2.R",
  "R/06_session_info.R",
  "R/07_validate_manuscript_counts.R"
)

for (f in steps) {
  cat("\n=============== ", f, " ===============\n", sep = "")
  status <- system2(rscript, f)
  if (!identical(status, 0L)) stop("analysis step failed: ", f)
}
