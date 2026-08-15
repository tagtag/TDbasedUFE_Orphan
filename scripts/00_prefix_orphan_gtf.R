## ---------------------------------------------------------------------------
## 00_prefix_orphan_gtf.R
## Add an explicit ORPHAN_ prefix to gene_id and transcript_id in the raw
## Ruiz-Orera human de novo GTF. This reproduces identifiers such as
## ORPHAN_hsa_00255579 and ORPHAN_XLOC_022250 in kallisto output and makes
## orphan features auditable by identifier as well as by an explicit ID list.
##
## Usage:
##   Rscript scripts/00_prefix_orphan_gtf.R input.gtf output.gtf
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript scripts/00_prefix_orphan_gtf.R input.gtf output.gtf")
}
infile <- args[1]
outfile <- args[2]
if (!file.exists(infile)) stop("input GTF not found: ", infile)

a <- readLines(infile, warn = FALSE)

prefix_attr <- function(x, key) {
  pat <- paste0('(', key, ' +")([^"]+)(")')
  m <- regexec(pat, x)
  hit <- regmatches(x, m)
  out <- x
  ok <- lengths(hit) == 4L
  for (i in which(ok)) {
    val <- hit[[i]][3]
    if (!startsWith(val, "ORPHAN_")) {
      repl <- paste0(hit[[i]][2], "ORPHAN_", val, hit[[i]][4])
      out[i] <- sub(pat, repl, out[i])
    }
  }
  out
}

idx <- !startsWith(a, "#") & nzchar(a)
a[idx] <- prefix_attr(a[idx], "gene_id")
a[idx] <- prefix_attr(a[idx], "transcript_id")
writeLines(a, outfile)
cat("prefixed orphan GTF written to ", outfile, "\n", sep = "")
