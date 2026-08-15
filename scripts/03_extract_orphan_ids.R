## Extract explicit orphan transcript and gene-locus ID lists from the GTF.
source("scripts/gtf_attr.R")

args <- commandArgs(trailingOnly = TRUE)
gtf_file <- if (length(args) >= 1L) args[1] else "data/reference/hsa_denovo.ORPHAN_prefixed.gtf"
out_dir <- if (length(args) >= 2L) args[2] else "metadata"
if (!file.exists(gtf_file)) stop("orphan GTF not found: ", gtf_file)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, quote = "", fill = TRUE)
if (ncol(gtf) < 9L) stop("not a valid GTF: fewer than 9 columns")
attr <- gtf[[9]]

tx <- unique(gtf_attribute(attr, "transcript_id"))
gid <- unique(gtf_attribute(attr, "gene_id"))
tx <- sort(tx[!is.na(tx) & nzchar(tx)])
gid <- sort(gid[!is.na(gid) & nzchar(gid)])

## If a raw, unprefixed GTF was supplied, normalize the exported IDs to the
## ORPHAN_ convention used by the reference-building script.
tx <- ifelse(startsWith(tx, "ORPHAN_"), tx, paste0("ORPHAN_", tx))
gid <- ifelse(startsWith(gid, "ORPHAN_"), gid, paste0("ORPHAN_", gid))

cat("orphan transcript ids: ", length(tx), " (manuscript: 2190)\n", sep = "")
cat("orphan gene ids:       ", length(gid), " (manuscript: 1226)\n", sep = "")

writeLines(tx, file.path(out_dir, "orphan_transcript_ids.txt"))
writeLines(gid, file.path(out_dir, "orphan_gene_ids.txt"))

if (length(tx) != 2190L) warning("orphan transcript count differs from 2190")
if (length(gid) != 1226L) warning("orphan gene-locus count differs from 1226")
