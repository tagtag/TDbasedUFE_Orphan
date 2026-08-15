## Convenience wrapper for one abundance file.
##
## Examples:
##   Rscript scripts/04a_aggregate_one_abundance.R \
##     data/reference/combined_GRCh37.gtf \
##     sample/abundance.csv sample/abundance_gene.csv
##
##   Rscript scripts/04a_aggregate_one_abundance.R \
##     data/reference/combined_GRCh37.gtf \
##     sample/abundance.tsv sample/abundance_gene.tsv
##
## This wrapper delegates to 04_aggregate_gene_level.R, whose single-file mode
## performs the actual transcript -> gene-locus aggregation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(paste(
    "usage: Rscript scripts/04a_aggregate_one_abundance.R",
    "COMBINED_GTF INPUT_ABUNDANCE [OUTPUT_GENE_ABUNDANCE]"
  ))
}
cmd <- c("scripts/04_aggregate_gene_level.R", args)
status <- system2("Rscript", cmd)
if (status != 0L) quit(status = status)
