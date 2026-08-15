## ---------------------------------------------------------------------------
## gtf_attr.R  --  shared GTF attribute parser
## Sourced by scripts/03_extract_orphan_ids.R and scripts/04_aggregate_gene_level.R
## ---------------------------------------------------------------------------

## Pull one attribute value out of GTF column 9.
## The key is anchored to the start of the field or to a preceding "; ", so that
## a key is not matched inside a longer attribute name (e.g. "gene_id" must not
## match "ref_gene_id").
gtf_attribute <- function(attr_col, key) {
  pat <- paste0('(^|; *)', key, ' +"[^"]*"')
  m <- regmatches(attr_col, regexpr(pat, attr_col))
  out <- rep(NA_character_, length(attr_col))
  ok <- nzchar(m)
  v <- sub(paste0('^(; *)?', key, ' +"'), "", m[ok])
  v <- sub('"$', "", v)
  out[ok] <- v
  out
}
