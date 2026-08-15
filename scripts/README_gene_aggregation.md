# Transcript abundance to gene-locus abundance

This repository includes the conversion used for the manuscript gene-locus
sensitivity analysis. The transcript-to-gene mapping is read from the same
combined GTF used to build the kallisto index.

For CSV input:

```bash
Rscript scripts/04_aggregate_gene_level.R \
  data/reference/combined_GRCh37.gtf \
  sample/abundance.csv \
  sample/abundance_gene.csv
```

For kallisto TSV input:

```bash
Rscript scripts/04_aggregate_gene_level.R \
  data/reference/combined_GRCh37.gtf \
  sample/abundance.tsv \
  sample/abundance_gene.tsv
```

The input must contain `target_id` and `tpm`. TPM values of all transcripts
sharing a `gene_id` are summed. If `est_counts` is present, it is summed too.
The output also records `n_transcripts` for auditability. Orphan loci keep an
`ORPHAN_` prefix.

The script also accepts arbitrary filenames, so `abundan.csv` can be written
as `abndan_gene.csv` simply by giving those names as the second and third
arguments.
