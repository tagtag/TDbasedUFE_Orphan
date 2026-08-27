# Disease-context-dependent directional dysregulation in human orphan genes
Y.-H. Taguchi  *,Turki Turki

https://doi.org/10.20944/preprints202608.1930.v1

Reproducible code for the patient-level analysis of human orphan-gene expression
in four paired RNA-seq cohorts.

The public-data pipeline is deliberately **manifest-free at the starting point**:
it does not require a hand-edited SRA run list, disease/normal table, patient-pair
table, or clinical override file.  `scripts/00_build_geo_patient_metadata.R`
queries the public GEO/ENA records, applies the cohort-specific rules documented
below, obtains the exact FASTQ URLs/MD5 checksums, constructs the matched-patient
sample sheets, and builds `GEO_patient_metadata_combined.tsv`.

| GEO | Disease | Matched pairs used in patient analysis |
|---|---|---:|
| GSE244679 | Psoriasis | 24 |
| GSE127165 | Laryngeal squamous cell carcinoma | 57 |
| GSE144269 | Hepatocellular carcinoma | 70 |
| GSE40419 | Lung adenocarcinoma | 69 |

The combined GRCh37 reference contains 198,507 transcripts, including 2,190
human de novo/orphan transcripts.  Transcript-level abundance is also collapsed
to 58,962 gene loci, including 1,226 orphan loci.

## What is reproduced

| Manuscript output | Script |
|---|---|
| Public GEO/ENA/SRA/FASTQ manifest and pairing | `scripts/00_build_geo_patient_metadata.R` |
| Kallisto quantification from exact ENA FASTQ URLs | `scripts/02_run_kallisto.sh` |
| Transcript -> gene-locus abundance | `scripts/04_aggregate_gene_level.R` |
| Table 1 | `R/01_table1_relative_expression.R` |
| Tables 2-4 | `R/02_patient_level_analysis.R` |
| Table 5, expression-matched controls | `R/03_table5_control_sets.R` |
| Clinical analyses / Table 7 | `R/04_clinical_association.R` |
| Figure 2 | `R/05_figure2.R` |
| Manuscript-count checks | `R/07_validate_manuscript_counts.R` |

## Repository layout

```text
scripts/
  00_build_geo_patient_metadata.R   # GEO -> ENA runs/FASTQ + pairing + clinical TSV
  make_orphan_GEO_metadata.R        # backward-compatible wrapper
  00_prefix_orphan_gtf.R
  01_build_reference.sh
  02_run_kallisto.sh                # uses exact URLs in metadata/sra_runs.csv
  03_extract_orphan_ids.R
  04_aggregate_gene_level.R
  04a_aggregate_one_abundance.R
R/
  config.R
  00_functions.R
  01_table1_relative_expression.R
  02_patient_level_analysis.R
  03_table5_control_sets.R
  04_clinical_association.R
  05_figure2.R
  06_session_info.R
  07_validate_manuscript_counts.R
  run_analysis.R
metadata/
  sra_runs.csv                       # generated; exact FASTQ URLs and MD5
  sample_sheet_GSE*.csv              # generated matched pairs
  orphan_delta_input.tsv             # optional compatibility input
metadata_output/                     # generated audit tables
results/
figures/
tests/
  smoke_test.R
run_all.sh
```

## Requirements

- R >= 4.0
- Bioconductor: `GEOquery`, `Biobase`
- CRAN: `jsonlite`, `readxl`, `xml2`, `BayesFactor`
- `kallisto`, `gffread`, `samtools`, `wget`

The metadata builder itself needs `GEOquery`, `Biobase`, `jsonlite`, `readxl`
and `xml2`.  `BayesFactor` is used by the PASI analysis.


## 1. Public metadata, exact FASTQs, disease/normal labels and patient pairing

Run:

```bash
Rscript scripts/00_build_geo_patient_metadata.R
```

The script reconstructs the dataset from public records, without a local
selection file.

### GSE244679

GEO contains 24 matched lesional/adjacent-normal pairs.  The sample title holds
the patient number and condition, e.g. `[71N]` / `[71D]`:

```text
N -> normal (adjacent normal skin)
D -> disease (psoriatic lesion)
```

The numeric value inside the brackets is used as the patient ID.  Each GSM's
GEO SRA relation is resolved through ENA to its exact `run_accession`,
`library_layout`, FASTQ URL(s), and MD5 checksum(s).

### GSE127165

GEO contains 57 LSCC and 57 paired adjacent normal mucosa samples.  Titles are
paired by subject number:

```text
LSCC_<subject> -> disease
ANM_<subject>  -> normal
```

GEO supplies subject number, tumour stage, tobacco smoking and alcohol
consumption.  Patient age, which is not exposed in the GEO sample
characteristics, is retrieved automatically from the primary paper's public
**Additional file 1, Table S1 (Clinical features of 57 LSCC samples for RNA
sequencing)**.  The script first uses the Figshare article API (article
`12414305`) and has the publisher-hosted ZIP as a fallback.  It identifies the
RNA-sequenced Table-S1 records by matching the published subject numbers to the
57 GEO subject IDs; in strict mode all 57 ages must be recovered or the script
stops rather than guessing.

### GSE144269

GEO itself defines 140 RNA-seq samples from 70 matched HCC tumour/non-tumour
pairs.  Titles are parsed as:

```text
pat <n> tumor     -> disease
pat <n> non-tumor -> normal
```

The patient phenotype deposited under controlled access is not required for the
manuscript's clinical analysis of this cohort.

### GSE40419

This cohort needs a run-level reconstruction rather than simply treating every
GSM as one downloaded library.  GEO states that the raw study is EBI-SRA
`ERP001058` / BioProject `PRJEB3132` and lists tumour titles such as `LC_S9`
and matched adjacent-normal titles such as `LC_S9_nor`.

The code therefore:

1. downloads the complete ENA `read_run` report for `ERP001058`;
2. retains only run rows for which ENA supplies a complete direct FASTQ set
   consistent with `library_layout` (two FASTQs for paired-end, one for
   single-end);
3. requires this public rule to reproduce the manuscript's **162 screened
   libraries**;
4. matches ENA `sample_alias` exactly to the GEO sample title;
5. classifies a matched title ending in `_nor` as normal and its base title as
   the patient ID; a matched title without `_nor` is disease;
6. requires the screened set to reproduce 69 disease, 76 normal and 17
   unclassified libraries;
7. retains the intersection of disease and normal patient IDs, which must give
   69 complete pairs; the other 24 screened libraries are written to an
   exclusion audit table.

No `KALLISTO_ROOT` lookup and no run-selection override are used to decide these
runs.  If the current public archive does not reproduce the manuscript counts,
the builder stops and exposes the mismatch instead of silently selecting a
replacement run.

### Generated metadata files

The builder writes:

```text
metadata/sra_runs.csv
metadata/sample_sheet_GSE244679.csv
metadata/sample_sheet_GSE127165.csv
metadata/sample_sheet_GSE144269.csv
metadata/sample_sheet_GSE40419.csv
metadata_output/GEO_patient_metadata_combined.tsv
metadata_output/GEO_sample_metadata_long.tsv
metadata_output/GEO_SRA_run_mapping.tsv
metadata_output/FASTQ_selection_audit.tsv
metadata_output/ERP001058_read_run.tsv
metadata_output/GSE40419_excluded_libraries.tsv
metadata_output/GSE127165_supplement_clinical.tsv
metadata_output/public_source_provenance.tsv
```

`sra_runs.csv` contains the exact files selected for download:

```text
run_accession,cohort,gsm,patient_id,condition,library_layout,
fastq_1,fastq_2,md5_1,md5_2
```

The unclassified/unmatched GSE40419 rows are kept in `sra_runs.csv` because the
library-level Table 1 uses all 162 screened libraries.  They are absent from the
matched-pair sample sheet used by Tables 2-5.

## 2. Build the combined orphan + Ensembl GRCh37 reference

Place the Ruiz-Orera human de novo GTF and the Ensembl GRCh37.75 GTF in
`data/reference/`, then run:

```bash
bash scripts/01_build_reference.sh
Rscript scripts/03_extract_orphan_ids.R \
  data/reference/hsa_denovo.ORPHAN_prefixed.gtf
```

Orphan transcript/locus identifiers retain the `ORPHAN_` prefix, e.g.
`ORPHAN_hsa_00255579` and `ORPHAN_XLOC_022250`.

## 2b. Build kallisto index

Prior to the following process, kallisto index must be build. See reference_build sub directory

## 3. Download the exact FASTQs and run kallisto

After Step 1, do:

```bash
bash scripts/02_run_kallisto.sh 8
```

The script does **not** infer FTP paths from SRA run names.  It downloads the
exact `fastq_1` / `fastq_2` URLs selected by the R metadata builder and verifies
ENA MD5 values when available.  Outputs are written as:

```text
data/kallisto/<GSE>_<SRR-or-ERR>/abundance.tsv
```

By default each FASTQ is removed after successful quantification.  To keep raw
FASTQs:

```bash
KEEP_FASTQ=1 bash scripts/02_run_kallisto.sh 8
```

## 4. Create `abundance_gene` from `abundance`

The gene-locus conversion is included.  For one file:

```bash
Rscript scripts/04_aggregate_gene_level.R \
  data/reference/combined_GRCh37.gtf \
  path/to/abundance.csv \
  path/to/abundance_gene.csv
```

The same command accepts `abundance.tsv`.  TPM values of transcripts sharing
the same GTF `gene_id` are summed.  If `est_counts` is present it is also
summed.  The batch form automatically visits all kallisto directories:

```bash
Rscript scripts/04_aggregate_gene_level.R \
  data/reference/combined_GRCh37.gtf data/kallisto
```

## 5. Run Tables 1-5, clinical analyses and Figure 2

If the abundance files already exist:

```bash
bash run_all.sh
```

To force public metadata reconstruction first:

```bash
REBUILD_GEO_METADATA=1 bash run_all.sh
```

Or use:

```bash
Rscript R/run_analysis.R
```

If abundance directories are elsewhere:

```bash
KALLISTO_ROOT=/path/to/kallisto_outputs bash run_all.sh
```

## Statistical implementation notes

**Sign convention.** Patient-level difference is `normal - disease`; negative
means higher expression in disease.

**Zero-abundance analysis.** The primary/sensitivity analyses can retain all
features or retain features nonzero in either member of the matched pair;
selection occurs before standardization/ranking.

**Rank analysis.** Within-sample ranks use `rank(..., ties.method="average")`.

**BH adjustment.** The two one-sided patient-level tests are adjusted separately
across patients within each cohort; manuscript threshold is adjusted `P < 0.01`.

**Expression-matched controls.** Controls are constructed independently for each
patient by the paired mean TPM `(normal + disease)/2`, using an exact-zero
stratum and 20 positive-expression quantile strata.  Each set contains 2,190
non-orphan transcripts and the procedure is repeated 20 times with fixed seeds.
This tests whether low abundance alone explains the rank-based directional
pattern; it is not claimed to prove uniqueness under all possible matching
schemes.

**Clinical metadata.** `R/04_clinical_association.R` uses only
`metadata_output/GEO_patient_metadata_combined.tsv` and matches analysed
patients through `normal_run_accession`, not row order.

## Reproducibility safeguards

By default the metadata builder asserts the manuscript reconstruction:

```text
screened libraries: 48, 114, 140, 162
matched pairs:       24, 57, 70, 69
GSE40419 screened classification: 69 disease, 76 normal, 17 unclassified
```

If public archive metadata changes, the script stops and leaves audit files.
For exploratory inspection only, strict assertions can be disabled with
`STRICT_MANUSCRIPT_COUNTS=0`; this should not be used for manuscript
reproduction.

Finally run:

```bash
Rscript R/07_validate_manuscript_counts.R
```

to check the generated Tables 2-5 against the manuscript counts.
