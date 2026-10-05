# TDbasedUFE_Orphan

Analysis code for *Reduced Detection of Human Orphan Gene Transcripts in Disease Tissue Relative to Expression- and Detection-Matched Controls*, https://doi.org/10.20944/preprints202608.1930.v2
(revision of Preprints 202608.1930.v1, https://doi.org/10.20944/preprints202608.1930.v1).

Everything from the transcript-level `abundance.tsv` files onward is here, together
with the scripts that build the reference and reconstruct the run selection from
public GEO/ENA metadata. Raw FASTQ and intermediate files are not tracked
(see `.gitignore`).

`analysis_manifest.tsv` maps every table and figure of the manuscript to the script
that produces it.

---

## What changed in the revision

The first version tested each library, then each patient, against a nominal
threshold. Reviewers asked for an explicit interaction term, for a p value that
respects transcript correlation, and for the rank analysis to be primary. Carrying
those out changed the conclusion: **the patient-level claim does not survive its own
negative control and has been withdrawn**, and the claim is now a cohort-level one.

The scripts that produced the first version's Tables 1–5 and 7 and Figure 2 are kept
in `legacy_v1/`. They are superseded, not deleted: the preprint they support is
public and citable.

---

## Layout

```
R/config.R                       constants (seeds, strata, expected counts)
R/00_functions.R                 shared I/O and transforms
R/07_interaction_permutation.R   the engine: beta3, resampling reference, external2
R/06_session_info.R              writes results/sessionInfo.txt

diag_slope.R                     Table 2   slope, intercept, amplification factor
de_framework.R                   Table 5   voom / edgeR / DESeq2 + camera
rho_compare.R                    Table 5   per-patient log2-ratio test, rho-bar
diag_batch.R                     Table 1   batch-confounding diagnostics
batch_extra.R                    Table 4   same-condition null
diag_mismatch.R                  Table 4   orph_excess with the pairing broken
diag_mismatch_beta3.R            Table 3   beta3 with the pairing broken
diag_mismatch_slope.R            Methods   the intercept with the pairing broken
mismatch_summary.R               Methods   pools the above across cohorts
check_draw_noise.R               Methods   resampling noise, from existing CSVs
clinical_assoc.R                 Table 7   PASI, tumour stage, alcohol
clinical_sensitivity.R           Table 7   influence points, extra covariates
R/07_matching_validation.R       Table S1  matching variables on synthetic data
test_2d.R                        Methods   two-dimensional matching on synthetic data
make_fig1.R                      Figure 1
make_fig2.R                      Figure 2

run_verify.sh                    batch: Tables 2 and 5
run_mismatch.sh                  batch: the intercept negative control
run_batch_audit.sh               batch: Tables 1 and 4

scripts/                         metadata reconstruction, kallisto, gene aggregation
reference_build/                 combined GRCh37 reference and kallisto index
metadata/                        sample sheets and run accessions
legacy_v1/                       superseded scripts of the first version
```

Every analysis script expects the **repository root as the working directory**; they
`source("R/config.R")` and read `R/07_interaction_permutation.R` by that path.

---

## Reference

| | |
|---|---|
| Genome / annotation | Ensembl GRCh37 release 75 (`Homo_sapiens.GRCh37.75`) |
| Orphan catalogue | hominoid orphan GTF, figshare article 1604892 |
| Combined reference | 2,190 orphan + 196,317 other = **198,507 transcripts** |
| Transcript FASTA | `gffread` from the combined GTF |
| Index | `kallisto index` with default k-mer length (31) |

`reference_build/build_kallisto_reference.sh` downloads the sources, records
SHA-256 checksums, asserts the three counts above, and writes
`reference_build_summary.txt` including the `kallisto version` string.

## Quantification

```
kallisto quant -i combined.idx -o <out> -t <threads> <R1> <R2>          # paired
kallisto quant -i combined.idx -o <out> -t <threads> --single -l 200 -s 30 <R1>
```

No bootstrap replicates. Single-end fragment length and s.d. are the script
defaults (`SE_FRAGMENT_LENGTH=200`, `SE_FRAGMENT_SD=30`) and can be overridden by
environment variable.

---

## Analysis settings used for the published numbers

| switch | value in the manuscript | affects | note |
|---|---|---|---|
| `MATCH_ON` | `external2` | Tables 2, 3, 4 | other patients' normal mean (20 strata) x cross-patient detection frequency (4 strata), leave-one-out |
| `TRANSFORMS` | `rank` primary, `scale` sensitivity | **Table 3 only** | zero-rate quantities pass through neither |
| `N_PERM` | 1000 | Table 3 | empirical P cannot fall below 1/1001 |
| `N_STRATA` / `N_STRATA2` | 20 / 4 | Tables 2, 3, 4, 5 | exact zeros form a stratum of their own in each dimension, so up to 21 x 5 cells |
| `NSET` (`diag_slope.R`) | **200** | Table 2 | the script default is 20; the manuscript uses 200 |
| `NSET` (`de_framework.R`, `rho_compare.R`) | 20 | Table 5 | |
| `MIN_COUNT` | 10 (edgeR default); 5 and 1 as sensitivity | Table 5, Supplementary | `filterByExpr` |
| `ENGINE` | `voom` primary; `deseq2` on GSE244679 only | Table 5, Supplementary | DESeq2 is slow with many blocking coefficients |
| `PRIOR` (`rho_compare.R`) | 1 TPM | Table 5 | the offset in log2((TPM+1)/(TPM+1)) |
| `NMIS` | 20 (`diag_mismatch*.R`), 5 (`diag_mismatch_beta3.R`) | Tables 3, 4 | mismatched partners per patient |
| `ZERO_RULE` | `or` (default) | — | `and` was considered and dropped |
| `NORMALIZE` | `tpm` (default) | — | `tmm` was considered and dropped; TMM already runs inside `de_framework.R` |
| `MIN_EXT` | 0 (disabled) | — | considered and dropped |
| `BASE_SEED` | 1 | all | per-patient seed is `BASE_SEED + i` |

Results are **bit-identical across `N_CORES`**, because the seed is set per patient
rather than per worker (checked with `N_CORES=1` against `N_CORES=4`).

Output file names carry the settings (`_<zero_rule>`, `_<normalize>`, `_ext<value>`,
`_<match_on>`, `_mc<min_count>`) so that a sensitivity run never overwrites the
primary one.

---

## Reproducing the tables

```bash
# Table 1 and Table 4
bash run_batch_audit.sh

# Table 2 and Table 5
bash run_verify.sh

# Table 3
for CO in GSE244679 GSE127165 GSE144269; do
  N_CORES=12 N_PERM=200 TRANSFORMS=rank NMIS=5 \
  KALLISTO_ROOT=Revised/$CO METADATA_DIR=Revised/metadata COHORT=$CO \
  Rscript diag_mismatch_beta3.R
done
# and the variance row, from the main run
N_CORES=12 Rscript R/07_interaction_permutation.R

# the intercept negative control (Methods)
bash run_mismatch.sh

# Table 7
METADATA_DIR=Revised/metadata Rscript clinical_assoc.R
METADATA_DIR=Revised/metadata Rscript clinical_sensitivity.R

# Table S1 and the two-dimensional matching check
cd R && Rscript 07_matching_validation.R && cd ..
Rscript test_2d.R

# figures
Rscript make_fig1.R
Rscript make_fig2.R

# software provenance
Rscript R/06_session_info.R
```

The synthetic-data scripts (`07_matching_validation.R`, `test_2d.R`) need no data
and are deterministic.

---

## Software

Recorded at the time of the analysis:

| | |
|---|---|
| limma | 3.58.1 |
| edgeR | 4.0.16 |
| DESeq2 | 1.42.1 |
| R | see `results/sessionInfo.txt` |
| kallisto | see `reference/reference_build_summary.txt` |

`Rscript R/06_session_info.R` regenerates `results/sessionInfo.txt`.

---

## Things worth knowing before reading the code

- **`make_strata` puts exact zeros in a stratum of their own** in each dimension, so
  `external2` admits up to 21 x 5 cells rather than 20 x 4. Empty cells are dropped.
- **The reference distribution for beta3 is a resampling reference, not a permutation
  of the class label.** Two mutually disjoint control sets are compared with each
  other. This instantiates the null hypothesis directly.
- **The matching implementation differs between analyses, deliberately.** Where a
  control set is drawn per patient the strata exclude that patient's own sample;
  where one cohort-level control set serves every patient no exclusion is needed;
  the differential-expression analysis builds strata from counts per million within
  the filtered set. This is stated in the Methods.
- **`diag_slope.R` draws the `NSET` control sets one at a time.** Drawing them in one
  call exhausts the strata (shortfall 209 in GSE127165).
- **Do not average beta3 within a cohort.** Where the sign is split between patients
  a mean is driven towards zero by construction.
- `nohup N_CORES=12 bash x.sh` does not work; write `N_CORES=12 nohup bash x.sh`.

## Licence and citation

See `LICENSE` and `CITATION.cff`.
