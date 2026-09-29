#!/usr/bin/env bash
set -euo pipefail

# Analysis runner from transcript-level abundance files onward.
#
# Metadata/run selection can be reconstructed entirely from public GEO/ENA and
# the GSE127165 public supplementary table.  No hand-written run list,
# normal/disease mapping, or clinical override is used.  By default an existing
# validated metadata set is reused; set REBUILD_GEO_METADATA=1 to force a fresh
# public-data reconstruction.

KROOT=${KALLISTO_ROOT:-data/kallisto}
GTF=${COMBINED_GTF:-data/reference/combined_GRCh37.gtf}
META_OUT=${METADATA_OUTPUT_DIR:-metadata_output}
META_FILE=${GEO_PATIENT_METADATA_FILE:-${META_OUT}/GEO_patient_metadata_combined.tsv}
REBUILD=${REBUILD_GEO_METADATA:-0}

required_meta=(
  "$META_FILE"
  "metadata/sra_runs.csv"
  "metadata/sample_sheet_GSE244679.csv"
  "metadata/sample_sheet_GSE127165.csv"
  "metadata/sample_sheet_GSE144269.csv"
  "metadata/sample_sheet_GSE40419.csv"
)

metadata_ready=1
for f in "${required_meta[@]}"; do
  if [[ ! -s "$f" ]] || [[ $(wc -l < "$f") -le 1 ]]; then
    metadata_ready=0
  fi
done

if [[ "$REBUILD" == "1" || "$metadata_ready" == "0" ]]; then
  echo "=============== public GEO/ENA metadata reconstruction ==============="
  Rscript scripts/00_build_geo_patient_metadata.R
else
  echo "=============== metadata: reuse validated files ==============="
  echo "Set REBUILD_GEO_METADATA=1 to reconstruct them from public sources."
fi

echo "=============== transcript -> gene loci ==============="
Rscript scripts/04_aggregate_gene_level.R "$GTF" "$KROOT"

for f in \
  R/01_table1_relative_expression.R \
  R/02_patient_level_analysis.R \
  R/03_table5_control_sets.R \
  R/04_clinical_association.R \
  R/05_figure2.R \
  R/06_session_info.R \
  R/07_validate_manuscript_counts.R; do
  echo "=============== $f ==============="
  Rscript "$f"
done
