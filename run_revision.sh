#!/usr/bin/env bash
# Revised analysis, end to end, from transcript-level abundance files onward.
#
#   KALLISTO_ROOT=Revised METADATA_DIR=Revised/metadata N_CORES=12 bash run_revision.sh
#
# Each stage writes its own logs; see analysis_manifest.tsv for what produces what.
# The three run_*.sh scripts are the heavy stages and can be run separately.
set -euo pipefail

N_CORES=${N_CORES:-1}
METADATA_DIR=${METADATA_DIR:-Revised/metadata}
REVISED_ROOT=${KALLISTO_ROOT:-Revised}
COHORTS=${COHORTS:-"GSE244679 GSE127165 GSE144269"}

echo "=============== Table 1 and Table 4 ==============="
N_CORES=$N_CORES bash run_batch_audit.sh

echo "=============== Table 2 and Table 5 ==============="
N_CORES=$N_CORES bash run_verify.sh

echo "=============== Table 3 ==============="
for CO in $COHORTS; do
  echo "--- $CO ---"
  N_CORES=$N_CORES N_PERM=200 TRANSFORMS=rank NMIS=5 \
  KALLISTO_ROOT="$REVISED_ROOT/$CO" METADATA_DIR="$METADATA_DIR" COHORT="$CO" \
  Rscript diag_mismatch_beta3.R
done
echo "--- between-patient variance, from the main run (N_PERM=1000) ---"
N_CORES=$N_CORES Rscript R/07_interaction_permutation.R

echo "=============== intercept negative control (Methods) ==============="
N_CORES=$N_CORES bash run_mismatch.sh

echo "=============== Table 7 ==============="
METADATA_DIR="$METADATA_DIR" Rscript clinical_assoc.R
METADATA_DIR="$METADATA_DIR" Rscript clinical_sensitivity.R

echo "=============== Table S1 and the synthetic checks ==============="
( cd R && Rscript 07_matching_validation.R )
Rscript test_2d.R

echo "=============== figures ==============="
Rscript make_fig1.R
Rscript make_fig2.R

echo "=============== software provenance ==============="
Rscript R/06_session_info.R

echo
echo "Done. See analysis_manifest.tsv for the table-to-script mapping."
