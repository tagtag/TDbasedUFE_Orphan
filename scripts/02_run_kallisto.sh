#!/usr/bin/env bash
# Download and quantify every exact FASTQ selected by
# scripts/00_build_geo_patient_metadata.R.
#
# IMPORTANT: this script does not reconstruct an ENA URL from the run accession.
# It consumes the exact fastq_1 / fastq_2 URLs and MD5 values returned by ENA
# when metadata/sra_runs.csv was built.  This makes the downloaded files part of
# the auditable metadata-selection step rather than a second implicit decision.
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT_DIR"

THREADS=${1:-8}
REF_DIR=${REF_DIR:-data/reference}
FASTQ_DIR=${FASTQ_DIR:-data/fastq}
OUT_DIR=${OUT_DIR:-data/kallisto}
INDEX="$REF_DIR/combined_GRCh37.idx"
RUNS=${RUNS:-metadata/sra_runs.csv}
SE_FRAGMENT_LENGTH=${SE_FRAGMENT_LENGTH:-200}
SE_FRAGMENT_SD=${SE_FRAGMENT_SD:-30}
KEEP_FASTQ=${KEEP_FASTQ:-0}

mkdir -p "$FASTQ_DIR" "$OUT_DIR"
[ -f "$INDEX" ] || { echo "index not found: $INDEX" >&2; exit 1; }
[ -f "$RUNS" ] || { echo "run manifest not found: $RUNS" >&2; exit 1; }
command -v wget >/dev/null || { echo "wget is required" >&2; exit 1; }
command -v kallisto >/dev/null || { echo "kallisto is required" >&2; exit 1; }

header=$(head -n 1 "$RUNS" | tr -d '\r')
expected='run_accession,cohort,gsm,patient_id,condition,library_layout,fastq_1,fastq_2,md5_1,md5_2'
if [ "$header" != "$expected" ]; then
  echo "Unexpected schema in $RUNS" >&2
  echo "Expected: $expected" >&2
  echo "Observed: $header" >&2
  echo "Rebuild it with: Rscript scripts/00_build_geo_patient_metadata.R" >&2
  exit 1
fi

check_md5 () {
  local file=$1
  local expected_md5=${2:-}
  [ -z "$expected_md5" ] && return 0
  if command -v md5sum >/dev/null; then
    local observed
    observed=$(md5sum "$file" | awk '{print $1}')
    if [ "${observed,,}" != "${expected_md5,,}" ]; then
      echo "MD5 mismatch: $file" >&2
      echo "  expected: $expected_md5" >&2
      echo "  observed: $observed" >&2
      return 1
    fi
  fi
}

download_one () {
  local url=$1
  local dest=$2
  local md5=${3:-}
  [ -n "$url" ] || { echo "empty FASTQ URL for $dest" >&2; return 1; }
  if [ ! -s "$dest" ]; then
    wget -c -O "$dest" "$url"
  fi
  check_md5 "$dest" "$md5"
}

# write.csv(..., quote=FALSE) is used by the R metadata builder.  The fields
# written here (accessions, identifiers, URLs and MD5 strings) do not contain
# commas, so shell CSV parsing is deterministic.
tail -n +2 "$RUNS" | tr -d '\r' | \
while IFS=, read -r run cohort gsm patient condition layout fq1 fq2 md51 md52; do
  [ -n "$run" ] || continue
  out="$OUT_DIR/${cohort}_${run}"

  if [ -s "$out/abundance.tsv" ]; then
    echo "[skip] $cohort $run: abundance.tsv already exists"
    continue
  fi

  echo "[$cohort] $run  GSM=${gsm:-NA} patient=${patient:-NA} condition=${condition:-NA} layout=$layout"

  if [ "$layout" = "PAIRED" ]; then
    [ -n "$fq1" ] && [ -n "$fq2" ] || {
      echo "paired run $run does not have exactly two FASTQ URLs in $RUNS" >&2
      exit 1
    }
    f1="$FASTQ_DIR/${run}_1.fastq.gz"
    f2="$FASTQ_DIR/${run}_2.fastq.gz"
    download_one "$fq1" "$f1" "$md51"
    download_one "$fq2" "$f2" "$md52"
    kallisto quant -i "$INDEX" -o "$out" -t "$THREADS" "$f1" "$f2"
    if [ "$KEEP_FASTQ" != "1" ]; then rm -f "$f1" "$f2"; fi

  elif [ "$layout" = "SINGLE" ]; then
    [ -n "$fq1" ] || {
      echo "single-end run $run has no FASTQ URL in $RUNS" >&2
      exit 1
    }
    f1="$FASTQ_DIR/${run}.fastq.gz"
    download_one "$fq1" "$f1" "$md51"
    kallisto quant -i "$INDEX" -o "$out" -t "$THREADS" --single \
      -l "$SE_FRAGMENT_LENGTH" -s "$SE_FRAGMENT_SD" "$f1"
    if [ "$KEEP_FASTQ" != "1" ]; then rm -f "$f1"; fi

  else
    echo "unknown library_layout '$layout' for $run" >&2
    exit 1
  fi
done

echo "quantification complete"
