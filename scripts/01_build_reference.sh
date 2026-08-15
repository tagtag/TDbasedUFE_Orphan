#!/usr/bin/env bash
# Build the combined GRCh37 reference used by the analysis.
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT_DIR"

REF_DIR=${REF_DIR:-data/reference}
mkdir -p "$REF_DIR"

ENSEMBL_GTF="$REF_DIR/Homo_sapiens.GRCh37.75.gtf"
DENOVO_RAW="$REF_DIR/hsa_denovo.gtf"
DENOVO_PREFIXED="$REF_DIR/hsa_denovo.ORPHAN_prefixed.gtf"
GENOME_FA="$REF_DIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
COMBINED_GTF="$REF_DIR/combined_GRCh37.gtf"
TRANSCRIPT_FA="$REF_DIR/combined_GRCh37.transcripts.fa"
INDEX="$REF_DIR/combined_GRCh37.idx"

if [ ! -f "$ENSEMBL_GTF" ]; then
  wget -c https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -O "$ENSEMBL_GTF.gz"
  gunzip -kf "$ENSEMBL_GTF.gz"
fi

if [ ! -f "$GENOME_FA" ]; then
  wget -c https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -O "$GENOME_FA.gz"
  gunzip -kf "$GENOME_FA.gz"
fi

if [ ! -f "$DENOVO_RAW" ]; then
  echo "Place Ruiz-Orera et al. hsa_denovo.gtf at: $DENOVO_RAW"
  echo "Source: https://doi.org/10.6084/m9.figshare.1604892"
  exit 1
fi

Rscript scripts/00_prefix_orphan_gtf.R "$DENOVO_RAW" "$DENOVO_PREFIXED"
cat "$ENSEMBL_GTF" "$DENOVO_PREFIXED" > "$COMBINED_GTF"

samtools faidx "$GENOME_FA"
gffread "$COMBINED_GTF" -g "$GENOME_FA" -w "$TRANSCRIPT_FA"

N_TX=$(grep -c '^>' "$TRANSCRIPT_FA")
echo "transcripts in reference: $N_TX (manuscript: 198507)"
if [ "$N_TX" -ne 198507 ]; then
  echo "WARNING: reference size differs from manuscript; verify GTF versions." >&2
fi

kallisto index -i "$INDEX" "$TRANSCRIPT_FA"
echo "index written to $INDEX"
echo "next: Rscript scripts/03_extract_orphan_ids.R $DENOVO_PREFIXED"
