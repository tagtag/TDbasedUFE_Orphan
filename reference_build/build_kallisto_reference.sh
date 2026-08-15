#!/usr/bin/env bash
# Reproduce the GRCh37 combined Ensembl + human orphan transcript reference
# and build the kallisto index used in the orphan-gene RNA-seq analysis.
set -euo pipefail

OUTDIR=${1:-reference}
SOURCE_DIR="$OUTDIR/source"
WORK_DIR="$OUTDIR/work"
mkdir -p "$SOURCE_DIR" "$WORK_DIR"

ENSEMBL_GTF_GZ="$SOURCE_DIR/Homo_sapiens.GRCh37.75.gtf.gz"
ENSEMBL_GTF_RAW="$WORK_DIR/human.raw.gtf"
ENSEMBL_GTF_CLEAN="$WORK_DIR/human.clean.gtf"
GENOME_FA_GZ="$SOURCE_DIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
GENOME_FA="$OUTDIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
ORPHAN_RAW="$SOURCE_DIR/hsa_denovo.gtf"
ORPHAN_NOCHR="$WORK_DIR/orphan.nochr.gtf"
ORPHAN_CLEAN="$WORK_DIR/orphan.clean.gtf"
ORPHAN_PREFIXED="$WORK_DIR/orphan.prefixed.gtf"
HUMAN_FILTERED="$WORK_DIR/human.filtered.gtf"
COMBINED_RAW="$WORK_DIR/combined.gtf"
COMBINED_CLEAN="$OUTDIR/combined.clean.gtf"
DROPPED="$WORK_DIR/dropped_human_transcripts.tsv"
TX2GENE="$OUTDIR/tx2gene.tsv"
TRANSCRIPT_FA="$OUTDIR/combined.transcripts.fa"
KALLISTO_INDEX="$OUTDIR/combined.idx"
SUMMARY="$OUTDIR/reference_build_summary.txt"
CHECKSUMS="$OUTDIR/source_checksums.sha256"

ENSEMBL_GTF_URL="https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
GENOME_FA_URL="https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
FIGSHARE_API="https://api.figshare.com/v2/articles/1604892"

EXPECTED_ORPHAN=2190
EXPECTED_OTHER=196317
EXPECTED_TOTAL=198507

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: $1" >&2
    exit 1
  }
}

for cmd in curl gzip python3 gffread samtools kallisto awk grep sort wc sha256sum; do
  need_cmd "$cmd"
done

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "[source] already present: $out"
  else
    echo "[source] downloading: $url"
    curl -fL --retry 5 --retry-delay 3 "$url" -o "$out"
  fi
}

count_gtf_transcripts() {
  python3 - "$1" <<'PY'
import re, sys
p = sys.argv[1]
ids = set()
pat = re.compile(r'transcript_id\s+"([^"]+)"')
with open(p, encoding='utf-8') as fh:
    for line in fh:
        if not line or line.startswith('#'):
            continue
        m = pat.search(line)
        if m:
            ids.add(m.group(1))
print(len(ids))
PY
}

count_orphan_gtf_transcripts() {
  python3 - "$1" <<'PY'
import re, sys
p = sys.argv[1]
ids = set()
pat = re.compile(r'transcript_id\s+"([^"]+)"')
with open(p, encoding='utf-8') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        m = pat.search(line)
        if m and m.group(1).startswith('ORPHAN_'):
            ids.add(m.group(1))
print(len(ids))
PY
}

echo "== 1. Download Ensembl release 75 / GRCh37 sources =="
fetch "$ENSEMBL_GTF_URL" "$ENSEMBL_GTF_GZ"
fetch "$GENOME_FA_URL" "$GENOME_FA_GZ"

echo "== 2. Resolve and download the Ruiz-Orera hsa_denovo.gtf from Figshare =="
if [[ ! -s "$ORPHAN_RAW" ]]; then
  FIGSHARE_JSON=$(curl -fsSL --retry 5 --retry-delay 3 "$FIGSHARE_API")
  ORPHAN_URL=$(python3 -c '
import json, sys
obj=json.load(sys.stdin)
hits=[f.get("download_url") for f in obj.get("files", []) if f.get("name")=="hsa_denovo.gtf"]
if len(hits)!=1 or not hits[0]:
    raise SystemExit("Could not uniquely resolve hsa_denovo.gtf from Figshare article 1604892")
print(hits[0])
' <<< "$FIGSHARE_JSON")
  fetch "$ORPHAN_URL" "$ORPHAN_RAW"
else
  echo "[source] already present: $ORPHAN_RAW"
fi

printf '%s  %s\n' "$(sha256sum "$ENSEMBL_GTF_GZ" | awk '{print $1}')" "$ENSEMBL_GTF_GZ" > "$CHECKSUMS"
printf '%s  %s\n' "$(sha256sum "$GENOME_FA_GZ" | awk '{print $1}')" "$GENOME_FA_GZ" >> "$CHECKSUMS"
printf '%s  %s\n' "$(sha256sum "$ORPHAN_RAW" | awk '{print $1}')" "$ORPHAN_RAW" >> "$CHECKSUMS"

echo "== 3. Clean the Ensembl GTF =="
gzip -dc "$ENSEMBL_GTF_GZ" > "$ENSEMBL_GTF_RAW"
gffread -E "$ENSEMBL_GTF_RAW" -T -o "$ENSEMBL_GTF_CLEAN"

echo "== 4. Normalize orphan chromosome names and clean the orphan GTF =="
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     {
       sub(/^chr/, "", $1)
       if ($1=="M") $1="MT"
       print
     }' "$ORPHAN_RAW" > "$ORPHAN_NOCHR"
gffread -E "$ORPHAN_NOCHR" -T -o "$ORPHAN_CLEAN"

echo "== 5. Prefix orphan gene_id and transcript_id with ORPHAN_ =="
python3 - "$ORPHAN_CLEAN" "$ORPHAN_PREFIXED" <<'PY'
import re, sys
src, dst = sys.argv[1], sys.argv[2]
patterns = {
    "gene_id": re.compile(r'(gene_id\s+")([^"]+)(")'),
    "transcript_id": re.compile(r'(transcript_id\s+")([^"]+)(")'),
}
def pref(line, key):
    p = patterns[key]
    def repl(m):
        value = m.group(2)
        if not value.startswith("ORPHAN_"):
            value = "ORPHAN_" + value
        return m.group(1) + value + m.group(3)
    return p.sub(repl, line, count=1)
with open(src, encoding="utf-8") as fi, open(dst, "w", encoding="utf-8") as fo:
    for line in fi:
        if not line.startswith("#"):
            line = pref(line, "gene_id")
            line = pref(line, "transcript_id")
        fo.write(line)
PY

ORPHAN_N=$(count_orphan_gtf_transcripts "$ORPHAN_PREFIXED")
if [[ "$ORPHAN_N" -ne "$EXPECTED_ORPHAN" ]]; then
  echo "ERROR: expected $EXPECTED_ORPHAN orphan transcripts, found $ORPHAN_N" >&2
  exit 1
fi

echo "== 6. Remove exact Ensembl transcript models duplicated by orphan models =="
# Exact duplicate = same chromosome, strand, and sorted exon start/end intervals.
# The orphan annotation is retained; only the matching Ensembl transcript is removed.
python3 - "$ENSEMBL_GTF_CLEAN" "$ORPHAN_PREFIXED" "$HUMAN_FILTERED" "$DROPPED" <<'PY'
import re, sys
from collections import defaultdict
human, orphan, human_out, dropped_out = sys.argv[1:]
pat = re.compile(r'transcript_id\s+"([^"]+)"')

def exon_signatures(path):
    ex = defaultdict(list)
    chrom = {}
    strand = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            m = pat.search(f[8])
            if not m:
                continue
            tid = m.group(1)
            ex[tid].append((int(f[3]), int(f[4])))
            chrom[tid] = f[0]
            strand[tid] = f[6]
    sig = {}
    for tid, intervals in ex.items():
        sig[tid] = (chrom[tid], strand[tid], tuple(sorted(intervals)))
    return sig

hsig = exon_signatures(human)
osig = exon_signatures(orphan)
orph_sigs = set(osig.values())
drop = {tid for tid, sig in hsig.items() if sig in orph_sigs}

with open(dropped_out, 'w', encoding='utf-8') as fo:
    fo.write('transcript_id\n')
    for tid in sorted(drop):
        fo.write(tid + '\n')

with open(human, encoding='utf-8') as fi, open(human_out, 'w', encoding='utf-8') as fo:
    for line in fi:
        if line.startswith('#'):
            fo.write(line)
            continue
        m = pat.search(line)
        if m and m.group(1) in drop:
            continue
        fo.write(line)

print(f'exact exon-structure Ensembl transcripts removed: {len(drop)}')
PY

echo "== 7. Merge annotations and clean the combined GTF =="
cat "$HUMAN_FILTERED" "$ORPHAN_PREFIXED" > "$COMBINED_RAW"
gffread -E "$COMBINED_RAW" -T -o "$COMBINED_CLEAN"

echo "== 8. Build transcript-to-gene table from the final combined GTF =="
python3 - "$COMBINED_CLEAN" "$TX2GENE" <<'PY'
import re, sys
src, dst = sys.argv[1:]
gpat = re.compile(r'gene_id\s+"([^"]+)"')
tpat = re.compile(r'transcript_id\s+"([^"]+)"')
pairs = set()
with open(src, encoding='utf-8') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        gm = gpat.search(line)
        tm = tpat.search(line)
        if gm and tm:
            pairs.add((tm.group(1), gm.group(1)))
with open(dst, 'w', encoding='utf-8') as fo:
    fo.write('transcript_id\tgene_id\n')
    for tid, gid in sorted(pairs):
        fo.write(f'{tid}\t{gid}\n')
PY

echo "== 9. Prepare the GRCh37 primary-assembly genome =="
gzip -dc "$GENOME_FA_GZ" > "$GENOME_FA"
samtools faidx "$GENOME_FA"

echo "== 10. Extract transcript sequences =="
gffread -w "$TRANSCRIPT_FA" -g "$GENOME_FA" "$COMBINED_CLEAN"

echo "== 11. Verify reference composition =="
TOTAL_N=$(grep -c '^>' "$TRANSCRIPT_FA")
ORPHAN_FASTA_N=$(grep '^>' "$TRANSCRIPT_FA" | sed 's/^>//' | awk '{print $1}' | grep -c '^ORPHAN_' || true)
OTHER_N=$((TOTAL_N - ORPHAN_FASTA_N))

printf 'total_transcripts\t%s\n' "$TOTAL_N"
printf 'orphan_transcripts\t%s\n' "$ORPHAN_FASTA_N"
printf 'other_transcripts\t%s\n' "$OTHER_N"

if [[ "$ORPHAN_FASTA_N" -ne "$EXPECTED_ORPHAN" || "$OTHER_N" -ne "$EXPECTED_OTHER" || "$TOTAL_N" -ne "$EXPECTED_TOTAL" ]]; then
  echo "ERROR: final reference composition does not match the manuscript." >&2
  echo "Expected total/orphan/other = $EXPECTED_TOTAL/$EXPECTED_ORPHAN/$EXPECTED_OTHER" >&2
  echo "Observed total/orphan/other = $TOTAL_N/$ORPHAN_FASTA_N/$OTHER_N" >&2
  exit 1
fi

echo "== 12. Build kallisto index =="
kallisto index -i "$KALLISTO_INDEX" "$TRANSCRIPT_FA"

echo "== 13. Record build summary =="
{
  echo "Human orphan-gene kallisto reference build"
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "ensembl_release=75"
  echo "genome_assembly=GRCh37"
  echo "orphan_source=Ruiz-Orera_et_al_2015_figshare_1604892"
  echo "total_transcripts=$TOTAL_N"
  echo "orphan_transcripts=$ORPHAN_FASTA_N"
  echo "other_transcripts=$OTHER_N"
  echo "kallisto_index=$KALLISTO_INDEX"
  echo
  echo "Software versions"
  printf 'kallisto: '; kallisto version 2>&1 | head -1 || true
  printf 'gffread: '; gffread --version 2>&1 | head -1 || true
  printf 'samtools: '; samtools --version 2>&1 | head -1 || true
  printf 'python3: '; python3 --version 2>&1 | head -1 || true
} > "$SUMMARY"

echo
echo "Reference build completed successfully."
echo "Combined GTF : $COMBINED_CLEAN"
echo "Transcript FASTA: $TRANSCRIPT_FA"
echo "kallisto index: $KALLISTO_INDEX"
echo "tx2gene table: $TX2GENE"
echo "summary: $SUMMARY"
