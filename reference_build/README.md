# Building the human orphan-gene kallisto reference (GRCh37)

This directory reproduces the reference-construction procedure used for the
human orphan-gene RNA-seq analysis, from the operational definition of the
orphan transcript set through construction of the kallisto index.

## Operational definition of human orphan transcripts

In this study, **human orphan transcripts** are defined operationally as the
human *de novo* transcripts reported by Ruiz-Orera et al. (2015), rather than by
rerunning a new homology search. The deposited annotation is `hsa_denovo.gtf`,
which contains the exon coordinates of the human de novo transcript set.

Source:

- Ruiz-Orera J, Hernandez-Rodriguez J, Chiva C, et al. *Origins of De Novo
  Genes in Human and Chimpanzee*. PLOS Genetics 11(12):e1005721 (2015).
  DOI: `10.1371/journal.pgen.1005721`
- Figshare dataset: `hsa_denovo.gtf`, DOI: `10.6084/m9.figshare.1604892`

The analysis is tied to **GRCh37 / Ensembl release 75**, matching the genomic
coordinate system used for the project. The Ensembl annotation and the primary
assembly genome sequence are therefore downloaded from Ensembl release 75:

- `Homo_sapiens.GRCh37.75.gtf.gz`
- `Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz`

## Why preprocessing is required

The Ruiz-Orera GTF and Ensembl GRCh37.75 do not use identical sequence-name
conventions. The orphan annotation uses UCSC-style chromosome names such as
`chr10`, whereas the Ensembl GRCh37 primary assembly uses names such as `10`.
Before the annotations are combined, the orphan GTF is therefore normalized as
follows:

- the leading `chr` prefix is removed (`chr10` -> `10`);
- mitochondrial `chrM`/`M` is converted to `MT`.

Both the Ensembl and orphan GTF files are then cleaned with `gffread -E -T`.

## Preserving orphan identifiers

To make orphan features unambiguous in every downstream file, the script
prefixes the orphan `gene_id` and `transcript_id` values with `ORPHAN_`.
For example:

```text
transcript_id "hsa_00255579"  -> transcript_id "ORPHAN_hsa_00255579"
gene_id       "XLOC_022250"   -> gene_id       "ORPHAN_XLOC_022250"
```

Consequently, kallisto output can be queried directly for IDs beginning with
`ORPHAN_`, and gene-locus aggregation retains the same explicit designation.

## Combining the Ensembl and orphan annotations

The combined reference is not produced by blindly concatenating the two GTF
files. Before merging, transcript exon structures are compared exactly.

For each transcript, an exon-structure signature is defined by:

1. chromosome;
2. strand; and
3. the ordered list of exon start/end coordinates.

If an Ensembl transcript has exactly the same exon structure as an orphan
transcript, the Ensembl copy is removed and the orphan annotation is retained.
This prevents the same transcript model from entering the kallisto index under
both an Ensembl identifier and an orphan identifier. No other Ensembl
transcripts are removed by this step.

The script records any removed Ensembl transcript IDs in
`dropped_human_transcripts.tsv`.

After the merge, the combined GTF is cleaned once more with `gffread` and a
`tx2gene.tsv` table is generated from the final annotation.

## Transcript FASTA and kallisto index

The Ensembl GRCh37.75 primary-assembly FASTA is indexed with `samtools faidx`.
`gffread` then extracts the transcript sequences described by the final combined
GTF:

```bash
gffread -w combined.transcripts.fa \
        -g Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        combined.clean.gtf
```

The kallisto index is finally constructed as:

```bash
kallisto index -i combined.idx combined.transcripts.fa
```

For the reference used in the manuscript, the expected final transcript counts
are:

```text
all transcripts:     198,507
orphan transcripts:    2,190
other transcripts:   196,317
```

The build script checks these counts and stops if they are not reproduced.
This is intentional: a mismatch usually means that a different Ensembl release,
genome FASTA, orphan annotation, or preprocessing rule has been used.

## Requirements

The executable requires:

- Bash
- `curl`
- `gzip`
- `python3`
- `gffread`
- `samtools`
- `kallisto`
- standard Unix tools (`awk`, `grep`, `sort`, `wc`, `sha256sum`)

No R package is required for reference construction.

## Run

From this directory:

```bash
chmod +x build_kallisto_reference.sh
./build_kallisto_reference.sh
```

By default, all files are written under `reference/`. A different output
directory can be supplied as the first argument:

```bash
./build_kallisto_reference.sh /path/to/reference
```

## Main outputs

```text
reference/
  source/
    Homo_sapiens.GRCh37.75.gtf.gz
    Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
    hsa_denovo.gtf
  work/
    human.clean.gtf
    orphan.nochr.gtf
    orphan.clean.gtf
    orphan.prefixed.gtf
    human.filtered.gtf
    combined.gtf
    dropped_human_transcripts.tsv
  Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
  combined.clean.gtf
  combined.transcripts.fa
  combined.idx
  tx2gene.tsv
  reference_build_summary.txt
  source_checksums.sha256
```

`reference_build_summary.txt` records transcript counts and software versions,
and `source_checksums.sha256` records checksums of the downloaded source files.

## Reproducibility note

The Figshare file is resolved by querying the public Figshare API for article
`1604892` and selecting the file whose name is exactly `hsa_denovo.gtf`. The
script therefore does not rely on a hard-coded Figshare file-download ID.
