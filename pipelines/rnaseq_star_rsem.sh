#!/usr/bin/env bash
# ==============================================================================
# RNA-seq pipeline: STAR (align) -> RSEM (quant)
# Self-contained, GitHub-safe, no auto-download of data.
# ==============================================================================
# DATA AVAILABILITY
# - Raw RNA-seq FASTQs: GEO **GSE304507** (download manually).
#   Place FASTQs in ./raw/ as <SAMPLE>_1.fq.gz and <SAMPLE>_2.fq.gz
#
# WHAT THIS DOES
# 1) (Optional) build STAR/RSEM reference (if PREPARE_REFERENCE=true)
# 2) Map with STAR (sorted BAM + transcriptome BAM for RSEM)
# 3) Quantify with RSEM (strand-aware)
# 4) Merge per-sample .genes.results into:
#    - counts matrix:   results/counts/GeneExpressionMatrix.tsv
#    - TPM matrix:      results/counts/GeneExpressionMatrixTPM.tsv
#
# REQUIREMENTS in PATH: STAR, RSEM, samtools, awk, paste, sort, grep
# Recommended versions (used in our runs): STAR ≥2.7, RSEM ≥1.3.0, samtools ≥1.10
# ==============================================================================

set -euo pipefail

# -------------------------
# CONFIGURATION (EDIT ME)
# -------------------------
THREADS="${THREADS:-12}"

# Samples (base names only; no _1/_2 suffix)
SAMPLES=(
  HSA2_0mM_1
  HSA2_0mM_2
  HSA2_0mM_3
  HSA2_25mM_1
  HSA2_25mM_2
  HSA2_25mM_3
  HSA3_0mM_1
  HSA3_0mM_2
  HSA3_0mM_3
  HSA3_25mM_1
  HSA3_25mM_2
  HSA3_25mM_3
)

# I/O layout (relative to repo root)
RAW_DIR="raw"                 # FASTQs from GSE304507
REF_DIR="ref"                 # FASTA/GTF live here
INDEX_DIR="rsem-star-index"   # STAR/RSEM reference prefix directory
OUT_DIR="results"
STAR_DIR="${OUT_DIR}/star"
RSEM_DIR="${OUT_DIR}/rsem"
MATRIX_DIR="${OUT_DIR}/counts"
mkdir -p "$STAR_DIR" "$RSEM_DIR" "$MATRIX_DIR" "$REF_DIR" "$INDEX_DIR"

# Reference (edit to your build; example: Ensembl CanFam3.1 r104)
FASTA="${REF_DIR}/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa"
GTF="${REF_DIR}/Canis_lupus_familiaris.CanFam3.1.104.gtf"
RSEM_PREFIX="${INDEX_DIR}/Ensembl_CanFam3.1_r104"   # will create files like <prefix>.grp, <prefix>.n2g.idx, etc.

# If you want the script to build the reference locally, set true and ensure FASTA/GTF exist.
PREPARE_REFERENCE=false

# Library strandedness for RSEM: one of {none, forward, reverse}
RSEM_STRAND="reverse"

# STAR options (array, safe quoting)
STAR_ARGS=(
  --genomeLoad NoSharedMemory
  --outSAMtype BAM SortedByCoordinate
  --quantMode TranscriptomeSAM
  --runThreadN "${THREADS}"
  --outSAMattributes All
  --readFilesCommand zcat
  --outSAMstrandField intronMotif
)

# -------------------------
# Sanity checks
# -------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }; }
for b in STAR rsem-prepare-reference rsem-calculate-expression samtools awk paste sort grep; do need "$b"; done

echo "==> Samples: ${#SAMPLES[@]}"
echo "==> FASTQs: ${RAW_DIR}"
echo "==> Reference FASTA: ${FASTA}"
echo "==> Reference GTF:   ${GTF}"
echo "==> RSEM/STAR prefix: ${RSEM_PREFIX}"
echo "==> Threads: ${THREADS}"
echo "==> RSEM strandedness: ${RSEM_STRAND}"
echo

# -------------------------
# (Optional) Build reference
# -------------------------
if $PREPARE_REFERENCE; then
  [[ -s "$FASTA" && -s "$GTF" ]] || {
    echo "ERROR: FASTA/GTF not found in ${REF_DIR}. Download them manually (e.g., Ensembl)."
    echo "Example (commented):"
    echo "  # wget https://ftp.ensembl.org/pub/release-104/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz"
    echo "  # wget https://ftp.ensembl.org/pub/release-104/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.CanFam3.1.104.gtf.gz"
    echo "  # gunzip -c ...fa.gz  > ${FASTA}"
    echo "  # gunzip -c ...gtf.gz > ${GTF}"
    exit 1
  }
  echo "==> Building STAR/RSEM reference at ${RSEM_PREFIX}"
  rsem-prepare-reference --star -p "${THREADS}" --gtf "${GTF}" "${FASTA}" "${RSEM_PREFIX}"
  echo "Reference build complete."
fi

# Ensure index exists
[[ -s "${RSEM_PREFIX}.grp" ]] || { echo "ERROR: RSEM/STAR index not found at ${RSEM_PREFIX}. Set PREPARE_REFERENCE=true or point to an existing prefix."; exit 1; }

# -------------------------
# STAR mapping
# -------------------------
for s in "${SAMPLES[@]}"; do
  R1="${RAW_DIR}/${s}_1.fq.gz"
  R2="${RAW_DIR}/${s}_2.fq.gz"
  [[ -s "$R1" && -s "$R2" ]] || { echo "ERROR: Missing FASTQs for ${s} in ${RAW_DIR}"; exit 1; }

  OUTPFX="${STAR_DIR}/${s}."
  GENOMEDIR="$(dirname "${RSEM_PREFIX}")"   # STAR genomeDir is the directory containing the indices made by rsem-prepare-reference

  if [[ ! -s "${STAR_DIR}/${s}.Aligned.toTranscriptome.out.bam" ]]; then
    echo "[STAR] ${s}"
    STAR       "${STAR_ARGS[@]}"       --genomeDir "${GENOMEDIR}"       --readFilesIn "${R1}" "${R2}"       --outFileNamePrefix "${OUTPFX}"
    # Make sure BAM is indexed (STAR writes a sorted genomic BAM)
    samtools index -@ "${THREADS}" "${STAR_DIR}/${s}Aligned.sortedByCoord.out.bam"
  else
    echo "[SKIP] STAR outputs exist for ${s}"
  fi
done

# -------------------------
# RSEM quantification
# -------------------------
for s in "${SAMPLES[@]}"; do
  TBAM="${STAR_DIR}/${s}.Aligned.toTranscriptome.out.bam"
  [[ -s "$TBAM" ]] || { echo "ERROR: Missing transcriptome BAM for ${s}: ${TBAM}"; exit 1; }

  if [[ ! -s "${RSEM_DIR}/${s}.genes.results" ]]; then
    echo "[RSEM] ${s}"
    rsem-calculate-expression       --paired-end --alignments       --estimate-rspd       --strandedness "${RSEM_STRAND}"       --no-bam-output       -p "${THREADS}"       "${TBAM}"       "${RSEM_PREFIX}"       "${RSEM_DIR}/${s}"
  else
    echo "[SKIP] RSEM outputs exist for ${s}"
  fi
done

# -------------------------
# Merge matrices
# -------------------------
echo "[MERGE] counts matrix"
# rsem-generate-data-matrix expects to be run in the directory holding *.genes.results
pushd "${RSEM_DIR}" >/dev/null
rsem-generate-data-matrix $(printf "%s.genes.results " "${SAMPLES[@]}") > "${MATRIX_DIR}/GeneExpressionMatrix.tsv"
popd >/dev/null

echo "[MERGE] TPM matrix"
# Build TPM matrix without external tools; assumes identical gene ordering across samples (true for RSEM with same reference).
TMPDIR="$(mktemp -d)"
# Extract gene_id list from the first sample (skip header)
awk 'NR>1{print $1}' "${RSEM_DIR}/${SAMPLES[0]}.genes.results" > "${TMPDIR}/gene_id.txt"
# Verify order matches across all samples
for s in "${SAMPLES[@]}"; do
  awk 'NR>1{print $1}' "${RSEM_DIR}/${s}.genes.results" > "${TMPDIR}/${s}.ids"
  if ! cmp -s "${TMPDIR}/gene_id.txt" "${TMPDIR}/${s}.ids"; then
    echo "ERROR: Gene order mismatch in ${s}. Consider merging TPM in R instead." >&2
    rm -rf "${TMPDIR}"; exit 1
  fi
  # Extract TPM column (column named "TPM"); in RSEM .genes.results it's column 6.
  awk 'NR==1{next} {print $6}' "${RSEM_DIR}/${s}.genes.results" > "${TMPDIR}/${s}.tpm"
done
# Header
{
  printf "gene_id"
  for s in "${SAMPLES[@]}"; do printf "\t%s" "${s}"; done
  printf "\n"
} > "${MATRIX_DIR}/GeneExpressionMatrixTPM.tsv"
# Paste columns
paste "${TMPDIR}/gene_id.txt" $(printf "%s " $(for s in "${SAMPLES[@]}"; do echo "${TMPDIR}/${s}.tpm"; done)) >> "${MATRIX_DIR}/GeneExpressionMatrixTPM.tsv"
rm -rf "${TMPDIR}"

echo "✔ Done."
echo "Outputs:"
echo "  STAR BAMs:        ${STAR_DIR}/*Aligned.sortedByCoord.out.bam"
echo "  RSEM results:     ${RSEM_DIR}/*.genes.results"
echo "  Counts matrix:    ${MATRIX_DIR}/GeneExpressionMatrix.tsv"
echo "  TPM matrix:       ${MATRIX_DIR}/GeneExpressionMatrixTPM.tsv"
