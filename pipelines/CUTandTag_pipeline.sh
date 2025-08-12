#!/usr/bin/env bash
# ==============================================================================
# CUT&Tag analysis pipeline (end-to-end, GitHub-safe)
# ==============================================================================
# • Self-contained: relative paths; no hard-coded server/usernames.
# • Runs on GNU/Linux/macOS Bash ≥ 4 (tested with bowtie2 2.5, samtools 1.19,
#   bedtools 2.31, MACS2 2.2, deepTools 3.5).
#
# DATA AVAILABILITY
# - Raw CUT&Tag sequence data: GEO Series **GSE304509** (download manually).
# - This script assumes you have already downloaded FASTQs locally into $RAW_DIR.
#
# USAGE
#   1) Edit the CONFIGURATION section below.
#   2) Put paired-end FASTQs in $RAW_DIR named as: <SAMPLE>_1.fq.gz and <SAMPLE>_2.fq.gz
#   3) Run:
#        bash pipelines/CUTandTag_pipeline.sh
#   4) Outputs go to $RESULTS_DIR/{aligned,peaks,bigwig,plots}
#
# NOTES
# - Peak calling is optional here (MACS2 example shown), since many downstream
#   figures use coverage-based heatmaps/profiles (deepTools).
# - RPGC normalization requires a correct effective genome size for your build.
#   Set $EFFECTIVE_GENOME_SIZE accordingly (integer).
# - For identical color scaling across groups in heatmaps, set --zMin/--zMax
#   (commented near plotHeatmap calls).
# ==============================================================================

set -euo pipefail

# ----------------------------------
# CONFIGURATION
# ----------------------------------
# Project & IO
PROJECT="CUTandTag_HSA"
RAW_DIR="raw_fastq"                  # where your FASTQs live (downloaded from GSE304509)
RESULTS_DIR="results"                # outputs root
THREADS="${THREADS:-8}"             # or set via environment

# Reference genome (edit to your build)
# Example placeholders — replace with your actual files/paths:
BOWTIE2_INDEX="ref/canfam_index/canfam"   # bowtie2 index prefix (e.g., ref/canfam)
GENOME_FA="ref/canfam.fa"                  # genome FASTA (for bedtools genomecov, etc.)
ANNOT_GTF="ref/canfam.gtf"                 # gene annotation (GTF) for computeMatrix
EFFECTIVE_GENOME_SIZE="2400000000"         # EDIT: integer (e.g., dog build; set correctly!)

# peak calling settings (optional, MACS2)
DO_PEAKS=false
MACS2_QVAL="0.01"
MACS2_SHIFT="-75"   # conservative CUT&Tag-like parameters (edit as needed)
MACS2_EXTSIZE="150"

# plotting switches
DO_GENE_BODY=true     # scale-regions view
DO_TSS=true           # TSS-centered view

# Color map for heatmaps
COLORMAP="RdPu"

# Sample list (base names only; no _1/_2 suffix). EDIT to your samples.
# Example: Kla mark, HSA2 0mM/25mM with 2 reps each. Add/remove as needed.
SAMPLES=(
  "HSA2_Kla_0mM_rep1"
  "HSA2_Kla_0mM_rep2"
  "HSA2_Kla_25mM_rep1"
  "HSA2_Kla_25mM_rep2"
)

# ----------------------------------
# Sanity checks
# ----------------------------------
check_bin() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not in PATH"; exit 1; }; }
for b in bowtie2 samtools bedtools macs2 bamCoverage computeMatrix plotHeatmap plotProfile; do
  check_bin "$b"
done

mkdir -p "$RESULTS_DIR/aligned" "$RESULTS_DIR/peaks" "$RESULTS_DIR/bigwig" "$RESULTS_DIR/plots" "$RESULTS_DIR/qc"

echo "==> Project: $PROJECT"
echo "==> RAW_DIR: $RAW_DIR"
echo "==> RESULTS_DIR: $RESULTS_DIR"
echo "==> THREADS: $THREADS"
echo "==> Using index: $BOWTIE2_INDEX"
echo "==> Annotation: $ANNOT_GTF"
echo "==> Effective genome size: $EFFECTIVE_GENOME_SIZE"
echo "==> Samples: ${#SAMPLES[@]}"

# ==============================================================================
# 1) ALIGNMENT + FILTERING
# ==============================================================================
for S in "${SAMPLES[@]}"; do
  R1="$RAW_DIR/${S}_1.fq.gz"
  R2="$RAW_DIR/${S}_2.fq.gz"
  [[ -s "$R1" && -s "$R2" ]] || { echo "Missing FASTQs for $S"; exit 1; }

  BAM="$RESULTS_DIR/aligned/${S}.bam"
  BAM_FILT="$RESULTS_DIR/aligned/${S}.filtered.bam"
  BAM_SORT="$RESULTS_DIR/aligned/${S}.filtered.sorted.bam"

  if [[ ! -s "$BAM_SORT" ]]; then
    echo "[ALIGN] $S"
    bowtie2       --very-sensitive-local -p "$THREADS"       -x "$BOWTIE2_INDEX" -1 "$R1" -2 "$R2"     | samtools view -bS -     > "$BAM"

    echo "[FILTER] properly paired, MAPQ≥10"
    samtools view -b -f 3 -q 10 "$BAM" > "$BAM_FILT"

    echo "[SORT/INDEX]"
    samtools sort -@ "$THREADS" -o "$BAM_SORT" "$BAM_FILT"
    samtools index "$BAM_SORT"

    # optional: remove duplicates (CUT&Tag often low-dup; leave off by default)
    # samtools markdup -r -@ "$THREADS" "$BAM_SORT" "${BAM_SORT%.bam}.dedup.bam"
    # mv "${BAM_SORT%.bam}.dedup.bam" "$BAM_SORT"
    # samtools index "$BAM_SORT"

    # QC: fragment size distribution (optional)
    # bamPEFragmentSize -b "$BAM_SORT" --hist "$RESULTS_DIR/qc/${S}_fraglen.txt" -T "$S" --samplesLabel "$S" -p "$THREADS"
  else
    echo "[SKIP] $S already aligned and sorted"
  fi
done

# ==============================================================================
# 2) (OPTIONAL) PEAK CALLING with MACS2
#    Not required for computeMatrix heatmaps/profiles. Toggle via DO_PEAKS.
# ==============================================================================
if $DO_PEAKS; then
  for S in "${SAMPLES[@]}"; do
    BAM_SORT="$RESULTS_DIR/aligned/${S}.filtered.sorted.bam"
    [[ -s "$BAM_SORT" ]] || { echo "Missing BAM for $S"; exit 1; }
    OUTPFX="$RESULTS_DIR/peaks/${S}"
    echo "[PEAKS] MACS2 callpeak: $S"
    macs2 callpeak -t "$BAM_SORT"       -f BAMPE --keep-dup all       --shift "$MACS2_SHIFT" --extsize "$MACS2_EXTSIZE"       -n "$OUTPFX" -g "$EFFECTIVE_GENOME_SIZE" -q "$MACS2_QVAL"
  done
fi

# ==============================================================================
# 3) BIGWIG (RPGC) — coverage tracks for plotting
# ==============================================================================
for S in "${SAMPLES[@]}"; do
  BAM_SORT="$RESULTS_DIR/aligned/${S}.filtered.sorted.bam"
  BW="$RESULTS_DIR/bigwig/${S}_RPGC.bw"
  if [[ ! -s "$BW" ]]; then
    echo "[BIGWIG] $S -> RPGC bigWig"
    bamCoverage -b "$BAM_SORT" -o "$BW"       --normalizeUsing RPGC       --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE"       --binSize 25 --centerReads       -p "$THREADS"
    # Note: For paired-end CUT&Tag, bamCoverage uses fragments; --extendReads not needed.
  else
    echo "[SKIP] bigWig exists: $BW"
  fi
done

# Build arrays of bigWigs and labels for deepTools
declare -a BW_FILES=()
declare -a LABELS=()
for S in "${SAMPLES[@]}"; do
  BW_FILES+=("$RESULTS_DIR/bigwig/${S}_RPGC.bw")
  LABELS+=("$S")
done

# ==============================================================================
# 4) HEATMAP/PROFILE over GENE BODIES (scale-regions)
# ==============================================================================
if $DO_GENE_BODY; then
  echo "[MATRIX] scale-regions over gene bodies"
  computeMatrix scale-regions     -S "${BW_FILES[@]}"     -R "$ANNOT_GTF"     --regionBodyLength 36000 -b 3000 -a 3000     --skipZeros --missingDataAsZero     -p "$THREADS"     -o "$RESULTS_DIR/plots/gene_body_matrix.gz"

  echo "[PLOT] heatmap (gene bodies)"
  plotHeatmap     -m "$RESULTS_DIR/plots/gene_body_matrix.gz"     -out "$RESULTS_DIR/plots/gene_body_heatmap.png"     --colorMap "$COLORMAP"     --samplesLabel "${LABELS[@]}"     --plotTitle "Signal heatmap over gene bodies"
    # To enforce identical color scales across tracks, uncomment and set:
    # --zMin 0 --zMax 5
fi

# ==============================================================================
# 5) HEATMAP/PROFILE at TSS (reference-point)
# ==============================================================================
if $DO_TSS; then
  echo "[MATRIX] reference-point at TSS"
  computeMatrix reference-point     -S "${BW_FILES[@]}"     -R "$ANNOT_GTF"     --referencePoint TSS -b 3000 -a 3000     --skipZeros --missingDataAsZero     -p "$THREADS"     -o "$RESULTS_DIR/plots/TSS_matrix_refpoint.gz"

  echo "[PLOT] heatmap (TSS)"
  plotHeatmap     -m "$RESULTS_DIR/plots/TSS_matrix_refpoint.gz"     -out "$RESULTS_DIR/plots/TSS_heatmap.png"     --colorMap "$COLORMAP"     --samplesLabel "${LABELS[@]}"     --plotTitle "Signal heatmap at TSS"
    # For identical color scaling across conditions, uncomment & set:
    # --zMin 0 --zMax 5

  echo "[PLOT] average profile (TSS)"
  plotProfile     -m "$RESULTS_DIR/plots/TSS_matrix_refpoint.gz"     -out "$RESULTS_DIR/plots/TSS_profile.png"     --samplesLabel "${LABELS[@]}"     --plotTitle "Average profile at TSS"
fi

echo "✔ Done. Outputs in: $(realpath "$RESULTS_DIR")"
