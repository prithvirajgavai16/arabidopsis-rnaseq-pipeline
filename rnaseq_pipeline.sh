#!/bin/bash
# RNA-seq pipeline: HISAT2 + featureCounts for Arabidopsis thaliana (SRP439187)
# Single-end, 50bp, Illumina HiSeq 1500

set -euo pipefail

# ─── Configuration ───────────────────────────────────────────────────────────
THREADS=8
SRA_DIR="/home/prithvirajgavai/SRP439187_raw"
BASE_DIR="/home/prithvirajgavai/rnaseq_pipeline"
FASTQ_DIR="$BASE_DIR/fastq"
TRIM_DIR="$BASE_DIR/trimmed"
QC_DIR="$BASE_DIR/qc"
ALIGN_DIR="$BASE_DIR/aligned"
COUNTS_DIR="$BASE_DIR/counts"
LOG_DIR="$BASE_DIR/logs"
REF_DIR="$BASE_DIR/reference"
INDEX_DIR="$REF_DIR/hisat2_index"

GENOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz"
GENOME_FA="$REF_DIR/TAIR10.fa"
GTF="$REF_DIR/TAIR10.gtf"
INDEX_PREFIX="$INDEX_DIR/tair10"

SRA_TOOLS="/home/prithvirajgavai/sratoolkit.3.4.1-ubuntu64/bin"
export PATH="$SRA_TOOLS:$PATH"

ACCESSIONS=(
  SRR24710315 SRR24710316 SRR24710317 SRR24710318 SRR24710319
  SRR24710320 SRR24710321 SRR24710322 SRR24710323 SRR24710324
  SRR24710325 SRR24710326 SRR24710327
)

# ─── Helpers ─────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_DIR/pipeline.log"; }
check_tool() { command -v "$1" &>/dev/null || { echo "ERROR: $1 not found. Run: bash install_tools.sh"; exit 1; }; }

# ─── Step 0: Check tools ─────────────────────────────────────────────────────
check_tools() {
  log "Checking required tools..."
  for tool in hisat2 samtools fastqc trim_galore featureCounts multiqc fasterq-dump; do
    check_tool "$tool"
  done
  log "All tools found."
}

# ─── Step 1: Create directories ───────────────────────────────────────────────
setup_dirs() {
  mkdir -p "$LOG_DIR"
  log "Creating directory structure..."
  mkdir -p "$FASTQ_DIR" "$TRIM_DIR" "$QC_DIR/raw" "$QC_DIR/trimmed" \
           "$ALIGN_DIR" "$COUNTS_DIR" "$LOG_DIR" "$REF_DIR" "$INDEX_DIR"
}

# ─── Step 2: Download reference genome & annotation ──────────────────────────
download_reference() {
  if [[ -f "$GENOME_FA" && -f "$GTF" ]]; then
    log "Reference files already exist, skipping download."
    return
  fi
  log "Downloading TAIR10 genome and GTF annotation..."
  wget -q --show-progress -O "${GENOME_FA}.gz" "$GENOME_URL"
  wget -q --show-progress -O "${GTF}.gz" "$GTF_URL"
  gunzip "${GENOME_FA}.gz"
  gunzip "${GTF}.gz"
  log "Reference download complete."
}

# ─── Step 3: Build HISAT2 index ──────────────────────────────────────────────
build_index() {
  if [[ -f "${INDEX_PREFIX}.1.ht2" ]]; then
    log "HISAT2 index already exists, skipping."
    return
  fi
  log "Building HISAT2 index (this takes ~10 min)..."
  hisat2-build -p "$THREADS" "$GENOME_FA" "$INDEX_PREFIX" \
    > "$LOG_DIR/hisat2_build.log" 2>&1
  log "Index build complete."
}

# ─── Step 4: SRA to FASTQ ────────────────────────────────────────────────────
sra_to_fastq() {
  log "Converting SRA files to FASTQ..."
  for acc in "${ACCESSIONS[@]}"; do
    local fq="$FASTQ_DIR/${acc}.fastq.gz"
    if [[ -f "$fq" ]]; then
      log "  $acc: already converted, skipping."
      continue
    fi
    log "  Converting $acc..."
    fasterq-dump "$SRA_DIR/$acc" \
      --outdir "$FASTQ_DIR" \
      --threads "$THREADS" \
      --progress \
      2>> "$LOG_DIR/fasterq_${acc}.log"
    gzip "$FASTQ_DIR/${acc}.fastq"
    log "  $acc done."
  done
}

# ─── Step 5: FastQC on raw reads ─────────────────────────────────────────────
fastqc_raw() {
  log "Running FastQC on raw reads..."
  fastqc -t "$THREADS" -o "$QC_DIR/raw" "$FASTQ_DIR"/*.fastq.gz \
    >> "$LOG_DIR/fastqc_raw.log" 2>&1
  log "Raw FastQC done."
}

# ─── Step 6: Trim adapters ───────────────────────────────────────────────────
trim_reads() {
  log "Trimming adapters with Trim Galore..."
  for acc in "${ACCESSIONS[@]}"; do
    local trimmed="$TRIM_DIR/${acc}_trimmed.fq.gz"
    if [[ -f "$trimmed" ]]; then
      log "  $acc: already trimmed, skipping."
      continue
    fi
    log "  Trimming $acc..."
    trim_galore \
      --quality 20 \
      --length 20 \
      --cores 4 \
      --gzip \
      --fastqc \
      --fastqc_args "--outdir $QC_DIR/trimmed" \
      --output_dir "$TRIM_DIR" \
      "$FASTQ_DIR/${acc}.fastq.gz" \
      >> "$LOG_DIR/trim_${acc}.log" 2>&1
    log "  $acc trimmed."
  done
}

# ─── Step 7: Align with HISAT2 ───────────────────────────────────────────────
align_reads() {
  log "Aligning reads with HISAT2..."
  for acc in "${ACCESSIONS[@]}"; do
    local bam="$ALIGN_DIR/${acc}.sorted.bam"
    if [[ -f "$bam" ]]; then
      log "  $acc: already aligned, skipping."
      continue
    fi
    log "  Aligning $acc..."
    hisat2 \
      -x "$INDEX_PREFIX" \
      -U "$TRIM_DIR/${acc}_trimmed.fq.gz" \
      --dta \
      -p "$THREADS" \
      --rna-strandness R \
      --summary-file "$LOG_DIR/hisat2_${acc}.log" \
      2>> "$LOG_DIR/hisat2_${acc}.log" \
    | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"
    log "  $acc aligned: $(samtools flagstat "$bam" | grep 'mapped (' | head -1)"
  done
}

# ─── Step 8: Count reads with featureCounts ──────────────────────────────────
count_reads() {
  log "Counting reads with featureCounts..."
  local bam_files=()
  for acc in "${ACCESSIONS[@]}"; do
    bam_files+=("$ALIGN_DIR/${acc}.sorted.bam")
  done

  featureCounts \
    -T "$THREADS" \
    -a "$GTF" \
    -o "$COUNTS_DIR/counts_matrix.txt" \
    -t exon \
    -g gene_id \
    -s 2 \
    "${bam_files[@]}" \
    > "$LOG_DIR/featurecounts.log" 2>&1

  # Clean up column names to just accession IDs
  python3 - <<'PYEOF'
import pandas as pd, re, os

counts_file = "/home/prithvirajgavai/rnaseq_pipeline/counts/counts_matrix.txt"
df = pd.read_csv(counts_file, sep="\t", comment="#")
df.columns = [re.sub(r".*/(.+)\.sorted\.bam", r"\1", c) for c in df.columns]
df.to_csv(counts_file.replace(".txt", "_clean.txt"), sep="\t", index=False)
print(f"Count matrix saved: {df.shape[0]} genes x {df.shape[1]-5} samples")
PYEOF

  log "featureCounts done. Matrix: $COUNTS_DIR/counts_matrix_clean.txt"
}

# ─── Step 9: MultiQC ─────────────────────────────────────────────────────────
run_multiqc() {
  log "Running MultiQC..."
  multiqc "$BASE_DIR" \
    --outdir "$QC_DIR/multiqc" \
    --filename "SRP439187_multiqc_report" \
    --force \
    >> "$LOG_DIR/multiqc.log" 2>&1
  log "MultiQC report: $QC_DIR/multiqc/SRP439187_multiqc_report.html"
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
  log "=========================================="
  log " RNA-seq Pipeline: SRP439187"
  log " Arabidopsis thaliana | HISAT2 + featureCounts"
  log "=========================================="

  check_tools
  setup_dirs
  download_reference
  build_index
  sra_to_fastq
  fastqc_raw
  trim_reads
  align_reads
  count_reads
  run_multiqc

  log "=========================================="
  log " Pipeline complete!"
  log " Counts matrix : $COUNTS_DIR/counts_matrix_clean.txt"
  log " MultiQC report: $QC_DIR/multiqc/SRP439187_multiqc_report.html"
  log " Next step     : Rscript deseq2_analysis.R"
  log "=========================================="
}

# Run only the step passed as argument, or full pipeline
case "${1:-all}" in
  all)            main ;;
  check)          check_tools ;;
  dirs)           setup_dirs ;;
  reference)      setup_dirs && download_reference ;;
  index)          build_index ;;
  fastq)          sra_to_fastq ;;
  fastqc_raw)     fastqc_raw ;;
  trim)           trim_reads ;;
  align)          align_reads ;;
  counts)         count_reads ;;
  multiqc)        run_multiqc ;;
  *) echo "Usage: $0 [all|check|reference|index|fastq|fastqc_raw|trim|align|counts|multiqc]"; exit 1 ;;
esac
