#!/bin/bash
# Install all RNA-seq pipeline dependencies via conda

set -euo pipefail

echo "Installing RNA-seq pipeline tools..."

# Check if conda is available
if ! command -v conda &>/dev/null; then
  echo "conda not found. Installing Miniconda..."
  wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
  eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
  conda init bash
  echo "Miniconda installed. Run: source ~/.bashrc, then re-run this script."
  exit 0
fi

# Create dedicated environment
conda create -y -n rnaseq -c conda-forge -c bioconda \
  python=3.11 \
  hisat2 \
  samtools \
  fastqc \
  trim-galore \
  subread \
  multiqc \
  pandas \
  r-base \
  bioconductor-deseq2 \
  bioconductor-enhancedvolcano \
  r-pheatmap \
  r-ggplot2 \
  r-dplyr \
  r-readr

echo ""
echo "Installation complete!"
echo "Activate with: conda activate rnaseq"
echo "Then run:      bash rnaseq_pipeline.sh"
