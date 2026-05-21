# RNA-seq Pipeline: *Arabidopsis thaliana* Heat-Stress Study (SRP439187)

End-to-end RNA-seq pipeline for differential expression and functional enrichment analysis of a 2×2 factorial heat-stress experiment in *Arabidopsis thaliana*.

---

## Dataset

| Field | Detail |
|---|---|
| SRA Project | [SRP439187](https://www.ncbi.nlm.nih.gov/sra/?term=SRP439187) |
| Organism | *Arabidopsis thaliana* (TAIR10) |
| Design | 2×2 factorial — line (heat-selected / control) × treatment (heat / control) |
| Samples | 13 single-end 50 bp Illumina HiSeq 1500 runs (SRR24710315–SRR24710327) |

---

## Pipeline Overview

```
SRA → FASTQ → FastQC → Trim Galore → HISAT2 → featureCounts → DESeq2 → clusterProfiler
```

| Step | Tool | Script |
|---|---|---|
| 0. Dependency check | — | `rnaseq_pipeline.sh` |
| 1. SRA → FASTQ | fasterq-dump | `rnaseq_pipeline.sh` |
| 2. Quality control | FastQC + MultiQC | `rnaseq_pipeline.sh` |
| 3. Adapter trimming | Trim Galore | `rnaseq_pipeline.sh` |
| 4. Alignment | HISAT2 + samtools | `rnaseq_pipeline.sh` |
| 5. Read counting | featureCounts (subread) | `rnaseq_pipeline.sh` |
| 6. Differential expression | DESeq2 | `deseq2_analysis.R` |
| 7. Functional enrichment | clusterProfiler (GO + KEGG + GSEA) | `enrichment_analysis.R` |

---

## Requirements

Install all tools with the provided conda script:

```bash
bash install_tools.sh
conda activate rnaseq
```

**Tools installed:** HISAT2, samtools, FastQC, Trim Galore, featureCounts, MultiQC, Python 3, R, DESeq2, EnhancedVolcano, clusterProfiler, org.At.tair.db

---

## Usage

### Full pipeline (steps 1–5)

```bash
conda activate rnaseq
bash rnaseq_pipeline.sh          # run all steps
bash rnaseq_pipeline.sh align    # run a single step
```

Available step arguments: `check | reference | index | fastq | fastqc_raw | trim | align | counts | multiqc`

### Differential expression

```bash
Rscript deseq2_analysis.R
```

Outputs four contrasts to `deseq2_results/`:

| Label | Contrast |
|---|---|
| A | Heat effect in control line |
| B | Line effect under control conditions |
| C | Line × treatment interaction |
| D | Heat-selected+heat vs. control line+control (overall) |

### Functional enrichment

```bash
Rscript enrichment_analysis.R
```

Produces GO (BP/MF/CC), KEGG, and GSEA results in `enrichment_results/`.

---

## Repository Structure

```
rnaseq_pipeline/
├── rnaseq_pipeline.sh       # Main pipeline (HISAT2 → featureCounts)
├── deseq2_analysis.R        # DESeq2 differential expression
├── enrichment_analysis.R    # GO / KEGG / GSEA enrichment
├── install_tools.sh         # Conda environment setup
├── counts/                  # featureCounts output matrix
├── deseq2_results/          # DE tables + volcano / heatmap plots
└── enrichment_results/      # GO, KEGG, GSEA tables and plots
```

> **Not tracked (large data):** `fastq/`, `trimmed/`, `aligned/`, `reference/`, `qc/`, `logs/`

---

## Reference

- Genome & annotation: Ensembl Plants release 59, TAIR10
- URLs embedded in `rnaseq_pipeline.sh` for reproducibility
