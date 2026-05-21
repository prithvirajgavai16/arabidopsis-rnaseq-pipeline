#!/usr/bin/env Rscript
# DESeq2 differential expression analysis for SRP439187
# Arabidopsis thaliana — 2x2 factorial design
# Factor 1: line (heat_selected / control_line)
# Factor 2: treatment (heat / control)

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

BASE_DIR  <- "/home/prithvirajgavai/rnaseq_pipeline"
COUNTS_FILE <- file.path(BASE_DIR, "counts/counts_matrix_clean.txt")
OUT_DIR   <- file.path(BASE_DIR, "deseq2_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load count matrix ──────────────────────────────────────────────────────
cat("Loading count matrix...\n")
counts_raw <- read_tsv(COUNTS_FILE, comment = "#", show_col_types = FALSE)
gene_ids <- counts_raw$Geneid
count_matrix <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(count_matrix) <- gene_ids
storage.mode(count_matrix) <- "integer"
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(count_matrix), ncol(count_matrix)))

# ── 2. Sample metadata (2x2 factorial) ───────────────────────────────────────
# line:      heat_selected = evolved under heat stress | control_line = wild-type
# treatment: heat = grown under heat | control = grown under normal conditions
coldata <- data.frame(
  line = factor(c(
    "heat_selected", # SRR24710315
    "heat_selected", # SRR24710316
    "heat_selected", # SRR24710317
    "heat_selected", # SRR24710318
    "heat_selected", # SRR24710319
    "heat_selected", # SRR24710320
    "heat_selected", # SRR24710321
    "control_line",  # SRR24710322
    "control_line",  # SRR24710323
    "control_line",  # SRR24710324
    "control_line",  # SRR24710325
    "control_line",  # SRR24710326
    "control_line"   # SRR24710327
  ), levels = c("control_line", "heat_selected")),
  treatment = factor(c(
    "heat",    # SRR24710315
    "heat",    # SRR24710316
    "heat",    # SRR24710317
    "control", # SRR24710318
    "control", # SRR24710319
    "control", # SRR24710320
    "control", # SRR24710321
    "heat",    # SRR24710322
    "heat",    # SRR24710323
    "heat",    # SRR24710324
    "control", # SRR24710325
    "control", # SRR24710326
    "control"  # SRR24710327
  ), levels = c("control", "heat")),
  row.names = colnames(count_matrix)
)
# Combined group label for easy contrasts
coldata$group <- factor(paste(coldata$line, coldata$treatment, sep = "_"))
cat("\nSample metadata:\n")
print(coldata)

# ── 3. DESeqDataSet with full factorial design ────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ line + treatment + line:treatment
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat(sprintf("\nAfter pre-filtering: %d genes retained\n", nrow(dds)))

# ── 4. Run DESeq2 ─────────────────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("Size factors:\n"); print(sizeFactors(dds))

# ── 5. Extract results for 4 key contrasts ───────────────────────────────────
save_results <- function(res, label) {
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)
  sig <- df %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
  write_tsv(df,  file.path(OUT_DIR, paste0(label, "_all.tsv")))
  write_tsv(sig, file.path(OUT_DIR, paste0(label, "_significant.tsv")))
  cat(sprintf("\n[%s] DE genes (padj<0.05, |LFC|>=1): %d (up=%d, down=%d)\n",
    label, nrow(sig),
    sum(sig$log2FoldChange > 0, na.rm = TRUE),
    sum(sig$log2FoldChange < 0, na.rm = TRUE)))
  return(df)
}

# Contrast A: Effect of heat treatment in control line (CL_heat vs CL_control)
resA <- results(dds, contrast = c("treatment", "heat", "control"), alpha = 0.05)
dfA  <- save_results(resA, "A_heat_effect_in_control_line")

# Contrast B: Effect of line under control conditions (HS_ctrl vs CL_ctrl)
resB <- results(dds, contrast = c("line", "heat_selected", "control_line"), alpha = 0.05)
dfB  <- save_results(resB, "B_line_effect_under_control")

# Contrast C: Effect of heat in heat-selected line vs control line (interaction)
# interaction term = lineheat_selected.treatmentheat
resC <- results(dds, name = "lineheat_selected.treatmentheat", alpha = 0.05)
dfC  <- save_results(resC, "C_interaction_line_x_treatment")

# Contrast D: Heat-selected+heat vs Control+control (overall effect)
dds_grp <- dds
design(dds_grp) <- ~ group
dds_grp <- DESeq(dds_grp)
resD <- results(dds_grp,
  contrast = c("group", "heat_selected_heat", "control_line_control"),
  alpha = 0.05)
dfD <- save_results(resD, "D_HeatSelected_Heat_vs_ControlLine_Control")

# ── 6. PCA plot ───────────────────────────────────────────────────────────────
cat("\nGenerating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("line", "treatment"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2,
    color = treatment, shape = line, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "#2166ac", "heat" = "#d6604d")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA — SRP439187 | Arabidopsis thaliana",
          subtitle = "Shape = line type, Colour = growth condition") +
  theme_bw(base_size = 13)
ggsave(file.path(OUT_DIR, "PCA_plot.pdf"), p_pca, width = 9, height = 7)

# ── 7. Volcano plots for each contrast ───────────────────────────────────────
cat("Generating volcano plots...\n")
make_volcano <- function(df, title, filename) {
  pdf(file.path(OUT_DIR, filename), width = 10, height = 8)
  print(EnhancedVolcano(df,
    lab            = df$gene_id,
    x              = "log2FoldChange",
    y              = "padj",
    title          = title,
    subtitle       = "SRP439187 | Arabidopsis thaliana",
    pCutoff        = 0.05,
    FCcutoff       = 1.0,
    pointSize      = 2,
    labSize        = 3,
    colAlpha       = 0.6,
    legendPosition = "right"
  ))
  dev.off()
}
make_volcano(dfA, "Heat effect in Control line (Heat vs Control)",
             "volcano_A_heat_effect_control_line.pdf")
make_volcano(dfB, "Line effect under Control conditions (HS vs CL)",
             "volcano_B_line_effect_control_conditions.pdf")
make_volcano(dfC, "Interaction: line × treatment",
             "volcano_C_interaction.pdf")
make_volcano(dfD, "Heat-selected+Heat vs Control+Control",
             "volcano_D_overall.pdf")

# ── 8. Heatmap — top 50 DE genes from interaction contrast ───────────────────
cat("Generating heatmap...\n")
sig_interaction <- read_tsv(file.path(OUT_DIR, "C_interaction_line_x_treatment_significant.tsv"),
                             show_col_types = FALSE)
top50 <- head(sig_interaction$gene_id, 50)

if (length(top50) >= 2) {
  mat     <- assay(vsd)[top50, ]
  mat     <- mat - rowMeans(mat)
  ann_col <- data.frame(
    line      = coldata$line,
    treatment = coldata$treatment,
    row.names = rownames(coldata)
  )
  ann_colors <- list(
    treatment = c(control = "#2166ac", heat = "#d6604d"),
    line      = c(control_line = "#999999", heat_selected = "#e69f00")
  )
  pdf(file.path(OUT_DIR, "heatmap_top50_interaction_genes.pdf"), width = 11, height = 13)
  pheatmap(mat,
    annotation_col  = ann_col,
    annotation_colors = ann_colors,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    show_rownames   = TRUE,
    fontsize_row    = 7,
    main            = "Top 50 Interaction Genes (VST, row-centered)"
  )
  dev.off()
} else {
  cat("Too few interaction genes — plotting top 50 from contrast D instead.\n")
  sig_overall <- read_tsv(file.path(OUT_DIR, "D_HeatSelected_Heat_vs_ControlLine_Control_significant.tsv"),
                           show_col_types = FALSE)
  top50 <- head(sig_overall$gene_id, 50)
  if (length(top50) >= 2) {
    mat     <- assay(vsd)[top50, ]
    mat     <- mat - rowMeans(mat)
    ann_col <- data.frame(line = coldata$line, treatment = coldata$treatment,
                          row.names = rownames(coldata))
    pdf(file.path(OUT_DIR, "heatmap_top50_overall_DE_genes.pdf"), width = 11, height = 13)
    pheatmap(mat, annotation_col = ann_col, cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, fontsize_row = 7,
             main = "Top 50 DE Genes — HS+Heat vs CL+Control (VST, row-centered)")
    dev.off()
  }
}

cat("\n=== DESeq2 analysis complete ===\n")
cat("Outputs in:", OUT_DIR, "\n\n")
cat("  Contrast A — Heat effect in control line:\n")
cat("    A_heat_effect_in_control_line_all.tsv / _significant.tsv\n")
cat("    volcano_A_heat_effect_control_line.pdf\n\n")
cat("  Contrast B — Line effect under control conditions:\n")
cat("    B_line_effect_under_control_all.tsv / _significant.tsv\n")
cat("    volcano_B_line_effect_control_conditions.pdf\n\n")
cat("  Contrast C — Interaction (line x treatment):\n")
cat("    C_interaction_line_x_treatment_all.tsv / _significant.tsv\n")
cat("    volcano_C_interaction.pdf\n")
cat("    heatmap_top50_interaction_genes.pdf\n\n")
cat("  Contrast D — HS+Heat vs CL+Control (overall):\n")
cat("    D_HeatSelected_Heat_vs_ControlLine_Control_all.tsv / _significant.tsv\n")
cat("    volcano_D_overall.pdf\n\n")
cat("  PCA_plot.pdf\n")
