#!/usr/bin/env Rscript
# GO & KEGG enrichment analysis for SRP439187
# Arabidopsis thaliana — 2x2 factorial design

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(readr)
})

BASE_DIR <- "/home/prithvirajgavai/rnaseq_pipeline"
DE_DIR   <- file.path(BASE_DIR, "deseq2_results")
OUT_DIR  <- file.path(BASE_DIR, "enrichment_results")
for (d in c("", "GO", "KEGG", "GSEA")) dir.create(file.path(OUT_DIR, d), showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────
run_go <- function(gene_list, label, ont = "BP") {
  res <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  if (is.null(res) || nrow(res) == 0) {
    cat(sprintf("  [%s] No GO-%s terms enriched.\n", label, ont)); return(NULL)
  }
  cat(sprintf("  [%s] GO-%s: %d enriched terms\n", label, ont, nrow(res)))
  write_tsv(as.data.frame(res),
    file.path(OUT_DIR, "GO", sprintf("%s_GO_%s.tsv", label, ont)))
  return(res)
}

run_kegg <- function(gene_list, label) {
  # Convert TAIR → ENTREZID
  eg <- bitr(gene_list, fromType = "TAIR", toType = "ENTREZID",
             OrgDb = org.At.tair.db, drop = TRUE)
  if (nrow(eg) == 0) { cat(sprintf("  [%s] No ENTREZ IDs mapped.\n", label)); return(NULL) }
  cat(sprintf("  [%s] Mapped %d/%d genes to ENTREZID\n", label, nrow(eg), length(gene_list)))

  res <- tryCatch(
    enrichKEGG(
      gene          = eg$ENTREZID,
      organism      = "ath",
      keyType       = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    ),
    error = function(e) { cat(sprintf("  [%s] KEGG error: %s\n", label, conditionMessage(e))); NULL }
  )
  if (is.null(res) || nrow(res) == 0) {
    cat(sprintf("  [%s] No KEGG pathways enriched.\n", label)); return(NULL)
  }
  cat(sprintf("  [%s] KEGG: %d pathways enriched\n", label, nrow(res)))
  write_tsv(as.data.frame(res),
    file.path(OUT_DIR, "KEGG", sprintf("%s_KEGG.tsv", label)))
  return(res)
}

save_go_plots <- function(go_res, label, ont) {
  prefix <- file.path(OUT_DIR, "GO", sprintf("%s_GO_%s", label, ont))
  title  <- sprintf("Contrast %s | GO-%s", label, ont)

  tryCatch({
    pdf(paste0(prefix, "_dotplot.pdf"), width = 10, height = 8)
    print(dotplot(go_res, showCategory = 20) + ggtitle(title))
    dev.off()
  }, error = function(e) cat("  dotplot error:", conditionMessage(e), "\n"))

  tryCatch({
    pdf(paste0(prefix, "_barplot.pdf"), width = 10, height = 8)
    print(barplot(go_res, showCategory = 20) + ggtitle(title))
    dev.off()
  }, error = function(e) cat("  barplot error:", conditionMessage(e), "\n"))

  tryCatch({
    go_sim <- pairwise_termsim(go_res)
    pdf(paste0(prefix, "_emapplot.pdf"), width = 12, height = 10)
    print(emapplot(go_sim, showCategory = 30) + ggtitle(paste(title, "— enrichment map")))
    dev.off()
  }, error = function(e) cat("  emapplot error:", conditionMessage(e), "\n"))

  tryCatch({
    pdf(paste0(prefix, "_upsetplot.pdf"), width = 12, height = 6)
    print(upsetplot(go_res, n = 10))
    dev.off()
  }, error = function(e) cat("  upsetplot error:", conditionMessage(e), "\n"))
}

save_kegg_plots <- function(kegg_res, label) {
  prefix <- file.path(OUT_DIR, "KEGG", sprintf("%s_KEGG", label))
  title  <- sprintf("Contrast %s | KEGG Pathways", label)

  tryCatch({
    pdf(paste0(prefix, "_dotplot.pdf"), width = 10, height = 7)
    print(dotplot(kegg_res, showCategory = 20) + ggtitle(title))
    dev.off()
  }, error = function(e) cat("  KEGG dotplot error:", conditionMessage(e), "\n"))

  tryCatch({
    pdf(paste0(prefix, "_barplot.pdf"), width = 10, height = 7)
    print(barplot(kegg_res, showCategory = 20) + ggtitle(title))
    dev.off()
  }, error = function(e) cat("  KEGG barplot error:", conditionMessage(e), "\n"))
}

# ── GSEA (full ranked gene list) ──────────────────────────────────────────────
run_gsea <- function(all_df, label) {
  min_nonzero_p <- min(all_df$pvalue[all_df$pvalue > 0], na.rm = TRUE)
  max_stat <- -log10(min_nonzero_p) * 1.1  # cap for zero p-values

  ranked <- all_df %>%
    filter(!is.na(log2FoldChange), !is.na(pvalue), is.finite(log2FoldChange)) %>%
    mutate(
      pvalue_capped = pmax(pvalue, min_nonzero_p * 0.1),
      rank_stat     = pmin(pmax(-log10(pvalue_capped) * sign(log2FoldChange),
                                -max_stat), max_stat)
    ) %>%
    arrange(desc(rank_stat))

  gene_rank <- setNames(ranked$rank_stat, ranked$gene_id)
  gene_rank <- gene_rank[!duplicated(names(gene_rank)) & is.finite(gene_rank)]

  res <- tryCatch(
    gseGO(
      geneList      = gene_rank,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = "BP",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      eps           = 0
    ),
    error = function(e) { cat(sprintf("  [%s] GSEA error: %s\n", label, conditionMessage(e))); NULL }
  )

  if (is.null(res) || nrow(res) == 0) {
    cat(sprintf("  [%s] No GSEA GO-BP terms significant.\n", label)); return(NULL)
  }
  cat(sprintf("  [%s] GSEA GO-BP: %d terms\n", label, nrow(res)))
  write_tsv(as.data.frame(res),
    file.path(OUT_DIR, "GSEA", sprintf("%s_GSEA_GOBP.tsv", label)))

  tryCatch({
    pdf(file.path(OUT_DIR, "GSEA", sprintf("%s_GSEA_dotplot.pdf", label)),
        width = 10, height = 8)
    print(dotplot(res, showCategory = 20, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle(sprintf("Contrast %s | GSEA GO-BP", label)))
    dev.off()
  }, error = function(e) cat("  GSEA dotplot error:", conditionMessage(e), "\n"))

  tryCatch({
    pdf(file.path(OUT_DIR, "GSEA", sprintf("%s_GSEA_ridgeplot.pdf", label)),
        width = 10, height = 8)
    print(ridgeplot(res, showCategory = 20) +
          ggtitle(sprintf("Contrast %s | GSEA GO-BP — enrichment score distribution", label)))
    dev.off()
  }, error = function(e) cat("  GSEA ridgeplot error:", conditionMessage(e), "\n"))

  return(res)
}

# ── Cross-contrast comparison ─────────────────────────────────────────────────
compare_contrasts <- function(contrast_gene_lists) {
  cat("\nRunning compareCluster across all contrasts...\n")
  # Limit to reasonably sized gene lists for speed
  trimmed <- lapply(contrast_gene_lists, function(g) {
    if (length(g) > 2000) sample(g, 2000) else g
  })

  res <- tryCatch(
    compareCluster(
      geneClusters  = trimmed,
      fun           = "enrichGO",
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) { cat("  compareCluster error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) return(NULL)

  write_tsv(as.data.frame(res),
    file.path(OUT_DIR, "GO", "compareCluster_all_contrasts_GOBP.tsv"))

  tryCatch({
    pdf(file.path(OUT_DIR, "GO", "compareCluster_dotplot.pdf"), width = 14, height = 12)
    print(dotplot(res, showCategory = 10) +
          ggtitle("GO-BP enrichment across all 4 contrasts"))
    dev.off()
    cat(sprintf("  compareCluster: %d enriched terms\n", nrow(as.data.frame(res))))
  }, error = function(e) cat("  compareCluster dotplot error:", conditionMessage(e), "\n"))

  return(res)
}

# ── Load significant gene lists ───────────────────────────────────────────────
contrasts <- list(
  A = "A_heat_effect_in_control_line",
  B = "B_line_effect_under_control",
  C = "C_interaction_line_x_treatment",
  D = "D_HeatSelected_Heat_vs_ControlLine_Control"
)

all_sig_genes <- list()

for (key in names(contrasts)) {
  label <- contrasts[[key]]
  cat(sprintf("\n======== Contrast %s: %s ========\n", key, label))

  sig_df    <- read_tsv(file.path(DE_DIR, paste0(label, "_significant.tsv")),
                        show_col_types = FALSE)
  sig_genes <- sig_df$gene_id
  all_sig_genes[[key]] <- sig_genes
  cat(sprintf("  Significant genes: %d\n", length(sig_genes)))

  if (length(sig_genes) < 10) { cat("  Too few genes, skipping.\n"); next }

  for (ont in c("BP", "MF", "CC")) {
    go_res <- run_go(sig_genes, key, ont)
    if (!is.null(go_res)) save_go_plots(go_res, key, ont)
  }

  kegg_res <- run_kegg(sig_genes, key)
  if (!is.null(kegg_res)) save_kegg_plots(kegg_res, key)

  all_df <- read_tsv(file.path(DE_DIR, paste0(label, "_all.tsv")),
                     show_col_types = FALSE)
  run_gsea(all_df, key)
}

compare_contrasts(all_sig_genes)

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n\n=== Enrichment Analysis Summary ===\n")
for (key in names(contrasts)) {
  for (ont in c("BP", "MF", "CC")) {
    f <- file.path(OUT_DIR, "GO", sprintf("%s_GO_%s.tsv", key, ont))
    if (file.exists(f)) cat(sprintf("  Contrast %s GO-%-2s: %d terms\n", key, ont,
                                    nrow(read_tsv(f, show_col_types = FALSE))))
  }
  f <- file.path(OUT_DIR, "KEGG", sprintf("%s_KEGG.tsv", key))
  if (file.exists(f)) cat(sprintf("  Contrast %s KEGG  : %d pathways\n", key,
                                  nrow(read_tsv(f, show_col_types = FALSE))))
  f <- file.path(OUT_DIR, "GSEA", sprintf("%s_GSEA_GOBP.tsv", key))
  if (file.exists(f)) cat(sprintf("  Contrast %s GSEA  : %d GO-BP terms\n", key,
                                  nrow(read_tsv(f, show_col_types = FALSE))))
}
cat("\nOutputs in:", OUT_DIR, "\n")
