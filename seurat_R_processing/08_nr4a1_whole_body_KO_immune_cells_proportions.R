# =============================================================================
# T cell GOBP Pathway Scoring + Gene Expression Pipeline
# For Nr4a1 KO scRNA-seq analysis
# Built to extend tcell_scoring_pipeline_v7.R
#
# NEW FUNCTIONALITY:
#   1. Subset Seurat object to T cells only
#   2. GOBP pathway scoring (26 pathways) via msigdbr + AUCell
#   3. plot_expression_custom() — bar and violin plots for any gene
#   4. Gene loops: collaborator_genes + collaborator_genes_2
#      → bar plots, violin plots, summary CSVs
#   5. Pathway score bar + violin plots for all scored pathways
#   6. Pseudo-bulk (per-sample) pathway score CSV + t-test table
#      for publication-quality statistics
# =============================================================================

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(AUCell)
library(msigdbr)
library(tidyverse)

set.seed(42)

# =============================================================================
# 0. PATHS AND SETUP
# =============================================================================

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)

# Output directories
dir.create("plots/pathway_scores",                               showWarnings = FALSE, recursive = TRUE)
dir.create("plots/gene_expression/collaborator_genes",          showWarnings = FALSE, recursive = TRUE)
dir.create("plots/gene_expression/collaborator_genes_2",        showWarnings = FALSE, recursive = TRUE)

# Load Seurat object if not already in environment
if (!exists("data")) {
  message("Loading Seurat object...")
  data <- readRDS(file.path(base_dir,
                            "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4_exonic.rds"))
  DefaultAssay(data) <- "RNA"
  data$Genotype <- data$Condition
}

# Shared constants — must match v7 exactly
tcell_subtypes <- c("CD4+ T cells", "CD8+ T cells", "Treg cells",
                    "DN T cells", "γδ T cells", "Cyc. T cells")
geno_levels    <- c("Nr4a1 KO", "WT")
geno_colors    <- c("Nr4a1 KO" = "#F8766D", "WT" = "#00BFC4")

# Shared ggplot theme (mirrors v7 cust_theme)
cust_theme <- theme_classic() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.title.x     = element_blank(),
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 14),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 12)
  )

# =============================================================================
# 1. SUBSET TO T CELLS
# =============================================================================

message("Subsetting to T cells...")
tcell_obj <- subset(data, subset = CellType == "T cells")
DefaultAssay(tcell_obj) <- "RNA"

tcell_obj$Genotype    <- factor(tcell_obj$Genotype,    levels = geno_levels)
tcell_obj$SubCellType <- factor(tcell_obj$SubCellType, levels = tcell_subtypes)

message(sprintf("T cell object: %d cells, %d samples",
                ncol(tcell_obj),
                length(unique(tcell_obj$BatchID))))

# =============================================================================
# 2. GOBP PATHWAY GENE SETS (26 pathways)
# =============================================================================

message("Fetching GOBP gene sets from msigdbr (Mus musculus, C5/BP)...")
gobp_all <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

pathway_names <- c(
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_ACTIVATION_VIA_T_CELL_RECEPTOR_CONTACT_WITH_ANTIGEN_BOUND_TO_MHC_MOLECULE_ON_ANTIGEN_PRESENTING_CELL",
  "GOBP_ACTIVATED_T_CELL_PROLIFERATION",
  "GOBP_T_CELL_PROLIFERATION",
  "GOBP_T_CELL_PROLIFERATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION",
  "GOBP_T_CELL_DIFFERENTIATION",
  "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_T_CELL_HOMEOSTASIS",
  "GOBP_T_CELL_MEDIATED_IMMUNITY",
  "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
  "GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",
  "GOBP_T_CELL_CYTOKINE_PRODUCTION",
  "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_T_CELL_TOLERANCE_INDUCTION",
  "GOBP_T_CELL_APOPTOTIC_PROCESS",
  "GOBP_ACTIVATION_INDUCED_CELL_DEATH_OF_T_CELLS",
  "GOBP_T_CELL_MIGRATION",
  "GOBP_T_CELL_EXTRAVASATION",
  "GOBP_T_CELL_CHEMOTAXIS",
  "GOBP_T_CELL_SELECTION",
  "GOBP_T_CELL_LINEAGE_COMMITMENT",
  "GOBP_T_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

# Build list: pathway → intersected gene symbols
pathway_gene_sets <- lapply(pathway_names, function(pw) {
  genes <- gobp_all %>%
    filter(gs_name == pw) %>%
    pull(gene_symbol) %>%
    unique()
  intersect(genes, rownames(tcell_obj))
})
names(pathway_gene_sets) <- pathway_names

# Report sizes and filter
n_genes <- sapply(pathway_gene_sets, length)
message("\nGene set sizes (expressed in data):")
print(data.frame(
  pathway = gsub("GOBP_", "", names(n_genes)),
  n_genes = n_genes
), row.names = FALSE)

pathway_gene_sets_ok <- pathway_gene_sets[n_genes >= 3]
message(sprintf("\n%d of %d pathways have >= 20 genes and will be scored.",
                length(pathway_gene_sets_ok), length(pathway_gene_sets)))

# =============================================================================
# 3. AUCELL SCORING
# =============================================================================

message("\nBuilding AUCell rankings (may take several minutes)...")
cells_rankings <- AUCell_buildRankings(
  GetAssayData(tcell_obj, layer = "data"),   # Seurat v5 syntax
  verbose = FALSE,
  plot    = FALSE
)

message("Calculating AUC for all pathways...")
cells_AUC <- AUCell_calcAUC(
  pathway_gene_sets_ok,
  cells_rankings,
  aucMaxRank = nrow(cells_rankings) * 0.05,   # top 5% ranked genes
  verbose    = FALSE
)
# Store scores in Seurat metadata
auc_mat <- getAUC(cells_AUC)
for (pw in rownames(auc_mat)) {
  col_name <- paste0("score_", gsub("GOBP_", "", pw))
  tcell_obj[[col_name]] <- as.numeric(auc_mat[pw, ])
}
message("AUCell scores added to tcell_obj@meta.data.")


# =============================================================================
# 4. plot_expression_custom()
# =============================================================================
plot_expression_custom <- function(seurat_obj,
                                   gene,
                                   plot_type     = "bar",
                                   group_by      = "SubCellType",
                                   condition_col = "Genotype",
                                   conditions    = geno_levels,
                                   assay         = "RNA",
                                   layer         = "data",
                                   save_path     = NULL) {
  
  if (!gene %in% rownames(seurat_obj)) {
    warning("Gene not found: ", gene)
    return(NULL)
  }
  
  expr_vec <- as.numeric(
    GetAssayData(seurat_obj, assay = assay, layer = layer)[gene, ]
  )
  
  plot_df <- seurat_obj@meta.data %>%
    dplyr::mutate(
      expr      = expr_vec,
      group     = .data[[group_by]],
      condition = .data[[condition_col]]
    ) %>%
    dplyr::filter(group %in% tcell_subtypes, condition %in% conditions) %>%
    dplyr::mutate(
      group     = factor(group,     levels = tcell_subtypes),
      condition = factor(condition, levels = conditions)
    )
  
  if (nrow(plot_df) == 0) {
    warning("No cells after filtering for: ", gene)
    return(NULL)
  }
  
  # Helper: save with dpi=300 for raster formats, plain ggsave for vector
  save_plot <- function(plot, path, w, h) {
    if (grepl("\\.png$", path, ignore.case = TRUE)) {
      ggsave(path, plot, width = w, height = h, dpi = 300)
    } else {
      ggsave(path, plot, width = w, height = h)
    }
    message("  Saved plot: ", basename(path))
  }
  
  comparisons <- list(conditions)
  
  # --- BAR ---
  if (plot_type == "bar") {
    
    p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition)) +
      geom_bar(stat = "summary", fun = "mean",
               position = position_dodge(0.9), color = "black", width = 0.7) +
      geom_errorbar(stat = "summary", fun.data = mean_se,
                    position = position_dodge(0.9),
                    width = 0.25, color = "black") +
      geom_jitter(aes(fill = condition), shape = 21,
                  color = "black", stroke = 0.4,
                  width = 0.12, size = 1.8, alpha = 0.75) +
      stat_compare_means(method = "wilcox.test", label = "p.signif",
                         comparisons = comparisons) +
      facet_wrap(~ group, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = geno_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
      labs(title = bquote(italic(.(gene))),
           y     = "Mean normalized expression") +
      cust_theme + NoLegend()
    
    summary_data <- plot_df %>%
      dplyr::group_by(group, condition) %>%
      dplyr::summarise(
        Mean_expression = mean(expr, na.rm = TRUE),
        SE_expression   = sd(expr, na.rm = TRUE) / sqrt(dplyr::n()),
        n_cells         = dplyr::n(),
        .groups         = "drop"
      ) %>%
      dplyr::rename(SubCellType = group, Genotype = condition)
    
    if (!is.null(save_path)) {
      save_plot(p, save_path, w = 14, h = 4)
      
      csv_path <- sub("\\.[^.]+$", "_summary.csv", save_path)
      write.csv(summary_data, csv_path, row.names = FALSE)
      message("  Saved CSV:  ", basename(csv_path))
    }
    
    return(list(plot = p, summary_data = summary_data))
    
    # --- VIOLIN ---
  } else if (plot_type == "violin") {
    
    p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
      geom_boxplot(width = 0.10, outlier.size = 0.3,
                   color = "grey40", linewidth = 0.4) +
      stat_compare_means(method = "wilcox.test", label = "p.signif",
                         comparisons = comparisons) +
      facet_wrap(~ group, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = geno_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
      labs(title    = bquote(italic(.(gene)) ~ "(violin)"),
           subtitle = "Wilcoxon rank-sum; whiskers = 1.5x IQR",
           x        = "Genotype",
           y        = "Normalized expression") +
      cust_theme + NoLegend()
    
    if (!is.null(save_path)) {
      save_plot(p, save_path, w = 14, h = 4)
    }
    
    return(p)
    
  } else {
    stop("plot_type must be 'bar' or 'violin'")
  }
}

# =============================================================================
# 5. GENE EXPRESSION LOOPS
# =============================================================================

# --- collaborator_genes ---
collaborator_genes <- c(
  "Bcl2",  "Cd38",  "Cd69",  "Ctla4", "Ifng",  "Il10",  "Il12a",
  "Il15",  "Il17a", "Il21",  "Il23a", "Il27",  "Il2",   "Il2ra",
  "Il4",   "Il6",   "Il7",   "Mki67", "Nkg7",  "Pdcd1", "Prf1",
  "Stat5a","Tgfb1", "Tigit", "Tnf",   "Top2a", "Tox"
)

message("\n=== collaborator_genes: bar + violin + CSV ===")
for (gene in collaborator_genes) {
  message("  ", gene)
  out_dir <- "plots/gene_expression/collaborator_genes"
  
  bar_res <- plot_expression_custom(
    tcell_obj, gene, plot_type = "bar",
    save_path = file.path(out_dir, paste0(gene, "_BarPlot_by_SubCellType.png"))
  )
  plot_expression_custom(
    tcell_obj, gene, plot_type = "violin",
    save_path = file.path(out_dir, paste0(gene, "_Violin_by_SubCellType.png"))
  )
  if (!is.null(bar_res))
    write.csv(bar_res$summary_data,
              file.path(out_dir, paste0(gene, "_summary_statistics_tcells.csv")),
              row.names = FALSE)
}

plot_expression_custom(
  tcell_obj, "Gzmb", plot_type = "bar")
  
# --- collaborator_genes_2 ---
collaborator_genes_2 <- c(
  "Ccl3",   "Ccl4",   "Ccl5",   "Egr1",   "Egr2",   "Eomes",
  "Fos",    "Gzmk",   "Havcr2", "Icos",   "Irf4",   "Jun",
  "Klrd1",  "Klrk1",  "Lag3",   "Lamp1",  "Myc",    "Nfkbia",
  "Rab27a", "Stat1",  "Stat4",  "Tbx21",  "Tnfrsf4","Tnfrsf9",
  "Zeb2"
)

message("\n=== collaborator_genes_2: bar + violin + CSV ===")
for (gene in collaborator_genes_2) {
  message("  ", gene)
  out_dir <- "plots/gene_expression/collaborator_genes_2"
  
  bar_res <- plot_expression_custom(
    tcell_obj, gene, plot_type = "bar",
    save_path = file.path(out_dir, paste0(gene, "_BarPlot_by_SubCellType.png"))
  )
  plot_expression_custom(
    tcell_obj, gene, plot_type = "violin",
    save_path = file.path(out_dir, paste0(gene, "_Violin_by_SubCellType.png"))
  )
  if (!is.null(bar_res))
    write.csv(bar_res$summary_data,
              file.path(out_dir, paste0(gene, "_summary_statistics.csv")),
              row.names = FALSE)
}

# =============================================================================
# 6. PATHWAY SCORE PLOTS
# =============================================================================

#' Plot an AUCell pathway score (bar or violin) across T cell subtypes
#'
#' @param seurat_obj    Seurat object with score in metadata
#' @param score_col     Metadata column name for the score
#' @param pathway_label Short human-readable pathway name for plot title
#' @param plot_type     "bar" or "violin"
#' @param save_path     Optional SVG save path

plot_pathway_score <- function(seurat_obj, score_col, pathway_label,
                               plot_type = "bar", save_path = NULL) {
  
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    warning("Score column not found: ", score_col)
    return(NULL)
  }
  
  plot_df <- seurat_obj@meta.data %>%
    dplyr::filter(SubCellType %in% tcell_subtypes, Genotype %in% geno_levels) %>%
    dplyr::mutate(
      score       = .data[[score_col]],
      SubCellType = factor(SubCellType, levels = tcell_subtypes),
      Genotype    = factor(Genotype,    levels = geno_levels)
    )
  
  comparisons <- list(geno_levels)
  
  if (plot_type == "bar") {
    p <- ggplot(plot_df, aes(x = Genotype, y = score, fill = Genotype)) +
      geom_bar(stat = "summary", fun = "mean",
               position = position_dodge(0.9), color = "black", width = 0.7) +
      geom_errorbar(stat = "summary", fun.data = mean_se,
                    position = position_dodge(0.9),
                    width = 0.25, color = "black") +
      geom_jitter(aes(fill = Genotype), shape = 21,
                  color = "black", stroke = 0.4,
                  width = 0.12, size = 1.8, alpha = 0.75) +
      stat_compare_means(method = "wilcox.test", label = "p.signif",
                         comparisons = comparisons) +
      facet_wrap(~ SubCellType, nrow = 1, scales = "free_y") +
      scale_fill_manual(values  = geno_colors) +
      scale_color_manual(values = geno_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
      labs(title    = pathway_label,
           subtitle = "AUCell score per cell; Wilcoxon rank-sum",
           y        = "AUCell score") +
      cust_theme + NoLegend()
    
  } else {
    p <- ggplot(plot_df, aes(x = Genotype, y = score, fill = Genotype)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
      geom_boxplot(width = 0.10, outlier.size = 0.3,
                   color = "grey40", linewidth = 0.4) +
      stat_compare_means(method = "wilcox.test", label = "p.signif",
                         comparisons = comparisons) +
      facet_wrap(~ SubCellType, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = geno_colors) +
      labs(title    = paste0(pathway_label, " (violin)"),
           subtitle = "Per-cell AUCell distribution; whiskers = 1.5 × IQR",
           x        = "Genotype",
           y        = "AUCell score") +
      cust_theme + NoLegend()
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 16, height = 4)
    message("  Saved: ", basename(save_path))
  }
  return(p)
}

message("\n=== Pathway score plots (bar + violin) ===")
for (pw in names(pathway_gene_sets_ok)) {
  score_col   <- paste0("score_", gsub("GOBP_", "", pw))
  short_label <- gsub("GOBP_", "", pw) %>% gsub("_", " ", .) %>% tools::toTitleCase()
  base_name   <- gsub("GOBP_", "", pw)
  
  tryCatch({
    plot_pathway_score(
      tcell_obj, score_col, short_label, "bar",
      file.path("plots/pathway_scores", paste0("score_", base_name, "_barplot.png"))
    )
    plot_pathway_score(
      tcell_obj, score_col, short_label, "violin",
      file.path("plots/pathway_scores", paste0("score_", base_name, "_violin.png"))
    )
  }, error = function(e) warning("Error for ", pw, ": ", e$message))
}

# =============================================================================
# 7. PSEUDO-BULK PATHWAY SCORES (sample-level; use for publication stats)
# =============================================================================
message("\n=== Pseudo-bulk pathway scores (per sample × subtype) ===")
score_cols <- paste0("score_", gsub("GOBP_", "", names(pathway_gene_sets_ok)))
pseudobulk_df <- tcell_obj@meta.data %>%
  dplyr::filter(SubCellType %in% tcell_subtypes, Genotype %in% geno_levels) %>%
  dplyr::group_by(BatchID, Genotype, SubCellType) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

write.csv(pseudobulk_df, "pathway_pseudobulk_scores_per_sample.csv", row.names = FALSE)
message("Saved: pathway_pseudobulk_scores_per_sample.csv")

# Sample-level t-tests
ttest_rows <- list()
for (sc in score_cols) {
  for (ct in tcell_subtypes) {
    sub <- dplyr::filter(pseudobulk_df, SubCellType == ct)
    ko  <- sub[sub$Genotype == "Nr4a1 KO", sc, drop = TRUE]
    wt  <- sub[sub$Genotype == "WT",       sc, drop = TRUE]
    if (length(ko) < 2 || length(wt) < 2) next
    res <- tryCatch(t.test(ko, wt), error = function(e) NULL)
    if (is.null(res)) next
    ttest_rows[[length(ttest_rows) + 1]] <- data.frame(
      pathway     = gsub("^score_", "", sc),
      SubCellType = ct,
      mean_KO     = round(mean(ko, na.rm = TRUE), 5),
      mean_WT     = round(mean(wt, na.rm = TRUE), 5),
      p_value     = round(res$p.value, 5),
      direction   = ifelse(mean(ko) > mean(wt), "KO > WT", "KO < WT"),
      sig         = dplyr::case_when(
        res$p.value < 0.001 ~ "***",
        res$p.value < 0.01  ~ "**",
        res$p.value < 0.05  ~ "*",
        TRUE                ~ "ns"),
      stringsAsFactors = FALSE
    )
  }
}

if (length(ttest_rows) > 0) {
  ttest_df <- dplyr::bind_rows(ttest_rows) %>% dplyr::arrange(p_value)
  write.csv(ttest_df, "pathway_pseudobulk_ttest_results.csv", row.names = FALSE)
  message("Saved: pathway_pseudobulk_ttest_results.csv")
  sig_hits <- dplyr::filter(ttest_df, sig != "ns")
  if (nrow(sig_hits) > 0) {
    message("\nSignificant pathways (pseudobulk t-test, p < 0.05):")
    print(sig_hits, row.names = FALSE)
  } else {
    message("No significant pseudobulk results at p < 0.05.")
    message("(This is common with small n; per-cell results remain exploratory.)")
  }
}

# =============================================================================
# SUMMARY
# =============================================================================
message("\n", strrep("=", 65))
message("Pipeline complete.")
message(sprintf("  Pathways scored:              %d", length(pathway_gene_sets_ok)))
message(sprintf("  collaborator_genes processed: %d genes × 2 plots + CSV",
                length(collaborator_genes)))
message(sprintf("  collaborator_genes_2 processed: %d genes × 2 plots + CSV",
                length(collaborator_genes_2)))
message("  Output CSVs:")
message("    pathway_pseudobulk_scores_per_sample.csv")
message("    pathway_pseudobulk_ttest_results.csv")
message("    plots/gene_expression/*/gene_summary_statistics*.csv")
message(strrep("=", 65))
message("\nIMPORTANT: Per-cell Wilcoxon tests on plots are exploratory only.")
message("Use pathway_pseudobulk_ttest_results.csv for publication statistics.")